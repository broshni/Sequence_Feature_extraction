#Streamlit front end: upload the input files, pick a reference genome, and run every pipeline
#step in order, watching each one finish. Run with:  streamlit run app.py
import argparse
import io
import os
import re
import tempfile

import pandas as pd
from sf_columns import read_csv
import streamlit as st

import ui_theme as ui
from find_reference import find_reference_genomes
from run_pipeline import pipeline_steps, run_step, StepFailed

st.set_page_config(page_title="SF Table Builder", page_icon=":material/biotech:",
                   layout="wide", initial_sidebar_state="expanded")
ui.apply_styles()


def save_upload(upload, directory):
    """Write an uploaded file to disk and return its path."""
    path = os.path.join(directory, upload.name)
    with open(path, 'wb') as handle:
        handle.write(upload.getbuffer())
    return path


def settings_from_form(**values):
    """pipeline_steps() reads plain attributes, so a namespace is all it needs."""
    return argparse.Namespace(**values)


def column_problem(upload, required, what, hint=""):
    """A message when an optional upload lacks the columns the pipeline needs, else None.
    Checked at upload time: otherwise the mistake only surfaces as an error several steps in."""
    if upload is None:
        return None
    try:
        columns = set(read_csv(io.BytesIO(upload.getvalue()), nrows=0).columns)
    except Exception:
        return f"{upload.name} could not be read as a CSV file."
    missing_columns = [c for c in required if c not in columns]
    if not missing_columns:
        return None
    return (f"{upload.name} is not {what}: it has no {', '.join(repr(c) for c in missing_columns)} "
            f"column. {hint}").strip()


def format_group(group, common=''):
    """One short line for a candidate reference. `common` is the text every candidate's name
    starts with ('Influenza A virus (A/'); dropping it keeps the part that tells them apart
    visible in a narrow dropdown."""
    name = group['name'][len(common):] if common and len(group['name']) > len(common) else group['name']
    first = group['accessions'][0]
    if group['records'] > 1:
        return f"{name} · {group['records']} segments · {first}…"
    return f"{name} · {first}"


def candidate_labels(groups):
    """{label: group}, with the shared name prefix removed from every label."""
    names = [g['name'] for g in groups]
    common = os.path.commonprefix(names)
    common = common[:common.rfind(' ') + 1] if ' ' in common else ''
    # Only strip when the names really differ: candidates that share one name (Zika's two
    # alternative genomes) are told apart by accession, and stripping would leave "virus".
    if len(set(names)) < len(names) or any(len(n) <= len(common) for n in names):
        common = ''
    return {format_group(g, common): g for g in groups}


def find_references():
    """Button callback: look the virus up and, when there is only one candidate, pick it."""
    state = st.session_state
    state['ref_error'], state['ref_groups'], state['ref_taxon'] = None, [], None
    try:
        taxon, groups = find_reference_genomes(state['virus_query'].strip(), state['email'].strip())
    except Exception as exc:
        state['ref_error'] = str(exc)
        return
    state['ref_taxon'], state['ref_groups'] = taxon, groups
    state['ref_choice'] = None
    # never leave the previous virus's accession in place while a choice is pending
    state['reference_genome'] = ''
    if not groups:
        state['ref_error'] = ("No RefSeq reference genome was found for this virus. Try a more "
                              "specific name, or type a genome accession below.")
    elif len(groups) == 1:
        state['reference_genome'] = ' '.join(groups[0]['accessions'])


def choose_reference():
    """Selectbox callback: copy the chosen candidate's accessions into the accession field."""
    chosen = candidate_labels(st.session_state.get('ref_groups', [])).get(
        st.session_state.get('ref_choice'))
    if chosen:
        st.session_state['reference_genome'] = ' '.join(chosen['accessions'])


ui.masthead(
    "BV-BRC · Sequence Feature Variant Types",
    "SF Table Builder",
    "Combine UniProt annotations, IEDB epitopes and PubMed-mined features onto a common "
    "reference sequence. Each step runs in order, feeding the next.",
)

# ---------------------------------------------------------------- settings
st.session_state.setdefault('reference_genome', '')
with st.sidebar:
    st.markdown("## Virus")
    email = st.text_input("Email for NCBI", key='email', placeholder="you@example.com",
                          help="NCBI requires an email address with Entrez requests.")
    virus_query = st.text_input("Virus name or taxonomy ID", key='virus_query',
                                placeholder="e.g. measles virus, or 11234")
    st.button("Find reference genomes", on_click=find_references,
              disabled=not (email.strip() and virus_query.strip()),
              help="Looks the virus up at NCBI and lists its RefSeq reference genomes. "
                   "Needs your email address above.")
    taxon, groups = st.session_state.get('ref_taxon'), st.session_state.get('ref_groups', [])
    if st.session_state.get('ref_error'):
        st.caption(st.session_state['ref_error'])
    elif taxon:
        st.caption(f"{taxon['name']} · taxonomy ID {taxon['taxid']}")
    if len(groups) > 1:
        st.selectbox("Reference genome", list(candidate_labels(groups)), index=None,
                     key='ref_choice', on_change=choose_reference,
                     placeholder="Choose one",
                     help="A segmented virus is listed once per strain with all its "
                          "segments. Pick the reference your collaborators use.")
    reference_genome = st.text_input(
        "Reference genome accession(s)", key='reference_genome',
        placeholder="e.g. NC_001498.1",
        help="Filled in by the search above, or type your own. A segmented virus needs one "
             "accession per segment, separated by spaces. Reference proteins are read from "
             "these records' annotated coding sequences.")

    st.markdown("## Run settings")
    default_name = re.sub(r'[^A-Za-z0-9_-]', '', (
        (taxon['name'] if taxon else virus_query).split() or ['SF'])[0]) or 'SF'
    name = st.text_input("Name prefix", placeholder=default_name,
                         help="Names every output file, e.g. <prefix>_SF_Final.csv")
    virus_name = st.text_input("Virus name in sf_id", value="", placeholder="defaults to prefix",
                               help="Short name used in sf_id and sf_name, e.g. Measles_N_SF3.")
    outdir_root = st.text_input("Output folder", value=os.path.expanduser("~/Downloads/SF_runs"),
                                help="Each run gets its own subfolder named after the prefix.")
    row_delay = st.number_input("Mapping row delay (s)", 0.0, 5.0, 0.0, 0.1,
                                help="Raise if NCBI or UniProt rate-limits you.")

# The PubMed search term: what the user typed, or the resolved name when they gave a taxonomy ID.
pubmed_query = (taxon['name'] if taxon and virus_query.strip().isdigit() else virus_query).strip()

run_name = name.strip() or default_name
outdir = os.path.join(os.path.expanduser(outdir_root), run_name)

# ---------------------------------------------------------------- inputs
ui.section(1, "Source files")
left, right = st.columns(2, gap="large")
with left:
    uniprot_upload = st.file_uploader("UniProt taxonomic-search export · .tsv",
                                      type=['tsv', 'txt'])
with right:
    iedb_upload = st.file_uploader("IEDB epitope export · .xlsx", type=['xlsx'])

st.markdown("###### PubMed miner")
mine_pubmed = st.checkbox(
    "Mine PubMed for sequence features", value=True,
    help="Searches PubMed for papers on this virus that report mutations, residues or epitopes, "
         "reads their abstracts, and adds what it finds as a third source. Every row is "
         "labelled 'PubMed miner' and cites its paper by PMID.")
pubmed_upload = st.file_uploader(
    "Extra papers to mine · .csv with a PMID column (optional)", type=['csv'],
    help="A list of papers, such as a pubmed_results.csv, to mine as well as the automatic "
         "search. Only the PMID column is used.")

with st.expander("Optional sources"):
    aliases_upload = st.file_uploader(
        "Extra product aliases · .csv", type=['csv'],
        help="Not needed: product names are matched to the reference proteins automatically. "
             "Upload alias,product pairs only to add to or override what is found.")

upload_problems = [p for p in (
    column_problem(aliases_upload, ['alias', 'product'], "a product alias file",
                   "It needs 'alias' and 'product' columns."),
    column_problem(pubmed_upload, ['PMID'], "a list of PubMed papers",
                   "It needs a 'PMID' column."),
) if p]
for problem in upload_problems:
    ui.failed(problem)

ui.rule()

# ---------------------------------------------------------------- run
ui.section(2, "Run")
missing = []
if uniprot_upload is None:
    missing.append("the UniProt export")
if iedb_upload is None:
    missing.append("the IEDB export")
if not email.strip():
    missing.append("an email address for NCBI")
if not reference_genome.strip():
    missing.append("a reference genome (use Find reference genomes)")
if missing:
    ui.note("Still needed: " + ", ".join(missing) + ".")
elif mine_pubmed and not pubmed_query and pubmed_upload is None:
    ui.note("PubMed mining needs a virus name (in the sidebar) or an uploaded list of papers, "
            "so it will be skipped for this run.")

ui.plan([title for title, _, _ in pipeline_steps(settings_from_form(
    uniprot_tsv='', iedb_xlsx='', reference_genome=reference_genome.strip() or 'REFERENCE',
    outdir=outdir, name=run_name, virus_name=virus_name or None, aliases=aliases_upload,
    pubmed_virus=pubmed_query if (mine_pubmed and pubmed_query) else None,
    pmid_file='papers.csv' if pubmed_upload is not None else None,
    email='', row_delay='0'))])

start = st.button("Run pipeline", type="primary", disabled=bool(missing or upload_problems))

if start:
    os.makedirs(outdir, exist_ok=True)
    upload_dir = tempfile.mkdtemp(prefix="sf_uploads_")

    settings = settings_from_form(
        uniprot_tsv=save_upload(uniprot_upload, upload_dir),
        iedb_xlsx=save_upload(iedb_upload, upload_dir),
        reference_genome=reference_genome.strip(),
        outdir=outdir,
        name=run_name,
        virus_name=virus_name.strip() or None,
        aliases=save_upload(aliases_upload, upload_dir) if aliases_upload else None,
        pubmed_virus=pubmed_query if (mine_pubmed and pubmed_query) else None,
        pmid_file=save_upload(pubmed_upload, upload_dir) if pubmed_upload else None,
        email=email.strip(),
        row_delay=str(row_delay),
    )

    steps = pipeline_steps(settings)
    progress = st.progress(0.0, text="")
    finished = []

    for index, (title, command, output_file) in enumerate(steps, 1):
        progress.progress((index - 1) / len(steps),
                          text=f"Step {index} of {len(steps)} · {title}")
        with st.status(f"{index}. {title}", expanded=False) as status:
            lines = []
            try:
                run_step(index, title, command, log=lines.append)
            except StepFailed as exc:
                status.update(label=f"{index}. {title} — failed", state="error", expanded=True)
                st.code('\n'.join(lines) or str(exc), language='text')
                progress.empty()
                ui.failed(f"{exc}<br>Output from the earlier steps is still in "
                          f"<code>{outdir}</code>.")
                st.stop()
            st.code('\n'.join(lines), language='text')
            status.update(label=f"{index}. {title}", state="complete")
        finished.append((title, output_file))

    progress.empty()
    final_table = finished[-1][1]
    ui.done("Pipeline finished.", final_table)

    # ------------------------------------------------------------ results
    ui.rule()
    ui.section(3, "Results")

    if os.path.exists(final_table):
        table = read_csv(final_table)
        ui.figures([
            (f"{len(table):,}", "features"),
            (f"{table['sf_id'].nunique():,}" if 'sf_id' in table else "—", "sf_ids"),
            (f"{table['gene_ref'].nunique():,}" if 'gene_ref' in table else "—", "genes"),
            (str(len(finished)), "steps"),
        ])
        if 'Source' in table:
            # a feature several sources report is one row listing all of them in Source
            counts = {name: int(table['Source'].astype(str).str.contains(name, regex=False).sum())
                      for name in ('Uniprot', 'IEDB', 'PubMed miner')}
            st.caption("Rows reported by · " + " · ".join(
                f"{k} {v:,}" for k, v in counts.items() if v)
                + "  (a feature reported by several sources is kept once)")
        st.dataframe(table.head(250), use_container_width=True, height=420)

    ui.matching_report(outdir, settings.name)

    ui.rule()
    ui.section(5, "Files")
    produced = [(t, f) for t, f in finished if os.path.exists(f)]
    for row_start in range(0, len(produced), 3):
        for column, (title, output_file) in zip(st.columns(3, gap="small"),
                                                produced[row_start:row_start + 3]):
            with column:
                with open(output_file, 'rb') as handle:
                    st.download_button(os.path.basename(output_file), handle.read(),
                                       file_name=os.path.basename(output_file),
                                       mime='text/csv', key=output_file)
