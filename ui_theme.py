#presentation layer for app.py - the stylesheet and the small set of markup helpers that go
#with it. Kept out of app.py so the app file stays about the pipeline, not about CSS.
import os

import pandas as pd
from sf_columns import read_csv
import streamlit as st

# One restrained palette: warm neutral surfaces, ink text, a single near-black accent. No
# gradients or drop shadows - hairline borders and whitespace carry the hierarchy instead.
STYLES = """
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600&family=JetBrains+Mono:wght@400;500&display=swap');
  :root {
    --ink:        #1C1917;
    --ink-soft:   #44403C;
    --muted:      #78716C;
    --faint:      #A8A29E;
    --line:       #E7E5E4;
    --line-soft:  #F0EFEE;
    --surface:    #FFFFFF;
    --canvas:     #FAFAF9;
    --accent:     #1C1917;
    --ok:         #15803D;
    --ok-bg:      #F0FDF4;
    --warn:       #B45309;
    --err:        #B91C1C;
    --err-bg:     #FEF2F2;
    --radius:     10px;
  }

  html, body, [class*="css"], .stApp, button, input, textarea, select {
    font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
    -webkit-font-smoothing: antialiased;
  }
  .stApp { background: var(--canvas); color: var(--ink); }

  /* Streamlit's own chrome - the hamburger, footer and coloured top bar add nothing here.
     stToolbar is deliberately NOT hidden: the button that reopens a collapsed sidebar lives
     inside it, so display:none there leaves no way to get the sidebar back. The deploy badge
     and menu are already gone via toolbarMode="minimal" in .streamlit/config.toml. */
  #MainMenu, footer, [data-testid="stDecoration"],
  [data-testid="stStatusWidget"] { display: none !important; }
  [data-testid="stHeader"] { background: transparent; height: 3.1rem; }

  /* sidebar show/hide controls - small, quiet, but always reachable */
  [data-testid="stExpandSidebarButton"],
  [data-testid="stSidebarCollapseButton"] button,
  [data-testid="stSidebar"] [data-testid="stBaseButton-headerNoPadding"] {
    /* Streamlit only reveals the collapse control on hover, which makes it easy to collapse
       the sidebar and then not find the way back. Keep both controls always visible. */
    visibility: visible !important; opacity: 1 !important;
    color: var(--muted) !important; border: 1px solid var(--line) !important;
    background: var(--surface) !important; border-radius: 7px !important;
    width: 28px !important; height: 28px !important;
    transition: border-color .15s ease, color .15s ease;
  }
  [data-testid="stExpandSidebarButton"]:hover,
  [data-testid="stSidebar"] [data-testid="stBaseButton-headerNoPadding"]:hover {
    color: var(--ink) !important; border-color: var(--faint) !important;
  }
  [data-testid="stExpandSidebarButton"] span, 
  [data-testid="stSidebar"] [data-testid="stBaseButton-headerNoPadding"] span {
    font-size: 17px !important;
  }

  .block-container { padding: 2.75rem 3rem 5rem; max-width: 1080px; }

  /* ---------- type ---------- */
  h1, h2, h3 { letter-spacing: -0.021em; color: var(--ink); font-weight: 600; }
  p, li, label, .stMarkdown { color: var(--ink-soft); }
  code, pre, .stCode, [data-testid="stCodeBlock"] * {
    font-family: 'JetBrains Mono', ui-monospace, SFMono-Regular, Menlo, monospace !important;
    font-size: 0.78rem !important;
  }

  .masthead { margin-bottom: 2.25rem; }
  .masthead .eyebrow {
    font-size: 0.68rem; font-weight: 600; letter-spacing: 0.1em;
    text-transform: uppercase; color: var(--faint); margin-bottom: 0.6rem;
  }
  .masthead h1 { font-size: 1.95rem; margin: 0 0 0.45rem; line-height: 1.15; }
  .masthead .lede { color: var(--muted); font-size: 0.93rem; max-width: 60ch; margin: 0; }

  .rule { height: 1px; background: var(--line); margin: 2.25rem 0; border: 0; }

  .section-label {
    display: flex; align-items: baseline; gap: 0.7rem;
    font-size: 0.68rem; font-weight: 600; letter-spacing: 0.1em;
    text-transform: uppercase; color: var(--faint); margin: 0 0 1rem;
  }
  .section-label .num {
    font-variant-numeric: tabular-nums; color: var(--ink);
    border: 1px solid var(--line); border-radius: 50%;
    width: 1.35rem; height: 1.35rem; display: inline-flex;
    align-items: center; justify-content: center; font-size: 0.68rem;
  }

  /* ---------- sidebar ---------- */
  [data-testid="stSidebar"] { background: var(--surface); border-right: 1px solid var(--line); }
  [data-testid="stSidebar"] .block-container { padding: 2rem 1.35rem; }
  [data-testid="stSidebar"] h2 {
    font-size: 0.68rem !important; font-weight: 600; letter-spacing: 0.1em;
    text-transform: uppercase; color: var(--faint); margin-bottom: 1.25rem;
  }
  [data-testid="stSidebar"] label { font-size: 0.8rem !important; font-weight: 500; color: var(--ink-soft); }

  /* ---------- inputs ---------- */
  [data-testid="stTextInput"] input, [data-testid="stNumberInput"] input {
    border: 1px solid var(--line); border-radius: 8px; background: var(--surface);
    font-size: 0.85rem; padding: 0.5rem 0.7rem; color: var(--ink);
    transition: border-color .15s ease, box-shadow .15s ease;
  }
  [data-testid="stTextInput"] input:focus, [data-testid="stNumberInput"] input:focus {
    border-color: var(--ink); box-shadow: 0 0 0 3px rgba(28,25,23,.07); outline: none;
  }

  /* ---------- file uploader ---------- */
  [data-testid="stFileUploader"] section {
    border: 1px dashed #DCD9D6; border-radius: var(--radius);
    background: var(--surface); padding: 1.6rem 1.25rem;
    transition: border-color .15s ease, background .15s ease;
  }
  [data-testid="stFileUploaderDropzone"] { gap: 0.9rem; }
  [data-testid="stFileUploader"] section:hover { border-color: var(--faint); background: #FDFDFC; }
  [data-testid="stFileUploader"] label { font-size: 0.82rem !important; font-weight: 500; }
  [data-testid="stFileUploader"] small { color: var(--faint); font-size: 0.72rem; }
  [data-testid="stFileUploaderDropzoneInstructions"] span { font-size: 0.82rem; color: var(--muted); }
  [data-testid="stFileUploader"] section button {
    border: 1px solid var(--line) !important; background: var(--surface) !important;
    color: var(--ink) !important; border-radius: 7px !important;
    font-size: 0.78rem !important; font-weight: 500 !important; padding: 0.35rem 0.8rem !important;
  }

  /* ---------- buttons ---------- */
  .stButton > button {
    background: var(--accent); color: #fff; border: 1px solid var(--accent);
    border-radius: 8px; font-weight: 500; font-size: 0.87rem;
    padding: 0.62rem 1.6rem; letter-spacing: -0.006em; min-width: 170px;
    transition: opacity .15s ease, transform .05s ease;
  }
  .stButton > button:hover:not(:disabled) { opacity: .88; color: #fff; border-color: var(--accent); }
  .stButton > button:active:not(:disabled) { transform: translateY(1px); }
  .stButton > button:disabled {
    background: var(--surface); color: var(--faint); border-color: var(--line); cursor: not-allowed;
  }
  [data-testid="stDownloadButton"] > button {
    background: var(--surface); color: var(--ink); border: 1px solid var(--line);
    border-radius: 8px; font-size: 0.8rem; font-weight: 500; padding: 0.5rem 0.85rem;
  }
  [data-testid="stDownloadButton"] > button:hover { border-color: var(--ink); color: var(--ink); }

  /* ---------- status / expander ---------- */
  [data-testid="stExpander"], [data-testid="stStatus"] {
    border: 1px solid var(--line) !important; border-radius: var(--radius) !important;
    background: var(--surface); box-shadow: none !important;
  }
  [data-testid="stExpander"] summary, [data-testid="stStatus"] summary {
    font-size: 0.85rem; font-weight: 500; color: var(--ink-soft); padding: 0.65rem 0.9rem;
  }
  [data-testid="stCodeBlock"] { background: #FBFBFA !important; border: 1px solid var(--line-soft); border-radius: 8px; }

  /* ---------- progress ----------
     Target the track by its own testid. A '> div > div' selector here lands on the label's
     markdown container instead and squashes it, so the caption and the bar overlap. */
  [data-testid="stProgressBarTrack"] { background: var(--line-soft) !important; height: 3px !important; border-radius: 2px; }
  [data-testid="stProgressBarTrack"] > div { background: var(--ink) !important; border-radius: 2px; }
  [data-testid="stProgress"] p { font-size: 0.78rem; color: var(--muted); margin-bottom: 0.4rem; }
  [data-testid="stProgress"] { margin-bottom: 1rem; }

  /* ---------- callouts ---------- */
  .note {
    border: 1px solid var(--line); border-left: 2px solid var(--faint);
    background: var(--surface); border-radius: 0 8px 8px 0;
    padding: 0.8rem 1rem; font-size: 0.84rem; color: var(--muted);
  }
  .done {
    border: 1px solid #BBF7D0; border-left: 2px solid var(--ok); background: var(--ok-bg);
    border-radius: 0 8px 8px 0; padding: 0.9rem 1.1rem; font-size: 0.86rem; color: #14532D;
  }
  .done .path { font-family: 'JetBrains Mono', monospace; font-size: 0.78rem; color: #166534; }
  .failed {
    border: 1px solid #FECACA; border-left: 2px solid var(--err); background: var(--err-bg);
    border-radius: 0 8px 8px 0; padding: 0.9rem 1.1rem; font-size: 0.86rem; color: #7F1D1D;
  }

  /* ---------- summary figures ---------- */
  .figures { display: flex; gap: 2.5rem; margin: 0.25rem 0 1.25rem; }
  .figure .value { font-size: 1.5rem; font-weight: 600; color: var(--ink); font-variant-numeric: tabular-nums; }
  .figure .key {
    font-size: 0.68rem; letter-spacing: 0.08em; text-transform: uppercase;
    color: var(--faint); margin-top: 0.15rem;
  }

  /* ---------- run plan ---------- */
  ol.plan { list-style: none; margin: 0.25rem 0 0; padding: 0; }
  ol.plan li {
    display: flex; align-items: baseline; gap: 0.85rem;
    padding: 0.5rem 0; border-bottom: 1px solid var(--line-soft);
  }
  ol.plan li:last-child { border-bottom: 0; }
  ol.plan .n {
    font-family: 'JetBrains Mono', monospace; font-size: 0.72rem;
    color: var(--faint); font-variant-numeric: tabular-nums; min-width: 1.1rem;
  }
  ol.plan .t { font-size: 0.86rem; color: var(--ink-soft); }

  /* ---------- dataframe ---------- */
  [data-testid="stDataFrame"] { border: 1px solid var(--line); border-radius: var(--radius); }

  [data-testid="stSpinner"] > div { border-top-color: var(--ink) !important; }
</style>
"""


def apply_styles():
    # st.html, not st.markdown: current Streamlit strips <style> out of markdown, so the
    # stylesheet gets printed on the page as text instead of applied.
    st.html(STYLES)


def masthead(eyebrow, title, lede):
    st.markdown(
        f'<div class="masthead"><div class="eyebrow">{eyebrow}</div>'
        f'<h1>{title}</h1><p class="lede">{lede}</p></div>',
        unsafe_allow_html=True,
    )


def section(number, label):
    st.markdown(
        f'<div class="section-label"><span class="num">{number}</span>{label}</div>',
        unsafe_allow_html=True,
    )


def rule():
    st.markdown('<hr class="rule">', unsafe_allow_html=True)


def note(text):
    st.markdown(f'<div class="note">{text}</div>', unsafe_allow_html=True)


def done(message, path):
    st.markdown(
        f'<div class="done">{message}<br><span class="path">{path}</span></div>',
        unsafe_allow_html=True,
    )


def failed(message):
    st.markdown(f'<div class="failed">{message}</div>', unsafe_allow_html=True)


def plan(steps):
    """A quiet numbered list of the steps a run will perform, shown before it starts."""
    items = ''.join(
        f'<li><span class="n">{index}</span><span class="t">{title}</span></li>'
        for index, title in enumerate(steps, 1)
    )
    st.html(f'<ol class="plan">{items}</ol>')


def figures(pairs):
    """A row of value/label figures, e.g. [('886', 'rows'), ('7', 'steps')]."""
    cells = ''.join(
        f'<div class="figure"><div class="value">{value}</div><div class="key">{key}</div></div>'
        for value, key in pairs
    )
    st.markdown(f'<div class="figures">{cells}</div>', unsafe_allow_html=True)


def matching_report(outdir, name):
    """How well product names were matched to reference proteins, and which names were
    linked automatically or left unmatched. Reads the files suggest_aliases.py and
    Edit_SFTable.py wrote into outdir."""
    merged_path = os.path.join(outdir, f"{name}_SF_merged.csv")
    aliases_path = os.path.join(outdir, f"{name}_aliases.csv")
    unmatched_path = os.path.join(outdir, f"{name}_aliases_unresolved.csv")
    if os.path.exists(merged_path):
        rule()
        section(4, "Product matching")
        merged = read_csv(merged_path)
        matched = int(merged['genbank_accession'].notna().sum())
        figures([
            (f"{matched:,}", "rows matched"),
            (f"{len(merged) - matched:,}", "rows unmatched"),
            (f"{matched / len(merged):.1%}" if len(merged) else "—", "coverage"),
        ])
        st.caption("Each source row is matched to a reference protein by product name. Names "
                   "that differ between NCBI, UniProt and IEDB are linked automatically.")

        if os.path.exists(aliases_path):
            aliases = read_csv(aliases_path)
            with st.expander(f"Names linked automatically · {len(aliases)}"):
                if len(aliases):
                    st.dataframe(aliases, use_container_width=True, hide_index=True)
                    with open(aliases_path, 'rb') as handle:
                        st.download_button("Download aliases", handle.read(),
                                           file_name=os.path.basename(aliases_path),
                                           mime='text/csv', key='aliases_dl')
                else:
                    st.write("Every product name already matched the reference directly.")
        if os.path.exists(unmatched_path):
            unmatched = read_csv(unmatched_path)
            if len(unmatched):
                with st.expander(f"Not matched · {len(unmatched)}"):
                    st.caption("These rows are skipped by mapping. They are usually host proteins, "
                               "or proteins the reference genome does not annotate. If one is a "
                               "protein the reference does have under another name, add the pair "
                               "under Extra product aliases and run again.")
                    st.dataframe(unmatched, use_container_width=True, hide_index=True)

    organism_path = os.path.join(outdir, f"{name}_dropped_other_organism.csv")
    if os.path.exists(organism_path):
        dropped = read_csv(organism_path)
        with st.expander(f"Epitopes from other organisms, dropped · {len(dropped)}"):
            st.caption("IEDB exports also hold epitopes from related viruses and from host "
                       "proteins. Mapping them onto this virus's reference would give "
                       "coordinates that mean nothing there.")
            if len(dropped):
                counts = dropped['virus_type'].value_counts().rename_axis('organism')
                st.dataframe(counts.reset_index(name='rows'), use_container_width=True,
                             hide_index=True)

    mapped_path = os.path.join(outdir, f"{name}_SF_mapped.csv")
    if os.path.exists(mapped_path):
        mapped = read_csv(mapped_path)
        placeable = mapped[mapped['genbank_accession'].notna()]
        failed = int(placeable['segments'].isna().sum())
        if failed:
            with st.expander(f"Matched but not mappable · {failed}"):
                st.caption("These rows have a reference protein but could not be placed on it: "
                           "the source sequence could not be fetched, was too different from "
                           "the reference to be the same protein, or the feature lies in a "
                           "stretch the reference lacks. The run log has each reason.")
                st.dataframe(placeable[placeable['segments'].isna()][
                    [c for c in ('source_id', 'sf_type', 'source_sf_location', 'product')
                     if c in placeable.columns]].head(200),
                    use_container_width=True, hide_index=True)
