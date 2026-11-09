#run the whole SF pipeline end to end: every step in order, each step's output wired into the
#next, so a run needs only the input files and a virus name.
import argparse
import os
import subprocess
import sys

PYTHON = sys.executable
HERE = os.path.dirname(os.path.abspath(__file__))


class StepFailed(Exception):
    pass


def run_step(number, title, command, log=print):
    """Run one pipeline step. Raises StepFailed with the captured output if it fails."""
    log(f"\n=== Step {number}: {title}")
    log('    ' + ' '.join(str(c) for c in command))
    process = subprocess.run([str(c) for c in command], cwd=HERE,
                             capture_output=True, text=True)
    output = (process.stdout or '') + (process.stderr or '')
    for line in output.rstrip().splitlines():
        log('    ' + line)
    if process.returncode != 0:
        raise StepFailed(f"Step {number} ({title}) failed with exit code {process.returncode}")
    return output


def pipeline_steps(args):
    """The ordered steps, as (title, command, output_file) triples.

    Every output path is derived from outdir+name, so nothing has to be wired up by hand.
    """
    out, name = args.outdir, args.name
    path = lambda filename: os.path.join(out, filename)
    # a segmented virus has one accession per segment: space- or comma-separated
    accessions = args.reference_genome.replace(',', ' ').split()
    if not accessions:
        raise ValueError("No reference genome accession was given")

    steps = [
        ("UniProt SF table",
         [PYTHON, 'SF_Table_Uniprot.py', args.uniprot_tsv, '-o', out, '-n', name],
         path(f'{name}_Uniprot_SF.csv')),
        ("IEDB epitope SF table",
         [PYTHON, 'SF_Table_Epitope.py', args.iedb_xlsx, '-o', out, '-n', name],
         path(f'{name}_Epitope_SF.csv')),
    ]

    source_tables = [path(f'{name}_Uniprot_SF.csv'), path(f'{name}_Epitope_SF.csv')]
    findings = getattr(args, 'pubmed_findings', None)
    if findings:
        # the third source: a PubMed-miner findings export, converted to the SF table's shape
        steps.append(
            ("PubMed miner SF table",
             [PYTHON, 'SF_Table_Pubmed.py', findings, '-o', out, '-n', name],
             path(f'{name}_Pubmed_SF.csv')))
        source_tables.append(path(f'{name}_Pubmed_SF.csv'))

    steps.append(
        # after the source tables: it checks the genome is the same virus as the UniProt table
        ("Reference table from the reference genome",
         [PYTHON, 'build_reference.py'] + accessions + ['-o', out, '-n', name]
         + (['--aliases', args.aliases] if args.aliases else [])
         + ['--check-against', path(f'{name}_Uniprot_SF.csv'), '--email', args.email],
         path(f'{name}_reference.csv')))

    steps += [
        ("Match product names to reference proteins",
         [PYTHON, 'suggest_aliases.py', path(f'{name}_reference.csv')] + source_tables
         + ['-o', out, '-n', name],
         path(f'{name}_reference_aliased.csv')),
    ]

    merge_inputs = list(source_tables)

    steps += [
        ("Merge source tables and backfill accessions",
         [PYTHON, 'Edit_SFTable.py'] + merge_inputs
         + ['-o', out, '-n', name, '--reference-file', path(f'{name}_reference_aliased.csv')],
         path(f'{name}_SF_merged.csv')),
        ("Drop epitopes from other organisms",
         [PYTHON, 'filter_organism.py', path(f'{name}_SF_merged.csv'),
          '--reference', accessions[0], '-o', out, '-n', name, '--email', args.email],
         path(f'{name}_SF_filtered.csv')),
        ("Map features onto the reference",
         [PYTHON, 'mapping.py', path(f'{name}_SF_filtered.csv'), '-o', out, '-n', name,
          '--email', args.email, '--row-delay', args.row_delay],
         path(f'{name}_SF_mapped.csv')),
        ("Collapse duplicate features",
         [PYTHON, 'merge_same_SF_coordinates.py', path(f'{name}_SF_mapped.csv'),
          '-o', out, '-n', name, '--drop-unmapped'],
         path(f'{name}_dedup.csv')),
        ("Assign sf_id and sf_name",
         [PYTHON, 'build_sf_ids_names.py', path(f'{name}_dedup.csv'), '-o', out, '-n', name,
          '--virus-name', args.virus_name or name],
         path(f'{name}_SF_Final.csv')),
    ]
    return steps


def run_all(args, log=print, on_step=None):
    """Run every step in order. on_step(index, title, output_file) fires after each one."""
    os.makedirs(args.outdir, exist_ok=True)
    steps = pipeline_steps(args)
    for index, (title, command, output_file) in enumerate(steps, 1):
        run_step(index, title, command, log=log)
        if on_step:
            on_step(index, title, output_file)
    final = steps[-1][2]
    log(f"\nPipeline finished. Final SF table: {final}")
    return final


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run the whole SF pipeline end to end, wiring each step's output into "
                    "the next."
    )
    parser.add_argument("uniprot_tsv", help="UniProt taxonomic-search TSV export.")
    parser.add_argument("iedb_xlsx", help="IEDB epitope search export (.xlsx).")
    parser.add_argument("reference_genome",
                        help="Reference genome accession, e.g. 'NC_001498.1' for measles. A "
                             "segmented virus has one per segment: pass them comma-separated "
                             "(find_reference.py lists them for a virus name).")
    parser.add_argument("-o", "--outdir", required=True,
                        help="Directory for all outputs (created if missing).")
    parser.add_argument("-n", "--name", default="output", help="Prefix for every output file.")
    parser.add_argument("--virus-name", default=None,
                        help="Short virus name for sf_id/sf_name (defaults to --name).")
    parser.add_argument("--aliases", default=None,
                        help="Optional alias CSV (alias, product) to add to the automatically "
                             "suggested ones - only needed to override or extend them.")
    parser.add_argument("--pubmed-findings", default=None,
                        help="A PubMed-miner findings CSV to merge as a third source. Its rows "
                             "are labelled Source=PMID with the paper's PMID as source_id.")
    parser.add_argument("--email", default="rbhattac@jcvi.org",
                        help="Email address sent with NCBI Entrez requests.")
    parser.add_argument("--row-delay", default="0",
                        help="Seconds to pause between rows in the mapping step (default: 0).")
    return parser.parse_args()


def main():
    try:
        run_all(parse_args())
    except StepFailed as exc:
        print(f"\n{exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
