def test_filter_rows_judges_iedb_and_pubmed_rows():
    """A paper about Ebola also discusses Marburg and Reston, and may name only the genus."""
    df = pd.DataFrame([
        {"Source": "IEDB", "virus_type": "Zaire ebolavirus"},
        {"Source": "IEDB", "virus_type": "Marburg virus"},
        {"Source": "PMID", "virus_type": "Ebola virus"},
        {"Source": "PMID", "virus_type": "Reston virus"},
        {"Source": "PMID", "virus_type": "Ebolavirus"},       # the genus this virus sits in
        {"Source": "PMID", "virus_type": "Mystery agent"},
        {"Source": "Uniprot", "virus_type": "Marburg virus"},  # never judged
    ])
    lineages = {"Zaire ebolavirus": ["186538", "1570291", "10292"],
                "Ebola virus": ["1570291", "186538", "10292"],
                "Marburg virus": ["11269", "10293"],
                "Reston virus": ["186539", "10292"],
                "Ebolavirus": ["10292"]}
    taxids = {"Zaire ebolavirus": "186538", "Ebola virus": "1570291",
              "Marburg virus": "11269", "Reston virus": "186539", "Ebolavirus": "10292"}
    lookup = lambda n: {"lineage": lineages[n], "taxid": taxids[n]} if n in lineages else None

    kept, dropped, unknown = filter_rows(df, "186538", lookup,
                                         target_lineage=["186538", "10292"],
                                         target_words=name_words("Zaire ebolavirus"))

    assert list(kept["virus_type"]) == ["Zaire ebolavirus", "Ebola virus", "Ebolavirus",
                                        "Marburg virus"], list(kept["virus_type"])
    assert list(dropped["virus_type"]) == ["Marburg virus", "Reston virus", "Mystery agent"]
    assert unknown == ["Mystery agent"]

    # NCBI has no searchable 'Marburg virus' (now Orthomarburgvirus marburgense), so an
    # unresolvable name is judged on whether it shares a word with the virus being built
    target = name_words("Zaire ebolavirus")
    assert names_agree("Ebola virus", target) and names_agree("Ebolavirus", target)
    assert not names_agree("Marburg virus", target)


#self-checks for the reference-building and reference-matching logic - no network needed
import os
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation

from build_reference import cds_rows, add_aliases
from Edit_SFTable import check_and_update_genbank_accession
from mapping import map_feature, percent_identity, is_pair_feature
from suggest_aliases import suggest, gene_tokens
from filter_organism import filter_rows, name_words, names_agree, organism_agreement
from SF_Table_Pubmed import clean_position, convert, feature_length, sf_category, split_features
from sf_columns import PUBMED_SOURCE, publication_label, read_csv, write_csv
from build_sf_ids_names import add_pubmed_flag, name_virus_column
from merge_same_SF_coordinates import merge_duplicate_sf_rows
from find_reference import group_records


def _cds(**qualifiers):
    return SeqFeature(SimpleLocation(0, 9), type="CDS", qualifiers=qualifiers)


def test_cds_rows():
    record = SeqRecord(Seq("ATGAAATAG"), id="TEST.1", features=[
        _cds(gene=["N"], product=["nucleocapsid protein"], protein_id=["NP_000001.1"]),
        # no /product -> falls back to the gene name
        _cds(gene=["M"], protein_id=["NP_000002.1"]),
        # no /protein_id -> skipped, there is nothing to map onto
        _cds(gene=["X"], product=["hypothetical protein"]),
        # not a CDS -> skipped
        SeqFeature(SimpleLocation(0, 9), type="mat_peptide",
                   qualifiers={"product": ["fusion glycoprotein F1"], "protein_id": ["NP_000003.1"]}),
    ])

    rows = cds_rows(record)

    assert rows == [
        {"product": "nucleocapsid protein", "gene": "N", "protein_accession": "NP_000001.1"},
        {"product": "M", "gene": "M", "protein_accession": "NP_000002.1"},
    ], rows


def test_add_aliases_reuses_target_gene_and_accession():
    rows = [{"product": "C protein", "gene": "P/V/C", "protein_accession": "NP_000003.1"}]
    aliases = pd.DataFrame([{"alias": "protein c", "product": "C protein"}])
    aliases.to_csv("_aliases_tmp.csv", index=False)

    out = add_aliases(list(rows), "_aliases_tmp.csv")

    assert out[1] == {"product": "protein c", "gene": "P/V/C",
                      "protein_accession": "NP_000003.1"}, out
    os.remove("_aliases_tmp.csv")


def test_match_respects_word_boundaries():
    """'myelin basic protein' must NOT match the reference product 'C protein', and 'PB1' must
    not be found inside the different protein 'PB1-F2'."""
    pd.DataFrame([
        {"product": "C protein", "gene": "P/V/C", "protein_accession": "NP_000003.1"},
        {"product": "matrix protein", "gene": "M", "protein_accession": "NP_000004.1"},
        {"product": "PB1", "gene": "PB1", "protein_accession": "NP_000005.1"},
    ]).to_csv("_ref_tmp.csv", index=False)
    merged = pd.DataFrame([
        {"product": "Myelin basic protein"},
        {"product": "C protein [Measles morbillivirus]"},
        {"product": "Protein PB1-F2"},
    ])

    out = check_and_update_genbank_accession(merged, "_ref_tmp.csv")

    assert pd.isna(out.at[0, "genbank_accession"]), out.at[0, "genbank_accession"]
    assert out.at[1, "genbank_accession"] == "NP_000003.1", out.at[1, "genbank_accession"]
    assert pd.isna(out.at[2, "genbank_accession"]), out.at[2, "genbank_accession"]
    os.remove("_ref_tmp.csv")


def test_single_protein_reference_takes_every_row():
    """A polyprotein virus has one reference protein, and every protein is part of it."""
    pd.DataFrame([{"product": "polyprotein", "gene": "POLY", "protein_accession": "YP_1.1"}]
                 ).to_csv("_ref_tmp.csv", index=False)
    merged = pd.DataFrame([{"product": "Envelope protein E"}, {"product": "polyprotein"}])

    out = check_and_update_genbank_accession(merged, "_ref_tmp.csv")

    assert list(out["genbank_accession"]) == ["YP_1.1", "YP_1.1"], list(out["genbank_accession"])
    assert list(out["gene"]) == ["POLY", "POLY"]
    os.remove("_ref_tmp.csv")


def test_gene_symbol_NA_is_not_read_as_missing():
    """'NA' is influenza neuraminidase. pandas reads the text 'NA' as missing by default, which
    silently stripped the gene (and the reference match) from every neuraminidase row."""
    pd.DataFrame([
        {"product": "neuraminidase", "gene": "NA", "protein_accession": "NP_1"},
        {"product": "haemagglutinin", "gene": "HA", "protein_accession": "NP_2"},
    ]).to_csv("_ref_tmp.csv", index=False)
    merged = pd.DataFrame([{"product": "Neuraminidase (EC 3.2.1.18)"}])

    out = check_and_update_genbank_accession(merged, "_ref_tmp.csv")
    assert out.at[0, "gene"] == "NA" and out.at[0, "gene_ref"] == "NA", out.loc[0].to_dict()

    # a real blank is still missing, and the symbol survives a write/read round trip
    write_csv(pd.DataFrame([{"gene": "NA", "start": 1}, {"gene": "", "start": 2}]), "_rt_tmp.csv")
    back = read_csv("_rt_tmp.csv")
    assert back.at[0, "gene"] == "NA" and pd.isna(back.at[1, "gene"]), back.to_dict("records")
    os.remove("_ref_tmp.csv"); os.remove("_rt_tmp.csv")


def test_pubmed_rows_are_cited_by_pmid_when_merged_into_another_source():
    """A feature UniProt and the PubMed miner both report is kept once; the miner's paper must
    still be cited, as a clear PMID."""
    assert PUBMED_SOURCE == "PMID"
    assert publication_label(PUBMED_SOURCE, "40493685") == "PMID: 40493685"
    assert publication_label("PubMed", "40493685") == "PMID: 40493685"   # older tables
    assert publication_label("IEDB", "2322") == "IEDB:2322"

    df = pd.DataFrame([
        {"gene": "H", "segments": "187-187", "sf_category": "Sequence_Alterations",
         "Source": "Uniprot", "source_id": "P08362", "publication": "PMID: 111"},
        {"gene": "H", "segments": "187-187", "sf_category": "Sequence_Alterations",
         "Source": PUBMED_SOURCE, "source_id": "40493685", "publication": "PMID: 40493685"},
        {"gene": "H", "segments": "300-300", "sf_category": "Sequence_Alterations",
         "Source": PUBMED_SOURCE, "source_id": "39641619", "publication": "PMID: 39641619"},
    ])

    merged, count = merge_duplicate_sf_rows(df)

    assert count == 1 and len(merged) == 2, (count, len(merged))
    survivor = merged[merged["segments"] == "187-187"].iloc[0]
    assert survivor["publication"] == "PMID: 111, PMID: 40493685", survivor["publication"]
    # a row names exactly one source, and source_id is that source's own identifier
    assert survivor["Source"] == "Uniprot" and survivor["source_id"] == "P08362", survivor.to_dict()
    # a feature only the miner found keeps its own label and PMID
    only_pubmed = merged[merged["segments"] == "300-300"].iloc[0]
    assert (only_pubmed["Source"], only_pubmed["publication"]) == (PUBMED_SOURCE, "PMID: 39641619")


def test_a_paper_cited_by_two_sources_is_listed_once():
    df = pd.DataFrame([
        {"gene": "P", "segments": "172-172", "sf_category": "Sequence_Alterations",
         "Source": "Uniprot", "source_id": "P12579", "publication": "PMID: 111, 222"},
        {"gene": "P", "segments": "172-172", "sf_category": "Sequence_Alterations",
         "Source": PUBMED_SOURCE, "source_id": "222", "publication": "PMID: 222"},
        {"gene": "P", "segments": "172-172", "sf_category": "Sequence_Alterations",
         "Source": PUBMED_SOURCE, "source_id": "333", "publication": "PMID: 333"},
    ])
    merged, _ = merge_duplicate_sf_rows(df)
    assert merged.iloc[0]["publication"] == "PMID: 111, 222, PMID: 333", merged.iloc[0]["publication"]


def test_coordinate_ranges_survive_excel_and_round_trip():
    """'1-36' is what Excel turns into a date. Quoting alone does not stop it, so those cells
    get an ="1-36" guard; every field is quoted; and reading strips the guard again so the
    pipeline and other tools see the plain value."""
    frame = pd.DataFrame([{"segments": "1-36", "source_sf_location": "188-606", "start": 1,
                           "length": 36, "gene": "NA", "comments": "Chain, Nucleoprotein"}])
    write_csv(frame, "_x_tmp.csv")
    text = open("_x_tmp.csv").read()

    assert '"=""1-36"""' in text and '"=""188-606"""' in text, text
    assert text.splitlines()[0] == '"segments","source_sf_location","start","length","gene","comments"'

    back = read_csv("_x_tmp.csv")
    assert back.at[0, "segments"] == "1-36" and back.at[0, "source_sf_location"] == "188-606"
    assert back.at[0, "start"] == 1 and str(back["start"].dtype) == "int64"
    assert back.at[0, "gene"] == "NA"                     # still not read as missing
    assert back.at[0, "comments"] == "Chain, Nucleoprotein"   # a comma inside a field survives

    write_csv(frame, "_x_tmp.csv", excel_safe=False)
    assert '"1-36"' in open("_x_tmp.csv").read()
    os.remove("_x_tmp.csv")


def test_a_reference_from_another_virus_is_caught():
    """A measles reference used for a Zika table mapped 1 row of 1005 with no error."""
    lineages = {"Zika virus (strain Mr 766)": ["64320"], "Zika virus": ["64320"],
                "Measles virus (strain Edmonston)": ["11234"]}
    lookup = lambda name: {"lineage": lineages[name]} if name in lineages else None
    zika_table = {"Zika virus (strain Mr 766)": 990, "Zika virus": 10, "Unknown thing": 5}

    assert organism_agreement(zika_table, "64320", lookup) == 1.0       # right reference
    assert organism_agreement(zika_table, "11234", lookup) == 0.0       # measles reference
    assert organism_agreement({"Unknown thing": 3}, "64320", lookup) is None


def test_miner_findings_become_sf_rows_labelled_PMID():
    """A findings export row -> an SF row: Source is PMID, source_id is the numeric PMID, the
    quote becomes the comment, and the yes/no flag sits right after source_id."""
    findings = pd.DataFrame([
        {"pmid": 38783022, "virus": "Ebola virus", "protein": "VP35",
         "sequence_feature_name": "alpha-helix", "start_position": 26, "end_position": 36,
         "position": "26-36", "quote": "NPBP adopts two orthogonal helices.",
         "target_type": "domain", "length": 11},
        # discrete residues: 'position' is authoritative, start/end are only its first and last
        {"pmid": 38783022, "virus": "Ebola virus", "protein": "VP35",
         "sequence_feature_name": "mutation_site", "start_position": 222, "end_position": 222,
         "position": "222, 225, 248, 251", "quote": "Mutating these residues.",
         "target_type": "mutation", "length": 4},
    ])

    table = convert(findings)

    assert list(table["Source"]) == [PUBMED_SOURCE, PUBMED_SOURCE] == ["PMID", "PMID"]
    assert set(table["Evidence_code"]) == {"EXP"}
    assert list(table["source_id"]) == ["38783022", "38783022"]
    assert list(table["source_sf_location"]) == ["26-36", "222,225,248,251"]
    assert list(table["comments"]) == ["NPBP adopts two orthogonal helices.",
                                       "Mutating these residues."]
    assert list(table["sf_category"]) == ["Structural", "Sequence_Alterations"]
    assert list(table["publication"]) == ["PMID: 38783022"] * 2
    # the quote is kept as evidence in its own right, beside the comment
    assert list(table["evidence_sentence"]) == list(table["comments"])
    assert table["review_note"].str.startswith("coordinates assumed").all()

    # length is counted from the position when the export does not give one
    assert feature_length("222, 225, 248, 251", None) == 4
    assert feature_length("26-36", None) == 11
    assert clean_position(" 251,222 , 225 ") == "222,225,251"
    # a list that mixes a span with single residues keeps every residue the span covers
    assert clean_position("143, 29-31") == "29,30,31,143"
    # names and positions as parallel lists are several features, not one
    assert split_features("Long loop, L209, E211", "146-340, 209, 211") == [
        ("Long loop", "146-340"), ("L209", "209"), ("E211", "211")]
    # one name against several positions is a single feature spanning them
    assert split_features("interface", "21, 22, 24") == [("interface", "21, 22, 24")]
    assert sf_category("", "neutralizing epitope") == "Epitope"


def test_virus_type_is_one_clean_name_and_the_source_name_is_kept():
    df = pd.DataFrame([{"virus_type": "Ebola virus (strain Mayinga-76)", "gene": "VP35"},
                       {"virus_type": "Zaire ebolavirus", "gene": "NP"}])

    out = name_virus_column(df, "Ebola")

    assert list(out.columns) == ["virus_type", "virus_type_source", "gene"]
    assert list(out["virus_type"]) == ["Ebola", "Ebola"]
    assert list(out["virus_type_source"]) == ["Ebola virus (strain Mayinga-76)",
                                              "Zaire ebolavirus"]
    # running it twice must not stack the rename
    assert list(name_virus_column(out, "Ebola").columns) == list(out.columns)


def test_pubmed_flag_column_sits_after_source_id():
    df = pd.DataFrame([{"Source": "PMID", "source_id": "38783022", "gene": "VP35"},
                       {"Source": "Uniprot", "source_id": "Q05318", "gene": "VP35"}])

    out = add_pubmed_flag(df)

    assert list(out.columns) == ["Source", "source_id", "sourced from pubmed miner", "gene"]
    assert list(out["sourced from pubmed miner"]) == ["yes", "no"]
    # running it twice must not add the column twice
    assert list(add_pubmed_flag(out).columns) == list(out.columns)


def test_paper_text_keeps_its_characters_through_the_pipeline():
    """A quote from a paper holds en-dashes and Greek letters. Reading a table as latin-1 (as
    the mapping and id steps used to) turned every one into mojibake."""
    sentence = "NPBP adopts one \u03b1-helix, residues 20\u201348, binding NP."
    write_csv(pd.DataFrame([{"comments": sentence, "segments": "20-48"}]), "_e_tmp.csv")

    assert read_csv("_e_tmp.csv").at[0, "comments"] == sentence

    # a stray Windows-1252 byte among proper UTF-8 must not stop the read
    with open("_e_tmp.csv", "rb") as handle:
        patched = handle.read().replace(b"20\xe2\x80\x9348", b"20\x9348")
    with open("_e_tmp.csv", "wb") as handle:
        handle.write(patched)
    recovered = read_csv("_e_tmp.csv").at[0, "comments"]
    assert "\u03b1-helix" in recovered and "binding NP." in recovered, recovered
    os.remove("_e_tmp.csv")


def _alignment(ref, source):
    from Bio.Align import MultipleSeqAlignment
    return MultipleSeqAlignment([SeqRecord(Seq(ref), id="ref"), SeqRecord(Seq(source), id="src")])


def test_pair_features_map_two_residues_not_the_span():
    """'a-b' on a disulfide bond means residue a and residue b, not everything between."""
    alignment = _alignment("ACDEFGHIK", "ACDEFGHIK")
    assert is_pair_feature("Disulfide bond") and is_pair_feature(" Cross-link ")
    assert not is_pair_feature("Helix")

    span = map_feature(alignment, "2-7", is_pair=False)
    pair = map_feature(alignment, "2-7", is_pair=True)

    assert span[2:] == ("2-7", "CDEFGH", "CDEFGH"), span
    assert pair[2:] == ("2,7", "CH", "CH"), pair


def test_a_single_residue_pair_is_one_residue_not_two_copies():
    """UniProt can report a disulfide as '139-139'. That is one residue, not the same residue
    twice ('139,139' with sequence 'CC')."""
    alignment = _alignment("ACDEFGHIK", "ACDEFGHIK")

    pair = map_feature(alignment, "2-2", is_pair=True)
    listed = map_feature(alignment, "2,2,5")

    assert pair == (2, 2, "2-2", "C", "C"), pair
    assert listed[2:] == ("2,5", "CF", "CF"), listed


def test_map_feature_reports_reference_positions_not_alignment_columns():
    """A source with 3 extra N-terminal residues puts every alignment column 3 past its
    reference position. Reading columns as positions misplaces every feature after the gap."""
    alignment = _alignment("---ACDEFGHIK", "MKVACDEFGHIK")

    start, end, segments, source_seq, ref_seq = map_feature(alignment, "7-9")

    # source residues 7-9 are E,F,G = reference residues 4-6 (columns 7-9)
    assert (start, end, segments) == (4, 6, "4-6"), (start, end, segments)
    assert (source_seq, ref_seq) == ("EFG", "EFG")

    # a range that starts inside the insertion keeps only the residues the reference has
    start, end, _, source_seq, _ = map_feature(alignment, "2-5")
    assert (start, end, source_seq) == (1, 2, "AC"), (start, end, source_seq)

    # a feature wholly inside the insertion has no reference position at all
    try:
        map_feature(alignment, "1-3")
        assert False, "expected a feature inside an insertion to be rejected"
    except ValueError:
        pass


def test_percent_identity_tells_a_different_protein_from_a_strain():
    assert percent_identity(_alignment("ACDEFGHIKL", "ACDEFGHIKL")) == 1.0
    assert percent_identity(_alignment("ACDEFGHIKL", "ACDEFGHIKV")) == 0.9
    assert percent_identity(_alignment("ACDEFGHIKL", "WWWWWWWWWW")) == 0.0
    # gaps are not mismatches: only residues present in both are compared
    assert percent_identity(_alignment("ACDEF---KL", "ACDEFGHIKL")) == 1.0


def test_suggest_aliases_links_synonyms_and_leaves_host_proteins():
    reference = pd.DataFrame([
        {"product": "nucleocapsid protein", "gene": "N", "protein_accession": "NP_1"},
        {"product": "C protein", "gene": "P/V/C", "protein_accession": "NP_2"},
        {"product": "large polymerase protein", "gene": "L", "protein_accession": "NP_3"},
    ])
    rows = {
        # matches directly, and its bracketed synonym links the IEDB spelling below
        "Nucleoprotein (Nucleocapsid protein) (NP) (Protein N)": 5,
        "nucleoprotein [Measles morbillivirus]": 2,
        "RNA-directed RNA polymerase L (Protein L) [Includes: GTP phosphohydrolase]": 4,
        "Protein C": 1,
        # contains the letters 'c protein' but is a different protein entirely
        "Myelin basic protein": 1,
    }

    aliases, unresolved = suggest(reference, rows)
    linked = {a["alias"].lower(): a["product"] for a in aliases}

    assert linked["nucleoprotein"] == "nucleocapsid protein", linked
    assert linked["protein c"] == "C protein", linked
    assert any(v == "large polymerase protein" for v in linked.values()), linked
    assert unresolved == ["Myelin basic protein"], unresolved


def test_suggest_aliases_uses_the_primary_gene_symbol():
    reference = pd.DataFrame([
        {"product": "haemagglutinin", "gene": "HA", "protein_accession": "NP_1"},
        {"product": "polymerase PA", "gene": "PA", "protein_accession": "NP_2"},
        {"product": "PA-X protein", "gene": "PA-X", "protein_accession": "NP_3"},
    ])
    rows = {"Hemagglutinin": 4, "Polymerase acidic protein": 3}
    # spelling differs ('haemagglutinin'), and UniProt lists PA's gene as 'PA PA-X' because the
    # gene also makes PA-X - only the first symbol is the gene itself
    genes = {"Hemagglutinin": {"HA"}, "Polymerase acidic protein": gene_tokens("PA PA-X")}

    aliases, unresolved = suggest(reference, rows, genes)
    linked = {a["alias"]: a["product"] for a in aliases}

    assert gene_tokens("PA PA-X") == {"PA"} and gene_tokens("P/V/C") == {"P", "V", "C"}
    assert linked == {"Hemagglutinin": "haemagglutinin",
                      "Polymerase acidic protein": "polymerase PA"}, linked
    assert unresolved == []


def test_filter_rows_drops_other_organisms_only_from_iedb():
    df = pd.DataFrame([
        {"Source": "IEDB", "virus_type": "Zika virus"},
        {"Source": "IEDB", "virus_type": "Dengue virus"},
        {"Source": "IEDB", "virus_type": "Homo sapiens"},
        {"Source": "IEDB", "virus_type": "Mystery agent"},
        {"Source": "Uniprot", "virus_type": "Dengue virus (strain X)"},
    ])
    lineages = {"Zika virus": ["64320", "10239"], "Dengue virus": ["12637", "10239"],
                "Homo sapiens": ["9606"]}
    lookup = lambda name: {"lineage": lineages[name]} if name in lineages else None

    kept, dropped, unknown = filter_rows(df, "64320", lookup)

    assert list(kept["virus_type"]) == ["Zika virus", "Mystery agent", "Dengue virus (strain X)"]
    assert list(dropped["virus_type"]) == ["Dengue virus", "Homo sapiens"]
    assert unknown == ["Mystery agent"]   # unidentifiable is kept, not silently discarded


def test_find_reference_separates_segments_from_alternative_genomes():
    rec = lambda a, t, n: {"AccessionVersion": a, "TaxId": t, "Length": n}
    zika = group_records([rec("NC_1.1", 64320, 10807), rec("NC_2.1", 64320, 10795),
                          rec("NG_3.1", 64320, 500)], {"64320": "Zika virus"})
    flu = group_records([rec(f"NC_0{i}.1", 211044, n) for i, n in
                         enumerate([2341, 2341, 2233, 1778, 1565, 1413, 1027, 890])],
                        {"211044": "Influenza A virus"})

    assert [g["accessions"] for g in zika] == [["NC_1.1"], ["NC_2.1"]], zika
    assert len(flu) == 1 and flu[0]["records"] == 8, flu


if __name__ == "__main__":
    test_cds_rows()
    test_add_aliases_reuses_target_gene_and_accession()
    test_match_respects_word_boundaries()
    test_single_protein_reference_takes_every_row()
    test_gene_symbol_NA_is_not_read_as_missing()
    test_coordinate_ranges_survive_excel_and_round_trip()
    test_paper_text_keeps_its_characters_through_the_pipeline()
    test_a_reference_from_another_virus_is_caught()
    test_pubmed_rows_are_cited_by_pmid_when_merged_into_another_source()
    test_a_paper_cited_by_two_sources_is_listed_once()
    test_miner_findings_become_sf_rows_labelled_PMID()
    test_pubmed_flag_column_sits_after_source_id()
    test_virus_type_is_one_clean_name_and_the_source_name_is_kept()
    test_pair_features_map_two_residues_not_the_span()
    test_a_single_residue_pair_is_one_residue_not_two_copies()
    test_map_feature_reports_reference_positions_not_alignment_columns()
    test_percent_identity_tells_a_different_protein_from_a_strain()
    test_suggest_aliases_links_synonyms_and_leaves_host_proteins()
    test_suggest_aliases_uses_the_primary_gene_symbol()
    test_filter_rows_judges_iedb_and_pubmed_rows()
    test_find_reference_separates_segments_from_alternative_genomes()
    print("ok")
