"""
Centralized configuration.
"""

from os.path import dirname, realpath

root_dir = dirname(realpath(__file__)) + "/../"

imgthla_git = "https://github.com/ANHIG/IMGTHLA.git"
imgthla_dir = root_dir + "dat/IMGTHLA/"
hla_dat = imgthla_dir + "hla.dat"
hla_nom_g = imgthla_dir + "wmda/hla_nom_g.txt"
hla_nom_p = imgthla_dir + "wmda/hla_nom_p.txt"

ref_dir = root_dir + "dat/ref/"
allele_group_json = ref_dir + "allele_groups.json"
cdna_json = ref_dir + "cDNA.json"
cdna_single_json = ref_dir + "cDNA.single.json"
hla_transcripts_json = ref_dir + "hla_transcripts.json"
hla_convert_json = ref_dir + "hla.convert.json"
hla_fa = ref_dir + "hla.fasta"
hla_idx = ref_dir + "hla.idx"
hla_json = ref_dir + "hla.p.json"
partial_fa = ref_dir + "hla_partial.fasta"
partial_idx = ref_dir + "hla_partial.idx"
partial_json = ref_dir + "hla_partial.p.json"
ref_zip_archive = ref_dir + "customization_reference_fastas.tar.gz"
zipped_ref_files = {
    "GRCh38_chr6": ref_dir + "GRCh38.chr6.noHLA.fasta",
    "GRCh38": ref_dir + "GRCh38.all.noHLA.fasta",
    "dummy_HLA_fa": ref_dir + "GRCh38.chr6.HLA.fasta",
}

info_dir = root_dir + "dat/info/"
alt_decoys = info_dir + "decoys_alts.json"
hla_freq = info_dir + "hla_freq.tsv"
parameters_json = info_dir + "parameters.json"
