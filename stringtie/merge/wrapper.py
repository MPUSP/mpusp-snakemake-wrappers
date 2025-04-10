__author__ = "Michael Jahn"
__copyright__ = "Copyright 2025, Michael Jahn"
__email__ = "jahn@mpusp.mpg.de"
__license__ = "MIT"


import os
from snakemake.shell import shell

# define default options for this wrapper
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
output = snakemake.output.get("gff")

# get input files
list_gff = snakemake.input.get("gff")
ref_tx = snakemake.input.get("reference", "")


def concatenate_input(files):
    return " ".join(files) if isinstance(files, list) else files


# prepare input
gff = concatenate_input(list_gff)

# supply genome or transcription start sites
if ref_tx:
    ref_tx = "-G " + ref_tx

# run StringTie
shell(
    "stringtie --merge "  # merge mode
    "{extra} "  # additional options
    "{ref_tx} "  # optional reference transcriptome
    "-o {output} "  # output file
    "{gff} {log};"  # list of GFF files
)
