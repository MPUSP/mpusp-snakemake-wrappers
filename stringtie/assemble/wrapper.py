__author__ = "Michael Jahn"
__copyright__ = "Copyright 2025, Michael Jahn"
__email__ = "jahn@mpusp.mpg.de"
__license__ = "MIT"


import os
from snakemake.shell import shell

# define default options for this wrapper
threads = snakemake.threads
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
output = snakemake.output.get("gff")


# get input files
short_read_bam = snakemake.input.get("short_reads", "")
long_read_bam = snakemake.input.get("long_reads", "")
ref_tx = snakemake.input.get("reference", "")
ptf = snakemake.input.get("point_features", "")


def concatenate_input(files):
    return " ".join(files) if isinstance(files, list) else files


# depending on the input, use different options for StringTie assembly
if short_read_bam and not long_read_bam:
    short_reads = concatenate_input(short_read_bam)
    long_reads = ""
    input_type = "only short reads"
elif not short_read_bam and long_read_bam:
    short_reads = ""
    long_reads = "-L " + concatenate_input(long_read_bam)
    input_type = "only long reads"
elif short_read_bam and long_read_bam:
    short_reads = "--mix " + concatenate_input(short_read_bam)
    long_reads = concatenate_input(long_read_bam)
    input_type = "mixed short and long reads"
else:
    raise ValueError("No input reads provided for StringTie assembly")

# supply reference or point features
if ref_tx:
    ref_tx = "-G " + ref_tx

if ptf:
    ptf = "--ptf " + str(ptf)


# run StringTie
shell(
    "(echo 'Running StringTie assembly with {input_type}';"
    "stringtie "
    "-p {threads} "  # number of threads
    "{extra} "  # additional options
    "{ref_tx} "  # optional reference transcriptome
    "{ptf} "  # optional point features
    "-o {output} "  # output file
    "{short_reads} "  # optional short reads
    "{long_reads}) "  # optional long reads
    "{log}"
)
