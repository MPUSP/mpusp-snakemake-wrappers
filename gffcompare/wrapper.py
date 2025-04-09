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
output_prefix = os.path.commonprefix(snakemake.output).rstrip(".")


# get input files
gff = snakemake.input.get("gff")
ref = snakemake.input.get("reference", "")
fasta = snakemake.input.get("fasta", "")

# prepare input
if isinstance(gff, list):
    gff = " ".join(gff)

if ref:
    ref = f"-r {ref}"

if fasta:
    fasta = f"-s {fasta}"


# run gffcompare tool with the specified options
shell(
    "gffcompare"
    " {ref}"  # optional reference input
    " -p {threads}"  # number of threads
    " {extra}"  # additional options
    " {fasta}"  # optional fasta input
    " -o {output_prefix}"  # output prefix
    " {gff}"
    " {log}"
)
