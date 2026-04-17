__author__ = "Michael Jahn"
__copyright__ = "Copyright 2026, Michael Jahn"
__email__ = "jahn@mpusp.mpg.de"
__license__ = "MIT"


import os
from snakemake.shell import shell

# define default options for this wrapper
threads = snakemake.threads
input_type = snakemake.params.get("input_type", "contig")
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# get input file
fasta = snakemake.input.get("fasta")
if not fasta:
    raise ValueError("Input 'fasta' is required for rgi main.")

# define output prefix; RGI always creates <prefix>.txt and <prefix>.json
report_txt = snakemake.output.get("report_txt")
if not report_txt:
    raise ValueError("Output 'report_txt' is required for this wrapper.")

if input_type not in {"contig", "protein"}:
    raise ValueError("Parameter 'input_type' must be either 'contig' or 'protein'.")

if isinstance(fasta, list):
    if len(fasta) != 1:
        raise ValueError("Exactly one fasta input file is supported by rgi main.")
    fasta = fasta[0]

output_prefix = os.path.splitext(report_txt)[0]
if os.path.dirname(output_prefix):
    os.makedirs(os.path.dirname(output_prefix), exist_ok=True)

# run rgi main
shell(
    "rgi main"
    " --input_sequence {fasta}"
    " --output_file {output_prefix}"
    " --input_type {input_type}"
    " -n {threads}"
    " {extra}"  # additional options
    " {log}"
)
