__author__ = "Michael Jahn"
__copyright__ = "Copyright 2025, Michael Jahn"
__email__ = "jahn@mpusp.mpg.de"
__license__ = "MIT"


import os
from snakemake.shell import shell

# define default options for this wrapper
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# get input files
summary = snakemake.input.get("summary")
bam = snakemake.input.get("bam", "")


# prepare input
def concatenate_input(files):
    return " ".join(files) if isinstance(files, list) else files


summary = concatenate_input(summary)
bam = concatenate_input(bam) if bam else ""

# prepare output
report_html = snakemake.output.get("report_html")
report_json = snakemake.output.get("report_json")

# run pycioQC
shell(
    "pycoQC"
    " -f {summary}"  # required summary input
    " {bam}"  # optional bam input
    " {extra}"  # additional options
    " -o {report_html}"  # html output
    " -j {report_json}"  # json output
    " {log}"
)
