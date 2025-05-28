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
barcode = snakemake.input.get("barcode", "")


# prepare input
def concatenate_input(prefix, files):
    if files:
        files = " ".join(files) if isinstance(files, list) else files
        return f"--{prefix} " + files
    else:
        if prefix == "summary_file":
            raise ValueError("'summary' input is required for pycoQC.")
        return ""


summary = concatenate_input("summary_file", summary)
bam = concatenate_input("bam_file", bam)
barcode = concatenate_input("barcode_file", barcode)

# prepare output
report_html = snakemake.output.get("report_html")
report_json = snakemake.output.get("report_json")

# run pycioQC
shell(
    "pycoQC"
    " {summary}"  # required summary input
    " {bam}"  # optional bam input
    " {barcode}"  # optional barcode input
    " {extra}"  # additional options
    " -o {report_html}"  # html output
    " -j {report_json}"  # json output
    " {log}"
)
