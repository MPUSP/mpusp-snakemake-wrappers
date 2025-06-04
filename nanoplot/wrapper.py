__author__ = "Michael Jahn"
__copyright__ = "Copyright 2025, Michael Jahn"
__email__ = "jahn@mpusp.mpg.de"
__license__ = "MIT"


from os import path
from snakemake.shell import shell

# define default options for this wrapper
threads = snakemake.threads
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# get input files, only one input is allowed
input = {}
input_args = {
    "summary": "",
    "fastq": "",
    "fasta": "",
    "fastq_rich": "",
    "fastq_minimal": "",
    "bam": "",
    "ubam": "",
    "cram": "",
    "pickle": "",
    "feather": "",
    "arrow": "",
}

for arg in input_args:
    input_args[arg] = snakemake.input.get(arg, "")
    if input_args[arg]:
        input[arg] = input_args[arg]

# test if more than one input was provided
if len(input) != 1:
    raise ValueError(
        "Exactly one of the following inputs is allowed: "
        + ", ".join(input_args.keys())
        + f"\nYou have used {len(input)} inputs: "
        + ", ".join(input.keys())
    )

# format input string
for k, v in input.items():
    files = " ".join(v) if isinstance(v, list) else v
    input_formatted = f"--{k} " + files

# prepare output
report = snakemake.output.get("report")
output_dir = path.abspath(path.dirname(report))

# run nanoplot
shell(
    "NanoPlot"
    " {input_formatted}"  # input can be summary, bam, fasta, fastq, ...
    " {extra}"  # additional options
    " --outdir {output_dir}"  # output dir
    " --verbose"  # needed to redirect log to snakemake logfile
    " --threads {threads}" # number of threads used
    " {log};"
    "mv {output_dir}/NanoPlot-report.html {report}"  # move to desired output file name
)
