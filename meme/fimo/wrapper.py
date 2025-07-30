__author__ = "Michael Jahn"
__copyright__ = "Copyright 2025, Michael Jahn"
__email__ = "jahn@mpusp.mpg.de"
__license__ = "MIT"


from os import path
from snakemake.shell import shell

# input/output definitions
meme = snakemake.input["meme"]
fasta = snakemake.input["fasta"]
output_prefix = path.commonprefix(snakemake.output).rstrip(".")
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# check output file name
if not output_prefix.endswith("fimo"):
    raise ValueError("Output 'dir/filename' must be defined like <mydir>/fimo")

# run fimo search
shell(
    "(fimo"
    " {extra}"
    " -oc {output_prefix}"
    " {meme}"
    " {fasta}"
    " && "
    "mv {output_prefix}/* {output_prefix}/../ &&"
    "rm -rf {output_prefix}/) {log}"
)
