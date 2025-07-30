__author__ = "Michael Jahn"
__copyright__ = "Copyright 2025, Michael Jahn"
__email__ = "jahn@mpusp.mpg.de"
__license__ = "MIT"


from os import path
from snakemake.shell import shell

# input/output definitions
counts = snakemake.input["counts"]
meme = snakemake.output["meme"]
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# run nanoplot
shell("cat {counts} | matrix2meme {extra} > {meme} {log}")
