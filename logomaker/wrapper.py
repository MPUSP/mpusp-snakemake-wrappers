__author__ = "Michael Jahn"
__copyright__ = "Copyright 2025, Michael Jahn"
__email__ = "jahn@mpusp.mpg.de"
__license__ = "MIT"


import Bio.SeqIO
import logomaker
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from os import path
from snakemake.shell import shell

# define default options
fasta = snakemake.input.get("fasta", "")
seqlist = snakemake.input.get("seqlist", "")
pwm_file = snakemake.input.get("pwm", "")
pwm_type = snakemake.params.get("pwm_type", "weight")
background = snakemake.params.get("background", None)
ignore = snakemake.params.get("ignore", "")
figsize = snakemake.params.get("figsize", [10, 2.5])
extra = snakemake.params.get("extra", "")
output = snakemake.output
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# parse extra arguments into key-value pairs
logs = []
extra_args = {}
for arg in extra.split():
    k, v = arg.split("=", 1)
    if k in [
        "baseline_width",
        "shade_below",
        "fade_below",
        "vpad",
        "vsep",
        "alpha",
        "zorder",
    ]:
        v = float(v)
    if k in ["center_values", "flip_below", "fade_probabilities", "show_spines", "ax"]:
        v = v.lower() == "true"
    extra_args[k] = v

# process input files
if fasta:
    if not path.exists(fasta):
        raise FileNotFoundError(f"Input FASTA file not found: {fasta}")
    fasta_db = Bio.SeqIO.parse(fasta, "fasta")
    fasta_db = [str(i.seq) for i in fasta_db]
    title = fasta
    logs += [f"Using fasta input from file {fasta}."]
elif seqlist:
    # read input text file to list
    with open(seqlist) as f:
        fasta_db = f.readlines()
        fasta_db = [i.rstrip("\n") for i in fasta_db]
    title = seqlist
    logs += [f"Using list input from file {seqlist}."]
elif pwm_file:
    pwm = pd.read_csv(pwm_file).set_index("pos")
    title = pwm_file
    logs += [f"Using position weight matrix input from file {pwm_file}"]
else:
    raise ValueError(
        "Missing input files. "
        "Exactly one of the three inputs 'fasta', 'pwm' or 'seqlist' must be supplied."
    )

# if not supplied, create the PWM
if not pwm_file:
    # background is relative freq of A,C,G,T
    bg = np.array(background) if background else None
    pwm = logomaker.alignment_to_matrix(
        fasta_db, to_type=pwm_type, background=bg, characters_to_ignore=ignore
    )
else:
    pwm = logomaker.validate_matrix(pwm)

# plot logo
logo = logomaker.Logo(df=pwm, figsize=figsize, **extra_args)
logo.ax.set_title(title, fontsize=11, x=0.01, y=1.02, color="#999999")

# export logo files
for outfile in output:
    plt.savefig(outfile)
    logs += [f"Exported logo to file {outfile}."]

# saved logs
logs = "\n".join(logs)
shell("echo '{logs}' {log}")
