__author__ = "Michael Jahn"
__copyright__ = "Copyright 2025, Michael Jahn"
__email__ = "jahn@mpusp.mpg.de"
__license__ = "MIT"


import pandas as pd
from snakemake.shell import shell
from scipy import stats
from Bio import SeqIO
from Bio import motifs

# input/output definitions
fasta = snakemake.input["fasta"]
pwm = snakemake.output["pwm"]
subseq = snakemake.params.get("subseq", "")
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# read fasta file with sequences containing motif
try:
    with open(fasta) as handle:
        records = [r for r in SeqIO.parse(handle, "fasta")]
except:
    raise FileNotFoundError(f"Input file {fasta} not found.")


# check length of fasta file, min 10 records
if len(records) < 10:
    raise ValueError(f"Fasta file {fasta} has less than 10 records.")

# use a subsequence by coordinate (zero based start)
if subseq:
    for r in records:
        r.seq = r.seq[subseq[0] : subseq[1]]

# find the mode of the sequence length and dismiss non-matching items
seq_lengths = [len(r.seq) for r in records]
records = [r for r in records if len(r.seq) == stats.mode(seq_lengths).mode]

# create motif
motif = motifs.create([str(r.seq) for r in records], alphabet="ACGT")

# export counts and pwm matrix dicts to disk
pd.DataFrame(motif.pwm).to_csv(snakemake.output["pwm"], sep="\t", index=False)
pd.DataFrame(motif.counts).to_csv(
    snakemake.output["counts"], sep="\t", index=False, header=False
)

# log message
shell(
    f"echo 'Counts and PWM matrix exported for Fasta with {len(records)} records.' {log}"
)
