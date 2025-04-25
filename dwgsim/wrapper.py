__author__ = "Michael Jahn"
__copyright__ = "Copyright 2025, Michael Jahn"
__email__ = "jahn@mpusp.mpg.de"
__license__ = "MIT"


from os import path
from snakemake.shell import shell

# define default options
read_length = snakemake.params.get("read_length", 70)
read_number = snakemake.params.get("read_number", 1e5)
mutations = snakemake.params.get("mutations", 0.001)
random_reads = snakemake.params.get("random_reads", 0.05)
output_type = snakemake.params.get("output_type", 1)
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# get input file
fasta = snakemake.input.get("fasta")
if not fasta:
    raise ValueError("Input FASTA file not specified")
if not path.exists(fasta):
    raise FileNotFoundError(f"Input FASTA file not found: {fasta}")

# determine output based on '-o' option
if output_type == 3 - len(snakemake.output):
    output_prefix = path.commonprefix(snakemake.output)
else:
    raise IOError(
        f"For output type {output_type} exactly {3-output_type} output files must be defined"
    )
pattern_bwa = "<name>.bwa.read<number>.fastq.gz"
pattern_bfast = "<name>.bfast.fastq.gz"

if output_type == 0:
    # separate files for read 1 and read 2, and 1 interleaved file
    if output_prefix.endswith(".b"):
        output_prefix = output_prefix.removesuffix(".b")
    else:
        raise IOError(
            "Three output files must be defined:"
            + " two files named "
            + pattern_bwa
            + " and one file named "
            + pattern_bfast
        )
elif output_type == 1:
    # separate files for read 1 and read 2
    if output_prefix.endswith(".bwa.read"):
        output_prefix = output_prefix.removesuffix(".bwa.read")
    else:
        raise IOError("Two output files must be defined named " + pattern_bwa)
elif output_type == 2:
    # single interleaved file for both reads
    if output_prefix.endswith(".bfast.fastq.gz"):
        output_prefix = output_prefix.removesuffix(".bfast.fastq.gz")
    else:
        raise IOError("One output file must be defined named " + pattern_bfast)
else:
    raise ValueError("Parameter 'output_type' must be one of 0, 1, 2")

# run tool
shell(
    "dwgsim "
    " -1 {read_length}"  # length simulated read 1
    " -2 {read_length}"  # length simulated read 2
    " -N {read_number}"  # total number reads per sample
    " -r {mutations}"  # frequency of mutations
    " -y {random_reads}"  # frequency of random reads
    " -o {output_type}"  # output type, '1' for separate files per read
    " {extra}"  # additional options
    " {fasta}"  # input fasta file
    " {output_prefix}"  # output prefix
    " {log}"
)
