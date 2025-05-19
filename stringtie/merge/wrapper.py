__author__ = "Michael Jahn"
__copyright__ = "Copyright 2025, Michael Jahn"
__email__ = "jahn@mpusp.mpg.de"
__license__ = "MIT"


from os import path
from snakemake.shell import shell

# define default options for this wrapper
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
output = snakemake.output.get("gff")

# get input files
list_gff = snakemake.input.get("gff")
ref_tx = snakemake.input.get("reference", "")

# supply reference transcriptome
if ref_tx:
    ref_tx = "-G " + ref_tx

# handle edge case of single input file
if not isinstance(list_gff, list):
    list_gff = [list_gff]

# strand definitions
strands = {"plus": "+", "minus": "-", "empty": "."}

# separate files by strand and run StringTie merge
to_combine = []
for strand in strands:
    to_merge = []
    for raw_gff in list_gff:
        temp_gff = path.splitext(raw_gff)[0]
        shell(f"awk '$7 == \"{strands[strand]}\"' {raw_gff} > {temp_gff}_{strand}.gff")
        with open(f"{temp_gff}_{strand}.gff", "r") as f:
            lines = len(f.readlines())
            if lines:
                to_merge.append(f"{temp_gff}_{strand}.gff")
            else:
                shell(f"rm {temp_gff}_{strand}.gff")
    if to_merge:
        input = " ".join(to_merge)
        if len(to_merge) > 1:
            output_temp = path.commonprefix(to_merge)
        else:
            output_temp = path.splitext(to_merge[0])[0]
        shell(
            "(stringtie --merge "  # merge mode
            "{extra} "  # additional options
            "{ref_tx} "  # optional reference transcriptome
            "-o {output_temp}_merged_{strand}.gff "  # output file
            "{input} && "  # list of GFF files
            "rm {input}) "  # remove temp files
            "{log}"
        )
        to_combine.append(f"{output_temp}_merged_{strand}.gff")

# combine results from strand-specific runs
to_combine = " ".join(to_combine)
shell(
    "(cat {to_combine} | "  # strand-specific gffs to combine
    "grep -v '#' | "  # ignore comments
    "sort -k4,4 -n | "  # sort by start position, numerically
    "awk 'gsub(\"MSTRG.[0-9]+\", \"MSTRG.\" int((NR+1)/2))' > "  # re-number tx/exons
    "{output} && "  # output file
    "rm {to_combine}) "  # remove temp files
    "{log}"
)
