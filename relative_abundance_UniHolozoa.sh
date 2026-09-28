#!/bin/bash

set -euo pipefail

###################################################
###################################################

            # From a group of ASVs
    # to relative abundance per environment
         # compared to established unicellular
            # Holozoa lineages

###################################################
###################################################

## Computes total reads per environment for a group of ASVs from EukBank.
#
# Input files must follow the naming pattern ASVs_GroupX.list, where each
# file lists the ASVs belonging to one group. Reads are summed per
# environment for the group as a whole, and separately for each characterized
# unicellular Holozoa lineage, using EukBank taxonomy assignments.

# Note: Atreyea and Filasterea are treated as the same group.

# Note: designed for EukBank v1. Column indices assume the v1 
# file structure and may need updating for other versions.

# WARNING: If an ASV is assigned to both a .list group and a characterized
# unicellular Holozoa taxon, it is counted only under the .list group.

## Usage:
#   ./script.sh <input_folder> <output_folder>

## Arguments:
#   input_folder    Folder containing the .list files
#   output_folder   Directory where results will be written

# Output:
#   UniHolozoa_AllGroups_reads_per_env.tsv
#   Columns: taxon_or_group, environment, summed_reads
#
#   Intermediate files are also written to output_folder:
#     UniHolozoa_eukbank_18S_V4_ASVtotaxon.tsv
#     Groups_ASVtotaxon.tsv
#     UniHolozoa_AllGroups_ASVtotaxon.tsv
#     UniHolozoa_AllGroups_counts.tsv
#     UniHolozoa_AllGroups_counts-envplot.tsv

## Example:
#   ./script.sh EukBank/ASV_lists EukBank/Relative_abundance_UniHolozoa

# Require two arguments:
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_folder_with_lists> <output_folder>"
    exit 1
fi

# Folder containing ASVs_Group*.list files:
INPUT_DIR="$1"

# Output folder:
OUTDIR="$2"
mkdir -p "$OUTDIR"

# Hardcoded input files — verify these paths before running: 
ASVS_FILE="eukbank_18S_V4_asvs.tsv"
COUNTS_FILE="eukbank_18S_V4_counts.tsv"
SAMPLES_FILE="eukbank_18S_V4_samples.tsv"


# 1. Filter ASVs belonging to unicellular Holozoa lineages (ASV, taxon):
UNI_FILE="$OUTDIR/UniHolozoa_eukbank_18S_V4_ASVtotaxon.tsv"

awk -F'\t' '$10 ~ /Choanoflagellatea|Filasterea|Ichthyosporea|Corallochytrea|Tunicaraptor/ {print $1, $10}' OFS='\t' \
    "$ASVS_FILE" > "$UNI_FILE"


# 2. Build an ASV-to-group mapping file (ASV, GroupX) from the .list files:
GROUP_FILE="$OUTDIR/Groups_ASVtotaxon.tsv"
> "$GROUP_FILE"

for f in "$INPUT_DIR"/ASVs_Group*.list; do
    [ -e "$f" ] || continue
    GROUP=$(basename "$f" .list)
    GROUP=${GROUP#ASVs_}
    awk -v group="$GROUP" '{print $1 "\t" group}' "$f" >> "$GROUP_FILE"
done

# 3. Combine unicellular Holozoa taxonomy and ASV group assignments:
FINAL_ASV_FILE="$OUTDIR/UniHolozoa_AllGroups_ASVtotaxon.tsv"
cat "$GROUP_FILE" "$UNI_FILE" > "$FINAL_ASV_FILE"
echo "Combined ASV + Taxon/Group file:"
echo "$FINAL_ASV_FILE"

# 4. Generate per-ASV counts:
OUTPUT_COUNTS="$OUTDIR/UniHolozoa_AllGroups_counts.tsv"

awk -F'\t' -v OFS='\t' '
    # First file (GROUP_FILE): store ASV-to-group assignments
    FNR==NR {
        group_asv[$1] = $2
        next
    }

    # Second file (UNI_FILE): fill in Holozoa taxon only if not already assigned to a group
    FILENAME==uni_file {
        if (!($1 in group_asv)) {
            group_asv[$1] = $2
        }
        next
    }

    # Third file (COUNTS_FILE): print counts prefixed with the resolved taxon/group
    {
        taxon = group_asv[$1]
        if (taxon == "Tunicaraptor") taxon="Filasterea"
        if (taxon != "") print taxon, $0
    }

' "$GROUP_FILE" uni_file="$UNI_FILE" "$UNI_FILE" "$COUNTS_FILE" > "$OUTPUT_COUNTS"

echo "Per-ASV counts file generated (GROUP priority applied):"
echo "$OUTPUT_COUNTS"

# 5. Add environment column:
OUTPUT_ENV="$OUTDIR/UniHolozoa_AllGroups_counts-envplot.tsv"

awk -F'\t' -v OFS='\t' '
    NR==FNR { sample2env[$1]=$15; next }   # Map SampleID to environment
    {
        sample = $3
        env = sample2env[sample] ? sample2env[sample] : "Unknown"
        print $0, env
    }
' "$SAMPLES_FILE" "$OUTPUT_COUNTS" > "$OUTPUT_ENV"
echo "Environment column added:"
echo "$OUTPUT_ENV"

# 6. Sum reads per taxon/group per environment:
FINAL_SUM="$OUTDIR/UniHolozoa_AllGroups_reads_per_env.tsv"

awk -F'\t' -v OFS='\t' '
{
    # $1 = taxon/group
    # $5 = environment
    reads[$1,$5] += $4
}
END {
    # Output: taxon/group, environment, summed reads
    for (key in reads) {
        split(key, arr, SUBSEP)
        print arr[1], arr[2], reads[key]
    }
}
' "$OUTPUT_ENV" > "$FINAL_SUM"

echo "Final table of reads per taxon/group per environment generated:"
echo "$FINAL_SUM"