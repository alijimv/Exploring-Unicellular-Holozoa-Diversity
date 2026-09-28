#!/bin/bash

set -euo pipefail

###################################################
###################################################

            # From ASVs to environmental
            # distribution and biogeography

###################################################
###################################################

## Extracts sample-level metadata for a given list of ASVs from EukBank.
#
# For each ASV, the script identifies the samples in which it occurs
# and retrieves associated data, including read counts, environment,
# and geographic coordinates. The resulting table can be used
# downstream, e.g. to generate biogeographic maps in R.

# Note: designed for EukBank v1. Column indices (e.g. environment
# in $SAMPLES) assume the v1 file structure and may need updating
# for other versions.

## Usage:
# ./envdistribution.sh <files.list> ... <output_folder>

## Arguments:
#   files.list      One or more ASV list files (e.g. *.list)
#   output_folder   Directory where results will be written

## Output:
#   One TSV file per input group, named:
#     <basename>_reads_per_sample_nreads_env_lat_long_relativeabundance.tsv
#   Columns: sample, reads, nreads, environment, lat, lon, relative_abundance

## Example:
#   ./envdistribution.sh *.list EukBank/Biogeography


# Require at least two arguments:
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <files_list> ... <outdir>"
    exit 1
fi

# Use the last argument as the output directory:
OUTDIR="${!#}"

# Set all arguments except the last one as the list of files to process:
LIST_FILES=("${@:1:$#-1}")

# Check that all files exist:
for f in "${LIST_FILES[@]}"; do
    if [ ! -f "$f" ]; then
        echo "File not found: $f"
        exit 1
    fi
done

# Create the output directory if it doesn't exist:
mkdir -p "$OUTDIR"

# Hardcoded input files — verify these paths before running: 
COUNTS="eukbank_18S_V4_counts.tsv"
SAMPLES="eukbank_18S_V4_samples.tsv"
ENVLIST="$OUTDIR/envplot.list"

# Create a list of environments:
cut -f15 "$SAMPLES" | sort | uniq | tail -n 9 > "$ENVLIST"

# ---------------------------------------------------------------------
# Process a single group: full pipeline from ASV list to output table.
# ---------------------------------------------------------------------

process_group() { 
    local LIST="$1"
    local BASENAME
    BASENAME=$(basename "$LIST" .list)
    local GDIR="$OUTDIR/$BASENAME"
    mkdir -p "$GDIR"

    echo "Processing $BASENAME ..."

    # 1. Extract read counts per sample for the given ASVs, and append environment data:
    awk -F '\t' '
    NR==FNR { asv[$1]=1; next }
    FNR==1 {
        while ((getline < samples) > 0)
            env[$1]=$15
    }
    $1 in asv {
        printf "%s\t%s\n", $0, env[$2]
    }
    ' samples="$SAMPLES" "$LIST" "$COUNTS" \
    > "$GDIR/${BASENAME}_counts_env.tsv"

    # 2. Sum number of reads per environment for each ASV
    awk -F '\t' '
    {
        sum[$1,$4] += $3
        asv[$1]=1
        env[$4]=1
    }
    END {
        for (k in sum) {
            split(k,a,SUBSEP)
            printf "%s\t%s\t%d\n", a[1], a[2], sum[k]
        }
    }
    ' "$GDIR/${BASENAME}_counts_env.tsv" \
    > "$GDIR/${BASENAME}_ASV_env_long.tsv"

    # 3. Create the wide-format file (ASVs as rows, environments as columns):
    awk -F '\t' '
    {
        val[$1,$2]=$3
        asv[$1]=1
        env[$2]=1
    }
    END {
        printf "ASV"
        for (e in env) printf ",%s", e
        print ""
        for (a in asv) {
            printf "EukBank_%s", a
            for (e in env)
                printf ",%d", val[a,e]+0
            print ""
        }
    }
    ' "$GDIR/${BASENAME}_ASV_env_long.tsv" \
    > "$GDIR/${BASENAME}_ASV_env_WIDE.csv"

    # 4. Sum reads per sample and append sample metadata (nreads, environment,
    #    lat, lon, relative abundance):
    awk -F '\t' '
    NR==FNR {
        nreads[$1]=$2
        env[$1]=$15
        lat[$1]=$6
        lon[$1]=$7
        next
    }
    {
        reads[$2]+=$3
    }
    END {
        for (s in reads)
            printf "%s\t%d\t%s\t%s\t%s\t%s\t%f\n",
                   s, reads[s], nreads[s], env[s], lat[s], lon[s], reads[s]/nreads[s]
    }
    ' "$SAMPLES" "$GDIR/${BASENAME}_counts_env.tsv" \
    > "$GDIR/${BASENAME}_reads_per_sample_nreads_env_lat_long_relativeabundance.tsv"

    printf "Group %s complete\n" "$BASENAME"
}

# ------------------------------
# Process each group in parallel
# ------------------------------

for L in "${LIST_FILES[@]}"; do
    process_group "$L" &
done

wait  # Wait for all groups to finish

printf "Pipeline complete. Results in $OUTDIR"


