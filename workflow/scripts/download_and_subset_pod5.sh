#!/usr/bin/env bash

set -euo pipefail

# Parse command line options
while getopts ":u:n:o:t:h" opt; do
    case $opt in
        u) URL="$OPTARG";;
        n) N_READS="$OPTARG";;
        o) OUTDIR="$OPTARG";;
        t) THREADS="$OPTARG";;
        h) echo "Usage: download_and_subset_pod5.sh -u <URL> -n <N_READS> -o <OUTDIR> -t <THREADS>"; exit 0;;
        \?) echo "Invalid option -$OPTARG" >&2; exit 1;;
    esac
done

mkdir -p "$OUTDIR"

tmpdir=$(mktemp -d)
trap "rm -rf $tmpdir" EXIT

cd "$tmpdir" || exit 1

echo "Downloading pod5 files from $URL..." >&2

aria2c -x "$THREADS" --retry-wait=60 --max-tries=10 --stderr=true -o pod5.tar "$URL"

tar -xf pod5.tar

rm pod5.tar

echo "Download complete. Subsetting $N_READS reads..." >&2

# get the filepath for the file starting with sequencing_summary
summary_file=$(fd -ag1 'sequencing_summary*.txt')
target_mapping_file=read_ids.csv
# the target_mapping_file needs to have the first column as the target filename
# and the second column as the read id - https://github.com/nanoporetech/pod5-file-format/blob/master/python/pod5/README.md#target-mapping-csv
csvtk -t cut -f read_id "$summary_file" |
        sort --random-sort |
        head -n "$N_READS" |
        awk -v n=1 -v OFS=, '{print "batch" int((NR-1)/10)+n ".pod5", $0}' |
        sed '1itarget,read_id' >"$target_mapping_file"

echo "Subsetting the pod5 files..." >&2

pod5 subset --csv "$target_mapping_file" --recursive -o "$OUTDIR" -t "$THREADS" -f .

echo "Subset complete. Cleaning up..." >&2

rm -rf "$tmpdir"

mv "$target_mapping_file" "$OUTDIR"

echo "Done." >&2

