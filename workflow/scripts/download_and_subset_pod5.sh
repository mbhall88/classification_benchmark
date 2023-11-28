exec 2> "${snakemake_log[0]}"

URL="${snakemake_params[url]}"
N_READS="${snakemake_params[n_reads]}"
OUTDIR="${snakemake_output[pod5_dir]}"
mkdir -p "$OUTDIR"

tmpdir=$(mktemp -d)
trap "rm -rf $tmpdir" EXIT

cd $tmpdir || exit 1

echo "Downloading pod5 files from $URL..." >&2

curl -L --retry 10 --retry-delay 60 "$URL" | tar -xz

echo "Download complete. Subsetting $N_READS reads..." >&2

# get the filepath for the file starting with sequencing_summary
summary_file=$(fd -ag1 'sequencing_summary*.txt')
target_mapping_file=read_ids.csv
# the target_mapping_file needs to have the first column as the target filename 
# and the second column as the read id - https://github.com/nanoporetech/pod5-file-format/blob/master/python/pod5/README.md#target-mapping-csv
csvtk -t cut -f read_id "$summary_file" \
    | sort --random-sort \
    | head -n $N_READS \
    | awk -v n=1 -v OFS=, '{print "batch" int((NR-1)/10)+n ".pod5", $0}' \
    | sed '1itarget,read_id' > "$target_mapping_file"

echo "Subsetting the pod5 files..." >&2

pod5 subset --csv "$target_mapping_file" --recursive -o "$OUTDIR" -t "${snakemake[threads]}" -f .

echo "Subset complete. Cleaning up..." >&2

rm -rf "$tmpdir" "$target_mapping_file"

echo "Done." >&2