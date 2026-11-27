#!/usr/bin/env bash
set -euo pipefail

if (( $# < 3 )); then
    printf 'Usage: bash tests/run_smoke.sh BIN_DIR NEW_OUTPUT_DIR LAUNCHER [LAUNCHER_ARGS...]\n' >&2
    exit 2
fi
bin=$(readlink -f "$1")
directory=$2
shift 2
launcher=("$@")
mkdir "$directory"
directory=$(readlink -f "$directory")

run() {
    local name=$1
    shift
    printf 'smoke: %s\n' "$name"
    "${launcher[@]}" "$@" </dev/null >"$directory/$name.log" 2>&1
}
app() {
    local name=$1
    shift
    run "$name" "$bin/SWBAM" "$@" --work-memory 8G
}

run make "$bin/swbam-smoke-data" make "$directory"
app flagstat flagstat -i "$directory/input.bam"
app rewrite view -i "$directory/input.bam" -o "$directory/rewrite.bam" \
    --rank-body-backend auto --rank-body-memory-limit 1 --rank-body-temp-dir "$directory"
app filter view -i "$directory/input.bam" -o "$directory/filtered.bam" \
    --min-mapq 30 --rank-body-temp-dir "$directory"
app bam2sam view -i "$directory/input.bam" -o "$directory/text.sam" \
    --io-backend mpiio --rank-body-temp-dir "$directory"
app sam2bam view -i "$directory/text.sam" -o "$directory/from-sam.bam" \
    --rank-body-temp-dir "$directory"
app collate collate -i "$directory/input.bam" -o "$directory/names.bam" \
    -m 4G -n 66 --rank-body-temp-dir "$directory"
app collate-external collate -i "$directory/input.bam" -o "$directory/names-external.bam" \
    -m 192M -n 66 -T "$directory/collate-runs" --rank-body-temp-dir "$directory"
app fixmate fixmate -i "$directory/names.bam" -o "$directory/mates.bam" \
    -m --rank-body-temp-dir "$directory"
app sort sort -i "$directory/mates.bam" -o "$directory/sorted.bam" \
    -m 4G --rank-body-temp-dir "$directory"
app sort-external sort -i "$directory/mates.bam" -o "$directory/sorted-external.bam" \
    -m 192M -T "$directory/sort-runs" --rank-body-temp-dir "$directory"
app markdup markdup -i "$directory/sorted.bam" -o "$directory/marked.bam" \
    -m 4G -l 300 --rank-body-temp-dir "$directory"
app workflow dedup-workflow -i "$directory/input.bam" -o "$directory/workflow.bam" \
    -m 4G -n 66 -T "$directory/workflow" --rank-body-temp-dir "$directory"
run sdk-count "$bin/swbam-sdk-count" "$directory/input.bam" --work-memory 8G
run sdk-rewrite "$bin/swbam-sdk-bam2bam" "$directory/input.bam" \
    -o "$directory/sdk-rewrite.bam" --work-memory 8G
run verify "$bin/swbam-smoke-data" verify "$directory"
grep -q '20000 + 0 in total' "$directory/flagstat.log"
grep -q 'records=20000' "$directory/sdk-count.log"
grep -q 'mode=external' "$directory/collate-external.log"
grep -q 'mode=external' "$directory/sort-external.log"
grep 'SWBAM_SMOKE PASS' "$directory/verify.log"
