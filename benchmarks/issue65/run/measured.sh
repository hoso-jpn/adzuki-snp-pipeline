#!/usr/bin/env bash
# measured.sh <label> <docker run args...>: run one container detached and record
# wall time, exit code and the container cgroup's memory.peak (polled each second,
# so a lower bound within the last second of the run) to <label>.resource.tsv.
set -uo pipefail
label=$1; shift
start=$(date +%s.%N)
cid=$(docker run -d "$@") || exit 1
peak=0
while [ "$(docker inspect -f '{{.State.Running}}' "$cid")" = true ]; do
  f=/sys/fs/cgroup/system.slice/docker-$cid.scope/memory.peak
  if [ -r "$f" ]; then v=$(cat "$f" 2>/dev/null || echo 0); [ "${v:-0}" -gt "$peak" ] && peak=$v; fi
  sleep 1
done
code=$(docker wait "$cid")
end=$(date +%s.%N)
docker logs "$cid" > "$label.container.log" 2>&1
docker rm "$cid" > /dev/null
printf 'label\texit_code\twall_seconds\tmemory_peak_bytes\n%s\t%s\t%.1f\t%s\n' "$label" "$code" "$(echo "$end - $start" | bc)" "$peak" > "$label.resource.tsv"
exit "$code"
