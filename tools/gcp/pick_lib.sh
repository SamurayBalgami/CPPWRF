# Sourced by gcp_run.sh, gcp_matrix.sh and gcp_wrfref.sh (P, BK and REPO set): the VMs to try, cheapest first.
# vm_candidates MACHINE NEED HOURS DATA ZONES prints "MACHINE ZONE DISK_TYPE" lines: tools/gcp/gcp_pick.py's ranking
# of spot price plus data transfer (MACHINE empty: each Intel family's size with NEED vCPUs; ZONES empty: all its
# zones), or, when the price lookup fails, MACHINE (default c3-highcpu with NEED vCPUs) in ZONES or the us zones.
vm_candidates() {
  local m=$1 need=$2 hours=$3 data=$4 zones=$5 gb objs=() d c
  for d in ${data//,/ }; do objs+=("gs://$BK/$d"); done
  gb=$(gcloud storage du -s -c "${objs[@]}" --project "$P" 2>/dev/null | LC_ALL=C awk '$2 == "total" {printf "%.1f", $1 / 1e9}')
  if c=$(python3 "$REPO/tools/gcp/gcp_pick.py" --project "$P" --vcpus "$need" --hours "$hours" --data-gb "${gb:-0}" \
       ${m:+--machine "$m"} ${zones:+--zones "$zones"} --top 15 2> "$OUT/pick.log"); then
    awk '{print $1, $2, $3}' <<< "$c"; return 0
  fi
  echo "[vm_candidates] price lookup failed ($(tail -1 "$OUT/pick.log")): c3 in fixed zones" >&2
  if [ -z "$m" ]; then for c in 4 8 22 44 88 176; do m=c3-highcpu-$c; [ "$need" -le $c ] && break; done; fi
  zones=${zones:-us-central1-a,us-central1-b,us-central1-c,us-central1-f,us-east1-b,us-east1-c,us-east1-d,us-east4-a,us-east4-b,us-east4-c}
  case $m in n4-*|c4-*) d=hyperdisk-balanced ;; *) d=pd-balanced ;; esac
  for c in ${zones//,/ }; do echo "$m $c $d"; done
}
