#!/usr/bin/env python3
"""Rank spot machine/zone pairs by the expected cost of a VM job.

    gcp_pick.py --vcpus N [--hours H] [--data-gb G] [--machine TYPE] [--families c3,n4,...] [--zones Z,...] [--top K]

Prints "MACHINE ZONE DISK_TYPE VCPUS COST_USD" lines, cheapest first: the order in which gcp_run.sh,
gcp_matrix.sh and gcp_wrfref.sh try to create their VM (spot capacity and quotas vary, so the next line is
the fallback). The job needs N vCPUs for H hours; a machine with fewer vCPUs takes longer in proportion
(the smallest size with N vCPUs, else the family's largest). The cost is the spot price of the
machine's vCPUs and memory (Cloud Billing Catalog, current prices, cached for 6 hours in
~/.cache/cppwrf_gcp_prices.json) over the job's hours plus the transfer of G GB of data from the bucket's
region (us-central1) to the VM's region. With CPPWRF_GCP_PROJECT (or --project) set, a zone without the
machine type (listed once a week) and a region without the vCPU quota left for it are left out. Only Intel families: glibc's libmvec differs on AMD, so AMD runs
are not bit-identical to the Intel-made WRF references. --machine TYPE ranks only the zones for that type,
--zones only the zones given.
"""
import argparse, json, os, re, subprocess, sys, time, urllib.request

BUCKET_REGION = "us-central1"
REGIONS = {  # zones tried; europe costs more transfer, so it wins only when its spot price is far lower
    "us-central1": "abcf", "us-east1": "bcd", "us-east4": "abc", "us-west1": "abc",
    "europe-west1": "bcd", "europe-west4": "abc",
}
TRANSFER = {"us": 0.02, "europe": 0.05}  # USD per GB from a us-central1 bucket (same region: free)
# family: (vCPU sizes, GB per vCPU, catalog name, boot disk type, speed of a vCPU relative to c3)
FAMILIES = {
    "c3": ("highcpu", (4, 8, 22, 44, 88, 176), 2.0, "C3", "pd-balanced", 1.0),
    "n4": ("highcpu", (4, 8, 16, 32, 48, 64, 80), 2.0, "N4", "hyperdisk-balanced", 1.0),
    "c4": ("highcpu", (4, 8, 16, 32, 48, 96, 144, 192), 1.875, "C4", "hyperdisk-balanced", 1.1),
    "h3": ("standard", (88,), 4.0, "H3", "pd-balanced", 1.5),  # no SMT: a vCPU is a whole core
}
CACHE = os.path.expanduser("~/.cache/cppwrf_gcp_prices.json")
COMPUTE = "6F81-5844-456A"


def spot_prices():
    """{(CATALOG_NAME, region): {"Core": usd/h, "Ram": usd/GB-h}} for spot instances."""
    if os.path.exists(CACHE) and time.time() - os.path.getmtime(CACHE) < 6 * 3600:
        with open(CACHE) as f:
            return {tuple(k.split("|")): v for k, v in json.load(f).items()}
    tok = subprocess.run(["gcloud", "auth", "print-access-token"], capture_output=True, text=True,
                         check=True).stdout.strip()
    names = {v[3] for v in FAMILIES.values()}
    out, page = {}, ""
    while True:
        url = (f"https://cloudbilling.googleapis.com/v1/services/{COMPUTE}/skus?pageSize=5000&currencyCode=USD"
               + (f"&pageToken={page}" if page else ""))
        with urllib.request.urlopen(urllib.request.Request(url, headers={"Authorization": f"Bearer {tok}"}),
                                    timeout=60) as r:
            d = json.load(r)
        for s in d.get("skus", []):
            m = re.match(r"Spot Preemptible (\w+) Instance (Core|Ram) running", s["description"])
            if not m or m.group(1) not in names:
                continue
            u = s["pricingInfo"][0]["pricingExpression"]["tieredRates"][-1]["unitPrice"]
            for reg in s.get("serviceRegions", []):
                if reg in REGIONS:
                    out.setdefault((m.group(1), reg), {})[m.group(2)] = int(u["units"]) + u["nanos"] / 1e9
        page = d.get("nextPageToken", "")
        if not page:
            break
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    with open(CACHE, "w") as f:
        json.dump({"|".join(k): v for k, v in out.items()}, f)
    return out


def zones_having(names, project):
    """{machine type: set of zones} in REGIONS (cached for a week: types rarely change)."""
    cache = CACHE.replace("prices", "zones")
    if os.path.exists(cache) and time.time() - os.path.getmtime(cache) < 7 * 86400:
        with open(cache) as f:
            got = json.load(f)
        if all(n in got for n in names):
            return {n: set(z) for n, z in got.items()}
    regs = "|".join(REGIONS)
    out = subprocess.run(["gcloud", "compute", "machine-types", "list", "--project", project, "--filter",
                          f"name~'^({'|'.join(names)})$' AND zone~'^({regs})-'", "--format=value(name,zone)"],
                         capture_output=True, text=True, check=True).stdout
    got = {n: [] for n in names}
    for line in out.split("\n"):
        if line.strip():
            n, z = line.split()
            got[n].append(z.rsplit("/", 1)[-1])
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    with open(cache, "w") as f:
        json.dump(got, f)
    return {n: set(z) for n, z in got.items()}


def free_quota(project):
    """{(region, METRIC): limit - usage} now; a family over its quota there cannot be created."""
    out = subprocess.run(["gcloud", "compute", "regions", "list", "--project", project, "--format=json"],
                         capture_output=True, text=True, check=True).stdout
    return {(r["name"], q["metric"]): q["limit"] - q["usage"] for r in json.loads(out) if r["name"] in REGIONS
            for q in r.get("quotas", [])}


def machine(fam, need, only):
    """--machine's type, or the smallest size with N vCPUs (the largest if none has)."""
    kind, sizes, gb, *_ = FAMILIES[fam]
    if only:
        return only, int(only.rsplit("-", 1)[1]), gb
    v = min([v for v in sizes if v >= need] or [max(sizes)])
    return f"{fam}-{kind}-{v}", v, gb


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--vcpus", type=int, required=True)
    ap.add_argument("--hours", type=float, default=1.0)
    ap.add_argument("--data-gb", type=float, default=0.0)
    ap.add_argument("--machine", default="")
    ap.add_argument("--families", default=",".join(FAMILIES))
    ap.add_argument("--zones", type=lambda z: z.split(","), default=[])
    ap.add_argument("--top", type=int, default=12)
    ap.add_argument("--project", default=os.environ.get("CPPWRF_GCP_PROJECT", ""))
    a = ap.parse_args()
    fams = [a.machine.split("-")[0]] if a.machine else a.families.split(",")
    if any(f not in FAMILIES for f in fams):
        sys.exit(f"gcp_pick: families are {','.join(FAMILIES)} (Intel only)")
    prices = spot_prices()
    pick = {fam: machine(fam, a.vcpus, a.machine) for fam in fams}
    where = zones_having([m[0] for m in pick.values()], a.project) if a.project else {}
    quota = free_quota(a.project) if a.project else {}
    rows = []
    for fam in fams:
        _, _, _, cat, disk, speed = FAMILIES[fam]
        name, v, gb = pick[fam]
        hours = a.hours * a.vcpus / (min(v, a.vcpus) * speed)
        for reg, zs in REGIONS.items():
            p = prices.get((cat, reg), {})
            if "Core" not in p or "Ram" not in p:
                continue
            if quota.get((reg, f"{cat}_CPUS"), v) < v or quota.get((reg, "CPUS"), v) < v:
                continue  # the family's (or the region's) vCPU quota is used up there
            move = 0.0 if reg == BUCKET_REGION else TRANSFER[reg.split("-")[0]] * a.data_gb
            cost = v * (p["Core"] + gb * p["Ram"]) * hours + move
            rows += [(cost, name, f"{reg}-{z}", disk, v) for z in zs
                     if (not where or f"{reg}-{z}" in where.get(name, ())) and (not a.zones or f"{reg}-{z}" in a.zones)]
    if not rows:
        sys.exit("gcp_pick: no priced machine")
    rows.sort()
    for cost, name, zone, disk, v in rows if a.machine else rows[:a.top]:
        print(f"{name} {zone} {disk} {v} {cost:.2f}")


if __name__ == "__main__":
    main()
