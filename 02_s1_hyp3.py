"""Step 02: Sentinel-1 interferograms on demand via ASF HyP3 (free, Earthdata login).

Builds a small-baseline network inside each summer (each date paired with
the next 1-2 acquisitions), so MintPy can invert it to a time series.
Uses INSAR_ISCE_MULTI_BURST on the bursts listed by step 01 (s1_bursts.csv),
so frame-number drift between years does not matter. Only dates with the
full burst set are used. Cost comes from the public HyP3 /costs table
(10x2 looks = 40 m; 4-9 bursts = 5 credits/pair). Free quota 8,000/month.

    python 02_s1_hyp3.py plan     --path P [--years 2020 2025]
    python 02_s1_hyp3.py submit   --path P [--years ...]
    python 02_s1_hyp3.py status   --path P
    python 02_s1_hyp3.py download --path P
Start small (e.g. --years 2020 2025), check one product, then add years.
"""
import argparse
import zipfile
from itertools import islice

import pandas as pd
import requests

from config import DIRS, YEARS

LOOKS = "10x2"
COSTS_URL = "https://hyp3-api.asf.alaska.edu/costs"


def credits_per_pair(n_bursts, looks=LOOKS):
    table = requests.get(COSTS_URL, timeout=30).json()
    return table["INSAR_ISCE_MULTI_BURST"]["cost_table"][looks][str(n_bursts)]


def build_pairs(dates, n_next=2, annual=False):
    """dates: list of (date, granules) where granules is a tuple of burst
    names for that date. Returns list of (g1, g2, d1, d2).

    Inside each year: pair each scene with the next `n_next` scenes.
    annual=True also links consecutive summers (late-summer to late-summer)
    so the network is connected across years; C-band coherence over one year
    is often poor on tundra, decent on built-up/rocky ground.
    """
    df = pd.DataFrame(dates, columns=["date", "granule"])
    df["date"] = pd.to_datetime(df["date"])
    df = df.drop_duplicates("date").sort_values("date")
    pairs = []
    for _, g in df.groupby(df["date"].dt.year):
        g = g.reset_index(drop=True)
        for i in range(len(g)):
            for j in range(i + 1, min(i + 1 + n_next, len(g))):
                pairs.append((g.granule[i], g.granule[j], g.date[i], g.date[j]))
    if annual:
        years = sorted(df["date"].dt.year.unique())
        for y0, y1 in zip(years, years[1:]):
            a = df[df["date"].dt.year == y0].iloc[-1]
            target = a.date + pd.DateOffset(years=y1 - y0)
            b_all = df[df["date"].dt.year == y1]
            b = b_all.iloc[(b_all.date - target).abs().argsort().iloc[0]]
            if abs((b.date - target).days) <= 15:
                pairs.append((a.granule, b.granule, a.date, b.date))
    return pairs


def chunks(it, n):
    it = iter(it)
    while batch := list(islice(it, n)):
        yield batch


def job_name(path):
    return f"norilsk_p{path}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["plan", "submit", "status", "download"])
    ap.add_argument("--path", type=int, required=True)
    ap.add_argument("--years", type=int, nargs="+", default=YEARS)
    ap.add_argument("--n-next", type=int, default=2)
    ap.add_argument("--annual", action="store_true")
    a = ap.parse_args()

    if a.cmd in ("plan", "submit"):
        b = pd.read_csv(DIRS["s1"] / "s1_bursts.csv")
        b = b[b.burst_id.str.startswith(f"{a.path}_") & b.complete
              & b.date.str[:4].astype(int).isin(a.years)]
        if b.empty:
            raise SystemExit("No complete burst dates: run 01 with --path.")
        n_bursts = b.burst_id.nunique()
        by_date = (b.sort_values("burst_id").groupby("date").granule
                   .apply(tuple).reset_index().values.tolist())
        pairs = build_pairs(by_date, a.n_next, a.annual)
        plan = pd.DataFrame(pairs, columns=["ref", "sec", "ref_date", "sec_date"])
        plan["dt_days"] = (plan.sec_date - plan.ref_date).dt.days
        plan[["ref_date", "sec_date", "dt_days"]].to_csv(
            DIRS["hyp3"] / f"pairs_p{a.path}.csv", index=False)
        cpp = credits_per_pair(n_bursts)
        print(plan.groupby(plan.ref_date.dt.year).size().rename("pairs").to_string())
        print(f"\n{n_bursts} bursts/pair, {LOOKS} looks, {cpp} credits/pair")
        print(f"{len(plan)} pairs, {len(plan) * cpp} credits")
        if a.cmd == "plan":
            return

    import hyp3_sdk
    from auth import earthdata_credentials
    user, pwd = earthdata_credentials()
    hyp3 = hyp3_sdk.HyP3(username=user, password=pwd)
    name = job_name(a.path)

    if a.cmd == "submit":
        print(f"Credits available: {hyp3.check_credits()}")
        done = {(tuple(sorted(j.job_parameters["reference"])),
                 tuple(sorted(j.job_parameters["secondary"])))
                for j in hyp3.find_jobs(name=name)
                if j.job_type == "INSAR_ISCE_MULTI_BURST"}
        todo = [p for p in pairs
                if (tuple(sorted(p[0])), tuple(sorted(p[1]))) not in done]
        print(f"Already submitted: {len(pairs) - len(todo)}; new: {len(todo)}")
        if input("Submit? [y/N] ").lower() != "y":
            return
        prepared = [hyp3_sdk.HyP3.prepare_insar_isce_multi_burst_job(
            list(g1), list(g2), name=name, looks=LOOKS,
            apply_water_mask=True) for g1, g2, *_ in todo]
        for batch in chunks(prepared, 50):
            hyp3.submit_prepared_jobs(batch)
        print("Submitted. Processing takes ~1 h per batch; run 'status'.")

    elif a.cmd == "status":
        jobs = hyp3.find_jobs(name=name)
        print(pd.Series([j.status_code for j in jobs]).value_counts().to_string())

    elif a.cmd == "download":
        out = DIRS["hyp3"] / name
        out.mkdir(exist_ok=True)
        jobs = hyp3.find_jobs(name=name).filter_jobs(
            succeeded=True, pending=False, running=False, failed=False,
            include_expired=False)
        have = {p.stem for p in out.iterdir()}
        new = [j for j in jobs
               if j.files[0]["filename"].removesuffix(".zip") not in have]
        print(f"{len(jobs)} succeeded, {len(new)} to download -> {out}")
        for j in new:
            for z in j.download_files(out):
                with zipfile.ZipFile(z) as zf:
                    zf.extractall(out)
                z.unlink()


if __name__ == "__main__":
    main()
