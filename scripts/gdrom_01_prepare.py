from pathlib import Path
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent
DATA = PROJECT_ROOT / "data"

RES_PATH = DATA / "studied_reservoirs.txt"
PROPS_PATH = DATA / "Processed_data" / "reserv_properties.csv"
TS_PATH = DATA / "Processed_data" / "time_series_for_simulations.csv"
PDSI_DAILY_PATH = DATA / "release_policies" / "gdrom" / "fit" / "PDSI" / "reservoir_pdsi_daily.csv"
TRAINING_DIR = DATA / "release_policies" / "gdrom" / "fit" / "GDROM_training"

M3S_TO_MM3DAY = 86400 / 1e6
MM3_TO_ACFT = 810.714

for d in [TRAINING_DIR]:
    d.mkdir(parents=True, exist_ok=True)


def find_col(df, candidates, required=False, label="column"):
    for c in candidates:
        if c in df.columns:
            return c
    if required:
        raise ValueError(f"Could not find {label}. Tried {candidates}. Available: {list(df.columns)}")
    return None


from typing import List

def load_studied_reservoirs(path: Path) -> List[int]:
    if not path.exists():
        raise FileNotFoundError(f"Missing studied reservoirs file: {path}")

    df = pd.read_csv(path, sep=r"\s+", engine="python")

    col = find_col(
        df,
        ["studied_reservoirs", "GRAND_ID"],
        required=True,
        label="studied reservoir id column"
    )

    return df[col].dropna().astype(int).tolist()


def load_inputs():
    if not RES_PATH.exists():
        raise FileNotFoundError(f"Missing {RES_PATH}")
    if not PROPS_PATH.exists():
        raise FileNotFoundError(f"Missing {PROPS_PATH}")
    if not TS_PATH.exists():
        raise FileNotFoundError(f"Missing {TS_PATH}")
    if not PDSI_DAILY_PATH.exists():
        raise FileNotFoundError(
            f"Missing {PDSI_DAILY_PATH}. Copy the PDSI folder into the new project, or create reservoir_pdsi_daily.csv there."
        )

    studied = load_studied_reservoirs(RES_PATH)
    props = pd.read_csv(PROPS_PATH)
    ts = pd.read_csv(TS_PATH, parse_dates=["date"])
    pdsi_daily = pd.read_csv(PDSI_DAILY_PATH, parse_dates=["date"])

    if "GRAND_ID" not in pdsi_daily.columns or "pdsi" not in pdsi_daily.columns:
        raise ValueError(
            f"{PDSI_DAILY_PATH} must contain GRAND_ID, date, pdsi columns. Available: {list(pdsi_daily.columns)}"
        )

    return studied, props, ts, pdsi_daily


def build_metadata(studied, props):
    cap_col = find_col(props, ["Cap_mcm", "capacity_Mm3", "CAP_MCM"], required=True, label="capacity column")
    name_col = find_col(props, ["Name", "DAM_NAME"], required=False, label="name column")

    meta = props.loc[props["GRAND_ID"].isin(studied), ["GRAND_ID"] + ([name_col] if name_col else [])].copy()
    meta["STORAGE_CAP"] = props.loc[props["GRAND_ID"].isin(studied), cap_col].values * MM3_TO_ACFT
    if name_col is None:
        meta["Name"] = ""
    meta = meta.rename(columns={name_col: "Name"} if name_col else {})
    meta.to_csv(TRAINING_DIR / "reservoir_metadata.csv", index=False)
    return cap_col


def build_training_files(studied, props, ts, pdsi_daily, cap_col):
    skipped = []
    ok = []

    ts_cols = set(ts.columns)
    err_col = find_col(ts, ["daily_err_Mm3", "extr", "error_balance_Mm3"], required=False, label="error column")
    if err_col is None:
        raise ValueError(
            "Could not find a daily error column. Expected one of: daily_err_Mm3, extr, error_balance_Mm3"
        )

    for gid in studied:
        prop_row = props.loc[props["GRAND_ID"] == gid]
        if len(prop_row) == 0:
            skipped.append((gid, "not in reserv_properties"))
            continue

        smax_mm3 = float(prop_row[cap_col].iloc[0])
        smax_acft = smax_mm3 * MM3_TO_ACFT

        df = ts.loc[ts["GRAND_ID"] == gid].copy().sort_values("date").reset_index(drop=True)
        if len(df) < 365 * 2:
            skipped.append((gid, f"only {len(df)} days"))
            continue

        pdsi_res = pdsi_daily.loc[pdsi_daily["GRAND_ID"] == gid, ["date", "pdsi"]].copy()
        if len(pdsi_res) == 0:
            skipped.append((gid, "no PDSI series found"))
            continue

        df = df.merge(pdsi_res, on="date", how="left")
        df["pdsi"] = df["pdsi"].ffill().bfill().fillna(0)

        df["Inflow_raw_Mm3"] = df["inflow"] * M3S_TO_MM3DAY
        df["DailyErr_Mm3"] = df[err_col].fillna(0.0).astype(float)
        df["Inflow_net_Mm3"] = df["Inflow_raw_Mm3"] + df["DailyErr_Mm3"]

        df["Storage_Mm3"] = df["storage"].astype(float)
        df["Storage_prev_Mm3"] = df["Storage_Mm3"].shift(1)

        df["Release_Mm3"] = df["Storage_prev_Mm3"] + df["Inflow_net_Mm3"] - df["Storage_Mm3"]
        df["Release_Mm3"] = df["Release_Mm3"].clip(lower=0)

        out = pd.DataFrame({
            "Date": df["date"].dt.strftime("%Y-%m-%d"),
            "Year": df["date"].dt.year,
            "Inflow_raw": (df["Inflow_raw_Mm3"] / smax_mm3).round(6),
            "Inflow": (df["Inflow_net_Mm3"] / smax_mm3).round(6),
            "Storage": (df["Storage_prev_Mm3"] / smax_mm3).clip(0, 1).round(6),
            "Storage_obs": (df["Storage_Mm3"] / smax_mm3).clip(0, 1).round(6),
            "Release": (df["Release_Mm3"] / smax_mm3).round(6),
            "DailyErr": (df["DailyErr_Mm3"] / smax_mm3).round(6),
            "PDSI": df["pdsi"].round(3),
            "DOY": df["date"].dt.dayofyear,
        })

        out = out.iloc[1:].copy()
        out = out.dropna(subset=["Inflow", "Storage", "Storage_obs", "Release", "PDSI", "DOY"])
        out = out[np.isfinite(out["Inflow"]) & np.isfinite(out["Storage"]) & np.isfinite(out["Storage_obs"]) & np.isfinite(out["Release"])]

        if len(out) < 365:
            skipped.append((gid, f"only {len(out)} valid rows"))
            continue

        out.to_csv(TRAINING_DIR / f"{gid}.csv", index=False)
        ok.append(gid)

    print(f"Done: {len(ok)} reservoirs | Skipped: {len(skipped)}")
    for gid, reason in skipped:
        print(f"  GRAND_ID {gid}: {reason}")
    print(f"Training files written to: {TRAINING_DIR}")


if __name__ == "__main__":
    studied, props, ts, pdsi_daily = load_inputs()
    cap_col = build_metadata(studied, props)
    build_training_files(studied, props, ts, pdsi_daily, cap_col)
