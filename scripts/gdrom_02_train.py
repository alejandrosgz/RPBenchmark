import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import GridSearchCV
from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor, _tree

warnings.filterwarnings("ignore")

PROJECT_ROOT = Path(__file__).resolve().parent
DATA = PROJECT_ROOT / "data"
TRAINING_DIR = DATA / "release_policies" / "gdrom" / "fit" / "GDROM_training"
OUT_DIR = DATA / "release_policies" / "gdrom" / "fit" / "gdrom_fitted"
RESULTS_DIR = OUT_DIR / "results"
RULES_DIR = OUT_DIR / "rules"

GDROM_ENV = Path(r"C:\ASG\UCDavis\NSF_PROJECT_FOLDER\Datasets\GDROM_v2\data\contents\Scripts\Environment")
sys.path.insert(0, str(GDROM_ENV))

try:
    from hmmlearn import hmm_DT
    HMM_OK = True
    print("hmm_DT loaded OK")
except Exception as e:
    HMM_OK = False
    print(f"hmm_DT not available: {e}")

CALIB_END = 2004
VAL_END = 2010
TEST_END = 2015
ITS_VALUES = [6, 7, 8]
TEST_IDS = None

for d in [RESULTS_DIR, RULES_DIR / "modules", RULES_DIR / "module_conditions"]:
    d.mkdir(parents=True, exist_ok=True)


def load_training_file(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, parse_dates=["Date"])
    required = ["Date", "Year", "Inflow", "Storage", "Storage_obs", "Release", "PDSI", "DOY"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in {path.name}: {missing}")
    return df


def nrmse_sd(obs, sim):
    obs = np.asarray(obs, dtype=float)
    sim = np.asarray(sim, dtype=float)
    ok = np.isfinite(obs) & np.isfinite(sim)
    if ok.sum() < 2:
        return np.inf
    sd = np.std(obs[ok], ddof=1)
    if not np.isfinite(sd) or sd == 0:
        return np.inf
    rmse = np.sqrt(np.mean((obs[ok] - sim[ok]) ** 2))
    return float(rmse / sd)


def train_hmdt(train_data, impurity):
    ITS_min = 10 ** (-impurity)
    O = train_data["Release"].values.reshape(-1, 1)
    F = train_data[["Inflow", "Storage"]].values
    tree_model = DecisionTreeRegressor(
        min_impurity_decrease=ITS_min,
        random_state=41,
        max_depth=15,
        min_samples_leaf=10,
    )
    tree_model.fit(F, O)
    GHMM = hmm_DT.GaussianHMM
    models = []
    for k in range(1, 8):
        try:
            m = GHMM(tree_model, relax="all", n_components=k, verbose=False, n_iter=200, trials=10)
            m.fit(O, F, lengths=None)
            models.append((k, m.best_model))
        except Exception as e:
            print(f"      K={k} error: {e}")
    return models


def train_cart(train_data, hmdt_model, n_states):
    O = train_data["Release"].values.reshape(-1, 1)
    F = train_data[["Inflow", "Storage"]].values
    hmdt_model.transmat_[np.sum(hmdt_model.transmat_, axis=1) == 0] = 1 / n_states
    _, state_seq = hmdt_model.decode(O, F)
    td = train_data.copy()
    td["Module"] = state_seq
    X = td[["Inflow", "Storage", "PDSI", "DOY"]].values
    y = td["Module"].values
    params = {
        "max_depth": [4, 5, 6, 8],
        "min_samples_split": [5, 10, 15],
        "min_samples_leaf": [5, 10, 15],
    }
    clf = GridSearchCV(DecisionTreeClassifier(), params, cv=3, n_jobs=-1)
    clf.fit(X, y)
    return clf.best_estimator_


def simulate_hmdt(data, hmdt_model, ct_model):
    pre_R, pre_S = [], []
    stor = float(data["Storage"].iloc[0])
    for _, row in data.iterrows():
        inflow = float(row["Inflow"])
        pdsi = float(row["PDSI"])
        doy = float(row["DOY"])
        if ct_model is not None:
            mod_id = int(ct_model.predict([[inflow, stor, pdsi, doy]])[0])
        else:
            mod_id = 0
        rel = max(float(hmdt_model.model[mod_id].predict([[inflow, stor]])[0]), 0.0)
        stor = max(min(stor + inflow - rel, 1.0), 0.0)
        pre_R.append(rel)
        pre_S.append(stor)
    return np.array(pre_R), np.array(pre_S)


def simulate_tree(tree, data):
    pre_R, pre_S = [], []
    stor = float(data["Storage"].iloc[0])
    for _, row in data.iterrows():
        rel = max(float(tree.predict([[row["Inflow"], stor]])[0]), 0.0)
        stor = max(min(stor + row["Inflow"] - rel, 1.0), 0.0)
        pre_R.append(rel)
        pre_S.append(stor)
    return np.array(pre_R), np.array(pre_S)


def compute_metrics(obs_r, sim_r, obs_s, sim_s):
    def nse(o, s):
        d = np.sum((o - o.mean()) ** 2)
        return float(1 - np.sum((o - s) ** 2) / d) if d > 0 else np.nan

    def rmse(o, s):
        return float(np.sqrt(np.mean((o - s) ** 2)))

    def kge(o, s):
        if np.std(o) == 0 or np.std(s) == 0:
            return np.nan
        r = np.corrcoef(o, s)[0, 1]
        a = np.std(s) / np.std(o)
        b = np.mean(s) / np.mean(o) if np.mean(o) != 0 else np.nan
        return float(1 - np.sqrt((r - 1) ** 2 + (a - 1) ** 2 + (b - 1) ** 2))

    def pbias(o, s):
        return float(np.sum(o - s) / np.sum(o) * 100) if np.sum(o) != 0 else np.nan

    return {
        "NSE_out": nse(obs_r, sim_r),
        "RMSE_out": rmse(obs_r, sim_r),
        "KGE_out": kge(obs_r, sim_r),
        "PBIAS_out": pbias(obs_r, sim_r),
        "NRMSE_out": nrmse_sd(obs_r, sim_r),
        "NSE_sto": nse(obs_s, sim_s),
        "RMSE_sto": rmse(obs_s, sim_s),
        "KGE_sto": kge(obs_s, sim_s),
        "PBIAS_sto": pbias(obs_s, sim_s),
        "NRMSE_sto": nrmse_sd(obs_s, sim_s),
        "combined": 0.5 * nrmse_sd(obs_r, sim_r) + 0.5 * nrmse_sd(obs_s, sim_s),
    }


def select_best_hmdt(train_d, val_d):
    train_d = train_d.dropna(subset=["Inflow", "Storage", "Storage_obs", "Release", "PDSI", "DOY"]).reset_index(drop=True)
    val_d = val_d.dropna(subset=["Inflow", "Storage", "Storage_obs", "Release", "PDSI", "DOY"]).reset_index(drop=True)
    best_score = np.inf
    best_hmdt = None
    best_ct = None
    best_K = 1
    best_ITS = ITS_VALUES[0]

    for its in ITS_VALUES:
        print(f"    ITS={its}...", end=" ", flush=True)
        models = train_hmdt(train_d, its)
        print(f"{len(models)} models", flush=True)

        for actual_K, model in models:
            ct = train_cart(train_d, model, actual_K) if actual_K > 1 else None
            sim_r, sim_s = simulate_hmdt(val_d, model, ct)
            score = 0.5 * nrmse_sd(val_d["Release"].values, sim_r) + 0.5 * nrmse_sd(val_d["Storage_obs"].values, sim_s)
            if score < best_score:
                best_score = score
                best_hmdt = model
                best_ct = ct
                best_K = actual_K
                best_ITS = its

    return best_hmdt, best_ct, best_K, best_ITS, best_score


def fallback_regression_tree(train_d, val_d):
    train_d = train_d.dropna(subset=["Inflow", "Storage", "Storage_obs", "Release"]).reset_index(drop=True)
    val_d = val_d.dropna(subset=["Inflow", "Storage", "Storage_obs", "Release"]).reset_index(drop=True)
    best_score = np.inf
    best_tree = None
    for depth in [3, 4, 5, 6, 8]:
        for leaf in [10, 20, 30, 50]:
            t = DecisionTreeRegressor(max_depth=depth, min_samples_leaf=leaf, random_state=41)
            t.fit(train_d[["Inflow", "Storage"]].values, train_d["Release"].values)
            sim_r, sim_s = simulate_tree(t, val_d)
            score = 0.5 * nrmse_sd(val_d["Release"].values, sim_r) + 0.5 * nrmse_sd(val_d["Storage_obs"].values, sim_s)
            if score < best_score:
                best_score = score
                best_tree = t
    return best_tree, best_score


def tree_to_rules(tree_obj, smax_acft, feature_names, is_classifier):
    tree_ = tree_obj.tree_
    feats = [feature_names[i] if i != _tree.TREE_UNDEFINED else "?" for i in tree_.feature]
    paths = []

    def recurse(node, path):
        if tree_.feature[node] != _tree.TREE_UNDEFINED:
            nm = feats[node]
            thr = tree_.threshold[node]
            s = (f"{thr * smax_acft:.1f}" if nm in ["Inflow", "Storage"] else str(int(thr)) if nm == "DOY" else f"{thr:.2f}")
            recurse(tree_.children_left[node], path + [f"({nm} <= {s})"])
            recurse(tree_.children_right[node], path + [f"({nm} > {s})"])
        else:
            paths.append(path + [(tree_.value[node], tree_.n_node_samples[node])])

    recurse(0, [])
    rules = []
    for path in paths:
        cond = " and ".join(c for c in path[:-1])
        vals = path[-1][0][0]
        rule = f"if {cond} then "
        rule += (f"module: {int(np.argmax(vals))}" if is_classifier else f"Release: {float(vals[0]) * smax_acft:.2f}")
        rules.append(rule)
    return rules


def export_params(grand_id, hmdt_model, ct_model, smax_acft, is_fallback=False, fallback_tree=None):
    NINF = -1e15
    PINF = 1e15
    module_rows = []
    cart_rows = []

    def parse_leaves(tree_obj, smax, feat_names, is_classifier):
        tree_ = tree_obj.tree_
        feats = [feat_names[i] if i != _tree.TREE_UNDEFINED else None for i in tree_.feature]
        rows = []

        def recurse(node, bounds):
            if tree_.feature[node] != _tree.TREE_UNDEFINED:
                nm = feats[node]
                thr = tree_.threshold[node]
                if nm in ["Inflow", "Storage"]:
                    thr = thr * smax
                b_l = {**bounds, nm + "_hi": min(bounds.get(nm + "_hi", PINF), thr)}
                b_r = {**bounds, nm + "_lo": max(bounds.get(nm + "_lo", NINF), thr)}
                recurse(tree_.children_left[node], b_l)
                recurse(tree_.children_right[node], b_r)
            else:
                vals = tree_.value[node][0]
                row = dict(bounds)
                if is_classifier:
                    row["module_id"] = int(np.argmax(vals))
                else:
                    row["release_acft"] = float(vals[0]) * smax
                rows.append(row)

        recurse(0, {})
        return rows

    if is_fallback:
        for li, leaf in enumerate(parse_leaves(fallback_tree, smax_acft, ["Inflow", "Storage"], False)):
            module_rows.append({
                "GRAND_ID": grand_id,
                "module_id": 0,
                "node_id": li,
                "smax_acft": smax_acft,
                "Inflow_lo": leaf.get("Inflow_lo", NINF),
                "Inflow_hi": leaf.get("Inflow_hi", PINF),
                "Storage_lo": leaf.get("Storage_lo", NINF),
                "Storage_hi": leaf.get("Storage_hi", PINF),
                "release_acft": leaf["release_acft"],
            })
    else:
        for mod_idx, rt in enumerate(hmdt_model.model):
            for li, leaf in enumerate(parse_leaves(rt, smax_acft, ["Inflow", "Storage"], False)):
                module_rows.append({
                    "GRAND_ID": grand_id,
                    "module_id": mod_idx,
                    "node_id": li,
                    "smax_acft": smax_acft,
                    "Inflow_lo": leaf.get("Inflow_lo", NINF),
                    "Inflow_hi": leaf.get("Inflow_hi", PINF),
                    "Storage_lo": leaf.get("Storage_lo", NINF),
                    "Storage_hi": leaf.get("Storage_hi", PINF),
                    "release_acft": leaf["release_acft"],
                })

    if ct_model is not None and not is_fallback:
        for li, leaf in enumerate(parse_leaves(ct_model, smax_acft, ["Inflow", "Storage", "PDSI", "DOY"], True)):
            cart_rows.append({
                "GRAND_ID": grand_id,
                "node_id": li,
                "smax_acft": smax_acft,
                "Inflow_lo": leaf.get("Inflow_lo", NINF),
                "Inflow_hi": leaf.get("Inflow_hi", PINF),
                "Storage_lo": leaf.get("Storage_lo", NINF),
                "Storage_hi": leaf.get("Storage_hi", PINF),
                "PDSI_lo": leaf.get("PDSI_lo", NINF),
                "PDSI_hi": leaf.get("PDSI_hi", PINF),
                "DOY_lo": leaf.get("DOY_lo", NINF),
                "DOY_hi": leaf.get("DOY_hi", PINF),
                "module_id": leaf["module_id"],
            })
    else:
        cart_rows.append({
            "GRAND_ID": grand_id,
            "node_id": 0,
            "smax_acft": smax_acft,
            "Inflow_lo": NINF,
            "Inflow_hi": PINF,
            "Storage_lo": NINF,
            "Storage_hi": PINF,
            "PDSI_lo": NINF,
            "PDSI_hi": PINF,
            "DOY_lo": NINF,
            "DOY_hi": PINF,
            "module_id": 0,
        })

    pd.DataFrame(module_rows).to_csv(RESULTS_DIR / f"module_params_{grand_id}.csv", index=False)
    pd.DataFrame(cart_rows).to_csv(RESULTS_DIR / f"cart_params_{grand_id}.csv", index=False)


if __name__ == "__main__":
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    (RULES_DIR / "modules").mkdir(parents=True, exist_ok=True)
    (RULES_DIR / "module_conditions").mkdir(parents=True, exist_ok=True)

    meta_global = pd.read_csv(TRAINING_DIR / "reservoir_metadata.csv")

    if TEST_IDS is not None:
        files = [TRAINING_DIR / f"{gid}.csv" for gid in TEST_IDS if (TRAINING_DIR / f"{gid}.csv").exists()]
        print(f"Test mode: {[f.stem for f in files]}")
    else:
        files = sorted([f for f in TRAINING_DIR.glob("*.csv") if f.stem.isdigit()])
        print(f"Full mode: {len(files)} reservoirs")

    for i, f in enumerate(files):
        grand_id = int(f.stem)
        print(f"\n[{i + 1}/{len(files)}] GRAND_ID {grand_id}", flush=True)

        smax_row = meta_global[meta_global["GRAND_ID"] == grand_id]
        if len(smax_row) == 0:
            print("  SKIP: not in reservoir_metadata")
            continue
        smax_acft = float(smax_row["STORAGE_CAP"].values[0])

        if (RESULTS_DIR / f"metadata_{grand_id}.csv").exists():
            print("  Already done, skipping")
            continue

        data = load_training_file(f)
        train_d = data[data["Year"] <= CALIB_END].copy()
        val_d = data[(data["Year"] > CALIB_END) & (data["Year"] <= VAL_END)].copy()
        test_d = data[(data["Year"] > CALIB_END) & (data["Year"] <= TEST_END)].copy()
        if len(test_d) < 100:
            test_d = val_d.copy()

        print(f"  train={len(train_d)} val={len(val_d)} test={len(test_d)}")

        is_fallback = False
        ft = None

        if HMM_OK:
            try:
                print("  Training HMDT+CART...", flush=True)
                best_hmdt, best_ct, best_K, best_ITS, val_score = select_best_hmdt(train_d, val_d)
                if best_hmdt is None:
                    raise ValueError("HMDT returned None")
                sim_r, sim_s = simulate_hmdt(test_d, best_hmdt, best_ct)
                m = compute_metrics(test_d["Release"].values, sim_r, test_d["Storage_obs"].values, sim_s)
                method = "HMDT+CART"
                print(f"  K={best_K} ITS={best_ITS} NSE_out={m['NSE_out']:.3f} NSE_sto={m['NSE_sto']:.3f}")
                features = ["Inflow", "Storage", "PDSI", "DOY"]
                for idx, rt in enumerate(best_hmdt.model):
                    rules = tree_to_rules(rt, smax_acft, features, False)
                    (RULES_DIR / "modules" / f"{grand_id}_{idx}.txt").write_text("\n".join(rules) + "\n")
                if best_ct is not None:
                    rules = tree_to_rules(best_ct, smax_acft, features, True)
                    (RULES_DIR / "module_conditions" / f"{grand_id}.txt").write_text("\n".join(rules) + "\n")
                export_params(grand_id, best_hmdt, best_ct, smax_acft)
            except Exception as e:
                print(f"  HMDT failed ({e}), using fallback", flush=True)
                is_fallback = True
        else:
            is_fallback = True

        if is_fallback:
            print("  Training regression tree fallback...", flush=True)
            ft, val_score = fallback_regression_tree(train_d, val_d)
            sim_r, sim_s = simulate_tree(ft, test_d)
            m = compute_metrics(test_d["Release"].values, sim_r, test_d["Storage_obs"].values, sim_s)
            best_K = 1
            best_ITS = 0
            method = "RegressionTree"
            print(f"  NSE_out={m['NSE_out']:.3f} NSE_sto={m['NSE_sto']:.3f}")
            export_params(grand_id, None, None, smax_acft, is_fallback=True, fallback_tree=ft)

        meta_row = {
            "GRAND_ID": grand_id,
            "K": best_K,
            "ITS": best_ITS,
            "smax_acft": smax_acft,
            "method": method,
            "n_train": len(train_d),
            "n_val": len(val_d),
            "n_test": len(test_d),
            "val_combined": val_score,
            "test_NSE_out": m["NSE_out"],
            "test_KGE_out": m["KGE_out"],
            "test_RMSE_out": m["RMSE_out"],
            "test_PBIAS_out": m["PBIAS_out"],
            "test_NRMSE_out": m["NRMSE_out"],
            "test_NSE_sto": m["NSE_sto"],
            "test_KGE_sto": m["KGE_sto"],
            "test_RMSE_sto": m["RMSE_sto"],
            "test_PBIAS_sto": m["PBIAS_sto"],
            "test_NRMSE_sto": m["NRMSE_sto"],
        }
        pd.DataFrame([meta_row]).to_csv(RESULTS_DIR / f"metadata_{grand_id}.csv", index=False)

    print("\nCollecting results...")
    mod_files = sorted(RESULTS_DIR.glob("module_params_*.csv"))
    cart_files = sorted(RESULTS_DIR.glob("cart_params_*.csv"))
    meta_files = sorted(RESULTS_DIR.glob("metadata_*.csv"))

    if not mod_files or not cart_files or not meta_files:
        raise RuntimeError("No result files were generated. Check the training data and model fitting step.")

    pd.concat([pd.read_csv(f) for f in mod_files], ignore_index=True).to_csv(OUT_DIR / "gdrom_module_params.csv", index=False)
    pd.concat([pd.read_csv(f) for f in cart_files], ignore_index=True).to_csv(OUT_DIR / "gdrom_cart_params.csv", index=False)
    meta_all = pd.concat([pd.read_csv(f) for f in meta_files], ignore_index=True)
    meta_all.to_csv(OUT_DIR / "gdrom_metadata.csv", index=False)

    print(f"\nDone: {len(mod_files)} reservoirs")
    print(f"Methods:\n{meta_all['method'].value_counts().to_string()}")
    print(f"NSE_out: mean={meta_all.test_NSE_out.mean():.3f}  median={meta_all.test_NSE_out.median():.3f}")
    print(f"NSE_sto: mean={meta_all.test_NSE_sto.mean():.3f}  median={meta_all.test_NSE_sto.median():.3f}")
    print(f"Outputs: {OUT_DIR}")
