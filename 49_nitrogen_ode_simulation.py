"""
49_nitrogen_ode_simulation.py
=============================
Nitrogen mass-balance ODE for the ISS pepper rhizosphere. Community-scale,
literature-parameterised, closed-budget model.

STATE VARIABLES (NH4, NO2, NO3, B_plant as mg-N per litre of rhizosphere
pore water; D dimensionless in [0, 1]):
    NH4       — ammonium pool (fed by Osmocote; consumed by AOB and plant)
    NO2       — nitrite intermediate
    NO3       — nitrate pool (plant-available)
    B_plant   — nitrogen stored in standing plant biomass (harvested)
    D         — per-crop physiological damage state (accumulates with NH4+
                toxicity exposure within a crop cycle; RESETS to 0 at each
                harvest event because harvest = new seedling transplant;
                damage represents crop-level physiology, NOT substrate
                legacy). Scales plant uptake by (1 - D).

RATE EXPRESSIONS
    Osmocote supply:       dNH4/dt += R_SUPPLY
    AOB oxidation:         NIT_MAX_AOB * fA(t) * NH4 / (K_AOB + NH4)
    NOB oxidation:         NIT_MAX_NOB * fN(t) * NO2 / (K_NOB + NO2)
    Plant NO3 uptake:      V_PLANT_NO3 * NO3 / (K_PLANT + NO3) * inh * (1-D)
    Plant NH4 uptake:      V_PLANT_NH4 * NH4 / (K_PLANT_NH4 + NH4) * inh * (1-D)
    Within-crop damage:    dD/dt = K_DAMAGE * max(0, NH4 - NH4_TOXICITY)^2 * (1-D)
    Harvest (closure):     discrete event at day 120 and 240 — 95% of B_plant
                            leaves the system; B_plant residual retained;
                            D resets to 0 (new seedling transplant).
                            Timing reflects the Advanced Plant Habitat (APH)
                            NuMex Española Improved chile pepper cycle length
                            (~120 days per generation; red harvest at day 109,
                            green harvest at day 137 per Khodadad et al. 2026
                            Sci. Rep.; 120 d is used as the single illustrative
                            value).
                            Nitrifier state f_A(t), f_N(t) persists across
                            harvest (microbial community unaffected by crop
                            replacement).

Rates fA(t), fN(t) are community-level fractions of maximal nitrifier
activity; scenario-dependent (healthy = 1, ISS = 0, inoc at T_INOC = logistic
recovery).

LITERATURE-PARAMETERISED VALUES
    NIT_MAX_AOB  = 60  mg-N / L / day  (community potential; Norton & Stark 2011;
                                         Schimel & Bennett 2004)
    K_AOB        = 2   mg-N / L        (low Km; Stein 2019; AOB substrate affinity
                                         is high at low NH4+)
    NIT_MAX_NOB  = 80  mg-N / L / day  (typically > AOB, rarely limiting;
                                         Sundermeyer-Klinger 1984)
    K_NOB        = 0.5 mg-N / L        (Sundermeyer-Klinger 1984)
    R_SUPPLY     = 8   mg-N / L / day  (Osmocote 14-14-14, Veggie pillow
                                         scale, Massa 2017; Shaviv 2001)
    V_PLANT_NO3  = 15  mg-N / L / day  (calibrated so plant uptake exceeds
                                         Osmocote supply, allowing a closed
                                         mass budget with harvest)
    K_PLANT      = 2
    V_PLANT_NH4  = 3   (plants take NO3 preferentially in pepper;
                         ~80:20 NO3:NH4 preference, Britto & Kronzucker 2002)
    K_PLANT_NH4  = 5

    NH4 toxicity threshold = 50  mg-N / L (crop stress; Britto 2002)
    NH4 severe stress      = 100 mg-N / L

    Harvest days           = 120, 240 (APH pepper cycle length ~120 days per
                            generation; Khodadad et al. 2026 Sci. Rep.)
    Harvest fraction       = 0.95 of standing B_plant removed

Outputs (acta_astronautica/simulation/):
    ode_trajectories.csv      — full time course, all state variables
    ode_scenario_summary.csv  — day-to-toxicity, T_recovery, end-state
    ode_mass_balance.csv      — cumulative in/out per scenario (closure check)
    ode_trajectories.pdf/.png — 4-panel figure
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config

OUT_DIR = config.AA_DIR / "simulation"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# -------------------------------------------------------------------
# Parameters (documented in header; mg-N/L and days)
# -------------------------------------------------------------------
NIT_MAX_AOB = 60.0
K_AOB       = 2.0
NIT_MAX_NOB = 80.0
K_NOB       = 0.5

R_SUPPLY    = 8.0

V_PLANT_NO3 = 15.0       # calibrated to close mass balance under healthy regime
K_PLANT     = 2.0
V_PLANT_NH4 = 3.0
K_PLANT_NH4 = 5.0

NH4_TOXICITY    = 50.0
NH4_HIGH_STRESS = 100.0

# Hill-type ammonium-toxicity inhibition on plant uptake. Britto & Kronzucker
# 2002 review NH4+ toxicity mechanisms (ion imbalance, proton-pump dissipation,
# photosystem II impairment) but do not publish a single Hill coefficient; we
# use n_Hill = 3 as a moderate-cooperativity default and run a one-way
# sensitivity sweep over n_Hill = {1, 2, 3, 4} (see ode_sensitivity_hill.csv).
NH4_HILL_N = 3.0

# Per-crop physiological damage state (5th state variable, D ∈ [0, 1]):
# accumulates irreversibly *within a crop cycle* with toxic NH4+ exposure
# and scales plant uptake by (1 - D). Reset to 0 at each harvest event
# because harvest = new seedling transplant (see simulate() and header
# docstring). Damage does NOT carry across crop cycles in this model;
# substrate-level legacy damage is a separate, unmodelled effect flagged
# in the Limitations. Calibration: D reaches ~0.5 after ~10 days at
# [NH4] ≈ 100 mg-N/L (severe stress); K_DAMAGE is the main calibration
# parameter controlling the late-inoculation first-crop loss magnitude
# (see ode_sensitivity_kdamage.csv).
K_DAMAGE = 2.0e-5        # per (mg-N/L)^2 per day

DAYS             = 260    # extended past day-240 harvest so post-harvest state
                          # is captured; summary flags indicate day-240 vs
                          # day-260 readouts separately.
HARVEST_DAYS     = [120, 240]   # APH pepper cycle length ~120 d per generation
                                # (Khodadad et al. 2026 Sci. Rep.; red harvest
                                # day 109, green day 137; 120 d used as single
                                # illustrative value between the two harvest
                                # readouts of the same crop).
HARVEST_FRACTION = 0.95


# -------------------------------------------------------------------
# Scenario-specific nitrifier activity
# -------------------------------------------------------------------

def nitrifier_fraction(t, scenario, t_inoc, r_growth=0.25):
    if scenario == "healthy":
        return 1.0, 1.0
    if scenario == "iss":
        return 0.0, 0.0
    if scenario.startswith("inoc_"):
        if t < t_inoc:
            return 0.0, 0.0
        dt = t - t_inoc
        seed = 0.01
        frac = 1.0 * seed / (seed + (1.0 - seed) * np.exp(-r_growth * dt))
        return frac, frac
    raise ValueError(scenario)


# -------------------------------------------------------------------
# ODE right-hand side (no harvest term here — handled piecewise)
# -------------------------------------------------------------------

def rhs(t, y, scenario, t_inoc, hill=None, k_dmg=None):
    NH4, NO2, NO3, B_plant, D = y
    NH4 = max(NH4, 0.0)
    NO2 = max(NO2, 0.0)
    NO3 = max(NO3, 0.0)
    D   = min(max(D, 0.0), 1.0)

    hill = NH4_HILL_N if hill is None else hill
    k_dmg = K_DAMAGE if k_dmg is None else k_dmg

    fa, fb = nitrifier_fraction(t, scenario, t_inoc)

    v_amo = NIT_MAX_AOB * fa * NH4 / (K_AOB + NH4)
    v_nxr = NIT_MAX_NOB * fb * NO2 / (K_NOB + NO2)

    # Reversible Hill-type toxicity inhibition (acute; dissipates when NH4
    # drops).
    inh = 1.0 / (1.0 + (NH4 / NH4_TOXICITY) ** hill)
    # Irreversible damage reduces plant uptake capacity on top of the
    # reversible inhibition.
    cap = max(0.0, 1.0 - D)

    v_plant_no3 = V_PLANT_NO3 * NO3 / (K_PLANT + NO3) * inh * cap
    v_plant_nh4 = V_PLANT_NH4 * NH4 / (K_PLANT_NH4 + NH4) * inh * cap

    excess = max(0.0, NH4 - NH4_TOXICITY)
    dD = k_dmg * excess ** 2 * (1.0 - D)   # saturates at D = 1

    dNH4     = R_SUPPLY - v_amo - v_plant_nh4
    dNO2     = v_amo - v_nxr
    dNO3     = v_nxr - v_plant_no3
    dB_plant = v_plant_no3 + v_plant_nh4
    return [dNH4, dNO2, dNO3, dB_plant, dD]


# -------------------------------------------------------------------
# Piecewise integrator with harvest events
# -------------------------------------------------------------------

def simulate(scenario, t_inoc=None, hill=None, k_dmg=None):
    """Integrate ODE segment-by-segment, applying harvest at HARVEST_DAYS.

    Returns
    -------
    df : pd.DataFrame
        Time course with pre- and post-harvest rows at each harvest day.
        `phase_order` column (0 = pre-harvest, 1 = post-harvest) preserves
        within-day ordering.
    harvest_events : list[dict]
    """
    # y0: [NH4, NO2, NO3, B_plant, D]; D starts at 0 (no accumulated damage)
    y = [5.0, 0.1, 0.1, 0.0, 0.0]
    segments = []
    segment_bounds = sorted(set([0] + HARVEST_DAYS + [DAYS]))
    harvest_events = []
    cumulative_harvest = 0.0

    for seg_start, seg_end in zip(segment_bounds[:-1], segment_bounds[1:]):
        t_eval = np.arange(seg_start, seg_end + 1)
        sol = solve_ivp(rhs, (seg_start, seg_end), y,
                        args=(scenario, t_inoc, hill, k_dmg), t_eval=t_eval,
                        method="LSODA", rtol=1e-8, atol=1e-10)
        if not sol.success:
            raise RuntimeError(f"ODE integration failed at segment "
                               f"[{seg_start},{seg_end}] for scenario {scenario}: "
                               f"{sol.message}")
        if np.any(sol.y < -1e-6):
            print(f"  WARNING: negative state in {scenario} "
                  f"[{seg_start},{seg_end}]: min = {sol.y.min():.2e}")
        sol_y = np.clip(sol.y, a_min=0.0, a_max=None)

        seg_df = pd.DataFrame({
            "day":                sol.t,
            "NH4":                sol_y[0],
            "NO2":                sol_y[1],
            "NO3":                sol_y[2],
            "B_plant":            sol_y[3],
            "D":                  sol_y[4],
            "scenario":           scenario,
            "cumulative_harvest": cumulative_harvest,
            "phase":              "pre",
            "phase_order":        0,
        })
        segments.append(seg_df)

        # Harvest event: remove 95% of standing B_plant, RESET per-crop damage
        # D to 0 (harvest = new seedling transplant; previous crop's tissue
        # damage does not carry over to the next crop). Nitrifier community
        # state (f_A, f_N via t_inoc) persists across harvest.
        y = list(sol_y[:, -1])
        if seg_end in HARVEST_DAYS:
            bp_pre = y[3]
            d_pre  = y[4]
            harvested = HARVEST_FRACTION * bp_pre
            y[3] = bp_pre - harvested
            y[4] = 0.0                  # reset per-crop damage
            cumulative_harvest += harvested
            harvest_events.append({
                "day":              int(seg_end),
                "Bplant_pre":       round(bp_pre, 2),
                "harvested":        round(harvested, 2),
                "Bplant_post":      round(y[3], 2),
                "D_pre_harvest":    round(d_pre, 3),
                "D_post_reset":     0.0,
            })
            post_row = pd.DataFrame({
                "day":                [int(seg_end)],
                "NH4":                [y[0]],
                "NO2":                [y[1]],
                "NO3":                [y[2]],
                "B_plant":            [y[3]],
                "D":                  [y[4]],
                "scenario":           [scenario],
                "cumulative_harvest": [cumulative_harvest],
                "phase":              ["post"],
                "phase_order":        [1],
            })
            segments.append(post_row)

    df = pd.concat(segments, ignore_index=True)
    # Sort by (day, phase_order) so the post-harvest row is *after* the
    # pre-harvest row on the same day (string sort would put "post" before
    # "pre"). Drop only true duplicates (same day + same phase).
    df = df.sort_values(["day", "phase_order"], kind="stable").drop_duplicates(
        subset=["day", "phase_order"]).reset_index(drop=True)
    return df, harvest_events


# -------------------------------------------------------------------
# Metrics
# -------------------------------------------------------------------

def detect_failure(df: pd.DataFrame, harvest_events, t_inoc):
    d_pre = df[df["phase"] == "pre"]
    nh4 = d_pre["NH4"].values
    day = d_pre["day"].values

    day_tox    = int(day[nh4 >= NH4_TOXICITY].min()) \
                 if (nh4 >= NH4_TOXICITY).any()    else -1
    day_severe = int(day[nh4 >= NH4_HIGH_STRESS].min()) \
                 if (nh4 >= NH4_HIGH_STRESS).any() else -1

    day_recovery = -1
    if day_tox > 0:
        after = d_pre[d_pre["day"] > day_tox].copy()
        after_sub = after[after["NH4"] < NH4_TOXICITY]
        if len(after_sub):
            cand = int(after_sub["day"].iloc[0])
            tail = after[after["day"] >= cand]
            if (tail["NH4"] < NH4_TOXICITY).all():
                day_recovery = cand

    harvest_days = [h["day"] for h in harvest_events]
    harvest_amts = [h["harvested"] for h in harvest_events]
    damage_at_harvest = [h["D_pre_harvest"] for h in harvest_events]
    harvest_1 = harvest_amts[0] if len(harvest_amts) > 0 else 0.0
    harvest_2 = harvest_amts[1] if len(harvest_amts) > 1 else 0.0
    total_harvested = sum(harvest_amts)

    # Per-crop peak damage: D resets at each harvest. Crop 1 is the
    # [0, HARVEST_DAYS[0]] window (pre-rows only) and crop 2 is the
    # (HARVEST_DAYS[0], HARVEST_DAYS[1]] window. Crop 2 starts at D=0 thanks
    # to the harvest reset.
    h1, h2 = HARVEST_DAYS[0], HARVEST_DAYS[1]
    d_crop1_peak = float(d_pre[d_pre["day"] <= h1]["D"].max())
    crop2_rows   = d_pre[(d_pre["day"] > h1) & (d_pre["day"] <= h2)]
    d_crop2_peak = float(crop2_rows["D"].max()) if len(crop2_rows) else -1.0

    # Bplant state RIGHT AFTER the second harvest (end of second crop cycle)
    bp_after_2nd = -1.0
    if df["day"].max() >= h2:
        post_h2 = df[(df["day"] == h2) & (df["phase"] == "post")]
        if len(post_h2):
            bp_after_2nd = float(post_h2["B_plant"].iloc[0])

    # Regrowth between the second harvest and end-of-simulation:
    # not a third crop cycle; represents partial regrowth monitoring only.
    final_row = df.iloc[-1]

    return {
        "day_NH4_toxicity_50":         day_tox,
        "day_NH4_severe_100":          day_severe,
        "day_recovery":                day_recovery,
        "T_recovery_days":             (day_recovery - day_tox
                                        if day_recovery > 0 and day_tox > 0 else -1),
        "peak_NH4":                    round(float(nh4.max()), 2),
        "harvest_day_1":               harvest_days[0] if harvest_days else -1,
        "harvest_day_2":               harvest_days[1] if len(harvest_days) > 1 else -1,
        "harvest_1_mgN_L":             round(harvest_1, 2),
        "harvest_2_mgN_L":             round(harvest_2, 2),
        "total_harvested_N_2cycles":   round(total_harvested, 2),
        "Bplant_post_h2":              round(bp_after_2nd, 2),
        "Bplant_final_partial_regrowth": round(float(final_row["B_plant"]), 2),
        "D_pre_harvest_1":             round(damage_at_harvest[0], 3) if len(damage_at_harvest) else "NA",
        "D_pre_harvest_2":             round(damage_at_harvest[1], 3) if len(damage_at_harvest) > 1 else "NA",
        "D_peak_crop1":                round(d_crop1_peak, 3),
        "D_peak_crop2":                round(d_crop2_peak, 3),
        "end_NH4":                     round(float(final_row["NH4"]), 2),
        "end_NO3":                     round(float(final_row["NO3"]), 2),
    }


def check_mass_balance(df: pd.DataFrame, harvest_events):
    t_end = df["day"].iloc[-1]
    cum_in = R_SUPPLY * t_end
    final_row = df.iloc[-1]
    # D is dimensionless damage state, not N; exclude from pool
    current_pool = float(final_row["NH4"] + final_row["NO2"]
                         + final_row["NO3"] + final_row["B_plant"])
    initial_pool = 5.0 + 0.1 + 0.1 + 0.0
    cum_harvest = sum(h["harvested"] for h in harvest_events)
    imbalance = cum_in - (current_pool + cum_harvest - initial_pool)
    return {
        "cum_supply_final":   round(cum_in, 2),
        "current_N_pool":     round(current_pool, 2),
        "cum_harvest_2cycles": round(cum_harvest, 2),
        "n_harvest_events":   len(harvest_events),
        "initial_pool":       round(initial_pool, 2),
        "imbalance_mgN_L":    round(imbalance, 4),
    }


_SENS_SCENARIOS = [
    ("healthy",  None),
    ("iss",      None),
    ("inoc_d5",  5),
    ("inoc_d10", 10),
    ("inoc_d30", 30),
]


def hill_sensitivity_sweep():
    """One-way sensitivity on the Hill coefficient for NH4 toxicity
    inhibition. Runs each focus scenario at n_Hill = 1, 2, 3, 4."""
    rows = []
    for n_h in (1.0, 2.0, 3.0, 4.0):
        for scenario, t_inoc in _SENS_SCENARIOS:
            df, h_events = simulate(scenario, t_inoc=t_inoc, hill=n_h)
            summ = detect_failure(df, h_events, t_inoc)
            rows.append({
                "n_Hill":               n_h,
                "scenario":             scenario,
                "day_NH4_toxicity_50":  summ["day_NH4_toxicity_50"],
                "day_NH4_severe_100":   summ["day_NH4_severe_100"],
                "peak_NH4":             summ["peak_NH4"],
                "total_harvested_N":    summ["total_harvested_N_2cycles"],
                "harvest_1":            summ["harvest_1_mgN_L"],
                "harvest_2":            summ["harvest_2_mgN_L"],
                "D_peak_crop1":         summ["D_peak_crop1"],
                "D_peak_crop2":         summ["D_peak_crop2"],
            })
    return pd.DataFrame(rows)


def kdamage_sensitivity_sweep():
    """One-way sensitivity on the K_DAMAGE calibration constant
    (irreversible-damage accumulation rate). Multiplier grid chosen to span
    an order of magnitude around the baseline (2e-5)."""
    rows = []
    for mult in (0.5, 1.0, 2.0, 5.0):
        k_val = K_DAMAGE * mult
        for scenario, t_inoc in _SENS_SCENARIOS:
            df, h_events = simulate(scenario, t_inoc=t_inoc, k_dmg=k_val)
            summ = detect_failure(df, h_events, t_inoc)
            rows.append({
                "K_DAMAGE_multiplier":  mult,
                "K_DAMAGE_value":       k_val,
                "scenario":             scenario,
                "day_NH4_toxicity_50":  summ["day_NH4_toxicity_50"],
                "day_NH4_severe_100":   summ["day_NH4_severe_100"],
                "peak_NH4":             summ["peak_NH4"],
                "total_harvested_N":    summ["total_harvested_N_2cycles"],
                "harvest_1":            summ["harvest_1_mgN_L"],
                "harvest_2":            summ["harvest_2_mgN_L"],
                "D_peak_crop1":         summ["D_peak_crop1"],
                "D_peak_crop2":         summ["D_peak_crop2"],
            })
    return pd.DataFrame(rows)


# -------------------------------------------------------------------
# Main
# -------------------------------------------------------------------

def main():
    print("=" * 60)
    print("49: Nitrogen mass-balance ODE simulation (closed budget)")
    print("=" * 60)

    scenarios = [
        ("healthy",    None),
        ("iss",        None),
        ("inoc_d0",    0),
        ("inoc_d5",    5),
        ("inoc_d10",   10),
        ("inoc_d20",   20),
        ("inoc_d30",   30),
    ]

    all_traj, summary_rows, balance_rows, harvest_rows = [], [], [], []
    for s, t_inoc in scenarios:
        print(f"\n[{s}] simulating {DAYS} days (t_inoc={t_inoc})...")
        df, harvest_events = simulate(s, t_inoc=t_inoc)
        all_traj.append(df)

        for he in harvest_events:
            harvest_rows.append({"scenario": s, **he})

        summ = detect_failure(df, harvest_events, t_inoc)
        summ["scenario"] = s
        summ["t_inoc"] = t_inoc
        summary_rows.append(summ)
        print(f"  peak NH4 = {summ['peak_NH4']:.1f}; "
              f"tox@{summ['day_NH4_toxicity_50']}; "
              f"severe@{summ['day_NH4_severe_100']}; "
              f"recover@{summ['day_recovery']}; "
              f"Trec={summ['T_recovery_days']} d; "
              f"total_harvest_2cycles={summ['total_harvested_N_2cycles']:.0f} "
              f"(crop1:{summ['harvest_1_mgN_L']:.0f}, crop2:{summ['harvest_2_mgN_L']:.0f}); "
              f"D_peak crop1={summ['D_peak_crop1']:.3f}, crop2={summ['D_peak_crop2']:.3f}")

        bal = check_mass_balance(df, harvest_events)
        bal["scenario"] = s
        balance_rows.append(bal)
        print(f"  mass balance: in {bal['cum_supply_final']:.0f}, pool {bal['current_N_pool']:.1f},"
              f" harvest {bal['cum_harvest_2cycles']:.1f}, imbalance {bal['imbalance_mgN_L']:.4f}")

    traj = pd.concat(all_traj, ignore_index=True)
    traj.to_csv(OUT_DIR / "ode_trajectories.csv", index=False)

    summary = pd.DataFrame(summary_rows)
    summary = summary[["scenario", "t_inoc",
                       "day_NH4_toxicity_50", "day_NH4_severe_100",
                       "day_recovery", "T_recovery_days", "peak_NH4",
                       "harvest_day_1", "harvest_1_mgN_L",
                       "harvest_day_2", "harvest_2_mgN_L",
                       "total_harvested_N_2cycles", "Bplant_post_h2",
                       "Bplant_final_partial_regrowth",
                       "D_peak_crop1", "D_pre_harvest_1",
                       "D_peak_crop2", "D_pre_harvest_2",
                       "end_NH4", "end_NO3"]]
    summary.to_csv(OUT_DIR / "ode_scenario_summary.csv", index=False)
    print("\nSummary:")
    print(summary.to_string(index=False))

    balance = pd.DataFrame(balance_rows)
    balance.to_csv(OUT_DIR / "ode_mass_balance.csv", index=False)

    harvest_df = pd.DataFrame(harvest_rows)
    harvest_df.to_csv(OUT_DIR / "ode_harvest_events.csv", index=False)

    # One-way sensitivity on NH4_HILL_N
    print("\n[sensitivity] sweeping n_Hill = 1..4...")
    sens_hill = hill_sensitivity_sweep()
    sens_hill.to_csv(OUT_DIR / "ode_sensitivity_hill.csv", index=False)
    print(sens_hill.to_string(index=False))

    # One-way sensitivity on K_DAMAGE (main calibration parameter for
    # irreversible damage accumulation; Hill sensitivity alone is not
    # sufficient since early-vs-late inoculation divergence is governed by
    # D accumulation).
    print("\n[sensitivity] sweeping K_DAMAGE multiplier = 0.5, 1, 2, 5...")
    sens_kdmg = kdamage_sensitivity_sweep()
    sens_kdmg.to_csv(OUT_DIR / "ode_sensitivity_kdamage.csv", index=False)
    print(sens_kdmg.to_string(index=False))

    # ──── Plot ────
    focus = ["healthy", "iss", "inoc_d5", "inoc_d30"]
    fig, axes = plt.subplots(4, 1, figsize=(11, 12), sharex=True)
    color_map = {"healthy": "#2ecc71", "iss": "#e74c3c",
                 "inoc_d5": "#2980b9", "inoc_d30": "#8e44ad"}
    label_map = {
        "healthy":  "Healthy (AOB + NOB active)",
        "iss":      "ISS observed (AOB = NOB = 0)",
        "inoc_d5":  "Inoculation @ day 5",
        "inoc_d30": "Inoculation @ day 30",
    }

    for s in focus:
        d = traj[traj["scenario"] == s]
        axes[0].plot(d["day"], d["NH4"],     color=color_map[s], lw=2.2,
                     label=label_map[s])
        axes[1].plot(d["day"], d["NO2"],     color=color_map[s], lw=2.2)
        axes[2].plot(d["day"], d["NO3"],     color=color_map[s], lw=2.2)
        axes[3].plot(d["day"], d["B_plant"], color=color_map[s], lw=2.2)

    axes[0].axhline(NH4_TOXICITY,    color="orange", ls=":", lw=1.3)
    axes[0].axhline(NH4_HIGH_STRESS, color="red",    ls=":", lw=1.3)
    # Harvest markers on B_plant panel
    for hd in HARVEST_DAYS:
        axes[3].axvline(hd, color="black", ls="--", lw=1, alpha=0.55)
        axes[3].text(hd, axes[3].get_ylim()[1] * 0.02, f"harvest day {hd}",
                     rotation=90, fontsize=8, va="bottom", ha="right", color="black")
    for ti, col in [(5, "#2980b9"), (30, "#8e44ad")]:
        axes[0].axvline(ti, color=col, ls="-.", lw=1, alpha=0.55)

    axes[0].set_ylabel(r"NH$_4^+$" + "\n(mg-N/L)")
    axes[1].set_ylabel(r"NO$_2^-$" + "\n(mg-N/L)")
    axes[2].set_ylabel(r"NO$_3^-$" + "\n(mg-N/L)")
    axes[3].set_ylabel(f"Plant biomass N\n(standing pool,\nharvested at d{HARVEST_DAYS[0]}/d{HARVEST_DAYS[1]})")
    axes[3].set_xlabel("Simulation day")
    for ax in axes:
        ax.grid(alpha=0.3)

    # Text annotations for thresholds (so legend doesn't block data)
    x_anno = traj["day"].max() * 0.98
    axes[0].text(x_anno, NH4_TOXICITY * 1.05, "50 mg-N/L (toxicity)",
                 color="#d4780f", fontsize=8, ha="right", va="bottom")
    axes[0].text(x_anno, NH4_HIGH_STRESS * 1.05, "100 mg-N/L (severe)",
                 color="#c0392b", fontsize=8, ha="right", va="bottom")

    # Compact legend outside axis area
    axes[0].legend(loc="center left", bbox_to_anchor=(1.01, 0.5),
                   fontsize=9, frameon=True)

    fig.suptitle("Closed-budget N-mass-balance ODE simulation — "
                 f"R_supply={R_SUPPLY} mg-N/L/day; harvest at days {HARVEST_DAYS}",
                 fontsize=10, fontweight="bold", y=0.995)
    fig.subplots_adjust(left=0.10, right=0.78, top=0.95, bottom=0.05, hspace=0.12)

    for ext in ("pdf", "png"):
        fig.savefig(OUT_DIR / f"ode_trajectories.{ext}",
                    dpi=300 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  Saved: ode_trajectories.pdf/.png")
    print(f"\nDone. Outputs in {OUT_DIR}/")


if __name__ == "__main__":
    main()
