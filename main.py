#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Prova de conceito de desgaste/RUL para bucha oscilante de pá Kaplan.

Modelo-alvo: bucha de deslizamento (plain bushing) no pivô/trunnion da pá do
rotor, sujeita a movimento angular oscilatório. O cálculo de desgaste usa a
forma de Archard dV = K * W * s / H e profundidade = volume/(D*L).

IMPORTANTE: o dataset Kaplan não fornece material, dureza, geometria, carga
normal do mancal, ângulo/movimentos da pá nem desgaste/falha. A diferença de
pressões é apenas um proxy operacional. Portanto, os valores físicos informados
na linha de comando são hipóteses de cenário, não especificações da máquina.
"""
from __future__ import annotations

import argparse
import gzip
import math
import sys
import urllib.request
import zipfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

DATA_URL = "https://api.researchdata.se/dataset/2025-35/1/file/zip"
DEFAULT_ZIP = Path("kaplan_dataset.zip")
REQUIRED = [
    "Time", "p_closeAvg", "p_closeMin", "p_closeMax",
    "p_openAvg", "p_openMin", "p_openMax", "Glycerol",
]


def download_dataset(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Baixando dataset Kaplan oficial para {path} (~116 MB)...")
    urllib.request.urlretrieve(DATA_URL, path)


def find_csv_member(zf: zipfile.ZipFile, turbine: str) -> str:
    expected = f"{turbine.lower()}.csv"
    matches = [name for name in zf.namelist() if Path(name).name.lower() == expected]
    if not matches:
        raise FileNotFoundError(f"Não encontrei {expected} dentro do ZIP do dataset.")
    matches.sort(key=lambda name: ("/data/" not in name.lower(), len(name)))
    return matches[0]


def read_profile(zf: zipfile.ZipFile, turbine: str, stride: int) -> tuple[pd.DataFrame, dict[str, Any]]:
    member = find_csv_member(zf, turbine)
    with zf.open(member) as stream:
        df = pd.read_csv(stream)

    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        raise ValueError(f"{member}: colunas requeridas ausentes: {missing}")

    for col in REQUIRED:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    for col in ["t_oilAvg", "p_sysAvg"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=REQUIRED).sort_values("Time")
    df = df.drop_duplicates("Time", keep="last").reset_index(drop=True)
    if len(df) < 3:
        raise ValueError(f"Poucas observações válidas em {member}.")

    sampled = df.iloc[::stride].copy()
    if sampled.index[-1] != df.index[-1]:
        sampled = pd.concat([sampled, df.iloc[[-1]]], ignore_index=False)
    sampled = sampled.reset_index(drop=True)

    t_days = sampled["Time"].to_numpy(dtype=float)
    dt = np.diff(t_days) * 24.0
    positive = dt[np.isfinite(dt) & (dt > 0)]
    if len(positive) == 0:
        raise ValueError(f"Não foi possível calcular intervalo temporal em {member}.")
    median_dt = float(np.median(positive))
    dt_hours = np.r_[dt, median_dt]
    dt_hours = np.where(np.isfinite(dt_hours) & (dt_hours > 0), dt_hours, median_dt)

    p_close = sampled["p_closeAvg"].to_numpy(float)
    p_open = sampled["p_openAvg"].to_numpy(float)
    delta_p = np.abs(p_close - p_open)
    q95_delta_p = float(np.nanquantile(delta_p, 0.95))
    if not np.isfinite(q95_delta_p) or q95_delta_p <= 1e-12:
        load_ratio = np.ones(len(sampled), dtype=float)
    else:
        # Normalização necessária: unidades da pressão foram omitidas no dataset.
        load_ratio = np.clip(delta_p / q95_delta_p, 0.25, 2.0)

    close_span = np.maximum(
        sampled["p_closeMax"].to_numpy(float) - sampled["p_closeMin"].to_numpy(float), 0.0
    )
    open_span = np.maximum(
        sampled["p_openMax"].to_numpy(float) - sampled["p_openMin"].to_numpy(float), 0.0
    )
    span = np.sqrt(close_span**2 + open_span**2)
    q95_span = float(np.nanquantile(span, 0.95))
    if not np.isfinite(q95_span) or q95_span <= 1e-12:
        activity_proxy = np.ones(len(sampled), dtype=float)
    else:
        activity_proxy = 0.25 + 0.75 * np.clip(span / q95_span, 0.0, 2.0)

    keep = [c for c in [
        "Time", "p_closeAvg", "p_closeMin", "p_closeMax",
        "p_openAvg", "p_openMin", "p_openMax", "Glycerol", "t_oilAvg", "p_sysAvg",
    ] if c in sampled.columns]
    profile = sampled[keep].copy().rename(columns={"Time": "Time_dataset_days"})
    profile["abs_pressure_delta_dataset_units"] = delta_p
    profile["pressure_load_ratio_proxy"] = load_ratio
    profile["activity_proxy_from_pressure_span"] = activity_proxy
    profile["dt_hours"] = dt_hours

    metadata = {
        "source_member": member,
        "valid_input_rows": int(len(df)),
        "profile_rows_after_stride": int(len(profile)),
        "stride": int(stride),
        "dataset_span_days": float(t_days[-1] - t_days[0]),
        "median_sample_interval_min": float(median_dt * 60.0),
        "oil_temperature_present_but_not_used_for_wear": bool("t_oilAvg" in profile.columns),
        "system_pressure_present": bool("p_sysAvg" in profile.columns),
        "q95_pressure_delta_dataset_units": q95_delta_p,
    }
    return profile, metadata


def simulate_scenario(
    profile: pd.DataFrame,
    turbine: str,
    scenario_id: int,
    args: argparse.Namespace,
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    diameter_m = args.bearing_diameter_mm / 1000.0
    length_m = args.bearing_length_mm / 1000.0
    hardness_pa = args.hardness_mpa * 1e6
    limit_m = args.wear_limit_um * 1e-6
    theta_rad = math.radians(args.theta_amplitude_deg)
    radius_m = diameter_m / 2.0
    projected_area_m2 = diameter_m * length_m

    dt = profile["dt_hours"].to_numpy(float)
    load_proxy_n = args.load_ref_n * profile["pressure_load_ratio_proxy"].to_numpy(float)
    activity = profile["activity_proxy_from_pressure_span"].to_numpy(float)
    # Ciclo completo: de -theta a +theta e retorno; distância angular total 4*theta.
    cycles_assumed = args.cycles_per_hour * activity * dt
    distance_m = cycles_assumed * (4.0 * radius_m * theta_rad)

    # V/K = W*s/H; h/K = (V/K)/(D*L).
    depth_per_k_m = load_proxy_n * distance_m / (hardness_pa * projected_area_m2)
    has_glycerol = profile["Glycerol"].to_numpy(float) >= 0.5
    glycerol_factor = np.where(has_glycerol, args.glycerol_factor, 1.0)
    depth_per_k_m = depth_per_k_m * glycerol_factor

    if args.k_archard is not None:
        k_reference = args.k_archard
        calibration = "K_archard fornecido pelo usuário"
    else:
        # Calibração de demonstração: escolhe K para atingir o limite em uma vida
        # nominal assumida, sob a taxa mediana do próprio sinal da turbina.
        median_rate_per_k_m_per_hour = float(np.median(depth_per_k_m / dt))
        if median_rate_per_k_m_per_hour <= 0 or not np.isfinite(median_rate_per_k_m_per_hour):
            raise ValueError("Taxa-base inválida; verifique parâmetros e entradas.")
        k_reference = limit_m / (median_rate_per_k_m_per_hour * args.nominal_life_days * 24.0)
        calibration = "K calibrado para vida nominal hipotética; não é medição experimental"

    k_multiplier = float(rng.lognormal(mean=0.0, sigma=args.k_sigma))
    increments_m = k_reference * k_multiplier * depth_per_k_m

    n = len(profile)
    max_steps = n * args.max_cycles
    repeated_increment = np.tile(increments_m, args.max_cycles)
    repeated_dt = np.tile(dt, args.max_cycles)
    cumulative = np.cumsum(repeated_increment)
    crossing_count = int(np.searchsorted(cumulative, limit_m, side="left") + 1)
    reaches_limit = crossing_count <= max_steps
    n_rows = min(crossing_count, max_steps)

    selected_inc = repeated_increment[:n_rows]
    selected_dt = repeated_dt[:n_rows]
    wear_before_m = np.r_[0.0, np.cumsum(selected_inc[:-1])]
    elapsed_before_h = np.r_[0.0, np.cumsum(selected_dt[:-1])]
    source_idx = np.arange(n_rows) % n

    trajectory = profile.iloc[source_idx].reset_index(drop=True).copy()
    trajectory.insert(0, "turbine", turbine)
    trajectory.insert(1, "scenario_id", scenario_id)
    trajectory["sim_elapsed_hours"] = elapsed_before_h
    trajectory["normal_load_proxy_N"] = np.tile(load_proxy_n, args.max_cycles)[:n_rows]
    trajectory["oscillation_cycles_assumed_interval"] = np.tile(cycles_assumed, args.max_cycles)[:n_rows]
    trajectory["sliding_distance_assumed_m_interval"] = np.tile(distance_m, args.max_cycles)[:n_rows]
    trajectory["wear_depth_sim_um"] = wear_before_m * 1e6
    trajectory["wear_limit_assumed_um"] = args.wear_limit_um
    trajectory["K_archard_used"] = k_reference * k_multiplier
    trajectory["glycerol_wear_factor_assumed"] = np.tile(glycerol_factor, args.max_cycles)[:n_rows]
    trajectory["is_failure_event_sim"] = 0
    trajectory["synthetic_data_flag"] = True

    if reaches_limit:
        last = n_rows - 1
        remaining_depth = max(limit_m - wear_before_m[last], 0.0)
        fraction = float(np.clip(remaining_depth / selected_inc[last], 0.0, 1.0))
        failure_h = float(elapsed_before_h[last] + selected_dt[last] * fraction)
        trajectory["RUL_sim_hours"] = np.maximum(failure_h - elapsed_before_h, 0.0)
        event = trajectory.iloc[[-1]].copy()
        event["sim_elapsed_hours"] = failure_h
        event["wear_depth_sim_um"] = args.wear_limit_um
        event["RUL_sim_hours"] = 0.0
        event["is_failure_event_sim"] = 1
        event["oscillation_cycles_assumed_interval"] = 0.0
        event["sliding_distance_assumed_m_interval"] = 0.0
        trajectory = pd.concat([trajectory, event], ignore_index=True)
        censored = False
    else:
        failure_h = math.nan
        trajectory["RUL_sim_hours"] = np.nan
        censored = True

    summary = {
        "turbine": turbine,
        "scenario_id": scenario_id,
        "rows_output": int(len(trajectory)),
        "simulated_failure_time_hours": failure_h,
        "simulated_failure_time_days": failure_h / 24.0 if np.isfinite(failure_h) else math.nan,
        "right_censored_at_max_cycles": bool(censored),
        "K_archard_reference": k_reference,
        "K_multiplier_sampled": k_multiplier,
        "calibration_method": calibration,
        "nominal_life_days_assumed": args.nominal_life_days,
        "load_reference_at_q95_N_assumed": args.load_ref_n,
        "bearing_diameter_mm_assumed": args.bearing_diameter_mm,
        "bearing_length_mm_assumed": args.bearing_length_mm,
        "hardness_MPa_assumed": args.hardness_mpa,
        "theta_amplitude_deg_assumed": args.theta_amplitude_deg,
        "cycles_per_hour_assumed": args.cycles_per_hour,
        "wear_limit_um_assumed": args.wear_limit_um,
        "glycerol_factor_assumed": args.glycerol_factor,
        "rul_label_synthetic": True,
    }
    return trajectory, summary


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Simula desgaste/RUL de uma bucha oscilante de Kaplan por um modelo Archard simplificado."
    )
    p.add_argument("--zip", dest="zip_path", type=Path, default=DEFAULT_ZIP,
                   help="ZIP oficial do dataset; se não existir, será baixado.")
    p.add_argument("--out", type=Path, default=Path("saidas_bucha_kaplan"), help="Pasta de saída.")
    p.add_argument("--turbines", nargs="+", default=["Turbine_1", "Turbine_2"])
    p.add_argument("--stride", type=int, default=15, help="1 = minuto a minuto; 15 = aproximadamente a cada 15 min.")
    p.add_argument("--scenarios", type=int, default=3)
    p.add_argument("--max-cycles", type=int, default=5, help="Máximo de repetições da série operacional.")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--k-sigma", type=float, default=0.35, help="Dispersão lognormal entre coeficientes K dos cenários.")

    cal = p.add_mutually_exclusive_group(required=True)
    cal.add_argument("--nominal-life-days", type=float,
                     help="Escala hipotética para calibrar K quando não há K de ensaio.")
    cal.add_argument("--k-archard", type=float,
                     help="Coeficiente de desgaste Archard informado pelo usuário/ensaio.")

    # Parâmetros físicos de cenário. Defaults são exemplos explícitos, não dados Kaplan.
    p.add_argument("--load-ref-n", type=float, default=10000.0,
                   help="Carga normal hipotética em N associada ao q95 de |p_close-p_open|.")
    p.add_argument("--bearing-diameter-mm", type=float, default=150.0,
                   help="Diâmetro da bucha; exemplo ilustrativo, substituir por desenho real.")
    p.add_argument("--bearing-length-mm", type=float, default=100.0,
                   help="Comprimento axial da bucha; exemplo ilustrativo.")
    p.add_argument("--hardness-mpa", type=float, default=50.0,
                   help="Dureza efetiva do material em MPa; exemplo, obter ficha/ensaio do material.")
    p.add_argument("--theta-amplitude-deg", type=float, default=3.0,
                   help="Amplitude angular de meia-oscilação em graus; exemplo hipotético.")
    p.add_argument("--cycles-per-hour", type=float, default=2.0,
                   help="Ciclos completos de ida-volta por hora; exemplo hipotético, não medido no dataset.")
    p.add_argument("--wear-limit-um", type=float, default=500.0,
                   help="Limite de desgaste radial assumido em micrômetros; substituir pelo critério técnico.")
    p.add_argument("--glycerol-factor", type=float, default=1.0,
                   help="Multiplicador hipotético de K quando Glycerol=1; 1.0 é neutro.")
    return p


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    if args.stride < 1 or args.scenarios < 1 or args.max_cycles < 1:
        parser.error("stride, scenarios e max-cycles precisam ser >= 1.")
    positive = [args.load_ref_n, args.bearing_diameter_mm, args.bearing_length_mm,
                args.hardness_mpa, args.theta_amplitude_deg, args.cycles_per_hour, args.wear_limit_um]
    if any(v <= 0 for v in positive):
        parser.error("Carga, dimensões, dureza, movimento, ciclos e limite de desgaste devem ser > 0.")
    if args.k_sigma < 0:
        parser.error("k-sigma deve ser >= 0.")
    if args.k_archard is not None and args.k_archard <= 0:
        parser.error("k-archard deve ser > 0.")
    if args.nominal_life_days is not None and args.nominal_life_days <= 0:
        parser.error("nominal-life-days deve ser > 0.")
    if not (0 < args.glycerol_factor <= 2):
        parser.error("glycerol-factor deve estar entre 0 e 2.")

    if not args.zip_path.exists():
        try:
            download_dataset(args.zip_path)
        except Exception as exc:
            print(f"Falha ao baixar dataset: {exc}", file=sys.stderr)
            return 2

    args.out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    summaries: list[dict[str, Any]] = []
    print("\nAVISO DE MODELAGEM")
    print("- A bucha e os parâmetros físicos abaixo são um cenário hipotético, não especificações das turbinas do dataset.")
    print("- A pressão diferencial é normalizada e mapeada para uma carga assumida; não é carga do mancal medida.")
    print("- Glycerol=1 usa fator neutro por padrão; qualquer fator diferente de 1 é hipótese.")

    try:
        with zipfile.ZipFile(args.zip_path) as zf:
            for turbine in args.turbines:
                print(f"\nLendo {turbine}...")
                profile, metadata = read_profile(zf, turbine, args.stride)
                print(f"  {metadata['valid_input_rows']:,} pontos válidos; {len(profile):,} usados após stride.")
                output_file = args.out / f"rul_sintetico_bucha_{turbine}.csv.gz"
                if output_file.exists():
                    output_file.unlink()
                with gzip.open(output_file, "wt", encoding="utf-8", newline="") as gz:
                    for sid in range(1, args.scenarios + 1):
                        trajectory, summary = simulate_scenario(profile, turbine, sid, args, rng)
                        trajectory.to_csv(gz, index=False, header=(sid == 1))
                        summary.update(metadata)
                        summaries.append(summary)
                        t_days = summary["simulated_failure_time_days"]
                        t_text = f"{t_days:.1f} d" if np.isfinite(t_days) else "censurado"
                        print(f"  cenário {sid}: {len(trajectory):,} linhas; evento simulado={t_text}")
                print(f"  Saída: {output_file}")

        summary_file = args.out / "resumo_parametros_e_cenarios.csv"
        pd.DataFrame(summaries).to_csv(summary_file, index=False)
        print(f"\nResumo: {summary_file}")
        print("RUL_sim_hours é rótulo sintético calculado pelo modelo; não é RUL validado em turbina real.")
    except Exception as exc:
        print(f"Erro na simulação: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
