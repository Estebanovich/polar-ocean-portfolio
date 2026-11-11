"""
plot_flujo_boca_bt_bc.py
------------------------
Series de tiempo y tabla-resumen del flujo de energía NETO a través de la boca
de la bahía, separado en barotrópico (p_bt'·u_bt'·H) y baroclínico (Nash 2005).

    F_total(t) [W] = Σ_x F_sec(t, x) · dx     (+ = hacia la bahía, +y)

Lee npz_flux/{scen}_{run}_flujo.npz generados por energy_flux_IWBay.py.

Figuras (una por escenario):
    figs/flujo_boca/flujo_bt_bc_{scen}.png    F_bt(t) y F_bc(t), bahía vs sin-bahía
Tabla en stdout: ⟨F⟩ y max|F| por escenario/run/etapa para cada componente.

Uso:
    python plot_flujo_boca_bt_bc.py
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(_BASE))

from energy_flux_IWBay import SCENARIOS, DX, cargar_flujo, SAVE_DIR

FIGS_DIR = _BASE / "figs" / "flujo_boca"
FIGS_DIR.mkdir(parents=True, exist_ok=True)

RUNS   = ["bay", "nobay"]
ETAPAS = [1, 2]
RUN_STYLE = {
    "bay":   dict(color="#d62728", lw=1.4, label="con bahía"),
    "nobay": dict(color="#1f77b4", lw=1.2, ls="--", label="sin bahía"),
}


def _integrar(Fsec, dx_sec=None):
    """Σ_x F_sec(t,x)·dx(x) → F_total(t) [W]. NaN si alguna columna es NaN.

    dx_sec viene del npz (dxF real de cada columna): la malla es telescópica
    desde jun-2026, así que el DX=5000 constante subestimaba ~1% en la boca.
    """
    if Fsec is None:
        return None
    dx = DX if dx_sec is None else np.asarray(dx_sec)
    F = np.nansum(Fsec * dx, axis=1)
    F[np.any(np.isnan(Fsec), axis=1)] = np.nan
    return F


def _serie_concatenada(scen, run, key_sec):
    """Concatena las dos etapas de un componente (key_sec ∈ {'Fy_sec','Fy_bt_sec'})."""
    t_list, F_list = [], []
    for et in ETAPAS:
        try:
            d = cargar_flujo(scen, run, et, npz_dir=SAVE_DIR)
        except FileNotFoundError:
            continue
        F = _integrar(d.get(key_sec), d.get("dx_sec"))
        if F is None:
            continue
        t = d["t_horas"]
        if t_list and abs(t[0] - t_list[-1][-1]) < 1e-6:
            t = t[1:]; F = F[1:]
        t_list.append(t); F_list.append(F)
    if not t_list:
        return None, None
    return np.concatenate(t_list), np.concatenate(F_list)


def plot_scenario(scen):
    fig, (ax_bt, ax_bc) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    any_data = False

    for run in RUNS:
        t_bt, F_bt = _serie_concatenada(scen, run, "Fy_bt_sec")
        t_bc, F_bc = _serie_concatenada(scen, run, "Fy_sec")
        if F_bt is not None:
            ax_bt.plot(t_bt, F_bt / 1e3, **RUN_STYLE[run]); any_data = True
        if F_bc is not None:
            ax_bc.plot(t_bc, F_bc / 1e3, **RUN_STYLE[run]); any_data = True

    if not any_data:
        plt.close(fig); return

    for ax, ttl in ((ax_bt, "Barotrópico  F_bt = p_bt'·u_bt'·H"),
                    (ax_bc, "Baroclínico  F_bc = ∫p'·v' dz  (Nash 2005)")):
        ax.axhline(0, color="k", lw=0.6)
        ax.axvline(60, color="0.5", ls=":", lw=1.0)   # forzado → libre (~60 h)
        ax.set_ylabel("F a la boca [kW]")
        ax.set_title(ttl, fontsize=10)
        ax.grid(ls=":", lw=0.5)
        ax.legend(fontsize=8)
    ax_bc.set_xlabel("Tiempo (h)")
    fig.suptitle(f"Flujo de energía por la boca de la bahía — {scen}\n"
                 f"(+ = hacia la bahía)  •  línea punteada vertical = forzado→libre",
                 fontsize=12)
    fig.tight_layout()
    out = FIGS_DIR / f"flujo_bt_bc_{scen}.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] {out}")


def tabla_resumen():
    hdr = (f"{'Escenario':<11} {'Run':<7} {'Etapa':<6} "
           f"{'<F_bt>[kW]':>11} {'<F_bc>[kW]':>11} "
           f"{'max|F_bt|':>10} {'max|F_bc|':>10}")
    print(hdr); print("-" * len(hdr))
    for scen in SCENARIOS:
        for run in RUNS:
            for et in ETAPAS:
                try:
                    d = cargar_flujo(scen, run, et, npz_dir=SAVE_DIR)
                except FileNotFoundError:
                    continue
                Fbt = _integrar(d.get("Fy_bt_sec"), d.get("dx_sec"))
                Fbc = _integrar(d.get("Fy_sec"),    d.get("dx_sec"))
                if Fbt is None or Fbc is None:
                    continue
                print(f"{scen:<11} {run:<7} {et:<6} "
                      f"{np.nanmean(Fbt)/1e3:>11.2f} {np.nanmean(Fbc)/1e3:>11.2f} "
                      f"{np.nanmax(np.abs(Fbt))/1e3:>10.2f} {np.nanmax(np.abs(Fbc))/1e3:>10.2f}")


if __name__ == "__main__":
    print("Tabla-resumen del flujo en la boca [kW] (+ = hacia la bahía):\n")
    tabla_resumen()
    print()
    for scen in SCENARIOS:
        plot_scenario(scen)
    print("\nFiguras en:", FIGS_DIR)
