import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def make_meridional_wind_patch(
        *,                    # obliga a usar argumentos por nombre → mayor claridad
        nt: int = 180,
        end_time: float = 60,          # h
        max_tau: float = 0.10,         # N m⁻²
        sigma_t: float = 1.0,          # h (efold temporal)
        center_t: float = 2.5,         # h (centro temporal)
        nx: int = 272,
        ny: int = 320,
        cx: int | None = None,         # centro x (celdas). None → nx//2
        cy: int | None = None,         # centro y (celdas). None → ny//2
        radius: int = 60,              # radio del parche (celdas)
        sigma_xy: float | None = None, # efold espacial. None → radius/2
        run_desc: str = 'windpatch',
        out_dir: str | Path = '.',
        write_bin: bool = True,
        dtype: str = '>f8',            # big-endian float64
        make_plots: bool = False
    ):
    """
    Genera un forzamiento meridional en forma de parche gaussiano espacio-tiempo.

    Devuelve:
        merid_tau  : ndarray (nx, ny, nt)
        fname_bin  : ruta del archivo .bin (None si write_bin=False)
    """
    # ---------- PARTE 1: pulso temporal -----------------------------------
    t = np.linspace(0, end_time, nt)                          # h
    gauss_t = max_tau * np.exp(-(t - center_t)**2 / (2*sigma_t**2))
    gauss_t[t > center_t] = 0                                 # media onda

    # ---------- PARTE 2: gaussiana espacial --------------------------------
    cx = cx if cx is not None else nx // 2
    cy = cy if cy is not None else ny // 2
    sigma_xy = sigma_xy if sigma_xy is not None else radius / 2

    x = np.arange(nx)[:, None]        # (nx,1)
    y = np.arange(ny)[None, :]        # (1,ny)
    r2 = (x - cx)**2 + (y - cy)**2
    gauss_xy = np.exp(-r2 / (2 * sigma_xy**2))
    gauss_xy[r2 > radius**2] = 0      # recorte opcional

    # ---------- PARTE 3: combina espacio y tiempo --------------------------
    windstress = gauss_xy[:, :, None] * gauss_t[None, None, :]  # (nx,ny,nt)
    alpha = np.radians(90)
    merid_tau = windstress * np.sin(alpha)                      # τ_mer

    # ---------- PARTE 4: escritura binaria ---------------------------------
    fname_bin = None
    if write_bin:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        fname_bin = out_dir / f'merid_{run_desc}_{nx}x{ny}.bin'

        with open(fname_bin, 'wb') as f:
            merid_tau.transpose(2, 1, 0).astype(dtype).tofile(f)
        print(f'Archivo binario escrito: {fname_bin}')

    # ---------- PARTE 5 (opcional): gráficos rápidos -----------------------
    if make_plots:
        # campo espacial (instante 3)
        plt.figure(figsize=(4.5,4))
        plt.contourf(merid_tau[:, :, 2], 20, cmap='Greens')
        plt.title('τ_mer (instante 3)')
        plt.colorbar(label='N m$^{-2}$')
        plt.xlabel('x'); plt.ylabel('y')
        plt.tight_layout(); plt.show()

        # pulso temporal
        plt.figure()
        plt.plot(t, gauss_t, '-', color='0.25')
        plt.xlabel('tiempo [h]'); plt.ylabel(r'$\tau_{mer}$ [N m$^{-2}$]')
        plt.xlim(0, end_time); plt.tight_layout(); plt.show()

    return merid_tau, fname_bin
# ---------------------------------------------------------------------------
# --------------------------  EJEMPLO RÁPIDO  -------------------------------
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    # Parche de radio 40 celdas, centrado en (100, 150)
    _, fname = make_meridional_wind_patch(
        cx=100, cy=150,
        radius=40,
        run_desc='ejemplo',
        make_plots=True
    )
    # fname → ruta del .bin listo para MITgcm
