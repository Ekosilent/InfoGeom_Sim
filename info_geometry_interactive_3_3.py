# ============================================================
# Interactive Information Geometry - 3.3
# DST(Dirichlet) Poisson + profilo 1D + contorni + geodetiche
#
# pip install numpy scipy matplotlib   (matplotlib qualsiasi versione)
# Eseguire come .py con backend interattivo (NON inline).
# ============================================================

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from scipy.fft import dstn, idstn
from scipy.integrate import solve_ivp
from scipy.interpolate import RegularGridInterpolator

# ---- griglia ----
NX, NY = 300, 150
x = np.linspace(-10, 10, NX)
y = np.linspace(-5, 5, NY)
dx, dy = x[1] - x[0], y[1] - y[0]
X, Y = np.meshgrid(x, y)
MID = NY // 2                      # riga y ~ 0 per il profilo

# ---- stato ----
params = {"A": 1.0, "B": 0.5, "ALPHA": 1.0}
DEF_SOURCES = [
    {"x": 3.0, "y": 0.5, "amp": 1.0, "sigma": 0.6},
    {"x": -2.0, "y": -1.0, "amp": 0.6, "sigma": 1.2},
    {"x": 0.5, "y": 2.0, "amp": 0.4, "sigma": 0.8},
]
sources = [dict(s) for s in DEF_SOURCES]

# ---- autovalori Laplaciano 5-punti, BC Dirichlet (Phi=0 al bordo) ----
_ii = np.arange(1, NX + 1); _jj = np.arange(1, NY + 1)
_lx = (2.0 * np.cos(np.pi * _ii / (NX + 1)) - 2.0) / dx**2
_ly = (2.0 * np.cos(np.pi * _jj / (NY + 1)) - 2.0) / dy**2
LAM = _ly[:, None] + _lx[None, :]

# ============================================================
# CAMPI
# ============================================================

def build_density():
    rho = np.zeros_like(X)
    for s in sources:
        rho += s["amp"] * np.exp(-((X - s["x"])**2 + (Y - s["y"])**2) / s["sigma"]**2)
    return rho

def solve_field():
    rho = build_density()
    gy, gx = np.gradient(rho, dy, dx)
    I = params["A"] * rho + params["B"] * (gx**2 + gy**2)
    Phi = idstn(dstn(params["ALPHA"] * I, type=1) / LAM, type=1)
    return Phi, -2.0 * Phi

LAUNCH = [(-9.0, yy, 0.6, 0.0) for yy in (-2.5, -1.0, 0.5, 2.0)]

def integrate_geodesics(Phi):
    gy, gx = np.gradient(Phi, dy, dx)
    Fx = RegularGridInterpolator((y, x), -gx, bounds_error=False, fill_value=0.0)
    Fy = RegularGridInterpolator((y, x), -gy, bounds_error=False, fill_value=0.0)
    def rhs(t, s):
        xp, yp, vx, vy = s
        return [vx, vy, float(Fx((yp, xp))), float(Fy((yp, xp)))]
    out = []
    for ic in LAUNCH:
        sol = solve_ivp(rhs, (0.0, 45.0), ic, t_eval=np.linspace(0.0, 45.0, 900),
                        method="RK45", rtol=1e-6, atol=1e-6)
        out.append(sol.y)
    return out

# ============================================================
# FIGURA (mappa + profilo 1D + controlli)
# ============================================================

fig = plt.figure(figsize=(12, 8))
ax        = fig.add_axes([0.07, 0.46, 0.84, 0.47])   # mappa
cax       = fig.add_axes([0.925, 0.46, 0.015, 0.47]) # colorbar
ax_prof   = fig.add_axes([0.07, 0.34, 0.84, 0.09])   # profilo 1D

Phi, h00 = solve_field()
im = ax.imshow(h00, extent=[x.min(), x.max(), y.min(), y.max()],
               origin="lower", cmap="viridis", aspect="auto")
cbar = fig.colorbar(im, cax=cax); cbar.set_label("h00 = -2 Phi")
ax.set_title("Interactive Information Geometry 3.3"); ax.set_ylabel("y")

contours = [ax.contour(X, Y, h00, levels=15, colors="white", alpha=0.30, linewidths=0.5)]

(profile_line,) = ax_prof.plot(x, h00[MID, :], color="#7a3fb5", lw=1.5)
ax_prof.set_title("Profilo lungo y = 0", fontsize=8, pad=2)
ax_prof.set_ylabel("h00(y=0)", fontsize=8); ax_prof.set_xlabel("x", fontsize=8)
ax_prof.set_xlim(x.min(), x.max()); ax_prof.tick_params(labelsize=7); ax_prof.grid(alpha=0.2)

ax.text(0.01, 0.99, "Drag: sposta | Dx: aggiungi | DblClick: rimuovi",
        transform=ax.transAxes, fontsize=8, va="top", color="white",
        bbox=dict(boxstyle="round", facecolor="black", alpha=0.5))

points, geo_lines, show_geo = [], [], {"on": False}

def rebuild_points():
    for p in points: p.remove()
    points.clear()
    for s in sources:
        (p,) = ax.plot(s["x"], s["y"], "ro", ms=10, zorder=5)
        points.append(p)
rebuild_points()

# ============================================================
# UPDATE
# ============================================================

def recompute():
    global Phi, h00
    Phi, h00 = solve_field()
    im.set_data(h00)

def rescale():
    im.set_clim(h00.min(), h00.max())

def clear_contours():
    for c in contours:
        try:
            c.remove()                                   # matplotlib >= 3.8
        except Exception:
            for coll in getattr(c, "collections", []):   # fallback versioni vecchie
                coll.remove()
    contours.clear()

def update_overlays():
    clear_contours()
    contours.append(ax.contour(X, Y, h00, levels=15, colors="white", alpha=0.30, linewidths=0.5))
    profile_line.set_ydata(h00[MID, :])
    ax_prof.set_ylim(h00.min() - 0.1, h00.max() + 0.1)

def update_geodesics():
    for l in geo_lines: l.remove()
    geo_lines.clear()
    if show_geo["on"]:
        for g in integrate_geodesics(Phi):
            (ln,) = ax.plot(g[0], g[1], color="yellow", lw=1.3, alpha=0.85)
            geo_lines.append(ln)

def full_update():
    recompute(); rescale(); update_overlays(); update_geodesics()
    fig.canvas.draw_idle()

def light_update():
    # durante il drag: campo + profilo (vettoriali, veloci). Niente clim/contorni/geodetiche.
    recompute()
    profile_line.set_ydata(h00[MID, :])
    fig.canvas.draw_idle()

# ============================================================
# SLIDER
# ============================================================

sA = Slider(plt.axes([0.10, 0.22, 0.34, 0.03]), "A (rho)",      0.0, 3.0, valinit=params["A"])
sB = Slider(plt.axes([0.10, 0.17, 0.34, 0.03]), "B (|grad|^2)", 0.0, 3.0, valinit=params["B"])
sL = Slider(plt.axes([0.10, 0.12, 0.34, 0.03]), "alpha",        0.0, 3.0, valinit=params["ALPHA"])

def on_slider(_):
    params["A"], params["B"], params["ALPHA"] = sA.val, sB.val, sL.val
    full_update()
sA.on_changed(on_slider); sB.on_changed(on_slider); sL.on_changed(on_slider)

# ============================================================
# BUTTON
# ============================================================

bReset = Button(plt.axes([0.60, 0.18, 0.12, 0.05]), "Reset")
bGeo   = Button(plt.axes([0.74, 0.18, 0.12, 0.05]), "Geodetiche")

def on_reset(_):
    sources[:] = [dict(s) for s in DEF_SOURCES]
    rebuild_points()
    sA.set_val(1.0); sB.set_val(0.5); sL.set_val(1.0)
def on_geo(_):
    show_geo["on"] = not show_geo["on"]
    update_geodesics(); fig.canvas.draw_idle()
bReset.on_clicked(on_reset); bGeo.on_clicked(on_geo)

# ============================================================
# MOUSE
# ============================================================

selected = {"i": None}
def nearest(xm, ym):
    d = [np.hypot(xm - s["x"], ym - s["y"]) for s in sources]
    i = int(np.argmin(d)); return i, d[i]

def on_press(event):
    if event.inaxes != ax or event.xdata is None: return
    if event.button == 3:
        sources.append({"x": event.xdata, "y": event.ydata, "amp": 0.8, "sigma": 0.8})
        (p,) = ax.plot(event.xdata, event.ydata, "ro", ms=10, zorder=5); points.append(p)
        full_update(); return
    i, d = nearest(event.xdata, event.ydata)
    if event.dblclick and d < 0.6 and len(sources) > 1:
        sources.pop(i); points.pop(i).remove(); full_update(); return
    if d < 0.6: selected["i"] = i

def on_motion(event):
    i = selected["i"]
    if i is None or event.inaxes != ax or event.xdata is None: return
    sources[i]["x"], sources[i]["y"] = event.xdata, event.ydata
    points[i].set_data([event.xdata], [event.ydata])
    light_update()

def on_release(event):
    if selected["i"] is not None:
        selected["i"] = None; full_update()

fig.canvas.mpl_connect("button_press_event", on_press)
fig.canvas.mpl_connect("motion_notify_event", on_motion)
fig.canvas.mpl_connect("button_release_event", on_release)

plt.show()
