"""Force-based stop rules for steady RANS runs (simpleFoam), used by the 3D cases (run3d).

Residuals of separated cases often keep creeping long after the forces stand, and runs with a shedding wake
never reach a residual limit at all. The rules therefore look at the force coefficients:

Plateau: the run stops as soon as Cl(f), Cl(r) (relative to |Cl|) and Cd have stayed within PLATEAU_BAND over
the last PLATEAU_WINDOW iterations; it is then averaged over the last 200 iterations. Checked on the 13 runs of
the 2D Tyrrell validation: h/c 0.179 stops after 1200 instead of 4000 iterations, 0.134 after 1100 instead of
2108, Cl within 0.08 % and Cd within 0.15 % of the full runs (mesh scatter +-0.1 %). Oscillating runs never
reach the band. What decides is the drift per iteration (band / window): 100 iterations with 0.2 % stopped
while the forces still crept towards their final value, 0.34 % off.

Settled mean: periodic runs never reach the plateau. Once SETTLE_MIN iterations have run, a run stops when the
means over two consecutive windows differ by less than SETTLE_TOL (lift relative to |Cl|, drag relative to
Cd). The window is a whole number of shedding periods (dominant period of the rear lift, at least SETTLE_WIN
iterations), so the swing averages out. It is averaged over the last window. A run that reaches its iteration
limit unsettled is averaged over the last window as well and marked "unsteady".
SETTLE_TOL 0.5 % / 1 %: in a long reference run of a car section with a developing underbody flow the
downforce wandered slowly by ~1.5 % after it had been reached; 0.3 % was never met. The remaining uncertainty
of ~1 % is that wander, not a lack of iterations.

Coefficient rows as written by the forceCoeffs function object: Time Cd Cd(f) Cd(r) Cl Cl(f) Cl(r).
"""

from pathlib import Path

import numpy as np

PLATEAU_WINDOW = 150      # iterations
PLATEAU_BAND = 1e-3       # max - min over the window, relative to |Cl| (front and rear lift) and Cd
SETTLE_MIN = 1200         # no settled-mean check before this iteration
SETTLE_WIN = 400          # shortest averaging window, iterations (rounded up to whole shedding periods)
SETTLE_TOL = (0.005, 0.01)   # change of the window means that counts as settled: Cl(f), Cl(r) on |Cl|; Cd on Cd


def coeffs(case, name):
    """Rows (Time Cd Cd(f) Cd(r) Cl Cl(f) Cl(r)) of postProcessing/<name> written so far, or None."""
    f = Path(case) / f"postProcessing/{name}/0/coefficient.dat"
    if not f.exists():
        return None
    rows = []
    for line in f.read_text().splitlines():
        v = line.split()
        if line.startswith("#") or len(v) < 7:
            continue
        try:
            rows.append([float(x) for x in v[:7]])
        except ValueError:          # line being written
            pass
    return np.array(rows) if rows else None


def swing(a):
    """Largest max - min of Cl(f), Cl(r) (relative to |Cl|) and Cd (relative to Cd) over the rows a."""
    cl, cd = abs(a[:, 4].mean()), abs(a[:, 1].mean())
    return max(np.ptp(a[:, 5]) / cl, np.ptp(a[:, 6]) / cl, np.ptp(a[:, 1]) / cd)


def plateau(a):
    """Iteration at which the forces stand, or None."""
    if a is None or len(a) < PLATEAU_WINDOW:
        return None
    w = a[-PLATEAU_WINDOW:]
    return int(w[-1, 0]) if swing(w) < PLATEAU_BAND else None


def period(a):
    """Dominant period (iterations) of the rear lift over the last 2000 iterations of the rows a, or None."""
    t = a[:, 0]
    sel = t > t[-1] - 2000
    y, t = a[sel, 6], t[sel]
    if t[-1] - t[0] < 800:
        return None
    dt = np.median(np.diff(t))
    x = np.arange(len(y))
    y = y - np.polyval(np.polyfit(x, y, 1), x)              # without the drift
    f = np.abs(np.fft.rfft(y))
    k = np.arange(len(f))
    span = len(y) * dt
    ok = (k >= span / 800) & (k <= span / 40)               # periods of 40 ... 800 iterations
    if not ok.any():
        return None
    i = np.argmax(np.where(ok, f, 0))
    return span / i if f[i] > 3 * np.median(f[ok]) else None


def window(a):
    """Averaging window, iterations: whole shedding periods, at least SETTLE_WIN."""
    P = period(a)
    return SETTLE_WIN if P is None else int(round(np.ceil(SETTLE_WIN / P) * P))


def settled(a, min_it=SETTLE_MIN):
    """(iteration, window) once the mean of the forces has settled, or None."""
    if a is None or a[-1, 0] < min_it:
        return None
    w = window(a)
    t = a[:, 0]
    if t[-1] - t[0] < 2 * w:
        return None
    m1 = a[(t > t[-1] - 2 * w) & (t <= t[-1] - w)].mean(axis=0)
    m2 = a[t > t[-1] - w].mean(axis=0)
    cl, cd = abs(m2[4]), abs(m2[1])
    lift = max(abs(m2[5] - m1[5]), abs(m2[6] - m1[6])) / cl
    ok = lift < SETTLE_TOL[0] and abs(m2[1] - m1[1]) / cd < SETTLE_TOL[1]
    return (int(t[-1]), w) if ok else None
