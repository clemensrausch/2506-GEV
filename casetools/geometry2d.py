"""Parametric 2D body contours (x = flow direction, y = vertical, ground at y = 0).

A body is described as a closed, counter-clockwise list of named segments. Each segment is a polyline
(N x 2 array); consecutive segments share their end points. Segment names become OpenFOAM patch names
(body_<segment>). Corners between segments are kept sharp by the mesher. A body may consist of several
separate closed contours (a wing with a flap): `loops` lists the segment names of each contour in order;
without it all segments form one contour.

Wings in ground effect are built from tabulated profiles (both surfaces from the leading to the trailing
edge, normalised by the chord, suction side y < 0) with single_element() and multi_element(), or from a
NACA 4-digit section with naca4().
"""

from dataclasses import dataclass, field

import numpy as np
from scipy.interpolate import CubicSpline


@dataclass
class Body:
    segments: dict = field(default_factory=dict)  # name -> (N, 2) array, CCW order
    params: dict = field(default_factory=dict)
    loops: list = None                            # [[segment name, ...], ...] per closed contour

    def contours(self):
        """Segment-name lists of the closed contours."""
        return self.loops or [list(self.segments)]

    def loop_points(self, loop):
        """One closed contour as one array (duplicate junction points removed)."""
        return np.vstack([self.segments[n][:-1] for n in loop])

    def points(self):
        """All contour points (for bounding boxes)."""
        return np.vstack([self.loop_points(l) for l in self.contours()])

    def check(self, min_ground_clearance=1e-3):
        for loop in self.contours():
            segs = [self.segments[n] for n in loop]
            for a, b in zip(segs, segs[1:] + segs[:1]):
                if not np.allclose(a[-1], b[0]):
                    raise ValueError(f"contour {loop} not closed between segments")
            p = self.loop_points(loop)
            x, y = p[:, 0], p[:, 1]
            area = 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)
            if area <= 0:
                raise ValueError(f"contour {loop} must be counter-clockwise")
        y_min = self.points()[:, 1].min()
        if y_min < min_ground_clearance:
            raise ValueError(f"body too close to ground: y_min = {y_min:.4g}")
        return self


def rectangle(length=1.044, height=0.288, ride_height=0.05, x0=0.0):
    """Floating rectangle; defaults are the Ahmed body's length, height and clearance."""
    x1, y0, y1 = x0 + length, ride_height, ride_height + height
    corners = {
        "bottom": [(x0, y0), (x1, y0)],
        "rear":   [(x1, y0), (x1, y1)],
        "top":    [(x1, y1), (x0, y1)],
        "front":  [(x0, y1), (x0, y0)],
    }
    body = Body(
        {k: np.array(v, dtype=float) for k, v in corners.items()},
        dict(length=length, height=height, ride_height=ride_height, x0=x0,
             ref_length=length, ref_area_per_span=height),
    )
    return body.check()


def rotate(p, angle_deg, pivot):
    """Rotate points counter-clockwise by angle_deg about pivot."""
    a = np.radians(angle_deg)
    R = np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
    return (p - pivot) @ R.T + pivot


def resample(suction, pressure, n_side):
    """Smooth resampling of a tabulated profile (both sides LE -> TE, common LE point): one parametric
    spline over TE_p -> LE -> TE_s (C2 through the LE), cosine clustering at LE and TE."""
    raw = np.vstack([pressure[::-1], suction[1:]])
    s = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(raw, axis=0), axis=1))])
    spl = CubicSpline(s, raw)
    s_le = s[len(pressure) - 1]
    beta = np.linspace(0, np.pi, n_side)
    s_p = s_le * (1 - np.cos(beta)) / 2               # TE_p .. LE
    s_s = s_le + (s[-1] - s_le) * (1 - np.cos(beta)) / 2
    pressure_s, suction_s = spl(s_p)[::-1], spl(s_s)  # both LE -> TE again
    pressure_s[0] = suction_s[0] = raw[len(pressure) - 1]  # exact common LE point
    return suction_s, pressure_s


def single_element(suction, pressure, h_over_c, alpha=0.0, chord=1.0, x_le=0.0, n_side=200, alpha_table=0.0,
                   params=None):
    """Inverted single element wing above the ground.

    suction, pressure: tabulated surfaces (LE -> TE, normalised by the chord, LE at the origin, suction y < 0),
    given at the incidence alpha_table. The ground height h is the distance from the lowest point of the suction
    surface to the ground at alpha = 0; the incidence is then set by rotating about the quarter chord. Positive
    alpha raises the trailing edge (more downforce). Segments suction, base (blunt trailing edge), pressure.
    """
    suction_s, pressure_s = resample(suction, pressure, n_side)
    le = np.array([0.0, 0.0])
    te_mid = 0.5 * (suction[-1] + pressure[-1])
    quarter = le + 0.25 * (te_mid - le)

    y_low = rotate(suction_s, -alpha_table, quarter)[:, 1].min()
    shift = np.array([0.0, h_over_c - y_low])

    def final(p):
        p0 = rotate(p, -alpha_table, quarter) + shift
        return rotate(p0, alpha, quarter + shift) * chord + np.array([x_le, 0.0])

    suc, pre = final(suction_s), final(pressure_s)
    segments = {
        "suction":  suc,                               # LE -> TE (bottom)
        "base":     np.array([suc[-1], pre[-1]]),      # blunt trailing edge
        "pressure": pre[::-1],                         # TE -> LE (top)
    }
    body = Body(segments, dict(
        h_over_c=h_over_c, alpha=alpha, chord=chord,
        quarter_chord=(rotate(quarter, -alpha_table, quarter) + shift) * chord + [x_le, 0],
        ref_length=chord, ref_area_per_span=chord,
        x0=float(min(suc[:, 0].min(), pre[:, 0].min())),
        **(params or {}),
    ))
    return body.check(min_ground_clearance=1e-4)


def multi_element(elements, h_over_c, alpha=0.0, chord=1.0, x_le=0.0, alpha_table=0.0):
    """Inverted multi-element wing above the ground (e.g. main element and flap, assembled as tested).

    elements: {name: (suction, pressure, n_side)} with surfaces LE -> TE, normalised by the combined chord,
    the first element's LE at the origin. h is the distance from the lowest point of the first element's
    suction surface to the ground at alpha = 0; the incidence rotates the whole wing about the first element's
    quarter chord. Segments <name>_suction, <name>_base, <name>_pressure; force groups per element.
    """
    els = {}
    for name, (suc, pre, n) in elements.items():
        els[name] = resample(suc, pre, n) + ((suc[-1] + pre[-1]) / 2,)
    first = next(iter(els))
    le = np.array([0.0, 0.0])
    quarter = le + 0.25 * (els[first][2] - le)
    y_low = rotate(els[first][0], -alpha_table, quarter)[:, 1].min()
    shift = np.array([0.0, h_over_c - y_low])

    def final(p):
        p0 = rotate(p, -alpha_table, quarter) + shift
        return rotate(p0, alpha, quarter + shift) * chord + np.array([x_le, 0.0])

    segments, loops = {}, []
    for name, (suc, pre, _) in els.items():
        suc, pre = final(suc), final(pre)
        segments.update({f"{name}_suction": suc, f"{name}_base": np.array([suc[-1], pre[-1]]),
                         f"{name}_pressure": pre[::-1]})
        loops.append([f"{name}_suction", f"{name}_base", f"{name}_pressure"])
    pts = np.vstack(list(segments.values()))
    body = Body(segments, dict(
        h_over_c=h_over_c, alpha=alpha, chord=chord,
        quarter_chord=(rotate(quarter, -alpha_table, quarter) + shift) * chord + [x_le, 0],
        ref_length=chord, ref_area_per_span=chord, x0=float(pts[:, 0].min()),
        force_groups={name: [f"body_{name}_.*"] for name in els}, cofr_forces=(x_le, 0.0),
    ), loops=loops)
    return body.check(min_ground_clearance=1e-4)


def naca4(code="4412", n=161):
    """NACA 4-digit section with its standard blunt trailing edge, inverted for downforce.

    Returns (suction, pressure), both LE -> TE, normalised by the chord, LE at the origin, suction side y < 0,
    in the format single_element() takes."""
    m, p, t = int(code[0]) / 100, int(code[1]) / 10, int(code[2:]) / 100
    x = 0.5 * (1 - np.cos(np.linspace(0, np.pi, n)))
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2 + 0.2843 * x ** 3 - 0.1015 * x ** 4)
    if m:
        yc = np.where(x < p, m / p ** 2 * (2 * p * x - x ** 2), m / (1 - p) ** 2 * (1 - 2 * p + 2 * p * x - x ** 2))
        dyc = np.where(x < p, 2 * m / p ** 2 * (p - x), 2 * m / (1 - p) ** 2 * (p - x))
    else:
        yc = dyc = np.zeros_like(x)
    th = np.arctan(dyc)
    upper = np.c_[x - yt * np.sin(th), yc + yt * np.cos(th)]
    lower = np.c_[x + yt * np.sin(th), yc - yt * np.cos(th)]
    flip = np.array([1.0, -1.0])               # inverted: the cambered upper surface faces the ground
    return upper * flip, lower * flip


def plot(body, ax=None):
    import matplotlib.pyplot as plt

    ax = ax or plt.gca()
    for name, s in body.segments.items():
        ax.plot(s[:, 0], s[:, 1], "-o", ms=2, label=name)
    ax.axhline(0, color="k", lw=2)
    ax.set_aspect("equal")
    ax.legend()
    return ax
