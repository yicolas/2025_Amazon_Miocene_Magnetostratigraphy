#!/usr/bin/env python3
"""Generate the geometric accessory-symbol SVGs (viewBox -10 -10 20 20, SVG y-down).

Hand-drawn symbols (bivalve, gastropod, concretion, fining/coarsening_upward) are edited
directly; everything built from circles/ellipses/rotations lives here. Re-run after edits:
    python symbols/_make_symbols.py
Only absolute M/L/C/Z path commands are emitted (strat_symbols.py parses nothing else).
Path ids: "fill" = symbol body (colour can be overridden per use), others drawn in order.
"""
import os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
K = 0.5523

def f(x): return f"{x:.2f}".rstrip("0").rstrip(".")
def ell(cx, cy, rx, ry, rot=0.0):
    """closed ellipse as 4 cubics, optional rotation (deg)."""
    pts = [(rx, 0), (rx, K*ry), (K*rx, ry), (0, ry), (-K*rx, ry), (-rx, K*ry), (-rx, 0),
           (-rx, -K*ry), (-K*rx, -ry), (0, -ry), (K*rx, -ry), (rx, -K*ry), (rx, 0)]
    t = np.radians(rot); R = np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
    P = [R @ np.array(p) + (cx, cy) for p in pts]
    d = f"M {f(P[0][0])} {f(P[0][1])}"
    for i in range(1, 13, 3):
        d += " C " + " ".join(f"{f(P[i+j][0])} {f(P[i+j][1])}" for j in range(3))
    return d + " Z"

def rot_path(d, deg):
    """rotate every coordinate pair of an M/L/C/Z path about the origin."""
    t = np.radians(deg); c, s = np.cos(t), np.sin(t); out = []; tok = d.split(); i = 0
    while i < len(tok):
        if tok[i].isalpha(): out.append(tok[i]); i += 1
        else:
            x, y = float(tok[i]), float(tok[i+1]); out += [f(c*x - s*y), f(s*x + c*y)]; i += 2
    return " ".join(out)

def svg(name, title, parts):
    body = "\n".join(
        f'  <path id="{pid}" fill="{fc}" stroke="{ec}" stroke-width="{lw}" '
        f'stroke-linecap="round" stroke-linejoin="round"\n        d="{d}"/>' for pid, fc, ec, lw, d in parts)
    txt = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="-10 -10 20 20" width="80" height="80">\n'
           f'  <title>{title}</title>\n{body}\n</svg>\n')
    open(os.path.join(HERE, f"{name}.svg"), "w").write(txt); print("wrote", name)

# ---------------- physical sedimentary structures ----------------
svg("trough_xbed", "Trough cross-bedding (nested concave-up foresets)", [
 ("fill", "none", "#000", 0.9,
  "M -9.5 -6.5 C -6 3 6 3 9.5 -6.5 M -6.5 -6.5 C -4.2 -0.6 4.2 -0.6 6.5 -6.5 M -3.4 -6.5 C -2.2 -3.8 2.2 -3.8 3.4 -6.5 "
  "M -9.5 1.2 C -6.5 9.5 6.5 9.5 9.5 1.2 M -6 3.4 C -4 7.2 4 7.2 6 3.4"),
 ("rule", "none", "#000", 1.1, "M -9.5 -6.5 L 9.5 -6.5")])

svg("planar_xbed", "Planar/tabular cross-bedding, foresets tangential at the base", [
 ("fill", "none", "#000", 0.9,
  "M -9 -6 C -7.4 -1 -6.4 4 -3.6 6 M -5 -6 C -3.4 -1 -2.4 4 0.4 6 M -1 -6 C 0.6 -1 1.6 4 4.4 6 M 3 -6 C 4.6 -1 5.6 4 8.4 6"),
 ("rule", "none", "#000", 1.1, "M -9.5 -6 L 9.5 -6 M -9.5 6 L 9.5 6")])

svg("epsilon_xstrat", "Epsilon cross-stratification (point-bar lateral accretion, IHS)", [
 ("fill", "none", "#000", 0.9,
  "M -9.5 -5.2 C -4 -5.2 -1 -1 3.5 1.6 C 6 3 8 3.2 9.5 3.2 "
  "M -9.5 -1.4 C -4.8 -1.4 -2 2.2 2 4.6 C 4.4 6 7 6.6 9.5 6.6 "
  "M -9.5 -8.6 C -3.8 -8.6 -0.4 -4.6 4.4 -1.8 C 6.6 -0.6 8.4 -0.2 9.5 -0.2"),
 ("drape", "none", "#6d4c2f", 1.3, "M -6.6 -4.9 L -4.6 -4.4 M -3.8 -0.6 L -1.8 0.6 M -0.2 -4.4 L 1.4 -3.2"),
 ("rule", "none", "#000", 1.1, "M -9.5 8.6 L 9.5 8.6")])

svg("planar_lamination", "Planar-parallel lamination", [
 ("fill", "none", "#000", 0.9, "M -9 -6 L 9 -6 M -9 -2 L 9 -2 M -9 2 L 9 2 M -9 6 L 9 6")])

svg("hcs", "Hummocky cross-stratification (storm beds)", [
 ("fill", "none", "#000", 0.9,
  "M -9.5 -3 C -7 -8 -2.5 -8 0 -4 C 2.5 0 7 0 9.5 -4.5 "
  "M -9.5 1 C -7 -3.6 -2.6 -3.6 0 0 C 2.6 3.8 7 3.8 9.5 -0.4 "
  "M -9.5 5 C -7 0.8 -2.6 0.8 0 4 C 2.6 7.4 7 7.4 9.5 3.6")])

svg("ripples", "Ripple lamination (oscillation / starved current ripples)", [
 ("fill", "none", "#000", 0.9,
  "M -9.5 -1 C -8 -4.6 -6.6 -4.6 -5 -1 C -3.4 -4.6 -2 -4.6 -0.5 -1 C 1 -4.6 2.4 -4.6 4 -1 C 5.6 -4.6 7 -4.6 8.5 -1 "
  "M -9.5 3.6 C -8 0 -6.6 0 -5 3.6 C -3.4 0 -2 0 -0.5 3.6 C 1 0 2.4 0 4 3.6 C 5.6 0 7 0 8.5 3.6"),
 ("rule", "none", "#000", 1.0, "M -9.5 7 L 9.5 7")])

svg("pebbles", "Conglomeratic lenses / pebble stringers", [
 ("fill", "#d9c7a0", "#000", 0.8,
  " ".join([ell(-6.2, 1.5, 2.6, 1.9, 15), ell(-1.2, 2.8, 2.1, 1.6, -20), ell(3.6, 1.2, 2.8, 2.0, 5),
            ell(7.6, 3.4, 1.5, 1.2, 30), ell(-3.4, -2.2, 1.6, 1.2, -10), ell(1.6, -2.2, 1.4, 1.1, 25)])),
 ("rule", "none", "#000", 0.9, "M -9.5 5.6 C -5 7.2 5 7.2 9.5 5.6")])

svg("mottling", "Pedogenic mottling (palaeosol, varicoloured mudstone)", [
 ("fill", "#c0664a", "#5a2a1c", 0.6,
  "M -7.6 -4.8 C -5.4 -7.2 -2.6 -5.8 -3.2 -3.6 C -3.8 -1.4 -6.4 -0.6 -7.8 -2 C -9 -3.2 -8.6 -3.8 -7.6 -4.8 Z "
  "M 1.2 -6.6 C 3.6 -8 6.8 -6.8 6.2 -4.4 C 5.6 -2.4 2.4 -2.2 1.2 -3.4 C 0.2 -4.4 0.2 -5.8 1.2 -6.6 Z "
  "M -3.2 2.2 C -1 0.6 2.6 1.2 2.4 3.6 C 2.2 6 -1.6 7 -3.4 5.4 C -4.8 4.2 -4.6 3.2 -3.2 2.2 Z "
  "M 5.8 1.6 C 7.4 0.8 9.4 2 8.8 3.8 C 8.2 5.4 6 5.6 5.2 4.4 C 4.6 3.4 4.8 2.2 5.8 1.6 Z"),
 ("mottle2", "#8fae7a", "#3e5a33", 0.6,
  "M -8 5.8 C -6.8 4.6 -4.8 5.4 -5.2 6.8 C -5.6 8.2 -7.8 8.4 -8.4 7.2 C -8.8 6.6 -8.6 6.2 -8 5.8 Z "
  "M 7.2 -1.6 C 8.2 -2.6 9.8 -1.8 9.4 -0.6 C 9 0.6 7.6 0.8 7 0 C 6.6 -0.6 6.8 -1.2 7.2 -1.6 Z")])

svg("roots", "Root traces / rhizoliths (descending from a rooted horizon)", [
 ("fill", "none", "#000", 0.9,
  "M 0 -8.6 L 0 -2 C 0 1.6 -1.6 4.8 -3.2 8.4 M 0 -3.4 C 1.6 -0.6 3.8 1.6 5.6 5.6 "
  "M 0 -0.6 C -1.8 0.6 -4.6 1.6 -7 4 M 3.2 0.6 C 4.6 0.2 6.2 0.6 7.6 1.8 "
  "M -1.8 4.2 C -3.6 4.8 -5 6.2 -5.8 7.8 M 4.2 3.4 L 3.6 7 M -4.8 2.4 L -6.6 1.6"),
 ("rule", "none", "#000", 1.4, "M -8 -8.6 L 8 -8.6")])

svg("bioturbation", "Bioturbation / burrows (Planolites, Teichichnus, Palaeophycus)", [
 ("fill", "none", "#000", 1.0,
  "M -8.6 -5.4 C -8.6 -0.6 -3.4 -0.6 -3.4 -5.4 "
  "M 1.2 -6.6 C 3.2 -4.2 5.6 -8.6 8.4 -5.8 "
  "M -6 5.6 C -4 2.6 -1.6 7.6 0.6 4.8 C 2.6 2.2 4.6 6.8 7.4 3.6 "
  "M 3.4 -1.2 C 3.4 1.2 6.8 1.2 6.8 -1.2")])

svg("firmground", "Glossifungites firmground: Thalassinoides burrows descending from a discontinuity", [
 ("fill", "none", "#5a4632", 2.6,
  "M -5 -7 L -5 -1 C -5 1.6 -7.2 3.6 -8 7.2 M -5 -0.4 C -3.4 1.8 -2.4 3.4 -2.4 7.4 "
  "M 4.4 -7 L 4.4 2 C 4.4 4 6.2 5.2 7.6 7"),
 ("rule", "none", "#000", 1.6, "M -9.5 -7.2 C -6 -8.4 -3 -6.2 0 -7.4 C 3 -8.6 6 -6.2 9.5 -7.4")])

svg("tuff_ash", "Volcanic ash / tuff level", [
 ("fill", "none", "#6a4fa3", 1.0,
  "M -8 -4 L -6.5 -1 L -5 -4 M -2 -4 L -0.5 -1 L 1 -4 M 4 -4 L 5.5 -1 L 7 -4 "
  "M -5 2 L -3.5 5 L -2 2 M 1 2 L 2.5 5 L 4 2")])

# ---------------- body fossils ----------------
bone = ("M -6 -1.4 L 6 -1.4 C 6.2 -3.6 8.6 -4.8 9.6 -3.2 C 10.3 -2.1 9.4 -0.6 8.6 0 "
        "C 9.4 0.6 10.3 2.1 9.6 3.2 C 8.6 4.8 6.2 3.6 6 1.4 L -6 1.4 "
        "C -6.2 3.6 -8.6 4.8 -9.6 3.2 C -10.3 2.1 -9.4 0.6 -8.6 0 "
        "C -9.4 -0.6 -10.3 -2.1 -9.6 -3.2 C -8.6 -4.8 -6.2 -3.6 -6 -1.4 Z")
svg("bone", "Vertebrate bone / tooth (generic)", [("fill", "#efe6cf", "#000", 0.9, rot_path(bone, -32))])

svg("plant", "Plant remains / carbonised leaves", [
 ("fill", "#7a8f4a", "#1f2a10", 0.8, rot_path("M 0 -9 C 5 -5 5 4 0 8 C -5 4 -5 -5 0 -9 Z", 30)),
 ("detail", "none", "#1f2a10", 0.7, rot_path("M 0 -8 L 0 9.4 M 0 -3 L 2.6 -5.4 M 0 -3 L -2.6 -5.4 "
                                               "M 0 1 L 3 -1.4 M 0 1 L -3 -1.4 M 0 4.6 L 2.2 2.8 M 0 4.6 L -2.2 2.8", 30))])

# Purussaurus: dorsal view of a broad-snouted caimanine skull, snout up.
skull = ("M 0 -9.8 C 2.4 -9.8 3.6 -9 3.6 -7.4 C 3.6 -6.2 3.1 -5.6 3.2 -4.6 C 3.4 -2.4 4.4 0.2 5.2 2.4 C 5.8 4.2 6.3 6.4 6.4 9.2 C 5.4 9 4.4 8.4 3.7 8.4 C 3.2 8.4 3.1 9.4 2.6 9.4 L -2.6 9.4 C -3.1 9.4 -3.2 8.4 -3.7 8.4 C -4.4 8.4 -5.4 9 -6.4 9.2 C -6.3 6.4 -5.8 4.2 -5.2 2.4 C -4.4 0.2 -3.4 -2.4 -3.2 -4.6 C -3.1 -5.6 -3.6 -6.2 -3.6 -7.4 C -3.6 -9 -2.4 -9.8 0 -9.8 Z")
svg("purussaurus_skull", "Purussaurus (caimanine crocodylian) skull, dorsal view, snout to the right", [
 ("teeth", "#fbf7ea", "#000", 0.35, rot_path("M 3.15 -9 L 4.65 -8.3 L 3.15 -7.6 Z M -3.15 -9 L -4.65 -8.3 L -3.15 -7.6 Z M 2.9 -6.4 L 4.4 -5.7 L 2.9 -5 Z M -2.9 -6.4 L -4.4 -5.7 L -2.9 -5 Z M 3.05 -3.9 L 4.55 -3.2 L 3.05 -2.5 Z M -3.05 -3.9 L -4.55 -3.2 L -3.05 -2.5 Z M 3.7 -1.5 L 5.2 -0.8 L 3.7 -0.1 Z M -3.7 -1.5 L -5.2 -0.8 L -3.7 -0.1 Z M 4.45 0.9 L 5.95 1.6 L 4.45 2.3 Z M -4.45 0.9 L -5.95 1.6 L -4.45 2.3 Z", 90)),
 ("fill", "#8e44ad", "#000", 0.6, rot_path(skull, 90)),
 ("openings", "#1d1024", "none", 0.0, rot_path(
  " ".join([ell(-1.95, 3.9, 1.0, 1.35, -15), ell(1.95, 3.9, 1.0, 1.35, 15),    # orbits
            ell(-1.35, 7.3, 0.6, 0.55), ell(1.35, 7.3, 0.6, 0.55),              # supratemporal fenestrae
            ell(0, -7.9, 0.95, 0.8)]), 90)),                                    # external naris
 ("detail", "none", "#000", 0.45, rot_path("M -2.9 6.1 C -1 5.8 1 5.8 2.9 6.1 M 0 -6.2 L 0 1.4", 90))])

# Purussaurus side-profile silhouette (original generic caimanine outline; snout to the right, ground at y=3)
croc = ("M -10 -0.35 C -7.4 -0.9 -5.2 -1.6 -3.2 -2 C -1 -2.45 1.6 -2.45 3.6 -2.05 "
        "C 4.6 -1.85 5.3 -2 5.9 -2.35 C 6.4 -2.65 6.9 -2.6 7.2 -2.2 C 7.5 -1.8 8.1 -1.65 8.9 -1.5 "
        "C 9.5 -1.4 9.95 -1.25 10 -0.85 C 10.05 -0.45 9.7 -0.3 9.3 -0.3 C 8.3 -0.3 7.2 -0.15 6.2 0.15 "
        "C 5.5 0.35 5 0.55 4.5 0.8 "
        "C 4.9 1.4 5.5 1.7 5.7 2.3 C 5.8 2.7 5.9 2.9 6.4 3 L 6.5 3.3 L 4.8 3.3 C 4.6 2.8 4.1 2.3 3.4 1.35 "
        "C 1.6 1.55 -0.3 1.55 -1.6 1.3 "
        "C -2.1 1.9 -2.2 2.5 -1.6 2.95 L -0.9 3.05 L -1 3.3 L -2.7 3.3 C -3.2 2.5 -3.4 1.7 -3.6 1.05 "
        "C -5.6 0.65 -8 0.15 -10 -0.35 Z")
scutes = " ".join(f"M {x-0.35:g} {y+0.12:g} L {x:g} {y-0.45:g} L {x+0.35:g} {y+0.12:g} Z"
                  for x, y in [(-6.5,-1.2),(-5.2,-1.58),(-3.9,-1.9),(-2.6,-2.15),(-1.3,-2.35),(0,-2.43),(1.3,-2.42),(2.6,-2.26)])
svg("purussaurus", "Purussaurus, side-profile silhouette (original drawing), snout to the right", [
 ("fill", "#8e44ad", "#000", 0.35, croc + " " + scutes),
 ("eye", "#ffffff", "none", 0.0, ell(6.55, -2.05, 0.3, 0.24)),
 ("jaw", "none", "#1d1024", 0.3, "M 9.85 -0.75 C 8.8 -0.8 7.6 -0.7 6.5 -0.35")])

svg("fish", "Fish remains / fish bed", [
 ("fill", "#9fb8c8", "#000", 0.8,
  "M -7.5 0 C -4.5 -4.6 3 -4.8 6.4 -0.6 L 9.6 -3.6 L 8.6 0 L 9.6 3.6 L 6.4 0.6 C 3 4.8 -4.5 4.6 -7.5 0 Z"),
 ("eye", "#000", "none", 0.0, ell(-4.6, -0.8, 0.6, 0.6)),
 ("detail", "none", "#000", 0.6, "M -2.6 -2.6 C -1.8 -1 -1.8 1 -2.6 2.6")])
