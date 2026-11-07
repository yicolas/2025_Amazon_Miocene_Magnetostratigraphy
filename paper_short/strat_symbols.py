"""Stratigraphic accessory symbols for Fig. 2, loaded from the SVGs in ./symbols.

Each SVG (viewBox -10 -10 20 20) holds one or more <path>s:
  id="fill"   the symbol body (its colour can be overridden per use via fc=)
  other ids  extra parts (ribs, rings, openings, stipple) drawn in document order
Paths use absolute M/L/C/Z only, so no svg-parsing dependency is needed.
Edit the SVGs to restyle; the figure and legend pick the changes up.

  place(ax, x, y, "bivalve")                     point symbol, sized in points (no distortion)
  interval(ax, x0, x1, y0, y1, "fining_upward")  symbol stretched over a data-coord box
  handle("concretion", "label")                  legend handle (use with HANDLER_MAP)
"""
import os, re, xml.etree.ElementTree as ET
import numpy as np
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerTuple
from matplotlib.markers import MarkerStyle

SYMDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "symbols")
_NS = "{http://www.w3.org/2000/svg}"
_TOK = re.compile(r"[MLCZmlcz]|-?\d*\.?\d+(?:[eE][-+]?\d+)?")
_CACHE = {}

def _parse_d(d):
    verts, codes, cmd, nums = [], [], None, []
    def flush():
        n = {"M": 2, "L": 2, "C": 6}.get(cmd)
        for i in range(0, len(nums), n or 1):
            if cmd == "M":
                verts.append(nums[i:i+2]); codes.append(Path.MOVETO if i == 0 else Path.LINETO)
            elif cmd == "L":
                verts.append(nums[i:i+2]); codes.append(Path.LINETO)
            elif cmd == "C":
                c = nums[i:i+6]; verts.extend([c[0:2], c[2:4], c[4:6]]); codes.extend([Path.CURVE4]*3)
    for t in _TOK.findall(d):
        if t.isalpha():
            if cmd: flush()
            cmd, nums = t.upper(), []
            if cmd == "Z":
                verts.append(verts[-1] if verts else [0, 0]); codes.append(Path.CLOSEPOLY); cmd = None
        else:
            nums.append(float(t))
    if cmd: flush()
    v = np.array(verts, float); v[:, 1] *= -1          # SVG y-down -> matplotlib y-up
    return v, codes

def load(name):
    """Ordered list of path parts: dict(id, v, c, fc, ec, lw). Drawn in document order."""
    if name in _CACHE: return _CACHE[name]
    root = ET.parse(os.path.join(SYMDIR, f"{name}.svg")).getroot()
    out = []
    for k, p in enumerate(root.iter(f"{_NS}path")):
        v, c = _parse_d(p.get("d"))
        out.append(dict(id=p.get("id", f"p{k}"), v=v, c=c, fc=p.get("fill", "none"),
                        ec=p.get("stroke", "none"), lw=float(p.get("stroke-width", 0.8))))
    _CACHE[name] = out; return out

def _part(name, pid):
    return next((q for q in load(name) if q["id"] == pid), None)

def _marker(part):
    # two lone MOVETOs pin the extent to the full viewBox, so matplotlib's per-marker
    # normalisation scales every part identically (keeps them registered)
    v = np.vstack([part["v"], [[-10, -10], [10, 10]]])
    return MarkerStyle(Path(v, list(part["c"]) + [Path.MOVETO, Path.MOVETO]),
                       capstyle="round", joinstyle="round")

def _lines(name, x, y, size, zorder, fc=None):
    """fc overrides the colour of the part with id="fill" (the symbol body)."""
    arts = []
    for k, q in enumerate(load(name)):
        face = fc if (fc and q["id"] == "fill") else q["fc"]
        arts.append(Line2D(np.atleast_1d(x), np.atleast_1d(y), ls="none", marker=_marker(q), ms=size,
                           mfc=face, mec=(q["ec"] if q["ec"] != "none" else face),
                           mew=(q["lw"]*size/20 if q["ec"] != "none" else 0), zorder=zorder+0.01*k))
    return arts

def place(ax, x, y, name, size=9, zorder=6, fc=None):
    for a in _lines(name, x, y, size, zorder, fc): ax.add_line(a)

def interval(ax, x0, x1, y0, y1, name, zorder=6, fc=None, lw=0.6):
    """Stretch the symbol's outline bounding box onto the data box [x0,x1]x[y0,y1]."""
    f = _part(name, "fill")
    lo, hi = f["v"].min(0), f["v"].max(0)
    def m(v): return np.c_[x0+(v[:, 0]-lo[0])/(hi[0]-lo[0])*(x1-x0), y0+(v[:, 1]-lo[1])/(hi[1]-lo[1])*(y1-y0)]
    ax.add_patch(PathPatch(Path(m(f["v"]), f["c"]), fc=fc or f["fc"], ec=f["ec"], lw=lw,
                           joinstyle="miter", zorder=zorder))
    d = _part(name, "detail")
    if d is not None:   # stipple dots: one per sub-path start
        pts = m(d["v"])[np.array(d["c"]) == Path.MOVETO]
        ax.plot(pts[:, 0], pts[:, 1], ls="none", marker="o", ms=1.1, mfc="k", mec="none", zorder=zorder+0.1)

def handle(name, label, size=10, fc=None):
    arts = _lines(name, [0], [0], size, 1, fc)
    return tuple(arts) if len(arts) > 1 else arts[0], label

HANDLER_MAP = {tuple: HandlerTuple(ndivide=1, pad=0)}
