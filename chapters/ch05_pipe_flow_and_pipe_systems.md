---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# 5. Pipe Flow and Pipe Systems

```{admonition} Core ideas
:class: tip
1. The head lost to wall friction in a pipe is given by the **Darcy–Weisbach equation**, $h_f = f\,(L/D)\,V^2/2g$. All the physics is in the **friction factor** $f$.
2. In turbulent flow $f$ depends on the **Reynolds number** and the **relative roughness** $k_s/D$. The **Colebrook–White equation** and the **Moody diagram** give it; the **Swamee–Jain** formula gives it explicitly.
3. Fittings, bends and valves cause **minor losses** $K\,V^2/2g$. In short pipes they can be larger than the friction loss.
4. Pipe systems are solved with two rules: **continuity at every junction** and **one head at every node**. In series the losses add; in parallel the losses are equal; in a network the losses around each loop sum to zero (**Hardy Cross**).
5. Low pressure limits **siphons** and pump **suction** lines: the absolute pressure must stay well above the vapour pressure.
6. Stopping the flow suddenly causes **water hammer**, a pressure rise $\Delta p = \rho c \Delta V$ that can burst a pipe.
```

## 5.1 Head loss in pipes: the Darcy–Weisbach equation

Between two sections 1 and 2 of a pipe, the energy equation of Chapter 3 includes a **head loss** $h_L$:

$$
\frac{p_1}{\gamma} + \frac{V_1^2}{2g} + z_1
= \frac{p_2}{\gamma} + \frac{V_2^2}{2g} + z_2 + h_L .
$$

(We take the kinetic energy correction factor $\alpha \approx 1$, which is accurate for turbulent pipe flow.) The head loss has two parts: the **friction loss** $h_f$ along straight pipe and the **minor losses** $h_m$ at fittings, so $h_L = h_f + \sum h_m$.

In steady uniform flow, Section 4.2 showed that the wall shear stress is $\tau_0 = \gamma R S_f$, with $R = D/4$ and $S_f = h_f/L$. Dimensional analysis (Chapter 9) shows that $\tau_0$ must be proportional to $\rho V^2$. We define the **Darcy friction factor** $f$ by

$$
\tau_0 = \frac{f}{8}\, \rho V^2 .
$$

Equating the two expressions for $\tau_0$ gives the **Darcy–Weisbach equation**:

$$
h_f = f\, \frac{L}{D}\, \frac{V^2}{2g} .
$$

It is dimensionally homogeneous and valid for any fluid, laminar or turbulent flow, and any pipe material. In laminar flow, $f = 64/Re$ (Chapter 4). It also links $f$ to the shear velocity: $u_* = V\sqrt{f/8}$.

The **energy grade line** (EGL) is the line of total head $p/\gamma + V^2/2g + z$ along a pipe. The **hydraulic grade line** (HGL) lies $V^2/2g$ below it and shows the level to which water would rise in a piezometer. In a uniform pipe both lines fall with slope $S_f$; at a fitting they drop suddenly. Wherever the pipe lies **above** the HGL, the pressure is below atmospheric.

## 5.2 The friction factor

For turbulent flow, Nikuradse (1933) measured $f$ in pipes roughened with uniform sand. Commercial pipes have irregular roughness, and their behaviour in the transitional regime is smoother than in Nikuradse's tests. Colebrook (1939), using his experiments with White, proposed a single equation that covers smooth, transitional and fully rough commercial pipes. This is the **Colebrook–White equation**:

$$
\frac{1}{\sqrt{f}} = -2\log_{10}\!\left(\frac{k_s}{3.7D} + \frac{2.51}{Re\sqrt{f}}\right).
$$

When $k_s = 0$ it gives the smooth-pipe law; when $Re \to \infty$ it gives the fully rough law $1/\sqrt{f} = 2\log_{10}(3.7D/k_s)$, which follows from the rough log law of Section 4.5.

The equation is **implicit**: $f$ appears on both sides. Write $x = 1/\sqrt{f}$. The fixed-point iteration

$$
x_{k+1} = -2\log_{10}\!\left(\frac{k_s}{3.7D} + \frac{2.51\, x_k}{Re}\right)
$$

converges in a few steps from any reasonable start (for example $x_0 = 8$). Swamee and Jain (1976) fitted an **explicit** formula. For $10^{-6} \le k_s/D \le 10^{-2}$ and $5000 \le Re \le 10^8$ it is usually within 1% of Colebrook–White, and within 3% at worst (rough pipes at low $Re$):

$$
f = \frac{0.25}{\left[\log_{10}\!\left(\dfrac{k_s}{3.7D} + \dfrac{5.74}{Re^{0.9}}\right)\right]^2} .
$$

**The Moody diagram.** Moody (1944) plotted the Colebrook–White equation as a chart of $f$ against $Re$ for many values of $k_s/D$. Figure 5.1 is generated in the same way. Read it from right to left: at high $Re$ every curve becomes **horizontal** (fully rough flow, $f$ depends only on $k_s/D$); at lower $Re$ the curves bend up towards the smooth-pipe line (transitional flow); on the far left is the laminar line $f = 64/Re$.

```{code-cell} python
:tags: [hide-input]
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullLocator
from IPython.display import Markdown

plt.rcParams.update({"font.size": 11})
ink, ink2, c1 = "#0b0b0b", "#52514e", "#2a78d6"
c2, c3 = "#d9822b", "#3a9a5b"
rho, g, nu = 998.2, 9.81, 1.004e-6       # water at 20 °C
gamma = rho * g

def tidy(ax):
    ax.grid(True, which="major", color="#e4e3df"); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)

def friction_factor(Re, ks_D, n_iter=40):
    """Darcy f: 64/Re if Re < 2000, otherwise Colebrook-White by fixed-point iteration."""
    Re = np.asarray(Re, dtype=float)
    x = np.full_like(Re, 8.0)                      # x = 1/sqrt(f)
    for _ in range(n_iter):
        x = -2.0 * np.log10(ks_D / 3.7 + 2.51 * x / Re)
    f = np.where(Re < 2000, 64.0 / Re, 1.0 / x**2)
    return float(f) if f.ndim == 0 else f

def swamee_jain(Re, ks_D):
    return 0.25 / np.log10(ks_D / 3.7 + 5.74 / Re**0.9) ** 2

fig, ax = plt.subplots(figsize=(6.5, 4.6))
Re_t = np.logspace(np.log10(4000), 8, 300)
for e in [0.05, 0.01, 0.002, 5e-4, 1e-4, 1e-5, 0.0]:
    f = friction_factor(Re_t, e)
    ax.loglog(Re_t, f, color=c1 if e > 0 else ink, lw=1.6)
    ax.text(1.15e8, f[-1], "smooth" if e == 0 else f"{e:.5f}".rstrip("0"), va="center",
            color=c1 if e > 0 else ink, clip_on=False)
ax.text(1.15e8, 0.088, "$k_s/D$", color=c1, clip_on=False)
Re_l = np.logspace(np.log10(600), np.log10(4000), 50)
ax.loglog(Re_l[Re_l <= 2000], 64 / Re_l[Re_l <= 2000], color=c2, lw=2)
ax.loglog(Re_l[Re_l >= 2000], 64 / Re_l[Re_l >= 2000], color=c2, lw=1.2, ls=":")
ax.axvspan(2000, 4000, color="#f6f5f1")
ax.text(650, 0.0105, "laminar\n$f = 64/Re$", color=c2)
ax.text(3.0e6, 0.09, "fully rough", color=ink2)
ax.set_xlim(600, 1e8); ax.set_ylim(0.006, 0.1)
ticks = [0.006, 0.008, 0.01, 0.015, 0.02, 0.03, 0.04, 0.06, 0.08, 0.1]
ax.yaxis.set_major_locator(FixedLocator(ticks)); ax.yaxis.set_minor_locator(NullLocator())
ax.set_yticklabels([f"{t:g}" for t in ticks])
ax.set_xlabel("Reynolds number $Re = VD/\\nu$", color=ink2)
ax.set_ylabel("Friction factor $f$", color=ink2)
tidy(ax); fig.subplots_adjust(left=0.13, right=0.84, bottom=0.13, top=0.97); plt.show()
```

*Figure 5.1. Moody diagram generated from the Colebrook–White equation (Colebrook, 1939), with the laminar line $f = 64/Re$. Numbers on the right are the relative roughness $k_s/D$. The shaded band ($2000 < Re < 4000$) is the transition zone, where $f$ is uncertain.*

**Roughness of pipe materials.** To use Figure 5.1 we need the equivalent sand roughness $k_s$ of the pipe wall (Section 4.5). Table 5.1 lists typical values. They are rough guides: $k_s$ varies between manufacturers, and it can grow by a factor of ten as old metal pipes corrode or collect deposits.

| Material (new unless stated) | $k_s$ (mm) |
|---|---|
| Glass, drawn copper or brass | 0.0015 |
| PVC, polyethylene | 0.0015–0.007 |
| Commercial or welded steel | 0.045 |
| Asphalted cast iron | 0.12 |
| Galvanised iron | 0.15 |
| Cast iron | 0.26 |
| Concrete | 0.3–3.0 |
| Riveted steel | 0.9–9.0 |
| Old, corroded cast iron | 1–5 |

*Table 5.1. Typical equivalent sand roughness values for common pipe materials. For design, use values recommended by the pipe manufacturer or local standards.*

## 5.3 Three types of pipe problem

Every single-pipe friction problem relates five quantities: $Q$ (or $V$), $D$, $L$, $h_f$ and $k_s$ (with $\nu$ known). There are three classic types:

1. **Find the head loss** $h_f$, given $Q$, $D$, $L$, $k_s$. This is direct: compute $V$ and $Re$, then $f$, then $h_f$ (Example 5.1).
2. **Find the discharge** $Q$, given $h_f$, $D$, $L$, $k_s$. $Re$ is unknown because $V$ is unknown. Guess $f$ (0.02 is a good start), compute $V$ from Darcy–Weisbach, update $Re$ and $f$, and repeat. It usually converges in two or three steps (Example 5.2).
3. **Find the diameter** $D$, given $Q$, $h_f$, $L$, $k_s$. This is the design problem. Since $h_f$ falls steeply as $D$ grows (roughly as $D^{-5}$), a bisection on $D$ is simple and safe. The answer is then rounded **up** to the next commercial size (Example 5.3).

The $D^{-5}$ dependence is worth remembering. With $V = 4Q/\pi D^2$, Darcy–Weisbach becomes

$$
h_f = \frac{8 f L}{g \pi^2 D^5}\, Q^2 = r\, Q^2 ,
$$

where $r$ is the **pipe resistance coefficient**. A 10% smaller pipe gives about 70% more head loss for the same discharge.

## 5.4 The Hazen–Williams formula

Water-supply engineers often use the empirical **Hazen–Williams formula**. In SI units (m, m³ s⁻¹):

$$
h_f = \frac{10.67\, L\, Q^{1.852}}{C^{1.852} D^{4.87}} ,
$$

where $C$ is the Hazen–Williams coefficient: about 140–150 for plastic pipes, 130 for new cast iron or steel, and 100 or less for old cast iron. Like the Manning equation (Chapter 6), it is **not dimensionally homogeneous**: the constant 10.67 is valid only in SI units. It also ignores viscosity, so it applies only to water at ordinary temperatures in turbulent flow, in pipes larger than about 50 mm. The Darcy–Weisbach equation has no such limits and is preferred for new work. Example 5.1 compares the two.

## 5.5 Minor losses

At entrances, exits, bends, valves and changes of section, the flow separates (Section 4.6) and the eddies in the wake dissipate energy. Each such loss is written as

$$
h_m = K\, \frac{V^2}{2g} ,
$$

where $K$ is the **loss coefficient** and $V$ is normally the velocity in the pipe downstream of the fitting.

**Sudden expansion.** This is the one minor loss that can be derived. The momentum equation (Chapter 3) applied between the expansion and the point where the jet has spread over the larger pipe, combined with the energy equation, gives the **Borda–Carnot** loss:

$$
h_m = \frac{(V_1 - V_2)^2}{2g}
= \left(1 - \frac{A_1}{A_2}\right)^2 \frac{V_1^2}{2g},
$$

based on the upstream velocity $V_1$. When a pipe discharges into a reservoir, $A_2 \to \infty$ and $K = 1$: the **exit loss** is the whole velocity head.

| Fitting | $K$ |
|---|---|
| Entrance, sharp-edged (flush with wall) | 0.5 |
| Entrance, well rounded | 0.05 |
| Entrance, projecting into reservoir | 0.8–1.0 |
| Exit into a reservoir | 1.0 |
| Sudden contraction (small pipe velocity) | up to 0.5 |
| 90° bend, smooth (radius ≈ $2D$) | 0.2–0.3 |
| 90° elbow, standard or mitred | 0.9–1.1 |
| Gate valve, fully open | 0.15–0.2 |
| Globe valve, fully open | about 10 |
| Swing check valve | about 2 |

*Table 5.2. Approximate minor-loss coefficients for turbulent flow. Values vary with size and manufacturer; use manufacturers' data for design.*

A minor loss can be expressed as an **equivalent length** of straight pipe, $L_e = K D/f$. With $f = 0.02$, a fully open globe valve ($K = 10$) is equivalent to 500 diameters of pipe. Minor losses are "minor" in long pipelines but can dominate in short pipes inside pumping stations and treatment plants.

## 5.6 Pipes in series and in parallel

**Series.** Pipes of different diameter connected end to end carry the same discharge, and their losses add:

$$
H = \sum_i \left(f_i \frac{L_i}{D_i} + \sum K_i\right)\frac{V_i^2}{2g},
$$

where $H$ is the total head available (for example the difference between two reservoir levels).

**Parallel.** When a pipe splits into branches that join again, every branch connects the same two nodes. Each node has only one head, so **every branch has the same head loss**, and the discharges add: $Q = \sum Q_i$. For given total $Q$, guess the common head loss $h$, find each $Q_i = \sqrt{h/r_i}$, and adjust $h$ until the sum equals $Q$. Large branches take most of the flow, because $Q_i \propto D_i^{2.5}$ roughly.

## 5.7 Branching pipes

In a **branching** system, three or more pipes meet at a junction J and lead to reservoirs at different levels $z_i$. The classic case is the **three-reservoir problem**. The unknown is the head $H_J$ at the junction. For any trial value, each pipe carries

$$
Q_i = \operatorname{sign}(z_i - H_J)\sqrt{\frac{|z_i - H_J|}{r_i}},
$$

counted positive towards the junction. Continuity requires $\sum Q_i = 0$. If the sum is positive, $H_J$ is too low; if negative, too high. A bisection on $H_J$ between the lowest and highest reservoir levels always finds the answer. The highest reservoir always supplies water and the lowest always receives it; the middle one may do either, depending on $H_J$.

## 5.8 Pipe networks and the Hardy Cross method

Water distribution systems are **networks** with closed loops. Two rules hold:

1. **Continuity at each node**: inflow equals outflow plus demand.
2. **Energy around each loop**: going round any closed loop and returning to the start, the head losses (with sign) must sum to zero, $\sum h_f = 0$, because each node has one head.

Hardy Cross (1936) proposed a simple iterative method. Start with any set of discharges that satisfies continuity. With $h_f = r Q|Q|$ (sign follows the flow direction), a loop is generally out of balance. A uniform correction $\Delta Q$ added to every pipe in the loop keeps continuity and, to first order, removes the imbalance:

$$
\Delta Q = -\frac{\sum r Q|Q|}{2 \sum r |Q|} .
$$

For the Hazen–Williams form $h_f = rQ^{1.852}$, the 2 becomes 1.852. Corrections are applied to all loops; a pipe shared by two loops receives both. The process repeats until the corrections are negligible.

Figure 5.2 shows a small two-loop network: 100 L s⁻¹ enters at node A and is drawn off at C, D, E and F. We take a constant $f = 0.02$ for clarity (Colebrook–White could be used to update $f$ at each step). Table 5.3 shows how fast the corrections shrink.

```{code-cell} python
:tags: [hide-input]
nodes = {"A": (0, 1), "B": (1, 1), "C": (2, 1), "D": (2, 0), "E": (1, 0), "F": (0, 0)}
# name: (from, to, L m, D m, initial Q m3/s in from->to direction)
pipes = {"AB": ("A", "B", 600, 0.30, 0.050), "BC": ("B", "C", 600, 0.25, 0.030),
         "CD": ("C", "D", 400, 0.20, 0.000), "ED": ("E", "D", 600, 0.20, 0.030),
         "FE": ("F", "E", 600, 0.20, 0.030), "AF": ("A", "F", 400, 0.25, 0.050),
         "BE": ("B", "E", 400, 0.20, 0.020)}
demand = {"C": 0.03, "D": 0.03, "E": 0.02, "F": 0.02}
# loops, clockwise: +1 if pipe direction is clockwise
loops = {"I": [("AB", 1), ("BE", 1), ("FE", -1), ("AF", -1)],
         "II": [("BC", 1), ("CD", 1), ("ED", -1), ("BE", -1)]}
f_hc = 0.02
r = {p: 8 * f_hc * L / (g * np.pi**2 * D**5) for p, (_, _, L, D, _) in pipes.items()}
Q = {p: v[4] for p, v in pipes.items()}

history = []
for it in range(1, 21):
    dQ, resid = {}, {}
    for lp, members in loops.items():
        resid[lp] = sum(s * r[p] * Q[p] * abs(Q[p]) for p, s in members)
        dQ[lp] = -resid[lp] / (2 * sum(r[p] * abs(Q[p]) for p, s in members))
    for lp, members in loops.items():
        for p, s in members:
            Q[p] += s * dQ[lp]
    history.append((it, dQ["I"], dQ["II"], resid["I"], resid["II"]))
    if max(abs(v) for v in dQ.values()) < 1e-6:
        break

fig, ax = plt.subplots(figsize=(6.5, 3.6))
for p, (a, b, L, D, _) in pipes.items():
    (xa, ya), (xb, yb) = nodes[a], nodes[b]
    ax.plot([xa, xb], [ya, yb], color=ink2, lw=1.5 + 10 * D, solid_capstyle="round", zorder=1)
    if Q[p] < 0:
        (xa, ya), (xb, yb) = (xb, yb), (xa, ya)
    xm, ym = (xa + xb) / 2, (ya + yb) / 2
    ax.annotate("", xy=(xm + 0.08 * (xb - xa), ym + 0.08 * (yb - ya)),
                xytext=(xm - 0.08 * (xb - xa), ym - 0.08 * (yb - ya)),
                arrowprops=dict(arrowstyle="-|>", color=c1, lw=2), zorder=3)
    off = (0, 0.1) if ya == yb else (0.07, 0)
    ax.text(xm + off[0], ym + off[1], f"{abs(Q[p])*1e3:.1f}", color=c1,
            ha="left" if ya != yb else "center", va="center", fontsize=11)
for n, (x, y) in nodes.items():
    ax.scatter([x], [y], s=260, color="white", edgecolor=ink, zorder=4)
    ax.text(x, y, n, ha="center", va="center", color=ink, zorder=5)
ax.annotate("100 in", xy=(-0.07, 1.04), xytext=(-0.5, 1.25), color=c2,
            arrowprops=dict(arrowstyle="-|>", color=c2, lw=1.5))
for n, q in demand.items():
    x, y = nodes[n]
    dy = 0.3 if y == 1 else -0.3
    ax.annotate(f"{q*1e3:.0f} out", xy=(x + 0.05, y + dy), xytext=(x + 0.05, y + 0.06 * np.sign(dy)),
                color=c2, va="bottom" if dy > 0 else "top",
                arrowprops=dict(arrowstyle="-|>", color=c2, lw=1.5))
ax.text(0.5, 0.5, "loop I", ha="center", va="center", color=ink2)
ax.text(1.5, 0.5, "loop II", ha="center", va="center", color=ink2)
ax.set_xlim(-0.6, 2.4); ax.set_ylim(-0.45, 1.45); ax.set_aspect("equal"); ax.axis("off")
plt.tight_layout(); plt.show()
```

*Figure 5.2. A two-loop network solved by the Hardy Cross method. Arrows and blue numbers give the final pipe discharges in L s⁻¹; orange arrows show the inflow at A and the demands. Thicker lines are larger pipes (0.20–0.30 m).*

```{code-cell} python
:tags: [remove-input]
rows = [f"| {it} | {a*1e3:+.3f} | {b*1e3:+.3f} | {ra:+.4f} | {rb:+.4f} |"
        for it, a, b, ra, rb in history]
Markdown("| Iteration | ΔQ loop I (L s⁻¹) | ΔQ loop II (L s⁻¹) | Σh loop I (m) | Σh loop II (m) |\n"
         "|---|---|---|---|---|\n" + "\n".join(rows))
```

*Table 5.3. Hardy Cross iterations for the network of Figure 5.2. Σh is the loop imbalance before each correction. The imbalance shrinks quickly and is below 1 mm of head from the seventh iteration on.*

Hardy Cross is easy to follow by hand, which is why it is still taught. Modern network software (such as EPANET) solves all loop and node equations at once by Newton's method, following Todini and Pilati (1988), and handles thousands of pipes, pumps and tanks.

## 5.9 Siphons

A **siphon** is a pipe that rises above the upstream water level before falling to a lower outlet. It flows full because the outlet is below the inlet, but at the crest the pipe lies above the HGL, so the pressure there is **below atmospheric**. Chapter 1 showed that the absolute pressure can never fall below the vapour pressure. We now include the flow.

Take the upstream water surface as datum (point 1, $p = p_{atm}$, $V \approx 0$) and the crest C at height $z_C$ above it. Using **absolute** pressures, the energy equation gives

$$
\frac{p_C}{\gamma} = \frac{p_{atm}}{\gamma} - z_C - \frac{V^2}{2g} - h_{L,1C} ,
$$

where $h_{L,1C}$ is the head loss (entrance plus friction) between the intake and the crest. Setting $p_C = p_v$ gives the highest possible crest:

$$
z_{C,max} = \frac{p_{atm} - p_v}{\gamma} - \frac{V^2}{2g} - h_{L,1C} .
$$

With no flow, the last two terms vanish and $z_{C,max} \approx 10.1$ m, the result of Chapter 1. When water flows, the velocity head and the losses upstream of the crest reduce the limit, typically by 2–3 m. Dissolved air also comes out of solution well before the pressure reaches $p_v$, collects at the crest and can break the flow. For these reasons practical siphons are limited to a crest about **7–8 m** above the upstream water level, and designers often keep it lower. Figure 5.3 shows the grade lines for the siphon of Example 5.4.

```{code-cell} python
:tags: [hide-input]
# Siphon of Example 5.4: datum = upstream water surface
D_s, f_s, Ke, Kx, H_s = 0.20, 0.020, 0.5, 1.0, 6.0
s_pipe = np.array([0.0, 40.0, 120.0])          # distance along the pipe (m)
z_pipe = np.array([-3.0, 5.0, -9.0])           # pipe centreline (m)
L_s = s_pipe[-1]
hv = H_s / (Ke + f_s * L_s / D_s + Kx)         # velocity head V²/2g
s = np.linspace(0, L_s, 300)
EGL = -Ke * hv - f_s * s / D_s * hv
HGL = EGL - hv
p_lim = (101.3e3 - 2.34e3) / gamma             # (p_atm - p_v)/γ

fig, ax = plt.subplots(figsize=(6.5, 4.0))
ax.plot([-15, 0], [0, 0], color=c1, lw=3); ax.plot([L_s, L_s + 15], [-H_s, -H_s], color=c1, lw=3)
ax.plot(s_pipe, z_pipe, color=ink2, lw=6, solid_capstyle="round")
ax.plot(s, EGL, color=c3, lw=1.6, ls="--")
ax.plot(s, HGL, color=c1, lw=2)
ax.plot(s, HGL + p_lim, color=c2, lw=1.6, ls="-.")
hgl_C = -(Ke + f_s * 40 / D_s + 1) * hv     # HGL at the crest
ax.annotate("", xy=(40, 5), xytext=(40, hgl_C), arrowprops=dict(arrowstyle="<->", color=ink, lw=1.2))
ax.text(41.5, 1.0, f"{hgl_C - 5:.1f} m\n(gauge)", color=ink, va="center")
ax.text(75, HGL[np.argmin(abs(s - 75))] + p_lim + 0.6, "HGL + $(p_{atm}-p_v)/\\gamma$", color=c2)
ax.text(21, -0.8, "EGL", color=c3, va="bottom"); ax.text(21, -2.1, "HGL", color=c1, va="top")
ax.text(-14, 0.6, "upstream", color=c1); ax.text(L_s + 15, -H_s - 0.5, "downstream", color=c1, ha="right", va="top")
ax.set_xlabel("Distance along pipe (m)", color=ink2)
ax.set_ylabel("Elevation above upstream\nwater surface (m)", color=ink2)
ax.set_xlim(-15, L_s + 15); ax.set_ylim(-10, 11)
tidy(ax); plt.tight_layout(); plt.show()
```

*Figure 5.3. Profile of the siphon of Example 5.4 (thick grey line). The pipe crest lies above the HGL, so its pressure head is negative (gauge), as marked. The dash-dot line lies $(p_{atm}-p_v)/\gamma$ above the HGL: any part of the pipe above it would have a pressure below the vapour pressure.*

## 5.10 Pumps in pipelines

A pump adds head $H_p$ to the flow. The energy equation from a lower reservoir (level $z_1$) to a higher one ($z_2$) becomes

$$
H_p = (z_2 - z_1) + h_L .
$$

The right-hand side is the head the pipeline **needs**. Written as a function of $Q$, it is the **system curve**:

$$
H_{sys} = H_{st} + \left(f\frac{L}{D} + \sum K\right)\frac{Q^2}{2gA^2},
$$

where $H_{st} = z_2 - z_1$ is the **static lift**. The pump manufacturer supplies the **pump curve**, the head the pump **delivers** at each discharge, which falls as $Q$ rises. The pump operates where the two curves cross: the **operating point** (Figure 5.4). The power drawn from the motor is

$$
P = \frac{\gamma\, Q\, H_p}{\eta},
$$

where $\eta$ is the pump efficiency (typically 0.6–0.85). A good design places the operating point near the pump's best efficiency point. As a pipe ages and its roughness grows, the system curve becomes steeper and the discharge falls. Two identical pumps in **parallel** double the discharge at a given head; in **series** they double the head at a given discharge.

```{code-cell} python
:tags: [hide-input]
# Pipeline of Example 5.5
L_p, D_p, K_p, H_st = 1500.0, 0.30, 6.0, 20.0
A_p = np.pi * D_p**2 / 4
pump = lambda Q: 50.0 - 900.0 * Q**2           # pump curve, H in m, Q in m3/s

def H_system(Q, ks):
    V = Q / A_p
    f = friction_factor(V * D_p / nu, ks / D_p)
    return H_st + (f * L_p / D_p + K_p) * V**2 / (2 * g)

def operating_point(ks):
    lo, hi = 1e-4, 0.23
    for _ in range(60):                           # bisection on pump - system
        mid = (lo + hi) / 2
        if pump(mid) - H_system(mid, ks) > 0: lo = mid
        else: hi = mid
    return mid, pump(mid)

Qs = np.linspace(0.002, 0.2, 200)
fig, ax = plt.subplots(figsize=(6.5, 3.9))
ax.plot(Qs, pump(Qs), color=c2, lw=2.5)
for ks, lab, col, dy in [(0.26e-3, "new pipe", c1, -0.035), (2.0e-3, "old pipe", ink, -0.07)]:
    ax.plot(Qs, H_system(Qs, ks), color=col, lw=2)
    Qo, Ho = operating_point(ks)
    ax.scatter([Qo], [Ho], color=col, s=50, zorder=3)
    ax.text(Qo + dy, Ho + 1.5, f"{lab}\nQ = {Qo:.3f}", color=col, ha="right")
ax.axhline(H_st, color=ink2, lw=1, ls="--")
ax.text(0.005, H_st - 3, "static lift $H_{st}$", color=ink2)
ax.text(0.16, 30, "pump\ncurve", color=c2)
ax.set_xlabel("Discharge $Q$ (m³ s⁻¹)", color=ink2)
ax.set_ylabel("Head (m)", color=ink2)
ax.set_xlim(0, 0.2); ax.set_ylim(0, 55)
tidy(ax); plt.tight_layout(); plt.show()
```

*Figure 5.4. Pump curve ($H = 50 - 900Q^2$, a hypothetical pump) and system curves for the pipeline of Example 5.5 when new ($k_s$ = 0.26 mm) and after heavy corrosion ($k_s$ = 2 mm). Dots mark the operating points ($Q$ in m³ s⁻¹); the aged pipe delivers less water.*

**Cavitation and NPSH.** The lowest pressure in a pumping system is at the pump inlet, on the suction side of the impeller. If it falls to the vapour pressure, the water boils and the collapsing bubbles damage the impeller (cavitation, Chapter 1). The **net positive suction head available** is the absolute head at the pump inlet above the vapour pressure head:

$$
\text{NPSH}_A = \frac{p_{atm} - p_v}{\gamma} - z_s - h_{L,s} ,
$$

where $z_s$ is the height of the pump inlet above the suction water surface (negative if the pump is below it) and $h_{L,s}$ is the head loss in the suction pipe. The manufacturer states the **NPSH required**; a safe design needs $\text{NPSH}_A$ comfortably larger. This is why suction pipes are kept short and wide, and why large pumps are often set below the water level.

## 5.11 Water hammer

When a valve at the end of a long pipeline closes suddenly, the water next to it stops. Its momentum is destroyed by a sharp pressure rise, which travels up the pipe as a wave at speed $c$. Applying the momentum equation across the moving wave front gives the **Joukowsky equation** (Joukowsky, 1898):

$$
\Delta p = \rho\, c\, \Delta V, \qquad \Delta h = \frac{c\, \Delta V}{g},
$$

where $\Delta V$ is the change in velocity. In a rigid pipe $c = \sqrt{K/\rho} \approx$ 1480 m s⁻¹ (Chapter 1). A real pipe wall stretches, which slows the wave. For a thin-walled elastic pipe of diameter $D$, wall thickness $e$ and Young's modulus $E$:

$$
c = \sqrt{\frac{K/\rho}{1 + \dfrac{K D}{E e}}} .
$$

Table 5.4 shows that stopping a flow of only 1 m s⁻¹ raises the head by more than 100 m in metal pipes. Flexible plastic pipes have a much slower wave and a smaller surge.

```{code-cell} python
:tags: [remove-input]
K_w = 2.2e9
mats = [("Rigid pipe (limit)", None, None), ("Steel", 200, 50), ("Ductile iron", 170, 40),
        ("Concrete", 30, 12), ("PVC", 3.0, 20), ("HDPE", 0.9, 11)]
rows = []
for name, E, De in mats:
    c = np.sqrt(K_w / rho) if E is None else np.sqrt((K_w / rho) / (1 + K_w * De / (E * 1e9)))
    Es = "∞" if E is None else f"{E:g}"
    Ds = "–" if De is None else f"{De}"
    rows.append(f"| {name} | {Es} | {Ds} | {c:.0f} | {c * 1.0 / g:.0f} |")
Markdown("| Pipe | E (GPa) | D/e | c (m s⁻¹) | Δh for ΔV = 1 m s⁻¹ (m) |\n"
         "|---|---|---|---|---|\n" + "\n".join(rows))
```

*Table 5.4. Pressure-wave speed and Joukowsky head rise for water in pipes of typical materials and wall thickness ($K$ = 2.2 GPa). Values of $E$ for plastics are short-term values and vary with temperature.*

The wave reflects at the reservoir and returns to the valve after a time $2L/c$. If the valve closes in less than this time (**rapid closure**), the full Joukowsky rise occurs. Slower closure lets the reflected waves cancel part of the rise. Water hammer is controlled by closing valves slowly, by surge tanks and air vessels, and by pressure-relief valves. Pump trips cause the opposite problem: a pressure **drop** that can reach the vapour pressure and separate the water column.

## 5.12 Worked examples

The examples use the function `friction_factor` (Colebrook–White by fixed-point iteration) defined in the code of Figure 5.1.

**Example 5.1 — Head loss in a water main.** A new cast-iron main ($k_s$ = 0.26 mm) of diameter 0.30 m and length 1200 m carries 0.12 m³ s⁻¹ of water at 20 °C. Find the friction head loss using Colebrook–White and Swamee–Jain, and compare with Hazen–Williams ($C$ = 130).

```{code-cell} python
D, L, Q, ks = 0.30, 1200.0, 0.12, 0.26e-3
V = Q / (np.pi * D**2 / 4)
Re = V * D / nu
f_cw = friction_factor(Re, ks / D)
f_sj = swamee_jain(Re, ks / D)
hf = lambda f: f * L / D * V**2 / (2 * g)
hf_hw = 10.67 * L * Q**1.852 / (130**1.852 * D**4.87)
print(f"V = {V:.2f} m/s, Re = {Re:.3g}, ks/D = {ks/D:.2e}")
print(f"Colebrook-White: f = {f_cw:.4f}, h_f = {hf(f_cw):.2f} m")
print(f"Swamee-Jain:     f = {f_sj:.4f}, h_f = {hf(f_sj):.2f} m")
print(f"Hazen-Williams (C = 130):  h_f = {hf_hw:.2f} m")
```

Swamee–Jain agrees with Colebrook–White to within 1%. Hazen–Williams gives a similar answer, but only because $C$ = 130 happens to suit this pipe; its coefficient has no direct link to $k_s$.

**Example 5.2 — Discharge for a given head loss.** A commercial steel pipe ($k_s$ = 0.045 mm), 0.20 m in diameter and 800 m long, connects two reservoirs whose levels differ by 6.0 m. Neglecting minor losses, find the discharge.

```{code-cell} python
D, L, hf, ks = 0.20, 800.0, 6.0, 0.045e-3
f = 0.02                                       # first guess
for k in range(1, 6):
    V = np.sqrt(2 * g * hf * D / (f * L))      # Darcy-Weisbach solved for V
    Re = V * D / nu
    f = friction_factor(Re, ks / D)
    print(f"step {k}: V = {V:.4f} m/s, Re = {Re:.4g}, new f = {f:.5f}")
Q = V * np.pi * D**2 / 4
print(f"Q = {Q:.4f} m3/s = {Q*1e3:.1f} L/s")
```

The iteration settles after two or three steps, as expected for Type 2 problems. If the minor losses (entrance 0.5, exit 1.0) were included, they would add only 1.5 velocity heads to the $fL/D$ ≈ 67 velocity heads of friction.

**Example 5.3 — Choosing a pipe diameter.** A PVC pipeline ($k_s$ = 0.005 mm) 1500 m long must carry 0.05 m³ s⁻¹ with a friction loss of no more than 10 m. Find the minimum diameter and select a commercial size from 0.15, 0.20, 0.25 and 0.30 m.

```{code-cell} python
Q, L, hmax, ks = 0.05, 1500.0, 10.0, 0.005e-3
def head_loss(D):
    V = Q / (np.pi * D**2 / 4)
    return friction_factor(V * D / nu, ks / D) * L / D * V**2 / (2 * g)
lo, hi = 0.05, 1.0
for _ in range(50):                            # bisection: h_f decreases with D
    mid = (lo + hi) / 2
    if head_loss(mid) > hmax: lo = mid
    else: hi = mid
D_sel = min(d for d in [0.15, 0.20, 0.25, 0.30] if d >= hi)
print(f"Minimum diameter = {hi:.3f} m")
print(f"Selected D = {D_sel:.2f} m, h_f = {head_loss(D_sel):.2f} m")
```

The minimum diameter, 0.215 m, is just larger than 0.20 m, so the 0.25 m pipe must be chosen. Its head loss of 4.8 m is less than half the allowed 10 m. Choosing a pipe size is always a trade-off between the cost of the pipe and the cost of the head lost (often pumping energy).

**Example 5.4 — How high can the siphon crest be?** A siphon of diameter 0.20 m and length 120 m connects two reservoirs whose levels differ by 6.0 m. The crest is 40 m along the pipe. Take $f$ = 0.020, an entrance loss $K$ = 0.5 and an exit loss $K$ = 1.0. Find the discharge, the pressure at a crest 5.0 m above the upstream water level, and the highest possible crest. Use $p_{atm}$ = 101.3 kPa and $p_v$ = 2.34 kPa.

```{code-cell} python
D, L, L1, f, H = 0.20, 120.0, 40.0, 0.020, 6.0
Ke, Kx, zC = 0.5, 1.0, 5.0
p_atm, p_v = 101.3e3, 2.34e3
hv = H / (Ke + f * L / D + Kx)                 # velocity head V²/2g
V = np.sqrt(2 * g * hv)
Q = V * np.pi * D**2 / 4
hL_1C = (Ke + f * L1 / D) * hv                 # losses from intake to crest
pC_abs = p_atm / gamma - zC - hv - hL_1C       # absolute pressure head at crest
z_max = (p_atm - p_v) / gamma - hv - hL_1C
print(f"V = {V:.2f} m/s, Q = {Q:.4f} m3/s")
print(f"Crest pressure head: {pC_abs:.2f} m abs = {pC_abs - p_atm/gamma:.2f} m gauge")
print(f"Absolute crest pressure = {pC_abs*gamma/1e3:.1f} kPa")
print(f"Highest possible crest: {z_max:.2f} m above upstream level")
```

At a crest height of 5 m the absolute pressure is about 28 kPa, safely above the vapour pressure. The crest could not be higher than about 7.7 m, well below the static limit of 10.1 m, because the velocity head and the losses upstream of the crest use up 2.4 m of the available head. This matches the practical limit of 7–8 m quoted in Chapter 1.

**Example 5.5 — Pump operating point and NPSH.** The pump of Figure 5.4 ($H = 50 - 900Q^2$) lifts water through a static lift of 20 m along a new cast-iron pipe ($k_s$ = 0.26 mm, $D$ = 0.30 m, $L$ = 1500 m, $\sum K$ = 6). Find the operating point and the power drawn if the efficiency there is 78%. The pump inlet is 3.0 m above the sump water level and the suction-pipe loss is 0.5 m. Is the NPSH available enough if the pump needs 5.0 m?

```{code-cell} python
Qo, Ho = operating_point(0.26e-3)              # bisection, see Figure 5.4 code
eta = 0.78
P = gamma * Qo * Ho / eta
NPSH_A = (p_atm - p_v) / gamma - 3.0 - 0.5
print(f"Operating point: Q = {Qo:.4f} m3/s, H = {Ho:.1f} m")
print(f"Power = {P/1e3:.1f} kW")
print(f"NPSH available = {NPSH_A:.2f} m (required 5.0 m)")
```

The pump delivers about 124 L s⁻¹ against 36 m of head and draws about 56 kW. The NPSH available exceeds the 5.0 m required by about 1.6 m, an acceptable margin. Raising the pump by another 1.6 m, or a longer suction pipe with more loss, would use up the margin and lead to cavitation.

## 5.13 Common mistakes

- **Using the Fanning friction factor.** Some chemical engineering texts use $f_F = f/4$. The Darcy $f$ used here gives $f = 64/Re$ in laminar flow; check which one a chart uses.
- **Forgetting that $f$ depends on the flow.** $f$ changes with $Re$; in Type 2 and Type 3 problems it must be updated, not kept at the first guess.
- **Using $k_s$ in mm with $D$ in m.** The relative roughness $k_s/D$ must be dimensionless; convert both to the same unit.
- **Leaving out the exit loss.** A pipe that discharges into a reservoir loses its whole velocity head ($K = 1$).
- **Adding losses in parallel pipes.** In parallel the head losses are **equal** and the discharges add; in series the discharges are equal and the losses add.
- **Using gauge pressure for cavitation checks.** Siphon crests and pump suction must be checked with **absolute** pressure against $p_v$.

## 5.14 Exercises

1. Water at 20 °C flows at 0.03 m³ s⁻¹ through 500 m of 0.15 m galvanised-iron pipe. Find the head loss using Colebrook–White and the Moody diagram (Figure 5.1), and the error made by Swamee–Jain.
2. Two reservoirs with a level difference of 15 m are joined by 2 km of 0.40 m concrete pipe ($k_s$ = 1.0 mm), with a sharp-edged entrance, two smooth bends and an exit. Find the discharge, and the percentage of the head lost in minor losses.
3. A 0.30 m pipe 1000 m long is to be replaced by two parallel pipes of equal diameter with the same total discharge and head loss. Using a constant $f = 0.02$, find the diameter of the new pipes.
4. Repeat the Hardy Cross example of Figure 5.2 with the demand at D increased to 50 L s⁻¹ (the inflow at A rises to 120 L s⁻¹). Which pipe changes its flow direction, if any?
5. For the siphon of Example 5.4, find the highest crest if the pipe is 0.15 m in diameter. Explain why a smaller siphon can have a higher crest.
6. A 2 km steel pipeline ($D$ = 0.6 m, wall 8 mm) carries water at 1.8 m s⁻¹. A valve at the downstream end closes in 2 s. Is this a rapid closure? Estimate the maximum pressure rise in kPa and in metres of water.

---

*Figures and tables in this chapter are generated by the Python code shown. Text and figures: © Daeha Kim, CC BY-NC 4.0.*
