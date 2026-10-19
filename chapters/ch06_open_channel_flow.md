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

# 6. Open channel flow

```{admonition} Core ideas
:class: tip
1. An open channel has a **free surface** at atmospheric pressure. Gravity drives the flow, and the depth is not fixed by the walls: it is part of the solution.
2. In **uniform flow** the weight component along the slope balances bed friction. The **Manning equation** describes this balance, but it is empirical and works only in the units it was written for.
3. For a given discharge, the **specific energy** has a minimum at the **critical depth**. The **Froude number** $Fr = V/\sqrt{gD}$ divides flow into **subcritical** ($Fr < 1$) and **supercritical** ($Fr > 1$).
4. A **hydraulic jump** turns supercritical flow into subcritical flow and destroys energy. Because energy is lost, the jump is analysed with **momentum**, not energy.
5. In **gradually varied flow** the depth changes slowly along the channel according to $dy/dx = (S_0 - S_f)/(1 - Fr^2)$. Profiles start at a **control** and are computed step by step.
```

## 6.1 Channel geometry and types of flow

An **open channel** is any conduit in which water flows with a free surface: a river, an irrigation canal, a road gutter, or a sewer flowing part full. The pressure at the surface is atmospheric, so the flow is driven by gravity, not by a pressure difference. A channel whose cross-section and bed slope do not change along its length is called **prismatic**. Most canals are prismatic; natural rivers are not.

We describe a cross-section with a few geometric quantities, all functions of the **depth** $y$, measured vertically from the lowest point of the bed to the water surface:

- the **flow area** $A$, the area of the cross-section occupied by water;
- the **wetted perimeter** $P$, the length of the boundary in contact with water (the free surface is not included);
- the **top width** $B$, the width of the water surface;
- the **hydraulic radius** $R = A/P$;
- the **hydraulic depth** $D = A/B$, the mean depth of the section.

For a trapezoidal channel with bed width $b$ and side slopes of $m$ horizontal to 1 vertical,

$$
\begin{aligned}
A &= (b + m y)\,y, \\
P &= b + 2y\sqrt{1 + m^2}, \\
B &= b + 2 m y .
\end{aligned}
$$

A rectangular channel is the special case $m = 0$. In a very wide rectangular channel ($b \gg y$), $P \approx b$, so $R \approx D \approx y$.

```{code-cell} python
:tags: [hide-input]
import numpy as np
import matplotlib.pyplot as plt

ink, ink2, c1 = "#0b0b0b", "#52514e", "#2a78d6"
c2, c3 = "#d9822b", "#3a9a5b"
plt.rcParams.update({"font.size": 11})
g = 9.81

b, m, y = 3.0, 1.5, 1.6                 # schematic only
yb = y + 0.6                            # top of bank
fig, ax = plt.subplots(figsize=(6.5, 2.7))
bank = [(-m * yb - 0.8, yb), (-m * yb, yb), (0, 0), (b, 0), (b + m * yb, yb), (b + m * yb + 0.8, yb)]
ax.fill([p[0] for p in bank] + [b + m * yb + 0.8, -m * yb - 0.8],
        [p[1] for p in bank] + [-0.5, -0.5], color="#d9d8d3")
ax.fill([-m * y, 0, b, b + m * y], [y, 0, 0, y], color="#cfe1f6")
ax.plot([-m * y, 0, b, b + m * y], [y, 0, 0, y], color=c1, lw=3.5)
ax.plot([-m * y, b + m * y], [y, y], color=c1, lw=1.2)
ax.text(b / 2, y * 0.45, "$A$", ha="center", fontsize=16, color=ink)
ax.text(b / 2, -0.3, "wetted perimeter $P$", ha="center", va="top", color=c1)
ax.annotate("", xy=(-m * y, y + 0.25), xytext=(b + m * y, y + 0.25),
            arrowprops=dict(arrowstyle="<->", color=ink2))
ax.text(b / 2, y + 0.32, "top width $B$", ha="center", va="bottom", color=ink2)
ax.annotate("", xy=(b + 0.25, 0), xytext=(b + 0.25, y), arrowprops=dict(arrowstyle="<->", color=ink2))
ax.text(b + 0.35, y / 2, "$y$", va="center", color=ink2, fontsize=13)
ax.annotate("", xy=(0, 0.15), xytext=(b, 0.15), arrowprops=dict(arrowstyle="<->", color=ink2))
ax.text(b / 2, 0.22, "$b$", ha="center", va="bottom", color=ink2, fontsize=13)
ax.plot([-m * 0.8 - 1.2, -1.2, -1.2], [0.8, 0.8, 0], color=ink2, lw=1)
ax.plot([-m * 0.8 - 1.2, -1.2], [0.8, 0], color=ink2, lw=1)
ax.text(-1.05, 0.4, "1", color=ink2, va="center")
ax.text(-1.2 - m * 0.4, 0.88, "$m$", color=ink2, ha="center")
ax.set_xlim(-m * yb - 0.8, b + m * yb + 0.8); ax.set_ylim(-0.6, yb + 0.2)
ax.set_aspect("equal"); ax.axis("off")
plt.tight_layout(); plt.show()
```

*Figure 6.1. A trapezoidal channel section. The flow area $A$ is shaded, the wetted perimeter $P$ is the thick line along the bed and banks, and the top width $B$ is measured at the water surface. The side slope is $m$ horizontal to 1 vertical.*

Open channel flow is classified in three ways.

- **Steady or unsteady.** In steady flow the depth and velocity at a point do not change with time. A flood wave is unsteady; the flow in a canal fed at a constant rate is steady.
- **Uniform or varied.** In uniform flow the depth and velocity are the same at every section. In varied flow they change along the channel, either **gradually** (over hundreds of metres, as upstream of a weir) or **rapidly** (over a few depths, as in a hydraulic jump).
- **Subcritical or supercritical.** This depends on the Froude number, defined in §6.4.

The Reynolds number $Re = VR/\nu$ (Chapter 4) is almost always above $10^4$ in engineering channels, so the flow is **fully turbulent**.

## 6.2 Uniform flow: the Chézy and Manning equations

In uniform flow nothing accelerates, so the forces on a length $L$ of water must balance. The component of the weight along the bed, $\gamma A L \sin\theta \approx \gamma A L S_0$, where $S_0$ is the **bed slope**, is resisted by the mean shear stress $\tau_0$ on the wetted boundary, $\tau_0 P L$. Equating the two gives

$$
\tau_0 = \gamma R S_0 .
$$

In turbulent flow the boundary shear stress is proportional to $\rho V^2$. Writing $\tau_0 = \rho f V^2 / 8$, with the Darcy–Weisbach friction factor $f$ of Chapter 5, and solving for $V$ gives the **Chézy equation**,

$$
V = C \sqrt{R S_0}, \qquad C = \sqrt{8g/f},
$$

where $C$ is the **Chézy coefficient** (units $\text{m}^{1/2}\,\text{s}^{-1}$). The Chézy form is dimensionally sound, but $C$ is not constant: it depends on roughness and, weakly, on depth.

The formula used most in practice is the **Manning equation**, proposed by the Irish engineer Robert Manning in 1889:

$$
V = \frac{1}{n} R^{2/3} S_0^{1/2},
\qquad
Q = \frac{A}{n} R^{2/3} S_0^{1/2},
$$

where $n$ is **Manning's roughness coefficient**. Comparing with Chézy shows that $C = R^{1/6}/n$.

**Manning is empirical and unit-dependent.** As promised in Chapter 1, look at the dimensions. The left side is $[L T^{-1}]$ and $R^{2/3}$ is $[L^{2/3}]$. For the equation to be homogeneous, $n$ would need units of $\text{s m}^{-1/3}$, which have no physical meaning. In practice $n$ is quoted as a plain number, and the equation is then valid **only in SI units** ($V$ in m s⁻¹, $R$ in m). In US customary units the same $n$ value can be used only if the equation is rewritten as

$$
V = \frac{1.49}{n} R^{2/3} S_0^{1/2}
\quad (V \text{ in ft s}^{-1},\ R \text{ in ft}),
$$

where $1.49 \approx 3.281^{1/3}$ is the cube root of the number of feet in a metre. Never use an $n$ value in a formula written for different units without checking which form it belongs to.

Despite this weakness, the Manning equation is simple and well tested, and $n$ has been tabulated for almost every kind of surface. Table 6.1 gives typical values.

```{list-table}
:header-rows: 1
:widths: 50 20 30

* - Surface
  - Typical $n$
  - Usual range
* - Concrete, smooth finish
  - 0.012
  - 0.011–0.014
* - Concrete, rough or shotcrete
  - 0.016
  - 0.014–0.020
* - Brick or stone masonry
  - 0.017
  - 0.013–0.025
* - Earth canal, straight and clean
  - 0.022
  - 0.018–0.025
* - Earth canal, short grass
  - 0.027
  - 0.022–0.033
* - Natural river, clean and straight
  - 0.030
  - 0.025–0.035
* - Natural river, winding, pools, weeds
  - 0.045
  - 0.035–0.060
* - Floodplain with light brush and trees
  - 0.060
  - 0.040–0.080
```

*Table 6.1. Typical values of Manning's $n$ for use with the SI form of the equation. Values for natural channels vary widely; a site visit and photographs of calibrated reaches are the best guide.*

**Normal depth.** The depth at which a given discharge flows uniformly is the **normal depth** $y_n$. For a given channel, $n$, $S_0$ and $Q$, the Manning equation is one equation in one unknown, $y_n$, but it cannot be solved by hand for a trapezoid because $A$ and $R$ are nonlinear in $y$. We solve it numerically: define $F(y) = (A/n) R^{2/3} S_0^{1/2} - Q$ and find the root by bisection or Newton's method (Example 6.1). Because $Q$ increases steadily with $y$, there is exactly one root.

## 6.3 Efficient sections and compound channels

For a given area $A$, the discharge from the Manning equation is largest when the hydraulic radius is largest, that is, when the wetted perimeter is smallest. Such a section is called the **best hydraulic section**. It needs the least lining and, roughly, the least excavation for a given discharge.

For a rectangle, $P = A/y + 2y$. Setting $dP/dy = -A/y^2 + 2 = 0$ gives $A = 2y^2$, so

$$
b = 2y, \qquad R = y/2 .
$$

The best rectangle is twice as wide as it is deep. For a trapezoid the same procedure, carried out for both $y$ and $m$, gives half of a regular hexagon: side slopes at 60° to the horizontal ($m = 1/\sqrt{3}$), $b = 2y/\sqrt{3}$, and again $R = y/2$. A semicircle flowing full to the brim also has $R = y/2$, and it is the most efficient shape of all.

In practice the side slope is set by soil stability (earth canals often need $m$ = 1.5 to 3), so the best hydraulic section is a starting point for design, not the final answer.

**Compound channels.** A river in flood flows in a deep **main channel** and on shallow, rough **floodplains**. A single $R$ for the whole section gives poor results: when water first spills onto the floodplain, $P$ jumps while $A$ barely changes, and the computed discharge can even fall as the depth rises. Instead, divide the section with vertical lines at the bank tops, apply the Manning equation to each part with its own $n$ and $R$ (the dividing lines are not counted in $P$), and add the discharges.

## 6.4 Specific energy and critical flow

The energy head of open channel flow relative to the channel bed is the **specific energy**

$$
E = y + \frac{V^2}{2g} = y + \frac{Q^2}{2 g A^2}.
$$

It is the total head of Chapter 3 with the datum at the bed, so the pressure head and the elevation of the surface above the bed combine into $y$. For a rectangular channel of width $b$ we use the **discharge per unit width** $q = Q/b$ (m² s⁻¹), so $V = q/y$ and

$$
E = y + \frac{q^2}{2 g y^2}.
$$

For fixed $q$, plotting $E$ against $y$ gives the **specific energy curve** (Figure 6.2). As $y \to \infty$ the curve approaches the line $E = y$; as $y \to 0$ the velocity head grows without limit. In between there is a minimum. Setting $dE/dy = 1 - q^2/(g y^3) = 0$ gives the **critical depth**

$$
y_c = \left( \frac{q^2}{g} \right)^{1/3},
$$

at which $V_c = \sqrt{g y_c}$ and the minimum specific energy is $E_c = \tfrac{3}{2} y_c$.

```{code-cell} python
:tags: [hide-input]
fig, ax = plt.subplots(figsize=(6.0, 4.6))
yy = np.linspace(0.12, 3.0, 400)
for qq, col, lw in [(1.0, "#9fc3ec", 1.5), (2.0, c1, 2.2), (3.0, "#9fc3ec", 1.5)]:
    EE = yy + qq**2 / (2 * g * yy**2)
    ax.plot(EE, yy, color=col, lw=lw)
    ycq = (qq**2 / g) ** (1 / 3)
    ax.plot(1.5 * ycq, ycq, "o", color=col, ms=6, zorder=3)
    ax.text(1.5 * ycq - 0.08, ycq, f"$q$ = {qq:.0f}", color=c1, ha="right", va="center",
            bbox=dict(fc="white", ec="none", pad=1))
ax.plot([0, 3.0], [0, 3.0], color=ink2, lw=1, ls="--")
ax.text(2.55, 2.3, "$E = y$", color=ink2)
# alternative depths for q = 2, E = E1
q0, ya = 2.0, 1.5
E1 = ya + q0**2 / (2 * g * ya**2)
yb_alt = 0.4167
box = dict(fc="white", ec="none", pad=1)
ax.plot([E1, E1], [0, 3.0], color=c2, lw=1)
ax.plot([E1, E1], [ya, yb_alt], "s", color=c2, ms=6, zorder=3)
ax.text(E1 + 0.06, ya + 0.12, "$y_1$ (subcritical)", color=c2, va="center", bbox=box)
ax.text(E1 + 0.06, yb_alt, "$y_2$ (supercritical)", color=c2, va="center", bbox=box)
ax.set_xlim(0, 3.0); ax.set_ylim(0, 3.0)
ax.set_xlabel("Specific energy $E$ (m)", color=ink2)
ax.set_ylabel("Depth $y$ (m)", color=ink2)
ax.grid(True, color="#e4e3df"); ax.set_axisbelow(True)
for s in ("top", "right"): ax.spines[s].set_visible(False)
plt.tight_layout(); plt.show()
```

*Figure 6.2. Specific energy curves for a rectangular channel with $q$ = 1, 2 and 3 m² s⁻¹. Dots mark critical flow, where $E$ is a minimum. For $q$ = 2 m² s⁻¹ and $E$ = 1.59 m, two alternative depths are possible: 1.50 m (subcritical) and 0.42 m (supercritical).*

For any $E$ larger than $E_c$ there are two possible depths with the same energy, called **alternative depths**. The upper one is **subcritical** flow: deep and slow. The lower one is **supercritical** flow: shallow and fast. Which one actually occurs depends on the conditions upstream and downstream, not on $E$ alone.

**The Froude number.** For a channel of any shape, $dA = B\,dy$, so

$$
\frac{dE}{dy} = 1 - \frac{Q^2 B}{g A^3} = 1 - \frac{V^2}{g D} = 1 - Fr^2,
$$

where

$$
Fr = \frac{V}{\sqrt{g D}}
$$

is the **Froude number** and $D = A/B$ is the hydraulic depth. The minimum of $E$ ($dE/dy = 0$) occurs at $Fr = 1$, so **critical flow** is flow at $Fr = 1$. In a rectangle $D = y$, and the condition $Fr = 1$ gives back $y_c = (q^2/g)^{1/3}$. For other shapes we solve $Q^2 B/(g A^3) = 1$ numerically.

The quantity $\sqrt{gD}$ is the speed of a small surface wave in still water of depth $D$ (Chapter 7). The Froude number therefore compares the flow speed with the wave speed. In subcritical flow ($Fr < 1$) a disturbance can travel upstream, so the flow is **controlled from downstream**. In supercritical flow ($Fr > 1$) every disturbance is swept downstream, so the flow is **controlled from upstream**. This single fact explains much of the behaviour in the rest of this chapter.

A channel slope is called **mild** if the normal depth is subcritical ($y_n > y_c$) and **steep** if it is supercritical ($y_n < y_c$). The same channel can be mild for a small discharge and steep for a large one.

## 6.5 Transitions: humps and contractions

Consider a short, smooth rise of height $\Delta z$ in the bed of a rectangular channel. Losses over a short transition are small, so the total head is conserved, and the specific energy over the hump is

$$
E_2 = E_1 - \Delta z .
$$

On the specific energy curve (Figure 6.2) the point moves to the left at constant $q$. If the approach flow is **subcritical**, the depth **decreases** over the hump, and the water surface dips even though the bed rises. If it is **supercritical**, the depth **increases**. Both results follow directly from the shape of the curve.

The point cannot move further left than the critical point. The largest hump the flow can pass without changing the upstream conditions is

$$
\Delta z_{\max} = E_1 - E_c .
$$

If the hump is higher, the flow is **choked**: the flow over the hump stays critical, and the water upstream backs up until its specific energy is large enough, $E_1 = E_c + \Delta z$. In subcritical flow this raises the upstream water level.

A **width contraction** works in the same way: $q$ increases, and so does $E_c = 1.5\,(q^2/g)^{1/3}$; the flow chokes when $E_c$ exceeds $E_1$. Critical flow at a hump or contraction is a **control**, a section where the depth is fixed by the discharge alone. This is the principle of the broad-crested weir and the Venturi flume (§6.8).

## 6.6 The hydraulic jump

Water released under a sluice gate or at the foot of a spillway is usually supercritical. If the channel downstream is mild, the flow must return to subcritical depth. It does so suddenly, in a turbulent roller called a **hydraulic jump**. The depth before the jump, $y_1$, and the depth after it, $y_2$, are called **sequent depths** (or conjugate depths).

The jump loses energy in intense turbulence, and we do not know the loss in advance. The energy equation therefore cannot be used. Instead we use the momentum equation (Chapter 3). For a short jump on a horizontal bed, bed friction is negligible, and the only horizontal forces are the hydrostatic forces at the two ends. Per unit width of a rectangular channel,

$$
\frac{\gamma y_1^2}{2} - \frac{\gamma y_2^2}{2} = \rho q \,(V_2 - V_1).
$$

With $V = q/y$ this can be written as a conserved **specific force** (momentum function)

$$
M = \frac{q^2}{g y} + \frac{y^2}{2}, \qquad M_1 = M_2 .
$$

Dividing by $y_1^2$ and solving the resulting quadratic for $y_2/y_1$ gives the **sequent depth ratio** (the Bélanger equation)

$$
\frac{y_2}{y_1} = \frac{1}{2}\left(\sqrt{1 + 8 Fr_1^2} - 1\right),
$$

with $Fr_1 = V_1/\sqrt{g y_1}$. The head lost in the jump follows from the energy equation, now that $y_2$ is known:

$$
\Delta E = E_1 - E_2 = \frac{(y_2 - y_1)^3}{4 y_1 y_2}.
$$

The loss is positive only if $y_2 > y_1$. A jump from subcritical to supercritical flow would create energy, so it cannot happen.

```{code-cell} python
:tags: [hide-input]
Fr1 = np.linspace(1.0, 12, 300)
ratio = 0.5 * (np.sqrt(1 + 8 * Fr1**2) - 1)
loss = (ratio - 1)**3 / (4 * ratio) / (1 + Fr1**2 / 2)      # ΔE / E1

fig, axes = plt.subplots(1, 2, figsize=(6.5, 3.4))
for ax, yv, lab in [(axes[0], ratio, "$y_2 / y_1$"), (axes[1], loss * 100, "$\\Delta E / E_1$ (%)")]:
    ax.axvspan(4.5, 9, color="#e8f0fa")
    for f in (1.7, 2.5, 4.5, 9):
        ax.axvline(f, color=ink2, lw=0.8, ls=":")
    ax.plot(Fr1, yv, color=c1, lw=2)
    ax.set_xlabel("$Fr_1$", color=ink2)
    ax.set_title(lab, loc="left", color=ink)
    ax.set_xlim(1, 12)
    ax.grid(True, color="#e4e3df"); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
axes[1].text(6.75, 8, "steady\njump", ha="center", color=c1)
plt.tight_layout(); plt.show()
```

*Figure 6.3. Sequent depth ratio (left) and head loss as a fraction of the upstream specific energy (right) for a hydraulic jump in a horizontal rectangular channel. Dotted lines mark the jump types at $Fr_1$ = 1.7, 2.5, 4.5 and 9; the shaded band is the steady jump.*

The appearance of the jump depends on $Fr_1$. The usual classification, from laboratory studies by the US Bureau of Reclamation, is:

- $1 < Fr_1 < 1.7$: **undular jump**, a train of standing waves with little loss;
- $1.7 < Fr_1 < 2.5$: **weak jump**, small rollers, loss below about 20%;
- $2.5 < Fr_1 < 4.5$: **oscillating jump**, an unstable jet that sends waves far downstream (to be avoided in design);
- $4.5 < Fr_1 < 9$: **steady jump**, well balanced and insensitive to tailwater, loss 45–70%;
- $Fr_1 > 9$: **strong jump**, very rough and efficient, but hard on the structure.

A jump is roughly $6 y_2$ long. Engineers use jumps on purpose in **stilling basins** below spillways and gates, to dissipate energy before the water reaches an erodible river bed. The basin floor is set so that the tailwater depth in the river matches the sequent depth $y_2$; if the tailwater is too low, the jump is swept downstream.

## 6.7 Gradually varied flow

Upstream of a weir, or where the slope changes, the depth changes slowly over long distances. This is **gradually varied flow**. We assume that the streamlines are nearly parallel (hydrostatic pressure), and that the friction slope $S_f$ at each section equals the slope that uniform flow would need at the same depth and velocity:

$$
S_f = \frac{n^2 V^2}{R^{4/3}} .
$$

The total head is $H = z + y + V^2/2g$, where $z$ is the bed elevation. Along the channel, $dH/dx = -S_f$ and $dz/dx = -S_0$. The velocity head changes with depth as

$$
\frac{d}{dx}\left(\frac{Q^2}{2gA^2}\right) = -\frac{Q^2 B}{g A^3} \frac{dy}{dx} = -Fr^2 \frac{dy}{dx}.
$$

Adding the three terms gives the **gradually varied flow equation**

$$
\frac{dy}{dx} = \frac{S_0 - S_f}{1 - Fr^2} .
$$

The numerator compares the actual depth with the normal depth: $S_f < S_0$ when $y > y_n$. The denominator compares it with the critical depth: $Fr < 1$ when $y > y_c$. The lines $y = y_n$ and $y = y_c$ divide the channel into three **zones**: zone 1 above both lines, zone 2 between them, and zone 3 below both. Combining zone and slope type names each **water surface profile** (Table 6.2).

```{list-table}
:header-rows: 1
:widths: 12 22 12 54

* - Profile
  - Depth range
  - $dy/dx$
  - Typical location
* - M1
  - $y > y_n > y_c$
  - $+$
  - upstream of a dam or weir
* - M2
  - $y_n > y > y_c$
  - $-$
  - approaching a free overfall or a steep reach
* - M3
  - $y_n > y_c > y$
  - $+$
  - below a sluice gate, ending in a jump
* - S1
  - $y > y_c > y_n$
  - $+$
  - behind a weir on a steep slope, after a jump
* - S2
  - $y_c > y > y_n$
  - $-$
  - just downstream of a change from mild to steep
* - S3
  - $y_c > y_n > y$
  - $+$
  - below a sluice gate on a steep slope
```

*Table 6.2. Gradually varied flow profiles on mild (M) and steep (S) slopes. A positive $dy/dx$ means the depth increases in the flow direction.*

On a **horizontal** bed ($S_0 = 0$) there is no normal depth, and only H2 and H3 profiles exist; on an **adverse** slope ($S_0 < 0$) only A2 and A3. A **critical** slope ($y_n = y_c$) gives C1 and C3. All profiles approach the normal depth line asymptotically, and they meet the critical depth line almost vertically, where the gradually varied flow assumption breaks down (a jump or a free overfall).

**Controls and direction of computation.** A profile must start at a **control**, a section where the depth is known from the discharge: a weir, a gate, a free overfall, or a point of critical flow. Because subcritical flow is controlled from downstream, subcritical profiles (M1, M2, S1) are computed **upstream** from their control. Supercritical profiles (M3, S2, S3) are computed **downstream**.

**The direct step method.** Write the energy equation between two sections a distance $\Delta x$ apart, with section 1 upstream and section 2 downstream. With $\bar{S}_f$ the average friction slope of the two sections,

$$
E_1 + S_0 \Delta x = E_2 + \bar{S}_f \Delta x ,
$$

so

$$
\Delta x = \frac{E_2 - E_1}{S_0 - \bar{S}_f}.
$$

In a prismatic channel we choose a sequence of depths, compute $E$ and $S_f$ for each, and obtain the distance between them directly. No iteration is needed. (In a natural river, where the sections are fixed by survey, the **standard step method** is used instead: $\Delta x$ is known and the depth is found by iteration. It is the method in HEC-RAS.)

The code below applies the direct step method to a rectangular channel with $b$ = 10 m, $n$ = 0.030, $S_0$ = 0.0005 and $Q$ = 30 m³ s⁻¹. A weir raises the depth just upstream of it to 4.0 m. The computation stops when the depth is within 1% of the normal depth.

```{code-cell} python
Q, b, n, S0 = 30.0, 10.0, 0.030, 0.0005

def bisect(f, lo, hi, tol=1e-8):
    """Root of f between lo and hi by bisection."""
    flo = f(lo)
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if flo * f(mid) <= 0: hi = mid
        else: lo, flo = mid, f(mid)
    return 0.5 * (lo + hi)

def manning_Q(y):
    A, P = b * y, b + 2 * y
    return A / n * (A / P) ** (2 / 3) * S0 ** 0.5

def direct_step(depths):
    """Distances upstream (m) of the given depths; depths[0] is at the control."""
    y = np.asarray(depths, float)
    A, P = b * y, b + 2 * y
    V = Q / A
    E = y + V**2 / (2 * g)
    Sf = (n * V) ** 2 / (A / P) ** (4 / 3)
    dx = (E[:-1] - E[1:]) / (S0 - 0.5 * (Sf[:-1] + Sf[1:]))
    return np.concatenate([[0.0], np.cumsum(dx)]), V, E, Sf

yn = bisect(lambda y: manning_Q(y) - Q, 0.1, 10)
yc = ((Q / b) ** 2 / g) ** (1 / 3)
ys = np.linspace(4.0, 1.01 * yn, 400)
xs, _, _, _ = direct_step(ys)
print(f"Normal depth yn = {yn:.2f} m, critical depth yc = {yc:.2f} m (mild slope)")
print(f"Depth falls to 1.01 yn at {xs[-1]/1000:.1f} km upstream of the weir")
```

The backwater extends almost 9 km upstream, although the weir raises the water by only 1.25 m above normal depth. Long backwater curves are typical of mild slopes, and they decide how far upstream a dam affects flood levels, bridges and drainage outlets.

```{code-cell} python
:tags: [hide-input]
fig, ax = plt.subplots(figsize=(6.5, 3.6))
xk = -xs / 1000                          # km, flow from left to right
zb = S0 * xs                             # bed rises upstream
ax.fill_between(xk, zb - 0.6, zb, color="#d9d8d3")
ax.plot(xk, zb, color=ink2, lw=1.2)
ax.plot(xk, zb + ys, color=c1, lw=2.2, label="water surface (M1)")
ax.plot(xk, zb + yn, color=c3, lw=1.3, ls="--", label="normal depth line")
ax.plot(xk, zb + yc, color=c2, lw=1.3, ls=":", label="critical depth line")
ax.plot([0, 0], [-0.6, 2.6], color=ink, lw=4)
ax.text(-0.15, 0.6, "weir", ha="right", color=ink)
ax.text(-6.0, S0 * 6000 - 1.3, "bed", color=ink2)
ax.set_xlabel("Distance from weir (km)", color=ink2)
ax.set_ylabel("Elevation above weir bed (m)", color=ink2)
ax.set_xlim(-9.5, 0.4); ax.set_ylim(-1.0, 9.5)
ax.legend(frameon=False, loc="upper right")
ax.grid(True, color="#e4e3df"); ax.set_axisbelow(True)
for s in ("top", "right"): ax.spines[s].set_visible(False)
plt.tight_layout(); plt.show()
```

*Figure 6.4. M1 backwater profile upstream of a weir in a rectangular channel ($b$ = 10 m, $n$ = 0.030, $S_0$ = 0.0005, $Q$ = 30 m³ s⁻¹), computed by the direct step method. The vertical scale is strongly exaggerated. The water surface approaches the normal depth line asymptotically upstream and stays far above the critical depth line.*

```{code-cell} python
:tags: [remove-input]
from IPython.display import Markdown
yt = [4.0, 3.8, 3.6, 3.4, 3.2, 3.0, 2.8]
xt, Vt, Et, Sft = direct_step(yt)
rows = []
for i in range(len(yt)):
    dxs = "–" if i == 0 else f"{xt[i] - xt[i-1]:.0f}"
    rows.append(f"| {yt[i]:.2f} | {Vt[i]:.3f} | {Et[i]:.4f} | {Sft[i]*1e4:.3f} | {dxs} | {xt[i]:.0f} |")
Markdown("| $y$ (m) | $V$ (m s⁻¹) | $E$ (m) | $S_f$ (×10⁻⁴) | $\\Delta x$ (m) | $x$ (m) |\n"
         "|---|---|---|---|---|---|\n" + "\n".join(rows))
```

*Table 6.3. Direct step computation with depth steps of 0.2 m for the channel of Figure 6.4; $x$ is the distance upstream of the weir. The steps grow longer as the depth approaches normal depth, because $S_0 - \bar{S}_f$ becomes small. Smaller depth steps near $y_n$ give a more accurate profile.*

## 6.8 Flow measurement with weirs

A **weir** is an obstruction over which water flows. The flow over it acts as a control, so the upstream head $H$ above the crest fixes the discharge. Weirs are the simplest reliable way to measure discharge in a channel.

**Sharp-crested weirs.** A thin plate with a sharp upstream edge spans the channel. Applying the Bernoulli equation along streamlines over the crest, and ignoring the approach velocity, gives for a rectangular weir of crest length $b$

$$
Q = C_d \, \frac{2}{3} \sqrt{2g}\; b\, H^{3/2},
$$

where the **discharge coefficient** $C_d$ accounts for the contraction of the nappe and is about 0.62 for a full-width weir with small $H$ compared with the weir height. For small discharges a **V-notch** (triangular) weir with notch angle $\theta$ is more sensitive:

$$
Q = C_d \, \frac{8}{15} \sqrt{2g}\, \tan\frac{\theta}{2}\, H^{5/2},
$$

with $C_d \approx 0.58$ for a 90° notch. The head is measured about $4H$ upstream, beyond the drawdown over the crest.

**Broad-crested weirs.** If the crest is long enough for the streamlines to become parallel, the flow on the crest is critical (§6.5). With $H$ the head above the crest and the approach velocity ignored, $E = H$ and $y_c = \tfrac{2}{3} H$. Then $q = \sqrt{g y_c^3}$, so

$$
Q = \sqrt{g}\left(\tfrac{2}{3}\right)^{3/2} b\, H^{3/2} \approx 1.705\, b\, H^{3/2}.
$$

Real broad-crested weirs pass a little less, so a coefficient of about 0.85–0.95 is applied, depending on the shape of the upstream edge. They are robust and pass sediment well, so they are common in rivers and irrigation canals.

## 6.9 Worked examples

**Example 6.1 — Normal depth in an earth canal.** A trapezoidal earth canal has bed width 3.0 m, side slopes 2H:1V, $n$ = 0.025 and bed slope 0.0008. Find the normal depth for $Q$ = 12 m³ s⁻¹, the mean velocity and the Froude number, and decide whether the slope is mild or steep.

```{code-cell} python
bw, m, n1, S01, Q1 = 3.0, 2.0, 0.025, 0.0008, 12.0
A = lambda y: (bw + m * y) * y
P = lambda y: bw + 2 * y * np.sqrt(1 + m**2)
B = lambda y: bw + 2 * m * y

yn1 = bisect(lambda y: A(y) / n1 * (A(y) / P(y)) ** (2 / 3) * S01**0.5 - Q1, 0.01, 10)
V1 = Q1 / A(yn1)
Fr = V1 / np.sqrt(g * A(yn1) / B(yn1))
yc1 = bisect(lambda y: Q1**2 * B(y) / (g * A(y) ** 3) - 1, 0.01, 10)
print(f"Normal depth yn = {yn1:.2f} m, V = {V1:.2f} m/s, Fr = {Fr:.2f}")
print(f"Critical depth yc = {yc1:.2f} m  ->  {'mild' if yn1 > yc1 else 'steep'} slope")
```

The normal depth is 1.67 m and the flow is clearly subcritical ($Fr$ = 0.35), so the slope is mild. A velocity of 1.14 m s⁻¹ is acceptable for a grassed earth canal but close to the limit for fine sandy soils.

**Example 6.2 — Flow over a raised hump.** Water flows at 2.0 m² s⁻¹ in a wide rectangular channel at a depth of 1.50 m. A smooth hump 0.20 m high is built across the bed. Find the depth over the hump, the change in water level, and the highest hump that would not choke the flow.

```{code-cell} python
q, y1, dz = 2.0, 1.50, 0.20
Efun = lambda y: y + q**2 / (2 * g * y**2)
E1 = Efun(y1)
yc2 = (q**2 / g) ** (1 / 3)
Ec = 1.5 * yc2
y2 = bisect(lambda y: Efun(y) - (E1 - dz), yc2, E1)     # subcritical branch
print(f"Fr1 = {q / y1 / np.sqrt(g * y1):.2f}, E1 = {E1:.3f} m, yc = {yc2:.3f} m, Ec = {Ec:.3f} m")
print(f"Depth over hump = {y2:.3f} m; water surface drops by {y1 - (y2 + dz):.3f} m")
print(f"Largest hump without choking = {E1 - Ec:.3f} m")
```

The depth over the hump is 1.26 m. Because the approach flow is subcritical, the water surface falls by about 4 cm even though the bed rises by 20 cm. Any hump higher than 0.48 m would choke the flow and raise the water level upstream.

**Example 6.3 — Hydraulic jump below a sluice gate.** A sluice gate in a rectangular channel 5.0 m wide releases 15 m³ s⁻¹ at a depth of 0.30 m. Find the sequent depth, the type of jump, the head loss and the power dissipated.

```{code-cell} python
Qj, bj, yj1 = 15.0, 5.0, 0.30
qj = Qj / bj
Vj1 = qj / yj1
Frj = Vj1 / np.sqrt(g * yj1)
yj2 = 0.5 * yj1 * (np.sqrt(1 + 8 * Frj**2) - 1)
dE = (yj2 - yj1) ** 3 / (4 * yj1 * yj2)
Ej1 = yj1 + Vj1**2 / (2 * g)
power = 998.2 * g * Qj * dE
print(f"V1 = {Vj1:.1f} m/s, Fr1 = {Frj:.2f}")
print(f"Sequent depth y2 = {yj2:.2f} m")
print(f"Head loss = {dE:.2f} m ({100 * dE / Ej1:.0f}% of E1), power = {power / 1e3:.0f} kW")
```

With $Fr_1$ = 5.83 this is a steady jump. The depth rises from 0.30 m to 2.33 m, and more than half of the incoming energy, about 440 kW, is turned into turbulence and finally heat. A stilling basin needs a tailwater depth close to 2.3 m to hold the jump in place.

**Example 6.4 — Extent of the backwater.** For the weir and channel of §6.7, how far upstream is the water level raised by at least 0.30 m above its normal (no-weir) level?

```{code-cell} python
target = yn + 0.30
ys_fine = np.linspace(4.0, 1.001 * yn, 4000)
xs_fine, _, _, _ = direct_step(ys_fine)
x_target = np.interp(target, ys_fine[::-1], xs_fine[::-1])
print(f"Depth {target:.2f} m is reached {x_target/1000:.2f} km upstream of the weir")
```

The weir raises the water level by 0.30 m or more for about 4.1 km upstream. A drainage outlet or a low bridge within this distance would be affected by the weir, even though it is out of sight of the structure.

## 6.10 Common mistakes

- **Using Manning's $n$ in the wrong units.** The SI form has 1 in the numerator; the US form has 1.49. The value of $n$ is the same, but only if the formula matches the units.
- **Confusing $R$ and $D$.** The hydraulic radius $R = A/P$ goes into friction formulas; the hydraulic depth $D = A/B$ goes into the Froude number. They are equal only in a very wide rectangle.
- **Using the energy equation across a hydraulic jump.** The head loss is unknown, so the sequent depth must come from momentum. Energy is used afterwards to find the loss.
- **Choosing the wrong alternative depth.** For a given $E$ there are two depths. Decide first whether the flow is subcritical or supercritical, and search on the correct side of $y_c$.
- **Computing a profile in the wrong direction.** Subcritical profiles are computed upstream from a downstream control; supercritical profiles downstream from an upstream control.
- **Forgetting that a hump can lower the water surface.** In subcritical flow the depth over a hump falls by more than the bed rises.

## 6.11 Exercises

1. Show that the best hydraulic trapezoidal section has $R = y/2$, whatever its side slope $m$.
2. A concrete rectangular channel ($n$ = 0.014) 4 m wide on a slope of 0.002 carries 10 m³ s⁻¹. Find the normal and critical depths and classify the slope. Repeat for a slope of 0.02.
3. Water flows at 1.2 m depth and 1.5 m s⁻¹ in a rectangular channel 3 m wide. Find the minimum width of a contraction that does not choke the flow. What happens upstream if the width is made smaller still?
4. A hydraulic jump in a rectangular channel has sequent depths 0.25 m and 1.60 m. Find the discharge per unit width, $Fr_1$ and the relative head loss.
5. A broad-crested weir 6 m long carries a head of 0.45 m. Estimate the discharge, taking a coefficient of 0.90. What is the error in discharge if the head is misread by 1 cm?
6. Modify the direct step code of §6.7 to compute the M2 profile upstream of a free overfall in the same channel, starting at the critical depth at the brink. How far upstream does the drawdown reach 0.99 $y_n$?

---

*Figures and tables in this chapter are generated by the Python code shown. Text and figures: © Daeha Kim, CC BY-NC 4.0.*
