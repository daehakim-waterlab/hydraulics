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

# 2. Water at Rest

```{admonition} Core ideas
:class: tip
1. In a fluid at rest there is no shear stress, so pressure at a point **acts equally in all directions** and always acts **normal to a surface**.
2. Pressure in still water **increases linearly with depth**: $p = p_0 + \gamma h$. Equivalently, the **piezometric head** $p/\gamma + z$ is the same everywhere in a connected body of still water.
3. Engineers usually work with **gauge pressure** (relative to the atmosphere). Cavitation and gases need **absolute pressure**.
4. The force on a plane surface equals the **pressure at its centroid times its area**, but it acts at the **centre of pressure**, which is always below the centroid.
5. A curved surface is handled by components: the **horizontal force** acts on the vertical projection, and the **vertical force** equals the weight of water above (or the missing water above) the surface.
6. A floating body is in equilibrium when its weight equals the **buoyancy force**; it is stable when the **metacentre** lies above its centre of gravity.
```

## 2.1 Pressure at a point

**Hydrostatics** is the study of fluids at rest. Water at rest cannot carry shear stress, because any shear stress, however small, would make it deform continuously (Chapter 1). The only stress left is **pressure**, the normal force per unit area, $p = dF/dA$, measured in pascals (Pa = N m⁻²).

Consider a tiny wedge of still water. The pressure forces on its faces and its weight must balance. As the wedge shrinks to a point, the weight (proportional to volume) becomes negligible compared with the surface forces (proportional to area), and the balance shows that the pressure on every face is the same. This is **Pascal's law**: at a point in a fluid at rest, pressure has the same value in all directions. Pressure is therefore a scalar, and the force it produces on any surface acts **normal** to that surface.

## 2.2 Variation of pressure with depth

Take a small vertical column of still water with cross-section $dA$ and height $dz$, where $z$ is measured upward. The pressure on its bottom face is $p$ and on its top face $p + dp$. Its weight is $\gamma\, dA\, dz$. Vertical equilibrium gives

$$
p\, dA - (p + dp)\, dA - \gamma\, dA\, dz = 0,
$$

so that

$$
\frac{dp}{dz} = -\gamma .
$$

This is the **hydrostatic equation**. Pressure decreases as we go up and increases as we go down. In the horizontal direction there is no weight to balance, so pressure is constant on any horizontal plane in a connected body of still water.

Because water is practically incompressible, $\gamma$ is constant and the equation integrates directly. If $h$ is the depth below a free surface where the pressure is $p_0$,

$$
p = p_0 + \gamma h .
$$

Writing the same result with the elevation $z$,

$$
\frac{p}{\gamma} + z = \text{constant}.
$$

The quantity $p/\gamma$ is the **pressure head**, a pressure expressed as the height of a water column that would produce it. The sum $p/\gamma + z$ is the **piezometric head**. In still water it is the same at every point. In Chapter 3 we will add the velocity head to obtain the total head of a moving fluid.

**Absolute and gauge pressure.** **Absolute pressure** is measured from a perfect vacuum. **Gauge pressure** is measured from the local atmospheric pressure:

$$
p_{\text{abs}} = p_{\text{atm}} + p_{\text{gauge}}.
$$

A negative gauge pressure is called a **vacuum** or **suction**. Since the atmosphere acts on all sides of most structures, its effect cancels, and in this book a pressure is a gauge pressure unless we say otherwise. Standard atmospheric pressure, 101.3 kPa, equals a pressure head of about 10.3 m of water. The absolute pressure can never be negative, and in water it cannot fall below the vapour pressure without cavitation (Chapter 1).

## 2.3 Measuring pressure

The hydrostatic equation is the basis of all liquid-column instruments.

**Piezometer.** A **piezometer** is a simple vertical tube connected to the point of interest and open to the atmosphere. Water rises in the tube to a height $h = p/\gamma$ above the point. It is accurate and cheap, but impractical for high pressures (a 200 kPa pressure needs a tube over 20 m tall) and useless for gases or for pressures below atmospheric.

**U-tube manometer.** A **manometer** uses a second, heavier liquid, usually mercury ($S_m = 13.6$), in a U-shaped tube. To analyse any manometer, use one rule: start at a point of known pressure, **add** $\gamma \Delta h$ for each step down through a fluid, **subtract** $\gamma \Delta h$ for each step up, and stop at the point of interest. Points at the same level in the same continuous liquid have the same pressure.

**Differential manometer.** A **differential manometer** connects two points A and B, for example two sections of a water pipe, through a U-tube of mercury. If the mercury levels differ by a reading $R$, the rule above gives

$$
p_A - p_B = (\gamma_m - \gamma) R + \gamma (z_B - z_A),
$$

or, in heads,

$$
\left(\frac{p_A}{\gamma} + z_A\right) - \left(\frac{p_B}{\gamma} + z_B\right) = R\,(S_m - 1).
$$

The manometer therefore measures the difference in **piezometric head**, not simply the pressure difference. With mercury and water, each millimetre of reading corresponds to 12.6 mm of water head. We will use this result for the Venturi meter in Chapter 3.

**Bourdon gauge.** For high pressures the **Bourdon gauge** is used. It contains a curved tube of oval cross-section that tends to straighten when the internal pressure rises; a linkage turns this movement into a needle reading. It reads gauge pressure, because the outside of the tube is at atmospheric pressure. Modern **pressure transducers** use a strain gauge or a piezoelectric element on a thin diaphragm and give an electrical signal suitable for data loggers.

## 2.4 Force on a plane submerged surface

Gates, walls and tank sides must resist the pressure of water. Consider a plane surface of area $A$, inclined at angle $\theta$ to the horizontal, with one face in water. Let $y$ be the distance measured **along the plane** from the line where the plane meets the free surface, so that the depth of any point is $h = y \sin\theta$.

The force on a small strip $dA$ is $\gamma h\, dA = \gamma y \sin\theta\, dA$. Integrating over the surface,

$$
F = \gamma \sin\theta \int_A y\, dA = \gamma \sin\theta \, \bar{y} A ,
$$

where $\bar{y}$ is the distance to the **centroid** of the area. Since $\bar{h} = \bar{y}\sin\theta$ is the depth of the centroid,

$$
F = \gamma \bar{h} A .
$$

The resultant force equals the pressure at the centroid times the area, whatever the angle of the plane.

The resultant does **not** act at the centroid. Pressure grows with depth, so the lower part of the surface carries more load, and the resultant acts lower down, at the **centre of pressure**. Taking moments about the surface line,

$$
F\, y_{cp} = \gamma \sin\theta \int_A y^2\, dA = \gamma \sin\theta \, I_O ,
$$

where $I_O$ is the second moment of area about the surface line. The parallel-axis theorem, $I_O = I_G + A\bar{y}^2$, gives

$$
y_{cp} = \bar{y} + \frac{I_G}{\bar{y} A},
$$

where $I_G$ is the second moment of area about a horizontal axis through the centroid. The second term is always positive, so the centre of pressure always lies below the centroid. For shapes symmetric about a vertical line, the centre of pressure lies on that line. The table below lists $A$, the centroid and $I_G$ for common gate shapes.

```{list-table} Properties of common plane shapes (axis horizontal, through the centroid)
:header-rows: 1
:widths: 30 20 25 25

* - Shape
  - Area $A$
  - Centroid
  - $I_G$
* - Rectangle, width $b$, height $h$
  - $bh$
  - $h/2$ from top
  - $bh^3/12$
* - Triangle, base $b$, height $h$
  - $bh/2$
  - $h/3$ from base
  - $bh^3/36$
* - Circle, diameter $D$
  - $\pi D^2/4$
  - centre
  - $\pi D^4/64$
* - Semicircle, radius $r$
  - $\pi r^2/2$
  - $4r/3\pi$ from flat side
  - $0.1098\, r^4$
```

**Pressure prism.** For a rectangular surface it is often easier to think of the pressure distribution as a solid, the **pressure prism**. The force equals the volume of the prism (pressure × area), and it acts through the prism's centroid. For a vertical rectangle reaching the free surface, the pressure diagram is a triangle and the force acts at two-thirds of the depth. If the top of the rectangle is below the surface, the diagram is a trapezoid (Figure 2.1).

```{code-cell} python
:tags: [hide-input]
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 11})

ink, ink2, c1 = "#0b0b0b", "#52514e", "#2a78d6"
c2, c3 = "#d9822b", "#3a9a5b"
rho, g = 998.2, 9.81
gamma = rho * g

h1, h2, k = 0.5, 3.0, 0.55           # gate from depth h1 to h2 (m); arrow scale
hbar = (h1 + h2) / 2
hcp = hbar + (h2 - h1) ** 2 / (12 * hbar)

fig, ax = plt.subplots(figsize=(6.5, 4.0))
ax.fill_between([-2.6, 0], -3.6, 0, color="#e8f1fb")             # water
ax.plot([-2.6, 0], [0, 0], color=c1, lw=1.2)                       # free surface
ax.plot([0, 0], [0.4, -3.6], color=ink, lw=2.5)                    # wall
ax.plot([0, 0], [-h1, -h2], color=c2, lw=5, solid_capstyle="butt")  # gate
ax.text(0.1, -h1 + 0.05, "gate", color=c2, va="top")
ax.plot([0, -k * h2], [0, -h2], color=ink2, lw=0.8, ls="--")      # full triangle
xs = np.array([-k * h1, 0, 0, -k * h2])
ax.fill(xs, [-h1, -h1, -h2, -h2], color=c1, alpha=0.25)
for d in np.linspace(h1, h2, 7):
    ax.annotate("", xy=(0, -d), xytext=(-k * d, -d),
                arrowprops=dict(arrowstyle="-|>", color=c1, lw=1.1))
ax.annotate("", xy=(0.02, -hcp), xytext=(1.3, -hcp),
            arrowprops=dict(arrowstyle="<|-", color=ink, lw=2))
ax.text(1.35, -hcp, "  $F$ at centre\n  of pressure", va="center", color=ink)
ax.plot([0, 1.1], [-hbar, -hbar], color=ink2, lw=0.8, ls=":")
ax.text(1.35, -hbar + 0.32, "  centroid", va="center", color=ink2)
ax.text(-2.5, -1.9, r"$p = \gamma h$", color=c1, fontsize=13)
ax.text(-2.5, 0.12, "free surface", color=ink2)
ax.set_xlim(-2.7, 3.2); ax.set_ylim(-3.7, 0.5); ax.axis("off")
plt.tight_layout(); plt.show()
```

*Figure 2.1. Pressure on a vertical gate in a wall, from 0.5 m to 3.0 m depth. The pressure diagram over the gate is a trapezoid, and its resultant acts at the centre of pressure, below the centroid.*

The distance between the centre of pressure and the centroid, $I_G/(\bar{y}A)$, shrinks as the surface goes deeper. At great depth the pressure is nearly uniform over the surface and the two points almost coincide (Figure 2.2). For deep outlet gates in dams it is often enough to apply the force at the centroid; for shallow gates it is not.

```{code-cell} python
:tags: [hide-input]
ybar = np.linspace(1.0, 20, 300)               # centroid depth (m), vertical gates
rect = (2.0 ** 2 / 12) / ybar                  # rectangle 2 m high: I_G/A = h²/12
circ = (2.0 ** 2 / 16) / ybar                  # circle 2 m diameter: I_G/A = D²/16

fig, ax = plt.subplots(figsize=(6.5, 3.6))
ax.plot(ybar, rect, color=c1, lw=2, label="rectangle, 2 m high")
ax.plot(ybar, circ, color=c2, lw=2, label="circle, 2 m diameter")
ax.legend(frameon=False, loc="upper right")
ax.set_xlabel(r"Depth of centroid $\bar{h}$ (m)", color=ink2)
ax.set_ylabel(r"$h_{cp} - \bar{h}$ (m)", color=ink2)
ax.set_xlim(0, 20); ax.set_ylim(0, 0.35)
ax.grid(True, color="#e4e3df"); ax.set_axisbelow(True)
for s in ("top", "right"): ax.spines[s].set_visible(False)
plt.tight_layout(); plt.show()
```

*Figure 2.2. Distance of the centre of pressure below the centroid for two vertical gates. The distance falls in inverse proportion to depth, so for deep gates the force acts almost at the centroid.*

Table 2.1 shows the same effect in numbers for a vertical rectangular gate, 2 m wide and 2 m high, lowered to increasing depths. The force grows in proportion to the centroid depth, while the offset of the centre of pressure falls.

```{code-cell} python
:tags: [remove-input]
from IPython.display import Markdown
b, hg = 2.0, 2.0
rows = []
for top in [0, 1, 2, 5, 10, 20]:
    hb = top + hg / 2
    F = gamma * hb * b * hg
    hc = hb + (b * hg ** 3 / 12) / (hb * b * hg)
    rows.append(f"| {top} | {hb:.1f} | {F/1e3:.1f} | {hc:.3f} | {hc - hb:.3f} |")
Markdown("| Top depth (m) | $\\bar{h}$ (m) | $F$ (kN) | $h_{cp}$ (m) | $h_{cp}-\\bar{h}$ (m) |\n"
         "|---|---|---|---|---|\n" + "\n".join(rows))
```

*Table 2.1. Force and centre of pressure for a vertical gate 2 m wide and 2 m high, water at 20 °C. When the top is at the surface the force acts at two-thirds of the height.*

## 2.5 Force on a curved surface

On a curved surface the pressure acts normal to the surface at every point, so the elementary forces point in different directions and cannot simply be added. Instead we find the **horizontal** and **vertical components** separately, using the equilibrium of the water next to the surface.

**Horizontal component.** The horizontal force on a curved surface equals the force on its **vertical projection**, a plane surface. It acts at the centre of pressure of that projection:

$$
F_H = \gamma \bar{h}_{\text{proj}} A_{\text{proj}} .
$$

**Vertical component.** The vertical force equals the **weight of the water vertically above the surface**, up to the level of the free surface. It acts through the centroid of that volume of water:

$$
F_V = \gamma \mathcal{V}_{\text{above}} .
$$

If the water lies above the surface, $F_V$ acts downward. If the water lies below the surface, as under an overhanging gate, the volume between the surface and the free-surface level contains no water, but the vertical force is still equal to the weight of that **imaginary water**, and it acts upward.

The resultant is $F = \sqrt{F_H^2 + F_V^2}$, inclined at $\tan^{-1}(F_V/F_H)$ to the horizontal. For a **circular** surface every elementary pressure force is normal to the surface and so passes through the centre of the circle. The resultant must therefore pass through the centre too.

This property is the reason for the shape of the **radial gate**, or **Tainter gate**, used on spillway crests and in irrigation canals. The gate face is a circular arc whose centre is the trunnion pin about which the gate rotates. Since the water force passes through the pin, it creates no moment about it. The hoist only needs to overcome the weight of the gate and friction at the pin, which makes large gates easy to operate.

```{code-cell} python
:tags: [hide-input]
Rg = 4.0
phi = np.linspace(0, np.pi / 2, 100)
xa, za = -Rg * np.cos(phi), -Rg * np.sin(phi)
alpha = np.arctan((np.pi / 4) / 0.5)          # tan α = F_V/F_H = (πR²/4)/(R²/2)

fig, ax = plt.subplots(figsize=(6.5, 4.2))
ax.fill(np.r_[-6.5, xa, -6.5], np.r_[0, za, -Rg], color="#e8f1fb")     # real water
ax.fill(np.r_[0, xa], np.r_[0, za], facecolor="none", edgecolor=c1, hatch="///", lw=0)
ax.plot([-6.5, -Rg], [0, 0], color=c1, lw=1.2)
ax.plot([-6.5, 3.0], [-Rg, -Rg], color=ink, lw=2)                    # floor
ax.plot(xa, za, color=c2, lw=4)                                       # gate face
ax.plot([0, -Rg], [0, 0], color=ink2, lw=1); ax.plot([0, 0], [0, -Rg], color=ink2, lw=1)
ax.scatter([0], [0], color=ink, zorder=5, s=40)
ax.text(0.15, 0.15, "trunnion O", color=ink)
xp, zp = -Rg * np.cos(alpha), -Rg * np.sin(alpha)
ax.annotate("", xy=(xp, zp), xytext=(xp - 0.7 * np.cos(alpha), zp - 0.7 * np.sin(alpha)),
            arrowprops=dict(arrowstyle="-|>", color=ink, lw=2.2))
ax.plot([xp, 0], [zp, 0], color=ink, lw=1, ls="--")
ax.text(xp - 0.5, zp - 0.3, "$F$", color=ink, fontsize=13, ha="right")
ax.text(-1.65, -1.3, "imaginary\nwater", color=c1, ha="center")
ax.text(-5.9, -2.6, "water", color=c1)
# force triangle on the dry side
ox, oz, s = 1.0, -3.4, 1.2
ax.annotate("", xy=(ox + s, oz), xytext=(ox, oz), arrowprops=dict(arrowstyle="-|>", color=c1, lw=2))
ax.annotate("", xy=(ox + s, oz + s * np.pi / 2), xytext=(ox + s, oz),
            arrowprops=dict(arrowstyle="-|>", color=c3, lw=2))
ax.annotate("", xy=(ox + s, oz + s * np.pi / 2), xytext=(ox, oz),
            arrowprops=dict(arrowstyle="-|>", color=ink, lw=2))
ax.text(ox + s / 2, oz - 0.15, "$F_H$", color=c1, ha="center", va="top")
ax.text(ox + s + 0.1, oz + 0.9, "$F_V$", color=c3)
ax.set_xlim(-6.6, 3.2); ax.set_ylim(-4.3, 0.6); ax.set_aspect("equal"); ax.axis("off")
plt.tight_layout(); plt.show()
```

*Figure 2.3. A radial gate whose face is a quarter circle centred on the trunnion O, with the water level at the trunnion. The horizontal force acts on the vertical projection; the upward force equals the weight of the hatched imaginary water. The resultant passes through O.*

## 2.6 Buoyancy and stability of floating bodies

**Buoyancy.** A body immersed in water feels a pressure on its lower surface that is larger than the pressure on its upper surface. Applying the vertical-force rule of Section 2.5 to the top and bottom halves of the body, the net upward force equals the weight of the water displaced:

$$
F_B = \gamma \mathcal{V}_d .
$$

This is **Archimedes' principle**. The **buoyancy force** acts through the centroid of the displaced volume, called the **centre of buoyancy** B. A floating body sinks until the displaced water weighs as much as the body. The horizontal components of pressure cancel, so there is no net horizontal force on a body in still water.

**Stability of a submerged body.** A fully submerged body, such as a submarine or a sunken pipe segment, is stable when its centre of gravity G lies **below** B. If the body tilts, the weight and buoyancy form a couple that rotates it back.

**Stability of a floating body.** A floating body, such as a pontoon, barge or caisson being towed to site, can be stable even when G lies above B. When the body heels by a small angle $\phi$, the shape of the displaced volume changes: one side goes deeper and the other rises. The centre of buoyancy moves sideways to a new position B′. The vertical line through B′ crosses the original centre line at the **metacentre** M. If M lies above G, the weight and buoyancy create a **righting moment**; if M lies below G, they create an overturning moment.

The distance from B to M, for small angles, is

$$
BM = \frac{I}{\mathcal{V}_d},
$$

where $I$ is the second moment of area of the **waterplane** (the plan area cut by the water surface) about the axis of tilting, and $\mathcal{V}_d$ is the displaced volume. Measuring heights from the keel K (the bottom of the body), the **metacentric height** is

$$
GM = KB + BM - KG .
$$

The body is stable if $GM > 0$. For small heel angles the righting moment is $W \cdot GM \sin\phi$, where $W$ is the weight. A large $GM$ means a stiff body that rolls quickly back; ships usually have $GM$ between about 0.3 and 1.5 m, while pontoons and work barges are much stiffer.

For a rectangular pontoon of width $B$ and draft $d$, rolling about its long axis, $KB = d/2$ and $BM = B^2/(12d)$. Figure 2.4 shows how the height of the metacentre changes with draft. A lightly loaded pontoon has a high metacentre because its displaced volume is small. As the draft increases, $BM$ falls while $KB$ rises, so $KM$ reaches a minimum at $d = B/\sqrt{6}$. A narrow pontoon can become unstable at moderate draft if its load is stacked high.

```{code-cell} python
:tags: [hide-input]
d = np.linspace(0.3, 2.5, 300)
KG = 1.2
fig, ax = plt.subplots(figsize=(6.5, 3.8))
for Bw, col in [(3.0, c2), (4.0, c1), (5.0, c3)]:
    KM = d / 2 + Bw ** 2 / (12 * d)
    ax.plot(d, KM, color=col, lw=2)
    ax.text(2.52, KM[-1], f"$B$ = {Bw:.0f} m", color=col, va="center")
ax.axhline(KG, color=ink, lw=1.2, ls="--")
ax.text(0.32, KG - 0.1, "centre of gravity, $KG$ = 1.2 m", color=ink, va="top")
ax.set_xlabel("Draft $d$ (m)", color=ink2)
ax.set_ylabel("Metacentre height $KM$ (m)", color=ink2)
ax.set_xlim(0.3, 3.1); ax.set_ylim(0.5, 4.0)
ax.grid(True, color="#e4e3df"); ax.set_axisbelow(True)
for s in ("top", "right"): ax.spines[s].set_visible(False)
plt.tight_layout(); plt.show()
```

*Figure 2.4. Height of the metacentre above the keel for rectangular pontoons of three widths. The pontoon is stable where the curve lies above its centre of gravity; the 3 m pontoon is only marginally stable near 1.2 m draft.*

## 2.7 Relative equilibrium

Water can also be "at rest" relative to a moving container, when every particle moves with the same acceleration and there is no relative motion and no shear. Pressure then varies hydrostatically, but with an effective gravity.

**Linear acceleration.** If an open tank moves horizontally with constant acceleration $a$, the free surface tilts at an angle $\theta$ to the horizontal, with

$$
\tan\theta = \frac{a}{g},
$$

and the surface is higher at the rear. Pressure is still $\gamma h$, where $h$ is the vertical depth below the tilted surface. This matters for tankers and for water in vehicles that brake suddenly.

**Rotating vessel.** If a cylindrical vessel rotates about its vertical axis at angular speed $\omega$ (rad s⁻¹), the water eventually rotates with it as a rigid body. The free surface becomes a **paraboloid**:

$$
z_s(r) = z_0 + \frac{\omega^2 r^2}{2g},
$$

where $z_0$ is the elevation of the surface on the axis and $r$ is the radius. The volume of a paraboloid is half that of its enclosing cylinder, so the surface on the axis falls by as much as the surface at the wall rises. The same principle is used in centrifugal pumps, where rotation raises the pressure at the outer edge of the impeller.

## 2.8 Worked examples

**Example 2.1 — Differential manometer between two pipes.** A mercury manometer ($S_m = 13.6$) connects pipe A to pipe B, both carrying water. The centre of B is 0.5 m above the centre of A, and the mercury reading is $R$ = 0.30 m. Find $p_A - p_B$ and the difference in piezometric head.

```{code-cell} python
rho, g, Sm = 998.2, 9.81, 13.6
gamma = rho * g
R, zA, zB = 0.30, 0.0, 0.5
dp = (Sm - 1) * gamma * R + gamma * (zB - zA)    # Pa
dH = R * (Sm - 1)                                # m of water
print(f"pA - pB = {dp/1e3:.1f} kPa")
print(f"Difference in piezometric head = {dH:.2f} m of water")
```

The pressure at A is about 41.9 kPa higher. Of this, 4.9 kPa simply reflects the 0.5 m difference in elevation; only the 3.78 m difference in piezometric head would drive flow between A and B.

**Example 2.2 — Inclined gate held by a stop.** A rectangular gate 1.5 m wide and 2.0 m long is hinged along its top edge and rests against a stop at its bottom edge. The gate is inclined at 60° to the horizontal, and its hinge is 1.0 m below the free surface. Find the water force, the position of the centre of pressure and the force on the stop. Ignore the weight of the gate.

```{code-cell} python
b, L, theta, h_top = 1.5, 2.0, np.radians(60), 1.0
A = b * L
y_top = h_top / np.sin(theta)            # along the plane, from the surface line
ybar = y_top + L / 2
hbar = ybar * np.sin(theta)
F = gamma * hbar * A
IG = b * L ** 3 / 12
ycp = ybar + IG / (ybar * A)
P_stop = F * (ycp - y_top) / L           # moments about the hinge
print(f"Centroid depth = {hbar:.3f} m, force F = {F/1e3:.1f} kN")
print(f"Centre of pressure: {ycp - ybar:.3f} m below the centroid, "
      f"{ycp - y_top:.3f} m from the hinge (along the gate)")
print(f"Force on the stop = {P_stop/1e3:.1f} kN")
```

The resultant of 54.8 kN acts 1.155 m from the hinge, a little beyond mid-length, because the lower part of the gate carries more pressure. The stop therefore carries slightly more than half the water force: 31.6 kN.

**Example 2.3 — Radial gate.** A radial gate 6 m wide has a face that is a quarter circle of radius 4 m, centred on the trunnion. The water level is at the trunnion and the bottom of the arc rests on the sill directly below it, as in Figure 2.3. Find the components and resultant of the water force, and its direction.

```{code-cell} python
Rg, bw = 4.0, 6.0
FH = gamma * (Rg / 2) * (Rg * bw)          # vertical projection: Rg high, bw wide
FV = gamma * (np.pi * Rg ** 2 / 4) * bw    # weight of imaginary water (quarter circle)
F = np.hypot(FH, FV)
ang = np.degrees(np.arctan2(FV, FH))
print(f"F_H = {FH/1e3:.0f} kN, F_V = {FV/1e3:.0f} kN (upward)")
print(f"Resultant F = {F/1e3:.0f} kN at {ang:.1f}° above the horizontal, through the trunnion")
```

The vertical component is larger than the horizontal one because the imaginary water above the arc (a quarter circle) is larger than the triangular pressure diagram. The resultant of about 875 kN passes through the trunnion, so it produces no moment about the pin; the trunnion bearings, however, must be designed to carry it.

**Example 2.4 — Stability of a pontoon.** A rectangular steel pontoon is 8 m long, 4 m wide and 2 m high. With its load it has a mass of 25 t, and its centre of gravity is 0.9 m above the keel. It floats in fresh water at 20 °C. Find the draft and the metacentric height for rolling about the long axis.

```{code-cell} python
Lp, Bp, m, KG = 8.0, 4.0, 25e3, 0.9
Vd = m / rho                       # displaced volume (m³)
dr = Vd / (Lp * Bp)                # draft
KB = dr / 2
BM = (Lp * Bp ** 3 / 12) / Vd      # waterplane I about the long axis / Vd
GM = KB + BM - KG
print(f"Draft = {dr:.3f} m, KB = {KB:.3f} m, BM = {BM:.3f} m")
print(f"GM = {GM:.2f} m  -> {'stable' if GM > 0 else 'unstable'}")
```

The pontoon floats with a draft of 0.78 m and a metacentric height of 1.19 m, so it is clearly stable. Most of the stability comes from $BM$: a wide waterplane on a shallow draft gives a high metacentre.

## 2.9 Common mistakes

- **Using the depth of the centre of pressure to compute the force.** The force is $\gamma \bar{h} A$, with the depth of the **centroid**. The centre of pressure only tells you where the force acts.
- **Measuring $\bar{y}$ vertically on an inclined surface.** In $y_{cp} = \bar{y} + I_G/(\bar{y}A)$, both $\bar{y}$ and $y_{cp}$ are measured **along the plane** from the surface line.
- **Mixing gauge and absolute pressure.** A Bourdon gauge reads gauge pressure. Add atmospheric pressure before comparing with vapour pressure.
- **Forgetting the water column in a manometer.** In a water–mercury manometer the factor is $(S_m - 1)$, not $S_m$, because the water in the limbs also has weight.
- **Taking $F_V$ downward on an overhanging surface.** When the water is below a curved surface, the vertical force acts **upward**, even though the volume used to compute it contains no water.
- **Using the wrong axis for $I$ in stability.** $BM = I/\mathcal{V}_d$ uses the second moment of the **waterplane** about the tilting axis. Rolling about the long axis uses $LB^3/12$, not $BL^3/12$.

## 2.10 Exercises

1. A Bourdon gauge on a pipe reads 350 kPa. Find the absolute pressure and the pressure head in metres of water. How tall would a piezometer have to be?
2. An oil layer ($S = 0.85$) 1.2 m thick floats on 2.0 m of water in an open tank. Sketch the pressure distribution on the tank wall and find the pressure at the bottom.
3. A circular gate of diameter 1.2 m in the vertical face of a dam has its centre 6 m below the water surface. Find the water force and the distance of the centre of pressure below the centre of the gate.
4. Show that for a vertical rectangle whose top edge lies at the free surface, the centre of pressure is at two-thirds of the height. Repeat for a triangle with its base at the surface.
5. A quarter-circle gate of radius 2 m and width 3 m has water above it, filling the space between the arc and the free surface, which is level with the top of the arc. Find the horizontal and vertical components of the water force.
6. A wooden block ($S = 0.6$), 1.0 m long with a square cross-section 0.4 m × 0.4 m, floats in water with its long axis horizontal. Is it stable when floating with one face horizontal? Find $GM$.

---

*Figures and tables in this chapter are generated by the Python code shown. Text and figures: © Daeha Kim, CC BY-NC 4.0.*
