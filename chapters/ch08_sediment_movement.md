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

# 8. Sediment Movement

```{admonition} Core ideas
:class: tip
1. Flowing water moves sediment when the **bed shear stress** exceeds a critical value. Below it the bed is stable; above it the bed is mobile.
2. The key size property of a grain is its **diameter**, and the key motion property is its **fall velocity**, which balances weight against drag.
3. The **Shields parameter** $\theta$ compares the force of the flow on a grain with the grain's submerged weight. Motion starts when $\theta$ exceeds a critical value $\theta_{cr}$ of about 0.03–0.06 for sand and gravel.
4. Grains roll and hop along the bed as **bed load**, or are carried up into the flow by turbulence as **suspended load**. The **Rouse number** $w_s/(\kappa u_*)$ decides which mode dominates.
5. Transport formulas are **empirical** and uncertain by a factor of two or more. Use them to estimate orders of magnitude and to compare options, not to predict exact quantities.
```

## 8.1 Why sediment matters

Rivers carry not only water but also sand, gravel, silt and clay. A large river can move millions of metric tons of sediment each year. Where the flow slows down, for example in a reservoir or at a river mouth, the sediment settles. Where the flow speeds up, for example around a bridge pier or below a dam, the flow picks up sediment and the bed is lowered.

These processes create some of the most expensive problems in water engineering:

- **Reservoir sedimentation.** Sediment fills the storage volume of a dam. Many reservoirs lose 0.5–1% of their capacity every year, so their useful life is limited.
- **Scour at structures.** Strong vortices around bridge piers and abutments dig holes in the bed. Scour is a leading cause of bridge failure during floods.
- **Channel stability.** An irrigation canal or a flood channel must neither erode its bed and banks nor fill with deposits.
- **Morphology and ecology.** Bars, pools and riffles in a river are made by sediment transport, and they form the habitat of fish and invertebrates.

Sediment transport is a two-phase problem: water and solid grains interact. A complete theory does not exist, so engineers combine simple force balances with experimental results. This chapter presents the classical framework that most practical methods still use.

## 8.2 Sediment properties

**Grain size.** The most important property of a sediment grain is its size. For sand and gravel, size is measured by sieving, and the **sieve diameter** $D$ is the opening of the smallest sieve the grain passes through. Silt and clay are too fine to sieve and are sized by their settling speed. Geologists group sizes in the **Wentworth scale**, in which each class boundary is double the previous one.

```{list-table} Grain-size classes (Wentworth scale)
:header-rows: 1
:widths: 25 25 50

* - Class
  - Diameter $D$ (mm)
  - Typical behavior in rivers
* - Clay
  - below 0.004
  - Cohesive; stays in suspension for days
* - Silt
  - 0.004–0.0625
  - Weakly cohesive; carried in suspension
* - Sand
  - 0.0625–2
  - Moves in suspension and as bed load; forms ripples and dunes
* - Gravel (granules and pebbles)
  - 2–64
  - Moves mainly as bed load during floods
* - Cobbles
  - 64–256
  - Moves only in large floods
* - Boulders
  - above 256
  - Rarely moves; acts as roughness
```

Natural sediments contain a range of sizes. The **grain-size distribution** is a plot of the percentage by mass finer than each size. It is summarized by percentiles: $D_{50}$ is the **median size** (50% of the mass is finer), and $D_{84}$ and $D_{16}$ describe the coarse and fine parts. The ratio $\sigma_g = \sqrt{D_{84}/D_{16}}$ is the **geometric standard deviation**; a sediment with $\sigma_g < 1.3$ is called well sorted (nearly uniform). In this chapter $D$ means $D_{50}$ unless stated.

**Density.** Most natural sand and gravel is quartz or feldspar, with **density** $\rho_s \approx 2650$ kg m⁻³. The **relative density** of the grains is

$$
s = \frac{\rho_s}{\rho} \approx 2.65,
$$

and the grains' submerged weight per unit volume is $(\rho_s - \rho) g = (s-1)\rho g$. The factor $(s-1) \approx 1.65$ appears in almost every sediment formula.

**Porosity.** A deposited bed is not solid: water fills the spaces between grains. The **porosity** $\lambda_p$ is the volume of voids divided by the total volume. For sand and gravel beds $\lambda_p \approx 0.3$–0.4. When we convert a transported volume of grains into a volume of deposit, we must divide by $(1-\lambda_p)$.

**Angle of repose.** A dry pile of loose grains stands at a maximum slope called the **angle of repose** $\phi_r$. It is about 30° for rounded sand and up to about 40° for angular gravel and rock. The angle of repose controls the stability of channel banks (Section 8.9).

## 8.3 Fall velocity

A grain released in still water accelerates until the drag force balances its submerged weight. The constant speed it then reaches is the **fall velocity** (or settling velocity) $w_s$. For a sphere of diameter $D$, with drag coefficient $C_D$ (Chapter 4):

$$
(\rho_s - \rho)\, g\, \frac{\pi D^3}{6}
= C_D\, \frac{\pi D^2}{4}\, \frac{\rho w_s^2}{2}.
$$

Solving for $w_s$:

$$
w_s = \sqrt{\frac{4 (s-1) g D}{3 C_D}}.
$$

The drag coefficient depends on the **particle Reynolds number** $Re_p = w_s D/\nu$. For very small grains ($Re_p < 1$, roughly $D < 0.1$ mm) the flow around the grain is laminar and $C_D = 24/Re_p$. Substituting gives **Stokes' law**:

$$
w_s = \frac{(s-1)\, g\, D^2}{18\, \nu}.
$$

Stokes' law says that fall velocity grows with $D^2$. For larger grains the drag is dominated by the pressure difference in the wake, $C_D$ becomes nearly constant, and $w_s$ grows only with $\sqrt{D}$. Natural grains are not spheres, so they fall more slowly than spheres of the same sieve diameter. A convenient formula that covers all sizes is that of **Ferguson and Church (2004)**:

$$
w_s = \frac{(s-1)\, g\, D^2}{C_1 \nu + \sqrt{0.75\, C_2\, (s-1)\, g\, D^3}},
$$

with $C_1 = 18$ and $C_2 = 1.0$ for natural sand and gravel measured by sieving. For small $D$ the first term in the denominator dominates and the formula becomes Stokes' law. For large $D$ the second term dominates and the formula gives the constant-drag result with $C_D = C_2$.

```{code-cell} python
:tags: [hide-input]
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 11})

ink, ink2, c1 = "#0b0b0b", "#52514e", "#2a78d6"
c2, c3 = "#d9822b", "#3a9a5b"
grid = "#e4e3df"

rho, rho_s, g, nu = 998.2, 2650.0, 9.81, 1.004e-6
s = rho_s / rho

def w_stokes(D):
    return (s - 1) * g * D**2 / (18 * nu)

def w_fc(D, C1=18.0, C2=1.0):
    """Fall velocity, Ferguson and Church (2004)."""
    return (s - 1) * g * D**2 / (C1 * nu + np.sqrt(0.75 * C2 * (s - 1) * g * D**3))

def style(ax):
    ax.grid(True, which="major", color=grid); ax.set_axisbelow(True)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)

D = np.logspace(-5, -1, 300)                 # m
fig, ax = plt.subplots(figsize=(6.5, 4.0))
ax.loglog(D * 1e3, w_fc(D), color=c1, lw=2.2)
ax.loglog(D * 1e3, w_stokes(D), color=c2, lw=1.8, ls="--")
ax.text(2.5, 0.012, "Ferguson and Church\n(natural grains)", color=c1, fontsize=11)
ax.text(0.25, 3, "Stokes' law", color=c2, fontsize=11)
for b in [0.0625, 2.0]:
    ax.axvline(b, color=ink2, lw=0.8, ls=":")
for x, lab in [(0.012, "silt"), (0.33, "sand"), (9, "gravel")]:
    ax.text(x, 3e-5, lab, color=ink2, fontsize=11, ha="center")
ax.set_xlim(1e-2, 100); ax.set_ylim(1e-5, 10)
ax.set_xlabel("Grain diameter $D$ (mm)", color=ink2)
ax.set_ylabel("Fall velocity $w_s$ (m s⁻¹)", color=ink2)
style(ax)
plt.tight_layout(); plt.show()
```

*Figure 8.1. Fall velocity of quartz grains in water at 20 °C. Stokes' law is accurate only for silt and very fine sand; for coarser grains it overestimates the fall velocity by orders of magnitude.*

Fall velocity is a better measure of how a fine grain behaves in water than its diameter. A silt particle of 0.01 mm falls at about 0.1 mm s⁻¹, or 8 m per day, so even slow turbulence keeps it suspended. A 1 mm sand grain falls at about 0.13 m s⁻¹ and settles quickly once the flow slows down.

## 8.4 Threshold of motion

**Bed shear stress.** The flow drags on the bed with the **bed shear stress** $\tau_0$. In steady uniform flow the weight component of the water along the slope balances the friction on the boundary (Chapter 6), which gives

$$
\tau_0 = \rho g R S,
$$

where $R$ is the hydraulic radius and $S$ the bed slope. In a wide channel $R \approx h$, the flow depth. It is convenient to express the shear stress as a velocity, the **shear velocity**:

$$
u_* = \sqrt{\frac{\tau_0}{\rho}} = \sqrt{g R S}.
$$

$u_*$ is not a real velocity of the water. It is a measure of the intensity of turbulence near the bed, and it is typically 5–10% of the mean velocity $V$.

**The Shields parameter.** A grain on the bed is pushed by a fluid force proportional to $\tau_0 D^2$ and held back by its submerged weight, proportional to $(\rho_s - \rho) g D^3$. The ratio of these two is the **Shields parameter** (or Shields stress):

$$
\theta = \frac{\tau_0}{(\rho_s - \rho)\, g\, D}.
$$

The grain starts to move when $\theta$ exceeds a **critical Shields parameter** $\theta_{cr}$. In his experiments in 1936, Shields found that $\theta_{cr}$ is not a constant but depends on how the grain size compares with the thin viscous layer next to the bed. This is described by the **dimensionless grain size**

$$
D_* = D \left[ \frac{(s-1)\, g}{\nu^2} \right]^{1/3}.
$$

$D_*$ depends only on the grain and the fluid, not on the flow, so it can be calculated before we know $\tau_0$. For quartz in water at 20 °C, $D_* \approx 25$ per millimeter of grain diameter. The relation between $\theta_{cr}$ and $D_*$ is the **Shields curve**. A convenient explicit fit is that of **Soulsby and Whitehouse (1997)**:

$$
\theta_{cr} = \frac{0.30}{1 + 1.2 D_*} + 0.055 \left[ 1 - e^{-0.020 D_*} \right].
$$

```{code-cell} python
:tags: [hide-input]
def Dstar(D):
    return D * ((s - 1) * g / nu**2) ** (1 / 3)

def theta_cr(Ds):
    """Critical Shields parameter, Soulsby and Whitehouse (1997)."""
    return 0.30 / (1 + 1.2 * Ds) + 0.055 * (1 - np.exp(-0.020 * Ds))

Ds = np.logspace(0, 4, 300)
k = Dstar(1e-3)                              # D* per mm of diameter
fig, ax = plt.subplots(figsize=(6.5, 4.0))
ax.loglog(Ds, theta_cr(Ds), color=c1, lw=2.4)
ax.fill_between(Ds, theta_cr(Ds), 1, color=c1, alpha=0.08)
ax.text(30, 0.25, "grains move", color=c1, fontsize=12)
ax.text(30, 0.016, "bed stable", color=ink2, fontsize=12)
ax.set_xlim(1, 1e4); ax.set_ylim(0.01, 0.5)
ax.set_xlabel("Dimensionless grain size $D_*$", color=ink2)
ax.set_ylabel(r"Critical Shields parameter $\theta_{cr}$", color=ink2)
top = ax.secondary_xaxis("top", functions=(lambda x: x / k, lambda d: d * k))
top.set_xlabel("$D$ for quartz in water at 20 °C (mm)", color=ink2)
style(ax)
plt.tight_layout(); plt.show()
```

*Figure 8.2. The Shields curve for the threshold of motion, from the fit of Soulsby and Whitehouse (1997). Above the curve the flow moves the grains. For gravel ($D_* > 300$), $\theta_{cr}$ is nearly constant at about 0.055.*

Three features of the curve are worth remembering. First, $\theta_{cr}$ has a minimum of about 0.03 for medium to coarse sand ($D \approx 0.5$–1 mm). Second, for fine grains $\theta_{cr}$ rises, because small grains hide inside the viscous layer near the bed (very fine silt and clay are also cohesive, which the curve does not include). Third, for gravel $\theta_{cr}$ is nearly constant, so the **critical shear stress** $\tau_{cr} = \theta_{cr} (\rho_s - \rho) g D$ is simply proportional to $D$. A useful rule for gravel is $\tau_{cr}$ (Pa) $\approx 0.9 \times D$ (mm).

```{code-cell} python
:tags: [remove-input]
from IPython.display import Markdown
rows = []
for Dmm in [0.0625, 0.125, 0.25, 0.5, 1, 2, 8, 32, 128]:
    Dm = Dmm / 1e3
    th = theta_cr(Dstar(Dm))
    tc = th * (rho_s - rho) * g * Dm
    rows.append(f"| {Dmm:g} | {w_fc(Dm)*1e3:.3g} | {Dstar(Dm):.3g} | {th:.3f} | {tc:.3g} | {np.sqrt(tc/rho):.3f} |")
Markdown("| $D$ (mm) | $w_s$ (mm s⁻¹) | $D_*$ | $\\theta_{cr}$ | $\\tau_{cr}$ (Pa) | $u_{*cr}$ (m s⁻¹) |\n"
         "|---|---|---|---|---|---|\n" + "\n".join(rows))
```

*Table 8.1. Fall velocity (Ferguson and Church, 2004) and threshold of motion (Soulsby and Whitehouse, 1997) for quartz grains in water at 20 °C. The grain diameter changes by a factor of 2000 across the table, but $\theta_{cr}$ stays between 0.03 and 0.11.*

## 8.5 Bed forms

When sand starts to move, the bed rarely stays flat. It organizes itself into regular waves called **bed forms**. Bed forms matter because they change the flow resistance of the channel (the Manning $n$ of a sand-bed river can double when dunes form) and because they control how sediment moves.

As the flow becomes stronger, a sand bed passes through a sequence of forms, which is linked to the Froude number $Fr = V/\sqrt{g h}$ (Chapter 6):

- **Ripples.** Small, triangular forms a few centimeters high and less than about 0.6 m long. They form only in fine and medium sand ($D$ less than about 0.6 mm) at low shear stress. Their size depends on the grain size, not on the depth.
- **Dunes.** Larger triangular forms with a gentle upstream face and a steep downstream face, with flow separation behind the crest. Dune height scales with flow depth (often 10–30% of $h$). Both ripples and dunes move slowly downstream. They belong to the **lower regime**, with $Fr$ well below 1.
- **Plane bed.** At higher velocity the dunes are washed out and the bed becomes flat again, with intense transport. Resistance drops sharply. This is the **upper regime**.
- **Antidunes.** When $Fr$ approaches or exceeds 1, the bed forms waves that are in phase with surface waves. They can move upstream, and they often break.

In a field river during a flood, the bed may change from dunes to plane bed as discharge rises. The depth may then fall even though the discharge increases, because the resistance has dropped. This is one reason why stage–discharge curves of sand-bed rivers show loops and jumps.

## 8.6 Bed load

**Bed load** is the part of the transport in which grains stay close to the bed. They roll, slide or make short jumps (**saltation**), and their weight is carried by contact with the bed. Bed load dominates in gravel rivers.

The transport rate is expressed as a volume of grains per unit time per unit width of channel, $q_b$ (m² s⁻¹). It is made dimensionless with the **Einstein parameter**:

$$
\Phi = \frac{q_b}{\sqrt{(s-1)\, g\, D^3}}.
$$

Most bed-load formulas relate $\Phi$ to the excess of $\theta$ above its critical value. The best-known is that of **Meyer-Peter and Müller (1948)**, from flume experiments with uniform sand and gravel:

$$
\Phi = 8\, (\theta - 0.047)^{3/2}, \quad \theta > 0.047.
$$

Note that this formula has its own threshold, 0.047, which is slightly lower than the Shields value for gravel in Figure 8.2. The threshold of motion is not a sharp line: a few grains move well below it, and the exact value depends on how "start of motion" is defined. A re-analysis of the original data by **Wong and Parker (2006)** gave $\Phi = 3.97 (\theta - 0.0495)^{3/2}$, which predicts about half the transport. The difference between the two curves is a reminder of the uncertainty in all such formulas.

```{code-cell} python
:tags: [hide-input]
th = np.linspace(0, 0.3, 400)
Phi_mpm = 8 * np.clip(th - 0.047, 0, None) ** 1.5
Phi_wp = 3.97 * np.clip(th - 0.0495, 0, None) ** 1.5

fig, ax = plt.subplots(figsize=(6.5, 4.0))
ax.plot(th, Phi_mpm, color=c1, lw=2.4)
ax.plot(th, Phi_wp, color=c2, lw=1.8, ls="--")
ax.text(0.095, 0.62, "Meyer-Peter and\nMüller (1948)", color=c1, fontsize=11)
ax.text(0.205, 0.04, "Wong and\nParker (2006)", color=c2, fontsize=11)
ax.axvline(0.047, color=ink2, lw=0.8, ls=":")
ax.text(0.052, 0.95, "threshold\n$\\theta$ = 0.047", color=ink2, fontsize=11, va="top")
ax.set_xlim(0, 0.3); ax.set_ylim(0, 1.0)
ax.set_xlabel(r"Shields parameter $\theta$", color=ink2)
ax.set_ylabel(r"Einstein parameter $\Phi$", color=ink2)
style(ax)
plt.tight_layout(); plt.show()
```

*Figure 8.3. Dimensionless bed-load rate against Shields parameter. Transport is zero below the threshold and grows quickly above it: doubling $\theta$ from 0.1 to 0.2 increases $\Phi$ by a factor of about five.*

The steep rise of the curve has a practical meaning. Most bed load is moved in a few days of large floods each year, not by the everyday flow. The Meyer-Peter and Müller formula should be applied with the shear stress acting on the grains themselves. On a bed with dunes, part of $\tau_0$ is spent as form drag on the dunes and does not move grains, so using the total $\tau_0$ overestimates bed load.

## 8.7 Suspended load

When the shear velocity becomes comparable with the fall velocity, turbulent eddies can lift grains away from the bed faster than they settle. The grains are then carried in the body of the flow as **suspended load**. Most of the sediment in large sand-bed rivers moves this way.

In steady conditions the downward settling flux balances the upward turbulent mixing flux at every height $z$ above the bed:

$$
w_s\, c + \varepsilon_s \frac{dc}{dz} = 0,
$$

where $c$ is the volumetric concentration and $\varepsilon_s$ the **eddy diffusivity** of sediment. If the diffusivity has the parabolic form $\varepsilon_s = \kappa u_* z (1 - z/h)$, the equation can be integrated to give the **Rouse profile** (Rouse, 1937):

$$
\frac{c}{c_a} = \left[ \frac{h - z}{z} \cdot \frac{a}{h - a} \right]^{P},
$$

where $c_a$ is a reference concentration at a small height $a$ above the bed (often $a = 0.05h$), and

$$
P = \frac{w_s}{\kappa\, u_*}
$$

is the **Rouse number**, with **von Kármán constant** $\kappa = 0.4$. The Rouse number compares the settling tendency of the grains ($w_s$) with the mixing power of the turbulence ($u_*$).

```{code-cell} python
:tags: [hide-input]
a = 0.05
z = np.linspace(a, 0.999, 400)                # z/h
fig, ax = plt.subplots(figsize=(6.5, 4.2))
for P, zl, col in zip([0.2, 0.5, 1.0, 2.5], [0.85, 0.65, 0.5, 0.3], [c1, c3, c2, ink]):
    c = (((1 - z) / z) * (a / (1 - a))) ** P
    ax.semilogx(c, z, color=col, lw=2)
    cl = (((1 - zl) / zl) * (a / (1 - a))) ** P
    if P < 2:   # label to the left of the curve
        ax.text(cl / 1.3, zl, f"P = {P}", color=col, fontsize=11, ha="right",
                va="center", bbox=dict(fc="white", ec="none", pad=0.5))
    else:       # label above and to the right of the curve
        ax.text(cl * 1.6, zl, f"P = {P}", color=col, fontsize=11, ha="left", va="bottom")
ax.set_xlim(1e-4, 1.2); ax.set_ylim(0, 1)
ax.set_xlabel("Relative concentration $c/c_a$ (log scale)", color=ink2)
ax.set_ylabel("Height above bed $z/h$", color=ink2)
style(ax)
plt.tight_layout(); plt.show()
```

*Figure 8.4. Rouse concentration profiles with reference height $a = 0.05h$ (note the logarithmic concentration axis). For a small Rouse number the sediment is mixed through the whole depth; for $P$ above about 2.5 almost all of it stays near the bed.*

As a rough guide, grains with $P > 2.5$ move mainly as bed load; for $0.8 < P < 2.5$ a large part is in suspension; and for $P < 0.8$ the grains are spread over the whole depth. Very fine material (silt and clay) that is not found in significant amounts in the bed is called **wash load**. Its amount depends on the supply from the catchment, not on the transport capacity of the river, so it cannot be predicted from the flow alone.

## 8.8 Total load

For many engineering purposes we need the **total load**, bed load plus bed-material suspended load. Total-load formulas skip the split between the two modes and relate the total rate directly to the flow. A simple and widely used one for sand-bed rivers is that of **Engelund and Hansen (1967)**. In dimensional form, with $q_t$ the total volumetric transport rate per unit width:

$$
q_t = 0.05\, V^2\, \theta^{3/2} \sqrt{\frac{D}{(s-1)\, g}}.
$$

Here $V$ is the mean velocity and $\theta$ is based on the total bed shear stress $\rho g R S$. The formula has no threshold, so it should be used only when the bed is clearly mobile, and it was developed for sand with $D > 0.15$ mm. Like all transport formulas, its predictions can be wrong by a factor of two or more. Where possible, engineers check them against measured data from the river of interest.

The **sediment concentration** of a river is often quoted in mg L⁻¹ (equal to g m⁻³). A volumetric transport rate converts to a mean concentration as $C = \rho_s q_t / q$, where $q$ is the water discharge per unit width.

## 8.9 Engineering applications

**Stable channels and riprap.** A channel that must not erode is designed so that the bed shear stress stays below the critical value of the bed material. If the natural soil is too fine, the bed and banks are protected with large stones called **riprap**. A first estimate of the required stone size follows directly from the Shields criterion:

$$
D \geq \frac{F_s\, \tau_0}{\theta_{cr}\, (\rho_s - \rho)\, g},
$$

with $\theta_{cr} \approx 0.047$ for large stones and a safety factor $F_s$ of about 1.5. On a bank sloping at angle $\alpha$, gravity helps to pull the stone down the slope, so it moves more easily. Lane (1955) expressed this as a reduction of the critical shear stress by the factor

$$
K = \sqrt{1 - \frac{\sin^2 \alpha}{\sin^2 \phi_r}},
$$

where $\phi_r$ is the angle of repose. At the same time the shear stress on the banks of a trapezoidal channel is lower than on the bed, about $0.75\,\rho g h S$. A bank steeper than $\phi_r$ is unstable whatever the flow. Detailed design manuals add factors for bends, waves and stone shape.

**Local scour.** Around a bridge pier the flow forms a **horseshoe vortex** that scours a hole at the upstream face. The equilibrium depth of local scour in sand can be of the order of twice the pier width. Pier foundations must be placed below the deepest expected scour, including the general lowering of the bed during floods.

**Reservoir sedimentation.** A reservoir traps most of the sediment that enters it. The fraction trapped is the **trap efficiency**. Brune (1953) related it to the ratio of reservoir capacity to mean annual inflow: large reservoirs trap almost all the sediment, while small ones on large rivers trap less. Engineers use such estimates to predict reservoir life, and they plan flushing, sluicing or dredging to extend it.

## 8.10 Worked examples

**Example 8.1 — Settling basin.** Water containing fine sand with $D = 0.15$ mm enters a settling basin 3.0 m deep. How long does a grain take to fall from the surface to the floor? Compare the result with Stokes' law.

```{code-cell} python
D, depth = 0.15e-3, 3.0
ws = w_fc(D)
ws_st = w_stokes(D)
print(f"Ferguson-Church: w_s = {ws*1e3:.1f} mm/s, settling time = {depth/ws:.0f} s")
print(f"Stokes' law:     w_s = {ws_st*1e3:.1f} mm/s, settling time = {depth/ws_st:.0f} s")
print(f"Particle Reynolds number = {ws*D/nu:.1f}")
```

The grain takes about three minutes to reach the floor, so the water must stay in the basin at least this long (in practice longer, because turbulence slows settling). Stokes' law overestimates the fall velocity by more than a third here because the particle Reynolds number is above 1.

**Example 8.2 — Will the gravel move?** A wide gravel-bed river has slope $S = 0.002$ and median bed size $D = 30$ mm. Find the depth at which the bed starts to move, and check whether it moves at a flood depth of 2.0 m.

```{code-cell} python
D, S = 0.030, 0.002
th_c = theta_cr(Dstar(D))
tau_c = th_c * (rho_s - rho) * g * D
h_c = tau_c / (rho * g * S)                     # wide channel: R ≈ h
h = 2.0
tau0 = rho * g * h * S
theta = tau0 / ((rho_s - rho) * g * D)
print(f"D* = {Dstar(D):.0f}, theta_cr = {th_c:.3f}, tau_cr = {tau_c:.1f} Pa")
print(f"Threshold depth = {h_c:.2f} m")
print(f"At h = {h} m: tau0 = {tau0:.1f} Pa, theta = {theta:.3f}, theta/theta_cr = {theta/th_c:.2f}")
```

The bed starts to move when the depth exceeds about 1.4 m. At 2.0 m the Shields parameter is about 1.5 times the critical value, so the bed is mobile but not far above the threshold, which is typical of gravel rivers in flood.

**Example 8.3 — Bed load during the flood.** For the river in Example 8.2, at $h = 2.0$ m, estimate the bed-load transport rate with the Meyer-Peter and Müller formula. The active channel is 40 m wide. Give the result in kg s⁻¹ and in metric tons per day.

```{code-cell} python
B = 40.0
Phi = 8 * max(theta - 0.047, 0) ** 1.5
qb = Phi * np.sqrt((s - 1) * g * D**3)          # m^2/s, volume of grains
Gs = qb * rho_s * B                              # kg/s
print(f"Phi = {Phi:.4f}, q_b = {qb*1e3:.2f} x 10^-3 m^2/s")
print(f"Bed load = {Gs:.0f} kg/s  =  {Gs*86400/1e3:,.0f} metric tons/day")
```

About 110 kg of gravel passes every second, or about 9,400 metric tons per day. If the flood lasted three days at this level, the deposit volume would be about $3 \times 9400 / (2.65 \times 0.65) \approx 16{,}000$ m³ (using a porosity of 0.35). The Wong and Parker version of the formula would give about half this amount.

**Example 8.4 — Total load in a sand-bed river.** A wide sand-bed river with $D = 0.3$ mm flows 3.0 m deep on a slope of 0.0002 with mean velocity 1.2 m s⁻¹. Estimate the total load per meter width with the Engelund–Hansen formula, and the mean sediment concentration in mg L⁻¹.

```{code-cell} python
D, h, S, V = 0.3e-3, 3.0, 0.0002, 1.2
tau0 = rho * g * h * S
theta = tau0 / ((rho_s - rho) * g * D)
qt = 0.05 * V**2 * theta**1.5 * np.sqrt(D / ((s - 1) * g))
C = rho_s * qt / (V * h) * 1e3                  # kg/m^3 -> mg/L
print(f"tau0 = {tau0:.2f} Pa, theta = {theta:.2f} (theta_cr = {theta_cr(Dstar(D)):.3f})")
print(f"q_t = {qt*1e3:.3f} x 10^-3 m^2/s = {qt*rho_s:.2f} kg/s per meter width")
print(f"Mean concentration = {C:.0f} mg/L")
```

The Shields parameter is more than 30 times the critical value, so the bed is fully mobile and the formula applies. Each meter of the river carries about 1.1 kg of sand per second, at a mean concentration of about 300 mg L⁻¹.

**Example 8.5 — How is the sand carried?** In the river of Example 8.4, find the Rouse number of the 0.3 mm sand. With $a = 0.05h$, what is the concentration at mid-depth relative to $c_a$?

```{code-cell} python
u_star = np.sqrt(g * h * S)
P = w_fc(D) / (0.4 * u_star)
a, zm = 0.05 * h, 0.5 * h
ratio = (((h - zm) / zm) * (a / (h - a))) ** P
print(f"u* = {u_star:.3f} m/s, w_s = {w_fc(D)*1e3:.1f} mm/s, P = {P:.2f}")
print(f"c/c_a at mid-depth = {ratio:.3f}")
```

The Rouse number is close to 1, so a large part of the sand moves in suspension, but the concentration at mid-depth is only about 2% of that near the bed. Samples taken only near the water surface would seriously underestimate the load.

**Example 8.6 — Riprap for a flood channel.** A trapezoidal flood channel has depth 2.0 m, slope 0.003 and side slopes of 1 vertical to 2 horizontal. Size the riprap for the bed and for the banks, with $\theta_{cr} = 0.047$, $F_s = 1.5$ and angle of repose 40°.

```{code-cell} python
h, S, Fs, thc = 2.0, 0.003, 1.5, 0.047
tau_bed = rho * g * h * S
tau_bank = 0.75 * tau_bed
alpha, phi_r = np.arctan(1 / 2), np.radians(40)
K = np.sqrt(1 - np.sin(alpha)**2 / np.sin(phi_r)**2)
D_bed = Fs * tau_bed / (thc * (rho_s - rho) * g)
D_bank = Fs * tau_bank / (K * thc * (rho_s - rho) * g)
print(f"tau_bed = {tau_bed:.1f} Pa, tau_bank = {tau_bank:.1f} Pa, K = {K:.2f}")
print(f"Required D50: bed {D_bed*1e3:.0f} mm, banks {D_bank*1e3:.0f} mm")
```

The bank stones need to be about as large as the bed stones even though the shear stress on the banks is lower, because gravity helps to pull them down the slope. A designer would specify one riprap grading with $D_{50}$ of about 120 mm for the whole section.

## 8.11 Common mistakes

- **Using $\rho_s$ instead of $\rho_s - \rho$.** The grain's resistance to motion is its **submerged** weight. Using $\rho_s$ makes the Shields parameter about 40% too small.
- **Applying Stokes' law to sand.** Stokes' law is valid only for $Re_p < 1$, roughly $D < 0.1$ mm. For coarser grains use a general formula such as Ferguson and Church.
- **Confusing $u_*$ with $V$.** The shear velocity is only about 5–10% of the mean velocity. Comparing $V$ with $w_s$ instead of $u_*$ with $w_s$ gives the wrong transport mode.
- **Forgetting porosity.** A volume of transported grains occupies $1/(1-\lambda_p)$ times more space once deposited, about 1.5 times for sand.
- **Mixing thresholds and formulas.** Each transport formula was fitted with its own threshold. Use 0.047 with Meyer-Peter and Müller, not a value read from a different Shields curve.
- **Trusting transport predictions too much.** Different formulas for the same river often differ by a factor of two to ten. Always state the formula used and, where possible, compare with measurements.

## 8.12 Exercises

1. Show that the Shields parameter, the Einstein parameter and $D_*$ are all dimensionless.
2. Find the grain diameter for which the Stokes fall velocity is twice the Ferguson and Church value. What is the particle Reynolds number at that size?
3. A wide channel lined with gravel of $D = 12$ mm carries water 0.8 m deep. What is the steepest slope at which the bed remains stable?
4. Using the Meyer-Peter and Müller formula, show that the bed-load rate per unit width $q_b$ for $\theta \gg 0.047$ is proportional to $\tau_0^{3/2}$, independent of grain size. Explain why this is physically reasonable.
5. A river carries 0.1 mm sand with shear velocity 0.05 m s⁻¹. Calculate the Rouse number and plot the concentration profile $c/c_a$ for $a = 0.05h$. At what height is $c/c_a = 0.5$?
6. A reservoir has a storage volume of 50 × 10⁶ m³ and traps 400,000 metric tons of sediment per year. With a deposit porosity of 0.4 and $\rho_s = 2650$ kg m⁻³, how many years does it take to lose half its storage?

---

*Figures and tables in this chapter are generated by the Python code shown. Text and figures: © Daeha Kim, CC BY-NC 4.0.*
