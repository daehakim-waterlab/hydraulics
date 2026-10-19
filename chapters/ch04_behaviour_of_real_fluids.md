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

# 4. Behaviour of real fluids

```{admonition} Core ideas
:class: tip
1. Real flows are either **laminar** (smooth layers) or **turbulent** (random eddies). The **Reynolds number** $Re = VD/\nu$, the ratio of inertial to viscous forces, decides which.
2. Laminar pipe flow can be solved exactly: the velocity profile is a **parabola** and the friction factor is $f = 64/Re$.
3. Turbulence mixes momentum across the flow. The mean velocity profile is much **flatter** than in laminar flow, and friction is much larger.
4. Near a wall, viscosity acts only in a thin **boundary layer**. Inside it, the velocity follows a **logarithmic law** scaled by the **shear velocity** $u_* = \sqrt{\tau_0/\rho}$.
5. A wall is **hydraulically smooth** when its roughness is buried in the viscous sublayer and **hydraulically rough** when the roughness sticks out of it.
6. When the boundary layer meets rising pressure it can **separate**. Separation creates a wake, form drag and most of the energy loss at fittings and bridge piers.
```

## 4.1 Laminar and turbulent flow

In 1883 Osborne Reynolds let water flow from a large tank through a horizontal glass tube and injected a thin filament of dye at the tube entrance. At low velocity the dye stayed as a straight, sharp line along the whole tube. At higher velocity the line began to wave, and above a certain velocity it broke up suddenly and coloured the whole tube within a short distance.

Reynolds had found the two basic types of real flow. In **laminar flow** the fluid moves in smooth layers (*lamina* is Latin for a thin sheet), and neighbouring layers exchange momentum only through molecular viscosity. In **turbulent flow** the fluid moves in irregular eddies of many sizes, and the eddies carry momentum, heat and dye across the flow far faster than viscosity can.

Reynolds showed that the change does not depend on velocity alone, but on the dimensionless group now called the **Reynolds number**:

$$
Re = \frac{\rho V D}{\mu} = \frac{V D}{\nu},
$$

where $V$ is the mean velocity, $D$ the pipe diameter and $\nu$ the kinematic viscosity (Chapter 1). The Reynolds number measures the ratio of **inertial forces** to **viscous forces**. When $Re$ is small, viscosity damps any disturbance and the flow stays laminar. When $Re$ is large, disturbances grow and the flow becomes turbulent.

For flow in circular pipes, engineers use these limits:

- $Re < 2000$: laminar;
- $2000 < Re < 4000$: **transitional** (the flow switches between the two states);
- $Re > 4000$: turbulent.

The lower **critical Reynolds number** of about 2000 is reliable. The upper value is not: in very smooth pipes with a quiet inflow, laminar flow has been kept to $Re \approx 10^5$, but any small disturbance then makes it turbulent.

**Open channels.** For a channel, the natural length scale is the **hydraulic radius** $R = A/P$, the flow area divided by the wetted perimeter, and $Re = VR/\nu$. A full circular pipe has $R = D/4$, so the pipe limit of 2000 becomes about **500** in terms of $R$. Open-channel flow is laminar below about 500; in practice it is turbulent in almost every canal and river.

```{code-cell} python
:tags: [remove-input]
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import Markdown

plt.rcParams.update({"font.size": 11})
ink, ink2, c1 = "#0b0b0b", "#52514e", "#2a78d6"
c2, c3 = "#d9822b", "#3a9a5b"
nu = 1.004e-6   # m² s⁻¹, water at 20 °C

def tidy(ax):
    ax.grid(True, color="#e4e3df"); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)

# (name, length-scale label, length scale L in m, V in m/s, kind, laminar limit)
flows = [("Groundwater in sand", "grain $d$ = 0.5 mm", 0.5e-3, 1e-4, "grain", 1),
         ("Rain sheet flow on a road", "$R$ ≈ depth = 2 mm", 2e-3, 0.1, "channel", 500),
         ("Household water pipe", "$D$ = 15 mm", 0.015, 1.0, "pipe", 2000),
         ("Irrigation canal", "$R$ = 0.5 m", 0.5, 0.8, "channel", 500),
         ("City water main", "$D$ = 0.5 m", 0.5, 1.5, "pipe", 2000),
         ("Large river in flood", "$R$ = 3 m", 3.0, 2.5, "channel", 500)]
rows = []
for name, lab, L, V, kind, lim in flows:
    Re = V * L / nu
    regime = "laminar" if Re < lim else "turbulent"
    e = int(np.floor(np.log10(Re)))
    Re_s = f"{Re:.2f}" if Re < 1 else f"${Re / 10**e:.1f}\\times10^{{{e}}}$"
    rows.append(f"| {name} | {lab} | {V:g} | {Re_s} | {regime} |")
Markdown("| Flow | Length scale | V (m s⁻¹) | Re | Regime |\n"
         "|---|---|---|---|---|\n" + "\n".join(rows))
```

*Table 4.1. Typical Reynolds numbers for flows of water at 20 °C. The length scale is the pipe diameter, the hydraulic radius or the grain size, as stated. Only very slow or very thin flows are laminar.*

Table 4.1 shows that laminar flow of water is rare in civil engineering. It occurs in groundwater, in thin films and in small laboratory tubes, but the flow in water mains, sewers, canals and rivers is turbulent. We still study laminar flow carefully, because it can be solved exactly and because a thin laminar layer survives next to every wall even in turbulent flow.

## 4.2 Laminar flow in a circular pipe

Consider steady, **fully developed** laminar flow in a straight pipe of radius $r_0 = D/2$. Fully developed means that the velocity profile no longer changes along the pipe. Take a coaxial cylinder of fluid with radius $r$ and length $\Delta x$. The difference in piezometric head pushes it forward, and shear stress $\tau$ on its curved surface holds it back. Because the fluid does not accelerate, these forces balance:

$$
\gamma\, S_f\, \pi r^2 \Delta x = \tau \, 2\pi r \Delta x,
$$

where $S_f = -\,d(p/\gamma + z)/dx$ is the slope of the hydraulic grade line, also called the **friction slope**. Therefore

$$
\tau = \frac{\gamma S_f}{2}\, r .
$$

The shear stress is zero on the centreline and grows linearly to its wall value $\tau_0 = \gamma S_f r_0/2$. Writing $r_0/2 = D/4 = R$ gives the important result

$$
\tau_0 = \gamma R S_f ,
$$

which holds for **any** steady uniform flow, laminar or turbulent, in pipes and in channels, because it comes only from a force balance.

In laminar flow, Newton's law of viscosity gives $\tau = -\mu\, du/dr$ (the minus sign is there because $u$ decreases as $r$ increases). Combining and integrating with $u = 0$ at the wall:

$$
u(r) = \frac{\gamma S_f}{4\mu}\left(r_0^2 - r^2\right).
$$

The velocity profile is a **parabola**. The maximum velocity is on the centreline, and integrating $u$ over the cross-section shows that the mean velocity is exactly half of it, $V = u_{max}/2$. Solving for the head loss $h_f = S_f L$ over a length $L$ gives the **Hagen–Poiseuille equation**:

$$
h_f = \frac{32\, \mu L V}{\gamma D^2}.
$$

In laminar flow the head loss is proportional to $V$ and to viscosity, and does not depend on the wall roughness. In Chapter 5 we will write every pipe head loss in the **Darcy–Weisbach** form $h_f = f (L/D) V^2/2g$, where $f$ is the **friction factor**. Setting the two expressions equal gives

$$
f = \frac{64}{Re} \qquad (\text{laminar flow}).
$$

## 4.3 Turbulence

If we record the velocity at one point in a turbulent flow, the signal jumps up and down rapidly and randomly. We split the instantaneous velocity into a **time-averaged** part and a **fluctuation**:

$$
u = \bar{u} + u' ,
$$

where $\bar{u}$ is averaged over a time much longer than the eddy time scale, so the average of $u'$ is zero. In hydraulics "velocity" in a turbulent flow almost always means $\bar{u}$. The size of the fluctuations is measured by the **turbulence intensity**, the root-mean-square of $u'$ divided by $\bar{u}$; in pipes and rivers it is of the order of 10%.

Eddies carry fast fluid towards the wall and slow fluid away from it. This exchange transfers momentum and acts like an extra shear stress, the **Reynolds stress** $-\rho\, \overline{u'v'}$, where $v'$ is the velocity fluctuation normal to the wall. Boussinesq proposed writing the total shear stress in the same form as Newton's law:

$$
\tau = \rho \left(\nu + \nu_t\right) \frac{d\bar{u}}{dy},
$$

where $\nu_t$ is the **eddy viscosity**. Unlike $\nu$, the eddy viscosity is not a property of the fluid but of the flow: it is zero at a wall and can be thousands of times larger than $\nu$ in the core of a pipe or river.

Strong mixing makes the turbulent velocity profile much **flatter** than the laminar parabola (Figure 4.1). The velocity rises steeply in a thin layer next to the wall, which means a large wall shear stress. For this reason turbulent flow loses much more energy than laminar flow at the same mean velocity.

```{code-cell} python
:tags: [hide-input]
s = np.linspace(-1, 1, 401)                      # r / r0
u_lam = 2 * (1 - s**2)                           # u / V, laminar
n = 7                                            # power-law exponent 1/n
V_over_umax = 2 * n**2 / ((n + 1) * (2 * n + 1))
u_turb = (1 - abs(s)) ** (1 / n) / V_over_umax   # u / V, turbulent

fig, ax = plt.subplots(figsize=(6.5, 3.6))
ax.plot(u_lam, s, color=c1, lw=2)
ax.plot(u_turb, s, color=c2, lw=2)
ax.axvline(1, color=ink2, lw=1, ls="--")
ax.text(1.04, 0.85, "mean velocity $V$", color=ink2)
ax.text(1.55, 0.35, "laminar\n(parabola)", color=c1)
ax.text(1.4, -0.8, "turbulent\n(1/7 power law)", color=c2)
ax.set_xlabel("Velocity $u/V$", color=ink2)
ax.set_ylabel("Radial position $r/r_0$", color=ink2)
ax.set_xlim(0, 2.2); ax.set_ylim(-1, 1)
tidy(ax); plt.tight_layout(); plt.show()
```

*Figure 4.1. Velocity profiles across a pipe for the same mean velocity $V$. The laminar profile is a parabola with $u_{max} = 2V$. The turbulent profile, shown with the empirical 1/7 power law (a fair fit near $Re \approx 10^5$), is flat in the core with $u_{max} \approx 1.22V$.*

## 4.4 Boundary layers

When a uniform stream of velocity $U$ meets a thin flat plate, the fluid touching the plate stops (no-slip). Viscosity spreads this slowing outwards, so a layer of retarded fluid grows along the plate. This layer is the **boundary layer**. Its **thickness** $\delta$ is usually defined as the distance from the wall where the velocity reaches $0.99U$. Outside the boundary layer the velocity gradients are small, so viscous stresses are negligible and the flow behaves almost as an ideal fluid (Chapter 3).

Near the leading edge the boundary layer is laminar. Its thickness was found by Blasius (1908):

$$
\frac{\delta}{x} \approx \frac{5.0}{\sqrt{Re_x}}, \qquad Re_x = \frac{U x}{\nu},
$$

where $x$ is the distance from the leading edge. When $Re_x$ reaches about $5\times10^5$ the layer becomes turbulent, and it then grows faster. A widely used estimate for the turbulent layer, based on the 1/7 power law, is

$$
\frac{\delta}{x} \approx \frac{0.37}{Re_x^{1/5}} .
$$

Even under a turbulent boundary layer, a very thin layer next to a smooth wall remains dominated by viscosity. This is the **viscous sublayer** (sometimes called the laminar sublayer), described in Section 4.5.

```{code-cell} python
:tags: [hide-input]
U = 1.0                                   # m/s
x = np.linspace(1e-4, 3.0, 600)           # m
Rex = U * x / nu
x_tr = 5e5 * nu / U                       # transition point
d_lam = 5.0 * x / np.sqrt(Rex)
d_turb = 0.37 * x / Rex**0.2

fig, ax = plt.subplots(figsize=(6.5, 3.6))
lam = x <= x_tr
ax.plot(x[lam], d_lam[lam] * 1e3, color=c1, lw=2.5)
ax.plot(x[~lam], d_lam[~lam] * 1e3, color=c1, lw=1.2, ls=":")
ax.plot(x[~lam], d_turb[~lam] * 1e3, color=c2, lw=2.5)
ax.plot(x[lam], d_turb[lam] * 1e3, color=c2, lw=1.2, ls=":")
ax.axvline(x_tr, color=ink2, lw=1, ls="--")
ax.text(x_tr + 0.05, 52, f"transition\n$Re_x = 5×10^5$\n($x$ = {x_tr:.2f} m)", color=ink2, va="top")
ax.text(2.0, 12, "laminar", color=c1)
ax.text(2.0, 49, "turbulent", color=c2)
ax.set_xlabel("Distance from leading edge $x$ (m)", color=ink2)
ax.set_ylabel("Thickness $\\delta$ (mm)", color=ink2)
ax.set_xlim(0, 3); ax.set_ylim(0, 60)
tidy(ax); plt.tight_layout(); plt.show()
```

*Figure 4.2. Boundary-layer thickness along a smooth flat plate in water at 20 °C with $U$ = 1 m s⁻¹. Solid lines show the layer that actually exists; dotted lines extend each formula outside its range. After transition the layer thickens quickly.*

**Pipes.** At the entrance of a pipe a boundary layer grows from the wall in the same way. After some distance the layer fills the whole pipe and the flow becomes fully developed. This **entrance length** is roughly $0.06\,Re\,D$ in laminar flow, but only about 20–40 diameters in turbulent flow. Most pipes are much longer than this, so we normally treat the whole pipe as fully developed.

## 4.5 Velocity distribution near a wall

Close to a wall, the important quantities are the wall shear stress $\tau_0$, the density and the viscosity. They combine into a velocity scale called the **shear velocity** (or friction velocity):

$$
u_* = \sqrt{\frac{\tau_0}{\rho}} .
$$

It is not a velocity that you can measure with a current meter; it is a way of writing the wall stress in velocity units. With $\tau_0 = \gamma R S_f$ from Section 4.2, $u_* = \sqrt{g R S_f}$. In pipes and rivers $u_*$ is typically 3–10% of the mean velocity.

**Viscous sublayer.** Very near a smooth wall the eddies are damped and the stress is purely viscous, $\tau_0 = \mu\, du/dy$. The velocity is linear in $y$, the distance from the wall:

$$
\frac{u}{u_*} = \frac{u_* y}{\nu} .
$$

**Logarithmic layer.** Further from the wall, turbulent stress dominates. Prandtl assumed that the size of the eddies is proportional to the distance from the wall, which gives $\nu_t = \kappa u_* y$, with $\kappa \approx 0.40$ the **von Kármán constant**. With $\tau \approx \tau_0$ near the wall, the eddy-viscosity equation of Section 4.3 integrates to a logarithm. For a smooth wall the constants measured by Nikuradse in pipes give

$$
\frac{u}{u_*} = \frac{1}{\kappa}\ln\frac{u_* y}{\nu} + 5.5 .
$$

The two laws meet at $u_* y/\nu \approx 11.6$, which gives a convenient estimate of the viscous sublayer thickness:

$$
\delta' \approx \frac{11.6\, \nu}{u_*} .
$$

In reality the change is gradual, through a **buffer layer** between about $u_*y/\nu = 5$ and 30 (Figure 4.3).

```{code-cell} python
:tags: [hide-input]
kappa, B = 0.40, 5.5
up = np.linspace(0, 30, 3000)                    # u+ = u / u*
# Spalding (1961) single formula for the whole wall region
yp = up + np.exp(-kappa * B) * (np.exp(kappa * up) - 1 - kappa * up
                                - (kappa * up)**2 / 2 - (kappa * up)**3 / 6)
yy = np.logspace(-0.5, 4, 300)

fig, ax = plt.subplots(figsize=(6.5, 3.8))
ax.axvspan(0.3, 5, color="#eef3fb"); ax.axvspan(5, 30, color="#f6f5f1")
ax.semilogx(yy[yy < 15], yy[yy < 15], color=ink2, lw=1.2, ls="--")
ax.semilogx(yy[yy > 5], np.log(yy[yy > 5]) / kappa + B, color=ink2, lw=1.2, ls="--")
ax.semilogx(yp, up, color=c1, lw=2.5)
ax.text(0.45, 26, "viscous\nsublayer", color=ink2, va="top")
ax.text(6.5, 26, "buffer", color=ink2, va="top")
ax.text(300, 26, "log layer", color=ink2, va="top")
ax.text(0.9, 5.5, "$u^+ = y^+$", color=ink2)
ax.text(1000, 14, "$u^+ = 2.5\\,\\ln y^+ + 5.5$", color=ink2, ha="center")
ax.set_xlabel("$y^+ = u_* y / \\nu$", color=ink2)
ax.set_ylabel("$u^+ = u / u_*$", color=ink2)
ax.set_xlim(0.3, 1e4); ax.set_ylim(0, 28)
tidy(ax); plt.tight_layout(); plt.show()
```

*Figure 4.3. Mean velocity near a smooth wall in wall units. Dashed lines: the linear sublayer law and the log law. Solid line: the single composite formula of Spalding (1961), which joins them smoothly through the buffer layer.*

**Smooth and rough walls.** Real walls are rough. Nikuradse (1933) glued uniform sand grains of size $k_s$ to the inside of pipes, and engineers now describe any real roughness by its **equivalent sand roughness** $k_s$: the sand size that gives the same friction. What matters is how $k_s$ compares with the sublayer thickness, measured by the **roughness Reynolds number** $k_s u_*/\nu$:

- $k_s u_*/\nu < 5$: **hydraulically smooth**. The roughness is hidden inside the viscous sublayer and has no effect.
- $5 < k_s u_*/\nu < 70$: **transitional**. Both viscosity and roughness matter.
- $k_s u_*/\nu > 70$: **hydraulically rough** (fully rough). The roughness elements shed their own eddies, the sublayer is destroyed and viscosity no longer matters.

The same pipe can therefore be smooth at low velocity and rough at high velocity. Over a fully rough wall the log law takes the form

$$
\frac{u}{u_*} = \frac{1}{\kappa}\ln\frac{y}{k_s} + 8.5 .
$$

Viscosity has disappeared from this equation. This is why friction in fully rough pipes (Chapter 5) and in rivers with gravel beds depends only on relative roughness, not on the Reynolds number.

**Depth-averaged velocity in a wide river.** If the rough log law holds over the whole depth $h$ of a wide river, averaging it over the depth gives $V/u_* = 2.5\ln(h/k_s) + 6.0$. The point velocity equals this mean at $y = h/e \approx 0.37h$ above the bed. This explains the familiar stream-gauging rule: measure the velocity at **0.6 of the depth below the surface** to obtain the mean velocity in a vertical.

## 4.6 Flow separation

Outside the boundary layer, the pressure along a body follows the Bernoulli equation. Where the outer flow speeds up, pressure falls in the flow direction (a **favourable pressure gradient**). Where it slows down, pressure rises (an **adverse pressure gradient**). The slow fluid near the wall has little momentum. In an adverse gradient it can be brought to rest and pushed backwards. The boundary layer then leaves the wall: this is **flow separation**. Behind the separation point lies a **wake** of slow, recirculating eddies.

Separation occurs wherever the flow must slow down quickly: behind a bridge pier, at a sudden pipe expansion, on the inside of a sharp bend, downstream of a partly open valve and in a diffuser that widens too fast. The wake converts kinetic energy into turbulence and finally into heat. This is the origin of most **minor losses** in pipe fittings (Chapter 5). Separation zones also produce strong eddies at the bed, which cause local scour around piers and abutments.

A turbulent boundary layer resists separation better than a laminar one, because eddies keep bringing fast fluid down to the wall. This explains the surprising drag behaviour of a sphere in Section 4.7.

## 4.7 Drag and lift on immersed bodies

A body in a moving fluid feels a force. The component in the direction of the approaching flow is the **drag** $F_D$; the component normal to it is the **lift** $F_L$. Both are written with dimensionless coefficients:

$$
F_D = C_D\, A\, \frac{\rho U^2}{2}, \qquad
F_L = C_L\, A\, \frac{\rho U^2}{2},
$$

where $U$ is the approach velocity and $A$ is a reference area, normally the **frontal area** projected on a plane normal to the flow. Drag has two parts:

- **Friction drag** (skin friction) is the sum of wall shear stresses. It dominates for thin, **streamlined** bodies aligned with the flow.
- **Form drag** (pressure drag) comes from the low pressure in the wake behind a separated flow. It dominates for **bluff** bodies such as spheres, cylinders and piers.

Streamlining a body delays separation and shrinks the wake, which is why bridge piers are given rounded or pointed noses.

**Drag on a sphere.** For a sphere of diameter $d$, $C_D$ depends only on $Re = U d/\nu$ (Figure 4.4). At very low Reynolds numbers ($Re \lesssim 1$) inertia is negligible and Stokes (1851) showed that

$$
F_D = 3\pi \mu U d ,
$$

which is the same as $C_D = 24/Re$. Between about $Re = 10^3$ and $2\times10^5$, $C_D$ is nearly constant at 0.4–0.5. Near $Re \approx 3\times10^5$ the boundary layer becomes turbulent before it separates. Separation moves to the rear, the wake narrows, and $C_D$ drops suddenly to about 0.1–0.2. This **drag crisis** is the reason golf balls have dimples. The fall velocity of sediment grains, which balances drag against submerged weight, is treated in Chapter 8.

```{code-cell} python
:tags: [hide-input]
def CD_sphere(Re):
    """Drag coefficient of a smooth sphere, empirical fit of Morrison (2013)."""
    return (24 / Re
            + 2.6 * (Re / 5.0) / (1 + (Re / 5.0) ** 1.52)
            + 0.411 * (Re / 2.63e5) ** -7.94 / (1 + (Re / 2.63e5) ** -8.0)
            + 0.25 * (Re / 1e6) / (1 + Re / 1e6))

Re = np.logspace(-1, 6, 400)
fig, ax = plt.subplots(figsize=(6.5, 3.8))
ax.loglog(Re, CD_sphere(Re), color=c1, lw=2.5)
Rs = np.logspace(-1, 1.3, 50)
ax.loglog(Rs, 24 / Rs, color=c2, lw=1.5, ls="--")
ax.text(0.6, 120, "Stokes: $C_D = 24/Re$", color=c2)
ax.text(1.5e3, 0.75, "$C_D$ ≈ 0.4–0.5", color=ink2)
ax.annotate("drag crisis", xy=(3.5e5, 0.15), xytext=(2e3, 0.12), color=ink2,
            arrowprops=dict(arrowstyle="-|>", color=ink2))
ax.set_xlabel("Reynolds number $Re = Ud/\\nu$", color=ink2)
ax.set_ylabel("Drag coefficient $C_D$", color=ink2)
ax.set_xlim(0.1, 1e6); ax.set_ylim(0.05, 300)
tidy(ax); plt.tight_layout(); plt.show()
```

*Figure 4.4. Drag coefficient of a smooth sphere, from the empirical fit of Morrison (2013), which follows the standard experimental curve from $Re = 0.1$ to $10^6$. The dashed line is Stokes' law, valid only for $Re \lesssim 1$.*

| Body | $C_D$ |
|---|---|
| Sphere | 0.47 |
| Long circular cylinder, axis normal to flow | 1.2 |
| Long flat strip normal to flow | 2.0 |
| Square flat plate normal to flow | 1.1–1.2 |
| Cube, face normal to flow | 1.05 |
| Streamlined body (teardrop shape) | 0.04–0.1 |

*Table 4.2. Approximate drag coefficients for common shapes at $Re \approx 10^4$–$10^5$, based on the frontal area. Values for bodies with sharp edges hardly change with $Re$, because the separation points are fixed at the edges.*

**Lift.** Lift appears whenever the flow around a body is not symmetric, so that the pressure on one side is lower than on the other (the Bernoulli equation, Chapter 3). Hydrofoils and pump impeller blades are designed to produce lift. In rivers, lift on bed grains helps the flow to pick them up (Chapter 8).

## 4.8 Worked examples

**Example 4.1 — Laminar or turbulent?** Water at 20 °C flows at 0.2 L s⁻¹ through a pipe of 25 mm diameter. Is the flow laminar or turbulent? What is the largest discharge for which it would stay laminar?

```{code-cell} python
nu, D, Q = 1.004e-6, 0.025, 0.2e-3
A = np.pi * D**2 / 4
V = Q / A
Re = V * D / nu
Q_lam = 2000 * nu * A / D          # Re = 2000
print(f"V = {V:.3f} m/s, Re = {Re:,.0f}  -> turbulent")
print(f"Largest laminar discharge ≈ {Q_lam*1e3:.3f} L/s")
```

The flow is turbulent with $Re$ ≈ 10,000. It would be laminar only below about 0.04 L s⁻¹, a thin dribble from a tap. Water pipes in buildings are almost always turbulent.

**Example 4.2 — Head loss in a small tube.** Water at 20 °C flows at a mean velocity of 0.5 m s⁻¹ through a horizontal tube with a diameter of 2 mm and a length of 1 m. Find the head loss and the friction factor.

```{code-cell} python
mu, rho, g = 1.002e-3, 998.2, 9.81
D, L, V = 0.002, 1.0, 0.5
Re = V * D / nu
hf = 32 * mu * L * V / (rho * g * D**2)       # Hagen–Poiseuille
f = 64 / Re
print(f"Re = {Re:.0f} (laminar)")
print(f"h_f = {hf:.3f} m, f = {f:.4f}")
print(f"Check with Darcy–Weisbach: {f * L / D * V**2 / (2 * g):.3f} m")
```

$Re$ ≈ 1000, so the flow is laminar and the Hagen–Poiseuille equation applies. The loss of about 0.41 m in only 1 m of tube shows why small tubes and soil pores resist flow so strongly.

**Example 4.3 — Shear velocity and roughness regime.** In a 0.30 m diameter cast-iron water main, the hydraulic grade line falls 2.0 m over 500 m. Taking $k_s$ = 0.26 mm, find the wall shear stress, the shear velocity and the viscous sublayer thickness. Is the pipe hydraulically smooth or rough?

```{code-cell} python
D, L, dh, ks = 0.30, 500.0, 2.0, 0.26e-3
S = dh / L
R = D / 4
tau0 = rho * g * R * S
u_star = np.sqrt(tau0 / rho)
delta_v = 11.6 * nu / u_star
ks_plus = ks * u_star / nu
print(f"tau0 = {tau0:.2f} Pa, u* = {u_star:.4f} m/s")
print(f"sublayer thickness = {delta_v*1e3:.3f} mm")
print(f"ks u*/nu = {ks_plus:.1f}")
```

The roughness Reynolds number is about 14, between 5 and 70, so the pipe is in the **transitional** regime: the roughness (0.26 mm) is similar in size to the sublayer (about 0.2 mm). Both roughness and viscosity affect the friction factor here.

**Example 4.4 — Velocity profile in a gravel-bed river.** A wide river is 2.0 m deep on a slope of 0.0005. The bed is gravel with $k_s$ = 0.05 m. Using the rough-wall log law, find the shear velocity, the velocity 0.2 m above the bed and at the surface, and the depth-averaged velocity.

```{code-cell} python
h, S, ks, kappa = 2.0, 0.0005, 0.05, 0.40
u_star = np.sqrt(g * h * S)              # wide channel: R ≈ h
u = lambda y: u_star * (np.log(y / ks) / kappa + 8.5)
V_mean = u_star * (np.log(h / ks) / kappa + 6.0)
print(f"u* = {u_star:.3f} m/s")
print(f"u(0.2 m) = {u(0.2):.2f} m/s, u(surface) = {u(h):.2f} m/s")
print(f"Depth-averaged V = {V_mean:.2f} m/s; u(0.4 h) = {u(0.4*h):.2f} m/s")
```

The velocity at the surface is only about 50% larger than 0.2 m above the bed, because the logarithm grows slowly. The velocity at $0.4h$ above the bed (0.6 of the depth below the surface) is almost exactly the depth average, as Section 4.5 predicted.

**Example 4.5 — Drag on a towed sphere.** A spherical instrument housing of 100 mm diameter is towed through still water at 20 °C at 1.0 m s⁻¹. Estimate the drag force using the sphere fit of Figure 4.4.

```{code-cell} python
d, U = 0.10, 1.0
Re = U * d / nu
CD = CD_sphere(Re)                     # Morrison (2013) fit
A = np.pi * d**2 / 4
FD = CD * A * rho * U**2 / 2
print(f"Re = {Re:.3g}, C_D = {CD:.2f}, drag = {FD:.2f} N")
```

$Re$ ≈ 10⁵ lies on the flat part of the curve, so $C_D$ ≈ 0.43 and the drag is about 1.7 N. Doubling the speed would roughly quadruple the drag, until the drag crisis is reached.

## 4.9 Common mistakes

- **Using the wrong length in the Reynolds number.** For pipes use the diameter and the limit 2000; for channels use the hydraulic radius and the limit 500. Mixing them gives errors of a factor of four.
- **Using $f = 64/Re$ for turbulent flow.** It holds only for laminar flow ($Re < 2000$). Turbulent friction factors come from the Colebrook–White equation or the Moody diagram (Chapter 5).
- **Treating $u_*$ as the flow velocity.** The shear velocity is a measure of wall stress and is usually only 3–10% of the mean velocity.
- **Judging roughness by its absolute size.** Whether a wall is hydraulically smooth or rough depends on $k_s u_*/\nu$, so it changes with the flow.
- **Using Stokes' law outside its range.** $C_D = 24/Re$ is valid only for $Re \lesssim 1$. At $Re = 100$ it underestimates the drag by a factor of about four.
- **Using the wrong area in the drag equation.** For bluff bodies $C_D$ is based on the frontal area, not the surface area.

## 4.10 Exercises

1. Water at 10 °C flows in a 50 mm pipe. Find the range of mean velocity for which the flow is transitional ($2000 < Re < 4000$). Repeat for water at 30 °C, using Table 1.1.
2. Starting from $u(r) = \gamma S_f (r_0^2 - r^2)/(4\mu)$, show by integration that the mean velocity is $V = u_{max}/2$ and derive the Hagen–Poiseuille equation.
3. A river is 1.5 m deep and 60 m wide, and flows at 0.8 m s⁻¹. Compute its Reynolds number based on the hydraulic radius and state the flow regime. What depth of sheet flow at 0.05 m s⁻¹ would be laminar?
4. Water flows at $U$ = 0.3 m s⁻¹ over a smooth flat plate 2 m long. Where does the boundary layer become turbulent? Find its thickness at the end of the plate.
5. A concrete pipe of 1.2 m diameter has $k_s$ = 1.0 mm and a friction slope of 0.002. Find $\tau_0$, $u_*$ and the roughness Reynolds number. Is the pipe hydraulically smooth, transitional or rough?
6. A cylindrical bridge pier of 1.5 m diameter stands in water 4 m deep flowing at 2 m s⁻¹. Using $C_D$ = 1.2, estimate the drag force on the pier. Explain why a pier with a rounded nose and tail would feel less force.

---

*Figures and tables in this chapter are generated by the Python code shown. Text and figures: © Daeha Kim, CC BY-NC 4.0.*
