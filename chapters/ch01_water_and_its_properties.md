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

# 1. Water and Its Properties

```{admonition} Core ideas
:class: tip
1. A fluid is a substance that **deforms continuously under any shear stress**, however small. A solid resists shear by deforming a finite amount; a fluid keeps flowing.
2. Hydraulics treats water as a **continuum**: properties such as density and velocity are smooth functions of position, even though water is made of molecules.
3. Every physically meaningful equation is **dimensionally homogeneous**. Checking dimensions is the fastest way to catch an error.
4. Five properties control almost everything in this book: **density** (inertia and weight), **viscosity** (friction), **surface tension** (capillarity), **compressibility** (pressure waves) and **vapour pressure** (cavitation).
5. All five depend on **temperature**. For most engineering problems water can be treated as incompressible, but never as frictionless when it flows over long distances.
```

## 1.1 What hydraulics is about

The word *hydraulics* comes from the Greek *hydor* (water) and *aulos* (pipe). Today it means the application of fluid mechanics to water in engineering: water at rest in tanks and behind gates, water flowing in pipes and channels, waves on the surface, and the sediment that moving water carries. Hydraulics sits between pure fluid mechanics, which seeks general laws, and hydrology, which describes the movement of water through the landscape. Much of hydraulics is therefore a blend of **theory** (conservation of mass, energy and momentum) and **experiment** (friction factors, loss coefficients, sediment transport rates).

Throughout this book we treat water as a **continuum**. A cubic millimetre of water contains about $3\times10^{19}$ molecules, so any volume we care about is large enough that molecular fluctuations average out. Density, pressure and velocity can then be defined at a point and differentiated, which is what allows us to write conservation laws as equations.

## 1.2 Dimensions and units

Every physical quantity can be expressed in terms of a few **fundamental dimensions**. In hydraulics we use either the mass system $[M, L, T]$ or the force system $[F, L, T]$, with temperature $[\Theta]$ when needed. Newton's second law, $F = ma$, links the two:

$$
[F] = [M L T^{-2}] \qquad \Longleftrightarrow \qquad [M] = [F L^{-1} T^{2}].
$$

An equation is **dimensionally homogeneous** when every additive term has the same dimensions. For example, the discharge equation $Q = AV$ is homogeneous because

$$
[L^3 T^{-1}] = [L^2][L T^{-1}].
$$

Dimensional homogeneity is a necessary condition for a physical law. It is not sufficient, because a dimensionless constant can still be wrong, but an inhomogeneous equation is certainly wrong or, at best, an empirical formula that only works in one set of units. Several well-known hydraulic formulas (the Manning equation in Chapter 6 is the classic example) are of this empirical kind, and their coefficients change when the units change.

This book uses **SI units** throughout: kilogram (kg), metre (m), second (s), newton (N = kg m s⁻²) and pascal (Pa = N m⁻²). Older Korean and Japanese texts often use the kilogram-force (kgf), the weight of 1 kg under standard gravity:

$$
1\ \text{kgf} = 1\ \text{kg} \times 9.81\ \text{m s}^{-2} = 9.81\ \text{N}.
$$

Mixing kgf and N is one of the most common sources of errors by a factor of ten.

```{list-table} Common quantities in hydraulics
:header-rows: 1
:widths: 30 20 25 25

* - Quantity
  - Symbol
  - Dimensions
  - SI unit
* - Velocity
  - $V$
  - $L T^{-1}$
  - m s⁻¹
* - Discharge
  - $Q$
  - $L^3 T^{-1}$
  - m³ s⁻¹
* - Pressure, shear stress
  - $p,\ \tau$
  - $F L^{-2} = M L^{-1} T^{-2}$
  - Pa
* - Density
  - $\rho$
  - $M L^{-3}$
  - kg m⁻³
* - Specific weight
  - $\gamma$
  - $F L^{-3} = M L^{-2} T^{-2}$
  - N m⁻³
* - Dynamic viscosity
  - $\mu$
  - $F T L^{-2} = M L^{-1} T^{-1}$
  - Pa s
* - Kinematic viscosity
  - $\nu$
  - $L^2 T^{-1}$
  - m² s⁻¹
* - Surface tension
  - $\sigma$
  - $F L^{-1} = M T^{-2}$
  - N m⁻¹
* - Bulk modulus
  - $K$
  - $F L^{-2}$
  - Pa
```

## 1.3 Density, specific weight and specific gravity

**Density** is mass per unit volume, $\rho = m / \mathcal{V}$. Pure water reaches its maximum density of about 1000 kg m⁻³ at 4 °C; at 20 °C it is 998.2 kg m⁻³. **Specific weight** is weight per unit volume,

$$
\gamma = \rho g ,
$$

about 9.79 kN m⁻³ at 20 °C. **Specific gravity** is the ratio of a substance's density to that of water at 4 °C, so it is dimensionless (mercury: 13.6; seawater: about 1.025).

Density matters in two ways. It sets the **weight** of water, which drives hydrostatic pressure (Chapter 2), and the **inertia** of water, which appears in every momentum calculation (Chapter 3).

## 1.4 Viscosity

Consider water between two parallel plates a distance $h$ apart. The lower plate is fixed and the upper plate moves at speed $U$. Water sticks to both plates (the **no-slip condition**), so the velocity varies from 0 to $U$ across the gap, and the moving plate must be pushed with a force to keep it going. Experiments show that for water the shear stress is proportional to the rate of shear deformation:

$$
\tau = \mu \frac{du}{dy}.
$$

This is **Newton's law of viscosity**, and fluids that obey it, including water and air, are called **Newtonian**. The constant $\mu$ is the **dynamic viscosity**. Dividing by density gives the **kinematic viscosity**,

$$
\nu = \frac{\mu}{\rho},
$$

which is the form that appears in the Reynolds number (Chapter 4). For water at 20 °C, $\mu \approx 1.0\times10^{-3}$ Pa s and $\nu \approx 1.0\times10^{-6}$ m² s⁻¹.

```{code-cell} python
:tags: [hide-input]
import numpy as np
import matplotlib.pyplot as plt

ink, ink2, c1 = "#0b0b0b", "#52514e", "#2a78d6"
fig, ax = plt.subplots(figsize=(6.5, 3.2))
h, U = 1.0, 1.0
y = np.linspace(0, h, 9)
ax.fill_between([0, 1.6], -0.08, 0, color="#d9d8d3")
ax.fill_between([0, 1.6], h, h + 0.08, color="#d9d8d3")
ax.plot([0, 0], [0, h], color=ink2, lw=1)
for yi in y[1:]:
    ax.annotate("", xy=(U * yi / h, yi), xytext=(0, yi),
                arrowprops=dict(arrowstyle="-|>", color=c1, lw=1.4))
ax.plot(U * y / h, y, color=c1, lw=2)
ax.annotate("", xy=(1.55, h + 0.04), xytext=(1.2, h + 0.04),
            arrowprops=dict(arrowstyle="-|>", color=ink, lw=1.6))
ax.text(1.22, h + 0.13, "moving plate, speed $U$", color=ink, fontsize=10)
ax.text(0.05, -0.2, "fixed plate ($u = 0$)", color=ink2, fontsize=10)
ax.text(1.05, 0.45, r"$\tau = \mu \, \dfrac{du}{dy} = \mu \, \dfrac{U}{h}$", fontsize=13, color=ink)
ax.annotate("", xy=(-0.12, h), xytext=(-0.12, 0), arrowprops=dict(arrowstyle="<->", color=ink2))
ax.text(-0.2, h / 2, "$h$", color=ink2, ha="right", va="center", fontsize=12)
ax.set_xlim(-0.35, 1.9); ax.set_ylim(-0.3, 1.3); ax.axis("off")
plt.tight_layout(); plt.show()
```

*Figure 1.1. Flow between a fixed and a moving plate. The no-slip condition forces the velocity to vary linearly across the gap, and the shear stress is the same everywhere in the fluid.*

Viscosity is the property that turns flow into **energy loss**. It is the reason a pipeline needs a pump, and the reason a river needs a slope to keep flowing. Viscosity decreases quickly with temperature: water at 5 °C is about 50% more viscous than at 20 °C.

## 1.5 Surface tension and capillarity

Molecules at a free surface are pulled inward by their neighbours, so the surface behaves like a stretched membrane. **Surface tension** $\sigma$ is the force per unit length acting along the surface (about 0.073 N m⁻¹ for water at 20 °C). It is negligible in most hydraulic structures but dominant at small scales: in thin tubes, in soil pores and in small physical models (Chapter 9).

**Capillary rise.** Water in a narrow vertical tube of diameter $d$ rises until the weight of the raised column balances the vertical component of surface tension around the wetted perimeter. With contact angle $\theta$ between the water surface and the wall:

$$
\underbrace{\rho g \, \frac{\pi d^2}{4} \, h}_{\text{weight of column}}
= \underbrace{\sigma \, \pi d \cos\theta}_{\text{surface tension}}
\qquad \Longrightarrow \qquad
h = \frac{4 \sigma \cos\theta}{\rho g d}.
$$

The rise is inversely proportional to diameter. In a 1 mm glass tube ($\theta \approx 0$) water rises about 3 cm; in a soil pore of 10 µm it would rise about 3 m. This is the physics behind the **capillary fringe** above the water table.

```{code-cell} python
:tags: [hide-input]
rho, g, sigma = 998.2, 9.81, 0.0728
d = np.logspace(-6, -2, 200)            # m
h_cap = 4 * sigma / (rho * g * d)       # m, contact angle 0

fig, ax = plt.subplots(figsize=(6.5, 3.6))
ax.loglog(d * 1e3, h_cap, color=c1, lw=2)
for dd, lab in [(1e-3, "1 mm glass tube"), (1e-5, "10 µm soil pore")]:
    hh = 4 * sigma / (rho * g * dd)
    ax.scatter([dd * 1e3], [hh], color=c1, zorder=3, s=40)
    ax.annotate(f"{lab}: {hh:.2g} m", (dd * 1e3, hh), xytext=(8, 6),
                textcoords="offset points", fontsize=9, color=ink)
ax.set_xlabel("Tube or pore diameter $d$ (mm)", color=ink2)
ax.set_ylabel("Capillary rise $h$ (m)", color=ink2)
ax.grid(True, which="major", color="#e4e3df"); ax.set_axisbelow(True)
for s in ("top", "right"): ax.spines[s].set_visible(False)
plt.tight_layout(); plt.show()
```

*Figure 1.2. Capillary rise of water at 20 °C for a perfectly wetting surface ($\theta = 0$). Halving the diameter doubles the rise.*

## 1.6 Compressibility

When pressure on a volume of water increases by $\Delta p$, the volume decreases slightly. The **bulk modulus of elasticity** measures the resistance to this compression:

$$
K = -\frac{\Delta p}{\Delta \mathcal{V} / \mathcal{V}} = \frac{\Delta p}{\Delta \rho / \rho}.
$$

For water, $K \approx 2.2\times10^{9}$ Pa. Reducing the volume of water by just 1% requires about 22 MPa, the pressure at a depth of 2.2 km in the ocean. For this reason we treat water as **incompressible** in almost every problem in this book. The exception is sudden valve closure in a pipeline (water hammer), where compressibility controls the speed of the pressure wave, $c = \sqrt{K/\rho} \approx 1480$ m s⁻¹.

## 1.7 Vapour pressure and cavitation

Molecules continually escape from a liquid surface. The pressure exerted by the vapour in equilibrium with the liquid is the **vapour pressure** $p_v$, which rises steeply with temperature: about 2.3 kPa at 20 °C and 101.3 kPa (1 atm) at 100 °C, which is why water boils at 100 °C at sea level.

If the **absolute** pressure inside flowing water drops to $p_v$, the water boils locally even at room temperature. Vapour bubbles form in the low-pressure zone and collapse violently when they reach higher pressure downstream. This is **cavitation**, which erodes pump impellers, spillways and valves. It also sets a hard limit on how high a siphon can lift water (Chapter 5): the pressure at the crest cannot fall below $p_v$.

```{code-cell} python
:tags: [hide-input]
T = np.linspace(0, 40, 200)                       # °C
TK = T + 273.15
# Density: Tanaka et al. (2001), 0–40 °C
a1, a2, a3, a4, a5 = -3.983035, 301.797, 522528.9, 69.34881, 999.974950
rho_T = a5 * (1 - (T + a1) ** 2 * (T + a2) / (a3 * (T + a4)))
# Dynamic viscosity: Vogel equation
mu_T = 2.414e-5 * 10 ** (247.8 / (TK - 140.0))
nu_T = mu_T / rho_T
# Surface tension: IAPWS (1994)
tau = 1 - TK / 647.096
sigma_T = 0.2358 * tau ** 1.256 * (1 - 0.625 * tau)
# Saturation vapour pressure: Tetens/Magnus form
pv_T = 0.6108 * np.exp(17.27 * T / (T + 237.3))   # kPa

panels = [(rho_T, "Density ρ (kg m⁻³)"), (nu_T * 1e6, "Kinematic viscosity ν (10⁻⁶ m² s⁻¹)"),
          (sigma_T * 1e3, "Surface tension σ (10⁻³ N m⁻¹)"), (pv_T, "Vapour pressure $p_v$ (kPa)")]
fig, axes = plt.subplots(2, 2, figsize=(8, 5.6), sharex=True)
for ax, (y, lab) in zip(axes.flat, panels):
    ax.plot(T, y, color=c1, lw=2)
    ax.set_title(lab, loc="left", fontsize=10, color=ink)
    ax.grid(True, color="#e4e3df"); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
for ax in axes[1]: ax.set_xlabel("Temperature (°C)", color=ink2)
plt.tight_layout(); plt.show()
```

*Figure 1.3. Properties of pure water between 0 and 40 °C, computed from standard empirical formulas (density: Tanaka et al., 2001; viscosity: Vogel equation; surface tension: IAPWS; vapour pressure: Tetens formula). Note the density maximum near 4 °C and the steep drop in viscosity with temperature.*

```{code-cell} python
:tags: [remove-input]
from IPython.display import Markdown
rows = []
for Tc in [0, 5, 10, 15, 20, 25, 30, 40]:
    i = np.argmin(abs(T - Tc))
    rows.append(f"| {Tc} | {rho_T[i]:.1f} | {rho_T[i]*9.81/1e3:.3f} | {mu_T[i]*1e3:.3f} | {nu_T[i]*1e6:.3f} | {sigma_T[i]*1e3:.1f} | {pv_T[i]:.2f} |")
Markdown("| T (°C) | ρ (kg m⁻³) | γ (kN m⁻³) | μ (10⁻³ Pa s) | ν (10⁻⁶ m² s⁻¹) | σ (10⁻³ N m⁻¹) | $p_v$ (kPa) |\n"
         "|---|---|---|---|---|---|---|\n" + "\n".join(rows))
```

*Table 1.1. Properties of pure water, computed from the same formulas as Figure 1.3. Values agree with standard tables to within about 2% (the largest difference is viscosity near 0 °C).*

## 1.8 Worked examples

**Example 1.1 — Shear stress in a thin water film.** A plate of area 0.5 m² slides at 1.0 m s⁻¹ over a 1.0 mm film of water at 20 °C. Find the force needed.

```{code-cell} python
mu, A, U, h = 1.002e-3, 0.5, 1.0, 1.0e-3
tau = mu * U / h          # Pa
F = tau * A               # N
print(f"Shear stress = {tau:.2f} Pa, force = {F:.2f} N")
```

The force is only about half a newton. Viscous stresses in water are small, but over the large surface area of a long pipe they add up to large energy losses.

**Example 1.2 — Capillary rise in a soil.** Estimate the height of the capillary fringe above a water table in a silty soil whose effective pore diameter is 20 µm. Assume water at 15 °C and $\theta = 0$.

```{code-cell} python
sigma15, rho15, g, d = 0.0735, 999.1, 9.81, 20e-6
h = 4 * sigma15 / (rho15 * g * d)
print(f"Capillary rise ≈ {h:.2f} m")
```

About 1.5 m. Real soils have a range of pore sizes, so the fringe is a gradual transition rather than a sharp line, but the order of magnitude explains why fine soils stay wet far above the water table.

**Example 1.3 — How compressible is water?** What pressure increase reduces the volume of water by 0.5%?

```{code-cell} python
K = 2.2e9
dp = K * 0.005
print(f"Δp = {dp/1e6:.0f} MPa  (≈ {dp/(998*9.81):.0f} m of water)")
```

**Example 1.4 — The highest possible siphon.** Atmospheric pressure is 101.3 kPa and the water is at 20 °C. Ignoring friction and velocity head, how high above the upstream water surface can the crest of a siphon be before the water boils?

```{code-cell} python
p_atm, p_v, gamma = 101.3e3, 2.34e3, 998.2 * 9.81
z_max = (p_atm - p_v) / gamma
print(f"Maximum crest height ≈ {z_max:.1f} m")
```

About 10 m in theory. In practice friction, velocity head and dissolved air reduce the safe limit to roughly 7–8 m, as we will see in Chapter 5.

## 1.9 Common mistakes

- **Mixing kgf and N.** 1 kgf = 9.81 N, not 1 N. Convert old data to SI before calculating.
- **Confusing $\mu$ and $\nu$.** Dynamic viscosity (Pa s) appears in shear stress; kinematic viscosity (m² s⁻¹) appears in the Reynolds number. They differ by a factor of $\rho \approx 1000$.
- **Using radius instead of diameter** in the capillary-rise formula. With radius $r$ the formula is $h = 2\sigma\cos\theta/(\rho g r)$.
- **Comparing gauge pressure with vapour pressure.** Cavitation depends on **absolute** pressure.
- **Assuming water properties are constant.** Viscosity at 5 °C is 50% higher than at 20 °C, which changes friction losses noticeably.

## 1.10 Exercises

1. Show that the Bernoulli sum $p/\gamma + V^2/2g + z$ has dimensions of length.
2. Express dynamic viscosity in the $[F L T]$ system and in the $[M L T]$ system, and show that they are consistent.
3. A 2 mm gap between a rotating shaft and its bearing is filled with water at 10 °C. The shaft surface moves at 3 m s⁻¹ and the wetted area is 0.02 m². Find the friction force and the power lost.
4. A glass tube is used to read the water level in a tank. How small can its diameter be if the capillary error must stay below 2 mm?
5. A pipeline carries water at 30 °C. At what absolute pressure would cavitation begin? Convert your answer to metres of water.

---

*Figures and tables in this chapter are generated by the Python code shown. Text and figures: © Daeha Kim, CC BY-NC 4.0.*
