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

# 9. Similitude and Hydraulic Models

```{admonition} Core ideas
:class: tip
1. A physical relation between $n$ variables that involve $m$ fundamental dimensions can always be rewritten as a relation between **$n - m$ dimensionless groups** (the Buckingham $\Pi$ theorem).
2. Dimensional analysis **reduces the number of experiments** and shows how to present results so that they apply at any scale. It cannot give the form of the function; that still needs experiment or theory.
3. The most important dimensionless groups in hydraulics are **ratios of forces**: Reynolds (inertia/viscosity), Froude (inertia/gravity), Weber (inertia/surface tension), Cauchy (inertia/elasticity) and Euler (pressure/inertia).
4. A model is **dynamically similar** to its prototype when all relevant force ratios are equal. With water in both, Froude and Reynolds similarity cannot be achieved together, so free-surface models follow the **Froude law**.
5. The forces that are not scaled correctly cause **scale effects**. Viscosity and surface tension become too strong in small models, so every model has a minimum size.
```

## 9.1 Why dimensional analysis

Many problems in hydraulics are too complex to solve from first principles. The flow over a spillway, the loss in a pipe bend and the force of waves on a breakwater all depend on turbulence and geometry in ways that equations alone cannot predict. Engineers solve such problems by **experiment**, either in the laboratory or in a scaled **physical model** of the real structure, the **prototype**.

Experiments raise two questions. Which variables should we measure, and how should we plot the results so that we need as few tests as possible? And how do we transfer a result measured in a 1 m long model to a 100 m long spillway? **Dimensional analysis** answers both. It is based on one idea from Chapter 1 (Section 1.2): every physically meaningful equation is **dimensionally homogeneous**. Each additive term has the same dimensions, so the equation can always be written with dimensionless groups only.

In this chapter we use the $[M, L, T]$ system. Recall that force is $[MLT^{-2}]$, pressure and shear stress are $[ML^{-1}T^{-2}]$, dynamic viscosity is $[ML^{-1}T^{-1}]$ and surface tension is $[MT^{-2}]$.

## 9.2 Rayleigh's method

The simplest method, due to Lord Rayleigh, assumes that the dependent variable is a product of powers of the independent variables. We then choose the exponents so that the dimensions balance.

**Discharge over a weir.** Suppose the discharge per unit width $q$ over a broad weir depends only on the head $H$ above the crest and on gravity $g$. We write

$$
q = C\, g^{a} H^{b},
$$

where $C$ is a dimensionless constant. In dimensions:

$$
[L^2 T^{-1}] = [L T^{-2}]^{a} [L]^{b}.
$$

Matching powers of $T$ gives $-1 = -2a$, so $a = 1/2$. Matching powers of $L$ gives $2 = a + b$, so $b = 3/2$. Therefore

$$
q = C \sqrt{g}\, H^{3/2}.
$$

Dimensional analysis has told us that $q$ grows with $H^{3/2}$, which is the form of every weir equation (Chapter 6). It cannot tell us the value of $C$; that comes from energy theory or from calibration.

**Drag on a sphere.** Now suppose the drag force $F$ on a sphere depends on the fluid density $\rho$, the velocity $V$, the diameter $D$ and the viscosity $\mu$. With $F = C \rho^{a} V^{b} D^{c} \mu^{d}$ there are four unknown exponents but only three dimensions, so we can solve for $a$, $b$, $c$ in terms of $d$. The result is

$$
F = C\, \rho V^2 D^2 \left( \frac{\rho V D}{\mu} \right)^{-d}.
$$

Since $d$ is unknown, the bracket can appear with any power, and the honest conclusion is

$$
\frac{F}{\rho V^2 D^2} = \phi\!\left( \frac{\rho V D}{\mu} \right) = \phi(Re).
$$

The drag coefficient is a function of the Reynolds number only, which is exactly what experiments show (Chapter 4). This example shows the limit of Rayleigh's method: when there are more variables than dimensions, the result contains an unknown function, and a more systematic procedure is useful.

## 9.3 The Buckingham Π theorem

The **Buckingham Π theorem** (Buckingham, 1914) states:

> If a physical relation involves $n$ variables and these variables contain $m$ fundamental dimensions, the relation can be rewritten as a relation between $n - m$ independent dimensionless groups, called $\Pi$ groups.

The groups are found by the **method of repeating variables**:

1. List all $n$ variables that affect the problem, including the dependent one. Write the dimensions of each.
2. Count the fundamental dimensions $m$ (usually 3: $M$, $L$, $T$).
3. Choose $m$ **repeating variables** that together contain all $m$ dimensions but do not form a dimensionless group among themselves. In hydraulics a good choice is one fluid property ($\rho$), one velocity ($V$) and one length ($D$ or $h$). Do not choose the dependent variable.
4. Form each $\Pi$ group as the product of one remaining variable with the repeating variables raised to unknown powers. Set the total powers of $M$, $L$ and $T$ to zero and solve.
5. Check that each group is dimensionless, and rearrange the groups into familiar forms if useful.

**Example: pressure drop in a pipe.** The pressure drop $\Delta p$ along a straight pipe depends on the density $\rho$, the mean velocity $V$, the diameter $D$, the viscosity $\mu$, the pipe length $L$ and the wall roughness height $k$. There are $n = 7$ variables and $m = 3$ dimensions, so we expect $7 - 3 = 4$ groups. We choose $\rho$, $V$ and $D$ as repeating variables.

*First group.* Write $\Pi_1 = \Delta p\, \rho^{a} V^{b} D^{c}$:

$$
\begin{aligned}
[M^0 L^0 T^0] &= [M L^{-1} T^{-2}] \\
&\quad \times [M L^{-3}]^{a}\,[L T^{-1}]^{b}\,[L]^{c}.
\end{aligned}
$$

The three balances are:

$$
\begin{aligned}
M &: \ 1 + a = 0 \ \Rightarrow\ a = -1,\\
T &: \ -2 - b = 0 \ \Rightarrow\ b = -2,\\
L &: \ -1 - 3a + b + c = 0 \ \Rightarrow\ c = 0.
\end{aligned}
$$

So $\Pi_1 = \Delta p / (\rho V^2)$, a form of the Euler number.

*Second group.* Write $\Pi_2 = \mu\, \rho^{a} V^{b} D^{c}$:

$$
\begin{aligned}
M &: \ 1 + a = 0 \ \Rightarrow\ a = -1,\\
T &: \ -1 - b = 0 \ \Rightarrow\ b = -1,\\
L &: \ -1 - 3a + b + c = 0 \ \Rightarrow\ c = -1.
\end{aligned}
$$

So $\Pi_2 = \mu/(\rho V D)$. We may invert any group, and $1/\Pi_2 = \rho V D/\mu = Re$ is the Reynolds number.

*Third and fourth groups.* $L$ and $k$ are lengths, so dividing by $D$ makes them dimensionless at once: $\Pi_3 = L/D$ and $\Pi_4 = k/D$.

The original relation $\Delta p = f(\rho, V, D, \mu, L, k)$ therefore becomes

$$
\frac{\Delta p}{\rho V^2} = \phi\!\left( Re,\ \frac{L}{D},\ \frac{k}{D} \right).
$$

We can go one step further with a physical argument. In a long pipe with fully developed flow, the pressure drop is proportional to the length, so $L/D$ must appear as a simple factor. Writing the remaining function as $f/2$:

$$
\Delta p = f\, \frac{L}{D}\, \frac{\rho V^2}{2},
\quad f = f\!\left( Re,\ \frac{k}{D} \right).
$$

This is the **Darcy–Weisbach equation** (Chapter 5). Dimensional analysis alone has shown that the friction factor depends on only two groups, which is why a single Moody diagram serves every pipe, fluid and velocity. Instead of testing six variables separately, we need to vary only two.

**Collapse of data.** The power of dimensionless groups is easiest to see in a plot. Figure 9.1 shows synthetic measurements on sharp-crested weirs of three different heights $P$. Each measured discharge is converted to a discharge coefficient $C_d$ (Chapter 6). Plotted against the head $H$ alone, the data form three separate lines. Plotted against the dimensionless ratio $H/P$, the same data fall onto one line, so a single curve describes weirs of any height.

```{code-cell} python
:tags: [hide-input]
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 11})

ink, ink2, c1 = "#0b0b0b", "#52514e", "#2a78d6"
c2, c3 = "#d9822b", "#3a9a5b"
grid = "#e4e3df"
rho, g, nu, sigma = 998.2, 9.81, 1.004e-6, 0.0728

def style(ax):
    ax.grid(True, which="major", color=grid); ax.set_axisbelow(True)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)

rng = np.random.default_rng(1)
fig, (axL, axR) = plt.subplots(1, 2, figsize=(6.5, 3.4), sharey=True)
for P, col in zip([0.10, 0.25, 0.60], [c1, c2, c3]):
    H = np.linspace(0.03, 0.15, 7)
    Cd = 0.611 + 0.08 * H / P                     # Rehbock-type relation
    q = (2 / 3) * Cd * np.sqrt(2 * g) * H**1.5 * (1 + 0.01 * rng.standard_normal(H.size))
    Cd_meas = q / ((2 / 3) * np.sqrt(2 * g) * H**1.5)   # from "measured" q and H
    axL.plot(H * 100, Cd_meas, "o", color=col, ms=5, label=f"P = {P*100:.0f} cm")
    axR.plot(H / P, Cd_meas, "o", color=col, ms=5)
x = np.linspace(0, 1.6, 50)
axR.plot(x, 0.611 + 0.08 * x, color=ink2, lw=1, zorder=0)
axL.set_xlabel("Head $H$ (cm)", color=ink2)
axL.set_ylabel("Discharge coefficient $C_d$", color=ink2)
axR.set_xlabel("$H/P$", color=ink2)
axL.legend(frameon=False, fontsize=11, loc="upper left", handletextpad=0.1)
style(axL); style(axR)
plt.tight_layout(); plt.show()
```

*Figure 9.1. Synthetic data for sharp-crested weirs of three heights $P$, generated from a Rehbock-type relation $C_d = 0.611 + 0.08H/P$ with 2% random scatter. Left: plotted against the head, the data form three lines. Right: plotted against $H/P$, they collapse onto one line.*

## 9.4 Important dimensionless numbers

In hydraulics the forces that can act on a fluid element are inertia, pressure, gravity, viscosity, surface tension and elasticity. Each standard dimensionless number is a ratio of two of these forces. For a flow with velocity scale $V$ and length scale $L$, the inertia force on an element is of order $\rho V^2 L^2$, the gravity force $\rho g L^3$, the viscous force $\mu V L$, the surface tension force $\sigma L$ and the elastic force $K L^2$. Dividing inertia by each of the others gives the numbers in the table below. (By tradition the Froude and Mach numbers are written as square roots of the force ratio.)

```{list-table} Dimensionless numbers in hydraulics
:header-rows: 1
:widths: 18 27 25 30

* - Number
  - Definition
  - Force ratio
  - Important in
* - Reynolds, $Re$
  - $\dfrac{\rho V L}{\mu} = \dfrac{V L}{\nu}$
  - inertia / viscous
  - Pipe flow, drag, boundary layers
* - Froude, $Fr$
  - $\dfrac{V}{\sqrt{g L}}$
  - (inertia / gravity)$^{1/2}$
  - Open channels, spillways, waves
* - Weber, $We$
  - $\dfrac{\rho V^2 L}{\sigma}$
  - inertia / surface tension
  - Small models, jets, droplets
* - Cauchy, $Ca$
  - $\dfrac{\rho V^2}{K}$
  - inertia / elastic
  - Water hammer
* - Mach, $Ma$
  - $\dfrac{V}{c} = \dfrac{V}{\sqrt{K/\rho}}$
  - (inertia / elastic)$^{1/2}$
  - Pressure waves; $Ma^2 = Ca$
* - Euler, $Eu$
  - $\dfrac{\Delta p}{\rho V^2}$
  - pressure / inertia
  - Losses, pressure coefficients, cavitation
```

A related group is the **cavitation number** $(p - p_v)/(\rho V^2/2)$, which compares the margin above vapour pressure (Chapter 1) with the dynamic pressure. In the Euler number the pressure is usually the dependent variable: once the other numbers are fixed, $Eu$ follows. This is why a loss coefficient or a drag coefficient measured in a model can be used directly in the prototype.

## 9.5 Similitude

A model gives reliable results only if it is **similar** to the prototype. Three levels of similarity are needed.

- **Geometric similarity.** All lengths in the model are a fixed fraction of the corresponding lengths in the prototype. We define the **length scale ratio** as prototype over model, $L_r = L_p/L_m$. A model at "scale 1:25" has $L_r = 25$. Angles are the same in model and prototype.
- **Kinematic similarity.** The velocities at corresponding points have the same direction and a fixed ratio $V_r = V_p/V_m$. The streamline patterns are then geometrically similar.
- **Dynamic similarity.** The forces at corresponding points have a fixed ratio, so the force polygons are similar. This requires every relevant dimensionless force ratio to have the same value in model and prototype.

Dynamic similarity is the condition that matters, because it guarantees the other two in a geometrically similar model. In practice it is impossible to make all force ratios equal at once, so the modeller identifies the **dominant forces** and scales those correctly.

## 9.6 Model laws

**The Froude law.** In free-surface flow, such as spillways, stilling basins, channels and waves, gravity and inertia dominate. Dynamic similarity then requires the same Froude number in model and prototype, $Fr_m = Fr_p$:

$$
\frac{V_m}{\sqrt{g L_m}} = \frac{V_p}{\sqrt{g L_p}}.
$$

Since $g$ is the same in both, the velocity ratio is

$$
V_r = L_r^{1/2}.
$$

All other ratios follow from the dimensions of each quantity. Time is length divided by velocity, so $t_r = L_r / L_r^{1/2} = L_r^{1/2}$. Discharge is velocity times area, so $Q_r = L_r^{1/2} L_r^{2} = L_r^{5/2}$. With the same fluid ($\rho_r = 1$), force scales as mass times acceleration, $F_r = L_r^{3} \times 1 = L_r^{3}$, and pressure as force per area, $p_r = L_r$. Power is force times velocity, $P_r = L_r^{7/2}$.

```{code-cell} python
:tags: [hide-input]
Lr = np.logspace(0, 2.3, 100)                  # 1 to 200
lines = [("velocity, time\n" + r"$L_r^{1/2}$", 0.5, c1),
         ("area  " + r"$L_r^{2}$", 2.0, ink2),
         ("discharge\n" + r"$L_r^{5/2}$", 2.5, c2),
         ("force  " + r"$L_r^{3}$", 3.0, c3),
         ("power  " + r"$L_r^{7/2}$", 3.5, ink)]
fig, ax = plt.subplots(figsize=(6.5, 4.4))
for lab, e, col in lines:
    ax.loglog(Lr, Lr**e, color=col, lw=2)
    ax.text(230, 200**e, lab, color=col, fontsize=11, ha="left", va="center")
ax.set_xlim(1, 200); ax.set_ylim(1, 1e9)
ax.set_xlabel("Length scale $L_r = L_p/L_m$", color=ink2)
ax.set_ylabel("Prototype / model ratio", color=ink2)
style(ax)
plt.tight_layout(); fig.subplots_adjust(right=0.74); plt.show()
```

*Figure 9.2. Scale ratios under the Froude law (same fluid in model and prototype) for velocity and time, area, discharge, force and power. At scale 1:50 a model velocity of 1 m s⁻¹ means about 7 m s⁻¹ in the prototype, but a model force of 1 N means 125 kN.*

```{code-cell} python
:tags: [remove-input]
from IPython.display import Markdown
qty = [("Length, depth", "$L_r$", 1), ("Area", "$L_r^{2}$", 2), ("Volume", "$L_r^{3}$", 3),
       ("Velocity", "$L_r^{1/2}$", 0.5), ("Time", "$L_r^{1/2}$", 0.5),
       ("Discharge", "$L_r^{5/2}$", 2.5), ("Force", "$L_r^{3}$", 3),
       ("Pressure", "$L_r$", 1), ("Power", "$L_r^{7/2}$", 3.5),
       ("Manning $n$", "$L_r^{1/6}$", 1 / 6)]
def fmt(x):
    return f"{x:,.0f}" if x >= 100 else f"{x:.3g}"
rows = [f"| {q} | {r} | {fmt(25**e)} | {fmt(50**e)} |" for q, r, e in qty]
Markdown("| Quantity | Ratio | 1:25 | 1:50 |\n|---|---|---|---|\n" + "\n".join(rows))
```

*Table 9.1. Froude-law scale ratios (prototype value divided by model value) for two common model scales, with the same fluid in model and prototype. The Manning $n$ ratio follows from the Manning equation (Chapter 6).*

**The Reynolds law.** In flows without a free surface, such as flow through closed pipes, valves and pumps, or around a submerged body far from the surface, gravity does not affect the flow pattern. The dominant ratio is then inertia to viscosity, and similarity requires $Re_m = Re_p$:

$$
\frac{V_m L_m}{\nu_m} = \frac{V_p L_p}{\nu_p}
\ \Rightarrow\
V_r = \frac{\nu_r}{L_r}.
$$

With the same fluid ($\nu_r = 1$), $V_r = 1/L_r$: a model at scale 1:10 must run ten times **faster** than the prototype. Pressures then scale as $\rho V^2$ and are 100 times larger in the model. This is why Reynolds-law models are often tested with a different fluid (for example air for a water prototype) or, more commonly, run at a Reynolds number that is lower than the prototype but still high enough that the result no longer depends on $Re$ (see Section 9.7).

**Why both laws cannot be satisfied.** Suppose both gravity and viscosity matter. The Froude law requires $V_r = L_r^{1/2}$ and the Reynolds law requires $V_r = \nu_r / L_r$. Both hold only if

$$
\nu_r = L_r^{3/2}.
$$

For a 1:10 model, the model fluid would need a kinematic viscosity 32 times smaller than water, about $3 \times 10^{-8}$ m² s⁻¹. No practical liquid is that thin. With water in both model and prototype ($\nu_r = 1$), a free-surface model can only satisfy the Froude law, and the Reynolds number in the model is $L_r^{3/2}$ times smaller than in the prototype. This unavoidable mismatch is the main source of scale effects.

## 9.7 Scale effects

**Scale effects** are the differences between model and prototype behaviour caused by forces that are not scaled correctly. In a Froude model with water, gravity and inertia are correct, but viscosity and surface tension are relatively too strong. The smaller the model, the larger the error.

**Viscosity.** In the prototype the flow is almost always fully turbulent and rough, so the friction factor and loss coefficients no longer depend on $Re$ (Chapter 5). If the model is large enough that its flow is also turbulent and rough, the coefficients are the same and the scale effect is small. If the model Reynolds number is too low, friction is relatively too high, losses are overestimated and jets and hydraulic jumps lose energy differently. Common guidelines ask for a model Reynolds number, based on depth, of at least about $10^4$, and higher (up to about $10^5$) when air entrainment or energy dissipation is studied. Figure 9.3 shows how quickly $Re_m$ falls as the scale becomes smaller.

```{code-cell} python
:tags: [hide-input]
Lr = np.logspace(0, np.log10(200), 200)
fig, ax = plt.subplots(figsize=(6.5, 4.0))
for qp, col in zip([1, 10, 50], [ink2, c1, c2]):
    Re_m = qp / nu / Lr**1.5                    # Re = V h / nu = q / nu
    ax.loglog(Lr, Re_m, color=col, lw=2)
    ax.text(1.15, qp / nu * 1.25, f"$q_p$ = {qp} m² s⁻¹", color=col, fontsize=11, va="bottom")
ax.axhspan(1e2, 1e4, color=c2, alpha=0.10)
ax.text(1.15, 2.2e3, "viscous scale effects likely", color=c2, fontsize=11)
ax.set_xlim(1, 200); ax.set_ylim(1e2, 3e8)
ax.set_xlabel("Length scale $L_r$", color=ink2)
ax.set_ylabel("Model Reynolds number $q_m/\\nu$", color=ink2)
style(ax)
plt.tight_layout(); plt.show()
```

*Figure 9.3. Model Reynolds number $Re_m = q_m/\nu$ (unit discharge over viscosity) for Froude-scaled models of flows with three prototype unit discharges. $Re_m$ falls as $L_r^{-3/2}$; below about $10^4$ viscous effects in the model become important.*

**Surface tension.** Surface tension affects flows with a strongly curved free surface: the nappe of a weir, small waves, vortices at intakes and the break-up of jets and air entrainment. Its importance is measured by the Weber number, which in a Froude model falls as $L_r^{2}$. A very thin nappe in a small model clings to the weir face instead of springing clear, and the discharge coefficient changes. As rules of thumb, model water depths and weir heads should exceed about 3 cm, and wave heights about 2 cm. Air entrainment is the hardest process to model, because it depends on viscosity and surface tension together; a small model usually entrains too little air.

**Other effects.** The model bed cannot always be made smooth enough. Under the Froude law Manning's $n$ scales as $L_r^{1/6}$ (Table 9.1), so a 1:50 model needs $n$ about 0.52 times the prototype value. A concrete spillway with $n = 0.014$ would need a model surface with $n \approx 0.0073$, which is close to the smoothest practical materials. Sediment is another difficulty: fine sand scaled by $L_r$ would become cohesive silt, so **movable-bed models** use light materials (such as plastic or crushed coal) chosen so that the Shields parameter of Chapter 8 is similar in model and prototype.

## 9.8 Distorted models

Rivers, estuaries and harbours are very wide compared with their depth. A river reach 10 km long and 200 m wide, but only 4 m deep, would be 4 cm deep at a scale of 1:100, too shallow for turbulent flow and accurate measurement. A long reach would also not fit in the laboratory at a larger scale. The solution is a **distorted model**, with a horizontal scale $X_r$ larger than the vertical scale $Y_r$. The ratio $X_r/Y_r$ is the **distortion**, typically 2–10.

Gravity waves and the Froude number depend on depth, so the Froude law is applied with the vertical scale:

$$
\begin{aligned}
V_r &= Y_r^{1/2},\\
Q_r &= X_r\, Y_r\, V_r = X_r\, Y_r^{3/2},\\
t_r &= X_r / V_r = X_r\, Y_r^{-1/2}.
\end{aligned}
$$

Slopes are exaggerated by the distortion, $S_r = Y_r/X_r$. For the friction in the model to match, the Manning equation $V = R^{2/3} S^{1/2}/n$ (with $R \approx h$ in a wide channel) requires

$$
n_r = \frac{Y_r^{2/3}}{X_r^{1/2}}.
$$

This ratio is often less than 1, so the model must be **rougher** than the prototype. Modellers add roughness with gravel, wire mesh or small blocks, and then **calibrate** the model: they adjust the roughness until the model reproduces water levels measured in the field. Distorted models reproduce water levels and discharges well, but not local three-dimensional flows such as those around piers or in sharp bends, because the shapes of eddies are distorted.

## 9.9 Hydraulic modelling practice

Building a physical model is a balance between accuracy and cost. The main steps are:

1. **Define the questions.** What must the model show: water levels, discharge capacity, pressures, scour, wave overtopping? This determines which forces dominate and therefore the model law.
2. **Choose the scale.** The model should be as large as the laboratory space, pump capacity and budget allow. Typical scales are 1:20 to 1:60 for spillways and stilling basins, 1:30 to 1:80 for coastal structures, and 1:50 to 1:200 (horizontal) for rivers. Check the model Reynolds and Weber numbers for scale effects.
3. **Build and calibrate.** Construct the model with accurate geometry, then compare it with any available prototype data (water levels, rating curves) and adjust roughness if needed.
4. **Test and interpret.** Run the design cases, convert the results to prototype values with the scale ratios, and judge which results may be affected by scale effects.

**Numerical models.** Today many questions are answered with **numerical models** that solve the flow equations on a computer, from one-dimensional river models to three-dimensional turbulence models. They have no scale effects and can test many options quickly. However, they depend on assumptions about turbulence, roughness and air entrainment, and they must be validated against data. Physical and numerical models are best used together: the physical model checks the critical details, and the numerical model extends the results to more cases. This combined approach is called **hybrid modelling**.

## 9.10 Worked examples

**Example 9.1 — Spillway model.** A spillway is studied in a Froude model at scale 1:40. The prototype design discharge is 1200 m³ s⁻¹. (a) What discharge must the model pump supply? (b) A velocity of 1.8 m s⁻¹ is measured at the toe of the model. What is the prototype velocity? (c) A gate in the model opens in 30 s. How long does it take in the prototype? (d) The force on a model training wall is 12 N. What is the prototype force?

```{code-cell} python
Lr = 40
Qp = 1200.0
Qm = Qp / Lr**2.5
Vp = 1.8 * Lr**0.5
tp = 30 * Lr**0.5
Fp = 12 * Lr**3
print(f"(a) Q_m = {Qm*1e3:.0f} L/s")
print(f"(b) V_p = {Vp:.1f} m/s")
print(f"(c) t_p = {tp:.0f} s = {tp/60:.1f} min")
print(f"(d) F_p = {Fp/1e3:.0f} kN")
```

The model needs about 119 L s⁻¹, a typical laboratory pump capacity. Forces grow very fast with scale: 12 N in the model become almost 770 kN in the prototype.

**Example 9.2 — Checking for scale effects.** In the spillway of Example 9.1 the crest is 60 m wide and the design head is 4.0 m. Check the model Reynolds number (based on unit discharge) and the Weber number (based on head and the mean crest velocity).

```{code-cell} python
B, Hp = 60.0, 4.0
qm = Qp / B / Lr**1.5           # m^2/s
Hm = Hp / Lr
Vm = qm / Hm
Re_m = qm / nu
We_m = rho * Vm**2 * Hm / sigma
print(f"Model head = {Hm*100:.0f} cm, q_m = {qm:.4f} m^2/s, V_m = {Vm:.2f} m/s")
print(f"Re_m = {Re_m:,.0f}, We_m = {We_m:.0f}")
```

The model head is 10 cm, well above the 3 cm rule, and the Reynolds number is about 79 000, above the usual limit of $10^4$. Discharge and water levels should be reliable; air entrainment in the jump would still be underestimated.

**Example 9.3 — Choosing the scale.** A stilling basin 30 m wide must pass 900 m³ s⁻¹. The laboratory flume is 1.0 m wide and the pump delivers at most 120 L s⁻¹. What is the largest possible model (smallest $L_r$)?

```{code-cell} python
Lr_width = 30.0 / 1.0
Lr_pump = (900.0 / 0.120) ** (1 / 2.5)
Lr_min = max(Lr_width, Lr_pump)
print(f"Width limit: L_r >= {Lr_width:.0f};  pump limit: L_r >= {Lr_pump:.1f}")
print(f"Choose L_r = 40: model Q = {900/40**2.5*1e3:.0f} L/s, model width = {30/40:.2f} m")
```

The pump, not the flume width, controls the scale here: $L_r$ must be at least 35.5. A round scale of 1:40 needs 89 L s⁻¹, which leaves a margin of about 25% on the pump.

**Example 9.4 — Distorted river model.** A river reach 10 km long, with mean width 200 m, depth 4.0 m, Manning $n = 0.030$ and flood discharge 800 m³ s⁻¹, is modelled with horizontal scale 1:250 and vertical scale 1:50. Find the model length, depth, discharge and required Manning $n$, the time scale, and the model Reynolds number $V_m h_m/\nu$.

```{code-cell} python
X, Y = 250.0, 50.0
Lm, hm = 10e3 / X, 4.0 / Y
Qm = 800.0 / (X * Y**1.5)
n_r = Y**(2 / 3) / X**0.5
n_m = 0.030 / n_r
t_r = X / Y**0.5
Vp = 800.0 / (200.0 * 4.0)
Vm = Vp / Y**0.5
print(f"Distortion = {X/Y:.0f}; model length = {Lm:.0f} m, depth = {hm*100:.0f} cm")
print(f"Model discharge = {Qm*1e3:.1f} L/s")
print(f"n_r = {n_r:.3f}, model n = {n_m:.3f}")
print(f"Time scale = {t_r:.1f}: 1 h in prototype = {3600/t_r/60:.1f} min in model")
print(f"Re_m = {Vm*hm/nu:,.0f}")
```

The model is 40 m long and 8 cm deep and needs only about 9 L s⁻¹. Its Manning $n$ must be about 0.035, rougher than the prototype, so artificial roughness is needed. A one-hour flood in the river lasts under two minutes in the model, and the Reynolds number of about 11 000 is just acceptable.

## 9.11 Common mistakes

- **Mixing up the direction of the scale ratio.** In this book $L_r = L_p/L_m$, so prototype values are model values times the ratio. Check that the prototype comes out larger.
- **Using the Froude law for everything.** Closed-conduit flow without a free surface is controlled by viscosity, not gravity. Use the Reynolds law, or check that $Re$ is high enough that it no longer matters.
- **Choosing the dependent variable as a repeating variable.** It must appear in only one $\Pi$ group, otherwise the result cannot be solved for it.
- **Leaving out a relevant variable.** Dimensional analysis cannot tell you what you forgot. If surface tension matters but $\sigma$ is not on the list, the Weber number never appears.
- **Ignoring scale effects in small models.** A model with 1 cm depths will show surface tension and laminar flow effects that do not exist in the prototype.
- **Treating dimensional analysis as a solution.** It gives the form of the relation, never the function or the constants. Those come from experiment.

## 9.12 Exercises

1. The period $T$ of small waves on deep water depends on the wavelength $\lambda$, gravity $g$ and density $\rho$. Use Rayleigh's method to find the form of the relation. Why does $\rho$ drop out?
2. The thrust $F$ of a ship's propeller depends on its diameter $D$, speed of rotation $N$, forward speed $V$, water density $\rho$ and viscosity $\mu$. Use the Buckingham $\Pi$ theorem with $\rho$, $N$ and $D$ as repeating variables to find the dimensionless groups.
3. The discharge $Q$ through a small orifice depends on the head $H$, the orifice diameter $d$, gravity $g$, density $\rho$, viscosity $\mu$ and surface tension $\sigma$. Find a set of $\Pi$ groups and name the familiar numbers among them.
4. A 1:25 model of a harbour is tested with waves of period 1.6 s and height 0.12 m. Find the prototype wave period and height. If the measured force on a model caisson is 85 N per metre length, what is the prototype force per metre length?
5. Show that, for a Froude model with the same fluid, the Weber number scales as $L_r^{2}$ and the Reynolds number as $L_r^{3/2}$. For a prototype weir head of 1.5 m, what is the largest scale ratio $L_r$ (the smallest model) that keeps the model head above 3 cm?
6. Design a distorted model of an estuary 30 km long and 2 km wide with mean depth 6 m, for a laboratory hall 60 m long. Choose the horizontal and vertical scales, and find the time scale for a tidal period of 12.4 h.

---

*Figures and tables in this chapter are generated by the Python code shown. Text and figures: © Daeha Kim, CC BY-NC 4.0.*
