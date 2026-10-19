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

# 3. Conservation Laws of Flow

```{admonition} Core ideas
:class: tip
1. Hydraulics describes flow by the velocity at fixed points (the **Eulerian** view) and usually reduces it to a **one-dimensional** flow with a mean velocity $V = Q/A$.
2. **Continuity** (conservation of mass): for incompressible water the discharge $Q = AV$ is the same at every section of a pipe or stream tube.
3. **Bernoulli's equation**: along a streamline in steady, frictionless, incompressible flow, the **total head** $p/\gamma + V^2/2g + z$ is constant. Real flows lose head through friction, which we add as a **head loss**.
4. The **energy grade line** and **hydraulic grade line** turn the energy equation into a picture: the gap between them is the velocity head, and the gap between the HGL and the pipe is the pressure head.
5. The **momentum equation** gives the forces that flowing water exerts on nozzles, bends, gates and vanes: the net force on a control volume equals the **momentum flux out minus the momentum flux in**.
```

## 3.1 Describing a flow

There are two ways to describe a moving fluid. In the **Lagrangian** view we follow individual particles and record their positions in time, as we would follow a float drifting down a river. In the **Eulerian** view we stand at fixed points and record the velocity and pressure of whatever fluid passes, as a current meter fixed in the river does. Hydraulics almost always uses the Eulerian view, because measurements and designs refer to fixed locations: a pipe section, a gauging station, the crest of a weir.

**Streamlines.** A **streamline** is a line that is everywhere tangent to the velocity at an instant. No flow crosses a streamline. A bundle of streamlines forms a **stream tube**, which behaves like a pipe with invisible walls. A **pathline** is the track of a single particle. In steady flow the streamline pattern does not change, and pathlines and streamlines coincide.

**Steady and unsteady, uniform and non-uniform.** A flow is **steady** if the velocity at every point does not change with time, and **unsteady** otherwise. A flow is **uniform** if the velocity does not change along the direction of flow, and **non-uniform** otherwise. These are independent properties:

- steady uniform flow: constant discharge in a long straight pipe of constant diameter;
- steady non-uniform flow: constant discharge through a tapering pipe or over a weir;
- unsteady flow: a pipe while a valve is closing, or a river during a flood.

Turbulent flow is never strictly steady, because the velocity at a point fluctuates rapidly. We call it steady when the **time-averaged** velocity is constant. The difference between laminar and turbulent flow is the subject of Chapter 4.

**One-dimensional flow.** The velocity in a pipe is zero at the wall and highest near the centre. For most engineering purposes we replace this profile by a single **mean velocity**,

$$
V = \frac{Q}{A} = \frac{1}{A}\int_A u \, dA ,
$$

where $u$ is the local velocity and $A$ the cross-sectional area. Treating the flow as varying only along its length, with one value of $V$, $p$ and $z$ at each section, is the **one-dimensional approximation**. It is the basis of pipe and channel hydraulics in this book. Section 3.7 shows how to correct it for a non-uniform velocity profile.

## 3.2 Continuity

Mass is neither created nor destroyed. For steady flow through a stream tube or pipe between section 1 and section 2, the mass flow rate in equals the mass flow rate out:

$$
\rho_1 A_1 V_1 = \rho_2 A_2 V_2 .
$$

For water, which we treat as incompressible, the densities cancel and we obtain the **continuity equation**:

$$
Q = A_1 V_1 = A_2 V_2 ,
$$

where the **discharge** $Q$ (m³ s⁻¹) is the volume passing a section per unit time. Where the pipe narrows, the water must speed up; where it widens, the water slows down. Halving the diameter of a pipe multiplies the velocity by four.

At a **junction** of pipes or channels, water cannot accumulate (in steady flow), so

$$
\sum Q_{\text{in}} = \sum Q_{\text{out}} .
$$

This simple statement is the first equation in the analysis of pipe networks (Chapter 5). In unsteady flow, the difference between inflow and outflow equals the rate of change of stored volume, $Q_{\text{in}} - Q_{\text{out}} = d\mathcal{V}/dt$, which we use for draining tanks in Section 3.5.

## 3.3 The Euler and Bernoulli equations

Consider a small cylindrical element of fluid along a streamline, with length $ds$ and cross-section $dA$. Let $z$ be elevation, measured upward. We assume the flow is **steady** and **frictionless**. The forces along the streamline are the net pressure force and the component of weight:

$$
-\frac{\partial p}{\partial s}\, ds\, dA - \gamma\, ds\, dA\, \frac{dz}{ds}.
$$

In steady flow a particle moving along the streamline accelerates because it moves to a place of different velocity, so its acceleration is $a_s = V\, dV/ds$. Newton's second law, divided by the volume $ds\, dA$, gives

$$
-\frac{dp}{ds} - \gamma \frac{dz}{ds} = \rho V \frac{dV}{ds}.
$$

This is **Euler's equation** along a streamline. It says that a particle speeds up when it moves towards lower pressure or downhill. For incompressible flow, dividing by $\gamma$ and integrating along the streamline gives

$$
\frac{p}{\gamma} + \frac{V^2}{2g} + z = H = \text{constant}.
$$

This is the **Bernoulli equation**, named after Daniel Bernoulli, whose *Hydrodynamica* (1738) first stated the principle; its modern form is due to Euler. Each term is an energy per unit weight of water and has dimensions of length, which is why it is called a **head**:

- $p/\gamma$ is the **pressure head** (the height to which water would rise in a piezometer);
- $V^2/2g$ is the **velocity head** (the height through which a body must fall freely to reach speed $V$);
- $z$ is the **elevation head** or datum head, measured from any convenient horizontal datum;
- $H$ is the **total head**.

The sum $p/\gamma + z$ is the piezometric head of Chapter 2. Bernoulli's equation shows that, without friction, a gain in velocity head is paid for by a fall in piezometric head, and vice versa.

**Assumptions.** Bernoulli's equation in this form holds only when the flow is (1) **steady**, (2) **incompressible**, (3) **frictionless**, (4) **along a single streamline**, and (5) has **no pump or turbine** between the two points. In practice it is applied between sections of a pipe or channel, using the mean velocity. This is accurate where the streamlines are straight and parallel, because then the pressure across the section is hydrostatic and $p/\gamma + z$ is the same at every point of the section. It is not accurate across sharply curved streamlines, such as at the crest of a spillway.

**The energy equation with losses.** Real water has viscosity, and friction converts some mechanical energy into heat. Between an upstream section 1 and a downstream section 2,

$$
\frac{p_1}{\gamma} + \frac{V_1^2}{2g} + z_1 = \frac{p_2}{\gamma} + \frac{V_2^2}{2g} + z_2 + h_L ,
$$

where $h_L$ is the **head loss**. A pump adds a head $h_p$ to the left-hand side; a turbine removes a head $h_t$ from it. How to estimate $h_L$ for pipes is the subject of Chapter 5. The **power** carried by a stream with total head $H$ (relative to the datum) is $P = \gamma Q H$, and the power a pump delivers to the water is $\gamma Q h_p$.

## 3.4 Energy and hydraulic grade lines

The energy equation is easiest to understand as a drawing along the pipe. The **energy grade line** (EGL) is a line at height $H = p/\gamma + V^2/2g + z$ above the datum. The **hydraulic grade line** (HGL) is a line at the piezometric head $p/\gamma + z$, which is the level to which water would rise in piezometers along the pipe. The two lines are separated by the velocity head $V^2/2g$.

A few rules make these lines easy to draw:

- The EGL can only **fall** in the direction of flow, except at a pump. Its slope is the rate of friction loss.
- At a reservoir surface, the EGL and HGL both lie at the water surface ($p = 0$, $V \approx 0$).
- Where the pipe diameter falls, $V^2/2g$ rises and the HGL drops further below the EGL.
- At a free outlet to the atmosphere, the HGL passes through the centre of the jet, since $p = 0$ there.
- Where the HGL falls **below the pipe**, the pressure is below atmospheric. If the absolute pressure there approaches the vapour pressure, cavitation occurs.

Figure 3.1 shows a reservoir feeding a horizontal pipe that contracts once and ends in a nozzle. Friction losses have been estimated with the methods of Chapter 5 (Darcy friction factor $f = 0.02$, entrance loss $0.5\,V^2/2g$, contraction loss $0.25\,V^2/2g$).

```{code-cell} python
:tags: [hide-input]
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 11})

ink, ink2, c1 = "#0b0b0b", "#52514e", "#2a78d6"
c2, c3 = "#d9822b", "#3a9a5b"
rho, g = 998.2, 9.81
gamma = rho * g

H0, f = 10.0, 0.02                          # reservoir level above pipe axis (m)
D1, L1, D2, L2, Dj = 0.20, 40.0, 0.15, 30.0, 0.075
a1, a2, aj = (np.pi * np.array([D1, D2, Dj]) ** 2 / 4)
# H0 = Vj²/2g + losses; write all in terms of Q
coef = (1 / aj**2 + (0.5 + f * L1 / D1) / a1**2 + (0.25 + f * L2 / D2) / a2**2) / (2 * g)
Q = np.sqrt(H0 / coef)
V1, V2, Vj = Q / a1, Q / a2, Q / aj
hv1, hv2, hvj = V1**2 / (2*g), V2**2 / (2*g), Vj**2 / (2*g)

E0 = H0 - 0.5 * hv1                        # after entrance
E1 = E0 - f * L1 / D1 * hv1                # end of pipe 1
E2 = E1 - 0.25 * hv2                       # after contraction
E3 = E2 - f * L2 / D2 * hv2                # at nozzle
xE = [-8, 0, 0, L1, L1, L1 + L2, L1 + L2 + 3]
yE = [H0, H0, E0, E1, E2, E3, E3]
xH = [-8, 0, 0, L1, L1, L1 + L2, L1 + L2, L1 + L2 + 3]
yH = [H0, H0, E0 - hv1, E1 - hv1, E2 - hv2, E3 - hv2, 0, 0]

fig, ax = plt.subplots(figsize=(6.5, 4.0))
ax.fill_between([-8, 0], -1.0, H0, color="#e8f1fb")
ax.plot([-8, 0], [H0, H0], color=c1, lw=1.2)
ax.plot([0, 0], [-1.0, H0 + 1], color=ink, lw=2)
ax.plot([0, L1], [0, 0], color=ink2, lw=7, solid_capstyle="butt")
ax.plot([L1, L1 + L2], [0, 0], color=ink2, lw=5, solid_capstyle="butt")
ax.plot([L1 + L2, L1 + L2 + 1], [0, 0], color=ink2, lw=2.5)
ax.plot(xE, yE, color=c2, lw=2)
ax.plot(xH, yH, color=c1, lw=2)
ax.text(18, E0 - 0.2, "EGL", color=c2, va="bottom", fontweight="bold")
ax.text(18, E0 - hv1 - 0.6, "HGL", color=c1, va="top", fontweight="bold")
xm = L1 + L2 - 4
hgl_xm = E2 - hv2 - f * (xm - L1) / D2 * hv2
ax.annotate("", xy=(xm, hgl_xm), xytext=(xm, 0),
            arrowprops=dict(arrowstyle="<->", color=ink2, lw=0.8))
ax.text(xm - 1, 1.6, r"$p/\gamma$", color=ink2, ha="right")
ax.annotate("", xy=(L1 + L2 + 2.5, E3), xytext=(L1 + L2 + 2.5, 0),
            arrowprops=dict(arrowstyle="<->", color=ink2, lw=0.8))
ax.text(L1 + L2 + 3.5, E3 / 2, r"$V_j^2/2g$", color=ink2, va="center")
ax.text(5, -0.4, "datum = pipe axis", color=ink2, va="top")
ax.set_xlim(-9, 82); ax.set_ylim(-1.5, 11.5)
ax.set_xlabel("Distance along pipe (m)", color=ink2)
ax.set_ylabel("Head (m)", color=ink2)
ax.grid(True, color="#e4e3df"); ax.set_axisbelow(True)
for s in ("top", "right"): ax.spines[s].set_visible(False)
plt.tight_layout(); plt.show()
```

*Figure 3.1. Energy and hydraulic grade lines for a reservoir, a 200 mm pipe, a 150 mm pipe and a 75 mm nozzle discharging to the atmosphere. The EGL falls by local losses at the entrance and the contraction and slopes down by friction; at the nozzle all of the remaining head becomes the velocity head of the jet.*

## 3.5 Applications of the Bernoulli equation

**Pitot tube.** A tube bent to face the flow brings the water in front of it to rest. The point at the tube mouth is a **stagnation point**, where all the velocity head has been converted to pressure head. Applying Bernoulli between a point upstream (velocity $V$, pressure $p$) and the stagnation point (velocity 0, pressure $p_s$) at the same level,

$$
\frac{p_s - p}{\gamma} = \frac{V^2}{2g}
\quad\Rightarrow\quad
V = \sqrt{2g\, \Delta h},
$$

where $\Delta h$ is the difference between the stagnation head and the static head. A **Pitot-static tube** measures both heads in one instrument, the static pressure through small holes in the side of the tube. It gives a point velocity, not the mean velocity.

**Venturi meter.** A **Venturi meter** is a short converging section, a throat, and a gradual diverging section built into a pipe. Between the inlet (section 1, area $A_1$) and the throat (section 2, area $A_2$), continuity and Bernoulli give

$$
\begin{aligned}
\Delta h &= \left(\frac{p_1}{\gamma}+z_1\right)-\left(\frac{p_2}{\gamma}+z_2\right) \\
&= \frac{V_2^2 - V_1^2}{2g}.
\end{aligned}
$$

Using $V_1 = V_2 A_2/A_1$ and solving for the discharge,

$$
Q = C_d \, A_2 \sqrt{\frac{2g\, \Delta h}{1 - (A_2/A_1)^2}} .
$$

The **discharge coefficient** $C_d$ corrects for the small friction loss between the sections; for a well-made Venturi meter it is about 0.97–0.99. The head difference $\Delta h$ is the difference in piezometric head, which is exactly what a differential manometer measures: $\Delta h = R(S_m - 1)$ (Chapter 2). So the meter reading does not depend on whether the pipe is horizontal or inclined. The diverging section recovers most of the pressure, so the overall loss across a Venturi meter is small. The throat pressure must stay well above the vapour pressure to avoid cavitation.

Because $\Delta h \propto Q^2$, a Venturi meter is sensitive at high flows but reads poorly at low flows: halving the discharge divides the reading by four (Figure 3.2, Table 3.1).

```{code-cell} python
:tags: [hide-input]
Cd, Dv = 0.98, 0.30
A1v = np.pi * Dv**2 / 4
dh = np.linspace(0, 5, 200)                 # piezometric head difference (m)
fig, ax = plt.subplots(figsize=(6.5, 3.8))
for dt, col in [(0.10, c2), (0.15, c1), (0.20, c3)]:
    A2v = np.pi * dt**2 / 4
    Qv = Cd * A2v * np.sqrt(2 * g * dh / (1 - (A2v / A1v)**2))
    ax.plot(dh, Qv * 1e3, color=col, lw=2)
    ax.text(5.1, Qv[-1] * 1e3, f"throat {dt*1e3:.0f} mm", color=col, va="center")
ax.set_xlabel(r"Head difference $\Delta h$ (m of water)", color=ink2)
ax.set_ylabel("Discharge $Q$ (L s⁻¹)", color=ink2)
ax.set_xlim(0, 6.6); ax.set_ylim(0, 360)
ax.grid(True, color="#e4e3df"); ax.set_axisbelow(True)
for s in ("top", "right"): ax.spines[s].set_visible(False)
plt.tight_layout(); plt.show()
```

*Figure 3.2. Calibration curves of Venturi meters in a 300 mm pipe, with $C_d = 0.98$. The discharge grows with the square root of the head difference, so a smaller throat gives a larger reading for the same flow.*

```{code-cell} python
:tags: [remove-input]
from IPython.display import Markdown
dt = 0.15
A2v = np.pi * dt**2 / 4
rows = []
for Qi in [0.02, 0.05, 0.08, 0.10, 0.15, 0.20]:
    V1i, V2i = Qi / A1v, Qi / A2v
    dhi = (Qi / (Cd * A2v))**2 * (1 - (A2v / A1v)**2) / (2 * g)
    rows.append(f"| {Qi:.2f} | {V1i:.2f} | {V2i:.2f} | {dhi:.3f} | {dhi/12.6*1e3:.0f} |")
Markdown("| $Q$ (m³ s⁻¹) | $V_1$ (m s⁻¹) | $V_2$ (m s⁻¹) | $\\Delta h$ (m) | $R$ (mm Hg) |\n"
         "|---|---|---|---|---|\n" + "\n".join(rows))
```

*Table 3.1. Readings of a Venturi meter with a 150 mm throat in a 300 mm pipe ($C_d = 0.98$). $R$ is the reading of a water–mercury manometer, $R = \Delta h/12.6$.*

**Orifices and jets.** An **orifice** is an opening in the wall or base of a tank through which water flows out as a jet. Applying Bernoulli from the free surface (at height $h$ above the orifice, $V \approx 0$, $p = 0$) to the jet ($p = 0$) gives the **Torricelli** velocity,

$$
V_{\text{ideal}} = \sqrt{2gh}.
$$

Two effects reduce the real discharge. First, the streamlines approaching a sharp-edged orifice cannot turn the corner instantly, so the jet keeps contracting for a short distance after it leaves the opening. The narrowest section is the **vena contracta**, where the streamlines are parallel and the pressure is atmospheric. The ratio of its area to the orifice area is the **coefficient of contraction** $C_c$, about 0.62 for a sharp-edged circular orifice. Second, a little friction reduces the velocity; the ratio of actual to ideal velocity is the **coefficient of velocity** $C_v$, about 0.97–0.99. The discharge is then

$$
Q = C_d\, a \sqrt{2gh}, \qquad C_d = C_c C_v ,
$$

where $a$ is the orifice area and $C_d \approx 0.60$–0.62 for a sharp-edged orifice. A well-rounded entrance removes most of the contraction and gives $C_d$ close to 0.97.

**Draining a tank.** When a tank of plan area $A_T$ drains through an orifice, the head $h$ falls with time and so does the outflow. Continuity gives $A_T\, dh/dt = -C_d a \sqrt{2gh}$. Separating variables and integrating from $h_1$ to $h_2$,

$$
t = \frac{2 A_T}{C_d\, a \sqrt{2g}} \left(\sqrt{h_1} - \sqrt{h_2}\right).
$$

The time to empty the tank is found with $h_2 = 0$. Since the outflow slows as the head falls, the last part of the water takes the longest to drain (Figure 3.3).

```{code-cell} python
:tags: [hide-input]
DT, h0, Cdo = 2.0, 3.0, 0.61
AT = np.pi * DT**2 / 4
fig, ax = plt.subplots(figsize=(6.5, 3.6))
for do, col in [(0.040, c2), (0.050, c1), (0.060, c3)]:
    ao = np.pi * do**2 / 4
    K = Cdo * ao * np.sqrt(2 * g) / (2 * AT)
    tE = np.sqrt(h0) / K
    t = np.linspace(0, tE, 200)
    h = (np.sqrt(h0) - K * t)**2
    ax.plot(t / 60, h, color=col, lw=2, label=f"{do*1e3:.0f} mm orifice")
ax.legend(frameon=False, loc="upper right")
ax.set_xlabel("Time (min)", color=ink2)
ax.set_ylabel("Water depth $h$ (m)", color=ink2)
ax.set_xlim(0, 60); ax.set_ylim(0, 3.2)
ax.grid(True, color="#e4e3df"); ax.set_axisbelow(True)
for s in ("top", "right"): ax.spines[s].set_visible(False)
plt.tight_layout(); plt.show()
```

*Figure 3.3. Depth against time for a cylindrical tank 2 m in diameter, initially 3 m deep, draining through sharp-edged orifices of three diameters ($C_d = 0.61$). The level falls fast at first and slowly near the end.*

## 3.6 The momentum equation

Bernoulli's equation gives velocities and pressures, but it cannot give the **forces** that flowing water exerts on a structure. For those we use Newton's second law in the form suited to a fluid: the **momentum equation**.

We draw a **control volume**, a fixed region of space through which water flows, for example the water inside a pipe bend between two sections. For steady flow, with water entering at section 1 and leaving at section 2 with uniform velocities,

$$
\sum \mathbf{F} = \rho Q \left(\mathbf{V}_2 - \mathbf{V}_1\right).
$$

The left side is the sum of all external forces **on the water** in the control volume. The right side is the rate at which momentum leaves the control volume minus the rate at which it enters, the **momentum flux**. Both sides are vectors, so we write the equation separately for each direction:

$$
\sum F_x = \rho Q \left(V_{2x} - V_{1x}\right),
$$

$$
\sum F_y = \rho Q \left(V_{2y} - V_{1y}\right).
$$

The external forces usually include: the **pressure forces** $p_1 A_1$ and $p_2 A_2$ on the end sections, each pushing **into** the control volume; the **weight** of the water, if the flow has a vertical component; and the unknown **force from the walls**, $R$, which is the force we want. The force of the water on the structure is equal and opposite to $R$. A good procedure is:

1. Sketch the control volume and choose $x$ and $y$ axes.
2. Find velocities from continuity and pressures from Bernoulli.
3. Show every force on the water, with a positive direction assumed for $R_x$ and $R_y$.
4. Write the momentum equation in each direction and solve.
5. Reverse the signs to get the force on the structure.

**Nozzle.** At the end of a hose or a fire nozzle, the water accelerates from $V_1$ to the jet velocity $V_2$, which requires a large pressure $p_1$ upstream. The nozzle pushes back on the water, so the water pushes the nozzle forward: the joint between the hose and the nozzle is in tension (Example 3.4).

**Pipe bend.** In a horizontal bend that turns the flow through an angle $\theta$, with the same diameter, velocity $V$ and pressure $p$ at both ends (no loss), the force of the water on the bend is

$$
F = 2\left(pA + \rho Q V\right)\sin\frac{\theta}{2},
$$

directed outward along the bisector of the bend. In most water mains the pressure term $pA$ is much larger than the momentum term $\rho QV$. These forces are resisted by **thrust blocks**, concrete blocks cast against undisturbed soil at bends, tees and dead ends (Figure 3.4, Table 3.2).

```{code-cell} python
:tags: [hide-input]
Db, Qb, pb = 0.30, 0.20, 100e3
Ab = np.pi * Db**2 / 4
Vb = Qb / Ab
th = np.radians(np.linspace(0, 180, 181))
Fp = pb * Ab                 # pressure part
Fm = rho * Qb * Vb           # momentum part
fig, ax = plt.subplots(figsize=(6.5, 3.6))
ax.plot(np.degrees(th), 2 * (Fp + Fm) * np.sin(th / 2) / 1e3, color=c1, lw=2)
ax.plot(np.degrees(th), 2 * Fm * np.sin(th / 2) / 1e3, color=c2, lw=2, ls="--")
ax.text(95, 12.0, "total force", color=c1)
ax.text(95, 1.9, r"momentum part $2\rho QV\sin(\theta/2)$", color=c2)
ax.set_xlabel(r"Bend angle $\theta$ (degrees)", color=ink2)
ax.set_ylabel("Force on bend (kN)", color=ink2)
ax.set_xlim(0, 180); ax.set_ylim(0, 17); ax.set_xticks(range(0, 181, 30))
ax.grid(True, color="#e4e3df"); ax.set_axisbelow(True)
for s in ("top", "right"): ax.spines[s].set_visible(False)
plt.tight_layout(); plt.show()
```

*Figure 3.4. Resultant force on a horizontal bend in a 300 mm pipe carrying 0.2 m³ s⁻¹ at 100 kPa gauge pressure. The pressure term dominates; the momentum flux contributes only about 7% of the total.*

```{code-cell} python
:tags: [remove-input]
rows = []
for d in [22.5, 45, 60, 90, 135, 180]:
    t = np.radians(d)
    Fx = (Fp + Fm) * (1 - np.cos(t))
    Fy = (Fp + Fm) * np.sin(t)
    rows.append(f"| {d:g} | {Fx/1e3:.2f} | {Fy/1e3:.2f} | {np.hypot(Fx, Fy)/1e3:.2f} |")
Markdown("| Bend angle (°) | $F_x$ (kN) | $F_y$ (kN) | $F$ (kN) |\n"
         "|---|---|---|---|\n" + "\n".join(rows))
```

*Table 3.2. Components of the force on the bend of Figure 3.4. The flow enters along $+x$; $F_x$ acts along the inflow direction and $F_y$ acts away from the side to which the flow turns.*

**Jets on plates and vanes.** A jet of velocity $V$ and discharge $Q$ that strikes a fixed flat plate normal to it spreads out sideways, losing all of its momentum in the original direction. The force on the plate is $F = \rho Q V$. A fixed curved **vane** that turns the jet through an angle $\theta$ without loss of speed receives a force $\rho Q V(1 - \cos\theta)$ in the jet direction, which reaches $2\rho Q V$ when the jet is turned right back ($\theta = 180°$). This is the principle of the Pelton wheel. If the plate or vane moves away from the jet at speed $u$, the relative velocity $V - u$ replaces $V$ in these expressions, and the water reaching it per second falls to $\rho a (V - u)$ for a single plate, where $a$ is the jet area.

The same equation is used in open channels to find the force on a sluice gate and, in Chapter 6, the depths before and after a hydraulic jump.

## 3.7 Kinetic energy and momentum correction coefficients

In Sections 3.3–3.6 we used the mean velocity $V$ as if every particle had that velocity. Because the velocity head and the momentum flux depend on $u^2$ and $u$ weighted by the flow, a non-uniform profile carries more energy and momentum than the mean velocity suggests. We correct for this with two coefficients:

$$
\alpha = \frac{1}{A}\int_A \left(\frac{u}{V}\right)^3 dA ,
$$

$$
\beta = \frac{1}{A}\int_A \left(\frac{u}{V}\right)^2 dA .
$$

The **kinetic energy correction coefficient** $\alpha$ multiplies $V^2/2g$ in the energy equation, and the **momentum correction coefficient** $\beta$ multiplies $\rho Q V$ in the momentum equation. For uniform flow both equal 1. For laminar flow in a pipe (a parabolic profile), $\alpha = 2$ and $\beta = 4/3$. For turbulent pipe flow, whose profile is much flatter, $\alpha \approx 1.03$–1.10 and $\beta \approx 1.01$–1.04, and in natural rivers with flood plains they can be considerably larger. Since turbulent flow is the usual case in hydraulic engineering, we normally take $\alpha = \beta = 1$.

## 3.8 Worked examples

**Example 3.1 — Pitot tube in a stream.** A Pitot tube is held in a river, facing the flow, at a point below the surface. The water in its vertical limb rises 0.12 m above the river surface. Find the velocity at that point.

```{code-cell} python
dh = 0.12
V = np.sqrt(2 * g * dh)
print(f"Point velocity V = {V:.2f} m/s")
```

The static head at the tube mouth equals the local depth below the surface, so the 0.12 m rise is the velocity head. A velocity of 1.53 m s⁻¹ is typical of a river in moderate flood.

**Example 3.2 — Draining a tank.** A cylindrical tank 2 m in diameter holds water 3 m deep. It drains through a sharp-edged orifice of 50 mm diameter in its base ($C_d = 0.61$). How long does it take for the level to fall to 1 m, and how long to empty?

```{code-cell} python
AT = np.pi * 2.0**2 / 4
ao = np.pi * 0.050**2 / 4
Cdo = 0.61
c = 2 * AT / (Cdo * ao * np.sqrt(2 * g))
t31 = c * (np.sqrt(3.0) - np.sqrt(1.0))
t10 = c * np.sqrt(1.0)
print(f"3 m -> 1 m: {t31:.0f} s = {t31/60:.1f} min")
print(f"1 m -> 0 m: {t10:.0f} s = {t10/60:.1f} min")
print(f"Total time to empty: {(t31 + t10)/60:.1f} min")
```

The first 2 m drain in 14.4 min, but the last metre takes 19.7 min, because the outflow falls with the square root of the head.

**Example 3.3 — Venturi meter.** A Venturi meter with a 150 mm throat is fitted in a 300 mm water main. A mercury manometer ($S_m = 13.6$) between the inlet and the throat reads $R$ = 0.25 m. Take $C_d = 0.98$. Find the discharge.

```{code-cell} python
Sm, R, Cd = 13.6, 0.25, 0.98
A1 = np.pi * 0.30**2 / 4
A2 = np.pi * 0.15**2 / 4
dh = R * (Sm - 1)                              # m of water
Q = Cd * A2 * np.sqrt(2 * g * dh / (1 - (A2 / A1)**2))
print(f"Head difference = {dh:.2f} m of water")
print(f"Q = {Q:.4f} m3/s = {Q*1e3:.0f} L/s")
```

The meter reads 3.15 m of water head and the discharge is about 141 L s⁻¹. The answer would be the same if the meter were inclined, because the manometer measures the difference in piezometric head.

**Example 3.4 — Force on a hose nozzle.** A hose of 80 mm internal diameter ends in a nozzle with a 25 mm outlet. The discharge is 10 L s⁻¹ and the jet leaves to the atmosphere. Ignoring losses and the weight of the water, find the pressure at the nozzle inlet and the force on the joint between hose and nozzle.

```{code-cell} python
Q, D1, D2 = 0.010, 0.080, 0.025
A1, A2 = np.pi * D1**2 / 4, np.pi * D2**2 / 4
V1, V2 = Q / A1, Q / A2
p1 = rho * (V2**2 - V1**2) / 2              # Bernoulli, p2 = 0, same level
Rx = rho * Q * (V2 - V1) - p1 * A1          # force of nozzle on water (+ = downstream)
print(f"V1 = {V1:.2f} m/s, V2 = {V2:.2f} m/s, p1 = {p1/1e3:.0f} kPa")
print(f"Force of nozzle on water = {Rx:.0f} N (negative = upstream)")
print(f"Force of water on nozzle = {-Rx:.0f} N downstream")
```

The water pushes the nozzle downstream with about 850 N, roughly the weight of an 86 kg person. The coupling must hold this force in tension, and the firefighter must resist the jet's reaction of $\rho Q V_2 \approx 203$ N.

**Example 3.5 — Reducing bend.** A horizontal 45° bend reduces a pipe from 400 mm to 300 mm diameter. It carries 0.30 m³ s⁻¹ of water and the gauge pressure at the inlet is 150 kPa. Ignoring losses, find the components and resultant of the force of the water on the bend. The flow enters along $+x$ and leaves at 45° towards $+y$.

```{code-cell} python
D1, D2, Q, p1, th = 0.40, 0.30, 0.30, 150e3, np.radians(45)
A1, A2 = np.pi * D1**2 / 4, np.pi * D2**2 / 4
V1, V2 = Q / A1, Q / A2
p2 = p1 + rho * (V1**2 - V2**2) / 2
# Momentum on the water: p1A1 - p2A2 cos(th) + Rx = rho Q (V2 cos(th) - V1)
Rx = rho * Q * (V2 * np.cos(th) - V1) - p1 * A1 + p2 * A2 * np.cos(th)
Ry = rho * Q * (V2 * np.sin(th)) + p2 * A2 * np.sin(th)
Fx, Fy = -Rx, -Ry                            # force of water on the bend
print(f"V1 = {V1:.2f} m/s, V2 = {V2:.2f} m/s, p2 = {p2/1e3:.1f} kPa")
print(f"Force on bend: Fx = {Fx/1e3:.2f} kN, Fy = {Fy/1e3:.2f} kN")
print(f"Resultant = {np.hypot(Fx, Fy)/1e3:.2f} kN at {np.degrees(np.arctan2(Fy, Fx)):.1f}° from +x")
```

The bend is pushed forward (along $+x$) and away from the direction of turning (towards $-y$), with a resultant of about 14.0 kN. As in Figure 3.4, almost all of it comes from the pressure terms, so the thrust block must be designed for the highest pressure the main can see, including test pressure.

## 3.9 Common mistakes

- **Applying Bernoulli across a pump, a turbine or a long pipe.** Bernoulli's equation assumes no energy is added or lost. Use the energy equation with $h_p$, $h_t$ and $h_L$ instead.
- **Mixing gauge and absolute pressure in one equation.** Either is fine, but use the same reference at both sections. A free jet or open water surface has $p = 0$ gauge.
- **Forgetting the pressure forces in the momentum equation.** In pipes, $p_1A_1$ and $p_2A_2$ are usually far larger than the momentum flux terms.
- **Confusing the force on the water with the force on the structure.** The momentum equation gives the force of the walls on the water; the force on the structure is its opposite.
- **Using the orifice area instead of the vena contracta area.** The jet is narrower than the orifice. Use $C_d$, not 1, unless the entrance is well rounded.
- **Reading a Venturi manometer as a pressure difference.** The manometer gives the difference in piezometric head, $R(S_m - 1)$, which already includes any change in elevation.

## 3.10 Exercises

1. Water flows through a pipe that branches into two pipes of 150 mm and 100 mm diameter. The main pipe, 200 mm in diameter, carries 60 L s⁻¹, and the velocity in the 150 mm branch is 2.0 m s⁻¹. Find the velocity in the 100 mm branch.
2. A pipe rises 4 m between section 1 (diameter 300 mm, pressure 200 kPa) and section 2 (diameter 200 mm). The discharge is 0.1 m³ s⁻¹ and the head loss between the sections is 1.5 m. Find the pressure at section 2 and sketch the EGL and HGL.
3. Show that the Venturi meter equation can be written as $Q = C_d A_1 A_2 \sqrt{2g\Delta h}/\sqrt{A_1^2 - A_2^2}$. What happens to the reading if $A_2$ approaches $A_1$?
4. A rectangular tank of plan area 4 m² is fed with a constant inflow of 5 L s⁻¹ and drains through an orifice with $C_d a = 0.002$ m². Find the steady water depth in the tank.
5. A horizontal jet of water, 40 mm in diameter, strikes a fixed vertical plate at 15 m s⁻¹. Find the force on the plate. What is the force if the plate moves away from the jet at 5 m s⁻¹?
6. For the 90° bend of Table 3.2, find the force if the pressure is raised to a test pressure of 1.0 MPa. By how much does the momentum flux contribute now?

---

*Figures and tables in this chapter are generated by the Python code shown. Text and figures: © Daeha Kim, CC BY-NC 4.0.*
