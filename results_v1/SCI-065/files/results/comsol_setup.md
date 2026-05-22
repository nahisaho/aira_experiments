# COMSOL Multiphysics Setup Guide for a Perfusion Brain Organoid Bioreactor

## 1. Model Scope
Use a **3D model** for direct comparison with the OpenFOAM case files, or a **2D axisymmetric model** for rapid design screening.

### Recommended physics interfaces
- **Laminar Flow (spf)** for hydrodynamics in the vessel and porous basket region.
- **Transport of Diluted Species (tds)** for nutrient or oxygen transport once the flow field is converged.

## 2. Geometry Setup
1. Create a cylindrical vessel:
   - Diameter: **80 mm**
   - Height: **120 mm**
2. Add axial inlet and outlet ports:
   - Port diameter: **6 mm**
   - Inlet centered on the bottom face.
   - Outlet centered on the top face.
3. Add the internal organoid basket as a porous cylinder:
   - Diameter: **60 mm**
   - Height: **80 mm**
   - Position the basket from **z = 20 mm** to **z = 100 mm**.
4. If a 2D axisymmetric model is used, build the geometry in the **r-z plane** with radius 40 mm and height 120 mm.

## 3. Material Properties
Define the culture medium as a Newtonian fluid:
- Density: **1007 kg/m^3**
- Dynamic viscosity: **0.001 Pa·s**

For the basket porous region use:
- Porosity: **0.4**
- Permeability: **1e-10 m^2**
- Forchheimer coefficient: **1e5 1/m**

## 4. Laminar Flow (spf) Settings
1. Assign the vessel domain to **Laminar Flow**.
2. Set the outer vessel wall and basket support surfaces to **No Slip**.
3. Define the inlet as either:
   - **Volumetric Flow Rate** boundary, or
   - **Normal inflow velocity** computed from port area.
4. Define the outlet as **Pressure, no viscous stress** with **0 Pa** gauge pressure.
5. Add the basket domain as a **Porous Medium** (Brinkman equations recommended).
   - Enter porosity, permeability, and inertial resistance values.
6. If nutrient transport is included later, couple `tds` to the converged velocity field from `spf`.

## 5. Mesh Strategy
Use a physics-controlled **finer** mesh as a starting point, then refine manually.

### Manual mesh recommendations
- Create a **boundary layer mesh** near vessel walls and the basket interface.
- Use at least **5 boundary layers** with growth factor **1.2**.
- Target first layer thickness: **0.15–0.25 mm**.
- Refine the porous basket and inlet/outlet regions more strongly than the bulk fluid.
- For a 2D axisymmetric model, aim for roughly **150–250 axial elements** and **60–100 radial elements**.

## 6. Study Sequence
### Study 1: Stationary
Use a stationary study first to obtain converged velocity and pressure fields for each flow condition.

### Study 2: Time-Dependent
Add a time-dependent study if you want to evaluate startup transients, pulsatile perfusion, or species wash-in/out.
- Suggested initial range: **0 to 600 s**
- Suggested output step: **1 to 5 s**

## 7. Parametric Sweep
Create a parametric sweep over inlet flow rate:
- **0.5 mL/min**
- **1.0 mL/min**
- **2.0 mL/min**
- **5.0 mL/min**

Useful implementation options:
- Sweep on **volumetric flow rate** directly, or
- Sweep on inlet velocity converted from port area.

## 8. Solver Notes
- Use a segregated flow solver for the stationary study if the fully coupled solver is slow.
- Enable consistent stabilization only if convergence becomes difficult.
- Tighten relative tolerance to **1e-4** or better for comparing shear stress trends.

## 9. Post-Processing
Generate the following outputs:
- **Velocity streamlines** through the vessel and basket.
- **Velocity magnitude contours** in axial slices.
- **Pressure contours** along the vessel height.
- **Wall shear stress contours** on the vessel wall and basket surface.
- **Line plots** of shear stress along the basket height.
- **Species concentration contours** if `tds` is included.

## 10. Validation Checks
Before trusting the results, verify:
- Mass conservation between inlet and outlet.
- Mesh independence of peak velocity and wall shear stress.
- Sensitivity of basket pressure drop to permeability and porosity.
- That predicted shear levels stay within acceptable organoid handling limits.
