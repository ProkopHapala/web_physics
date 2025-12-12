# [ppafm_web](cci:7://file:///home/prokop/git/web_physics/ppafm_web:0:0-0:0) overview (what is where, who does what)

## High-level architecture
- **CPU (JavaScript in [index.html](cci:7://file:///home/prokop/git/web_physics/ppafm_web/index.html:0:0-0:0))**:
  - Loads/parses molecule (`.xyz`)
  - Converts UI parameters into **shader uniforms**
  - Precomputes a few helper arrays (notably **Giessibl convolution weights**)
  - Triggers renders (single frame or animation loop)

- **GPU (fragment shader [PP_AFM_shader.glslf](cci:7://file:///home/prokop/git/web_physics/ppafm_web/PP_AFM_shader.glslf:0:0-0:0))**:
  - For **each pixel**: interprets that pixel as one lateral `(x,y)` probe position
  - Runs the **probe-particle relaxation + approach/oscillation loop**
  - Produces one scalar output per pixel (df / Fz / residual / iteration-count) mapped to grayscale

There are only 2 real “modules” here: [index.html](cci:7://file:///home/prokop/git/web_physics/ppafm_web/index.html:0:0-0:0) (UI + JS driver) and [PP_AFM_shader.glslf](cci:7://file:///home/prokop/git/web_physics/ppafm_web/PP_AFM_shader.glslf:0:0-0:0) (the simulation kernel).

# File responsibilities

## [ppafm_web/index.html](cci:7://file:///home/prokop/git/web_physics/ppafm_web/index.html:0:0-0:0) (everything on the CPU side)
### 1) UI + parameter collection
- HTML inputs define simulation parameters:
  - Imaging plane & view mapping: `inpZ`, `inpScale`, `inpCX`, `inpCY`
  - Probe parameters: `inpProbeQ`, `inpProbeR`, `inpProbeE`
  - Tip mechanics: `inpKLat`, `inpKRad`, `inpRTip`
  - Relaxation/integration knobs: `inpRelaxIters`, `inpDt`, `inpF2Conv`, `inpOscSteps`, `inpDz`, `inpOscAmp`, `inpPreRelax`, `inpRelaxSub`
  - Output selection: `inpRender` (df/Fz/etc), `inpAlgo` (declared, but see shader notes)

### 2) WebGL/Three.js setup
- Creates:
  - `THREE.WebGLRenderer`, `THREE.Scene`, `THREE.OrthographicCamera`
  - Fullscreen quad (`PlaneBufferGeometry(2,2)`) + `THREE.ShaderMaterial`
- Loads fragment shader source via `fetch('PP_AFM_shader.glslf')` and plugs it into `ShaderMaterial.fragmentShader`.

### 3) Uniforms definition and lifecycle
- Defines `uniforms` object (Three.js-style), e.g.:
  - `uResolution`, `uZPlane`, `uScale`, `uCenter`, `uContrast`, …
  - `uAtoms[256]` and `uREQK[256]` as JS arrays of `THREE.Vector4`
  - `uWeights[32]` as `Float32Array`
- Initializes `uAtoms/uREQK` with **256 dummy entries** to avoid undefined behavior in Three.js.

### 4) Parsing XYZ and preparing atom parameters ([applyXYZ](cci:1://file:///home/prokop/git/web_physics/ppafm_web/index.html:641:6-768:7))
[applyXYZ(text)](cci:1://file:///home/prokop/git/web_physics/ppafm_web/index.html:641:6-768:7) does:
- Parses XYZ supporting:
  - With comment line: `N`, comment, then `N` atoms
  - Without comment line: `N`, then `N` atoms
- Packs per-atom data into:
  - `uAtoms[i] = vec4(x,y,z,q_atom)` (q optional, defaults 0)
- Computes molecule bounds `molBounds` for auto-centering and top-Z detection.
- Applies optional **Z shift** (`chkShiftZ`): shifts so top-most atom is at `z=0`.

**Most important preparation step:** building `uREQK`:
- For each atom, JS computes a mixed interaction parameter vector:
  - `R0 = Rpp + Rii`
  - `E0 = sqrt(Epp * Eii)` (LB mixing)
  - `Qij = qProbe * q_atom`
  - `K_default = 1.5` (currently constant for all)
- Stores this in `uREQK[i] = vec4(R0, E0, Qij, K)`
- Pads both arrays out to exactly 256 entries (shader loops over max 256).

So: **the shader never needs element types**; it only sees already-mixed parameters.

### 5) Convolution weights ([calculateGiessiblWeights](cci:1://file:///home/prokop/git/web_physics/ppafm_web/index.html:504:8-523:9))
- [calculateGiessiblWeights(nSteps, dz)](cci:1://file:///home/prokop/git/web_physics/ppafm_web/index.html:504:8-523:9) precomputes a discrete weight vector used to convert `Fz(z)` samples into a `df`-like signal (Giessibl-style).
- Writes these into `uniforms.uWeights` (size 32, so `uOscSteps` is clamped to 32 in JS).

### 6) Updating the simulation parameters ([updateViewFromInputs](cci:1://file:///home/prokop/git/web_physics/ppafm_web/index.html:379:8-502:9))
- Reads all numeric inputs, sanitizes constraints (positive stiffness, dt, clamp steps).
- Computes actual `uZPlane` as:
  - `zTop + userZ` where `zTop` depends on `chkShiftZ` and molecule bounds
- Computes `uCenter` either:
  - molecule center + user offsets (if `chkCenterXY`)
  - or direct user coordinates
- Converts stiffness from N/m to eV/Å² by dividing by `16.02176634`
- Finally calls [renderOnce()](cci:1://file:///home/prokop/git/web_physics/ppafm_web/index.html:607:6-612:7).

### 7) Rendering control
- [renderOnce()](cci:1://file:///home/prokop/git/web_physics/ppafm_web/index.html:607:6-612:7) calls `renderer.render(scene, camera)`
- [animate()](cci:1://file:///home/prokop/git/web_physics/ppafm_web/index.html:614:6-622:7) calls [renderOnce()](cci:1://file:///home/prokop/git/web_physics/ppafm_web/index.html:607:6-612:7) in a `requestAnimationFrame` loop if `chkAnimate` is enabled.

## [ppafm_web/PP_AFM_shader.glslf](cci:7://file:///home/prokop/git/web_physics/ppafm_web/PP_AFM_shader.glslf:0:0-0:0) (everything on the GPU side)
This is the “simulation kernel”. It runs once per pixel.

### Inputs (uniforms)
- View mapping: `uResolution`, `uScale`, `uCenter`
- Scan height: `uZPlane`
- Atom data:
  - `uNumAtoms`
  - `uAtoms[i].xyz` positions, `.w` charge (not used directly in main kernel here)
  - `uREQK[i] = (R0,E0,Q,K)` used for Morse + Coulomb
- Tip model & relaxation:
  - `uKLat`, `uKRad`, `uRtip`
  - `uPreRelax`, `uRelaxSub`, `uOscSteps`, `uDz`, `uOscAmp`
  - convergence threshold `uF2Conv`
- Output: `uRenderMode`, `uContrast`
- Notes:
  - Relaxation uses a fixed-point / preconditioned update of the probe-particle position based on total force and an approximate inverse stiffness.
  - In high repulsion the sample potential curvature can dominate the effective stiffness, so naive preconditioned steps may overshoot unless stabilized.

### Relaxation stability notes (high repulsion)
- The core update uses a *preconditioned* step proposal `dN = Ftot * invK`, where `invK ~ (1/k_lat, 1/k_lat, 1/k_rad)`.
- This assumes the dominant curvature comes from the probe springs (`k_lat/k_rad`). In highly repulsive regions, the sample Hessian can dominate, so `dN` can become too large and cause divergence/artefacts.
- The issue is stronger for small `k_lat` because `invK.xy ~ 1/k_lat` amplifies lateral steps.
- Increasing `uRelaxSub` can make it worse because it applies an unstable update repeatedly.

To stabilize this, the shader supports a simple trust-region clamp (and a few experimental modes) controlled by `uAlgo` and a few step parameters.

**Recommended default (found empirically):**
- `Algo = 1` (trust region)
- `StepXY ≈ 0.10 Å`

### Relaxation algorithm modes (`uAlgo`)
- `uAlgo = 0`: current (legacy) update, no stabilization
- `uAlgo = 1`: trust-region clamp of the preconditioned step (caps lateral/vertical step sizes)
- `uAlgo = 2`: switch between preconditioned step and gradient descent step based on step magnitude threshold
- `uAlgo = 3`: gradient descent (bounded step) baseline

### Stabilization parameters
- `uStepMaxXY` (Å): trust region radius for lateral step per sub-iteration (used in `uAlgo>=1`)
- `uStepMaxZ`  (Å): trust region clamp for z step per sub-iteration
- `uStepSwitch` (Å): threshold on `|Ftot*invK|` for switching to GD in `uAlgo=2`
- `uGDStep` (Å): gradient descent step length (used in `uAlgo>=2`)

### Core physics routines
- `getMorseCoulombForce(pos, apos, REQK)`:
  - Morse (in a shifted squared form) + Coulomb term from `Q`
  - Returns force vector in eV/Å
- `computeSampleForce(pos)`:
  - Sums forces from all atoms up to `uNumAtoms`
- `computeTipForce(dpos)`:
  - Lateral Hooke spring: `Fxy = -uKLat * dxy`
  - Radial spring to enforce `|dpos| ≈ uRtip`: `-uKRad*(l-uRtip)` along `dpos`

### Main per-pixel algorithm (`main`)
For each pixel:
1. Map pixel to lateral position:
   - `uv = gl_FragCoord/uResolution → [-1,1]`
   - `posXY = uCenter + uv*uScale`
2. Initialize anchor + probe position:
   - `zStart = uZPlane + uOscAmp`
   - `anchor = (x,y,zStart)`
   - `pos = anchor`, then `pos.z -= uRtip` (equilibrium bond offset)
3. Pre-relax (far field):
   - up to `uPreRelax` iterations, update:
     - `Ftot = Fsamp(pos) + Ftip(pos-anchor)`
     - stop if `|Ftot|² < uF2Conv`
     - update `pos += Ftot * invK` where `invK ~ 1/k` (a fixed-point-like relaxation)
4. Approach/oscillation loop (`step < uOscSteps`):
   - For each step:
     - do `uRelaxSub` relaxation sub-iterations at current Z (same update scheme)
     - measure final `FsampFinal.z`
     - accumulate `dfAccum += Fz * uWeights[step]`
     - then advance down: `anchor.z -= uDz` and `pos.z -= uDz`
5. Output mapping:
   - `uRenderMode==0`: grayscale from `dfAccum`
   - `==1`: grayscale from relaxed `Fz`
   - `==2`: grayscale from residual force norm
   - `==3`: grayscale from iteration count fraction

So the shader is effectively doing a **“scan” image in one pass**: each fragment is one lateral grid point; the vertical relax+approach is done inside that fragment.

# Data flow summary (CPU → GPU)
- **XYZ text/file** → [applyXYZ()](cci:1://file:///home/prokop/git/web_physics/ppafm_web/index.html:641:6-768:7):
  - computes `uAtoms[256]`, `uREQK[256]`, `uNumAtoms`, `molBounds`
- **UI params** → [updateViewFromInputs()](cci:1://file:///home/prokop/git/web_physics/ppafm_web/index.html:379:8-502:9):
  - computes `uZPlane`, `uCenter`, `uScale`, stiffness, relax knobs, `uWeights`
- **Render**:
  - GPU runs fragment shader for every pixel, producing grayscale scalar field.

# Notable “boundaries” / good places to improve
- **Physics kernel changes** mostly go into [PP_AFM_shader.glslf](cci:7://file:///home/prokop/git/web_physics/ppafm_web/PP_AFM_shader.glslf:0:0-0:0).
- **Changing interaction parametrization** (element tables, mixing rules, per-element K, etc.) is currently in [applyXYZ()](cci:1://file:///home/prokop/git/web_physics/ppafm_web/index.html:641:6-768:7) (CPU side) and only affects the shader through `uREQK`.
- **Performance / resolution / multi-pass** is constrained by “one fragment = one full relaxation + approach” (expensive but simple).
