# web_physics

Interactive browser demos of my physics and chemistry codes, rewritten in JavaScript/WebGL so they run entirely client‑side (no installation needed). Each directory here is a small web app that showcases a stripped‑down version of the underlying C++/Python tools.

- **Live demos:** open the subfolders in a browser (or use the main `index.html` front page).
- **Full codes:** see linked GitHub projects for production‑grade simulators and scripts.

it is hosted here: **[https://prokophapala.github.io/web_physics/](https://prokophapala.github.io/web_physics/)**

Below, the demos are grouped into three themes.

## 1. AFM / STM Imaging and Pauli Master Equation

Browser ports of my scanning probe microscopy (SPM) and transport models: Probe‑Particle AFM/STM and a Pauli master‑equation dI/dV solver.

- **ppafm_web**, **ppafm_web_inline**  
  Interactive AFM image simulators based on the Probe‑Particle Model (PPM). They visualize frequency‑shift contrast of functionalized tips over molecules.
  - Upstream code: [PPAFM](https://github.com/Probe-Particle/ppafm)
  - Key refs:  
    - Hapala, P. et al., "Mechanism of high-resolution STM/AFM imaging with functionalized tips", *Phys. Rev. B* **90**, 085421 (2014). DOI: [10.1103/PhysRevB.90.085421](https://doi.org/10.1103/PhysRevB.90.085421)  
    - Hapala, P. et al., "Origin of High-Resolution IETS-STM Images of Organic Molecules with Functionalized Tips", *Phys. Rev. Lett.* **113**, 226101 (2014). DOI: [10.1103/PhysRevLett.113.226101](https://doi.org/10.1103/PhysRevLett.113.226101)  
    - Oinonen, N. et al., "Advancing Scanning Probe Microscopy Simulations: A Decade of Development in Probe-Particle Models", *Comput. Phys. Commun.* **305**, 109341 (2024). DOI: [10.1016/j.cpc.2024.109341](https://doi.org/10.1016/j.cpc.2024.109341)

- **ppstm_web**  
  STM / IETS simulator for flexible tips, using Chen’s approximation and PPM‑based tip models.
  - Upstream code: [PPSTM](https://github.com/Probe-Particle/PPSTM)
  - Key ref: Krejčí, O. et al., "Principles and simulations of high-resolution STM imaging with a flexible tip apex", *Phys. Rev. B* **95**, 045407 (2017). DOI: [10.1103/PhysRevB.95.045407](https://doi.org/10.1103/PhysRevB.95.045407)

- **pauli_web**  
  Pauli master‑equation (PME) demo for quantum transport and dI/dV maps in molecular assemblies.
  - Paper: Li, C. et al., "Negative differential conductance in triangular molecular assemblies", *arXiv:2508.05575* (2025). DOI: [10.48550/arXiv.2508.05575](https://doi.org/10.48550/arXiv.2508.05575)  
  - Related package: [QmeQ](https://github.com/gedaskir/qmeq), Kiršanskas, G. et al., "QmeQ 1.0: An open-source Python package for calculations of transport through quantum dot devices", *Comput. Phys. Commun.* **221**, 317–342 (2017). DOI: [10.1016/j.cpc.2017.07.024](https://doi.org/10.1016/j.cpc.2017.07.024)

---

## 2. Molecular Modeling and Engineering (GridFF, FireCore, MolGUI)

Lightweight molecular modeling tools, focusing on molecules on rigid substrates and CAD‑like molecular editing.

- **molgui_web**  
  Web molecular editor and viewer using force‑field parameters (MMParams) and WebGL rendering. Intended as a front‑end for on‑surface self‑assembly studies.
  - Related engines:  
    - [FireCore](https://github.com/ProkopHapala/FireCore) – QM/MM and on‑surface chemistry framework  
    - Grid‑accelerated force field (GridFF) for molecules on rigid substrates
  - Key refs:  
    - Mal, I. et al., "GridFF: Efficient Simulation of Organic Molecules on Rigid Substrates", *J. Chem. Theory Comput.* **21**, 12214–12226 (2025). DOI: [10.1021/acs.jctc.5c01223](https://doi.org/10.1021/acs.jctc.5c01223)  
    - Hapala, P. et al., "Polymer templates for nanofabrication", *ACS Nano* (2024). DOI: [10.1021/acsnano.3c10575](https://doi.org/10.1021/acsnano.3c10575)

---

## 3. Space‑Warfare Game and General Physics Sandbox

Gameplay‑oriented physics demos derived from my C++ **SimpleSimulationEngine**: orbital warfare, spacecraft design, plasma/MHD tests, and 3D mesh viewers.

- **mhd_demo**  
  Simple magnetohydrodynamics (MHD) testbed: fluid/plasma solver with interactive visualization.

- **spacecraft_editor**  
  Spacecraft design playground: parametric hulls, trusses, engine blocks and 3D previews, inspired by orbital‑warfare simulations.

- **OBJ_Loader**  
  Minimal OBJ mesh viewer used for quick inspection of 3D assets (e.g. turrets, ship parts) in the browser.

Underlying engine and game prototypes:

- [SimpleSimulationEngine](https://github.com/ProkopHapala/SimpleSimulationEngine) – C++ framework (with Python/Lua bindings) for physics, numerics and OpenGL visualization.
  - Orbital combat concept: [OrbitalWar wiki](https://github.com/ProkopHapala/SimpleSimulationEngine/wiki/OrbitalWar)
  - Background reading: [space_warfare encyclopedia](https://github.com/ProkopHapala/SimpleSimulationEngine/tree/master/encyclopedia/space_warfare)

---

## Repository Layout (this repo)

- **index.html** – simple front page linking to individual web demos.
- **ppafm_web/**, **ppafm_web_inline/** – AFM imaging demos (Probe‑Particle Model).
- **ppstm_web/** – STM / IETS demo.
- **pauli_web/** – Pauli master‑equation transport demo.
- **molgui_web/** – molecular editor / viewer.
- **mhd_demo/** – MHD / plasma demo.
- **spacecraft_editor/** – spacecraft design and visualization demo.
- **OBJ_Loader/** – OBJ model viewer.
- **common_js/**, **common_resources/** – shared JS utilities, shaders and data.

For background on the research and full‑scale codes, see `Supplement_info.md` or visit my GitHub profile: <https://github.com/ProkopHapala>.

---

## About the author

This project is developed by **Prokop Hapala**, computational physicist at the surface-science group of the Institute of Physics of the Czech Academy of Sciences in Prague. His work focuses on high‑resolution scanning probe microscopy, on‑surface chemistry, and computer‑aided design of molecular assemblies within the GAČR JUNIOR STAR project **CADTARSIS** (*Computer Aided Design of Templated Assembling, Replication and Synthesis on Ionic Substrates*).

In CADTARSIS he develops photosensitive polymer templates that encode structural information and drive deterministic bottom‑up assembly of molecular components, aiming at scalable nanofabrication of molecular electronic and photonic devices below the resolution limits of photolithography, as described in
"Computational Design of Photosensitive Polymer Templates To Drive Molecular Nanofabrication" (*ACS Nano* 2024, DOI [10.1021/acsnano.3c10575](https://doi.org/10.1021/acsnano.3c10575)) and on the
[CADTARSIS project page](https://www.fzu.cz/en/research/divisions-and-departments/division-2/department-23/research/computational-design-nanofabrication-molecular-computers).

- FZU profile: [Prokop Hapala](https://www.fzu.cz/en/people/ing-prokop-hapala-phd)  
- ORCID: [0000-0003-4807-0326](https://orcid.org/0000-0003-4807-0326)  
- GitHub: [github.com/ProkopHapala](https://github.com/ProkopHapala)