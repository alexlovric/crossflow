# Crossflow - FEM Multi-Physics Framework

## Dependencies

Crossflow requires the following dependencies:
- CMake 3.10 or higher
- C++11 compatible compiler
- Eigen3 (included in the repository)
- VTK 9.0 or higher (for visualization)

## Installation

1. **Install system dependencies**

   On Ubuntu/Debian:
   ```bash
   sudo apt-get update
   sudo apt-get install -y \
       build-essential \
       cmake \
       libvtk9-dev
   ```

   On Fedora:
   ```bash
   sudo dnf install -y \
       gcc-c++ \
       cmake \
       vtk-devel
   ```

2. **Build Crossflow**

   ```bash
   mkdir build && cd build
   cmake .. -DCMAKE_BUILD_TYPE=Release
   make -j$(nproc)
   ```

3. **Build and run examples**

   ```bash
   # Build all examples
   make -j$(nproc)
   
   # Run a specific example (e.g., poisson)
   cd examples/poisson
   ./crossflow_poisson
   ```

## Overview

Efficient FEM code capable of solving fluid dynamic problems for laminar viscous Newtonian fluids. The solvers included are:

* Standard incompressible
  * For coupled/decoupled velocity-pressure systems.
  * Decoupling based on a unique generalised-$\alpha$ projection scheme developed by [Lovrić et al. 2018](https://www.sciencedirect.com/science/article/pii/S0045782518302494?via%3Dihub), for segregating velocity and pressure fields. 

The framework allows for easy implementation of new linear/nonlinear PDEs and uses a Newton-Raphson nonlinear solver.

## Governing models

### Standard incompressible projection model

The **first** step of the generalised-$\alpha$ projection method with end-of-step elimination requires to compute the intermediate velocity $\tilde{\bm{u}}^{n+1}$ from

$$
    \rho\left(1 - \frac{\alpha_m}{\gamma} \right) \dot{\bm{u}}^n + \frac{\rho\alpha_m}{\gamma \Delta t} \left(\tilde{\bm{u}}^{n+1} - \tilde{\bm{u}}^n\right) + \rho\left(\tilde{\bm{u}}^{*,n+\alpha_f} \cdot \nabla \right)\tilde{\bm{u}}^{n+\alpha_f} 
    
    \\[5pt]+ \mu\Delta\tilde{\bm{u}}^{n+\alpha_f} + \nabla \left( \left(\alpha_f + \delta \right)p^n + \left(1 - \alpha_f - \delta \right) p^{n-1}\right) = \bm{f}^{n+\alpha_f}
    
    \\[12pt]\tilde{\bm{u}}^{n+1} \rvert_{\Gamma_g} = \bm{u}^{n+1}_{D}
$$

with

$$
	\tilde{\bm{u}}^{*,n+\alpha_f} = \alpha_f \tilde{\bm{u}}^{*,n+1} + \left(1-\alpha_f \right) \tilde{\bm{u}}^n,
	\\\tilde{\bm{u}}^{n+\alpha_f} = \alpha_f\tilde{\bm{u}}^{n+1} + \left(1 - \alpha_f \right)\tilde{\bm{u}}^n,
	\\\bm{f}^{n+\alpha_f} = \alpha_f\bm{f}^{n+1} + \left(1 - \alpha_f \right)\bm{f}^n,
$$

where the parameters $\alpha_m$, $\alpha_f$ and $\gamma$ are taken as

$$
	\alpha_m = \frac{1}{2} \ \frac{3 - \rho_\infty^h}{1 + \rho_\infty^h}
	\text{,}\qquad
	\alpha_f = \frac{1}{1 + \rho_\infty^h}	
	\text{,}\qquad
	\gamma = \frac{1}{2} + \alpha_m - \alpha_f.
$$

The convective velocity extrapolation is taken as

$$
    \tilde{\bm{u}}^{*,n+1} = 2\tilde{\bm{u}}^n - \tilde{\bm{u}}^{n-1}.
$$

In the **second** step the relations

$$
	\frac{\rho\alpha_m}{\gamma \Delta t} \left(\bm{u}^{n+1} - \tilde{\bm{u}}^{n+1} \right) + \nabla \left(\alpha_f p^{n+1} + \left(1 - \alpha_f - \delta \right) p^{n} \right) = \bm{0},
	\\\nabla \cdot \bm{u}^{n+1} = 0,
	\\\bm{u}^{n+1}\cdot\bm{n} \rvert_\Gamma = 0,
$$

give the following Poisson equation for the pressure $p^{n+1}$

$$
	-\frac{\rho \alpha_m}{\gamma \Delta t}\nabla\cdot\tilde{\bm{u}}^{n+1} + \Delta \left(\alpha_f p^{n+1} + \left(1 - \alpha_f - \delta \right) p^{n} \right) = 0.
$$

Finally the acceleration is computed from

$$
	\rho\bm{u}^{n+1} = \frac{\rho}{\gamma \Delta t} \left(\tilde{\bm{u}}^{n+1} - \tilde{\bm{u}}^n \right) 
    - \rho\frac{1 - \gamma}{\gamma} \dot{\bm{u}}^n 
	\\- \frac{1}{\alpha_m}\nabla \left(\alpha_f p^{n+1} + (1 - 2\alpha_f - \delta) p^{n}  - (1 - \alpha_f - \delta) p^{n-1}\right).
$$

## Code Directories
    ext  - external libraries
    src  - source files
    inc  - header files
    examples - example files

## Licence

*Yet to be written*


## Acknowledgements 

*Yet to be written*
