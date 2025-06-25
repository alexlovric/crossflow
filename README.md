# Crossflow - Multiphysics FEM Framework
[![Build Status](https://github.com/alov/crossflow/actions/workflows/build-and-test.yml/badge.svg)](https://github.com/alov/crossflow/actions)

Crossflow is a finite element method (FEM) framework designed for solving multiphysics problems. It provides a flexible and easily extendible platform for simulating fluid dynamics, structural mechanics, and coupled physics phenomena. Easily switch out the solver for something beefier if you want.

![Vortex Shedding](docs/vortex_shedding.gif)

<sup>Simulation of vortex shedding around a cylinder (from cylinder example).</sup>

## Dependencies

Crossflow requires the following dependencies:
- CMake 3.10 or higher
- C++17 compatible compiler
- Eigen3 (>=3.3)
- VTK (>=9.0)
- Google Test (for tests)

### Installation on Ubuntu/Debian

```bash
sudo apt-get update
sudo apt-get install -y \
    build-essential \
    cmake \
    libeigen3-dev \
    libvtk9-dev \
    libgtest-dev \
    libgmock-dev
```

### Installation on Fedora

```bash
sudo dnf install -y \
    gcc-c++ \
    cmake \
    eigen3-devel \
    vtk-devel \
    gtest-devel \
    gmock-devel
```

## Building

```bash
git clone https://github.com/yourusername/crossflow.git
cd crossflow
mkdir build
cd build
cmake ..
make
```

### Installing Crossflow as a library

To install Crossflow system-wide:

```bash
sudo make install
```

This will install the headers and configuration files to your system's default include directories.

## Using Crossflow in Your Project

To use Crossflow in your own project, create a `CMakeLists.txt` like this:

```cmake
project(YourProject)

# Find Crossflow
cmake_minimum_required(VERSION 3.10)
find_package(crossflow REQUIRED)

# Create your executable
add_executable(your_program main.cpp)

target_link_libraries(your_program PRIVATE crossflow)
```

## Running Examples

### Poisson Equation

```bash
cd build
make run_poisson
```

### Cylinder Flow

```bash
cd build
make run_cylinder
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
