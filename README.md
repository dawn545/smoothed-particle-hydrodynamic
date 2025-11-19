# smoothed-particle-hydrodynamic
# SPH Fluid Simulation

A 2D fluid simulation using the Smoothed Particle Hydrodynamics (SPH) method, rendered with OpenGL.

## Prerequisites

To run this on **Windows**, you'll need:

- **MinGW** or **Visual Studio** for compiling C++ code.
- **GLUT** for OpenGL windowing and input handling.
- **Eigen3** for linear algebra operations.

### Install Dependencies

1. **MinGW** (for C++ compilation) - Download and install MinGW from [here](http://www.mingw.org/).
2. **GLUT**:
   - Download [freeglut](http://freeglut.sourceforge.net/) and set up the library.
   - Put the `freeglut.dll` in your system PATH or the same folder as your executable.
3. **Eigen3**:
   - Download the latest version of Eigen from [Eigen's website](https://eigen.tuxfamily.org/dox/GettingStarted.html).
   - Extract the files and place them in your project directory.

## Building the Project

1. Clone the repository or download the source files.
2. Open a command prompt or terminal in the project directory.
3. Compile the project using **MinGW** (or any compatible C++ compiler):

   ```bash
   g++ main.cpp SPHsolver.cpp renderer.cpp -o2 sph_simulation -lfreeglut -lopengl32 -lglu32
