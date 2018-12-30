## Tyrolian Computational Hydrodynamics

http://tycho-cfd.at/

TYCHO is a multidimensional (1D/2D and 3D) compressible hydrodynamics code written in C and parallelized with OpenMP. A Lagrangian remap version of the Piecewise Parabolic Method developed by Paul Woodward and Phil Colella (1984) is applied. The code is based on the freely available VH-1 Package.

A short video-tutorial on how to set-up a simulation for a fan-grill from a STL-file and how to use PARAVIEW/VISIT for the visualisation can be found here.
(The STL file was taken from http://www.thingiverse.com/thing:74964 - licensed under the Creative Commons license)

Find some examples on my channel at https://www.youtube.com/user/kapfology/videos?view=0


# TYCHO offers you:

..* Obstacle - Gas interaction
..* Heat-diffusion within obstacles
..*Thermal-exchange between gas and obstacles
..*Viscosity (Sutherland's law)
..*Gravity as a constant background field
..*Inclusion of marker-fields, which are advected with the velocity-field of the gas
..*Wind-emitters
..*Up to six different boundary conditions (easy to implement your own)
..*Can write four different file types (including VTK)
..*A graphical user interface to set up all simulation parameters
..*An easy 2D initial-condition generator to extract boundaries from Pixelgraphics
..*An easy 3D initial-condition generator to extract boundaries from STL and Point Data (such as 3D scans)
..*Sound-simulations for investigating dB-Maps with obstacle absorbing and reflecting sound-waves

Images of tychoGUI, tychoBCGEN and tychoBCG3D

TYCHO for Windows (v1.3.1) - the complete simulation package for Windows with installer

TYCHO (v1.3) - the simulation package for Linux (configure/make/make install)

tychoGUI (v1.0) - a graphical user-interface to make interaction with TYCHO easier for Linux (qmake-qt4/make)

tychoBCG3D (v0.2) - a 3D boundary condition generator for STL- and point data (e.g. 3D scans)

tychoBCGEN (v0.5) - a 2D boundary condition generator from pixel graphics for Linux (qmake-qt4/make)

The simulation package is focused on gas-obstacle interaction experiments and has special routines for obstacles in wind-streams and advection of marker fields. In addition momenta and their direction on obstacle-surfaces, thermal diffusion, and viscosity can be studied within in your simulation. Gravity is included as a constant background-field and a stratified atmosphere can be set up, if needed.

TYCHO is freely available to everyone. You are welcome to download it and do whatever you want with it. I would appreciate, however, if you would acknowledge the package in your publications/work and if you would send me information for what purpose you use the code. Keep in mind that this code does not come with any guarantee.

Here are some examples of TYCHO simulations. For some of them initial conditions can be downloaded.

Wind through a model city
Sound in a flat
Wind in Town
Pressure gradient in boiler
3D toy model simulation enterprise
Scramjet at MACH 2 temperature
Turbulent airflow around obstacles 
Soundspeaker simulation
Windtunnel temperature	Sydney Opera House	Kelvin Helmhotz Instability	Wind through rocky region
Gas in a Pipe
Thermal exchange of gas and a fractal shaped obstacle
Car in windtunnel
Airfoil

Windtunnel velocity
 
