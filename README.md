# Hele-Shaw



INTRODUCTION



Pattern formation in low $$Re$$ fluids, where its highly viscous and instabilities in these type when explored in 2 phase fluids, results in interesting finger like branches when a fluid is injected in a closely spaced plates also known as Hele Shaw cell. Pattern in these pressure driven flow appear as a result of fluid with lesser viscosity displacing a fluid with higher viscosity resulting in viscous fingering. This is as a result of "Saffman-Taylor instability" where the interface between these fluid becomes unstable due to the contrast of viscosity forming fingers.



These pattern are mostly studied in porous medium for example to understand the flow in spacings and in general these are determined by two factors i.e. the contrast in viscosity and miscibility.



MODELLING



Since this cell consist of closely spaced plates, when a fluid is injected into a fluid with different viscosity, branches appears which are inherently two dimensional and the velocity follows a parabolic profile like in a Poiseuille flow.



$$\\mathbf{u} = \\frac{- b^2}{12 \\mu} \\nabla p$$



where $$b$$ is the gap between the two plates and $$\\mu$$ is the viscosity of the fluid. This is a porous medium analogous to Darcy's law. Additionally with the flow purely governed by the Poisson equation of pressure $$\\nabla^2 p = 0$$ where at the interface the pressure across is set by surface tension. i.e. $$P\_{inside} - P\_{outside} = \\sigma \\kappa$$, where $$\\sigma$$ is the surface tension and $$\\kappa$$ is the interface curvature.



Following the footprint of the reference 1, it becomes quite possible to evaluate this dynamic interface using collocated grid points and to avoid adaptative grid refinement either in VOF or advance schemes like the level set methods or the front tracking methods. Since these is a complicated topological problem to handle in case of 3D. A rather subtle alternative to this is the "Phase field method", where the model is build on fluids free energy given as



$$f = \\frac{1}{q}\\alpha |\\nabla C|^q + \\beta \\psi(C)$$



Note here that the first term is the gradient of free energy, helpful to track the interface where as the second term represent the bulk fluid energy. C here is the measure of the phase and $$\\psi(C)$$ models the fluids immiscibility and has minima corresponding to the fluids two stable phase. Interface separating these two phases are $$O(\\sqrt{\\alpha / \\beta})$$ in width and have a surface tension of $$O(\\sqrt{\\alpha \\beta})$$ and the model depends on the coefficient $$q$$, $$\\alpha$$ and $$\\beta$$.



with parameters $$q=2$$, $$\\alpha \\sim O(\\epsilon)$$ and $$\\beta \\sim O(1/\\epsilon)$$, gives the phase field model.



Cahn-Hillard in their approach to lay a time dependent solution, gave an approximation of interfacial diffusion fluxes being proportional to chemical potential gradients and thus making it possible to solve the dynamics of this interface on the same grid while being to solve a transport for scalar $$c$$ also known as the Cahn-Hillard equation.



A general formulation to understand this system is given as



$$\\nabla \\cdot \\mathbf{u} = 0$$

$$\\nabla p = \\frac{-12\\eta}{b^2} \\mathbf{u} - \\epsilon \\rho \\nabla \\cdot \\left( \\left(\\nabla C\\right) + \\left(\\nabla C\\right)^T \\right)$$

$$\\rho \\left( \\partial\_t C + \\mathbf{u} \\cdot \\nabla C \\right) = \\alpha \\nabla ^2 \\mu$$



here $$\\epsilon$$ is coefficient of capillarity, $$\\alpha$$ is the coefficient of mobility and $$\\mu$$ is the chemical potential of the phases.



In diffuse-Interface framework, the viscosity of the fluid mixture ($$\\eta$$) is assume to be related to c with an exponential contrast constant $$R\_v = ln(\\eta\_2/\\eta\_1)$$



$$\\eta(c) = \\eta\_1 exp(R\_v(1-c))$$



where the injection strength is given as $$Q = \\pi(R\_c^2 - R\_0^2)/t\_c$$



Following the expression for phase potential $$\\mu$$ and the profile of the Helmholtz free energy $$f\_0$$ for a partially miscible interface is provided as 



$$\\mu(c) = \\partial\_c f\_0 - \\epsilon \\nabla^2 c$$

$$f\_0 = (c-c\_{s1})^2(c-c\_{s2})^2f^\*$$



where $$c\_{s1}$$ is the miscibility of the binary fluid and $$c\_{s2}$$ being its complementary.



Futhermore, a nondimensional form could be given as 



$$\\nabla \\cdot \\mathbf{u} = 0$$

$$\\nabla p = \\eta \\mathbf{u} - \\frac{C}{I} \\nabla \\cdot \\left( \\left(\\nabla C\\right) + \\left(\\nabla C\\right)^T \\right)$$

$$\\rho \\left( \\partial\_t C + \\mathbf{u} \\cdot \\nabla C \\right) = \\frac{1}{Pe} \\nabla ^2 \\mu$$



and



$$\\eta(c) = \\eta\_1 exp(R\_v(1-c))$$

$$\\mu(c) = \\partial\_c f - \\epsilon \\nabla^2 c$$

$$f = f\_0/f^\* = (c-c\_{s1})^2(c-c\_{s2})^2$$



where the nondimensional number are



Pe: Pelect number

$$Pe = \\frac{4 \\rho R\_c^2}{\\alpha f^\*t\_c}$$



Rv: contrast ratio

$$ln(\\eta\_1/\\eta\_2)$$



C: interface thickness

$$C = \\frac{\\epsilon}{4R\_c^2f^\*}$$



I: Injection parameter

$$I = \\frac{48 \\eta\_1 R\_c^2}{\\rho b^2 f^\* t\_c}$$



A rather interesting approach to deal with this problem is the vorticity stream function formulation where explicit linking for pressure is avoided and a coupled set of contraction appears in the expression of vorticity.



Following are the set



$$u = \\frac{\\partial \\psi}{\\partial y}; \\quad v = -\\frac{\\partial \\psi}{\\partial x}$$

$$\\omega = \\nabla \\times \\mathbf{u}$$

$$\\nabla^2 \\psi = -\\omega$$



where the vorticity is given as

$$\\omega = -R\_v \\left(u\\frac{\\partial c}{\\partial y} - v\\frac{\\partial c}{\\partial x}\\right) + \\frac{C}{\\eta I}\\left\[\\frac{\\partial c}{\\partial x} \\left( \\frac{\\partial^3 c}{\\partial x^2 \\partial y} + \\frac{\\partial^3 c}{\\partial y^3} \\right) - \\frac{\\partial c}{\\partial y} \\left( \\frac{\\partial^3 c}{\\partial x \\partial y^2} + \\frac{\\partial^3 c}{\\partial x^3} \\right)\\right]$$



$$\\mathbf{u} = u\_{rotational} + u\_{potential}$$



where $$u\_{rotational}$$ is obtained by solving the streamfunction and the potential is modelled by a smooth gaussian injection at the centre of the cell as 

$$u\_{potential} = \\frac{Q}{2\\pi r}\\left\[1 - exp\\left(\\frac{-4r^2}{R\_0^2}\\right)\\right]\\hat{r}$$



where $$\\hat{r}$$ is the unit normal. here the injection rate is modified as $$Q = \\pi(1-4R\_0^2)/4$$



Please see the references 1 and 2 for detailed understanding of the model.



FUNCTIONALITY



Numerical Methods

Spatial Discretization



Grid: Uniform Cartesian grid with Nx × Ny points

Domain: \[-Lx/2, Lx/2] × \[-Ly/2, Ly/2]

Grid spacing: dx = Lx/(Nx-1), dy = Ly/(Ny-1)



Finite Difference Schemes



Gradients: Central differences in interior, forward/backward at boundaries

Laplacian: 5-point stencil with ghost point method for boundary conditions

Boundary conditions: Ghost point method for Neumann, direct enforcement for Dirichlet



Time Integration



Scheme: Explicit Euler method

Stability: CFL-limited time step

Update: $$c^(n+1) = c^n + dt \* (advection + diffusion)$$



Poisson Solver



Method: Jacobi iteration

Convergence: $$L\_{\\infty}$$ norm with tolerance 1e-6

Maximum iterations: 1000 (configurable)



CahnHilliardViscousFingeringSolver

├── \_\_init\_\_(...)                    # Initialize solver

├── Physics Methods

│   ├── free\_energy\_derivative(c)    # ∂f/∂c computation

│   ├── viscosity(c)                 # η(c) = exp(Rv(1-c))

│   └── injection\_velocity()         # Radial injection field

├── Numerical Methods  

│   ├── apply\_neumann\_bc(field)      # ∂field/∂n = 0

│   ├── apply\_dirichlet\_bc(field)    # field = 0 on boundaries

│   ├── gradient\_x(field)            # ∂field/∂x

│   ├── gradient\_y(field)            # ∂field/∂y

│   ├── laplacian(field, bc\_type)    # ∇²field

│   └── solve\_poisson\_jacobi(rhs)    # ∇²φ = rhs

├── Evolution Methods

│   ├── compute\_chemical\_potential() # μ = ∂f/∂c - C∇²c

│   ├── compute\_vorticity\_source()   # ω calculation

│   ├── compute\_streamfunction\_...() # ψ and velocity

│   └── time\_step()                  # Advance one time step

├── Simulation Control

│   ├── initialize\_concentration()   # Set initial condition

│   └── run\_simulation(t\_final)      # Main simulation loop

└── Visualization

&nbsp;   ├── plot\_results()               # Static plots

&nbsp;   ├── plot\_snapshots(results)      # Export image sequence

&nbsp;   ├── animate\_field(...)           # Create field animations

&nbsp;   └── animate\_velocity(...)        # Create velocity animations



REFERENCES



\[1] Yuka F. Deki et al Numerical simulation of effects of phase separation on viscous fingering in radial Hele-Shaw flows. Journal of Fluid Mechanics, Vol. 1003, A12, 2025.



\[3] https://visualpde.com/nonlinear-physics/cahn-hilliard.html



\[4] https://web.tuat.ac.jp/~yamanaka/pcoms2019/Cahn-Hilliard-2d.html

