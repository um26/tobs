from ngsolve import GridFunction, LinearForm, BilinearForm
from ngsolve.la import BaseMatrix

def solve(bilinearForm :LinearForm, 
          linearForm : BilinearForm
          ) -> tuple[GridFunction, BaseMatrix]:
    """ Solve the linear system associated to the given bilinear and linear forms."""
    K = bilinearForm.Assemble().mat
    f = linearForm.Assemble().vec
    fes = bilinearForm.space
    sol = GridFunction(fes)
    Kinv = K.Inverse(freedofs = fes.FreeDofs())
    sol.vec.data = Kinv * f
    return sol, Kinv


from ngsolve import CF, cos, sin, pi
from ngsolve.webgui import Draw

def DrawMaterial(rho, theta, scene = None):
    """ plot the fibers along the parallel direction """
    mesh = rho.space.mesh
    if scene is None :
        scene = Draw(rho*CF((cos(theta+pi/2),sin(theta+pi/2))), mesh, vectors = { "grid_size":20}, 
                     settings = {"Objects" : { "Wireframe" : False } } , min=0, max=1)
        return scene
    else:
        scene.Redraw(rho*CF((cos(theta+pi/2),sin(theta+pi/2))), mesh, vectors = { "grid_size":20},
                     settings = {"Objects" : { "Wireframe" : False } }, min=0, max=1)


from ngsolve import CF, grad

def curl(v):
    """ 2D curl operator | scalar -> vector """
    R = CF(((0,1),(-1,0)), dims = (2,2))  # Rotation matrix of angle -pi/2
    return R * grad(v)

from ngsolve import H1, BilinearForm, LinearForm, dx

def solve_state(nu : GridFunction, j : float = 1e6):
    """ Solve the state equation for a given reluctivity field """
    mesh = nu.space.mesh
    fes = H1(mesh, order = 1, dirichlet = "dOmega")
    a, v = fes.TnT()    # define the trial and test functions 
    bf = BilinearForm(curl(v) * (nu * curl(a)) * dx)  
    lf = LinearForm(j * v * dx("Pp") - j * v * dx("Pm"))
    return solve(bf, lf)


#%% Post processing

from ngsolve import Integrate

def flux(sol : GridFunction, Ns : float = 100, Lz : float = 1,
         positive_coil = "Sp", negative_coil = "Sm"):
    """ Compute the flux in the secondary coil """
    mesh = sol.space.mesh
    Sp = Integrate(1, mesh, definedon = mesh.Materials(positive_coil), order = 1)
    Sm = Integrate(1, mesh, definedon = mesh.Materials(negative_coil), order = 1)
    aSp = Integrate(sol, mesh, definedon = mesh.Materials(positive_coil))
    aSm = Integrate(sol, mesh, definedon = mesh.Materials(negative_coil))
    return Ns * Lz * (aSp/Sp - aSm/Sm)