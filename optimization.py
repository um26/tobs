from ngsolve import GridFunction, LinearForm
from ngsolve.la import BaseMatrix

def solve_adjoint(state : GridFunction,     # state solution
                  rho : GridFunction,       # density
                  Kinv : BaseMatrix,        # decomposition of the state matrix
                  df : callable             # directional derivative of f, with signature df(state, rho, v) where v is the test function
                  ) -> GridFunction:
    """ Solve the linear system associated to the given bilinear and linear forms."""
    fes = state.space
    v = fes.TestFunction()
    f = LinearForm(df(state, rho, v)).Assemble().vec
    adjoint = GridFunction(fes)
    adjoint.vec.data = -1 * Kinv.T * f
    return adjoint

def project(x : GridFunction, 
            x_min : float = 0, 
            x_max : float = 1):
    x.vec.data = np.maximum(x_min, np.minimum(x.vec, x_max))
    return x


import numpy as np
from copy import copy
from ngsolve import Integrate, Norm
from ngsolve.webgui import Draw

def gradient_descent(state : callable,           # a_rho : GridFunction, Kinv : BaseMatrix = state(rho : GridFunction)
                     objective : callable,       # obj : float = objective(rho : GridFunction, a_rho : GridFunction)
                     d_objective : callable,     # riesz_repr : GridFunction = d_objective(a_rho: GridFunction, rho: GridFunction, Kinv : BaseMatrix)
                     x0 : GridFunction, 
                     descent : callable = lambda x: -x, # takes d_objective GridFunction as input, returns descent GridFunction
                     step : float = 1,
                     x_min : float =  0, 
                     x_max : float =  1,
                     maxit : int = 100,
                     tol : float = 1e-6,
                     step_min : float = 1e-6,
                     step_max : float = np.inf,
                     fac_increase_step : float= 1.2,
                     fac_decrease_step : float= 0.5,
                     draw : bool = True,
                     verbose : int = 2):
    
    
    # Initialization
    x_list = [x0]
    x_accepted = x0
    F_list = []
    dF_list = []
    crit_list = [1]
    characteristic_function = GridFunction(x0.space)
    characteristic_function.Set(1)
    surface = Integrate(characteristic_function, x0.space.mesh)

    if draw:
        scene = Draw(x_list[-1], x0.space.mesh,
                         min = x_min, max = x_max,
                         settings = {"Objects" : {"Wireframe" : False}, 
                                     "Colormap" : {"ncolors" : 32}})
    
    iter = 0

    while 1:
        iter += 1

        # 1) Compute physical state
        sol, Kinv = state(x_list[-1])
    
        # 2) Compute objective function
        F_list.append(objective(x_list[-1], sol))

        # 3) Update step & compute stop criterion
        if len(F_list) == 1:            # initial descent
            dF_list.append(d_objective(sol, x_list[-1], Kinv))
            x_accepted, F_accepted, dF_accepted = x_list[-1], F_list[-1], dF_list[-1]
            dx_crit = GridFunction(x0.space)
            dx_crit.Set(x_list[-1] - step_min * dF_list[-1])
            crit_ref = Integrate(Norm(project(dx_crit, x_min, x_max) - x_list[-1]), x0.space.mesh) / (step_min*surface)
            if verbose >= 2: print(f"it {iter} ✅| obj = {F_list[-1] :.3e} | {step = :.2e} | crit = {crit_list[-1] :.2e}")
            
        elif F_list[-1] < F_accepted:   # accept the step
            dF_list.append(d_objective(sol, x_list[-1], Kinv))
            x_accepted, F_accepted, dF_accepted = x_list[-1], F_list[-1], dF_list[-1]
            dx_crit.Set(x_list[-1] - step_min * dF_list[-1])
            crit_list.append(Integrate(Norm(project(dx_crit, x_min, x_max) - x_list[-1]), x0.space.mesh) / (step_min*surface) / crit_ref)
            step = min(step_max, fac_increase_step*step)
            if verbose >= 2: print(f"it {iter} ✅| obj = {F_list[-1] :.3e} | {step = :.2e} | crit = {crit_list[-1] :.2e}")

        else:                           # reject the step
            step *= fac_decrease_step
            crit_list.append(crit_list[-1])
            if verbose >= 2:  print(f"it {iter} ❌| obj = {F_list[-1] :.3e} | {step = :.2e} | crit = {crit_list[-1] :.2e}")

        # 4) Check convergence
        if crit_list[-1] < tol: # stop if converged
            print(f"Converged! Relative norm of projected gradient = {crit_list[-1] :.2e} < {tol :.2e}")
            break

        if iter >= maxit:
            print("Maximum number of iterations reached!")
            break

        if step < step_min:
            print("Step lower than minimum step size!")
            break

        # 5) Update variable
        x_test = copy(x_accepted)
        x_test.vec.data += step * descent(dF_accepted).vec
        x_list.append(project(x_test, x_min, x_max))

        if draw:
            scene.Redraw(x_list[-1], x0.space.mesh,
                         min = x_min, max = x_max,
                         settings = {"Objects" : {"Wireframe" : False}, 
                                     "Colormap" : {"ncolors" : 32}})
            
    results = {"solution" :  x_list,
               "objective" : F_list,
               "gradient" :  dF_list,
               "criterion" : crit_list,
               "state" : sol}
        
    return results


from utils.solver import DrawMaterial

def modulo(x : GridFunction, 
            mod : float = 2*np.pi):
    x.vec.data = np.mod(x.vec, mod)
    return x

def gradient_descent2(state : callable,          # a_rho_theta : GridFunction, Kinv : BaseMatrix = state(rho : GridFunction, theta : GridFunction)
                     objective : callable,       # obj : float = objective(rho : GridFunction, theta: GridFunction, a_rho_theta : GridFunction)
                     drho_objective : callable,     # riesz_repr : GridFunction = d_objective(a_rhotheta : GridFunction, rho: GridFunction, theta: GridFunction, Kinv : BaseMatrix)
                     dtheta_objective : callable,     # riesz_repr : GridFunction = d_objective(a_rhotheta : GridFunction, rho: GridFunction, theta: GridFunction, Kinv : BaseMatrix)
                     rho0 : GridFunction, 
                     theta0 : GridFunction, 
                     descent_rho : callable = lambda x: -x, # takes d_objective GridFunction as input, returns descent GridFunction
                     descent_theta : callable = lambda x: -x, # takes d_objective GridFunction as input, returns descent GridFunction
                     step : float = 1,
                     rho_min : float =  0, 
                     rho_max : float =  1,
                     maxit : int = 100,
                     step_min : float = 1e-6,
                     step_max : float = np.inf,
                     fac_increase_step : float= 1.2,
                     fac_decrease_step : float= 0.5,
                     draw : bool = True,
                     verbose : int = 2):
    
    
    # Initialization
    rho_list, theta_list = [rho0], [theta0]
    rho_accepted, theta_accepted = rho0, theta0
    F_list = []
    dF_drho_list, dF_dtheta_list = [], []

    if draw:
        scene = DrawMaterial(rho_list[-1], theta_list[-1])
    
    iter = 0

    while 1:
        iter += 1

        # 1) Compute physical state
        sol, Kinv = state(rho_list[-1], theta_list[-1])
    
        # 2) Compute objective function
        F_list.append(objective(rho_list[-1], theta_list[-1], sol))

        # 3) Update step & compute stop criterion
        if len(F_list) == 1:            # initial descent
            dF_drho_list.append(drho_objective(sol, rho_list[-1], theta_list[-1], Kinv))
            dF_dtheta_list.append(dtheta_objective(sol, rho_list[-1], theta_list[-1], Kinv))
            rho_accepted, theta_accepted, F_accepted = rho_list[-1], theta_list[-1], F_list[-1]
            drhoF_accepted, dthetaF_accepted= dF_drho_list[-1], dF_dtheta_list[-1]
            if verbose >= 2: print(f"it {iter} ✅| obj = {F_list[-1] :.3e} | {step = :.2e}")
            
        elif F_list[-1] < F_accepted:   # accept the step
            dF_drho_list.append(drho_objective(sol, rho_list[-1], theta_list[-1], Kinv))
            dF_dtheta_list.append(dtheta_objective(sol, rho_list[-1], theta_list[-1], Kinv))
            rho_accepted, theta_accepted, F_accepted = rho_list[-1], theta_list[-1], F_list[-1]
            drhoF_accepted, dthetaF_accepted= dF_drho_list[-1], dF_dtheta_list[-1]
            step = min(step_max, fac_increase_step*step)
            if verbose >= 2: print(f"it {iter} ✅| obj = {F_list[-1] :.3e} | {step = :.2e}")

        else:                           # reject the step
            step *= fac_decrease_step
            if verbose >= 2:  print(f"it {iter} ❌| obj = {F_list[-1] :.3e} | {step = :.2e}")

        # 4) Check convergence
        if iter >= maxit:
            print("Maximum number of iterations reached!")
            break

        if step < step_min:
            print("Step lower than minimum step size!")
            break

        # 5) Update variable
        rho_test = copy(rho_accepted)
        rho_test.vec.data  += step * descent_rho(drhoF_accepted).vec
        rho_list.append(project(rho_test, rho_min, rho_max))

        theta_test = copy(theta_accepted)
        theta_test.vec.data  += step * descent_theta(dthetaF_accepted).vec

        theta_list.append(modulo(theta_test, 2*np.pi))

        if draw:
            DrawMaterial(rho_list[-1], theta_list[-1], scene)
            
            
    results = {"rho" :  rho_list,
               "theta" : theta_list,
               "objective" : F_list,
               "state" : sol}
        
    return results



