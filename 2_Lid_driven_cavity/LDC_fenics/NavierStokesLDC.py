"""
This demo program solves the steady incompressible Navier-Stokes equations
for lid-driven cavity problem using Taylor-Hood elements.
Author: Praveen. C
www   : http://math.tifrbng.res.in/~praveen

Adjustments made by Felix
"""

# python3.6 LDC_nonlinear.py

################################################################################
### Import
################################################################################
#from __future__ import division
#from __future__ import print_function
from fenics import *
#from dolfin import *
#from scipy.sparse import *
import scipy.io as sciio
#import scipy.io
import numpy as np
import sys
import scipy.io
from scipy import sparse
#import petsc4py
from numpy import intc
import time
#from scipy.io import savemat
#from scipy.io import loadmat

set_log_level(LogLevel.WARNING)
#Record simulation time
# start = time.time()

# load mesh size
#content = sciio.loadmat('./nx.mat')
#nx = int(content['nx'])

def Navier_Stokes_LDC(u_top, nu, nx):
    # Inputs:
    # u_top - velocity of top plate
    # nu - viscosity of fluid
    # sample_i - sample number
    # nx - number of cells along edge, ie 4 or 32

    # Outputs: x and y velocities, and magnitude
    # In previous versions the output has been pressure along the top surface.

    ny = nx

    n_points = nx+1

    # Define Mesh
    mesh = UnitSquareMesh(nx, ny)
    x = mesh.coordinates()
    # skew mesh towards walls (as in Oasis)
    x[:] = (x - 0.5) * 2
    x[:] = 0.5*(np.cos(np.pi*(x-1) / 2) + 1)

    # nu varies linearly from nu_0 to nu_1
    #nu = Expression('nu_0 + (nu_1-nu_0)*x[1]', nu_1= nu_1, nu_0=nu_0, degree=1)

    # nu sigmoid depression centered at vortex.
    # vortex center is 0.6215, 0.7357
    # sigmoid that is 1 at origin and 0 at about -0.75 is defined as:
    # nu_new = (1+s)*nu_old where s = 1/(1+exp(-(15x+5))) where x is distance from vortex.

    # sig_tight = 30;
    # sig_center = 0.25;
    # vort_x = 0.6684;
    # vort_y = 0.7357;
    #
    # nu = Expression('nu_0+(nu_1-nu_0)/(1+exp(sig_tight*((sqrt(pow(x[0]-vort_x,2)+pow(x[1]-vort_y,2)))-sig_center)))', nu_0=nu_0,nu_1=nu_1,sig_tight=sig_tight,vort_x=vort_x, vort_y=vort_y, sig_center=sig_center, degree=2)

    V1 = VectorElement("Lagrange",mesh.ufl_cell(),2)
    P1 = FiniteElement("Lagrange",mesh.ufl_cell(),1)

    # This gives the right size. But order is wacky.
    W = FunctionSpace(mesh,V1*P1)

    # load pertubation
    # content = sciio.loadmat('./delta_v_vec.mat')
    # delta_v = content['delta_v']
    # content = sciio.loadmat('./delta_nu_vec.mat')
    # delta_nu = content['delta_nu']
    #
    # delta_v = float(delta_v)
    # delta_nu = float(delta_nu)

    # Define test functions
    (v,q) = TestFunctions(W)

    # Define trial functions
    w     = Function(W)
    # This change seemed to break things...
    #(u,p) = TrialFunctions(W)
    (u,p) = (as_vector((w[0], w[1])), w[2])

    # Set parameter values
    Re = 100
    #nu = 1.0/Re
    nu = Constant(nu)


    # Define boundary conditions

    #u_top = 1    # Velocity of top plate.
    u_top = Constant(u_top)
    noslip  = DirichletBC(W.sub(0), (0, 0), "x[0] < DOLFIN_EPS || x[0] > 1.0 - DOLFIN_EPS || x[1] < DOLFIN_EPS")
    lid  = DirichletBC(W.sub(0), (u_top,0), "x[1] > 1.0 - DOLFIN_EPS")
    pref = DirichletBC(W.sub(1), 0, "x[0] < DOLFIN_EPS && x[1] < DOLFIN_EPS", "pointwise")

    bcs = [noslip, lid, pref]


    for bc in bcs:
        bc.apply(w.vector())

    # Tentative velocity step
    F =   inner(grad(u)*u, v)*dx \
        + nu*inner(grad(u), grad(v))*dx \
        - div(v)*p*dx \
        - q*div(u)*dx

    dw = TrialFunction(W)
    dF = derivative(F, w, dw)

    #if True:
    nsproblem = NonlinearVariationalProblem(F, w, bcs, dF)
    solver = NonlinearVariationalSolver(nsproblem)

    prm = solver.parameters
    prm['newton_solver']['absolute_tolerance'] = 1E-6
    prm['newton_solver']['relative_tolerance'] = 1E-6
    prm['newton_solver']['maximum_iterations'] = 1000
    prm['newton_solver']['relaxation_parameter'] = 1.0

    # umfpack (Unsymmetric MultiFrontal sparse LU factorization)
    #petsc PETSc built in LU solver
    #cg Conjugate Gradient method - no good
    # try mumps..
    # prm['newton_solver']['linear_solver'] = "umfpack"
    prm['newton_solver']['linear_solver'] = "mumps"

    a = solver.solve()


    (u,p) = w.split(True)

    #####################################################################
    ### QoI computation
    #####################################################################
    # 0 u_mat, 1 P_mat 2, u_vec mid, 3, p_vec mid, 4, p_vec top
    #
    # y_8 = x[x[:,:][:,0]==0.5]
    # x_8 = x[x[:,:][:,1]==0.5]

    # y_64 = x[x[:,:][:,0]==0.0]
    # x_64 = x[x[:,:][:,1]==0.0]

    # If I do interpolate:
    # Load 32 coordinates and interpolate x and y velocites mid slices of domain.

    # timing test
    #####################################################################
    # Load 64 coordinates and interpolate u and p slices

    content = sciio.loadmat('./LDC_data/x_64.mat')
    x_64 = content['x_64']


    u_y_array = np.zeros(x_64.shape[0])
    # u_x_array = np.zeros(x_64.shape[0])
    for i_coords in range(x_64.shape[0]):
        u_y_array[i_coords] = u(x_64[i_coords][0],0.5)[1]

    u_x_array = np.zeros(x_64.shape[0])
    # u_x_array = np.zeros(x_64.shape[0])
    for i_coords in range(x_64.shape[0]):
        u_x_array[i_coords] = u(0.5,x_64[i_coords][0])[0]

    p_array_mid = np.zeros(x_64.shape[0])
    for i_coords in range(x_64.shape[0]):
        p_array_mid[i_coords] = p(x_64[i_coords][0],0.5)


    p_array_vert = np.zeros(x_64.shape[0])
    for i_coords in range(x_64.shape[0]):
        p_array_vert[i_coords] = p(0.5,x_64[i_coords][0])

    p_array_base = np.zeros(x_64.shape[0])
    for i_coords in range(x_64.shape[0]):
        p_array_base[i_coords] = p(x_64[i_coords][0],0)


    ######################
    p_field = np.zeros(x_64.shape[0]**2)
    u_x_field = np.zeros(x_64.shape[0]**2)
    u_y_field = np.zeros(x_64.shape[0]**2)

    for i_x in range(x_64.shape[0]):
        for i_y in range(x_64.shape[0]):
            p_field[i_x*x_64.shape[0]+i_y] = p(x_64[i_x][0],x_64[i_y][0])
            u_x_field[i_x*x_64.shape[0]+i_y] = u(x_64[i_x][0],x_64[i_y][0])[0]
            u_y_field[i_x*x_64.shape[0]+i_y] = u(x_64[i_x][0],x_64[i_y][0])[1]

    ######################

    return u_y_array, u_x_array, p_array_mid, p_array_vert, p_array_base, p_field, u_x_field, u_y_field

    ## for Saby
    #return u_coords_array, u_x_array, u_y_array, u_mag #, p_array
