"""
This code is to generate solution for a 2D poisson equation with stochastic lognormal diffusivity on an L-shape domain.

"""

#from dolfin import *
from fenics import *
from scipy.sparse import *
import scipy.io as sciio
import scipy.io
import numpy as np
import sys
from numpy import intc
from scipy.io import loadmat
from scipy.io import savemat

#parameters["linear_algebra_backend"] = "uBLAS"
parameters["linear_algebra_backend"] = "Eigen"

# Read mesh
mesh = Mesh('./mesh/mesh.xml')
#plot(mesh, interactive=True)


# Create mesh function over cell facets (for boundary subdomains) and define boundaries
my_eps = .0001
boundary_parts = MeshFunction("size_t", mesh, mesh.topology().dim()-1)

# Mark upper boundary facets as subdomain 0
class UpperDirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and abs(x[1]-1.0) < my_eps

Gamma_D_U = UpperDirichletBoundary()
Gamma_D_U.mark(boundary_parts, 0)


# Mark right boundary facets as subdomain 1
class RightDirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and abs(x[0]-1.0) < my_eps

Gamma_D_R = RightDirichletBoundary()
Gamma_D_R.mark(boundary_parts, 1)


# Mark left boundary facets as subdomain 2
class LeftNeumannBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and abs(x[0]+1.0) < .01 #my_eps

Gamma_N_L = LeftNeumannBoundary()
Gamma_N_L.mark(boundary_parts, 2)

# Define the boundary integral variables
#ds = ds[boundary_parts]
ds = Measure('ds',domain=mesh, subdomain_data=boundary_parts)

# Load Young's modulus
content = loadmat('./fenics_inputs/Youngs.mat')
young = content['Youngs']


content = loadmat('./fenics_inputs/inputs.mat')
nu = float(content['nu'])
mode_qoi = float(content['mode_qoi'])
theta = float(content['theta'])
delta_q = float(content['delta_q'])

# Define the traction force on the left boundary (Neumann)
# use theta in trapezoid. This will possibly work...
# Currently does not work. Fix this. Shouldn't be hard.
#load_profile = ('1.0+ x[1]*tan(theta)', '0.0')
#load_profile = (1.0,0.0)
# Test some other stuff
# f = Constant(Expression(load_profile))
# maybe this does something?
#f1 = Constant((1.0,0.0))


f = Expression(('1.0*(1 + delta_q) + x[1]*tan(theta)', '0.0'),degree=1, delta_q=delta_q, theta=theta)



# Define function space and basis functions
V = VectorFunctionSpace(mesh, "CG", 1)
W = FunctionSpace(mesh, "CG", 1)
u = TrialFunction(V)
v = TestFunction(V)

E = Function(W)

V_0, V_1 = V.split()

# Define Dirichlet boundaries
bcs = [DirichletBC(V_1, Constant(0.0), boundary_parts, 0),
       DirichletBC(V_0, Constant(0.0), boundary_parts, 1)]




E = Function(W)

# Define stress
def sigma(u):
    return 2.0*mu*sym(grad(u)) + lmbda*tr(sym(grad(u)))*Identity(u.geometric_dimension())


# Setup dof map
# can I change this from V to W?? V is  a vector function space, W a function space...
#dof_map = dof_to_vertex_map(V)
dof_map = dof_to_vertex_map(W)

dof_coords = W.tabulate_dof_coordinates().reshape((W.dim(),-1))
#scipy.io.savemat('dof_coords.mat', mdict={'dof_coords':dof_coords})


# Create files for storing solution
#ufile = File("solution.pvd")
#umfile = open("solution_middle.txt","w")

d = mesh.geometry().dim()
#Vsig = FunctionSpace(mesh, 'CG', 1)
# Vsig = FunctionSpace(mesh, 'P', 1)


for i in range(len(young)):
    #print(i)
    u = TrialFunction(V)
    v = TestFunction(V)
    E = Function(W)
    for j in range(len(E.vector())):
        #print(j)
       	E.vector()[j] = np.array(np.mat(young[i,:]))[0,dof_map[j]]
        #print E.vector()[j]
        #IndexError: index 146 is out of bounds for axis 1 with size 79.
        # 910 is out of bounds for size 465
        # E has size 465, dof_map[j] is 910...
    #print 'done with reading E'
    #E_proj = project(E,W)
    #plot(E_proj,interactive=True)
    #print 'done with plotting E'
    # set nu earlier as variable
    #nu = 0.3
    #nu = 0.15
    mu = E / (2.0*(1.0 + nu))
    lmbda = E*nu / ((1.0 + nu)*(1.0 - 2.0*nu))
    a = inner(sigma(u), grad(v))*dx
    L = inner(f, v)*ds(2)
    u = Function(V)
    solve(a == L, u, bcs)
    #print 'done with solving'
    # Extract components of solution
    u_0, u_1 = u.split(True)
    #plot(u_0,interactive=True)
    np.set_printoptions(threshold='nan')
    fname = 'L_data/solution.' + str(i) + '.txt'
    #Put u in array form
    #u_0_nodal_values = u_0.vector()
    #u_0_array = u_0_nodal_values.array()
    #File("u_pvd.pvd") << u_0

    # hardcode coordinates? Nah. Use dof_map / find coordinates.
    if mode_qoi == 0:
        indices = np.where(np.logical_and(abs(dof_coords[:, 0]) < 0.000001, dof_coords[:, 1] > -0.000001))[0]
        # indices 4 and 5 (and likely more are in wrong order.)
        # dof_y = dof_coords[:,1]
        # dof_y[indices]
        u_0_array = u_0.vector()[indices]
    elif mode_qoi == 1:
        u_0_array=u_0.vector()
    elif mode_qoi == 2:
        #https://fenicsproject.org/pub/tutorial/html/._ftut1008.html
        s = sigma(u) - (1./3)*tr(sigma(u))*Identity(d)
        von_Mises_1 = sqrt(3./2*inner(s, s))
        von_Mises = project(von_Mises_1, W)
        u_0_array=von_Mises.vector()

    np.savetxt(fname,u_0_array)

    # to save pvd...
    # if i == 0:
    #     # File("u_0_coarse.pvd") << u_0
    #     # File("u_1_coarse.pvd") << u_1
    #     # File("sig_coarse.pvd") << von_Mises
    #     File("u_0_fine.pvd") << u_0
    #     File("u_1_fine.pvd") << u_1
    #     File("sig_fine.pvd") << von_Mises
