"""
This code is to generate solution for a 2D cantilever beam with composite cross section.

"""
# currently exiting on segmentation fault...
# whwere does this happen? Run parts of the code to find out


#from dolfin import *
from fenics import *
from mshr import *
from scipy.sparse import *
import scipy.io as sciio
import scipy.io
import numpy as np
import sys
from numpy import intc
#parameters["linear_algebra_backend"] = "uBLAS"


#if not has_cgal():
#        print "DOLFIN must be compiled with CGAL to run this demo."
#        exit(0)

r = Rectangle(Point(0.0, 0.0), Point(50.0, 5.4))#, diagonal="right")
c1 = Circle(Point(5,2.7),1.5,64)
c2 = Circle(Point(15,2.7),1.5,64)
c3 = Circle(Point(25,2.7),1.5,64)
c4 = Circle(Point(35,2.7),1.5,64)
c5 = Circle(Point(45,2.7),1.5,64)

g2d = r - c1 - c2 - c3 - c4 - c5
mesh = generate_mesh(g2d, 70)

for i in range(2):
    #cell_markers = CellFunction("bool", mesh,0)
    cell_markers = MeshFunction("bool", mesh,0)
    cell_markers.set_all(False)
    for cell in cells(mesh):
        p = cell.midpoint()
        if (p[1] < 0.4) or (p[1] > 5.0):
            cell_markers[cell] = True
    mesh = refine(mesh, cell_markers)

for i in range(1):
    cell_markers = CellFunction("bool", mesh)
    cell_markers.set_all(False)
    for cell in cells(mesh):
        p = cell.midpoint()
        if p.distance(Point(5,2.7)) < 2.5 or \
           p.distance(Point(15,2.7)) < 2.5 or \
           p.distance(Point(25,2.7)) < 2.5 or \
           p.distance(Point(35,2.7)) < 2.5 or \
           p.distance(Point(45,2.7)) < 2.5:
           cell_markers[cell] = True
    mesh = refine(mesh, cell_markers)


#print("mesh.hmax() = {}  mesh.hmin() = {}".format(mesh.hmax(), mesh.hmin()))

# Create mesh function over cell facets (for boundary subdomains) and define boundaries
my_eps = 1E-5
#boundary_parts = MeshFunction("size_t", mesh, mesh.topology().dim()-1)

print("Made it to 62")

# Mark upper boundary facets as subdomain 0
class UpperTractionBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and abs(x[1]-5.4) < my_eps

def ClampedBoundary(x, on_boundary):
    return on_boundary and x[0] < my_eps


# Define the boundary integral variables
#ds = ds(subdomain_data=boundary_parts)

# Define the traction force on the left boundary (Neumann)
f = Constant((0.0,q_samp))


# Define function space and basis functions
V = VectorFunctionSpace(mesh, "CG", 1)
W = FunctionSpace(mesh, "CG", 1)
u = TrialFunction(V)
v = TestFunction(V)

V_0, V_1 = V.split()

# Define Dirichlet BC
bc = DirichletBC(V, Constant((0, 0)), ClampedBoundary)

# Define traction BC
uppertractionboundary = UpperTractionBoundary()
boundaries = FacetFunction('size_t', mesh, 0)
uppertractionboundary.mark(boundaries, 1)

# Define outer surface measure aware of Dirichlet and Neumann boundaries
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)


# Define stress
def sigma(u):
    return 2.0*mu*sym(grad(u)) + lmbda*tr(sym(grad(u)))*Identity(u.geometric_dimension())


# Create files for storing solution
utfile = open("solution_top.txt","w")

# Define appropriate function spaces
u = TrialFunction(V)
v = TestFunction(V)

# Define Young's moudlus

class E(Expression):
    def set_E_values(self, E1, E2, E3):
        self.E_1 = E1
        self.E_2 = E2
        self.E_3 = E3
    def eval(self, value, x):
        "Set value[0] to value at point x"
        if x[1] > 5.2 - my_eps:
            value[0] = self.E_1
        if (x[1] < .2 + my_eps):
            value[0] = self.E_2
        if (x[1] >= .2 + my_eps) and (x[1] <= 5.2 - my_eps):
            value[0] = self.E_3

# Initialize kappa
young = E(degree=0)
young.set_E_values(E1_samp,E2_samp,E3_samp)

# Define other properties
nu = 0.3
mu = young / (2.0*(1.0 + nu))
lmbda = young*nu / ((1.0 + nu)*(1.0 - 2.0*nu))

# Define the variational form
a = inner(sigma(u), grad(v))*dx
L = inner(f, v)*ds(1)
u = Function(V)
solve(a == L, u, bc)
#print 'done with solving'

# Extract components of solution
u_0, u_1 = u.split(True)

#Put u in array form
u_1_nodal_values = u_1.vector()
u_1_array = u_1_nodal_values.array()
#np.savetxt(fname,u_1_array)


#vtkfile = File('solution.pvd')
#vtkfile << u

coordinates = V.tabulate_dof_coordinates()

for i in range(len(u_1_array)):
    #print (i, coordinates[4*i], coordinates[4*i+1], u_1_array[i])
    if coordinates[4*i+1]> 5.4 - my_eps:   # fixed line y = 5.2
        utfile.write(str(u_1_array[i])+ '\n') #(coordinates[4*i], coordinates[4*i+1], u_1_array[i])
