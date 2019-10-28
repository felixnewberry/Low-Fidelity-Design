from fenics import *
from mshr import *
import numpy as np
#import dolfin as dfn
import scipy.io as sciio
import scipy.io
import scipy.sparse as sparse
from scipy.io import savemat
from scipy.io import loadmat
#import os
import time

#set_log_level(WARNING)
#set_log_active(False)

################################################################################
### LDC Approach

#Record simulation time
start = time.time()

from NavierStokesLDC import Navier_Stokes_LDC

# Input parameters
nu = Constant(0.01)
u_lid = Constant(1.0)
nx_vec = np.array((4,8,16,32,64,128,256))
#nx_vec = np.array((128,256))

# old system. Now obtain all QoI at once.
QoI = float(1)
# 0 is u_mat - u and v
# 1 is P_mat
# 2 is u_vec mid - just u...
# 3 is p_vec mid
# 4 is p_vec top
# 5 is p_vec vert
# 6 is p_vec base

# % also check vertical velocity u and v seperately.
# I should do one run with all these qoi output? Maybe...
# look at convergence first,
# also consider increasing range of random variables.
# interpolate solution to 256 ..
# For now interpolate to 32.


u_matrix = np.zeros((len(nx_vec),65))
# 65, 4225, 8450
# else:
#     u_matrix = np.zeros((200,33))


for i in range(len(nx_vec)): #200
    sample_i = i;

    #sample_i = int(Nan_vec[0,i])
    nx = int(nx_vec[i])
    print(nx)
    u_array_samp = Navier_Stokes_LDC(u_lid, nu, sample_i, nx, QoI)

    u_matrix[i,:] = u_array_samp

# Save QoI
scipy.io.savemat('./u_meshes/u_mesh_u_vert.mat', mdict={'u_matrix':u_matrix})

end = time.time()
print(end-start)
