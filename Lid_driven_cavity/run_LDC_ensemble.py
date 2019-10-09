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

# Load Random Variables
content = loadmat('./LDC_data/nu_vec.mat')
nu_vec = content['nu_vec']
#mu_vec = float(mu_vec)
content = loadmat('./LDC_data/u_top_vec.mat')
u_lid_vec = content['u_top_vec']

# content = loadmat('./sample_i.mat')
# sample_i = int(content['sample_i'])

content = loadmat('./LDC_data/inputs_vec.mat')
nx = int(content['nx'])
delta_u = float(content['delta_u'])
delta_nu = float(content['delta_nu'])
QoI = float(content['QoI'])

# # 32
# u_x_matrix = np.zeros((10,1089))
# u_y_matrix = np.zeros((10,1089))
# u_mag_matrix = np.zeros((10,1089))

# 32
# u_field_matrix = np.zeros((10,2178))

# 8
# u_x_matrix = np.zeros((10,81))
# u_y_matrix = np.zeros((10,81))
# u_mag_matrix = np.zeros((10,81))

# p_field_matrix = np.zeros((200,1089))
# p_line_matrix = np.zeros((200,33))

# # 8
# u_field_matrix = np.zeros((200,578))
# p_field_matrix = np.zeros((200,81))

# # 6
# u_field_matrix = np.zeros((200,338))
# p_field_matrix = np.zeros((200,49))

# # 4
# u_field_matrix = np.zeros((200,162))
# p_field_matrix =

# u_matrix = np.zeros((200,66))

if QoI == 0:
    u_matrix = np.zeros((200,2178))
elif QoI == 1:
    u_matrix = np.zeros((200,1089))
elif QoI == 2:
    u_matrix = np.zeros((200,66))
elif QoI == 3:
    u_matrix = np.zeros((200,33))
else:
    u_matrix = np.zeros((200,33))

# Apply deltas to u and nu % want them to be multiplicative. :/
nu_vec = nu_vec*(1+delta_nu)
u_lid_vec = u_lid_vec*(1+delta_u)

print(QoI)

for i in range(200): #200
    sample_i = i;
    #sample_i = int(Nan_vec[0,i])
    u_lid = float(u_lid_vec[sample_i])
    u_lid = Constant(u_lid)

    nu = float(nu_vec[sample_i])
    nu = Constant(nu)

    u_array_samp = Navier_Stokes_LDC(u_lid, nu, sample_i, nx, QoI)

    u_matrix[i,:] = u_array_samp

# Save QoI
scipy.io.savemat('./u_meshes/u_matrix.mat', mdict={'u_matrix':u_matrix})

end = time.time()
print(end-start)
