import sys
sys.path.append('..')
from tricubic_interpolation import Tricubic_Interpolation
import matplotlib.pyplot as plt
import numpy as np
import sympy

x,y,z = sympy.symbols('x y z')
f = 0.0001*sympy.exp(-z)*(0.03*x**4-1.5*x**3)*y**3+10*x*y*z
lamf = sympy.lambdify((x,y,z), f, modules='numpy')

x0 = 0.
y0 = 0.
z0 = 0.
dx = 0.5
dy = 0.5
dz = 0.3
discard_x = 3
discard_y = 5
discard_z = 1

### create array
Nx = 100
Ny = 100
Nz = 100

A = np.empty([Nx,Ny,Nz])
for i in range(Nx):
    xi=dx*i + x0
    for j in range(Ny):
        yi=dy*j + y0
        for k in range(Nz):
            zi=dz*k + z0
            A[i,j,k] = lamf(xi, yi, zi) 

#default values for x,y,z when they are not variable
x_obs = 3.
y_obs = 6.
z_obs = 2.

ip = Tricubic_Interpolation(A, x0, y0, z0, dx, dy, dz, discard_x, discard_y, discard_z)

print('f     = %f'%ip.val(x_obs,y_obs,z_obs))
print('df/dx = %f'%ip.ddx(x_obs,y_obs,z_obs))
print('df/dy = %f'%ip.ddy(x_obs,y_obs,z_obs))
print('df/dz = %f'%ip.ddz(x_obs,y_obs,z_obs))
