import sys
sys.path.append('../..')
#from TricubicInterpolation.pyTricubic import Tricubic_Interpolation
from TricubicInterpolation.cTricubic import Tricubic_Interpolation
import matplotlib.pyplot as plt
import numpy as np
import sympy
#plt.style.use('kostas')

plt.close('all')

x,y,z = sympy.symbols('x y z')
f = 0.0001*sympy.exp(-z)*(0.03*x**4-1.5*x**3)*y**3+10*x*y*z
f = z**2*(0.03*x**1-1.5*x**3)*y**3*sympy.sin(x*y/10)
#f = z**2*(x**3)*y**3

# symbolic derivatives
dfdx = sympy.diff(f,x)
dfdy = sympy.diff(f,y)
dfdz = sympy.diff(f,z)
dfdxdy = sympy.diff(dfdx,y)
dfdxdz = sympy.diff(dfdx,z)
dfdydz = sympy.diff(dfdy,z)
dfdxdydz = sympy.diff(dfdxdy,z)

# convert symbolic functions to numerical functions
lamf = sympy.lambdify((x,y,z), f, modules='numpy')
lamdfdx = sympy.lambdify((x,y,z), dfdx, modules='numpy')
lamdfdy = sympy.lambdify((x,y,z), dfdy, modules='numpy')
lamdfdz = sympy.lambdify((x,y,z), dfdz, modules='numpy')
lamdfdxdy = sympy.lambdify((x,y,z), dfdxdy, modules='numpy')
lamdfdxdz = sympy.lambdify((x,y,z), dfdxdz, modules='numpy')
lamdfdydz = sympy.lambdify((x,y,z), dfdydz, modules='numpy')
lamdfdxdydz = sympy.lambdify((x,y,z), dfdxdydz, modules='numpy')

exact_derivatives = True 

x0 = 0.
y0 = 0.
z0 = 0.
dx = 0.5/2
dy = 0.5/2
dz = 0.3/2
discard_x = 1
discard_y = 5
discard_z = 1
Nx = 100
Ny = 40
Nz = 30


if exact_derivatives:
    A = np.empty([Nx,Ny,Nz,8])
    for i in range(Nx):
        xi = dx*i + x0
        for j in range(Ny):
            yi = dy*j + y0
            for k in range(Nz):
                zi = dz*k + z0
                A[i,j,k,0] = lamf(xi, yi, zi) 
                A[i,j,k,1] = lamdfdx(xi, yi, zi)
                A[i,j,k,2] = lamdfdy(xi, yi, zi)
                A[i,j,k,3] = lamdfdz(xi, yi, zi)
                A[i,j,k,4] = lamdfdxdy(xi, yi, zi)
                A[i,j,k,5] = lamdfdxdz(xi, yi, zi)
                A[i,j,k,6] = lamdfdydz(xi, yi, zi)
                A[i,j,k,7] = lamdfdxdydz(xi, yi, zi)
else:
    A = np.empty([Nx,Ny,Nz])
    for i in range(Nx):
        xi = dx*i + x0
        for j in range(Ny):
            yi = dy*j + y0
            for k in range(Nz):
                zi = dz*k + z0
                A[i,j,k] = lamf(xi, yi, zi) 

#default values of x,y,z when they are not variable
x_obs = 3.
y_obs = 6.
z_obs = 2.

if exact_derivatives:
    ip = Tricubic_Interpolation(A, x0, y0, z0, dx, dy, dz, discard_x, discard_y, discard_z, method='Exact')
else:
    ip = Tricubic_Interpolation(A, x0, y0, z0, dx, dy, dz, discard_x, discard_y, discard_z, method='FD')

X = np.linspace(x0+(discard_x)*dx,x0+(Nx-2-discard_x)*dx,3000)
Y = np.linspace(y0+(discard_y)*dy,y0+(Ny-2-discard_y)*dy,3000)
Z = np.linspace(z0+(discard_z)*dz,z0+(Nz-2-discard_z)*dz,3000)

fig=plt.figure(1,[18,12])
fig.suptitle('red: tricubic interpolation, black: true')

ax1=fig.add_subplot(3,4,1)
ax1.plot(X, lamf(X,y_obs,z_obs),'k')
ax1.plot(X, np.array([ip.val(i,y_obs,z_obs) for i in X]),'r')
ax1.set_xlabel('x')
ax1.set_title('f')

ax2=fig.add_subplot(3,4,2)
ax2.plot(X, lamdfdx(X,y_obs,z_obs),'k')
ax2.plot(X, np.array([ip.ddx(i,y_obs,z_obs) for i in X]),'r')
ax2.set_xlabel('x')
ax2.set_title('df/dx')

ax3=fig.add_subplot(3,4,3)
ax3.plot(X, lamdfdy(X,y_obs,z_obs),'k')
ax3.plot(X, np.array([ip.ddy(i,y_obs,z_obs) for i in X]),'r')
ax3.set_xlabel('x')
ax3.set_title('df/dy')

ax4=fig.add_subplot(3,4,4)
ax4.plot(X, lamdfdz(X,y_obs,z_obs),'k')
ax4.plot(X, np.array([ip.ddz(i,y_obs,z_obs) for i in X]),'r')
ax4.set_xlabel('x')
ax4.set_title('df/dz')


ax5=fig.add_subplot(3,4,5)
ax5.plot(Y, lamf(x_obs,Y,z_obs),'k')
ax5.plot(Y, np.array([ip.val(x_obs,i,z_obs) for i in Y]),'r')
ax5.set_xlabel('y',labelpad=-4)

ax6=fig.add_subplot(3,4,6)
ax6.plot(Y, lamdfdx(x_obs, Y, z_obs),'k')
ax6.plot(Y, np.array([ip.ddx(x_obs,i,z_obs) for i in Y]),'r')
ax6.set_xlabel('y',labelpad=-4)

ax7=fig.add_subplot(3,4,7)
ax7.plot(Y, lamdfdy(x_obs,Y,z_obs),'k')
ax7.plot(Y, np.array([ip.ddy(x_obs,i,z_obs) for i in Y]),'r')
ax7.set_xlabel('y',labelpad=-4)

ax8=fig.add_subplot(3,4,8)
ax8.plot(Y, lamdfdz(x_obs,Y,z_obs),'k')
ax8.plot(Y, np.array([ip.ddz(x_obs,i,z_obs) for i in Y]),'r')
ax8.set_xlabel('y',labelpad=-4)


ax9=fig.add_subplot(3,4,9)
ax9.plot(Z, lamf(x_obs,y_obs, Z),'k')
ax9.plot(Z, np.array([ip.val(x_obs,y_obs, i) for i in Z]),'r')
ax9.set_xlabel('z')

ax10=fig.add_subplot(3,4,10)
ax10.plot(Z, lamdfdx(x_obs, y_obs, Z),'k')
ax10.plot(Z, np.array([ip.ddx(x_obs,y_obs, i) for i in Z]),'r')
ax10.set_xlabel('z')

ax11=fig.add_subplot(3,4,11)
ax11.plot(Z, lamdfdy(x_obs,y_obs, Z),'k')
ax11.plot(Z, np.array([ip.ddy(x_obs,y_obs, i) for i in Z]),'r')
ax11.set_xlabel('z')

ax12=fig.add_subplot(3,4,12)
ax12.plot(Z, lamdfdz(x_obs,y_obs, Z),'k')
ax12.plot(Z, np.array([ip.ddz(x_obs,y_obs, i) for i in Z]),'r')
ax12.set_xlabel('z')

plt.show()
