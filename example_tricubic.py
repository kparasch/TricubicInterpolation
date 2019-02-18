from tricubic_interpolation import Tricubic_Interpolation
import matplotlib.pyplot as plt
import numpy as np
import sympy
#plt.style.use('kostas')

x,y,z = sympy.symbols('x y z')
f = 0.0001*sympy.exp(-z)*(0.03*x**4-1.5*x**3)*y**3
dfdx=sympy.diff(f,x)
dfdy=sympy.diff(f,y)
dfdz=sympy.diff(f,z)
lamf = sympy.lambdify((x,y,z), f, modules='numpy')
lamdfdx = sympy.lambdify((x,y,z), dfdx, modules='numpy')
lamdfdy = sympy.lambdify((x,y,z), dfdy, modules='numpy')
lamdfdz = sympy.lambdify((x,y,z), dfdz, modules='numpy')

x0 = 0.
y0 = 0.
z0 = 0.
dx = 0.513
dy = 0.1
dz = 0.36
discard_x = 1
discard_y = 1
discard_z = 1
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

x_obs = 3.
y_obs = 2.
z_obs = 2.

ip = Tricubic_Interpolation()
ip.discard_points(discard_x, discard_y, discard_z)
ip.set_steps(dx, dy, dz)
ip.set_origin(x0, y0, z0)
ip.initialize(A)

ix = int(x0)
iy = int(y0)
iz = int(z0)


X = np.linspace(x0+dx,x0+(Nx-3)*dx,3000)
Y = np.linspace(y0+dy,y0+(Ny-3)*dy,3000)
Z = np.linspace(z0+dz,z0+(Nz-3)*dz,3000)

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
ax11.plot(Z, np.array([ip.ddy(x_obs,z_obs, i) for i in Z]),'r')
ax11.set_xlabel('z')

ax12=fig.add_subplot(3,4,12)
ax12.plot(Z, lamdfdz(x_obs,y_obs, Z),'k')
ax12.plot(Z, np.array([ip.ddz(x_obs,y_obs, i) for i in Z]),'r')
ax12.set_xlabel('z')

plt.show(False)
input()
