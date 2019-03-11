import sys
sys.path.append('..')

from tricubic_interpolation import Tricubic_Interpolation
from tricubic_interpolation2 import Tricubic_Interpolation2
import matplotlib.pyplot as plt
import numpy as np
import sympy
#plt.style.use('kostas')

def run_test():
    x,y,z = sympy.symbols('x y z')
    upper_int = 9
    #print('One random number: %d'%np.random.randint(upper_int))
    coefs_true = np.random.randint(upper_int, size=64)
    f = sum([ coefs_true[i + 4*j + 16*k] * x**i * y**j * z**k for i in range(4) for j in range(4) for k in range(4)])
    dfdx=sympy.diff(f,x)
    dfdy=sympy.diff(f,y)
    dfdz=sympy.diff(f,z)
    dfdxdy=sympy.diff(dfdx,y)
    dfdxdz=sympy.diff(dfdx,z)
    dfdydz=sympy.diff(dfdy,z)
    dfdxdydz=sympy.diff(dfdxdy,z)
    lamf = sympy.lambdify((x,y,z), f, modules='numpy')
    lamdfdx = sympy.lambdify((x,y,z), dfdx, modules='numpy')
    lamdfdy = sympy.lambdify((x,y,z), dfdy, modules='numpy')
    lamdfdz = sympy.lambdify((x,y,z), dfdz, modules='numpy')
    lamdfdxdy = sympy.lambdify((x,y,z), dfdxdy, modules='numpy')
    lamdfdxdz = sympy.lambdify((x,y,z), dfdxdz, modules='numpy')
    lamdfdydz = sympy.lambdify((x,y,z), dfdydz, modules='numpy')
    lamdfdxdydz = sympy.lambdify((x,y,z), dfdxdydz, modules='numpy')

    x0 = 0
    y0 = 0
    z0 = 0
    dx = 1#2
    dy = 1#3
    dz = 1#7
    discard_x = 1#1
    discard_y = 1#5
    discard_z = 1#1
    Nx = 30
    Ny = 30
    Nz = 30

    A = np.empty([Nx,Ny,Nz])
    B = np.empty([Nx,Ny,Nz,8])
    for i in range(Nx):
        xi=dx*i + x0
        for j in range(Ny):
            yi=dy*j + y0
            for k in range(Nz):
                zi=dz*k + z0
                A[i,j,k] = lamf(xi, yi, zi) 
                B[i,j,k,0] = lamf(xi, yi, zi) 
                B[i,j,k,1] = lamdfdx(xi, yi, zi)
                B[i,j,k,2] = lamdfdy(xi, yi, zi)
                B[i,j,k,3] = lamdfdz(xi, yi, zi)
                B[i,j,k,4] = lamdfdxdy(xi, yi, zi)
                B[i,j,k,5] = lamdfdxdz(xi, yi, zi)
                B[i,j,k,6] = lamdfdydz(xi, yi, zi)
                B[i,j,k,7] = lamdfdxdydz(xi, yi, zi)
    
    #default values of x,y,z when they are not variable
    x_obs = 3.
    y_obs = 6.
    z_obs = 2.
    
    ip = Tricubic_Interpolation(A, x0, y0, z0, dx, dy, dz, discard_x, discard_y, discard_z)
    ip2 = Tricubic_Interpolation(B, x0, y0, z0, dx, dy, dz, discard_x, discard_y, discard_z, 'Exact')
    
    passed = True
    n = 1
    for test in range(n):
        xv = np.random.rand()*(Nx-3-2*discard_x)*dx + x0 +discard_x*dx
        yv = np.random.rand()*(Ny-3-2*discard_y)*dy + y0 +discard_y*dy
        zv = np.random.rand()*(Nz-3-2*discard_z)*dz + z0 +discard_z*dz
        ix1,iy1,iz1 = ip.coords_to_indices(xv,yv,zv)
        ix2,iy2,iz2 = ip2.coords_to_indices(xv,yv,zv)
        coefs1 = ip.get_coefs(ip.construct_b(ix1, iy1, iz1))
        coefs2 = ip2.get_coefs(ip2.construct_b(ix2, iy2, iz2))
        fxyz = sympy.simplify(sum([ int(coefs2[i + 4*j + 16*k]) * ((x-ix2)/dx)**i * ((y-iy2)/dy)**j * ((z-iz2)/dz)**k for i in range(4) for j in range(4) for k in range(4)]))
#        passed = passed and np.array_equal(coefs_true, coefs2)
        print(fxyz)
        print(coefs_true)
        #print(coefs2)
        
        #fxyz = f
        coefs_test = np.zeros([64])
        for Terms_fyz in sympy.Poly(fxyz,x).all_terms():
            ii = Terms_fyz[0][0]
            fyz = Terms_fyz[1]
            for Terms_fz in sympy.Poly(fyz,y).all_terms():
                jj = Terms_fz[0][0]
                fz = Terms_fz[1]
                for Terms_f in sympy.Poly(fz,z).all_terms():
                    kk = Terms_f[0][0]
                    coef = Terms_f[1]
                    coefs_test[ii + 4 * jj + 16 * kk] = int(coef)
                
        passed = passed and np.array_equal(coefs_true, coefs_test)
        print(coefs_test)
    return passed

n_tests = 300
passed_flag = True
for i in range(n_tests):
    test_result = run_test()
    passed_flag = passed_flag and test_result
    if test_result:
        passfail = 'passed.'
    else:
        passfail = 'failed.'
    print('Random test %d/%d has '%(i,n_tests)+ passfail) 

if passed_flag:
    print('All tests have passed.')
else:
    print('Test did not pass.')

sys.exit(1)
#X = np.linspace(x0+(discard_x)*dx,x0+(Nx-2-discard_x)*dx,3000)
#Y = np.linspace(y0+(discard_y)*dy,y0+(Ny-2-discard_y)*dy,3000)
#Z = np.linspace(z0+(discard_z)*dz,z0+(Nz-2-discard_z)*dz,3000)
