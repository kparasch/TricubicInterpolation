import sys
sys.path.append('..')

from tricubic_interpolation import Tricubic_Interpolation
import matplotlib.pyplot as plt
import numpy as np
import sympy
#plt.style.use('kostas')

def run_test(debug=False):
    x,y,z = sympy.symbols('x y z')
    upper_int = 9
    #print('One random number: %d'%np.random.randint(upper_int))
    coefs_true = np.random.randint(-upper_int,upper_int, size=64)
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
    dx = 1
    dy = 1
    dz = 1
    discard_x = 1
    discard_y = 1
    discard_z = 1
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
    
    ip = Tricubic_Interpolation(B, x0, y0, z0, dx, dy, dz, discard_x, discard_y, discard_z, 'Exact')
    
    passed = True
    n = 1
    for test in range(n):
        xv = np.random.rand()*(Nx-3-2*discard_x)*dx + x0 +discard_x*dx
        yv = np.random.rand()*(Ny-3-2*discard_y)*dy + y0 +discard_y*dy
        zv = np.random.rand()*(Nz-3-2*discard_z)*dz + z0 +discard_z*dz
        ix,iy,iz = ip.coords_to_indices(xv,yv,zv)
        coefs = ip.get_coefs(ip.construct_b(ix, iy, iz))
        fxyz = sympy.simplify(sum([ (coefs[i + 4*j + 16*k]) * ((x-ix)/dx)**i * ((y-iy)/dy)**j * ((z-iz)/dz)**k for i in range(4) for j in range(4) for k in range(4)]))
        if debug:
            print('Output (transformed) polynomial:')
            print(fxyz)
            print('Input coefficients:')
            print(coefs_true)
        
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
                    coefs_test[ii + 4 * jj + 16 * kk] = coef#int(coef)
                
        passed = passed and np.array_equal(coefs_true, coefs_test)
        if debug:
            print('Output coefficients:')
            print(coefs_test)
    return passed


n_tests = 5
debug = False


passed_flag = True
for i in range(n_tests):
    test_result = run_test(debug=debug)
    passed_flag = passed_flag and test_result
    if test_result:
        passfail = 'passed.'
    else:
        passfail = 'failed.'
    print('Random test %d/%d has '%(i,n_tests)+ passfail) 

if passed_flag:
    print('All tests have passed.')
    sys.exit(0)
else:
    print('Test did not pass.')
    sys.exit(1)

