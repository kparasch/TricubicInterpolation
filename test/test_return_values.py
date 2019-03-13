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

    x0 = -11
    y0 = -32
    z0 = -19
    dx = 2.
    dy = 6.1
    dz = 4.2
    discard_x = 1
    discard_y = 1
    discard_z = 1
    Nx = 10
    Ny = 10
    Nz = 10

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
    
    ip = Tricubic_Interpolation(B, x0, y0, z0, dx, dy, dz, discard_x, discard_y, discard_z, 'Exact')
    
    passed = True
    n = 10
    for test in range(n):
        xv = np.random.rand()*(Nx-3-2*discard_x)*dx + x0 +discard_x*dx
        yv = np.random.rand()*(Ny-3-2*discard_y)*dy + y0 +discard_y*dy
        zv = np.random.rand()*(Nz-3-2*discard_z)*dz + z0 +discard_z*dz
        
        output_true  = np.array([ lamf(xv,yv,zv),
                                  lamdfdx(xv,yv,zv),
                                  lamdfdy(xv,yv,zv),
                                  lamdfdz(xv,yv,zv)
                                ])

        output_test = np.array([ ip.val(xv,yv,zv),
                                 ip.ddx(xv,yv,zv),
                                 ip.ddy(xv,yv,zv),
                                 ip.ddz(xv,yv,zv)
                               ])

        if debug:       
            for i in range(len(output_true)):
                print('%f %f'%(output_true[i], output_test[i])) 
                
        # np.allclose(a,b, rtol, atol) checks if: absolute(a - b) <= (atol + rtol * absolute(b))
        passed = passed and np.allclose(output_test, output_true, rtol=1.e-10, atol=1.e-10)

    return passed

n_tests = 15
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

