import numpy as np
import matplotlib.pyplot as plt
from tricubic_matrix import tricubicMat

def finite_differences(A, d0, d1, d2):
    n0 = A.shape[0] - 2 * d0
    n1 = A.shape[1] - 2 * d1
    n2 = A.shape[2] - 2 * d2

    v = np.empty([n0,n1,n2,8],dtype=np.float64)

    #f#############################################################
    v[:,:,:,0] =         A[d0  :n0+d0  ,d1  :n1+d1  ,d2  :n2+d2  ] 
    ###############################################################

    #df/dx0#########################################################
    v[:,:,:,1] =   0.5*(+A[d0+1:n0+d0+1,d1  :n1+d1  ,d2  :n2+d2  ] 
                        -A[d0-1:n0+d0-1,d1  :n1+d1  ,d2  :n2+d2  ])
    ###############################################################

    #df/dx1########################################################
    v[:,:,:,2] =   0.5*(+A[d0  :n0+d0  ,d1+1:n1+d1+1,d2  :n2+d2  ] 
                        -A[d0  :n0+d0  ,d1-1:n1+d1-1,d2  :n2+d2  ])
    ###############################################################
    
    #df/dx2########################################################
    v[:,:,:,3] =   0.5*(+A[d0  :n0+d0  ,d1  :n1+d1  ,d2+1:n2+d2+1] 
                        -A[d0  :n0+d0  ,d1  :n1+d1  ,d2-1:n2+d2-1])
    ###############################################################

    #df/dx0dx1#####################################################
    v[:,:,:,4] =  0.25*(+A[d0+1:n0+d0+1,d1+1:n1+d1+1,d2  :n2+d2  ] 
                        -A[d0+1:n0+d0+1,d1-1:n1+d1-1,d2  :n2+d2  ] 
                        -A[d0-1:n0+d0-1,d1+1:n1+d1+1,d2  :n2+d2  ]
                        +A[d0-1:n0+d0-1,d1-1:n1+d1-1,d2  :n2+d2  ])
    ###############################################################

    #df/dx0dx2#####################################################
    v[:,:,:,5] =  0.25*(+A[d0+1:n0+d0+1,d1  :n1+d1  ,d2+1:n2+d2+1] 
                        -A[d0+1:n0+d0+1,d1  :n1+d1  ,d2-1:n2+d2-1] 
                        -A[d0-1:n0+d0-1,d1  :n1+d1  ,d2+1:n2+d2+1]
                        +A[d0-1:n0+d0-1,d1  :n1+d1  ,d2-1:n2+d2-1])
    ###############################################################

    #df/dx1dx2#####################################################
    v[:,:,:,6] =  0.25*(+A[d0  :n0+d0  ,d1+1:n1+d1+1,d2+1:n2+d2+1] 
                        -A[d0  :n0+d0  ,d1-1:n1+d1-1,d2+1:n2+d2+1] 
                        -A[d0  :n0+d0  ,d1+1:n1+d1+1,d2-1:n2+d2-1]
                        +A[d0  :n0+d0  ,d1-1:n1+d1-1,d2-1:n2+d2-1])
    ###############################################################

    #df/dx1dx2#####################################################
    v[:,:,:,6] = 0.125*(+A[d0+1:n0+d0+1,d1+1:n1+d1+1,d2+1:n2+d2+1] 
                        -A[d0-1:n0+d0-1,d1+1:n1+d1+1,d2+1:n2+d2+1] 
                        -A[d0+1:n0+d0+1,d1-1:n1+d1-1,d2+1:n2+d2+1] 
                        +A[d0-1:n0+d0-1,d1-1:n1+d1-1,d2+1:n2+d2+1]
                        -A[d0+1:n0+d0+1,d1+1:n1+d1+1,d2-1:n2+d2-1] 
                        +A[d0-1:n0+d0-1,d1+1:n1+d1+1,d2-1:n2+d2-1] 
                        +A[d0+1:n0+d0+1,d1-1:n1+d1-1,d2+1:n2+d2+1]
                        -A[d0-1:n0+d0-1,d1-1:n1+d1-1,d2-1:n2+d2-1])
    ###############################################################

    return v

def tricubic_coefs(v):
    
    n0 = v.shape[0] - 1
    n1 = v.shape[1] - 1
    n2 = v.shape[2] - 1
    
    coefs = np.empty([n0,n1,n2,64], dtype= np.float64)
    b = np.empty([64], dtype = np.float64)
    print('Calculating Coefficients.')
    for i in range(n0):
        if i%10 == 0: print('%d/%d'%(i,n0))
        for j in range(n1):
            for k in range(n2):
                for l in range(8):
                    b[8*l+0]= v[i  ,j  ,k  ,l]
                    b[8*l+1]= v[i+1,j  ,k  ,l]
                    b[8*l+2]= v[i  ,j+1,k  ,l]
                    b[8*l+3]= v[i+1,j+1,k  ,l]
                    b[8*l+4]= v[i  ,j  ,k+1,l]
                    b[8*l+5]= v[i+1,j  ,k+1,l]
                    b[8*l+6]= v[i  ,j+1,k+1,l]
                    b[8*l+7]= v[i+1,j+1,k+1,l]
                coefs[i,j,k,:] = np.matmul(tricubicMat, b)
    print('Coefficients are calculated.')

    return coefs

def fval(x,y,z,coefs):
    dz=1
    result = 0
    i = 0

    for k in range(4):
        dy = 1
        for j in range(4):
            result += dy*dz * ( coefs[i] + x * ( coefs[i+1] + x * ( coefs[i+2] + x * coefs[i+3])))
            i += 4
            dy *= y
        dz *= z
    return result


A = np.empty([10,10,10], np.float64)
for i in range(10):
    for j in range(10):
        for k in range(10):
            A[i,j,k] = np.cos(np.power(0.1*i*j*k,1)+0.1)+10*np.exp(-i*j*k/300.)

x0 = 3.
y0 = 4.
z0 = 4.

v = finite_differences(A,1,1,1)
coefs = tricubic_coefs(v)

ix = int(x0)
iy = int(y0)
iz = int(z0)

print(A[ix,iy,iz])
print(fval(x0-int(ix),y0-int(iy),z0-int(iz),coefs[ix-1,iy-1,iz-1]))

tx = np.linspace(1,6,300)
dx = tx[1]-tx[0]
ty = np.cos(np.power(0.1*tx*iy*iz,1)+0.1)+10*np.exp(-tx*iy*iz/300.)
tyi = np.empty_like(tx)
for i,xx in enumerate(tx):
    tyi[i]= fval(xx-int(xx),y0-iy,z0-iz,coefs[int(xx)-1,iy-1,iz-1])
plt.plot(A[:,iy,iz],'ko')
plt.plot(tx,ty,'r-')
plt.plot(tx,tyi,'g-')
plt.show(False)
input()







