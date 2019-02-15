import numpy as np
import matplotlib.pyplot as plt
from tricubic_matrix import tricubicMat

class Tricubic_Interpolation:
    def __init__(self):
        self.discard0 = 1
        self.discard1 = 1
        self.discard2 = 1
        
        self.x0 = 0. 
        self.y0 = 0. 
        self.z0 = 0.

        self.dx = 1.
        self.dy = 1.
        self.dz = 1.

    def initialize(self, A):
        if A.shape[0]- 2*self.discard0 < 2:
            raise Exception('n0 < 2: Interpolating array is too small along the first dimension (after discards).')
        if A.shape[1]- 2*self.discard1 < 2:
            raise Exception('n1 < 2: Interpolating array is too small along the second dimension (after discards).')
        if A.shape[2]- 2*self.discard2 < 2:
            raise Exception('n2 < 2: Interpolating array is too small along the third dimension (after discards).')

        v = self.finite_differences(A, self.discard0, self.discard1, self.discard2)
        self.coefs = self.tricubic_coefs(v)

    def set_origin(self, x0, y0, z0):
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0

    def set_steps(self, dx, dy, dz):
        self.dx = dx
        self.dy = dy
        self.dz = dz

    def discard_points(self, discard0, discard1, discard2):
        self.discard0 = discard0
        self.discard1 = discard1
        self.discard2 = discard2

    def finite_differences(self, A, d0, d1, d2):
        n0 = A.shape[0] - 2 * d0
        n1 = A.shape[1] - 2 * d1
        n2 = A.shape[2] - 2 * d2
    
        v = np.empty([n0,n1,n2,8],dtype=np.float64)
    
        #f#############################################################
        v[:,:,:,0] =         A[d0  :n0+d0  ,d1  :n1+d1  ,d2  :n2+d2  ] 
        ###############################################################
    
        #df/dx0########################################################
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

    def tricubic_coefs(self, v):
        
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


    def ddx(self, x, y, z):

        fx = (x - self.x0)/self.dx
        fy = (y - self.y0)/self.dy
        fz = (z - self.z0)/self.dz

        ix = int(fx)
        iy = int(fy)
        iz = int(fz)
        
        if ix < 1 or ix > self.coefs.shape[0]:
            raise Exception('Position is outside bounding box (first dimension)')

        if iy < 1 or iy > self.coefs.shape[1]:
            raise Exception('Position is outside bounding box (second dimension)')

        if iz < 1 or iz > self.coefs.shape[2]:
            raise Exception('Position is outside bounding box (third dimension)')

        coefs = self.coefs[ix-1,iy-1,iz-1,:]

        x1 = fx - ix
        y1 = fy - iy
        z1 = fz - iz

        dz=1
        result = 0
        i = 0
    
        xv=np.empty([4],dtype=np.float64)
        yv=np.empty([4],dtype=np.float64)
        zv=np.empty([4],dtype=np.float64)
        xv[0]=0
        yv[0]=1
        zv[0]=1
        xv[1]=1
        yv[1]=y1
        zv[1]=z1
        xv[2]=2*x1*xv[1]
        yv[2]=y1*yv[1]
        zv[2]=z1*zv[1]
        xv[3]=1.5*x1*xv[2]
        yv[3]=y1*yv[2]
        zv[3]=z1*zv[2]
        res = np.matmul(coefs.reshape(4,4,4),xv)
        res = np.matmul(res,yv)
        res = np.matmul(res,zv)
        return res/self.dx

    def ddy(self, x, y, z):

        fx = (x - self.x0)/self.dx
        fy = (y - self.y0)/self.dy
        fz = (z - self.z0)/self.dz

        ix = int(fx)
        iy = int(fy)
        iz = int(fz)
        
        if ix < 1 or ix > self.coefs.shape[0]:
            raise Exception('Position is outside bounding box (first dimension)')

        if iy < 1 or iy > self.coefs.shape[1]:
            raise Exception('Position is outside bounding box (second dimension)')

        if iz < 1 or iz > self.coefs.shape[2]:
            raise Exception('Position is outside bounding box (third dimension)')

        coefs = self.coefs[ix-1,iy-1,iz-1,:]

        x1 = fx - ix
        y1 = fy - iy
        z1 = fz - iz

        dz=1
        result = 0
        i = 0
    
        xv=np.empty([4],dtype=np.float64)
        yv=np.empty([4],dtype=np.float64)
        zv=np.empty([4],dtype=np.float64)
        xv[0]=1
        yv[0]=0
        zv[0]=1
        xv[1]=x1
        yv[1]=1
        zv[1]=z1
        xv[2]=x1*xv[1]
        yv[2]=2*y1*yv[1]
        zv[2]=z1*zv[1]
        xv[3]=x1*xv[2]
        yv[3]=1.5*y1*yv[2]
        zv[3]=z1*zv[2]
        res = np.matmul(coefs.reshape(4,4,4),xv)
        res = np.matmul(res,yv)
        res = np.matmul(res,zv)
        return res/self.dy

    def ddz(self, x, y, z):

        fx = (x - self.x0)/self.dx
        fy = (y - self.y0)/self.dy
        fz = (z - self.z0)/self.dz

        ix = int(fx)
        iy = int(fy)
        iz = int(fz)
        
        if ix < 1 or ix > self.coefs.shape[0]:
            raise Exception('Position is outside bounding box (first dimension)')

        if iy < 1 or iy > self.coefs.shape[1]:
            raise Exception('Position is outside bounding box (second dimension)')

        if iz < 1 or iz > self.coefs.shape[2]:
            raise Exception('Position is outside bounding box (third dimension)')

        coefs = self.coefs[ix-1,iy-1,iz-1,:]

        x1 = fx - ix
        y1 = fy - iy
        z1 = fz - iz

        dz=1
        result = 0
        i = 0
    
        xv=np.empty([4],dtype=np.float64)
        yv=np.empty([4],dtype=np.float64)
        zv=np.empty([4],dtype=np.float64)
        xv[0]=1
        yv[0]=1
        zv[0]=0
        xv[1]=x1
        yv[1]=y1
        zv[1]=1
        xv[2]=x1*xv[1]
        yv[2]=y1*yv[1]
        zv[2]=2*z1*zv[1]
        xv[3]=x1*xv[2]
        yv[3]=y1*yv[2]
        zv[3]=1.5*z1*zv[2]
        res = np.matmul(coefs.reshape(4,4,4),xv)
        res = np.matmul(res,yv)
        res = np.matmul(res,zv)
        return res/self.dz

    def val(self, x, y, z):

        fx = (x - self.x0)/self.dx
        fy = (y - self.y0)/self.dy
        fz = (z - self.z0)/self.dz

        ix = int(fx)
        iy = int(fy)
        iz = int(fz)
        
        if ix < 1 or ix > self.coefs.shape[0]:
            raise Exception('Position is outside bounding box (first dimension)')

        if iy < 1 or iy > self.coefs.shape[1]:
            raise Exception('Position is outside bounding box (second dimension)')

        if iz < 1 or iz > self.coefs.shape[2]:
            raise Exception('Position is outside bounding box (third dimension)')

        coefs = self.coefs[ix-1,iy-1,iz-1,:]

        x1 = fx - ix
        y1 = fy - iy
        z1 = fz - iz

        dz=1
        result = 0
        i = 0
    
        xv=np.empty([4],dtype=np.float64)
        yv=np.empty([4],dtype=np.float64)
        zv=np.empty([4],dtype=np.float64)
        xv[0]=1
        yv[0]=1
        zv[0]=1
        xv[1]=x1
        yv[1]=y1
        zv[1]=z1
        xv[2]=x1*xv[1]
        yv[2]=y1*yv[1]
        zv[2]=z1*zv[1]
        xv[3]=x1*xv[2]
        yv[3]=y1*yv[2]
        zv[3]=z1*zv[2]
        res = np.matmul(coefs.reshape(4,4,4),xv)
        res = np.matmul(res,yv)
        res = np.matmul(res,zv)
        return res

    def val2(self, x, y, z):

        fx = (x - self.x0)/self.dx
        fy = (y - self.y0)/self.dy
        fz = (z - self.z0)/self.dz

        ix = int(fx)
        iy = int(fy)
        iz = int(fz)
        
        coefs = self.coefs[ix-1,iy-1,iz-1,:]

        x1 = fx - ix
        y1 = fy - iy
        z1 = fz - iz

        dz=1
        result = 0
        i = 0
    
        xv=np.empty([4],dtype=np.float64)
        yv=np.empty([4],dtype=np.float64)
        zv=np.empty([4],dtype=np.float64)
        xv[0]=1
        yv[0]=1
        zv[0]=1
        xv[1]=x1
        yv[1]=y1
        zv[1]=z1
        xv[2]=x1*xv[1]
        yv[2]=y1*yv[1]
        zv[2]=z1*zv[1]
        xv[3]=x1*xv[2]
        yv[3]=y1*yv[2]
        zv[3]=z1*zv[2]
        X = np.tensordot(np.tensordot(zv,yv,axes=0),xv,axes=0).reshape(64)
        return np.tensordot(coefs,X,axes=1)

    def val3(self, x, y, z):

        fx = (x - self.x0)/self.dx
        fy = (y - self.y0)/self.dy
        fz = (z - self.z0)/self.dz

        ix = int(fx)
        iy = int(fy)
        iz = int(fz)
        
        coefs = self.coefs[ix-1,iy-1,iz-1,:]

        x1 = fx - ix
        y1 = fy - iy
        z1 = fz - iz

        dz=1
        result = 0
        i = 0
    
        for k in range(4):
            dy = 1
            for j in range(4):
                result += dy*dz * ( coefs[i] + x1 * ( coefs[i+1] + x1 * ( coefs[i+2] + x1 * coefs[i+3])))
                i += 4
                dy *= y1
            dz *= z1
        return result

#cc=0.050
#
#Nx=100
#A = np.empty([Nx,10,10])
#xa = np.empty([Nx])
#for i in range(Nx):
#    for j in range(10):
#        for k in range(10):
#            x=cc*i
#            y=j
#            z=k
#            #A[i,j,k] = 5*x**3*y**3+x**3*z**3 
#            A[i,j,k] = np.cos(np.power(0.1*x*y*z,1)+0.1)+10*np.exp(-x*y*z/300.)
#    xa[i]=x
#B= np.empty_like(xa)
#B[0]=0
#B[9]=0
#x0 = 3.
#y0 = 4.
#z0 = 4.
#
#ip = Tricubic_Interpolation()
#ip.discard0 = 1
#ip.dx = cc
#ip.initialize(A)
#coefs = ip.coefs
##v = finite_differences(A,1,1,1)
##coefs = tricubic_coefs(v)
#
#ix = int(x0)
#iy = int(y0)
#iz = int(z0)
#
#print(A[ix,iy,iz])
#for i in range(1,Nx-1):
#    B[i]=(A[i+1,iy,iz]-A[i-1,iy,iz])/cc/2.
#print(ip.val(x0,y0,z0))
#
#tx = np.linspace(0.2,3,3000)
#ty = np.cos(np.power(0.1*tx*iy*iz,1)+0.1)+10*np.exp(-tx*iy*iz/300.)
#tyi = np.empty_like(tx)
#tyid = np.empty_like(tx)
#for i,xx in enumerate(tx):
#    tyi[i]= ip.val(xx,y0,z0)
#    tyid[i]= ip.ddx(xx,y0+0.5,z0)
#fig1 = plt.figure(1)
#ax1 = fig1.add_subplot(111)
#ax1.plot(xa,A[:,iy,iz],'ko')
#ax1.plot(tx,ty,'r-')
#ax1.plot(tx,tyi,'g-')
#fig2 = plt.figure(2)
#ax2 = fig2.add_subplot(111)
#ax2.plot(xa,B[:],'ko')
#ax2.plot(tx,tyid,'g-')
#plt.show(False)
#input()
#






