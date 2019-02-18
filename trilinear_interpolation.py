import numpy as np
from tricubic_matrix import tricubicMat

class Trilinear_Interpolation:
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

        self.v = self.finite_differences(A, self.discard0, self.discard1, self.discard2)

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
    
        v = np.empty([n0,n1,n2,4],dtype=np.float64)
    
        #f#############################################################
        v[:,:,:,0] =         A[d0  :n0+d0  ,d1  :n1+d1  ,d2  :n2+d2  ] 
        ###############################################################
    
        #df/dx0########################################################
        v[:,:,:,1] =   1./self.dx*0.5*(+A[d0+1:n0+d0+1,d1  :n1+d1  ,d2  :n2+d2  ] 
                                 -A[d0-1:n0+d0-1,d1  :n1+d1  ,d2  :n2+d2  ])
        ###############################################################
    
        #df/dx1########################################################
        v[:,:,:,2] =   1./self.dy*0.5*(+A[d0  :n0+d0  ,d1+1:n1+d1+1,d2  :n2+d2  ] 
                                 -A[d0  :n0+d0  ,d1-1:n1+d1-1,d2  :n2+d2  ])
        ###############################################################
        
        #df/dx2########################################################
        v[:,:,:,3] =   1./self.dz*0.5*(+A[d0  :n0+d0  ,d1  :n1+d1  ,d2+1:n2+d2+1] 
                                 -A[d0  :n0+d0  ,d1  :n1+d1  ,d2-1:n2+d2-1])
        ###############################################################
    
        return v

    def interp(self, x, y, z, k):
        fx = (x - self.x0)/self.dx
        fy = (y - self.y0)/self.dy
        fz = (z - self.z0)/self.dz

        ix = int(fx)
        iy = int(fy)
        iz = int(fz)
        
        if ix < 1 or ix > self.v.shape[0]-1:
            raise Exception('Position is outside bounding box (first dimension)', x)

        if iy < 1 or iy > self.v.shape[1]-1:
            raise Exception('Position is outside bounding box (second dimension)', iy)

        if iz < 1 or iz > self.v.shape[2]-1:
            raise Exception('Position is outside bounding box (third dimension)', iz)

        x1 = fx - ix
        y1 = fy - iy
        z1 = fz - iz

        x2 = 1-x1
        y2 = 1-y1
        z2 = 1-z1

        ix -= 1
        iy -= 1
        iz -= 1

        c00=self.v[ix,iy,iz,k]*x2     + self.v[ix+1,iy,iz,k]*x1
        c01=self.v[ix,iy,iz+1,k]*x2   + self.v[ix+1,iy,iz+1,k]*x1
        c10=self.v[ix,iy+1,iz,k]*x2   + self.v[ix+1,iy+1,iz,k]*x1
        c11=self.v[ix,iy+1,iz+1,k]*x2 + self.v[ix+1,iy+1,iz+1,k]*x1
        c0 = c00*y2 + c10*y1
        c1 = c01*y2 + c11*y1
        c = c0*z2+c1*z1
        return c

    def val(self, x, y, z):
        return self.interp(x,y,z,0)

    def ddx(self, x, y, z):
        return self.interp(x,y,z,1)

    def ddy(self, x, y, z):
        return self.interp(x,y,z,2)

    def ddz(self, x, y, z):
        return self.interp(x,y,z,3)


