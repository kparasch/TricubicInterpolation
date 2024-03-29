import numpy as np

class Trilinear_Interpolation:
    def __init__(self, A, x0=0., y0=0., z0=0., dx=1., dy=1., dz=1., discardx=1, discardy=1, discardz=1, method='FD'):
        self.discardx = discardx
        self.discardy = discardy
        self.discardz = discardz

        self.dx = dx
        self.dy = dy
        self.dz = dz

        self.x0 = x0 + (discardx-1)*dx
        self.y0 = y0 + (discardy-1)*dx
        self.z0 = z0 + (discardz-1)*dx

        if A.shape[0]- 2*self.discardx < 2:
            raise Exception('n0 < 2: Interpolating array is too small along the first dimension (after discards).')
        if A.shape[1]- 2*self.discardy < 2:
            raise Exception('n1 < 2: Interpolating array is too small along the second dimension (after discards).')
        if A.shape[2]- 2*self.discardz < 2:
            raise Exception('n2 < 2: Interpolating array is too small along the third dimension (after discards).')

        self.v = self.finite_differences(A, self.discardx, self.discardy, self.discardz)

    def finite_differences(self, A, d0, d1, d2):
        n0 = A.shape[0] - 2 * d0
        n1 = A.shape[1] - 2 * d1
        n2 = A.shape[2] - 2 * d2
    
        v = np.empty([n0,n1,n2,4],dtype=float)
    
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
            raise Exception('Position is outside bounding box (first dimension)', ix)

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


