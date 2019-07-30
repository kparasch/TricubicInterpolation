import numpy as np
from tricubic_matrix import tricubicMat


class Tricubic_Interpolation(object):
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

        if method=='Exact':
            #print('Using exact derivatives.')
            if len(A.shape) != 4:
                raise Exception('Input array should be 4-dimensional when using exact derivatives method. It\'s not.')
            self.construct_b = self.exact_diff
        elif method=='FD':
            #print('Using finite difference approximation for derivatives.')
            if len(A.shape) != 3:
                raise Exception('Input array should be 3-dimensional when using finite differences method. It\'s not.')
            self.construct_b = self.finite_diff
        else:
            raise ValueError('Invalid method: %s'%method)


        self.A = A[:,:,:]
        if self.discardx > 1:
            self.A = self.A[self.discardx-1:-(self.discardx-1),:,:]
        if self.discardy > 1:
            self.A = self.A[:,self.discardy-1:-(self.discardy-1),:]
        if self.discardz > 1:
            self.A = self.A[:,:,self.discardz-1:-(self.discardz-1)]
        
        self.ix_bound_up = self.A.shape[0] - 3  #    a -1 because counting starts from 0,
        self.iy_bound_up = self.A.shape[1] - 3  #    another -1 because one more point is needed for finite differences,
        self.iz_bound_up = self.A.shape[2] - 3  #    and a last -1 because the bound corresponds to the bound 
                                                #    for the lower index inclusive

        self.ix_bound_low = 1
        self.iy_bound_low = 1
        self.iz_bound_low = 1

        if self.ix_bound_up <= self.ix_bound_low:
            raise Exception('Interpolating array is too small along the first dimension (after discards).')
        if self.iy_bound_up <= self.iy_bound_low:
            raise Exception('Interpolating array is too small along the second dimension (after discards).')
        if self.iz_bound_up <= self.iz_bound_low:
            raise Exception('Interpolating array is too small along the third dimension (after discards).')

        #if A.shape[0]- 2*self.discardx < 2:
        #    raise Exception('n0 < 2: Interpolating array is too small along the first dimension (after discards).')
        #if A.shape[1]- 2*self.discardy < 2:
        #    raise Exception('n1 < 2: Interpolating array is too small along the second dimension (after discards).')
        #if A.shape[2]- 2*self.discardz < 2:
        #    raise Exception('n2 < 2: Interpolating array is too small along the third dimension (after discards).')


    def exact_diff(self, ix, iy, iz): 
        
        #scaling factors for derivatives
        scale = [1, self.dx, self.dy, self.dz, self.dx*self.dy,
                 self.dx*self.dz, self.dy*self.dz,
                 self.dx*self.dy*self.dz
                ]

        b = np.empty([64],dtype=np.float64)
        
        for l in range(8):
            b[8*l+0] = self.A[ix  ,iy  ,iz  ,l]*scale[l]
            b[8*l+1] = self.A[ix+1,iy  ,iz  ,l]*scale[l]
            b[8*l+2] = self.A[ix  ,iy+1,iz  ,l]*scale[l]
            b[8*l+3] = self.A[ix+1,iy+1,iz  ,l]*scale[l]
            b[8*l+4] = self.A[ix  ,iy  ,iz+1,l]*scale[l]
            b[8*l+5] = self.A[ix+1,iy  ,iz+1,l]*scale[l]
            b[8*l+6] = self.A[ix  ,iy+1,iz+1,l]*scale[l]
            b[8*l+7] = self.A[ix+1,iy+1,iz+1,l]*scale[l]

        return b


    def finite_diff(self, ix, iy, iz):
        b = np.empty([64],dtype=np.float64)
        
        b[ 0] = self.A[ix  ,iy  ,iz  ]
        b[ 1] = self.A[ix+1,iy  ,iz  ]
        b[ 2] = self.A[ix  ,iy+1,iz  ]
        b[ 3] = self.A[ix+1,iy+1,iz  ]
        b[ 4] = self.A[ix  ,iy  ,iz+1]
        b[ 5] = self.A[ix+1,iy  ,iz+1]
        b[ 6] = self.A[ix  ,iy+1,iz+1]
        b[ 7] = self.A[ix+1,iy+1,iz+1]

        b[ 8] = 0.5*(self.A[ix+1,iy  ,iz  ] - self.A[ix-1,iy  ,iz  ]) 
        b[ 9] = 0.5*(self.A[ix+2,iy  ,iz  ] - self.A[ix  ,iy  ,iz  ])
        b[10] = 0.5*(self.A[ix+1,iy+1,iz  ] - self.A[ix-1,iy+1,iz  ])
        b[11] = 0.5*(self.A[ix+2,iy+1,iz  ] - self.A[ix  ,iy+1,iz  ])
        b[12] = 0.5*(self.A[ix+1,iy  ,iz+1] - self.A[ix-1,iy  ,iz+1])
        b[13] = 0.5*(self.A[ix+2,iy  ,iz+1] - self.A[ix  ,iy  ,iz+1])
        b[14] = 0.5*(self.A[ix+1,iy+1,iz+1] - self.A[ix-1,iy+1,iz+1])
        b[15] = 0.5*(self.A[ix+2,iy+1,iz+1] - self.A[ix  ,iy+1,iz+1])

        b[16] = 0.5*(self.A[ix  ,iy+1,iz  ] - self.A[ix  ,iy-1,iz  ]) 
        b[17] = 0.5*(self.A[ix+1,iy+1,iz  ] - self.A[ix+1,iy-1,iz  ])
        b[18] = 0.5*(self.A[ix  ,iy+2,iz  ] - self.A[ix  ,iy  ,iz  ])
        b[19] = 0.5*(self.A[ix+1,iy+2,iz  ] - self.A[ix+1,iy  ,iz  ])
        b[20] = 0.5*(self.A[ix  ,iy+1,iz+1] - self.A[ix  ,iy-1,iz+1])
        b[21] = 0.5*(self.A[ix+1,iy+1,iz+1] - self.A[ix+1,iy-1,iz+1])
        b[22] = 0.5*(self.A[ix  ,iy+2,iz+1] - self.A[ix  ,iy  ,iz+1])
        b[23] = 0.5*(self.A[ix+1,iy+2,iz+1] - self.A[ix+1,iy  ,iz+1])

        b[24] = 0.5*(self.A[ix  ,iy  ,iz+1] - self.A[ix  ,iy  ,iz-1]) 
        b[25] = 0.5*(self.A[ix+1,iy  ,iz+1] - self.A[ix+1,iy  ,iz-1])
        b[26] = 0.5*(self.A[ix  ,iy+1,iz+1] - self.A[ix  ,iy+1,iz-1])
        b[27] = 0.5*(self.A[ix+1,iy+1,iz+1] - self.A[ix+1,iy+1,iz-1])
        b[28] = 0.5*(self.A[ix  ,iy  ,iz+2] - self.A[ix  ,iy  ,iz  ])
        b[29] = 0.5*(self.A[ix+1,iy  ,iz+2] - self.A[ix+1,iy  ,iz  ])
        b[30] = 0.5*(self.A[ix  ,iy+1,iz+2] - self.A[ix  ,iy+1,iz  ])
        b[31] = 0.5*(self.A[ix+1,iy+1,iz+2] - self.A[ix+1,iy+1,iz  ])

        b[32] = 0.25*(self.A[ix+1,iy+1,iz  ] - self.A[ix-1,iy+1,iz  ] - self.A[ix+1,iy-1,iz  ] + self.A[ix-1,iy-1,iz  ])
        b[33] = 0.25*(self.A[ix+2,iy+1,iz  ] - self.A[ix  ,iy+1,iz  ] - self.A[ix+2,iy-1,iz  ] + self.A[ix  ,iy-1,iz  ])
        b[34] = 0.25*(self.A[ix+1,iy+2,iz  ] - self.A[ix-1,iy+2,iz  ] - self.A[ix+1,iy  ,iz  ] + self.A[ix-1,iy  ,iz  ])
        b[35] = 0.25*(self.A[ix+2,iy+2,iz  ] - self.A[ix  ,iy+2,iz  ] - self.A[ix+2,iy  ,iz  ] + self.A[ix  ,iy  ,iz  ])
        b[36] = 0.25*(self.A[ix+1,iy+1,iz+1] - self.A[ix-1,iy+1,iz+1] - self.A[ix+1,iy-1,iz+1] + self.A[ix-1,iy-1,iz+1])
        b[37] = 0.25*(self.A[ix+2,iy+1,iz+1] - self.A[ix  ,iy+1,iz+1] - self.A[ix+2,iy-1,iz+1] + self.A[ix  ,iy-1,iz+1])
        b[38] = 0.25*(self.A[ix+1,iy+2,iz+1] - self.A[ix-1,iy+2,iz+1] - self.A[ix+1,iy  ,iz+1] + self.A[ix-1,iy  ,iz+1])
        b[39] = 0.25*(self.A[ix+2,iy+2,iz+1] - self.A[ix  ,iy+2,iz+1] - self.A[ix+2,iy  ,iz+1] + self.A[ix  ,iy  ,iz+1])

        b[40] = 0.25*(self.A[ix+1,iy  ,iz+1] - self.A[ix-1,iy  ,iz+1] - self.A[ix+1,iy  ,iz-1] + self.A[ix-1,iy  ,iz-1])
        b[41] = 0.25*(self.A[ix+2,iy  ,iz+1] - self.A[ix  ,iy  ,iz+1] - self.A[ix+2,iy  ,iz-1] + self.A[ix  ,iy  ,iz-1])
        b[42] = 0.25*(self.A[ix+1,iy+1,iz+1] - self.A[ix-1,iy+1,iz+1] - self.A[ix+1,iy+1,iz-1] + self.A[ix-1,iy+1,iz-1])
        b[43] = 0.25*(self.A[ix+2,iy+1,iz+1] - self.A[ix  ,iy+1,iz+1] - self.A[ix+2,iy+1,iz-1] + self.A[ix  ,iy+1,iz-1])
        b[44] = 0.25*(self.A[ix+1,iy  ,iz+2] - self.A[ix-1,iy  ,iz+2] - self.A[ix+1,iy  ,iz  ] + self.A[ix-1,iy  ,iz  ])
        b[45] = 0.25*(self.A[ix+2,iy  ,iz+2] - self.A[ix  ,iy  ,iz+2] - self.A[ix+2,iy  ,iz  ] + self.A[ix  ,iy  ,iz  ])
        b[46] = 0.25*(self.A[ix+1,iy+1,iz+2] - self.A[ix-1,iy+1,iz+2] - self.A[ix+1,iy+1,iz  ] + self.A[ix-1,iy+1,iz  ])
        b[47] = 0.25*(self.A[ix+2,iy+1,iz+2] - self.A[ix  ,iy+1,iz+2] - self.A[ix+2,iy+1,iz  ] + self.A[ix  ,iy+1,iz  ])

        b[48] = 0.25*(self.A[ix  ,iy+1,iz+1] - self.A[ix  ,iy-1,iz+1] - self.A[ix  ,iy+1,iz-1] + self.A[ix  ,iy-1,iz-1])
        b[49] = 0.25*(self.A[ix+1,iy+1,iz+1] - self.A[ix+1,iy-1,iz+1] - self.A[ix+1,iy+1,iz-1] + self.A[ix+1,iy-1,iz-1])
        b[50] = 0.25*(self.A[ix  ,iy+2,iz+1] - self.A[ix  ,iy  ,iz+1] - self.A[ix  ,iy+2,iz-1] + self.A[ix  ,iy  ,iz-1])
        b[51] = 0.25*(self.A[ix+1,iy+2,iz+1] - self.A[ix+1,iy  ,iz+1] - self.A[ix+1,iy+2,iz-1] + self.A[ix+1,iy  ,iz-1])
        b[52] = 0.25*(self.A[ix  ,iy+1,iz+2] - self.A[ix  ,iy-1,iz+2] - self.A[ix  ,iy+1,iz  ] + self.A[ix  ,iy-1,iz  ])
        b[53] = 0.25*(self.A[ix+1,iy+1,iz+2] - self.A[ix+1,iy-1,iz+2] - self.A[ix+1,iy+1,iz  ] + self.A[ix+1,iy-1,iz  ])
        b[54] = 0.25*(self.A[ix  ,iy+2,iz+2] - self.A[ix  ,iy  ,iz+2] - self.A[ix  ,iy+2,iz  ] + self.A[ix  ,iy  ,iz  ])
        b[55] = 0.25*(self.A[ix+1,iy+2,iz+2] - self.A[ix+1,iy  ,iz+2] - self.A[ix+1,iy+2,iz  ] + self.A[ix+1,iy  ,iz  ])

        b[56] = 0.125*(self.A[ix+1,iy+1,iz+1] - self.A[ix-1,iy+1,iz+1] - self.A[ix+1,iy-1,iz+1] + self.A[ix-1,iy-1,iz+1] - self.A[ix+1,iy+1,iz-1] + self.A[ix-1,iy+1,iz-1] + self.A[ix+1,iy-1,iz-1] - self.A[ix-1,iy-1,iz-1])
        b[57] = 0.125*(self.A[ix+2,iy+1,iz+1] - self.A[ix  ,iy+1,iz+1] - self.A[ix+2,iy-1,iz+1] + self.A[ix  ,iy-1,iz+1] - self.A[ix+2,iy+1,iz-1] + self.A[ix  ,iy+1,iz-1] + self.A[ix+2,iy-1,iz-1] - self.A[ix  ,iy-1,iz-1])
        b[58] = 0.125*(self.A[ix+1,iy+2,iz+1] - self.A[ix-1,iy+2,iz+1] - self.A[ix+1,iy  ,iz+1] + self.A[ix-1,iy  ,iz+1] - self.A[ix+1,iy+2,iz-1] + self.A[ix-1,iy+2,iz-1] + self.A[ix+1,iy  ,iz-1] - self.A[ix-1,iy  ,iz-1])
        b[59] = 0.125*(self.A[ix+2,iy+2,iz+1] - self.A[ix  ,iy+2,iz+1] - self.A[ix+2,iy  ,iz+1] + self.A[ix  ,iy  ,iz+1] - self.A[ix+2,iy+2,iz-1] + self.A[ix  ,iy+2,iz-1] + self.A[ix+2,iy  ,iz-1] - self.A[ix  ,iy  ,iz-1])
        b[60] = 0.125*(self.A[ix+1,iy+1,iz+2] - self.A[ix-1,iy+1,iz+2] - self.A[ix+1,iy-1,iz+2] + self.A[ix-1,iy-1,iz+2] - self.A[ix+1,iy+1,iz  ] + self.A[ix-1,iy+1,iz  ] + self.A[ix+1,iy-1,iz  ] - self.A[ix-1,iy-1,iz  ])
        b[61] = 0.125*(self.A[ix+2,iy+1,iz+2] - self.A[ix  ,iy+1,iz+2] - self.A[ix+2,iy-1,iz+2] + self.A[ix  ,iy-1,iz+2] - self.A[ix+2,iy+1,iz  ] + self.A[ix  ,iy+1,iz  ] + self.A[ix+2,iy-1,iz  ] - self.A[ix  ,iy-1,iz  ])
        b[62] = 0.125*(self.A[ix+1,iy+2,iz+2] - self.A[ix-1,iy+2,iz+2] - self.A[ix+1,iy  ,iz+2] + self.A[ix-1,iy  ,iz+2] - self.A[ix+1,iy+2,iz  ] + self.A[ix-1,iy+2,iz  ] + self.A[ix+1,iy  ,iz  ] - self.A[ix-1,iy  ,iz  ])
        b[63] = 0.125*(self.A[ix+2,iy+2,iz+2] - self.A[ix  ,iy+2,iz+2] - self.A[ix+2,iy  ,iz+2] + self.A[ix  ,iy  ,iz+2] - self.A[ix+2,iy+2,iz  ] + self.A[ix  ,iy+2,iz  ] + self.A[ix+2,iy  ,iz  ] - self.A[ix  ,iy  ,iz  ])

        return b


    def coords_to_indices(self, x, y, z):
        
        fx = (x - self.x0)/self.dx
        fy = (y - self.y0)/self.dy
        fz = (z - self.z0)/self.dz

        ix = int(fx)
        iy = int(fy)
        iz = int(fz)

        return ix, iy, iz

    
    def is_inside_box(self, x, y, z):
        fx = (x - self.x0)/self.dx
        fy = (y - self.y0)/self.dy
        fz = (z - self.z0)/self.dz

        ix = int(fx)
        iy = int(fy)
        iz = int(fz)

        x1 = fx - ix
        y1 = fy - iy
        z1 = fz - iz

        inside_box = True
        if   ix < self.ix_bound_low or ix > self.ix_bound_up:
            inside_box = False
        elif iy < self.iy_bound_low or iy > self.iy_bound_up:
            inside_box = False
        elif iz < self.iz_bound_low or iz > self.iz_bound_up:
            inside_box = False

        return inside_box

    def coords_to_indices_and_floats(self, x, y, z):
        
        fx = (x - self.x0)/self.dx
        fy = (y - self.y0)/self.dy
        fz = (z - self.z0)/self.dz

        ix = int(fx)
        iy = int(fy)
        iz = int(fz)

        x1 = fx - ix
        y1 = fy - iy
        z1 = fz - iz

        inside_box = True
        if   ix < self.ix_bound_low or ix > self.ix_bound_up:
            inside_box = False
        elif iy < self.iy_bound_low or iy > self.iy_bound_up:
            inside_box = False
        elif iz < self.iz_bound_low or iz > self.iz_bound_up:
            inside_box = False

        if not inside_box:
            print('***WARNING: Coordinates outside bounding box.***')
            #raise RuntimeWarning('Coordinates outside bounding box.\n\t    (x0,y0,z0) = (%f,%f,%f) \n\t input (x,y,z) = (%f,%f,%f) '%(self.x0,self.y0,self.z0,x,y,z))

        return ix, iy, iz, x1, y1, z1, inside_box
    
    def get_coefs(self,b):
        return np.matmul(tricubicMat, b)

    def val(self, x, y, z):
        
        ix, iy, iz, x1, y1, z1, inside_box = self.coords_to_indices_and_floats(x, y, z)

        if not inside_box:
            return 0

        b = self.construct_b(ix,iy,iz)
        coefs = np.matmul(tricubicMat, b)

        res=0
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    res += coefs[i+4*j+16*k]*x1**i*y1**j*z1**k
        return res


    def ddx(self, x, y, z):

        ix, iy, iz, x1, y1, z1, inside_box = self.coords_to_indices_and_floats(x, y, z)
        
        if not inside_box:
            return 0

        b = self.construct_b(ix,iy,iz)
        coefs = np.matmul(tricubicMat, b)
        
        res=0
        for i in range(1,4):
            for j in range(4):
                for k in range(4):
                    res += i*coefs[i+4*j+16*k]*x1**(i-1)*y1**j*z1**k
        return res/self.dx


    def ddy(self, x, y, z):

        ix, iy, iz, x1, y1, z1, inside_box = self.coords_to_indices_and_floats(x, y, z)
        
        if not inside_box:
            return 0

        b = self.construct_b(ix,iy,iz)
        coefs = np.matmul(tricubicMat, b)

        res=0
        for i in range(4):
            for j in range(1,4):
                for k in range(4):
                    res += j*coefs[i+4*j+16*k]*x1**i*y1**(j-1)*z1**k
        return res/self.dy


    def ddz(self, x, y, z):

        ix, iy, iz, x1, y1, z1, inside_box = self.coords_to_indices_and_floats(x, y, z)
        
        if not inside_box:
            return 0

        b = self.construct_b(ix,iy,iz)
        coefs = np.matmul(tricubicMat, b)

        res=0
        for i in range(4):
            for j in range(4):
                for k in range(1,4):
                    res += k*coefs[i+4*j+16*k]*x1**i*y1**j*z1**(k-1)
        return res/self.dz


    def ddxdy(self, x, y, z):

        ix, iy, iz, x1, y1, z1, inside_box = self.coords_to_indices_and_floats(x, y, z)
        
        if not inside_box:
            return 0

        b = self.construct_b(ix,iy,iz)
        coefs = np.matmul(tricubicMat, b)
        
        res=0
        for i in range(1,4):
            for j in range(1,4):
                for k in range(4):
                    res += i*j*coefs[i+4*j+16*k]*x1**(i-1)*y1**(j-1)*z1**k
        return res/self.dx/self.dy


    def ddxdz(self, x, y, z):

        ix, iy, iz, x1, y1, z1, inside_box = self.coords_to_indices_and_floats(x, y, z)
        
        if not inside_box:
            return 0

        b = self.construct_b(ix,iy,iz)
        coefs = np.matmul(tricubicMat, b)
        
        res=0
        for i in range(1,4):
            for j in range(4):
                for k in range(1,4):
                    res += i*k*coefs[i+4*j+16*k]*x1**(i-1)*y1**j*z1**(k-1)
        return res/self.dx/self.dz


    def ddydz(self, x, y, z):

        ix, iy, iz, x1, y1, z1, inside_box = self.coords_to_indices_and_floats(x, y, z)
        
        if not inside_box:
            return 0

        b = self.construct_b(ix,iy,iz)
        coefs = np.matmul(tricubicMat, b)
        
        res=0
        for i in range(4):
            for j in range(1,4):
                for k in range(1,4):
                    res += j*k*coefs[i+4*j+16*k]*x1**i*y1**(j-1)*z1**(k-1)
        return res/self.dy/self.dz


    def ddxdydz(self, x, y, z):

        ix, iy, iz, x1, y1, z1, inside_box = self.coords_to_indices_and_floats(x, y, z)
        
        if not inside_box:
            return 0

        b = self.construct_b(ix,iy,iz)
        coefs = np.matmul(tricubicMat, b)
        
        res=0
        for i in range(1,4):
            for j in range(1,4):
                for k in range(1,4):
                    res += i*j*k*coefs[i+4*j+16*k]*x1**(i-1)*y1**(j-1)*z1**(k-1)
        return res/self.dx/self.dy/self.dz


    def ddx2(self, x, y, z):

        ix, iy, iz, x1, y1, z1, inside_box = self.coords_to_indices_and_floats(x, y, z)
        
        if not inside_box:
            return 0

        b = self.construct_b(ix,iy,iz)
        coefs = np.matmul(tricubicMat, b)
        
        res=0
        for i in range(2,4):
            for j in range(4):
                for k in range(4):
                    res += (i-1)*i*coefs[i+4*j+16*k]*x1**(i-2)*y1**j*z1**k
        return res/self.dx/self.dx


    def ddy2(self, x, y, z):

        ix, iy, iz, x1, y1, z1, inside_box = self.coords_to_indices_and_floats(x, y, z)
        
        if not inside_box:
            return 0

        b = self.construct_b(ix,iy,iz)
        coefs = np.matmul(tricubicMat, b)

        res=0
        for i in range(4):
            for j in range(2,4):
                for k in range(4):
                    res += (j-1)*j*coefs[i+4*j+16*k]*x1**i*y1**(j-2)*z1**k
        return res/self.dy/self.dy


    def ddz2(self, x, y, z):

        ix, iy, iz, x1, y1, z1, inside_box = self.coords_to_indices_and_floats(x, y, z)
        
        if not inside_box:
            return 0

        b = self.construct_b(ix,iy,iz)
        coefs = np.matmul(tricubicMat, b)

        res=0
        for i in range(4):
            for j in range(4):
                for k in range(2,4):
                    res += (k-1)*k*coefs[i+4*j+16*k]*x1**i*y1**j*z1**(k-2)
        return res/self.dz/self.dz

