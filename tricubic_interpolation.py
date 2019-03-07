import numpy as np
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
        #self.x0 += (self.discard0-1)*self.dx
        #self.y0 += (self.discard1-1)*self.dy
        #self.z0 += (self.discard2-1)*self.dz
        self.A = A

        #if A.shape[0]- 2*self.discard0 < 2:
        #    raise Exception('n0 < 2: Interpolating array is too small along the first dimension (after discards).')
        #if A.shape[1]- 2*self.discard1 < 2:
        #    raise Exception('n1 < 2: Interpolating array is too small along the second dimension (after discards).')
        #if A.shape[2]- 2*self.discard2 < 2:
        #    raise Exception('n2 < 2: Interpolating array is too small along the third dimension (after discards).')

        #v = self.finite_differences(A, self.discard0, self.discard1, self.discard2)
        #self.coefs = self.tricubic_coefs(v)

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

    def val(self, x, y, z):

        fx = (x - self.x0)/self.dx
        fy = (y - self.y0)/self.dy
        fz = (z - self.z0)/self.dz

        ix = int(fx)
        iy = int(fy)
        iz = int(fz)
        
        if ix < 1 or ix > self.A.shape[0]:
            return 0
 #           raise Exception('Position is outside bounding box (first dimension): %d'%ix)

        if iy < 1 or iy > self.A.shape[1]:
            return 0
 #           raise Exception('Position is outside bounding box (second dimension): %d'%iy)

        if iz < 1 or iz > self.A.shape[2]:
            return 0
 #           raise Exception('Position is outside bounding box (third dimension): %d'%iz)

        b = self.finite_diff(ix,iy,iz)
        coefs = np.matmul(tricubicMat, b)

        x1 = fx - ix
        y1 = fy - iy
        z1 = fz - iz
        
        res=0
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    res += coefs[i+4*j+16*k]*x1**i*y1**j*z1**k
        return res


    def ddx(self, x, y, z):

        fx = (x - self.x0)/self.dx
        fy = (y - self.y0)/self.dy
        fz = (z - self.z0)/self.dz

        ix = int(fx)
        iy = int(fy)
        iz = int(fz)
        
        if ix < 1 or ix > self.A.shape[0]:
            return 0
 #           raise Exception('Position is outside bounding box (first dimension): %d'%ix)

        if iy < 1 or iy > self.A.shape[1]:
            return 0
 #           raise Exception('Position is outside bounding box (second dimension): %d'%iy)

        if iz < 1 or iz > self.A.shape[2]:
            return 0
 #           raise Exception('Position is outside bounding box (third dimension): %d'%iz)

        b = self.finite_diff(ix,iy,iz)
        coefs = np.matmul(tricubicMat, b)

        x1 = fx - ix
        y1 = fy - iy
        z1 = fz - iz
        
        res=0
        for i in range(1,4):
            for j in range(4):
                for k in range(4):
                    res += i*coefs[i+4*j+16*k]*x1**(i-1)*y1**j*z1**k
        return res/self.dx


    def ddy(self, x, y, z):

        fx = (x - self.x0)/self.dx
        fy = (y - self.y0)/self.dy
        fz = (z - self.z0)/self.dz

        ix = int(fx)
        iy = int(fy)
        iz = int(fz)
        
        if ix < 1 or ix > self.A.shape[0]:
            return 0
 #           raise Exception('Position is outside bounding box (first dimension): %d'%ix)

        if iy < 1 or iy > self.A.shape[1]:
            return 0
 #           raise Exception('Position is outside bounding box (second dimension): %d'%iy)

        if iz < 1 or iz > self.A.shape[2]:
            return 0
 #           raise Exception('Position is outside bounding box (third dimension): %d'%iz)

        b = self.finite_diff(ix,iy,iz)
        coefs = np.matmul(tricubicMat, b)

        x1 = fx - ix
        y1 = fy - iy
        z1 = fz - iz
        
        res=0
        for i in range(4):
            for j in range(1,4):
                for k in range(4):
                    res += j*coefs[i+4*j+16*k]*x1**i*y1**(j-1)*z1**k
        return res/self.dy


    def ddz(self, x, y, z):

        fx = (x - self.x0)/self.dx
        fy = (y - self.y0)/self.dy
        fz = (z - self.z0)/self.dz

        ix = int(fx)
        iy = int(fy)
        iz = int(fz)
        
        if ix < 1 or ix > self.A.shape[0]:
            return 0
 #           raise Exception('Position is outside bounding box (first dimension): %d'%ix)

        if iy < 1 or iy > self.A.shape[1]:
            return 0
 #           raise Exception('Position is outside bounding box (second dimension): %d'%iy)

        if iz < 1 or iz > self.A.shape[2]:
            return 0
 #           raise Exception('Position is outside bounding box (third dimension): %d'%iz)

        b = self.finite_diff(ix,iy,iz)
        coefs = np.matmul(tricubicMat, b)

        x1 = fx - ix
        y1 = fy - iy
        z1 = fz - iz
        
        res=0
        for i in range(4):
            for j in range(4):
                for k in range(1,4):
                    res += k*coefs[i+4*j+16*k]*x1**i*y1**j*z1**(k-1)
        return res/self.dz

