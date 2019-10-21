import sys
import numpy as np
if sys.version_info[0] < 3:
    from .Tricubic2_c import *
else:
    from .Tricubic_c import *


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
            if discardx != 1 or discardy != 1 or discardz != 1:
                raise Exception('Discards should be equal to 1 when using Exact derivatives method.')
            self.imethod = 2
            self.A = A
            self.A = np.ascontiguousarray(self.A)
            
            self.ix_bound_up = self.A.shape[0] - 2  #    a -1 because counting starts from 0,
            self.iy_bound_up = self.A.shape[1] - 2  #    
            self.iz_bound_up = self.A.shape[2] - 2  #    and a last -1 because the bound corresponds to the bound 
                                                    #    for the lower index inclusive

            self.ix_bound_low = 0
            self.iy_bound_low = 0
            self.iz_bound_low = 0
        elif method=='FD':
            #print('Using finite difference approximation for derivatives.')
            if len(A.shape) != 3:
                raise Exception('Input array should be 3-dimensional when using finite differences method. It\'s not.')
            self.imethod = 1
            self.A = A
            if self.discardx > 1:
                self.A = self.A[self.discardx-1:-(self.discardx-1),:,:]
            if self.discardy > 1:
                self.A = self.A[:,self.discardy-1:-(self.discardy-1),:]
            if self.discardz > 1:
                self.A = self.A[:,:,self.discardz-1:-(self.discardz-1)]
            self.A = np.ascontiguousarray(self.A)
            
            self.ix_bound_up = self.A.shape[0] - 3  #    a -1 because counting starts from 0,
            self.iy_bound_up = self.A.shape[1] - 3  #    another -1 because one more point is needed for finite differences,
            self.iz_bound_up = self.A.shape[2] - 3  #    and a last -1 because the bound corresponds to the bound 
                                                    #    for the lower index inclusive

            self.ix_bound_low = 1
            self.iy_bound_low = 1
            self.iz_bound_low = 1
        elif method=='Exact-Mirror2':
            #print('Using exact derivatives.')
            if len(A.shape) != 4:
                raise Exception('Input array should be 4-dimensional when using Exact-Mirror2 method. It\'s not.')
            if discardx != 1 or discardy != 1 or discardz != 1:
                raise Exception('Discards should be equal to 1 when using Exact-Mirror2 method.')
            self.imethod = 3
            self.A = A
            self.A = np.ascontiguousarray(self.A)
            
            self.ix_bound_up = self.A.shape[0] - 2  #    a -1 because counting starts from 0,
            self.iy_bound_up = self.A.shape[1] - 2  #    
            self.iz_bound_up = self.A.shape[2] - 2  #    and a last -1 because the bound corresponds to the bound 
                                                    #    for the lower index inclusive

            self.ix_bound_low = 0
            self.iy_bound_low = 0
            self.iz_bound_low = 0
        else:
            raise ValueError('Invalid method: %s'%method)



        if self.ix_bound_up < self.ix_bound_low:
            raise Exception('Interpolating array is too small along the first dimension (after discards).')
        if self.iy_bound_up < self.iy_bound_low:
            raise Exception('Interpolating array is too small along the second dimension (after discards).')
        if self.iz_bound_up < self.iz_bound_low:
            raise Exception('Interpolating array is too small along the third dimension (after discards).')

    def construct_b(self, x, y, z):
        return tricubic_py_get_b(self.A, x, y, z, self.x0, self.y0, self.z0, self.dx, self.dy, self.dz, self.ix_bound_low, self.ix_bound_up, self.iy_bound_low, self.iy_bound_up, self.iz_bound_low, self.iz_bound_up, self.imethod)

    def get_coefs(self, b):
        return tricubic_py_get_coefs(b)


    def coords_to_indices(self, x, y, z):
        return tricubic_py_coords_to_indices(self.A, x, y, z, self.x0, self.y0, self.z0, self.dx, self.dy, self.dz, self.ix_bound_low, self.ix_bound_up, self.iy_bound_low, self.iy_bound_up, self.iz_bound_low, self.iz_bound_up, self.imethod)
    
    def is_inside_box(self, x, y, z):
        return tricubic_py_is_inside_box(self.A, x, y, z, self.x0, self.y0, self.z0, self.dx, self.dy, self.dz, self.ix_bound_low, self.ix_bound_up, self.iy_bound_low, self.iy_bound_up, self.iz_bound_low, self.iz_bound_up, self.imethod)

    def val(self, x, y, z):
        return tricubic_get_val(self.A, x, y, z, self.x0, self.y0, self.z0, self.dx, self.dy, self.dz, self.ix_bound_low, self.ix_bound_up, self.iy_bound_low, self.iy_bound_up, self.iz_bound_low, self.iz_bound_up, self.imethod)

    def ddx(self, x, y, z):
        return tricubic_get_ddx(self.A, x, y, z, self.x0, self.y0, self.z0, self.dx, self.dy, self.dz, self.ix_bound_low, self.ix_bound_up, self.iy_bound_low, self.iy_bound_up, self.iz_bound_low, self.iz_bound_up, self.imethod)

    def ddy(self, x, y, z):
        return tricubic_get_ddy(self.A, x, y, z, self.x0, self.y0, self.z0, self.dx, self.dy, self.dz, self.ix_bound_low, self.ix_bound_up, self.iy_bound_low, self.iy_bound_up, self.iz_bound_low, self.iz_bound_up, self.imethod)

    def ddz(self, x, y, z):
        return tricubic_get_ddz(self.A, x, y, z, self.x0, self.y0, self.z0, self.dx, self.dy, self.dz, self.ix_bound_low, self.ix_bound_up, self.iy_bound_low, self.iy_bound_up, self.iz_bound_low, self.iz_bound_up, self.imethod)

    def ddxdy(self, x, y, z):
        return tricubic_get_ddxdy(self.A, x, y, z, self.x0, self.y0, self.z0, self.dx, self.dy, self.dz, self.ix_bound_low, self.ix_bound_up, self.iy_bound_low, self.iy_bound_up, self.iz_bound_low, self.iz_bound_up, self.imethod)

    def ddxdz(self, x, y, z):
        return tricubic_get_ddxdz(self.A, x, y, z, self.x0, self.y0, self.z0, self.dx, self.dy, self.dz, self.ix_bound_low, self.ix_bound_up, self.iy_bound_low, self.iy_bound_up, self.iz_bound_low, self.iz_bound_up, self.imethod)

    def ddydz(self, x, y, z):
        return tricubic_get_ddydz(self.A, x, y, z, self.x0, self.y0, self.z0, self.dx, self.dy, self.dz, self.ix_bound_low, self.ix_bound_up, self.iy_bound_low, self.iy_bound_up, self.iz_bound_low, self.iz_bound_up, self.imethod)

    def ddxdydz(self, x, y, z):
        return tricubic_get_ddxdydz(self.A, x, y, z, self.x0, self.y0, self.z0, self.dx, self.dy, self.dz, self.ix_bound_low, self.ix_bound_up, self.iy_bound_low, self.iy_bound_up, self.iz_bound_low, self.iz_bound_up, self.imethod)

    def ddx2(self, x, y, z):
        return tricubic_get_ddx2(self.A, x, y, z, self.x0, self.y0, self.z0, self.dx, self.dy, self.dz, self.ix_bound_low, self.ix_bound_up, self.iy_bound_low, self.iy_bound_up, self.iz_bound_low, self.iz_bound_up, self.imethod)

    def ddy2(self, x, y, z):
        return tricubic_get_ddy2(self.A, x, y, z, self.x0, self.y0, self.z0, self.dx, self.dy, self.dz, self.ix_bound_low, self.ix_bound_up, self.iy_bound_low, self.iy_bound_up, self.iz_bound_low, self.iz_bound_up, self.imethod)

    def ddz2(self, x, y, z):
        return tricubic_get_ddz2(self.A, x, y, z, self.x0, self.y0, self.z0, self.dx, self.dy, self.dz, self.ix_bound_low, self.ix_bound_up, self.iy_bound_low, self.iy_bound_up, self.iz_bound_low, self.iz_bound_up, self.imethod)

    def kick(self, x, y, z):
        return tricubic_get_kick(self.A, x, y, z, self.x0, self.y0, self.z0, self.dx, self.dy, self.dz, self.ix_bound_low, self.ix_bound_up, self.iy_bound_low, self.iy_bound_up, self.iz_bound_low, self.iz_bound_up, self.imethod)
