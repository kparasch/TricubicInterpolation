#include "pywrap.h"
#include <stdio.h>

#define TRICUBIC_PROTOTYPE_GET_MACRO                                           \
    double dx, dy, dz;                                                         \
    double x, y, z;                                                            \
    double x0, y0, z0;                                                         \
                                                                               \
    int ix_bound_low, ix_bound_up;                                             \
    int iy_bound_low, iy_bound_up;                                             \
    int iz_bound_low, iz_bound_up;                                             \
                                                                               \
    int method;                                                                \
    PyArrayObject *py_A = NULL;                                                \
    if(!PyArg_ParseTuple(args, "Odddddddddiiiiiii", &py_A, &x,                 \
                         &y, &z, &x0, &y0, &z0, &dx, &dy, &dz, &ix_bound_low,  \
                         &ix_bound_up, &iy_bound_low, &iy_bound_up,            \
                         &iz_bound_low, &iz_bound_up, &method                  \
                        )                                                      \
      )                                                                        \
        return NULL;                                                           \
                                                                               \
    if(py_A == NULL)                                                           \
        return NULL;                                                           \
                                                                               \
    double* c_A = (double*)PyArray_DATA((PyArrayObject*)py_A);                 \
                                                                               \
    npy_intp* shape = PyArray_DIMS((PyArrayObject*)py_A);                      \
    int shape1 = (int)shape[1];                                                \
    int shape2 = (int)shape[2];                                                \
                                                                               \
    int sx = 1;                                                                \
    int sy = 1;                                                                \
                                                                               \
    if(method == 3)                                                            \
    {                                                                          \
        if(x < 0)                                                              \
            sx = -1;                                                           \
        if(y < 0)                                                              \
            sy = -1;                                                           \
    }                                                                          \
                                                                               \
    int ix, iy, iz;                                                            \
    double xn, yn, zn;                                                         \
    int is_inside = tricubic_coords_to_indices_and_floats(sx*x, sy*y, z,       \
                                                          x0, y0, z0,          \
                                                          dx, dy, dz, &ix, &iy,\
                                                          &iz, &xn, &yn, &zn,  \
                                                          ix_bound_low,        \
                                                          ix_bound_up,         \
                                                          iy_bound_low,        \
                                                          iy_bound_up,         \
                                                          iz_bound_low,        \
                                                          iz_bound_up          \
                                                         );                    \
                                                                               \
    if(!is_inside)                                                             \
    {                                                                          \
        printf("***WARNING: Coordinates outside bounding box.***\n");          \
        return Py_BuildValue("d", 0.);                                         \
    }                                                                          \
                                                                               \
    double* b = NULL;                                                          \
    if(method == 1)                                                            \
        b = tricubic_finite_diff(shape1, shape2, c_A, ix, iy, iz);             \
    else if(method == 2 || method == 3)                                                       \
    {                                                                          \
        int shape3 = (int)shape[3];                                            \
        b = tricubic_exact_diff(shape1, shape2, shape3, c_A, ix, iy, iz, dx,   \
                                dy, dz                                         \
                               );                                              \
    }                                                                          \
    else                                                                       \
    {                                                                          \
        printf("Method not recognized");                                       \
        Py_INCREF(Py_None);                                                    \
        return NULL;                                                           \
    }                                                                          \
                                                                               \
    double* coefs = tricubic_get_coefs(b);                                     \
                                                                               \
    double xni[4], ynj[4], znk[4];                                             \
    tricubic_xyz_powers(xn, yn, zn, xni, ynj, znk);                            
                                            
static PyObject* tricubic_get_val(PyObject* self, PyObject* args)
{
    TRICUBIC_PROTOTYPE_GET_MACRO

    double val = tricubic_eval(coefs, xni, ynj, znk);

    free(coefs);
    free(b);

    return Py_BuildValue("d", val);
}

static PyObject* tricubic_get_ddx(PyObject* self, PyObject* args)
{
    TRICUBIC_PROTOTYPE_GET_MACRO

    double val = sx*tricubic_ddx(coefs, xni, ynj, znk, dx);

    free(coefs);
    free(b);

    return Py_BuildValue("d", val);
}

static PyObject* tricubic_get_ddy(PyObject* self, PyObject* args)
{
    TRICUBIC_PROTOTYPE_GET_MACRO

    double val = sy*tricubic_ddy(coefs, xni, ynj, znk, dy);

    free(coefs);
    free(b);

    return Py_BuildValue("d", val);
}

static PyObject* tricubic_get_ddz(PyObject* self, PyObject* args)
{
    TRICUBIC_PROTOTYPE_GET_MACRO

    double val = tricubic_ddz(coefs, xni, ynj, znk, dz);

    free(coefs);
    free(b);

    return Py_BuildValue("d", val);
}

static PyObject* tricubic_get_ddxdy(PyObject* self, PyObject* args)
{
    TRICUBIC_PROTOTYPE_GET_MACRO

    double val = (sx*sy)*tricubic_ddxdy(coefs, xni, ynj, znk, dx, dy);

    free(coefs);
    free(b);

    return Py_BuildValue("d", val);
}

static PyObject* tricubic_get_ddxdz(PyObject* self, PyObject* args)
{
    TRICUBIC_PROTOTYPE_GET_MACRO

    double val = sx*tricubic_ddxdz(coefs, xni, ynj, znk, dx, dz);

    free(coefs);
    free(b);

    return Py_BuildValue("d", val);
}

static PyObject* tricubic_get_ddydz(PyObject* self, PyObject* args)
{
    TRICUBIC_PROTOTYPE_GET_MACRO

    double val = sy*tricubic_ddydz(coefs, xni, ynj, znk, dy, dz);

    free(coefs);
    free(b);

    return Py_BuildValue("d", val);
}

static PyObject* tricubic_get_ddxdydz(PyObject* self, PyObject* args)
{
    TRICUBIC_PROTOTYPE_GET_MACRO

    double val = (sx*sy)*tricubic_ddxdydz(coefs, xni, ynj, znk, dx, dy, dz);

    free(coefs);
    free(b);

    return Py_BuildValue("d", val);
}

static PyObject* tricubic_get_ddx2(PyObject* self, PyObject* args)
{
    TRICUBIC_PROTOTYPE_GET_MACRO

    double val = tricubic_ddx2(coefs, xni, ynj, znk, dx);

    free(coefs);
    free(b);

    return Py_BuildValue("d", val);
}

static PyObject* tricubic_get_ddy2(PyObject* self, PyObject* args)
{
    TRICUBIC_PROTOTYPE_GET_MACRO

    double val = tricubic_ddy2(coefs, xni, ynj, znk, dy);

    free(coefs);
    free(b);

    return Py_BuildValue("d", val);
}

static PyObject* tricubic_get_ddz2(PyObject* self, PyObject* args)
{
    TRICUBIC_PROTOTYPE_GET_MACRO

    double val = tricubic_ddz2(coefs, xni, ynj, znk, dz);

    free(coefs);
    free(b);

    return Py_BuildValue("d", val);
}

static PyObject* tricubic_get_kick(PyObject* self, PyObject* args)
{
    double dx, dy, dz;                                                         
    double x, y, z;                                                            
    double x0, y0, z0;                                                         
                                                                               
    int ix_bound_low, ix_bound_up;                                             
    int iy_bound_low, iy_bound_up;                                             
    int iz_bound_low, iz_bound_up;                                             
                                                                               
    int method;                                                                
    PyArrayObject *py_A = NULL;                                                
    if(!PyArg_ParseTuple(args, "Odddddddddiiiiiii", &py_A, &x,                 
                         &y, &z, &x0, &y0, &z0, &dx, &dy, &dz, &ix_bound_low,  
                         &ix_bound_up, &iy_bound_low, &iy_bound_up,            
                         &iz_bound_low, &iz_bound_up, &method                  
                        )                                                      
      )                                                                        
        return NULL;                                                           
                                                                               
    if(py_A == NULL)                                                           
        return NULL;                                                           
                                                                               
    double* c_A = (double*)PyArray_DATA((PyArrayObject*)py_A);                 
                                                                               
    npy_intp* shape = PyArray_DIMS((PyArrayObject*)py_A);                      
    int shape1 = (int)shape[1];                                                
    int shape2 = (int)shape[2];                                                
                                                                               
    int sx = 1;                                                                
    int sy = 1;                                                                
                                                                               
    if(method == 3)                                                            
    {                                                                          
        if(x < 0)                                                              
            sx = -1;                                                           
        if(y < 0)                                                              
            sy = -1;                                                           
    }                                                                          
                                                                               
    int ix, iy, iz;                                                            
    double xn, yn, zn;                                                         
    int is_inside = tricubic_coords_to_indices_and_floats(sx*x, sy*y, z,       
                                                          x0, y0, z0, 
                                                          dx, dy, dz, &ix, &iy,
                                                          &iz, &xn, &yn, &zn,  
                                                          ix_bound_low,        
                                                          ix_bound_up,         
                                                          iy_bound_low,        
                                                          iy_bound_up,         
                                                          iz_bound_low,        
                                                          iz_bound_up          
                                                         );                    
                                                                               
    if(!is_inside)                                                             
    {                                                                          
        printf("***WARNING: Coordinates outside bounding box.***\n");          
        return Py_BuildValue("ddd", 0., 0., 0.);                               
    }                                                                          
                                                                               
    double* b = NULL;                                                          
    if(method == 1)                                                            
        b = tricubic_finite_diff(shape1, shape2, c_A, ix, iy, iz);             
    else if(method == 2 || method == 3)                                                       
    {                                                                          
        int shape3 = (int)shape[3];                                            
        b = tricubic_exact_diff(shape1, shape2, shape3, c_A, ix, iy, iz, dx,   
                                dy, dz                                         
                               );                                              
    }                                                                          
    else                                                                       
    {                                                                          
        printf("Method not recognized");                                       
        Py_INCREF(Py_None);                                                    
        return NULL;                                                           
    }                                                                          
                                                                               
    double* coefs = tricubic_get_coefs(b);                                     
                                                                               
    double xni[4], ynj[4], znk[4];                                             
    tricubic_xyz_powers(xn, yn, zn, xni, ynj, znk);                            

    double xkick = -sx*tricubic_ddx(coefs, xni, ynj, znk, dx);
    double ykick = -sy*tricubic_ddy(coefs, xni, ynj, znk, dy);
    double zkick = -tricubic_ddz(coefs, xni, ynj, znk, dz);

    free(coefs);
    free(b);

    return Py_BuildValue("ddd", xkick, ykick, zkick);
}

static PyObject* tricubic_py_coords_to_indices(PyObject* self, PyObject* args)
{
    TRICUBIC_PROTOTYPE_GET_MACRO

    free(coefs);
    free(b);

    return Py_BuildValue("iii", ix, iy, iz);
}

static PyObject* tricubic_py_get_b(PyObject* self, PyObject* args)
{
    double dx, dy, dz;                                                            

    int ix, iy, iz;                                                            
    int ix_bound_low, ix_bound_up;                                             
    int iy_bound_low, iy_bound_up;                                             
    int iz_bound_low, iz_bound_up;                                             
                                                                               
    int method;                                                                
    PyArrayObject *py_A = NULL;                                                
    if(!PyArg_ParseTuple(args, "Oiiidddiiiiiii", &py_A, &ix,                 
                         &iy, &iz, &dx, &dy, &dz, &ix_bound_low,  
                         &ix_bound_up, &iy_bound_low, &iy_bound_up,            
                         &iz_bound_low, &iz_bound_up, &method                  
                        )                                                      
      )                                                                        
        return NULL;                                                           
                                                                               
    if(py_A == NULL)                                                           
        return NULL;                                                           
                                                                               
    double* c_A = (double*)PyArray_DATA((PyArrayObject*)py_A);                 
                                                                               
    npy_intp* shape = PyArray_DIMS((PyArrayObject*)py_A);                      
    int shape1 = (int)shape[1];                                                
    int shape2 = (int)shape[2];                                                
                                                                               
    int is_inside = 1;
    if( ix < ix_bound_low || ix > ix_bound_up )
        is_inside = 0;
    else if( iy < iy_bound_low || iy > iy_bound_up )
        is_inside = 0;
    else if( iz < iz_bound_low || iz > iz_bound_up )
        is_inside = 0;

    if(!is_inside)                                                             
    {                                                                          
        printf("***WARNING: Coordinates outside bounding box.***\n");          
        return Py_BuildValue("d", 0.);                                         
    }                                                                          
                                                                               
    double* b = NULL;                                                          
    if(method == 1)                                                            
        b = tricubic_finite_diff(shape1, shape2, c_A, ix, iy, iz);             
    else if(method == 2)                                                       
    {                                                                          
        int shape3 = (int)shape[3];                                            
        b = tricubic_exact_diff(shape1, shape2, shape3, c_A, ix, iy, iz, dx,   
                                dy, dz                                         
                               );                                              
    }                                                                          
    else                                                                       
    {                                                                          
        printf("Method not recognized");                                       
        Py_INCREF(Py_None);                                                    
        return NULL;                                                           
    }                                                                          
                                                                               
    npy_intp dims[1];
    dims[0] = 64;

    PyObject* b_numpy = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT64, (void*) b);
    PyArray_ENABLEFLAGS((PyArrayObject *) b_numpy, NPY_ARRAY_OWNDATA);

    return Py_BuildValue("N", b_numpy);
}

static PyObject* tricubic_py_get_coefs(PyObject* self, PyObject* args)
{
    
    PyObject *arg1 = NULL;
    if(!PyArg_ParseTuple(args, "O!", &PyArray_Type, &arg1))
         return NULL;

    PyObject *py_b = NULL;                                                     
    py_b = PyArray_FROM_OTF(arg1, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);            
                                                                               
    if(py_b == NULL)                                                           
        return NULL;                                                           
                                                                               
    double* c_b = (double*)PyArray_DATA((PyArrayObject*)py_b);                 
    double* coefs = tricubic_get_coefs(c_b);
    npy_intp dims[1];
    dims[0] = 64;

    PyObject* coefs_numpy = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT64, (void*) coefs);
    PyArray_ENABLEFLAGS((PyArrayObject *) coefs_numpy, NPY_ARRAY_OWNDATA);

    Py_DECREF(py_b);

    return Py_BuildValue("N", coefs_numpy);
}

static PyObject* tricubic_py_is_inside_box(PyObject* self, PyObject* args)
{
    double dx, dy, dz;                                                         
    double x, y, z;                                                            
    double x0, y0, z0;                                                         
                                                                               
    int ix_bound_low, ix_bound_up;                                             
    int iy_bound_low, iy_bound_up;                                             
    int iz_bound_low, iz_bound_up;                                             
                                                                               
    int method;                                                                
    PyArrayObject *py_A = NULL;                                                
    if(!PyArg_ParseTuple(args, "Odddddddddiiiiiii", &py_A, &x,                 
                         &y, &z, &x0, &y0, &z0, &dx, &dy, &dz, &ix_bound_low,  
                         &ix_bound_up, &iy_bound_low, &iy_bound_up,            
                         &iz_bound_low, &iz_bound_up, &method                  
                        )                                                      
      )                                                                        
        return NULL;                                                           
                                                                               
    if(py_A == NULL)                                                           
        return NULL;                                                           
                                                                               
    // double* c_A = (double*)PyArray_DATA((PyArrayObject*)py_A);                 
    //                                                                            
    // npy_intp* shape = PyArray_DIMS((PyArrayObject*)py_A);                      
    // int shape1 = (int)shape[1];                                                
    // int shape2 = (int)shape[2];                                                
                                                                               
    int ix, iy, iz;                                                            
    double xn, yn, zn;                                                         
    int is_inside = tricubic_coords_to_indices_and_floats(x, y, z, x0, y0, z0, 
                                                          dx, dy, dz, &ix, &iy,
                                                          &iz, &xn, &yn, &zn,  
                                                          ix_bound_low,        
                                                          ix_bound_up,         
                                                          iy_bound_low,        
                                                          iy_bound_up,         
                                                          iz_bound_low,        
                                                          iz_bound_up          
                                                         );                    
                                                                               
    return Py_BuildValue("i", is_inside);
}

static PyMethodDef TricubicMethods[] = 
{
    {"tricubic_get_val", tricubic_get_val, METH_VARARGS, "returns interpolated value"},
    {"tricubic_get_ddx", tricubic_get_ddx, METH_VARARGS, "returns interpolated first derivative with respect to x"},
    {"tricubic_get_ddy", tricubic_get_ddy, METH_VARARGS, "returns interpolated first derivative with respect to y"},
    {"tricubic_get_ddz", tricubic_get_ddz, METH_VARARGS, "returns interpolated first derivative with respect to z"},
    {"tricubic_get_ddxdy", tricubic_get_ddxdy, METH_VARARGS, "returns interpolated first derivative with respect to x and y"},
    {"tricubic_get_ddxdz", tricubic_get_ddxdz, METH_VARARGS, "returns interpolated first derivative with respect to x and z"},
    {"tricubic_get_ddydz", tricubic_get_ddydz, METH_VARARGS, "returns interpolated first derivative with respect to y and z"},
    {"tricubic_get_ddxdydz", tricubic_get_ddxdydz, METH_VARARGS, "returns interpolated first derivative with respect to x, y and z"},
    {"tricubic_get_ddx2", tricubic_get_ddx2, METH_VARARGS, "returns interpolated second derivative with respect to x"},
    {"tricubic_get_ddy2", tricubic_get_ddy2, METH_VARARGS, "returns interpolated second derivative with respect to y"},
    {"tricubic_get_ddz2", tricubic_get_ddz2, METH_VARARGS, "returns interpolated second derivative with respect to z"},
    {"tricubic_get_kick", tricubic_get_kick, METH_VARARGS, "returns interpolated kicks"},
    {"tricubic_py_coords_to_indices", tricubic_py_coords_to_indices, METH_VARARGS, "returns indices."},
    {"tricubic_py_get_b", tricubic_py_get_b, METH_VARARGS, "returns b."},
    {"tricubic_py_get_coefs", tricubic_py_get_coefs, METH_VARARGS, "returns coefs."},
    {"tricubic_py_is_inside_box", tricubic_py_is_inside_box, METH_VARARGS, "returns if position is inside bounding box."},
    {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef cModTricubic_c = 
{
    PyModuleDef_HEAD_INIT,
    "Tricubic_c", "Tricubic interpolation library.",
    -1,
    TricubicMethods
};

PyMODINIT_FUNC PyInit_Tricubic_c(void)
{
    import_array();
    return PyModule_Create(&cModTricubic_c);
}

#else

PyMODINIT_FUNC initTricubic2_c(void)
{
    (void) Py_InitModule("Tricubic2_c", TricubicMethods);
    import_array();
}

#endif
