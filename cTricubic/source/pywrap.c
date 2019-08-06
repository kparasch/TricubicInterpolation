#include "pywrap.h"

#define TRICUBIC_PROTOTYPE_GET_MACRO                                           \
                                                                               \
    double dx, dy, dz;                                                         \
    double x, y, z;                                                            \
    double x0, y0, z0;                                                         \
                                                                               \
    int ix_bound_low, ix_bound_up;                                             \
    int iy_bound_low, iy_bound_up;                                             \
    int iz_bound_low, iz_bound_up;                                             \
                                                                               \
    int method;                                                                \
    PyObject *arg1 = NULL;                                                     \
    if(!PyArg_ParseTuple(args, "O!dddddddddiiiiiii", &PyArray_Type, &arg1, &x, \
                         &y, &z, &x0, &y0, &z0, &dx, &dy, &dz, &ix_bound_low,  \
                         &ix_bound_up, &iy_bound_low, &iy_bound_up,            \
                         &iz_bound_low, &iz_bound_up, &method                  \
                        )                                                      \
      )                                                                        \
        return NULL;                                                           \
                                                                               \
    PyObject *py_A = NULL;                                                     \
    py_A = PyArray_FROM_OTF(arg1, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);            \
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
    int ix, iy, iz;                                                            \
    double xn, yn, zn;                                                         \
    int is_inside = tricubic_coords_to_indices_and_floats(x, y, z, x0, y0, z0, \
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
    else if(method == 2)                                                       \
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
    Py_DECREF(py_A);

    return Py_BuildValue("d", val);
}

static PyObject* tricubic_get_ddx(PyObject* self, PyObject* args)
{
    TRICUBIC_PROTOTYPE_GET_MACRO

    double val = tricubic_ddx(coefs, xni, ynj, znk, dx);

    free(coefs);
    free(b);
    Py_DECREF(py_A);

    return Py_BuildValue("d", val);
}

static PyObject* tricubic_get_ddy(PyObject* self, PyObject* args)
{
    TRICUBIC_PROTOTYPE_GET_MACRO

    double val = tricubic_ddy(coefs, xni, ynj, znk, dy);

    free(coefs);
    free(b);
    Py_DECREF(py_A);

    return Py_BuildValue("d", val);
}

static PyObject* tricubic_get_ddz(PyObject* self, PyObject* args)
{
    TRICUBIC_PROTOTYPE_GET_MACRO

    double val = tricubic_ddz(coefs, xni, ynj, znk, dz);

    free(coefs);
    free(b);
    Py_DECREF(py_A);

    return Py_BuildValue("d", val);
}

static PyObject* tricubic_get_ddxdy(PyObject* self, PyObject* args)
{
    TRICUBIC_PROTOTYPE_GET_MACRO

    double val = tricubic_ddxdy(coefs, xni, ynj, znk, dx, dy);

    free(coefs);
    free(b);
    Py_DECREF(py_A);

    return Py_BuildValue("d", val);
}

static PyObject* tricubic_get_ddxdz(PyObject* self, PyObject* args)
{
    TRICUBIC_PROTOTYPE_GET_MACRO

    double val = tricubic_ddxdz(coefs, xni, ynj, znk, dx, dz);

    free(coefs);
    free(b);
    Py_DECREF(py_A);

    return Py_BuildValue("d", val);
}

static PyObject* tricubic_get_ddydz(PyObject* self, PyObject* args)
{
    TRICUBIC_PROTOTYPE_GET_MACRO

    double val = tricubic_ddydz(coefs, xni, ynj, znk, dy, dz);

    free(coefs);
    free(b);
    Py_DECREF(py_A);

    return Py_BuildValue("d", val);
}

static PyObject* tricubic_get_ddxdydz(PyObject* self, PyObject* args)
{
    TRICUBIC_PROTOTYPE_GET_MACRO

    double val = tricubic_ddxdydz(coefs, xni, ynj, znk, dx, dy, dz);

    free(coefs);
    free(b);
    Py_DECREF(py_A);

    return Py_BuildValue("d", val);
}

static PyObject* tricubic_py_coords_to_indices(PyObject* self, PyObject* args)
{
    TRICUBIC_PROTOTYPE_GET_MACRO

    free(coefs);
    free(b);
    Py_DECREF(py_A);

    return Py_BuildValue("iii", ix, iy, iz);
}

static PyObject* tricubic_py_get_b(PyObject* self, PyObject* args)
{
    TRICUBIC_PROTOTYPE_GET_MACRO
    
    npy_intp dims[1];
    dims[0] = 64;

    PyObject* b_numpy = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT64, (void*) b);
    PyArray_ENABLEFLAGS((PyArrayObject *) b_numpy, NPY_ARRAY_OWNDATA);

    free(coefs);
    Py_DECREF(py_A);

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
    {"tricubic_py_coords_to_indices", tricubic_py_coords_to_indices, METH_VARARGS, "returns indices."},
    {"tricubic_py_get_b", tricubic_py_get_b, METH_VARARGS, "returns b."},
    {"tricubic_py_get_coefs", tricubic_py_get_coefs, METH_VARARGS, "returns coefs."},
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
