#include "pywrap.h"

static PyObject* get_val(PyObject* self, PyObject* args)
{

    int ix, iy, iz;
    double xn, yn, zn;
    double dx, dy, dz;
    int method;
    PyObject *arg1 = NULL;
    if(!PyArg_ParseTuple(args, "O!iiiddddddi", &PyArray_Type, &arg1, &ix, &iy, &iz, &xn, &yn, &zn, &dx, &dy, &dz, &method))
        return NULL;

    PyObject *py_A = NULL;
    py_A = PyArray_FROM_OTF(arg1, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);

    if(py_A == NULL)
        return NULL;

    double* c_A = (double*)PyArray_DATA((PyArrayObject*)py_A);
    // int nd = (int)PyArray_NDIM((PyArrayObject*)py_A);
       
    npy_intp* shape = PyArray_DIMS((PyArrayObject*)py_A);

    int shape1 = (int)shape[1];
    int shape2 = (int)shape[2];

    double* b = NULL;
    if(method == 1)
        b = finite_diff(shape1, shape2, c_A, ix, iy, iz);
    else if(method == 2)
    {
        int shape3 = (int)shape[3];
        b = exact_diff(shape1, shape2, shape3, c_A, ix, iy, iz, dx, dy, dz);
    }
    else
    {
        printf("Method not recognized");
        Py_INCREF(Py_None);
        return NULL;
    }
    
    double* coefs = get_coefs(b);
    
    double xni[4], ynj[4], znk[4];
    xyz_powers(xn, yn, zn, xni, ynj, znk);

    double val = eval(coefs, xni, ynj, znk);

    free(coefs);
    free(b);
    Py_DECREF(py_A);

    return Py_BuildValue("d", val);
}

static PyObject* get_ddx(PyObject* self, PyObject* args)
{

    int ix, iy, iz;
    double xn, yn, zn;
    double dx, dy, dz;
    int method;
    PyObject *arg1 = NULL;
    if(!PyArg_ParseTuple(args, "O!iiiddddddi", &PyArray_Type, &arg1, &ix, &iy, &iz, &xn, &yn, &zn, &dx, &dy, &dz, &method))
        return NULL;

    PyObject *py_A = NULL;
    py_A = PyArray_FROM_OTF(arg1, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);

    if(py_A == NULL)
        return NULL;

    double* c_A = (double*)PyArray_DATA((PyArrayObject*)py_A);
    // int nd = (int)PyArray_NDIM((PyArrayObject*)py_A);
       
    npy_intp* shape = PyArray_DIMS((PyArrayObject*)py_A);

    int shape1 = (int)shape[1];
    int shape2 = (int)shape[2];

    double* b = NULL;
    if(method == 1)
        b = finite_diff(shape1, shape2, c_A, ix, iy, iz);
    else if(method == 2)
    {
        int shape3 = (int)shape[3];
        b = exact_diff(shape1, shape2, shape3, c_A, ix, iy, iz, dx, dy, dz);
    }
    else
    {
        printf("Method not recognized");
        Py_INCREF(Py_None);
        return NULL;
    }
    
    double* coefs = get_coefs(b);
    
    double xni[4], ynj[4], znk[4];
    xyz_powers(xn, yn, zn, xni, ynj, znk);

    double ddx_value = ddx(coefs, xni, ynj, znk, dx);

    free(coefs);
    free(b);
    Py_DECREF(py_A);

    return Py_BuildValue("d", ddx_value);
}
static PyMethodDef TricubicMethods[] = 
{
    {"get_val", get_val, METH_VARARGS, "returns interpolated value"},
    {"get_ddx", get_ddx, METH_VARARGS, "returns interpolated first derivative with respect to x"},
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
