#include "pywrap.h"

static PyObject* get_val(PyObject* self, PyObject* args)
{

    int ix, iy, iz;
    double xn, yn, zn;
    double dx, dy, dz;
    int method;
    PyObject *arg1 = NULL;
    if(!PyArg_ParseTuple(args, "0!iiiddddddi", &PyArray_Type, &arg1, &ix, &iy, &iz, &xn, &yn, &zn, &dx, &dy, &dz, &method))
        return NULL;

    PyObject *py_A = NULL;
    py_A = PyArray_FROM_OTF(arg1, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);

    if(py_A == NULL)
        return NULL;

    double* c_A = (double*)PyArray_DATA((PyArrayObject*)py_A);
    int nd = (int)PyArray_NDIM((PyArrayObject*)py_A);
       
    int* shape = (int*)PyArray_DIMS((PyArrayObject*)py_A);

    if(method == 1)
        double* b = finite_diff(shape[1], shape[2], c_A, ix, iy, iz);
    else if(method == 2)
        double* b = exact_diff(shape[1], shape[2], shape[3], c_A, ix, iy, iz, dx, dy, dz);
    else
    {
        printf("Method not recognized");
        Py_INCREF(Py_None);
        return NULL;
    }
    
    double* coefs = get_coefs(b);
    
    double val = eval(coefs, xn, yn, zn);

    free(coefs);
    free(b);
    Py_DECREF(py_A);

    return Py_BuildValue("d", val);
}

static PyMethodDef TricubicMethods[] = 
{
    {"get_val", get_val, METH_VARARGS, "returns interpolated value"},
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
