#ifndef __TRICUBIC_PYWRAP_H__
#define __TRICUBIC_PYWRAP_H__

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERISON

#include <Python.h>
#include <numpy/arrayobject.h>

#include "evaluate.h"
#include "coefs.h"
#include "derivs.h"
#include "coords.h"

static PyObject* tricubic_get_val(PyObject* self, PyObject* args);
static PyObject* tricubic_get_ddx(PyObject* self, PyObject* args);

#endif
