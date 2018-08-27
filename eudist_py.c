#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"
// https://docs.scipy.org/doc/numpy/reference/c-api.array.html

#include "eudist.c"

/*wrappted eudist*/
static PyObject* eudist_func(PyObject* self, PyObject* args)
{

  import_array();
  PyArrayObject *pyA, *pyB;

  /*  parse the input, from python float to c double */
  if (!PyArg_ParseTuple(args, "OO", &pyA, &pyB))
    return NULL;

  size_t NA = (size_t) PyArray_Size((PyObject *) pyA);
  size_t NB = (size_t) PyArray_Size((PyObject *) pyB);

  if(NA==0 || NB==0)
  {
    PyErr_SetString(PyExc_RuntimeError, "Arguments have to be numpy arrays.");
    return NULL;
  }

 if(!(NB==3))
  {
    PyErr_SetString(PyExc_RuntimeError, "Voxel size has to be given by three numbers\n");
    return NULL;
  }

  double * A = (double*) PyArray_DATA(pyA);
  double * B = (double*) PyArray_DATA(pyB);

  

  /*  construct the output */
  return Py_BuildValue("f", 3);
}

/*  define functions in module */
static PyMethodDef eudistMethods[] =
{
  {"eudist", eudist_func, METH_VARARGS, "euclidean distance function"},
  {NULL, NULL, 0, NULL}
};

/* Python version 3*/
static struct PyModuleDef cModPyDem =
{
  PyModuleDef_HEAD_INIT,
  "eudist_module", "Some documentation",
  -1,
  eudistMethods
};

  PyMODINIT_FUNC
PyInit_eudist(void)
{
  return PyModule_Create(&cModPyDem);
}
