#include <Python.h>
#include <eudist.c>

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

 if(!(NA==NB))
  {
    PyErr_SetString(PyExc_RuntimeError, "Arrays have to have the same number of elements\n");
    return NULL;
  }

  double * A = (double*) PyArray_DATA(pyA);
  double * B = (double*) PyArray_DATA(pyB);

   size_t N = NA;

  if(N%2 ==1)
  {
    PyErr_SetString(PyExc_RuntimeError, "The array has to have an even number of elements\n");
    return NULL;
  }

  N = N/2;

  // printf("N: %zu\n", N);

  /* if the above function returns -1, an appropriate Python exception will
   * have been set, and the function simply returns NULL
   */

  /* call cos from libm */
  double tau = eudist2(A, B, N);

  /*  construct the output from cos, from c double to python float */
  return Py_BuildValue("f", tau);
}

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

//  printf("%zu\n", NA);

  double * A = (double*) PyArray_DATA(pyA);
  double * B = (double*) PyArray_DATA(pyB);
    

  if(!(NA==NB))
  {
    PyErr_SetString(PyExc_RuntimeError, "Arrays have to have the same number of elements");
    return NULL;
  }

  size_t N = NA;

  if(N%2 ==1)
  {
    PyErr_SetString(PyExc_RuntimeError, "The array has to have an even number of elements\n");
    return NULL;
  }

  N = N/2;

  // printf("N: %zu\n", N);

  /* if the above function returns -1, an appropriate Python exception will
   * have been set, and the function simply returns NULL
   */

  /* call cos from libm */
  double tau = eudist(A, B, N);

  /*  construct the output from cos, from c double to python float */
  return Py_BuildValue("f", tau);
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
