#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"
// https://docs.scipy.org/doc/numpy/reference/c-api.array.html

#include "eudist.h"

/*wrappted eudist*/
static PyObject* eudist_func(PyObject* self, PyObject* args)
{

    import_array();
    PyArrayObject *pyB, *pySz;

    /*  parse the input, from python float to c double */
    if (!PyArg_ParseTuple(args, "OO", &pyB, &pySz))
        return NULL;

    size_t nSz = (size_t) PyArray_Size((PyObject *) pySz);

    if(nSz!=3)
    {
        PyErr_SetString(PyExc_RuntimeError, "Voxel size has to be specified by three numbers");
        return NULL;
    }

    int        ndim     = PyArray_NDIM(pyB);
    npy_intp*  dims     = PyArray_DIMS(pyB);
    int        typenum  = PyArray_TYPE(pyB);

    printf("Typenum: %d\n", typenum);

    if(typenum != NPY_DOUBLE)
    {
        PyErr_SetString(PyExc_RuntimeError, "Mask has to be of type double");
        return NULL;
    }


    // Build output array for distances, D
    PyArrayObject *pyD;
    pyD = (PyArrayObject *) PyArray_NewFromDescr(&PyArray_Type,
                                                 PyArray_DescrFromType(NPY_DOUBLE), ndim, dims, NULL, NULL, 1, NULL);

    double * D = (double *) PyArray_GetPtr(pyD,(npy_intp[]){0, 0, 0} );
    double * B = (double *) PyArray_GetPtr(pyB,(npy_intp[]){0, 0, 0} );

    double * Sz = (double*) PyArray_DATA(pySz); // Voxel size

    size_t M = dims[0];
    size_t N = dims[1];
    size_t P = 1;
    if(ndim > 2)
    {
        P = dims[2];
    }

    printf("pyD: %d x %d x %d\n", M, N, P);
    printf("sz: %f %f %f\n", Sz[0], Sz[1], Sz[2]);

    edt(B, D,
        M, N, P,
        Sz[0], Sz[1], Sz[1],
        0);


    /*  construct the output */
    return pyD;
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
