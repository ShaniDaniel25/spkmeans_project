#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>
#include "spkmeans.h"


void convertVectors(PyObject *pyVecs, int n, int m);
void renormalizeMat(double **U, int n, int k);
void cMatToPy(double **vecs, int n, int m);
double **vectors;
PyObject *pyResMat;

static PyObject* spk_capi(PyObject *self, PyObject *args){
    PyObject *pyVectors;
    PyObject *initialIndices;
    int goal;
    int k;
    int n;
    int m;
    int i;
    int j;
    double **W;
    double **D;
    double **L;
    double **JEig;
    double **U;
    double **resMat;
    double **J;
    double **centroids;
    double *eigenVals;
    int *initIndices;
    int centIndex;

    if(!PyArg_ParseTuple(args, "OiiiiO", &pyVectors, &goal, &k, &n, &m, &initialIndices)){
        return NULL;
    }

    convertVectors(pyVectors, n, m);
    
    switch (goal)
    {
    case 0: //spk
        W = WadjacencyMatrix(vectors, n, m);
        D = computeD(W, n);
        L = LnormMatrix(W, D, n);
        JEig = jacobi(L, n);

        J = callCreateMatrix(n, n);
        for(i = 0; i < n; i++){
            for(j = 0; j < n; j ++){
                J[i][j] = JEig[i][j];
            }
        }

        eigenVals = calloc(n, sizeof(double));
        for(i = 0; i < n; i++){
            eigenVals[i] = JEig[n][i];
        }

        
        if(k == 0){
            k = eigengapHeuristic(eigenVals, n);
        }


        U = computeU(J, eigenVals, n, k);
        renormalizeMat(U, n, k);

        resMat = callCreateMatrix(1, n*k + 1);
        resMat[0][0] = k;

        for(i = 1; i < n*k + 1; i++){
            resMat[0][i] = U[0][i-1];
        }

        cMatToPy(resMat, 1, n*k + 1);
        free(resMat[0]);
        free(resMat);
        free(W[0]);
        free(W);
        free(D[0]);
        free(D);
        free(L[0]);
        free(L);
        free(J[0]);
        free(J);
        free(U[0]);
        free(U);
        free(JEig[0]);
        free(JEig);
        free(eigenVals);

        break;

    case 1: //wam
        resMat = WadjacencyMatrix(vectors, n, m);
        cMatToPy(resMat, n, n);
        free(resMat[0]);
        free(resMat);
        break;

    case 2: //ddg
        W = WadjacencyMatrix(vectors, n, m);
        resMat = computeD(W, n);
        cMatToPy(resMat, n, n);
        free(W[0]);
        free(W);
        free(resMat[0]);
        free(resMat);
        break;

    case 3: //lnorm
        W = WadjacencyMatrix(vectors, n, m);
        D = computeD(W, n);
        resMat = LnormMatrix(W, D, n);
        cMatToPy(resMat, n, n);
        free(W[0]);
        free(W);
        free(D[0]);
        free(D);
        free(resMat[0]);
        free(resMat);
        break;

    case 4: //jacobi
        resMat = jacobi(vectors, n);
        cMatToPy(resMat, n+1, n);
        free(resMat[0]);
        free(resMat);
        break;

    case 5:
        initIndices = calloc(k, sizeof(int));
        for(i = 0; i < k; i++){
            centIndex = (int) PyLong_AsLong(PySequence_GetItem(initialIndices,i));
            initIndices[i] = centIndex;
        }

        if((centroids = callKmeans(k, n, initIndices, vectors)) != NULL){
            cMatToPy(centroids, k, k);
        }

        free(centroids[0]);
        free(centroids);
        free(initIndices);
        break;

    default:
        break;
    }

    free(vectors[0]);
    free(vectors);
    return pyResMat;
}


static PyMethodDef capiMethods[] = {
    {"execute_spk", // python function's name
        (PyCFunction) spk_capi,
        METH_VARARGS,
        PyDoc_STR("skmeans algorithm implementation")},
    {NULL, NULL, 0, NULL} // end, sentinel for python
};

static struct PyModuleDef moduledef = { // defining module
    PyModuleDef_HEAD_INIT,
    "spkmeans",
    NULL,
    -1,
    capiMethods
};

PyMODINIT_FUNC
PyInit_spkmeans(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}

void convertVectors(PyObject *pyVecs, int n, int m){
    int i;
    vectors = callCreateMatrix(n, m);
    for(i = 0; i < n*m; i++){
        vectors[0][i] = PyFloat_AsDouble(PySequence_GetItem(pyVecs, i));
    }
}

void renormalizeMat(double **U, int n, int k){
    int i;
    int j;
    double NormRow;

    for(i = 0; i < n; i++){
        NormRow = 0;
        for(j = 0; j < k; j++){
            NormRow += pow(U[i][j], 2);
        }
        NormRow = sqrt(NormRow);

        for(j = 0; j < k; j++){
            if(NormRow == 0){
                U[i][j] = 0;
            }

            else{
                U[i][j] = U[i][j] / NormRow;
            }
        }
    }
}

void cMatToPy(double **vecs, int n, int m){ 
    int i;
    pyResMat = PyList_New(0);
    for(i = 0; i < n*m; i++){
        PyList_Append(pyResMat, Py_BuildValue("d", vecs[0][i]));
    }
}
