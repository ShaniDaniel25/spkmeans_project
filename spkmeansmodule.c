#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>
#include "spkmeans.h"

void convertVectors(PyObject *pyVecs, int n, int m);
void renormalizeMat(double **U, int n, int k);
PyObject *cMatToPy(double **centroids, int n, int m);
double **vectors;

static PyObject* spk_capi(PyObject *self, PyObject *args){
    PyObject *pyVectors;
    PyObject *pyResMat;
    int goal;
    int k;
    int n;
    int m;
    int i;
    double **W;
    double **D;
    double **L;
    double **J;
    double **U;
    double **resMat;

    if(!PyArg_ParseTuple(args, "Oiiii", &pyVectors, &goal, &k, &n, &m)){
        return NULL;
    }

    convertVectors(pyVectors, n, m);
    
    switch (goal)
    {
    case 0: //spk
        W = WadjacencyMatrix(vectors, n, m);
        D = computeD(W, n);
        L = LnormMatrix(W, D, n);
        J = jacobi(L, n);
        
        if(k == 0){
            k = eigengapHeuristic(L, J, n);
        }

        U = computeU(L, J, n, k);
        renormalizeMat(U, n, k);
        
        resMat = callCreateMatrix(1, n*k + 1);
        resMat[0][0] = k;

        for(i = 1; i < n*k + 1; i++){
            resMat[0][i] = U[0][i-1];
        }

        pyResMat = cMatToPy(resMat, 1, n*k + 1);
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

        break;

    case 1: //wam
        resMat = WadjacencyMatrix(vectors, n, m);
        pyResMat  = cMatToPy(resMat, n, n);
        free(resMat[0]);
        free(resMat);
        break;

    case 2: //ddg
        W = WadjacencyMatrix(vectors, n, m);
        resMat = computeD(W, n);
        pyResMat  = cMatToPy(resMat, n, n);
        free(W[0]);
        free(W);
        free(resMat[0]);
        free(resMat);
        break;

    case 3: //lnorm
        W = WadjacencyMatrix(vectors, n, m);
        D = computeD(W, n);
        resMat = LnormMatrix(W, D, n);
        pyResMat  = cMatToPy(resMat, n, n);
        free(W[0]);
        free(W);
        free(D[0]);
        free(D);
        free(resMat[0]);
        free(resMat);
        break;

    case 4: //jacobi
        J = jacobi(vectors, n);
        resMat = outputJacobi(vectors, J, n);
        pyResMat  = cMatToPy(resMat, n+1, n);
        free(J[0]);
        free(J);
        free(resMat[0]);
        free(resMat);
        break;
    default:
        break;
    }

    free(vectors[0]);
    free(vectors);
    return pyResMat;
}


static PyObject* kmeans_capi(PyObject *self, PyObject *args){
    double **centroids;
    double *initIndices;
    PyObject *pyVectors;
    PyObject *pyCentroids;
    PyObject *initialIndices;
    int K;
    int n;
    int i;
    int maxIter;
    int centIndex;

    if(!PyArg_ParseTuple(args, "iidOOi", &K, &maxIter, &pyVectors, &initialIndices, &n)){
        return NULL;
    }

    initIndices = calloc(K, sizeof(int));
    for(i = 0; i < K; i++){
        centIndex = (int) PyLong_AsLong(PySequence_GetItem(initialIndices,i));
        initIndices[i] = centIndex;
    }

    convertVectors(pyVectors, n, K);

    if((centroids = callKmeans(K, n, initIndices, vectors)) != NULL){
        pyCentroids = cMatToPy(centroids, K, K);
    }

    free(centroids[0]);
    free(centroids);
    free(vectors[0]);
    free(vectors);

    return pyCentroids;
}


static PyMethodDef capiMethods[] = {
    {"execute_spk", // executes spk goals
        (PyCFunction) spk_capi,
        METH_VARARGS,
        PyDoc_STR("goals algorithm implementation")},

    {"kmeans", // executes kmeans
        (PyCFunction) kmeans_capi,
        METH_VARARGS,
        PyDoc_STR("kmeans algorithm implementation")},

    {NULL, NULL, 0, NULL} // end, sentinel for python
};

static struct PyModuleDef moduledef = { // defining module
    PyModuleDef_HEAD_INIT,
    "spkmeansmodule",
    NULL,
    -1,
    capiMethods
};

PyMODINIT_FUNC
PyInit_spkmeansmodule(void)
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
            U[i][j] = U[i][j] / NormRow;
        }
    }
}

PyObject *cMatToPy(double **centroids, int n, int m){
    PyObject *finalCentroids;
    int i;
    finalCentroids = PyList_New(0);
    for(i = 0; i < n*m; i++){
        PyList_Append(finalCentroids, Py_BuildValue("d", centroids[0][i]));
    }
    return finalCentroids;
}
