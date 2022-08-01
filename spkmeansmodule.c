#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>

void convertVectors(PyObject *pyVecs, int n, int m);
void renormalizeMat(double **U, int n, int k);
double **vectors;

static PyObject* spk_capi(PyObject *self, PyObject *args){
    PyObject *pyVectors;
    int goal;
    int k;
    int n;
    int m;
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
        resMat = U;
        free(W[0]);
        free(W);
        free(D[0]);
        free(D);
        free(L[0]);
        free(L);
        free(J[0]);
        free(J);
        break;

    case 1: //wam
        resMat = WadjacencyMatrix(vectors, n, m);
        break;

    case 2: //ddg
        W = WadjacencyMatrix(vectors, n, m);
        resMat = computeD(W, n);
        free(W[0]);
        free(W);
        break;

    case 3: //lnorm
        W = WadjacencyMatrix(vectors, n, m);
        D = computeD(W, n);
        resMat = LnormMatrix(W, D, n);
        free(W[0]);
        free(W);
        free(D[0]);
        free(D);
        break;

    case 4: //jacobi
        J = jacobi(vectors, n);
        resMat = outputJacobi(vectors, J, n);
        free(J[0]);
        free(J);
        break;
    default:
        break;
    }

    free(vectors[0]);
    free(vectors);
    return resMat;
}


// static PyObject* kmeans_capi(PyObject *self, PyObject *args){
//     PyObject *pyVectors;
//     int K;
//     int maxIter;

//     if(!PyArg_ParseTuple(args, "iidOOii", &K, &maxIter, &epsilon, &pyVectors, &initialIndices, &vecLen, &numOfVecs)){
//         return NULL;
//     }

//     convertVectors(pyVectors);

//     if(kMeans(K, maxIter) == 0){
//         centroidsToPy(K);
//     }
//     free(vecList[0]);
//     free(vecList);
//     free(centroids[0]);
//     free(centroids);
//     free(sumClusters[0]);
//     free(sumClusters);
//     return finalCentroids;
// }

// static PyMethodDef capiMethods[] = {
//     {"execute_spk", // executes spk goals
//         (PyCFunction) spk_capi,
//         METH_VARARGS,
//         PyDoc_STR("goals algorithm implementation")},

//         {"kmeans", // executes kmeans
//         (PyCFunction) kmeans_capi,
//         METH_VARARGS,
//         PyDoc_STR("kmeans algorithm implementation")},

//     {NULL, NULL, 0, NULL} // end, sentinel for python
// };

// static struct PyModuleDef moduledef = { // defining module
//     PyModuleDef_HEAD_INIT,
//     "spkmeansmodule",
//     NULL,
//     -1,
//     capiMethods
// };

// PyMODINIT_FUNC
// PyInit_spkmeansmodule(void)
// {
//     PyObject *m;
//     m = PyModule_Create(&moduledef);
//     if (!m) {
//         return NULL;
//     }
//     return m;
// }

void convertVectors(PyObject *pyVecs, int n, int m){
    int i;
    vectors = createMatrix(n, m);
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
