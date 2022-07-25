#define PY_SSIZE_T_CLEAN
// #include <Python.h>
#include <math.h>

void convertVectors(PyObject *pyVecs, int n, int m);
double **vectors;

static PyObject* spk_capi(PyObject *self, PyObject *args){
    PyObject *pyVectors;
    int goal;
    int kZero;
    int n;
    int m;
    double **W;
    double **D;
    double **resMat;

    if(!PyArg_ParseTuple(args, "Oiiii", &pyVectors, &goal, &kZero, &n, &m)){
        return NULL;
    }

    convertVectors(pyVectors, n, m);
    
    switch (goal)
    {
    case 0: //spk
        /* code */
        break;
    case 1: //wam
        resMat = WadjacencyMatrix(vectors, n, m);
        break;
    case 2: //ddg
        W = WadjacencyMatrix(vectors, n, m);
        break;
    case 3:
        /* code */
        break;
    case 4:
        /* code */
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