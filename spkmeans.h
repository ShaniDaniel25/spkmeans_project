#ifndef SPKMEANS_H
#define SPKMEANS_H

double **callCreateMatrix(int rows, int columns);
double **WadjacencyMatrix(double **datapoints, int n, int m);
double **LnormMatrix(double **matrixW, double **matrixD, int n);
double **jacobi(double **matrix, int n);
double **outputJacobi(double **matrix, double **V, int n);
int eigengapHeuristic(double **matrixA, double **matrixV, int n);
double **computeD(double **matrixW, int n);
double **computeU(double **matrix, double **matrixJ, int n, int k);
double **callKmeans(int k, int size, int *initialIndices, double **vectors);

#endif /* #ifndef SPKMEANS_H */