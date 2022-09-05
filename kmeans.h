#ifndef KMEANS_H
#define KMEANS_H

double **kMeans(int k, int size, int *initialIndices, double **vectors);
double **createMatrix(int rows, int columns);

#endif /* #ifndef KMEANS_H */