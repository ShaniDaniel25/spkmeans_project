#include <math.h>

double **kMeans(int k, int size, int *initialIndices, double **vectors);
double computeSquareDistance(double vector1[], double vector2[], int vecLen);
int computeCentroid(double** sumCluster, int* numOfVectors, int k, int vecLen);
double **createMatrix(int rows, int columns);
void initCentroids(double **vecs, int *indices, int len);
double **centroids;
int const MAX_ITER = 300;


double **kMeans(int k, int size, int *initialIndices, double **vectors){
    int i1;
    int j1;
    int i2;
    int j2;
    int p;
    int iterCount = 0;
    int reachedEpsilon = 0;
    int *vecCnt;
    double **sumClusters;
    vecCnt = calloc(k, sizeof(int));
    centroids = createMatrix(k, k);
    sumClusters = createMatrix(k, k);
    
    if(centroids == NULL || sumClusters == NULL){
        return 1;
    }

    initCentroids(vectors, initialIndices, k);
    
    while(iterCount < MAX_ITER && reachedEpsilon == 0){
        reachedEpsilon = 1;
        
        for (i1 = 0; i1 < size; i1++){
            
            int clusterIndex = 0;
            double minDist = computeSquareDistance(vectors[i1], centroids[0], k);
            
            for(j1 = 1; j1 < k; j1++){
                double dist = computeSquareDistance(vectors[i1], centroids[j1], k);
                if(dist < minDist){
                    minDist = dist;
                    clusterIndex = j1;
                }
            }

            for(p = 0; p < k; p++){
                sumClusters[clusterIndex][p] += vectors[i1][p];
            }
            vecCnt[clusterIndex] += 1;
        }
        
        if (computeCentroid(sumClusters, vecCnt, k, k) == 0){
            return 1;
        }
        
        for(i2 = 0; i2 < k; i2++){
            if(sqrt(computeSquareDistance(centroids[i2], sumClusters[i2], k)) > 0){
                reachedEpsilon = 0;
            }
            
            for(j2 = 0; j2 < k; j2++){
                centroids[i2][j2] = sumClusters[i2][j2];
                sumClusters[i2][j2] = 0;
            }
        }
        iterCount++;
    }
    free(vecCnt);
    free(sumClusters[0]);
    free(sumClusters);
    return centroids;
}

double **createMatrix(int rows, int columns){
    double *p;
    double **a;
    int i;
    p = calloc(rows*columns, sizeof(double));
    a = calloc(rows, sizeof(double *));
    
    for(i = 0; i < rows; i++){
        a[i] = p + i*columns;
    }
    return a;
}

int computeCentroid(double** sumCluster, int* numOfVectors, int k, int vecLen){ 
    int i;
    int j;
    for(i = 0; i < k; i++){
        int divisor = numOfVectors[i];
        if(divisor == 0){
            return 0;
        }
        for(j = 0; j < vecLen; j++){
            sumCluster[i][j] = sumCluster[i][j] / divisor;
            numOfVectors[i] = 0;
        }
    }
    return 1;
}

double computeSquareDistance(double vector1[], double vector2[], int vecLen){
    double res = 0;
    int i = 0;
    for (; i < vecLen; i++){
        double tmp = vector1[i] - vector2[i];
        res += pow(tmp, 2);
    }
    return res;
}

void initCentroids(double **vecs, int *indices, int len){
    int i;
    int j;
    int idx;

    for(i = 0; i < len; i++){
        idx = indices[i];
        for(j = 0; j < len; j++){
            centroids[i][j] = vecs[idx][j];
        }
    }
}