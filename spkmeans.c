#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "kmeans.h"

double **callCreateMatrix(int rows, int columns);
double **initializeVectors(const char inputFile[]);
double computeDist(double *vec1, double *vec2);
double **WadjacencyMatrix(double **datapoints, int n, int m);
double **matrixMult(double **matrix1, double **matrix2);
double **subMats(double **matrix1, double **matrix2);
double **computeD(double **matrixW, int n);
double **LnormMatrix(double **matrixW, double **matrixD, int n);
double calculateT(double **matrixA, int maxi, int maxj);
void nextMatrix(double **matrixA, double c, double s, int maxi, int maxj);
double computeOff(double **matrix);
double **jacobi(double **matrix, int n);
int eigengapHeuristic(double *eigVals, int n);
double **computeU(double **matrixJ, double *eigVals, int n, int k);
int *sortedIndices(int n);
int cmpfunc(const void *a, const void *b);
double absVal(double d);
double **callKmeans(int k, int size, int *initialIndices, double **vectors);
int isTXTorCSV(const char filename[]);

int vecLen;
int numOfVecs;
double *eigenVals;


double **callCreateMatrix(int rows, int columns){
    return createMatrix(rows, columns);
}


double **initializeVectors(const char inputFile[]){
    FILE *vectorsFile;
    double **vecList;
    vectorsFile = fopen(inputFile, "r");
    
    if(vectorsFile != NULL){
        char c;
        double d;
        int i;
        int j;
        vecLen = 1;
        numOfVecs = 0;

        while ( (c = fgetc(vectorsFile)) != EOF) {
            if (c == ',' && numOfVecs == 0){
                vecLen += 1;
            }
            if (c == '\n'){
                numOfVecs += 1;
            }
        }
        fclose(vectorsFile);

        vecList = createMatrix(numOfVecs, vecLen);
        vectorsFile = fopen(inputFile, "r");
        for(i = 0; i < numOfVecs; i++){
            for(j = 0; j < vecLen - 1; j++){
                fscanf(vectorsFile, "%lf %c", &d, &c);
                vecList[i][j] = d;
            }
            fscanf(vectorsFile, "%lf", &d);
            vecList[i][vecLen - 1] = d;
        }
        fclose(vectorsFile);  
    }

    else{
        vecList = NULL;
    }
    return vecList;
}


double computeDist(double *vec1, double *vec2){
    int i;
    double dist;
    dist = 0;
    for(i = 0; i < vecLen; i++){
        dist += pow(vec1[i] - vec2[i], 2);
    }
    return sqrt(dist);
}


double **WadjacencyMatrix(double **datapoints, int n, int m){
int i;
int j;
double weight;
double **retMat;
numOfVecs = n;
vecLen = m;

retMat = createMatrix(numOfVecs, numOfVecs);
for(i = 0; i < numOfVecs; i++){
    retMat[i][i] = 0;
    for(j = i + 1; j < numOfVecs; j++){
        weight = exp((-1 * computeDist(datapoints[i], datapoints[j]))/ 2);
        retMat[i][j] = weight;
        retMat[j][i] = weight;
        }
    }

    return retMat;
}


double **matrixMult(double **matrix1, double **matrix2){
    double **retMat;
    double elem;
    int i;
    int j;
    int k;

    retMat = createMatrix(numOfVecs, numOfVecs);
    for(i = 0; i < numOfVecs; i++){
        for(j = 0; j < numOfVecs; j++){
            elem = 0;
            for(k = 0; k < numOfVecs; k++){
                elem += matrix1[i][k] * matrix2[k][j];
            }
            retMat[i][j] = elem;
        }
    }
    return retMat;
}


double **subMats(double **matrix1, double **matrix2){
    int i;
    int j;
    double **retMat;
    retMat = createMatrix(numOfVecs, numOfVecs);
    for(i = 0; i < numOfVecs; i++){
        for(j = 0; j < numOfVecs; j++){
            retMat[i][j] = matrix1[i][j] - matrix2[i][j];
        }
    }
    return retMat;
}

double **computeD(double **matrixW, int n){
    double **D;
    double d;
    int i;
    int j;
    numOfVecs = n;

    D = createMatrix(numOfVecs, numOfVecs);
    for(i = 0; i < numOfVecs; i++){
        d = 0;
        for(j = 0 ; j < numOfVecs; j++){
            d += matrixW[i][j];
            D[i][j] = 0;
        }
        D[i][i] = d;
    }
    return D;
}

double **LnormMatrix(double **matrixW, double **matrixD, int n){
    double **I;
    double **L;
    double **tmp1;
    double **tmp2; 
    int i;
    int j;
    numOfVecs = n;
    
    I = createMatrix(numOfVecs, numOfVecs);
    for(i = 0; i < numOfVecs; i++){
        for(j = i ; j < numOfVecs; j++){
            I[i][j] = 0;
            I[j][i] = 0;
        }
        I[i][i] = 1;
    }
    
    for(i = 0; i < numOfVecs; i++){
        matrixD[i][i] = 1 / (sqrt(matrixD[i][i]));
    }

    tmp1 = matrixMult(matrixW, matrixD);
    tmp2 = matrixMult(matrixD ,tmp1);
    L = subMats(I, tmp2);

    free(I[0]);
    free(I);
    free(tmp1[0]);
    free(tmp1);
    free(tmp2[0]);
    free(tmp2);

    return L;
}


double calculateT(double **matrixA, int maxi, int maxj){
    double theta;
    double t;
    int sign;

    theta = (matrixA[maxj][maxj] - matrixA[maxi][maxi]) / (2 * matrixA[maxi][maxj]);
    if(theta >= 0){
        sign = 1;
    } else {
        sign = -1;
    }

    t = sign / (absVal(theta) + sqrt(pow(theta, 2) + 1));
    return t;
}


void nextMatrix(double **matrixA, double c, double s, int maxi, int maxj){
    int r;
    double ari;
    double arj;
    double aii;
    double ajj;
    double aij;
    
    for(r = 0; r < numOfVecs; r++){
        if(r != maxi && r != maxj){
            ari = matrixA[r][maxi];
            arj = matrixA[r][maxj];
            matrixA[r][maxi] = c * ari - s * arj;
            matrixA[maxi][r] = matrixA[r][maxi];
            matrixA[r][maxj] = c * arj + s * ari;
            matrixA[maxj][r] = matrixA[r][maxj];
        }
    }

    aii = matrixA[maxi][maxi];
    ajj = matrixA[maxj][maxj];
    aij = matrixA[maxi][maxj];
    matrixA[maxi][maxi] = pow(c, 2) * aii + pow(s, 2) * ajj - 2 * s * c * aij;
    matrixA[maxj][maxj] = pow(s, 2) * aii + pow(c, 2) * ajj + 2 * s * c * aij;
    matrixA[maxi][maxj] = 0;
    matrixA[maxj][maxi] = 0;
}


double computeOff(double **matrix){
    int i;
    int j;
    double off;
    off = 0;

    for(i = 0; i < numOfVecs; i++){
        for(j = 0; j < numOfVecs; j++){
            if( i != j){
                off += pow(matrix[i][j], 2);
            }
        }
    }
    return off;
}


double **jacobi(double **matrix, int n){
    int i;
    int j;
    double t;
    double c;
    double s;
    double **P;
    double **V;
    double **tmp;
    double **copyMatrix;
    double **resMat;
    int maxi;
    int maxj;
    double EPS;
    int MAXROT;
    int contLoop;
    int numRot;
    double offA;
    double offNextA;
    numRot = 0;
    maxi = -1;
    maxj = -1;
    contLoop = 1;
    numOfVecs = n;
    V = createMatrix(numOfVecs, numOfVecs);
    copyMatrix = createMatrix(numOfVecs, numOfVecs);

    for(i = 0; i < numOfVecs; i++){
        for(j = 0; j < numOfVecs; j++){
            copyMatrix[i][j] = matrix[i][j];
        }
    }

    for(i = 0; i < numOfVecs; i++){
            for(j = 0; j < numOfVecs; j++){
                if(i == j){
                    V[i][j] = 1;
                } else {
                    V[i][j] = 0;
                }
            }
        }

    EPS = 0.00001;
    MAXROT = 100;

    offA = computeOff(copyMatrix);
    if(offA == 0){
        contLoop = 0;
    }

    while(contLoop == 1){
        contLoop = 0;
        for(i = 0; i < numOfVecs; i++){
            for(j = i + 1; j < numOfVecs; j++){
                if( maxj == -1 || (absVal(copyMatrix[i][j]) > absVal(copyMatrix[maxi][maxj]))){
                    maxi = i;
                    maxj = j;
                }

            }
        }

        t = calculateT(copyMatrix, maxi, maxj);
        c = 1 / (sqrt(pow(t, 2) + 1));
        s = t * c;

        offA = computeOff(copyMatrix);
        nextMatrix(copyMatrix, c, s, maxi, maxj);
        offNextA = computeOff(copyMatrix);
        
        P = createMatrix(numOfVecs, numOfVecs);
        for(i = 0; i < numOfVecs; i++){
            for(j = 0; j < numOfVecs; j++){
                if(i == j){
                    P[i][j] = 1;
                } else {
                    P[i][j] = 0;
                }
            }
        }
        P[maxi][maxi] = c;
        P[maxj][maxj] = c;
        P[maxi][maxj] = s;
        P[maxj][maxi] = -s;

        tmp = matrixMult(V, P);

        free(P[0]);
        free(P);
        
        for(i = 0; i < numOfVecs; i++){
            for(j = 0; j < numOfVecs; j++){
                V[i][j] = tmp[i][j];
            }
        }

        free(tmp[0]);
        free(tmp);

        numRot++;
        if((offA - offNextA > EPS) && (numRot < MAXROT)){
            contLoop = 1;
        }
    }

    eigenVals = calloc(numOfVecs, sizeof(double));
    for(i = 0; i < numOfVecs; i++){
        eigenVals[i] = copyMatrix[i][i];
    }

    free(copyMatrix[0]);
    free(copyMatrix);

    resMat = createMatrix(n+1, n);

    for(i = 0; i < n; i++){
        resMat[n][i] = eigenVals[i];
    }

    for(i = 0; i < n; i ++){
        for(j = 0; j < n; j++){
            resMat[i][j] = V[i][j];
        }
    }

    free(V[0]);
    free(V);
    free(eigenVals);

    return resMat;
}


int eigengapHeuristic(double *eigVals, int n){
    int i;
    double gap;
    double delta;
    int k;
    int *indices;
    numOfVecs = n;

    eigenVals = calloc(numOfVecs, sizeof(double));
    for(i = 0; i < n; i++){
        eigenVals[i] = eigVals[i];
    }

    indices = sortedIndices(numOfVecs);
    gap = 0;
    k = 0;

    for (i = 0; i < n/2 ; i++){
        delta = eigenVals[indices[i]] - eigenVals[indices[i+1]];
        if(delta > gap){
            gap = delta;
            k = i;
        }
    }

    free(indices);
    free(eigenVals);
    return k + 1;
}


double **computeU(double **matrixJ, double *eigVals, int n, int k){
    int i;
    int j;
    int *indices;
    double **U;
    numOfVecs = n;
    U = createMatrix(n, k);

    eigenVals = calloc(numOfVecs, sizeof(double));
    for(i = 0; i < n; i++){
        eigenVals[i] = eigVals[i];
    }

    indices = sortedIndices(numOfVecs);

    for(i = 0; i < n; i++){
        for(j = 0; j < k; j++){
            U[i][j] = matrixJ[i][indices[j]];
        }    
    }

    free(indices);
    free(eigenVals);

    return U;
}


int *sortedIndices(int n){
    int i;
    int *indices;
    numOfVecs = n;
    indices = calloc(numOfVecs, sizeof(int));

    for (i = 0; i < numOfVecs; i++){
        indices[i] = i;
    }
    
    qsort(indices, numOfVecs, sizeof(int), cmpfunc);
    return indices;
}


int cmpfunc(const void *a, const void *b){
    int idx1;
    int idx2;
    idx1 = *(int*) a;
    idx2 = *(int*) b;
    if (eigenVals[idx1] < eigenVals[idx2]){
        return 1;
    }
    else if(eigenVals[idx1] == eigenVals[idx2]){
        return 0;
    }
    return -1;
}


double absVal(double d){
    if(d > 0) return d;
    else return -d;
}

 double **callKmeans(int k, int size, int *initialIndices, double **vectors){
     return kMeans(k, size, initialIndices, vectors);
 }

int isTXTorCSV(const char filename[]){
    int len = -1;
    int isTXT = 1;
    int isCSV = 1;
    for(; filename[++len] != '\0';){
    }
    if(len < 4){
        return 0;
    }
    if(filename[len-4] != '.' || filename[len-3] != 't' || filename[len-2] != 'x' || filename[len-1] != 't'){
        isTXT = 0;
    }
    if(filename[len-4] != '.' || filename[len-3] != 'c' || filename[len-2] != 's' || filename[len-1] != 'v'){
        isCSV = 0;
    }

    if(isTXT == 1 || isCSV == 1){
        return 1;
    }
    return 0;
}


int main(int argc, char const *argv[])
{
    int n;
    int m;
    int i;
    int j;
    const char *goal;
    const char *fileName;
    double **datapoints;
    double **res;
    double **W;
    double **D;

    if(argc != 3){
        printf("Invalid Input!\n");
        return 1;
    }
    goal = argv[1];
    fileName = argv[2];

    if(isTXTorCSV(fileName) == 0){
        printf("Invalid Input!\n");
        return 1;
    }

    datapoints = initializeVectors(fileName);
    if(datapoints == NULL){
        printf("Invalid Input!\n");
        return 1;
    }
    n = numOfVecs;
    m = numOfVecs;
    
    if(strcmp(goal, "wam") == 0){
        res = WadjacencyMatrix(datapoints, numOfVecs, vecLen);
    }
    else if(strcmp(goal, "ddg") == 0){
        W = WadjacencyMatrix(datapoints, numOfVecs, vecLen);
        res = computeD(W, numOfVecs);
        free(W[0]);
        free(W);
    }
    else if(strcmp(goal, "lnorm") == 0){
        W = WadjacencyMatrix(datapoints, numOfVecs, vecLen);
        D = computeD(W, numOfVecs);
        res = LnormMatrix(W, D, numOfVecs);
        free(W[0]);
        free(W);
        free(D[0]);
        free(D);
    }
    else if(strcmp(goal, "jacobi") == 0){
        res = jacobi(datapoints, numOfVecs);

        for(i = 0; i < numOfVecs - 1; i++){
          printf("%.4f,", res[numOfVecs][i]);   
        }
        printf("%.4f\n", res[numOfVecs][numOfVecs - 1]);
    }
    
    else{
        free(datapoints[0]);
        free(datapoints);
        printf("Invalid Input!\n");
        return 1;
    }

    for(i = 0; i < n; i++){
        for(j = 0; j < m - 1; j++){
            printf("%.4f,", res[i][j]);
        }
        printf("%.4f\n", res[i][m - 1]);
    }

    free(datapoints[0]);
    free(datapoints);
    free(res[0]);
    free(res);
    return 0;
}
