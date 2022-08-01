#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
// #include <kmeans.h>

double **createMatrix(int rows, int columns);
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
void eigenValues(double **matrixA, double **matrixV);
double **outputJacobi(double **matrix, double **V, int n);
int eigengapHeuristic(double **matrixA, double **matrixV, int n);
double **computeU(double **matrix, double **matrixJ, int n, int k);
int *sortedIndices(int n);
int cmpfunc(const void *a, const void *b);
double absVal(double d);
double **callKmeans(int k, int size, int *initialIndices, double **vectors);
int isTXTorCSV(const char filename[]);

int vecLen;
int numOfVecs;
double *eigenVals;


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


double **initializeVectors(const char inputFile[]){
    int cols;
    int rows;

    FILE *vectorsFile;
    vectorsFile = fopen(inputFile, "r");
    double **vecList;

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
    numOfVecs = n;
    double **D;
    double d;
    int i;
    int j;

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
    double d;
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

    tmp1 =  matrixMult(matrixW, matrixD);
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
        free(V[0]);
        free(V);
        
        V = tmp;

        numRot++;
        if((offA - offNextA > EPS) && (numRot < MAXROT)){
            contLoop = 1;
        }
    }
    free(copyMatrix[0]);
    free(copyMatrix);

    return V;
}


void eigenValues(double **matrixA, double **matrixV){
    int i;
    int j;
    double val;
    double **eigenMat;
    eigenMat = matrixMult(matrixA, matrixV);

    eigenVals = calloc(numOfVecs, sizeof(double));
    for(i = 0; i < numOfVecs; i++){
        for (j = 0; j < numOfVecs; j++){
            if(matrixV[j][i] != 0){ 
                val = eigenMat[j][i] / matrixV[j][i];
                if(val != 0){
                    eigenVals[i] = val;
                    break;
                }
            }
        }
    }
    free(eigenMat[0]);
    free(eigenMat);
}


double **outputJacobi(double **matrix, double **V, int n){
    double **res;
    int i;
    int j;
    numOfVecs = n;

    eigenValues(matrix, V);
    res = createMatrix(numOfVecs + 1, numOfVecs);

    for(i = 0; i < numOfVecs; i++){
            res[0][i] = eigenVals[i];
    }
        
    for (i = 1; i < numOfVecs + 1 ; i++){
        for(j = 0; j < numOfVecs; j++){
                res[i][j] = V[i - 1][j];
        }
    }
    free(eigenVals);
    return res;
}


int eigengapHeuristic(double **matrixA, double **matrixV, int n){
    int i;
    int j;
    double gap;
    double delta;
    int k;
    int *indices;
    numOfVecs = n;

    eigenValues(matrixA, matrixV);

    indices = sortedIndices(numOfVecs);

    gap = 0;
    k = 0;
    for (i = 0; i < numOfVecs/2 ; i++){
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


double **computeU(double **matrix, double **matrixJ, int n, int k){
    int i;
    int j;
    int *indices;
    double **U;
    numOfVecs = n;
    U = createMatrix(n, k);

    eigenValues(matrix, matrixJ);

    indices = sortedIndices(numOfVecs);

    for(i = 0; i < numOfVecs; i++){
        for(j = 0; j < k; j++){
            U[i][j] = matrixJ[i][indices[j]];
        }
        
    }
    free(indices);
    free(eigenVals);
}


int *sortedIndices(int n){
    int i;
    int *indices;
    indices = calloc(numOfVecs, sizeof(int));
    numOfVecs = n;

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
    const char *goal;
    const char *fileName;
    double **datapoints;
    double **res;
    double **W;
    double **D;
    double **J;

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
        J = jacobi(datapoints, numOfVecs);
        res = outputJacobi(datapoints, J, numOfVecs);
        n = n + 1;
        free(J[0]);
        free(J);
    }
    else{
        printf("Invalid Input!\n");
        return 1;
    }

    for(int i = 0; i < n; i++){
        for(int j = 0; j < m - 1; j++){
            printf("%.4f, ", res[i][j]);
        printf("%.4f\n", res[i][numOfVecs - 1]);
        }
    }
    return 0;
}

