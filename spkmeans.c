#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <kmeans.h>

double **createMatrix(int rows, int columns);
double **initializeVectors(const char inputFile[]);
double computeDist(double *vec1, double *vec2);
double **WadjacencyMatrix(double **datapoints);
double **matrixMult(double **matrix1, double **matrix2);
double **subMats(double **matrix1, double **matrix2);
double **LnormMatrix(double **matrix);
double calculateT(double **matrixA, int maxi, int maxj);
void nextMatrix(double **matrixA, double c, double s, int maxi, int maxj);
double computeOff(double **matrix);
double **jacobi(double **matrix);
double *eigenValues(double **matrixA, double **matrixV);
int eigengapHeuristic(double **matrixA, double **matrixV);
int cmpfunc(const void *a, const void *b);
double absVal(double d);
double **callKmeans(int k, int size, int *initialIndices, double **vectors);

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


double **WadjacencyMatrix(double **datapoints){
int i;
int j;
double weight;
double **retMat;

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


double **LnormMatrix(double **matrix){
    double **I;
    double **D;
    double **L;
    double **tmp1;
    double **tmp2; 
    int i;
    int j;
    double d;
    
    I = createMatrix(numOfVecs, numOfVecs);
    for(i = 0; i < numOfVecs; i++){
        for(j = i ; j < numOfVecs; j++){
            I[i][j] = 0;
            I[j][i] = 0;
        }
        I[i][i] = 1;
    }
    
    D = createMatrix(numOfVecs, numOfVecs);
    for(i = 0; i < numOfVecs; i++){
        d = 0;
        for(j = 0 ; j < numOfVecs; j++){
            d += matrix[i][j];
            D[i][j] = 0;
        }
        D[i][i] = 1 / (sqrt(d));
    }

    tmp1 =  matrixMult(matrix, D);
    tmp2 = matrixMult(D,tmp1);
    L = subMats(I, tmp2);

    free(I[0]);
    free(I);
    free(D[0]);
    free(D);
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


double **jacobi(double **matrix){
    int i;
    int j;
    double t;
    double c;
    double s;
    double **P;
    double **V;
    double **tmp;
    int maxi;
    int maxj;
    double EPS;
    int MAXROT;
    int contLoop;
    int numRot = 0;
    double offA;
    double offNextA;
    maxi = -1;
    maxj = -1;
    contLoop = 1;
    V = createMatrix(numOfVecs, numOfVecs);

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

    while(contLoop==1){
        contLoop = 0;
        for(i = 0; i < numOfVecs; i++){
            for(j = i + 1; j < numOfVecs; j++){
                if( maxj == -1 || (absVal(matrix[i][j]) > absVal(matrix[maxi][maxj]))){
                    maxi = i;
                    maxj = j;
                }

            }
        }

        t = calculateT(matrix, maxi, maxj);
        c = 1 / (sqrt(pow(t, 2) + 1));
        s = t * c;

        offA = computeOff(matrix);
        nextMatrix(matrix, c, s, maxi, maxj);
        offNextA = computeOff(matrix);
        
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
    return V;
}

double *eigenValues(double **matrixA, double **matrixV){
    int i;
    int j;
    double val;
    double **eigenMat;
    eigenMat = matrixMult(matrixA, matrixV);

    for(i = 0; i < numOfVecs; i++){
        eigenVals = calloc(numOfVecs, sizeof(double));

        for (j = 0; i < numOfVecs; j++){
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

    return eigenVals;
}


int eigengapHeuristic(double **matrixA, double **matrixV){
    int i;
    int j;
    double gap;
    double delta;
    int k;
    int *indices;
    indices = calloc(numOfVecs, sizeof(int));

    eigenValues(matrixA, matrixV);
    
    for (i = 0; i < numOfVecs; i++){
        indices[i] = i;
    }
    
    qsort(indices, numOfVecs, sizeof(int), cmpfunc);

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
    return k + 1;
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

// TODO: ma shetzarich
int main(int argc, char const *argv[])
{
    double **mat = initializeVectors("input_1.txt");
    double **adj = WadjacencyMatrix(mat);
    double **lnorm = LnormMatrix(adj);
    double **jac = jacobi(lnorm);
    int k = eigengapHeuristic(lnorm, jac);
    for(int i = 0; i < numOfVecs; i++){
        for(int j = 0; j < numOfVecs  - 1; j++){
            printf("%.12f, ", jac[i][j]);
        }
        printf("%.12f\n", jac[i][numOfVecs - 1]);
    }
    printf("%d", k);
    return 0;
}

