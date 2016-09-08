#include "/opt/lib/mpi/intel/14.0.2/mvapich2/1.9/include/mpi.h"
#include "iindexx.c"
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#define NOT_NULL_ELEMENTS 1000
#define N 1000
#define NMAX_A 5000
#define NMAX_B 5000
#define NMAX_C 1000000
#define THRESH 0.0



float ** initBookMat(){
    float ** inputMat=malloc(sizeof *inputMat * (5 +1) );
    inputMat[1]=malloc(sizeof *inputMat[1] * (5+1));
    inputMat[1][1]=3;
    inputMat[1][2]=0;
    inputMat[1][3]=1;
    inputMat[1][4]=0;
    inputMat[1][5]=0;
    inputMat[2]=malloc(sizeof *inputMat[2] * (5+1));
    inputMat[2][1]=0;
    inputMat[2][2]=4;
    inputMat[2][3]=0;
    inputMat[2][4]=0;
    inputMat[2][5]=0;
    inputMat[3]=malloc(sizeof *inputMat[3] * (5+1));
    inputMat[3][1]=0;
    inputMat[3][2]=7;
    inputMat[3][3]=5;
    inputMat[3][4]=9;
    inputMat[3][5]=0;
    inputMat[4]=malloc(sizeof *inputMat[4] * (5+1));
    inputMat[4][1]=0;
    inputMat[4][2]=0;
    inputMat[4][3]=0;
    inputMat[4][4]=0;
    inputMat[4][5]=2;
    inputMat[5]=malloc(sizeof *inputMat[5] * (5+1));
    inputMat[5][1]=0;
    inputMat[5][2]=0;
    inputMat[5][3]=0;
    inputMat[5][4]=6;
    inputMat[5][5]=5;

    /*
     * Print matrice libro
     * */
//    printf("**********Matrice A***********\n");
//    int i;
//    for (i = 1; i <= 5; i++){
//        int j;
//        for ( j = 1; j <= 5; j++)
//        {
//            printf(" |%f| ",inputMat[i][j]);
//        }
//        printf("\n");
//    }

    return inputMat;
}
float ** initRandMat(int rows, int cols){
    float **inputMat2=malloc(sizeof *inputMat2 * (rows +1) );
    int i;
    for (i = 1; i <= rows; i++)
    {
        /*
         * per ogni cella della prima riga alloco la colonna corrispondente,
         * in pratica espando lo spazio di memoria occupato dalla cella
         * lo spazio occupato dalla cella passa da sizeof float a sizeof float * col
         * */
        inputMat2[i] = malloc(sizeof *inputMat2[i] * (cols+1));
        if (inputMat2[i])
        {
            size_t j;
            for (j = 1; j <= cols; j++)
            {
                inputMat2[i][j] = 0;
                if(i==j){
                    inputMat2[i][j] = 1;
                }
            }
        }
    }

    /*
     * Set sparse values in random way
     * */
    srand(time(NULL));
    for(i=0; i< NOT_NULL_ELEMENTS;i++){

        int rows_index = rand()%(rows)+1;
        int col_index = rand()%(cols)+1;
        int rand_value = rand()%1000;
        inputMat2[rows_index][col_index]=rand_value;
    }

    /*
     * Print matrice 2
     * */
//    printf("**********Matrice B***********\n");
//    for (i = 1; i <= rows; i++){
//        int j;
//        for (j = 1; j <= cols; j++)
//        {
//            printf(" |%f| ",inputMat2[i][j]);
//        }
//        printf("\n");
//    }

    return inputMat2;
}
/*
 * Converts matrix a, represented with classic structure, to row-indexed structure
 * */
void sprsin(float **a, int n, float thresh, unsigned long nmax, float sa[], unsigned long ija[]) {
    int i, j;
    unsigned long k;
    //Store diagonal elements.

    for (j = 1; j <= n; j++) {
        sa[j] = a[j][j];
    }

    //Index to 1st row off-diagonal element, if any.
    ija[1] = n + 2;

    k = n + 1;
    //Loop over rows.
    for (i = 1; i <= n; i++) {
        //Loop over columns.
        for (j = 1; j <= n; j++) {
            if (fabs(a[i][j]) >= thresh && i != j) {
                if (++k > nmax) nrerror("sprsin: nmax too small");
                //Store off-diagonal elements and their columns.
                sa[k] = a[i][j];
                ija[k] = j;
            }
        }
        //As each row is completed, store index to next.
        ija[i + 1] = k + 1;
    }
}
/*Construct the transpose of a sparse square matrix, from row-index sparse storage arrays sa and
        ija into arrays sb and ijb .*/

/*
 * sprstm
 * */
/*
 * Notice that, according to the storage rules, the value of N
(namely 5) is ija[1]-2, and the length of each array is ija[ija[1]-1]-1, namely 11.
The diagonal element in row i is sa[i], and the off-diagonal elements in that row are in
sa[k] where k loops from ija[i] to ija[i+1]-1, if the upper limit is greater or equal to
the lower one (as in C’s for loops).
 * */
int sprstm(float sa[], unsigned long ija[], float sb[], unsigned long ijb[],
           float thresh, unsigned long nmax, float sc[], unsigned long ijc[]) {
    unsigned long i, ijma, ijmb, j, k, ma, mb, mbb;
    float sum;
    if (ija[1] != ijb[1]) nrerror("sprstm: sizes do not match");//Check if matrices have the same size
    ijc[1] = k = ija[1];//c size is the same as input matrices size
    printf("\nSPRSTM Serial sums: \n");
    int cntStep=0;
    for (i = 1; i <= ija[1] - 2; i++) {//iterate from 1 to N of A
        //Loop over rows of A,
        for (j = 1; j <= ijb[1] - 2; j++) {//iterate from 1 to N of B
            //and rows of B.
            cntStep++;
            if (i == j) {
                sum = sa[i] * sb[j];
            }
            else {
                sum = 0.0e0;
            }



            mb = ijb[j];
            for (ma = ija[i]; ma <= ija[i + 1] - 1; ma++) {//loop over elements of row i of A
                ijma = ija[ma];
                if (ijma == j) sum += sa[ma] * sb[j];
                else {
                    while (mb < ijb[j + 1]) {//loop over elements of row i of B
                        ijmb = ijb[mb];
                        if (ijmb == i) {
                            sum += sa[i] * sb[mb++];
                            continue;
                        } else if (ijmb < ijma) {
                            mb++;
                            continue;
                        } else if (ijmb == ijma) {
                            sum += sa[ma] * sb[mb++];
                            continue;
                        }
                        break;
                    }
                }
            }
            for (mbb = mb; mbb <= ijb[j + 1] - 1; mbb++) {
                //Exhaust the remainder of B’s row.
                if (ijb[mbb] == i) sum += sa[i] * sb[mbb];
            }

            //printf("|%d-%f|",cntStep,sum);

            if (i == j)
                sc[i] = sum;
                //Where to put the answer...
            else if (fabs(sum) > thresh) {
                if (k > nmax) nrerror("sprstm: nmax too small");
                sc[k] = sum;
                ijc[k++] = j;
            }
        }
        ijc[i + 1] = k;
    }
    return k;
}



int main(int argc, char *argv[]) {


    int cntSend=0;
    int cntRcv=0;
    int size, rank, rc;
    MPI_Status status;
    rc = MPI_INIT (&argc, &argv);
    double start = MPI_Wtime();
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &size);

    //float **matA=initBookMat();
    //float **matB=initBookMat();
    float **matA=initRandMat(N,N);
    float **matB=initRandMat(N,N);
    /*Convert matrices in row-indexed sparse storage mode*/
    float sA[NMAX_A];
    unsigned long ijA[NMAX_A];
    sprsin(matA, N, 0.1, NMAX_A-1, sA, ijA);
    float sB[NMAX_B];
    unsigned long ijB[NMAX_B];
    sprsin(matB, N, 0.1, NMAX_B-1, sB, ijB);

    if (rank == 0) {
        unsigned long i , j;
        /*********************************************************************
        ********************** Execute serial SPRSTM**************************
        ****************************************************************** */
        float *sCSer=malloc(sizeof *sCSer * (NMAX_C) );
        unsigned long *ijCSer=malloc(sizeof *ijCSer * (NMAX_C) );

        double ticSer = MPI_Wtime();
        int kSer=sprstm( sA,ijA, sB, ijB, THRESH, NMAX_C-1, sCSer, ijCSer);
        double tocSer = MPI_Wtime();

        printf("\nSerial SPRSTM, Elapsed time: %f seconds\n", (tocSer - ticSer) );
        printf("\nResult sc:\n");
        for(i=1;i<=kSer-1;i++){
            if(i == ijCSer[1]-1){
                printf("|nill|");
                continue;
            }
            printf("|%f|",sCSer[i]);
        }
        printf("\nResult ijc:\n");
        for(i=1;i<=kSer-1;i++){
            printf("|%d|",ijCSer[i]);
        }
        /*********************************************************************
        *********************END Execute serial SPRSTM************************
        ****************************************************************** */

        /*********************************************************************
        ********************** Execute parallel SPRSTM**************************
        ****************************************************************** */
        double tic = MPI_Wtime();

        /*Distribute workload*/
        int rowsPerNode=(ijA[1]-2)/size;
        for(i=1; i<size;i++){//iterate over processors
            unsigned long startIdx=((i-1)*rowsPerNode)+1;
            unsigned long endIdx=startIdx+rowsPerNode;
            if(endIdx>ijA[1]-2){
                endIdx=ijA[1]-2;
            }
            MPI_Send(&startIdx,  1, MPI_UNSIGNED_LONG, i, 0, MPI_COMM_WORLD);
            MPI_Send(&endIdx,  1, MPI_UNSIGNED_LONG, i, 0, MPI_COMM_WORLD);
            //printf("Processor: %d, start: %d, end: %d\n", i, startIdx, endIdx);
        }

        /*Receive & compose result*/
        float sc[NMAX_C];
        unsigned long ijc[NMAX_C];
        int k=ijA[1];
        ijc[1]=k;
        int isNmaxTooSmall=0;
        unsigned long currentI=1;
        for(i=1; i<size;i++) {//iterate over processors
            unsigned long cntSums;
            MPI_Recv(&cntSums, 1, MPI_UNSIGNED_LONG, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            float sums[cntSums];
            unsigned long sumsI[cntSums];
            unsigned long sumsJ[cntSums];
            MPI_Recv(&sums , cntSums, MPI_FLOAT, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            MPI_Recv(&sumsI, cntSums, MPI_UNSIGNED_LONG, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            MPI_Recv(&sumsJ, cntSums, MPI_UNSIGNED_LONG, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            if(isNmaxTooSmall==1)continue;

            /*compose result*/
            for(j=1; j<cntSums;j++){
                if(sumsI[j]==sumsJ[j]){
                    sc[sumsI[j]] = sums[j];
                }
                else{
                    if(k>NMAX_C){
                        isNmaxTooSmall=1;
                        break;
                    }
                    sc[k] = sums[j];
                    ijc[k] = sumsJ[j];
                    k++;
                }

                if(sumsI[j+1] > currentI){//every "i" for loop step
                    ijc[currentI + 1] = k;
                    currentI++;
                }
            }

        }

        double toc = MPI_Wtime();

        /*Print Parallel SPRSTM Result*/
        if(isNmaxTooSmall==1){
            printf("NMAX_C too small");
        }
        else{
            printf("\nParallel SPRSTM, Elapsed time: %f seconds\n", (toc - tic));
            printf("\nResult sc:\n");
            for(i=1;i<=k-1;i++){
                if(i == ijc[1]-1){
                    printf("|nill|");
                    continue;
                }
                printf("|%f|",sc[i]);
            }
            printf("\nResult ijc:\n");
            for(i=1;i<=k-1;i++){
                printf("|%d|",ijc[i]);
            }
        }

        /*********************************************************************
        ********************END Execute parallel SPRSTM***********************
        ****************************************************************** */
    }
    else {
        unsigned long startIdx, endIdx, cntSc, cntIjc,i,ijma, ijmb, j, ma, mb, mbb, k;
        MPI_Recv(&startIdx, 1, MPI_UNSIGNED_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv(&endIdx, 1, MPI_UNSIGNED_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        unsigned long rowsAToCompute=endIdx-startIdx+1;
        unsigned long resultSize=rowsAToCompute*ijB[1] - 2;
        float sum, sums[resultSize];
        unsigned long sumsI[resultSize], sumsJ[resultSize];
        unsigned long cntSums=1;
        for (i = startIdx; i <= endIdx; i++) {
            //Loop over rows of A,
            for (j = 1; j <= ijB[1] - 2; j++) {//iterate from 1 to N of B
                //and rows of B.

                if (i == j) {
                    sum = sA[i] * sB[j];
                }
                else {
                    sum = 0.0e0;
                }



                mb = ijB[j];
                for (ma = ijA[i]; ma <= ijA[i + 1] - 1; ma++) {//loop over elements of row i of A
                    ijma = ijA[ma];
                    if (ijma == j) sum += sA[ma] * sB[j];
                    else {
                        while (mb < ijB[j + 1]) {//loop over elements of row i of B
                            ijmb = ijB[mb];
                            if (ijmb == i) {
                                sum += sA[i] * sB[mb++];
                                continue;
                            } else if (ijmb < ijma) {
                                mb++;
                                continue;
                            } else if (ijmb == ijma) {
                                sum += sA[ma] * sB[mb++];
                                continue;
                            }
                            break;
                        }
                    }
                }
                for (mbb = mb; mbb <= ijB[j + 1] - 1; mbb++) {
                    //Exhaust the remainder of B’s row.
                    if (ijB[mbb] == i) sum += sA[i] * sB[mbb];
                }

                if (i == j || fabs(sum) > THRESH){
                    sums[cntSums]=sum;
                    sumsI[cntSums]=i;
                    sumsJ[cntSums]=j;
                    cntSums++;
                }

            }
        }
        /*
         * Send result
         * */

        MPI_Send(&cntSums, 1, MPI_UNSIGNED_LONG, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&sums,  cntSums, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&sumsI, cntSums, MPI_UNSIGNED_LONG, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&sumsJ, cntSums, MPI_UNSIGNED_LONG, 0, 0, MPI_COMM_WORLD);
        //printf("Slave: %d, cntSums: %d\n", rank, cntSums);

    }
    double end = MPI_Wtime();
    double elapsed = end - start;
    //printf("\nTotal elapsed time of node with rank %d: %f seconds\n", rank,elapsed );
    MPI_FINALIZE();
}