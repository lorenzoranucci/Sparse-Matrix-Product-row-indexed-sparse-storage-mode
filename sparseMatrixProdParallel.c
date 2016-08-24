#include "/opt/lib/mpi/intel/14.0.2/mvapich2/1.9/include/mpi.h"
#include "iindexx.c"
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#define FINALIZE 300
#define NOT_NULL_ELEMENTS 9
#define NMAX_A 12
#define NMAX_B 12
#define NMAX_C 20
#define TRESH 0.0

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
    printf("**********Matrice A***********\n");
    int i;
    for (i = 1; i <= 5; i++){
        int j;
        for ( j = 1; j <= 5; j++)
        {
            printf(" |%f| ",inputMat[i][j]);
        }
        printf("\n");
    }

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
    printf("**********Matrice B***********\n");
    for (i = 1; i <= rows; i++){
        int j;
        for (j = 1; j <= cols; j++)
        {
            printf(" |%f| ",inputMat2[i][j]);
        }
        printf("\n");
    }

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
void sprstp(float sa[], unsigned long ija[], float sb[], unsigned long ijb[]) {
    unsigned long j, jl, jm, jp, ju, k, m, n2, noff, inc, iv;
    float v;
    n2 = ija[1];
    for (j = 1; j <= n2 - 2; j++){
        sb[j] = sa[j];
    }
    iindexx(ija[n2 - 1] - ija[1], (long *) &ija[n2 - 1], &ijb[n2 - 1]);
    jp = 0;
    for (k = ija[1]; k <= ija[n2 - 1] - 1; k++) {
        m = ijb[k] + n2 - 1;
        sb[k] = sa[m];
        for (j = jp + 1; j <= ija[m]; j++) ijb[j] = k;
        jp = ija[m];
        jl = 1;
        ju = n2 - 1;
        while (ju - jl > 1) {
            jm = (ju + jl) / 2;
            if (ija[jm] > m) ju = jm; else jl = jm;
        }
        ijb[k] = jl;


    }
    for (j = jp + 1; j < n2; j++) ijb[j] = ija[n2 - 1];
    for (j = 1; j <= n2 - 2; j++) {
        jl = ijb[j + 1] - ijb[j];
        noff = ijb[j] - 1;
        inc = 1;
        do {
            inc *= 3;
            inc++;
        } while (inc <= jl);
        do {
            inc /= 3;
            for (k = noff + inc + 1; k <= noff + jl; k++) {
                iv = ijb[k];
                v = sb[k];
                m = k;
                while (ijb[m - inc] > iv) {
                    ijb[m] = ijb[m - inc];
                    sb[m] = sb[m - inc];
                    m -= inc;
                    if (m - noff <= inc) break;
                }
                ijb[m] = iv;
                sb[m] = v;
            }
        } while (inc > 1);
    }
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

    if (rank == 0) {
        float **matA=initBookMat();
        //float **matB=initRandMat(6,6);
        float **matB=initBookMat();
        /*Convert matrices in row-indexed sparse storage mode*/
        float sA[NMAX_A];
        long ijA[NMAX_A];
        sprsin(matA, 5, 0.1, NMAX_A-1, sA, ijA);
        printf("Sprsin #1 executed\n");
        float sB[NMAX_B];
        long ijB[NMAX_B];
        sprsin(matB, 5, 0.1, NMAX_B-1, sB, ijB);
        printf("Sprsin #2 executed\n");
        /*Transpose matrix B*/
        float sBt[NMAX_B];
        long ijBt[NMAX_B];
        sprstp(sB,ijB, sBt, ijBt);
        printf("Sprstp executed\n");



        /*Distribute the workload among the slaves*/

        /*Send sA, ijA, sB, ijB to all the salves*/
        int i;
        int j;
        for (i=1; i<size; i++) {
            MPI_Send(&sA, NMAX_A, MPI_FLOAT, i,0, MPI_COMM_WORLD);
            cntSend++;
            MPI_Send(&ijA,NMAX_A, MPI_LONG, i,0, MPI_COMM_WORLD);
            cntSend++;

            MPI_Send(&sBt,  NMAX_B, MPI_FLOAT, i,0, MPI_COMM_WORLD);
            cntSend++;
            MPI_Send(&ijBt, NMAX_B, MPI_LONG, i,0, MPI_COMM_WORLD);
            cntSend++;
        }
        printf("I'm the master, sA, ijA, sBt, ijBt sent \n");

        /*
         * Iterate over rows of A
         *  Iterate over rows of B
         *      -receive when freeProcessors == 0
         *      -send(i, j, stepCounter) to next slave
         *      -increment stepCounter and decrement freeProcessor
         *
         * */
        int freeProcessors=size-1;
        int stepCounter=1;
        /*TODO improve waste of memory*/
        float sumBuffer[(ijA[1] - 2)*(ijA[1] - 2)+1];
        for (i = 1; i <= ijA[1] - 2; i++) {//iterate from 1 to N of A
            //Loop over rows of A,
            for (j = 1; j <= ijB[1] - 2; j++) {//iterate from 1 to N of B
                //and rows of B.
                int dest=freeProcessors;
                if(freeProcessors==0){
                    //receive sc, ijc, stepCounter
                    float sum;
                    int tmpIndex;
                    MPI_Recv(&sum,  1, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    dest=status.MPI_SOURCE;
                    MPI_Recv(&tmpIndex,  1, MPI_INT, dest, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    printf("I'm the master, %d° job received \n",tmpIndex);
                    sumBuffer[tmpIndex]=sum;
                    freeProcessors++;
                }

                MPI_Send(&i,  1, MPI_INT, dest, 0, MPI_COMM_WORLD);
                MPI_Send(&j,  1, MPI_INT, dest, 0, MPI_COMM_WORLD);
                MPI_Send(&stepCounter,  1, MPI_INT, dest, 0, MPI_COMM_WORLD);
                printf("I'm the master, %d° job sent \n",stepCounter);
                freeProcessors--;
                stepCounter++;
            }
        }

        /*
         * Receive last size-1 results and Send FINALIZE
         *
         * */
        for (i=1; i<size; i++) {
            float sum;
            int tmpIndex;
            MPI_Recv(&sum,  1, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            int dest=status.MPI_SOURCE;
            MPI_Recv(&tmpIndex,  1, MPI_INT, dest, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            printf("I'm the master, %d° job received \n",tmpIndex);
            sumBuffer[tmpIndex]=sum;
            //send finalize
            MPI_Send(&tmpIndex,  1, MPI_INT, dest, FINALIZE, MPI_COMM_WORLD);
            printf("I'm the master, FINALIZE sent to slave #%d \n",dest);
        }

        /*
         * Compose result
         * */
        float sC[NMAX_C];
        long ijC[NMAX_C];

        int k = ijA[1];
        int bufferElementIndex=1;//stepCounter
        for (i = 1; i <= ijA[1] - 2; i++) {//iterate from 1 to N of A
            //Loop over rows of A,
            for (j = 1; j <= ijB[1] - 2; j++) {//iterate from 1 to N of B
                if(i==j){//diagonal element
                    sC[i] = sumBuffer[bufferElementIndex];
                }
                else if (fabs(sumBuffer[bufferElementIndex]) > TRESH) {
                    if (k > NMAX_C) {
                        printf("sprstm: NMAX_C too small");
                    }
                    else{
                        sC[k] = sumBuffer[bufferElementIndex];
                        ijC[k++] = j;
                    }
                }
                bufferElementIndex++;
            }
            ijC[i + 1] = k;
        }

        printf("\nResult sc:\n");
        for(i=1;i<=NMAX_C;i++){
            printf("|%f|",sC[i]);
        }
        printf("\nBResult ijc:\n");
        for(i=1;i<=NMAX_C;i++){
            printf("|%d|",ijC[i]);
        }
    }
    else {
        printf("I'm the slave #%d \n",rank);
        /*Receive sA, ijA, sB, ijB*/
        float sA[12];
        long ijA[12];
        float sBt[20];
        long ijBt[20];
        MPI_Recv(sA,  NMAX_A, MPI_FLOAT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv(ijA, NMAX_A, MPI_LONG,  0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv(sBt,  NMAX_B, MPI_FLOAT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv(ijBt, NMAX_B, MPI_LONG,  0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        printf("I'm the slave #%d: sA, ijA, sBt, ijBt received \n",rank);

        /*
         * Receive (i, j, counter) until receive FINALIZE
         *  Given i, j compute sc, ijc
         *  Send sc, ijc, stepCounter
         * */
        int cntLocalJobs=0;
        while(1) {
            int i, j, stepCounter;
            MPI_Recv(&i, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            if (status.MPI_TAG == FINALIZE) {
                printf("I'm the slave %d. FINALIZE Received\n", rank);
                break;
            }
            MPI_Recv(&j, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            MPI_Recv(&stepCounter, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            cntLocalJobs++;
            printf("I'm the slave #%d: %d° job received \n",rank,cntLocalJobs);





            /*********************************************************************
             *********************** Compute product******************************
             ****************************************************************** */
            unsigned long  ijma, ijmb,  k, ma, mb, mbb;
            float sum;
            if (i == j) {
                sum = sA[i] * sBt[j];
            }
            else {
                sum = 0.0e0;
            }
            mb = ijBt[j];
            for (ma = ijA[i]; ma <= ijA[i + 1] - 1; ma++) {//loop over elements of row i of A
                ijma = ijA[ma];
                if (ijma == j) sum += sA[ma] * sBt[j];
                else {
                    while (mb < ijBt[j + 1]) {//loop over elements of row i of B
                        ijmb = ijBt[mb];
                        if (ijmb == i) {
                            sum += sA[i] * sBt[mb++];
                            continue;
                        } else if (ijmb < ijma) {
                            mb++;
                            continue;
                        } else if (ijmb == ijma) {
                            sum += sA[ma] * sBt[mb++];
                            continue;
                        }
                        break;
                    }
                }
            }
            for (mbb = mb; mbb <= ijBt[j + 1] - 1; mbb++) {
                //Exhaust the remainder of B’s row.
                if (ijBt[mbb] == i) sum += sA[i] * sBt[mbb];
            }

            /*********************************************************************
             **********************End Compute product****************************
             ****************************************************************** */

            /*
             * Send product
             * */

            MPI_Send(&sum,  1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&stepCounter,  1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            printf("I'm the slave #%d: %d° job computed and submitted \n",rank,cntLocalJobs);
        }
    }
    double end = MPI_Wtime();
    MPI_FINALIZE();
    double elapsed = end - start;
    printf("Tempo trascorso: %f secondi\n", elapsed );
    printf("Send: %d \n", cntSend );
    printf("Rcv: %d \n", cntRcv );
}