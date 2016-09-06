#include "/opt/lib/mpi/intel/14.0.2/mvapich2/1.9/include/mpi.h"
#include "iindexx.c"
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#define FINALIZE 300
#define NOT_NULL_ELEMENTS 10
#define N 5
#define NMAX_A 12
#define NMAX_B 20
#define NMAX_C 50
#define TRESH 0.0

struct node {
    float value;
    unsigned long i;
    unsigned long j;
    unsigned long position;
    struct node *next;
};

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
            /*
             * La moltiplicazione degli elementi sulla diagonale tra di loro avviene solamente se i==j
             * */
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

            printf("|%d-%f|",cntStep,sum);

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

void insertSumBufferList(float value, int i, int j, struct node **head,  int Nb){
    int position=((i-1)*Nb)+j;

    struct node *newNodePtr=(struct node *) malloc(sizeof(struct node));
    newNodePtr->i=i;
    newNodePtr->j=j;
    newNodePtr->position=position;
    newNodePtr->value=value;
    newNodePtr->next=0;
    if(*head==0){//insert first element
        *head=newNodePtr;
    }
    else{
        struct node *currentPtr=*head;
        struct node *prevPtr=0;
        while(currentPtr!=0){
            if(currentPtr->position > newNodePtr->position &&
               (prevPtr==0 || prevPtr->position < newNodePtr->position)){
                if(prevPtr==0){//insert first position
                    *head=newNodePtr;
                }
                else{//insert between prev and current
                    prevPtr->next=newNodePtr;
                }
                newNodePtr->next=currentPtr;
                break;

            }
            else if(currentPtr->next==0){//insert last position
                currentPtr->next=newNodePtr;
                break;

            }
            else{
                prevPtr=currentPtr;
                currentPtr=currentPtr->next;
            }

        }
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

    float **matA=initBookMat();
    float **matB=initBookMat();
//        float **matA=initRandMat(N,N);
//        float **matB=initRandMat(N,N);
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

        clock_t ticSer = clock();
        int kSer=sprstm( sA,ijA, sB, ijB, TRESH, NMAX_C-1, sCSer, ijCSer);
        clock_t tocSer = clock();

        printf("\nSerial SPRSTM, Elapsed time: %f seconds\n", (double)(tocSer - ticSer) / CLOCKS_PER_SEC);
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
        clock_t tic = clock();

        int freeProcessors=size-1;
        struct node * bufferListHeadPtr=0;
        for (i = 1; i <= ijA[1] - 2; i++) {//iterate from 1 to N of A
            //Loop over rows of A,
            for (j = 1; j <= ijB[1] - 2; j++) {//iterate from 1 to N of B
                //and rows of B.
                int dest=freeProcessors;
                if(freeProcessors==0){
                    float sum;
                    unsigned long tmpI, tmpJ;
                    MPI_Recv(&sum,  1, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    dest=status.MPI_SOURCE;
                    MPI_Recv(&tmpI,  1, MPI_UNSIGNED_LONG, dest, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    MPI_Recv(&tmpJ,  1, MPI_UNSIGNED_LONG, dest, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    //printf("I'm the master, (%f,%d, %d) job received \n",sum,tmpI, tmpJ);
                    insertSumBufferList(sum,tmpI,tmpJ,&bufferListHeadPtr,ijB[1] - 2);
                    freeProcessors++;
                }

                MPI_Send(&i,  1, MPI_UNSIGNED_LONG, dest, 0, MPI_COMM_WORLD);
                MPI_Send(&j,  1, MPI_UNSIGNED_LONG, dest, 0, MPI_COMM_WORLD);
                freeProcessors--;
            }
        }

        /*
         * Receive last size-1 results and Send FINALIZE
         *
         * */
        for (i=1; i<size; i++) {
            float sum;
            unsigned long tmpI, tmpJ;
            MPI_Recv(&sum,  1, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            int dest=status.MPI_SOURCE;
            MPI_Recv(&tmpI,  1, MPI_UNSIGNED_LONG, dest, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            MPI_Recv(&tmpJ,  1, MPI_UNSIGNED_LONG, dest, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            //printf("I'm the master, (%d, %d) job received \n",tmpI, tmpJ);
            insertSumBufferList(sum,tmpI,tmpJ,&bufferListHeadPtr,ijB[1] - 2);
            //send finalize
            MPI_Send(&tmpI,  1, MPI_UNSIGNED_LONG, dest, FINALIZE, MPI_COMM_WORLD);
            //printf("I'm the master, FINALIZE sent to slave #%d \n",dest);
        }

        /*
         * Compose result
         * */
        float sC[NMAX_C];
        unsigned long ijC[NMAX_C];

        unsigned long k = ijA[1];
        ijC[1]=k;
        struct node *currentNodePtr=bufferListHeadPtr;
        int currentI=1;
        printf("\nSPRSTM Parallel sums: \n");
        while(currentNodePtr != 0){
            printf("|%d-%f|",currentNodePtr->position,currentNodePtr->value);


            if(currentNodePtr->i == currentNodePtr->j){//diagonal element
                sC[currentNodePtr->i] = currentNodePtr->value;
            }
            else if (fabs(currentNodePtr->value) > TRESH) {
                if (k > NMAX_C) {
                    printf("sprstm: NMAX_C too small");
                }
                else{
                    sC[k] = currentNodePtr->value;
                    ijC[k++] = currentNodePtr->j;
                }
            }
            currentNodePtr=currentNodePtr->next;

            if(currentNodePtr==0 || currentNodePtr->i > currentI){//every "i" for loop step
                ijC[currentI + 1] = k;
                currentI++;
            }
        }
        clock_t toc = clock();
        printf("\nParallel SPRSTM, Elapsed time: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
        printf("\nResult sc:\n");
        for(i=1;i<=k-1;i++){
            if(i == ijC[1]-1){
                printf("|nill|");
                continue;
            }
            printf("|%f|",sC[i]);
        }
        printf("\nResult ijc:\n");
        for(i=1;i<=k-1;i++){
            printf("|%d|",ijC[i]);
        }
        /*********************************************************************
        ********************END Execute parallel SPRSTM***********************
        ****************************************************************** */
    }
    else {
        //printf("I'm the slave #%d \n",rank);
        while(1) {
            //printf("I'm the slave %d. Receiving...\n", rank);
            unsigned long i,ijma, ijmb, j, ma, mb, mbb;
            MPI_Recv(&i, 1, MPI_UNSIGNED_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            if (status.MPI_TAG == FINALIZE) {
                //printf("I'm the slave %d. FINALIZE Received\n", rank);
                break;
            }
            MPI_Recv(&j, 1, MPI_UNSIGNED_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            //printf("I'm the slave %d. Received (%d,%d)\n", rank, i, j);

            /*********************************************************************
             *********************** Compute sum******************************
             ****************************************************************** */
            float sum;
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

            /*********************************************************************
             **********************End Compute sum****************************
             ****************************************************************** */

            /*
             * Send product
             * */
            //printf("I'm the slave %d. Computed (%f,%d,%d)\n", rank,sum, i, j);
            MPI_Send(&sum,  1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&i,  1, MPI_UNSIGNED_LONG, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&j,  1, MPI_UNSIGNED_LONG, 0, 0, MPI_COMM_WORLD);
        }
    }
    double end = MPI_Wtime();
    double elapsed = end - start;
    //printf("\nTotal elapsed time of node with rank %d: %f seconds\n", rank,elapsed );
    MPI_FINALIZE();
}