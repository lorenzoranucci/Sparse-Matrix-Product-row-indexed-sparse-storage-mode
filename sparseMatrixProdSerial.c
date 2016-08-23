#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "iindexx.c"
#define SIZE 5
#define NOT_NULL_ELEMENTS 9



/*
 * Converts matrix a, represented with classic structure, to row-indexed structure
 * */
void sprsin(float **a, int n, float thresh, unsigned long nmax, float sa[],
            unsigned long ija[]) {
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
    for (i = 1; i <= ija[1] - 2; i++) {//iterate from 1 to N of A
        //Loop over rows of A,
        for (j = 1; j <= ijb[1] - 2; j++) {//iterate from 1 to N of B
            //and rows of B.

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
}

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
    for (int i = 1; i <= 5; i++){
        for (int j = 1; j <= 5; j++)
        {
            printf(" |%f| ",inputMat[i][j]);
        }
        printf("\n");
    }

    return inputMat;
}
float ** initRandMat(int rows, int cols){
    float **inputMat2=malloc(sizeof *inputMat2 * (rows +1) );
    for (int i = 1; i <= rows; i++)
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
    for(int i=0; i< NOT_NULL_ELEMENTS;i++){

        int rows_index = rand()%(rows)+1;
        int col_index = rand()%(cols)+1;
        int rand_value = rand()%1000;
        inputMat2[rows_index][col_index]=rand_value;
    }

    /*
     * Print matrice 2
     * */
    printf("**********Matrice B***********\n");
    for (int i = 1; i <= rows; i++){
        for (int j = 1; j <= cols; j++)
        {
            printf(" |%f| ",inputMat2[i][j]);
        }
        printf("\n");
    }

    return inputMat2;
}


int main(int argc, char *argv[]) {
    /*
     * RICORDA che questo bellissimo libro conta a partire da 1 e non da 0, che figata
     * Comportati di conseguenza
     * */

    /*Alloca Matrice del libro*/
    float **inputMat=initBookMat();

    /*
     * Converti matrice libro
     * */
    float sa[12];
    long ija[12];
    sprsin(inputMat, 5, 0.1, 11, sa, ija);

    /*
     * Alloca matrice 2
     * */
    float **inputMat2=initRandMat(6,6);

    /*
     * Converti matrice2
     * */
    float sa2[20];
    long ija2[20];
    sprsin(inputMat2, 5, 0.1, 19, sa2, ija2);


    /*
     * Trasposta matrice 2
     * */
    float sb2[20];
    long ijb2[20];
    sprstp(sa2,ija2, sb2, ijb2);

    /*
     * Moltiplica A * B^T
     * */
    float sc[100];
    long ijc[100];
    sprstm( sa,ija, sb2, ijb2, 0.1, 99, sc, ijc);

}

