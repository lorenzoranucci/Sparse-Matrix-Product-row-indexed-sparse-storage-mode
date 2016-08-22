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
int sprstm(float sa[], unsigned long ija[], float sb[], unsigned long ijb[],
           float thresh, unsigned long nmax, float sc[], unsigned long ijc[]) {
    unsigned long i, ijma, ijmb, j, k, ma, mb, mbb;
    float sum;
    if (ija[1] != ijb[1]) nrerror("sprstm: sizes do not match");//L'elemento 1 ha le size.... di cosa?
    ijc[1] = k = ija[1];//la size del risultato è uguale a quella di a e b
    for (i = 1; i <= ija[1] - 2; i++) {
        //Loop over rows of A,
        for (j = 1; j <= ijb[1] - 2; j++) {
            //and rows of B.
            if (i == j) sum = sa[i] * sb[j]; else sum = 0.0e0;
            mb = ijb[j];
            for (ma = ija[i]; ma <= ija[i + 1] - 1; ma++) {
                //Loop through elements in A’s row. Convoluted logic, following, accounts for the
                //various combinations of diagonal and off-diagonal elements.
                ijma = ija[ma];
                if (ijma == j) sum += sa[ma] * sb[j];
                else {
                    while (mb < ijb[j + 1]) {
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

int main(int argc, char *argv[]) {
    /*
     * RICORDA che questo bellissimo libro conta a partire da 1 e non da 0, che figata
     * Comportati di conseguenza
     * */
    int row=SIZE;
    int col=SIZE;

    /*Alloca Matrice del libro*/
    float **inputMat=malloc(sizeof *inputMat * (5 +1) );
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


    /*
     * Converti matrice libro
     * */
    float sa[12];
    long ija[12];
    sprsin(inputMat, 5, 0.1, 11, sa, ija);

    /*
     * Alloca matrice 2
     * */
    float **inputMat2=malloc(sizeof *inputMat2 * (row +1) );
    for (int i = 1; i <= row; i++)
    {
        /*
         * per ogni cella della prima riga alloco la colonna corrispondente,
         * in pratica espando lo spazio di memoria occupato dalla cella
         * lo spazio occupato dalla cella passa da sizeof float a sizeof float * col
         * */
        inputMat2[i] = malloc(sizeof *inputMat2[i] * (col+1));
        if (inputMat2[i])
        {
            size_t j;
            for (j = 1; j <= col; j++)
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

        int row_index = rand()%(row)+1;
        int col_index = rand()%(col)+1;
        int rand_value = rand()%1000;
        inputMat2[row_index][col_index]=rand_value;
    }

    /*
     * Print matrice 2
     * */
    printf("**********Matrice B***********\n");
    for (int i = 1; i <= row; i++){
        for (int j = 1; j <= col; j++)
        {
            printf(" |%f| ",inputMat2[i][j]);
        }
        printf("\n");
    }

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