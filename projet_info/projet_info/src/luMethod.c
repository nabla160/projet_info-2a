//
//  luMethod.c
//  projet_info
//
//  Created by Alban Laborde-Laulhé on 23/09/2025.
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "matrice.h"
#include "luMethod.h"

/* INPUT: A - array of pointers to rows of a square matrix having dimension N
 *        Tol - small tolerance number to detect failure when the matrix is near degenerate
 * OUTPUT: Matrix A is changed, it contains a copy of both matrices L-E and U as A=(L-E)+U such that P*A=L*U.
 *        The permutation matrix is not stored as a matrix, but in an integer vector P of size N+1
 *        containing column indexes where the permutation matrix has "1". The last element P[N]=S+N,
 *        where S is the number of row exchanges needed for determinant computation, det(P)=(-1)^S
 */
int LUPDecompose(Matrice A, double Tol, int *P) {
    if(A.n != A.m)
    {
        printf("Décomposition LU sur une matrice non carré !!!\n");
        exit(EXIT_FAILURE);
    }

    int i, j, k, imax;
    int N = A.n;
    double maxA, *ptr, absA;

    for (i = 0; i <= N; i++)
        P[i] = i; //Unit permutation matrix, P[N] initialized with N

    for (i = 0; i < N; i++) {
        maxA = 0.0;
        imax = i;

        for (k = i; k < N; k++)
            if ((absA = fabs(A.data[k][i])) > maxA) {
                maxA = absA;
                imax = k;
            }

        if (maxA < Tol) return 0; //failure, matrix is degenerate

        if (imax != i) {
            //pivoting P
            j = P[i];
            P[i] = P[imax];
            P[imax] = j;

            //pivoting rows of A
            ptr = A.data[i];
            A.data[i] = A.data[imax];
            A.data[imax] = ptr;

            //counting pivots starting from N (for determinant)
            P[N]++;
        }

        for (j = i + 1; j < N; j++) {
            A.data[j][i] /= A.data[i][i];

            for (k = i + 1; k < N; k++)
                A.data[j][k] -= A.data[j][i] * A.data[i][k];
        }
    }
    
    return 1;  //decomposition done
}


/* INPUT: A,P filled in LUPDecompose; N - dimension
 * OUTPUT: IA is the inverse of the initial matrix
 */
Matrice LUPInvert(Matrice A, int *P) {
    
    Matrice IA = copieMatrice(A);
    
    if(A.n != A.m)
    {
        printf("Inversion LU sur une matrice non carré !!!\n");
        exit(EXIT_FAILURE);
    }
    int N = A.n;
  
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
            IA.data[i][j] = P[i] == j ? 1.0 : 0.0;

            for (int k = 0; k < i; k++)
                IA.data[i][j] -= A.data[i][k] * IA.data[k][j];
        }

        for (int i = N - 1; i >= 0; i--) {
            for (int k = i + 1; k < N; k++)
                IA.data[i][j] -= A.data[i][k] * IA.data[k][j];

            IA.data[i][j] /= A.data[i][i];
        }
    }
    
    return IA;
}

/* INPUT: A,P filled in LUPDecompose; N - dimension.
 * OUTPUT: Function returns the determinant of the initial matrix
 */
double LUPDeterminant(Matrice A, int *P) {
    
    if(A.n != A.m)
    {
        printf("Calcul de determinant LU sur une matrice non carré !!!\n");
        exit(EXIT_FAILURE);
    }
    int N = A.n;

    double det = A.data[0][0];

    for (int i = 1; i < N; i++)
        det *= A.data[i][i];

    return (P[N] - N) % 2 == 0 ? det : -det;
}

