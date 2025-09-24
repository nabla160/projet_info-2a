//
//  main.c
//  projet_info
//
//  Created by Alban Laborde-Laulhé on 23/09/2025.
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../includes/matrice.h"
#include "../includes/luMethod.h"

int main(int argc, const char * argv[]) {
    Matrice M = createMatrice(2, 2);
    M.data[0][0] = 1; M.data[0][1] = 4; M.data[1][0] = 0; M.data[1][1] = 9;
    afficheMatrice(M);
    Matrice MLU = copieMatrice(M);
    
    printf("Et maintenant calculons le déterminant \n");
    int* P = malloc( MLU.n * sizeof(int));
    LUPDecompose(MLU, 0.1, P);
    double detM = LUPDeterminant(MLU, P);
    printf("%8.2f \n", detM);
    printf("Et si on calculait son inverse ? \n");
    Matrice IM = LUPInvert(MLU, P);
    afficheMatrice(IM);
    printf("Et maintenant : la matrice identité I3 : \n");
    afficheMatrice(Identite(3));
    
    
    return EXIT_SUCCESS;
}