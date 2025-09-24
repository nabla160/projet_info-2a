//
//  luMethod.h
//  projet_info
//
//  Created by Alban Laborde-Laulhé on 23/09/2025.
//

#ifndef LUMETHOD_H
#define LUMETHOD_H

#include "matrice.h"

int LUPDecompose(Matrice A, double Tol, int *P);
Matrice LUPInvert(Matrice A, int *P);
double LUPDeterminant(Matrice A, int *P);

#endif


