//
//  mactrice.h
//  projet_info
//
//  Created by Alban Laborde-Laulhé on 23/09/2025.
//

#ifndef MATRICE_H
#define MATRICE_H

//Définition d'une Matrice
typedef struct {
    int n;
    int m;
    double** data;
} Matrice;

Matrice createMatrice(int n, int m);
Matrice copieMatrice(Matrice M);
void pScaMat(Matrice M, double a);
double pScalaire(Matrice X);
Matrice addMat(Matrice M, Matrice N, double a);
double norme(Matrice M);
Matrice Identite(int n);
void afficheMatrice(Matrice mat);
void testMatrice(Matrice M);
Matrice transpose(Matrice M);
Matrice produit(Matrice M, Matrice N);

#endif

