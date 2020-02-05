//    AUTHOR : SERGIO CHAVES GARCÍA-MASCARAQUE
//    E-MAIL : SERGIOZTESKATE@GMAIL.COM
//    DICIEMBRE DE 2018, SWANSEA, WALES - MADRID, SPAIN


#ifndef LINALG_H
#define LINALG_H

#include <stdexcept>
#include "./loadData.h"

struct LUcont {         // Structure containing the LU decomposition
    double* lMat;
    double* uMat;
    unsigned rowSize;
};

struct LUcont LUDecompose( const double* inMat, const unsigned rowSize ) {
    /*
       Function that returns the LU decomposition using Doolittle's method.

        Args:
            const double*:
                Matrix to be decomposed, must be vectorised.
            const unsigned:
                Dimension of the matrix, must be a square matrix.

        Return:
            struct LUcont:
                Structure containing the LU decomposition of the initial
                matrix.
    */

    // Generate the pointers that will hold the LU decomposition
    double* uMat = new double[rowSize*rowSize];
    double* lMat = new double[rowSize*rowSize];

    // Fill the matrix lMat as identity, fill uMat as zeros
    for( unsigned i = 0; i < rowSize; i++ ) {
        for( unsigned j = 0; j < rowSize; j++ ) {
            uMat[i*rowSize+j] = 0.0;
            if( i == j )
                lMat[i*rowSize+j] = 1.0;
            else
                lMat[i*rowSize+j] = 0.0;
        }
    }

    // Generate the decomposition
    double auxVar;
    for( unsigned i = 0; i < rowSize; i++ ) {
        for( unsigned j = 0; j < rowSize; j++ ) {
            if( i <= j ) {
                auxVar = 0.0;
                for( unsigned k = 0; k < i; k++ ) {
                    auxVar += lMat[i*rowSize+k] * uMat[k*rowSize+j];
                }
                uMat[i*rowSize+j] = inMat[i*rowSize+j] - auxVar;
            }
            if( i > j ) {
                auxVar = 0.0;
                for( unsigned k = 0; k < j; k++ ) {
                    auxVar += lMat[i*rowSize+k] * uMat[k*rowSize+j];
                }
                lMat[i*rowSize+j] = ( 1 / uMat[j*rowSize+j] ) *
                                    ( inMat[i*rowSize+j] - auxVar );
            }
        }
    }

    struct LUcont LU = { lMat, uMat, rowSize };
    return LU;
}

double* solveLU( const double* inMat, const double* b, const unsigned rowSize ) {
    /*
       Function that returns the values that are solution of A * x = b.

        Args:
            const double*:
                Matrix representing A, must be vectorised.
            const double*:
                Vector representing b.
            const unsigned:
                Dimension of the matrix A and b, note that A has to be a
                square matrix.

        Return:
            double*:
                Pointer representing the vector x that fullfils A * x = b
    */

    // Get the LU decomposition of the matrix
    struct LUcont LU = LUDecompose( inMat, rowSize );

    // We must solve y and x using Ly = b, Ux = y
    double* xSol = new double[rowSize];
    double* ySol = new double[rowSize];

    for( unsigned i = 0; i < rowSize; i++ ) {
        xSol[i] = 0.0;
        ySol[i] = 0.0;
    }

    // First we solve Ly = b
    double auxVar;
    for( unsigned i = 0; i < rowSize; i++ ) {
        auxVar = 0.0;
        for( unsigned j = 0; j < i; j++ ) {
            auxVar += LU.lMat[i*rowSize+j] * ySol[j];
        }
        ySol[i] = ( 1 / LU.lMat[i*rowSize+i] ) * ( b[i] - auxVar );
    }

    // Then we solve for Ux = y
    unsigned auxIter = rowSize - 1;
    for( unsigned i = 0; i < rowSize; i++ ) {
        auxVar = 0.0;
        for( unsigned j = i - 1; j <= rowSize; j++ ) {
            auxVar += LU.uMat[auxIter*rowSize+j] * xSol[j];
        }
        xSol[auxIter] = ( 1 / LU.uMat[auxIter*rowSize+auxIter] ) * ( ySol[auxIter] - auxVar );
        auxIter -= 1;
    }

    return xSol;
}

double* auxSolveLU( const struct LUcont LU, const double* b ) {
    /*
       Function that returns the values that are solution of A * x = b given
       the LU decomposition of the matrix A. It is written to be called by
       another function.

        Args:
            const struct LUcont*:
                Structure containg the LU decomposition of A.
            const double*:
                Vector representing b.
        Return:
            double*:
                Pointer representing the vector x that fullfils A * x = b
    */

    // Get the dimension of the space
    unsigned rowSize = LU.rowSize;

    // We must solve y and x using Ly = b, Ux = y
    double* xSol = new double[rowSize];
    double* ySol = new double[rowSize];

    for( unsigned i = 0; i < rowSize; i++ ) {
        xSol[i] = 0.0;
        ySol[i] = 0.0;
    }

    // First we solve Ly = b
    double auxVar;
    for( unsigned i = 0; i < rowSize; i++ ) {
        auxVar = 0.0;
        for( unsigned j = 0; j < i; j++ ) {
            auxVar += LU.lMat[i*rowSize+j] * ySol[j];
        }
        ySol[i] = ( 1 / LU.lMat[i*rowSize+i] ) * ( b[i] - auxVar );
    }

    // Then we solve for Ux = y
    unsigned auxIter = rowSize - 1;
    for( unsigned i = 0; i < rowSize; i++ ) {
        auxVar = 0.0;
        for( unsigned j = i - 1; j <= rowSize; j++ ) {
            auxVar += LU.uMat[auxIter*rowSize+j] * xSol[j];
        }
        xSol[auxIter] = ( 1 / LU.uMat[auxIter*rowSize+auxIter] )
                        * ( ySol[auxIter] - auxVar );
        auxIter -= 1;
    }

    return xSol;
}

double* invertMat( const double* inMat, const unsigned rowSize ) {
    /*
       Function that returns the inverse matrix of the matrix A. The matrix
       given has to be non-singlular, det(A) != 0. It is based on the LU
       decomposition, it finds the matrix B such that, A * B = 1

        Args:
            const double*:
                Matrix representing A, must be vectorised.
            const unsigned:
                Dimnesion of the matrix A. A must be a square matrix.

        Return:
            double*:
                Pointer representing the matrix B such that A * B = 1. The
                matrix if vectorised. To access the element (i,j) just use
                B[i * rowSize + j]
    */

    // Get the LU decomposition of inMat
    struct LUcont LU = LUDecompose( inMat, rowSize );

    // We need to solve N equations
    double* invMat = new double[rowSize*rowSize];

    // Auxiliar pointer that will hold the solution
    double* auxSol = new double[rowSize];
    // Identity matrix vector in each case
    double* idVec = new double[rowSize];

    for( unsigned i = 0; i < rowSize; i++ ) {
        // Fill the solve vector automatically
        for( unsigned j = 0; j < rowSize; j++ ) {
            if( j != i )
                idVec[j] = 0.0;
            else
                idVec[j] = 1.0;
        }

        // Solve the equation
        auxSol = auxSolveLU( LU, idVec );
        for( unsigned k = 0; k < rowSize; k++ )
            invMat[k*rowSize+i] = auxSol[k];
    }

    return invMat;
}

#endif

