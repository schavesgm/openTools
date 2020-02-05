//------------------------------------------------------------------------------------------------//
//    AUTHOR : SERGIO CHAVES GARCÍA-MASCARAQUE
//    E-MAIL : SERGIOZTESKATE@GMAIL.COM
//    DICIEMBRE DE 2018, SWANSEA, WALES - MADRID, SPAIN
//------------------------------------------------------------------------------------------------//


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
     * Function that returns the LU decomposition using Doolittle's method
     *
     * args:
     * const double*    : matrix to be decompose, must be folded into an array
     * const unsigned   : dimensions of the matrix, must be a square matrix
     *
     * return:
     * struct LUcont    : structure containing lMat (L), uMat (U) and rowSize as the dimension
     *
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
                lMat[i*rowSize+j] = ( 1 / uMat[j*rowSize+j] ) * ( inMat[i*rowSize+j] - auxVar );
            }
        }
    }

    struct LUcont LU = { lMat, uMat, rowSize };
    return LU;
}

double* solveLU( const double* inMat, const double* b, const unsigned rowSize ) {
    /* 
     * Function that returns the values that are solution of A*x = b
     *
     * args:
     * const double*    : Pointer to array containing the elements of the matrix A
     * const double*    : Pointer to array containing the elements of the vector b
     * const unsigned   : Dimension of the space
     *
     * return:
     * double*          : Pointer containing the values that solve the system
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
     * Function that returns the value of the system Ax = b given the LU decompositon of A. This
     * routine is generated to recycle the LU decomposition while inverting the matrix. It is not
     * intended to be used while trying to solve the system
     *
     * args:
     * const struct LUcont  : LUcont structure containing the LU decomposition of a matrix A
     * const double*        : Pointer to array containing the elements of the vector b
     *
     * return:
     * double*              : solution of the system Ax = b
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
        xSol[auxIter] = ( 1 / LU.uMat[auxIter*rowSize+auxIter] ) * ( ySol[auxIter] - auxVar );
        auxIter -= 1;
    }
    
    return xSol;
}

double* invertMat( const double* inMat, const unsigned rowSize ) {
    /*
     * Function that calculates the inverse matrix of a given non-singular matrix whose determinant
     * is different from zero. It is based on the LU decomposition, implemented following Doolittle
     * method. We would like to obtain a matrix B such that A * B = 1
     *
     * args:
     * const double*    : Pointer to array containing the elements of the matrix A that we want to
     *                    invert. A must be a square matrix.
     * const unsigned   : Dimension of the space in which A lives
     *
     * return:
     * double*          : Pointer to array containing the elements of the inverse matrix A^{-1}
     */
    
    // Get the LU decomposition of inMat
    struct LUcont LU = LUDecompose( inMat, rowSize );

    // We need to solve N equations
    double* invMat = new double[rowSize*rowSize];

    double* auxSol = new double[rowSize];   // Auxiliar pointer that will hold the solution
    double* idVec = new double[rowSize];    // Identity matrix vector in each case

    for( unsigned i = 0; i < rowSize; i++ ) {
        // Fill the solve vector automatically
        for( unsigned j = 0; j < rowSize; j++ ) {
            if( j != i )
                idVec[j] = 0.0;
            else
                idVec[j] = 1.0;
        }
        
        // Solve the equation
        auxSol = auxSolveLU( LU, idVec );   // We use auxSolveLU as it recycles LU decomposition
        for( unsigned k = 0; k < rowSize; k++ )
            invMat[k*rowSize+i] = auxSol[k];
    }
    
    return invMat;
}
    



#endif

