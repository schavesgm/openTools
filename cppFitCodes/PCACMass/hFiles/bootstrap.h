//------------------------------------------------------------------------------------------------//
//    AUTHOR : SERGIO CHAVES GARCÍA-MASCARAQUE
//    E-MAIL : SERGIOZTESKATE@GMAIL.COM
//    DICIEMBRE DE 2018, SWANSEA, WALES - MADRID, SPAIN
//------------------------------------------------------------------------------------------------//


#ifndef BOOTSTRAP_H
#define BOOTSTRAP_H

#include <iostream>
#include <random>
#include <ctime>

#include "./loadData.h"     // We need the matrix structure

double meanValue( const double*, const unsigned );
double stdeValue( const double*, const unsigned );

struct myMat bootstrap( const myMat inFile, const unsigned timePoints,
                       const unsigned inPosFit, unsigned ouPosFit,
                       int seedAux ) {

    if( ouPosFit <= inPosFit )      // In order to avoid errors
        ouPosFit = timePoints;

    // Number of points contained in structure
    unsigned numPoints = inFile.rowSize;
    unsigned numEnsemble = numPoints / timePoints;
    unsigned pointsUsed = ouPosFit - inPosFit;

    double* resampPoints = new double[numEnsemble*pointsUsed];
    double* bootstPoints = new double[pointsUsed];

    // INITIALIZE THE RANDOM ENGINE
    typedef std::mt19937 MyRng;
    uint32_t seedVal = seedAux;
    MyRng rng;
    rng.seed(seedVal);
    std::uniform_int_distribution<int> dist(0, numEnsemble - 1 );

    // Generate the resampled vector
    int bootChoice, auxCol = 0;
    for( unsigned i = 0; i < numEnsemble; i++ ) {
        bootChoice = dist(rng);
        for( unsigned j = inPosFit; j < ouPosFit; j++ ) {
            resampPoints[i*pointsUsed+auxCol] = inFile.m_Matrix[bootChoice*timePoints+j];
            auxCol += 1;
        }
        auxCol = 0;
    }

    // Calculate the mean of the resampled vector at each time
    unsigned auxTime = 0;
    for( unsigned i = 0; i < pointsUsed; i++ ) {
        bootstPoints[auxTime] = 0.0;
        for( unsigned j = 0; j < numEnsemble; j++ ) {
            bootstPoints[auxTime] += resampPoints[j*pointsUsed+i] / numEnsemble;
        }
        //std::cout << bootstPoints[auxTime] << std::endl;
        auxTime += 1;
    }

    struct myMat Return = { bootstPoints, pointsUsed, inFile.colSize };
    delete resampPoints;    // Delete pointer containing the resampled data
    return Return;
}

struct myMat pcacMass( const myMat dataAP, const myMat dataPP,
                       const unsigned timePoints, int seedAux ) {

    // Number of points contained in structure
    unsigned numPoints = dataAP.rowSize;
    unsigned numEnsemble = numPoints / timePoints;

    // Buffers needed
    double* buffPCAC = new double[numEnsemble*(timePoints-1)];
    double* resPCAC = new double[timePoints - 1];

    // INITIALIZE THE RANDOM ENGINE
    std::mt19937 rng;
    rng.seed( seedAux );
    std::uniform_int_distribution<int> dist(0, numEnsemble - 1 );

    // Generate the resampled vector
    int pickAP, pickPP;
    for( unsigned nE = 0; nE < numEnsemble; nE++ ) {
        pickAP = dist(rng);
        pickPP = dist(rng);
        for( unsigned t = 0; t < timePoints; t++ ) {
            buffPCAC[t+nE*(timePoints-1)] = \
                fabs( ( dataAP.m_Matrix[pickAP*timePoints+t] - \
                        dataAP.m_Matrix[pickAP*timePoints+t+1] ) ) / \
                ( dataPP.m_Matrix[pickPP*timePoints+t] );
        }
    }

    // Calculate the mean of the resampled vector at each time
    for( unsigned t = 0; t < timePoints - 1; t++ ) {
        resPCAC[t] = 0.0;
        for( unsigned nE = 0; nE < numEnsemble; nE++ ) {
            resPCAC[t] += buffPCAC[nE*(timePoints-1)+t];
        }
        resPCAC[t] = resPCAC[t] / ( numEnsemble );
    }

    struct myMat Return = { resPCAC, timePoints - 1, dataAP.colSize };
    delete buffPCAC;    // Delete pointer containing the resampled data
    return Return;
}

double meanValue( const double* inPoint, const unsigned numPoints ) {
    /*
     * Function that calculates the mean value of the vector given.
     *
     * args:
     * double* inPoint    : input vector from which we want to calculate the mean value
     * unsigned numPoints : number of points contained in the vector.
     *
     * return:
     * double             : mean value of the pointer to array given
     */

    double retMean = 0.0;
    for( unsigned i = 0; i < numPoints; i++ )
        retMean += inPoint[i];

    return retMean / numPoints;
}

double stdeValue( const double* inPoint, const unsigned numPoints ) {
    /*
     * Function that calculates the standard deviation of the vector given.
     *
     * double* inPoint    : input vector from which we want to calculate the std deviation
     * unsigned numPoints : number of points contained in the vector
     *
     * return:
     * double             : standard deviation of the pointer to array given
     */
    double retStde = 0.0;
    double meanVal = meanValue( inPoint, numPoints );   // We need the mean value
    for( unsigned i = 0; i < numPoints; i++ )
        retStde += ( inPoint[i] - meanVal ) * ( inPoint[i] - meanVal );

    return std::sqrt( retStde / ( numPoints - 1 ) );
}

struct myMat covMat( struct myMat inVal ) {
    /*
     * Function that calculates the covariance matrix given the errors, i.e, it calculates the
     * matrix containing \Sigma = \sqrt{ \sigma_i * \sigma_j }
     *
     * args:
     * struct myMat       : Matrix containing the errors of the quantities
     *
     * return:
     * double*            : Pointer to array containing the elements of the matrix
     */

    // Get the dimension of the vector
    unsigned rowSize = inVal.rowSize;

    // Pointer that will hold the covariance matrix
    double* retMat = new double[rowSize*rowSize];

    for( unsigned i = 0; i < rowSize; i++ ) {
        for( unsigned j = 0; j < rowSize; j++ ) {
            retMat[i*rowSize+j] = std::sqrt( inVal.m_Matrix[i] * inVal.m_Matrix[j] );
        }
    }

    struct myMat Mat = { retMat, rowSize, rowSize };
    return Mat;
}

#endif
