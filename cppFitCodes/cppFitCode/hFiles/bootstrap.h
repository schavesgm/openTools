//    AUTHOR : SERGIO CHAVES GARCÍA-MASCARAQUE
//    E-MAIL : SERGIOZTESKATE@GMAIL.COM
//    DICIEMBRE DE 2018, SWANSEA, WALES - MADRID, SPAIN


#ifndef BOOTSTRAP_H
#define BOOTSTRAP_H

#include <iostream>
#include <random>
#include <ctime>

#include "./loadData.h"     // We need the matrix structure

double meanValue( const double*, const unsigned );

struct myMat pickAnalyzeBootstrap( const myMat inFile, const unsigned timePoints,
                              const unsigned inPosFit, unsigned ouPosFit,
                              int seedAux ) {
    /*
        Function to generate a bootstrapped resample from a given set of
        values.

        Args:
            const struct myMat:
                Input file to be resampled.
            const unsigned:
                Number of time points of the configuration.
            const unsigned:
                Initial time to resample. It must hold inPosFit < timePoints
            unsigned:
                Final time to resample.
            int:
                Seed to feed the random engine.

        Return:
            struct myMat:
                Structure containing the resampled estimation.
    */

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
        // We iterate from one time to another
        for( unsigned j = inPosFit; j < ouPosFit; j++ ) {
            resampPoints[i*pointsUsed+auxCol] =
                    inFile.m_Matrix[bootChoice*timePoints+j];
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
        auxTime += 1;
    }

    struct myMat Return = { bootstPoints, pointsUsed, inFile.colSize };
    delete resampPoints;    // Delete pointer containing the resampled data
    return Return;
}

double* bootStrap( const myMat inVec, const unsigned numBoots ) {
    /*
        Function to generate the central value of a data using bootstrap.

        Args:
            const struct myMat:
                Input file to be resampled.
            const unsigned:
                Number of bootstrap iterations.

        Return:
            double*:
                Pointer containing the resampled estimation.
    */

    unsigned numPoints = inVec.rowSize;

    double* resampAux = new double[numPoints];
    double* meanBoot = new double[numBoots];

    // INITIALIZE THE RANDOM ENGINE
    typedef std::mt19937 MyRng;
    uint32_t seedVal = std::time(0);
    MyRng rng;
    rng.seed(seedVal);
    std::uniform_int_distribution<int> dist(0,numPoints-1);

    for( unsigned i = 0; i < numBoots; i++ ) {
        for( unsigned j = 0; j < numPoints; j++ )
            resampAux[j] = inVec.m_Matrix[dist(rng)];   // Get a random resample

        meanBoot[i] = meanValue( resampAux, numPoints );
    }

    // RETURN THE ARRAY OF MEANS
    return meanBoot;
}

double meanValue( const double* inPoint, const unsigned numPoints ) {
    /*
       Function that calculates the average of a given vector.

        Args:
            const double*:
                Input pointer from which we would like to generate the average.
            const unsigned:
                Dimension of the vector.

        Return:
            double:
                Average of the vector.
    */

    double retMean = 0.0;
    for( unsigned i = 0; i < numPoints; i++ )
        retMean += inPoint[i];

    return retMean / numPoints;
}

double stdeValue( const double* inPoint, const unsigned numPoints ) {
    /*
       Function that calculates the standard deviation of a given vector.

        Args:
            const double*:
                Input pointer from which we would like to generate the average.
            const unsigned:
                Dimension of the vector.

        Return:
            double:
                Standard deviation of the vector.
    */

    double retStde = 0.0;
    double meanVal = meanValue( inPoint, numPoints );   // We need the mean value
    for( unsigned i = 0; i < numPoints; i++ )
        retStde += ( inPoint[i] - meanVal ) * ( inPoint[i] - meanVal );

    return std::sqrt( retStde / ( numPoints - 1 ) );
}

struct myMat covMat( struct myMat inVal ) {
    /*
       Function that calculates the covariance matrix of a vector, that is,

                    \Sigma_{ij} = \sqrt( y_i * y_j }.

        Args:
            struct myMat:
                Structure containing the data.

        Return:
            struct myMat:
                Structure containing the covariance matrix of the data provided..
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
