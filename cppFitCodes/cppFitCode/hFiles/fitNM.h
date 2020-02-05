//    AUTHOR : SERGIO CHAVES GARCÍA-MASCARAQUE
//    E-MAIL : SERGIOZTESKATE@GMAIL.COM
//    DICIEMBRE DE 2018, SWANSEA, WALES - MADRID, SPAIN


#ifndef SIMPLEX_H
#define SIMPLEX_H

#include <iostream>
#include <iomanip>
#include <random>
#include <ctime>
#include <cmath>

#include "./loadData.h"     // We need to load the matrix structure

// DEFINED FUNCTION NEEDED TO OPTIMIZE THE SIMPLEX
void orderSimplex( double* ,double* ,const unsigned );
void reordMat( double* ,double* ,const unsigned* ,const unsigned );
void centrMat( double*, double*, const unsigned );


struct GetVal {             // Structure to return various values from function
    double* optmPoint;
    double  optmEval;
};

template<class TargetFunc>
double chiSquare( TargetFunc f,
                  double* params,
                  myMat xData,
                  myMat yData,
                  myMat yError,
                  const unsigned matVec) {
    /*
        Function that calculates the chiSquare function to minimize. It uses three
        definitions:

                    chiSquare = \sum_i A * ( yDat_i - f( xDat_i, omega ) ) ^ 2
            * Error function:
                Error in the yAxis are not taken into account. A = 1
            * Maximum likelihood:
                The errors weight each value. A = 1 / diag( yError )
            * Maximum likelihood correlation:
                The chiSquare takes into account the correlation between data.
                A = 1 / corr( yError, yError )

        Args:
            TargetFunc:
                Function to be evaluated in the chiSquare. It must return a
                double and accept a pointer with the paramets and xData.
            double*:
                Parameters omega to evaluate the chiSquare to.
            struct myMat:
                Data in the x axis, it represents xDat
            struct myMat:
                Data in the y axis, it represents yDat, that is, the measurements.
            struct myMat:
                Errors of the measurements. It represents yError.
            const unsigned:
                Flag that tells the problem if you are using a vector or a matrix.
                It is only useful for the maximum likelihood values.

        Return:
            double:
                Value of the chiSquare for the experimental data provided and
                the function given with the parameters of the function being
                omega.

    */

    unsigned rowSize = yData.rowSize;
    unsigned dimParam = yData.colSize;

    double evalChiSq = 0.0;
    double auxVal = 0.0;

    if( matVec == 0  ) {    // Error function
        for( unsigned i = 0; i < rowSize; i++ ) {
            auxVal = f(params,xData.m_Matrix[i]) - yData.m_Matrix[i];
            evalChiSq += std::pow( auxVal, 2 );
        }
    }
    if( matVec == 1 ) { // Maximum likelihood without correlation
        for( unsigned i = 0; i < rowSize; i++ )  {
            auxVal = ( f(params,xData.m_Matrix[i]) - yData.m_Matrix[i] ) /
                    yError.m_Matrix[i];
            evalChiSq += std::pow( auxVal, 2 );
        }
    }
    if( matVec == 2 ) { // Maximum likelihood with correlation
        for( unsigned i = 0; i < rowSize; i++ ) {
            for( unsigned j = 0; j < rowSize; j++ ) {
                auxVal =  ( f(params,xData.m_Matrix[i]) - yData.m_Matrix[i] );
                auxVal *= yError.m_Matrix[i*rowSize+j];
                auxVal *= ( f(params,xData.m_Matrix[j]) - yData.m_Matrix[j] );
                evalChiSq += std::pow( auxVal, 2 );
            }
        }
    }
    return evalChiSq; // / ( rowSize - dimParam ) ;
}


template<class TargetFunc>
struct GetVal fitNM( TargetFunc f,
                     myMat xData,
                     myMat yData,
                     myMat yError,
                     double* initGuess,
                     const double* expGuess,
                     const int dimSpace,
                     int seedAux,
                     const unsigned matVec,
                     const int maxIter = 1e8) {
    /*
        Function that returns the fitted parameters for a given function
        knowing the experimental data. It implements a Nelder-Mead algorithm
        to fit the data.

        Args:
            TargetFunc:
                Target function to fit. It must return a double and be able to
                accpet a pointer with the parameters and a double xData.
            struct myMat:
                Structure containing the points in the x axis.
            struct myMat:
                Structure containing the results of the experiment that we
                would like to fit the function to.
            struct myMat:
                Structure containing the errors in the experiments.
            double*:
                Pointer containing the initial guess for the parameters.
            const double*:
                Pointer containing the volume in which we locate the initial
                simplex.
            const int:
                Dimension of the parameter space.
            int:
                Seed for the random engine.
            const unsigned:
                Flag to check if we are using a matrix or a vector inside the
                chiSquare calculation.
            const int:
                Maximum number of iterations in the calculation.

        Return:
            struct GetVal:
                Structure containing the fitted parameters that minimize the
                chiSquare function. It also returns the value of the chiSquare
                at the produced estimation of the parameters.

    */

    // INITIALIZE THE RANDOM ENGINE
    typedef std::mt19937 MyRng;
    // uint32_t seed_val = seedAux;
    MyRng rng;
    rng.seed(time(0) + seedAux);
    std::uniform_real_distribution<double> dist(-1,1);

    // DEFINE THE CONSTANTS OF THE ALGORITHM
    double alpha = 1.0;
    double gamma = 2.0;
    double rho   = 0.5;
    double sigma = 0.5;

    // DECLARE ALL THE VARIABLES AND POINTERS NEEDED
    unsigned dimSimp = dimSpace + 1;        // Dimension of the simplex
    double noImprThreshold = 1e-8;          // Non improvement threshold
    unsigned noImprovBreak = 100;           // Maximum iterations without improvement

    double lowestValue;                     // Best value achieved
    double refEval;                         // Image of the reflection point
    double expEval;                         // Image of the expansion point
    double conEval;                         // Image of the contraction point

    // AUXILIAR VARIABLES AND POINTERS NEEDED
    double* newPoint = new double[dimSpace]; // Needed to generate the initial simplex
    double* cenPoint = new double[dimSpace]; // Centroid point of the simplex
    double* xRef = new double[dimSpace];     // Reflection point of the simplex
    double* xExp = new double[dimSpace];     // Expansion point of the simplex
    double* xCon = new double[dimSpace];     // Contraction point of the simplex
    double* xRed = new double[dimSpace];     // Reduction point of the simplex

    double valFunc = chiSquare(f, initGuess, xData, yData, yError, matVec );

    // RETURN VARIABLES
    double* optmPoint = new double[dimSimp*dimSpace];

    // INITIATE THE ALGORITHM

    // GENERATE THE INITIAL SIMPLEX USING THE initGuess given
    double* imPoints = new double[dimSimp*dimSpace];
    double* imValues = new double[dimSimp];

    for( unsigned i = 0; i < dimSimp; i++ ) {
        // Substitute initGuess into the simplex
        if( i == 0 ) {
            imValues[i] = valFunc;
            for( unsigned j = 0; j < dimSpace; j++ )
                imPoints[i*dimSpace+j] = initGuess[j];
        }
        // Generate a simplex around initGuess
        else {
            for( unsigned j = 0; j < dimSpace; j++ ) {
                newPoint[j] = initGuess[j] + expGuess[j] * dist(rng);
                imPoints[i*dimSpace+j] = newPoint[j];
            }
            imValues[i] = chiSquare( f, newPoint, xData, yData, yError, matVec );
        }
    }

    unsigned noImprov = 0;
    unsigned iterStep = 0;
    while( 1 ) {

        // CALCULATE THE VALUE OF THE FUNCTION EACH STEP
        for( unsigned i = 0; i <  dimSimp; i++ ) {
            for( unsigned j = 0; j < dimSpace; j++ ) {
                newPoint[j] = imPoints[i*dimSpace+j];
            }
            imValues[i] = chiSquare( f, newPoint, xData, yData, yError, matVec );
        }

        // SORT VALUES USING THE FUNCTIONS CREATED
        orderSimplex( imPoints, imValues, dimSpace );

        // CONTROL THE ALGORITHM

        // Break if maximum iteration value is achieved
        if( maxIter and iterStep >= maxIter ) {
            for( unsigned i = 0; i < dimSpace; i++ ) {
                optmPoint[i] = imPoints[i];
            }
            GetVal r = { optmPoint, imValues[0] };
            return r;
        }
        iterStep += 1;

        // Control the non-improvement
        if( lowestValue < valFunc - noImprThreshold ) {
            noImprov = 0;
            valFunc = lowestValue;
        }
        else
            noImprov += 1;

        // Break if we do not improve results
        if( noImprov >= noImprovBreak ){
            for( unsigned i = 0; i < dimSpace; i++ )
                optmPoint[i] = imPoints[i];
            GetVal r = { optmPoint, imValues[0] };
            return r;
        }

        // CALCULATE THE CENTROID OF THE BEST N POINTS
        centrMat( cenPoint, imPoints, dimSpace );

        // REFLECTION PART OF THE ALGORITHM
        for( unsigned i = 0; i < dimSpace; i++ )
            xRef[i] = cenPoint[i] + alpha * ( cenPoint[i] -
                                imPoints[(dimSimp-1)*dimSpace+i] );
        refEval = chiSquare( f, xRef, xData, yData, yError, matVec );

        if( imValues[0] <= refEval and refEval < imValues[dimSimp-2] ) {
            for( unsigned i = 0; i < dimSpace; i++ )
                imPoints[(dimSimp-1)*dimSpace+i] = xRef[i];
            continue;
        }

        // EXPANSION PART OF THE ALGORITHM
        if( refEval < imValues[0] ) {
            for( unsigned i = 0; i < dimSpace; i++ )
                xExp[i] = cenPoint[i] + gamma * ( xRef[i] - cenPoint[i] );
            expEval = chiSquare( f, xRef, xData, yData, yError, matVec );

            if( expEval < refEval ) {
                for( unsigned i = 0; i < dimSpace; i++ )
                    imPoints[(dimSimp-1)*dimSpace+i] = xExp[i];
                continue;
            }
            else {
                for( unsigned i = 0; i < dimSpace; i++ )
                    imPoints[(dimSimp-1)*dimSpace+i] = xRef[i];
                continue;
            }
        }

        // CONTRACTION PART OF THE ALGORITHM
        for( unsigned i = 0; i < dimSpace; i++ )
            xCon[i] = cenPoint[i] + rho * ( imPoints[(dimSimp-1)*dimSpace+i] -
                      cenPoint[i] );
        conEval = chiSquare( f, xCon, xData, yData, yError, matVec );

        if( conEval < imValues[dimSimp-1] ) {
            for( unsigned i = 0; i < dimSpace; i++ )
                imPoints[(dimSimp-1)*dimSpace+i] = xCon[i];
            continue;
        }

        // SHRINK PART OF THE ALGORITHM
        for( unsigned i = 1; i < dimSimp; i++ ) {
            for( unsigned j = 0; j < dimSpace; j++ ) {
                imPoints[i*dimSpace+j] = imPoints[j] +
                            sigma * ( imPoints[i*dimSpace+j] - imPoints[j] );
            }
        }
    }   // END OF WHILE

} // END OF FUNCTION

void centrMat( double* cenPoint, double* inMat, const unsigned dimSpace ) {
    /*
       Function to generate the centroid of the simplex. The centroid is
       defined as the average of all the components of a column.

        Args:
            double*:
                Pointer that will contain the centroid of the simplex.
            double*:
                Input simplex coordinates to calculate the centroid.
            const unsigned:
                Dimension of the parameter space.
    */

    int dimSimp = dimSpace + 1;
    double centVal = 0.0;

    for( unsigned i = 0; i < dimSpace; i++ ) {
        for( unsigned j = 0; j < dimSimp - 1; j++ ) {
            centVal += inMat[j*dimSpace+i];
        }
        cenPoint[i] = centVal / ( dimSimp - 1 );
        centVal = 0.0;
    }
}

void reordMat( double* inMat, double* auxMat,
               const unsigned* auxColoc, const unsigned dimSpace ) {
    /*
        Function to reorder a matrix given a vector containing the order
        specified by row.

        Args:
            double*:
                Input matrix to be reordered.
            double*:
                Auxiliar matrix used to reorder the input matrix.
            const unsigned*:
                Array containing the new locations of each row.
            const unsigned:
                Dimension of the matrix, that is, number of rows of the
                matrix to be reordered.
    */

    int dimSimp = dimSpace + 1;
    for( unsigned i = 0; i < dimSimp; i++ ) {
        for( unsigned j = 0; j < dimSpace; j++ ) {
            auxMat[i*dimSpace+j] = inMat[auxColoc[i]*dimSpace+j];
        }
    }

    for( unsigned i = 0; i < dimSimp * dimSpace; i++ )
        inMat[i] = auxMat[i];
}

void orderSimplex( double* imPoints, double* evSimp, const unsigned dimSpace ) {
    /*
        Function to reorder a matrix given a metric. The sorting is done in
        ascendent order.

        Args:
            double*:
                Input pointer to be sorted ascendently.
            double*:
                Input pointer containing the metric that dictates the order
                to sort the input pointer.
            const unsigned:
                Dimension of the matrix, that is, number of rows of the
                matrix to be reordered.
    */

    int dimSimp = dimSpace + 1;

    unsigned* auxColoc = new unsigned[dimSimp];
    // Collocation auxiliar pointer
    for( unsigned i = 0; i < dimSimp; i++ )
        auxColoc[i] = i;

    // Ascendant order of evSimp
    double tempVarA;
    unsigned tempVarB;
    for( unsigned i = 0; i < dimSimp; i++ ) {
        for( unsigned j = 0; j < dimSimp - i - 1; j++ ) {
            if( evSimp[j] > evSimp[j+1] ) {
                // Order evSimp
                tempVarA = evSimp[j];
                evSimp[j] = evSimp[j+1];
                evSimp[j+1] = tempVarA;

                // The same for the initial collocation
                tempVarB = auxColoc[j];
                auxColoc[j] = auxColoc[j+1];
                auxColoc[j+1] = tempVarB;
            }
        }
    }

    double* auxMat = new double[dimSimp*dimSpace];
    reordMat( imPoints, auxMat, auxColoc, dimSpace );
}

#endif
