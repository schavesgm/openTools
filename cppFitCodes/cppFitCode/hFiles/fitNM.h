//------------------------------------------------------------------------------------------------//
//    AUTHOR : SERGIO CHAVES GARCÍA-MASCARAQUE
//    E-MAIL : SERGIOZTESKATE@GMAIL.COM
//    DICIEMBRE DE 2018, SWANSEA, WALES - MADRID, SPAIN
//------------------------------------------------------------------------------------------------//


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

//
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
     * Function that calculates the chiSquare function to minimize. It uses two definitions.
     * The error function, where errors in the y Axis are not taken into account and the maximum
     * likelihood definition, where the errors weight each value.
     *
     * TargetFunc f : target function to fit, must return a double and be able to accept a pointer
     *                with the parameters and a double xData
     * myMat xData  : object in myMat structure containing the points in the x Axis, must have 
     *                the number of rows included in the definition.
     * myMat yData  : object in myMat structure containing the points measured in the y Axis. Must
     *                have the number of rows included in the definition. The element colSize 
     *                contained in the definition of yData must have the dimension of the parameter
     *                space that we would like to fit to calculate the degrees of freedom of the 
     *                problem.
     * myMat yError : object in myMat structure containing the errors in the measurements. In the 
     *                case of number of columns given == 0 it selects error function to minimize.
     *                If the error provided is the covariance matrix it takes into account 
     *                correlations in the data.
     * unsigned     : flag that notifies when you are using a vector or a matrix. matVec = 0 means
     *                that you are passing a vector, matVec = 1 means that you are passing a matrix
     *
     * return double: value of the function at the points given
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
            auxVal = ( f(params,xData.m_Matrix[i]) - yData.m_Matrix[i] ) / yError.m_Matrix[i]; 
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
     * TargetFunc f : target function to fit, must return a double and be able to accept a pointer
     *                with the parameters and a double xData
     * myMat xData  : object in myMat structure containing the points in the x Axis, must have 
     *                the number of rows included in the definition
     * myMat yData  : object in myMat structure containing the points measured in the y Axis. Must
     *                have the number of rows included in the definition
     * myMat yError : object in myMat structure containing the errors in the measurements. In the 
     *                case of number of columns given == 0 it selects error function to minimize
     * double* iniGuess : pointer to double array that contains initial guess to the optimization
     * double* expGuess : pointer to double array containing the 'error' in each parameter space.
     *                     This means that we expect the value to be around that value
     * const int dimSpace  : integer containing the dimension of the parameter space
     * const int maxIter  : maximum number of iterations allowed in the calculation
     *
     * return struct    : Return a structure containing the double pointer with the points obtained
     *                    in the optimization and a doube containing the value of the function 
     *                    at the most optimized points.
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
    unsigned dimSimp = dimSpace + 1;             // Dimension of the simplex 
    double noImprThreshold = 1e-8;               // Non improvement threshold 
    unsigned noImprovBreak = 100;                // Maximum iterations without improvement
    
    double lowestValue;                     // Best value achieved
    double refEval;                         // Image of the reflection point
    double expEval;                         // Image of the expansion point
    double conEval;                         // Image of the contraction point

    // AUXILIAR VARIABLES AND POINTERS NEEDED
    double* newPoint = new double[dimSpace];        // Needed to generate the initial simplex
    double* cenPoint = new double[dimSpace];        // Centroid point of the simplex
    double* xRef = new double[dimSpace];            // Reflection point of the simplex
    double* xExp = new double[dimSpace];            // Expansion point of the simplex
    double* xCon = new double[dimSpace];            // Contraction point of the simplex
    double* xRed = new double[dimSpace];            // Reduction point of the simplex

    double valFunc = chiSquare(f, initGuess, xData, yData, yError, matVec );

    // RETURN VARIABLES
    double* optmPoint = new double[dimSimp*dimSpace];

    // INITIATE THE ALGORITHM
    
    // GENERATE THE INITIAL SIMPLEX USING THE initGuess given
    double* imPoints = new double[dimSimp*dimSpace];
    double* imValues = new double[dimSimp];

    for( unsigned i = 0; i < dimSimp; i++ ) { 
        if( i == 0 ) {                              // Substitute initGuess into the simplex
            imValues[i] = valFunc;
            for( unsigned j = 0; j < dimSpace; j++ ) 
                imPoints[i*dimSpace+j] = initGuess[j];
        }
        else {                                      // Generate a simplex around initGuess
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
        if( maxIter and iterStep >= maxIter ) {     // Break if maximum iteration value is achieved
            for( unsigned i = 0; i < dimSpace; i++ ) {
                optmPoint[i] = imPoints[i];
            }
            GetVal r = { optmPoint, imValues[0] };
            // std::cout << "Number of iterations: " << iterStep << std::endl;
            return r;
        }
        iterStep += 1;
        
        if( lowestValue < valFunc - noImprThreshold ) {     // Control the non-improvement
            noImprov = 0;
            valFunc = lowestValue;
        }
        else 
            noImprov += 1;
        
        if( noImprov >= noImprovBreak ){                    // Break if we do not improve results
            for( unsigned i = 0; i < dimSpace; i++ )
                optmPoint[i] = imPoints[i];
            GetVal r = { optmPoint, imValues[0] };
            // std::cout << "Number of iterations: " << iterStep << std::endl;
            return r;
        }
        
        // CALCULATE THE CENTROID OF THE BEST N POINTS
        centrMat( cenPoint, imPoints, dimSpace );
        
        // REFLECTION PART OF THE ALGORITHM
        for( unsigned i = 0; i < dimSpace; i++ ) 
            xRef[i] = cenPoint[i] + alpha * ( cenPoint[i] - imPoints[(dimSimp-1)*dimSpace+i] );
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
            xCon[i] = cenPoint[i] + rho * ( imPoints[(dimSimp-1)*dimSpace+i] - cenPoint[i] );
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
     * double* cenPoint : pointer to double array that will contain the centroid of the matrix
     * double* inMat    : input matrix to calculate the centroid. It calculates the mean value of
     *                    each column averaging over all the rows at fixed column.
     * int dimSpace     : integer containing the dimension of the parameter space
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

void reordMat( double* inMat, double* auxMat, const unsigned* auxColoc, const unsigned dimSpace ) {
    /* 
     * double* inMat      : pointer to array that contains the input matrix to reorder
     * double* auxMat     : auxiliar matrix used to reorder the input matrix
     * unsigned* auxColoc : array containing the new locations of each row
     * int dimSpace       : integer containing the dimension of the parameter space
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
     * double* imPoints   : input pointer containing the matrix that we want to sort ascendent
     * double* evSimp     : input pointer containing the evaluation of the simplex in the function
     * int dimSpace       : integer containing the dimension of the parameter space
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
