#include <iostream>
#include <string>
#include <ctime>
#include "mpi.h"

// Include header files
#include "./hFiles/loadData.h"
#include "./hFiles/readInput.h"
#include "./hFiles/bootstrap.h"

#include "./hFiles/obtBest.h"
#include "./hFiles/pcacMass.h"

// Include headers used in the main file
#include <random>

int main( int argc, char* argv[] ) {

    // Initialize the MPI engine
    MPI_Init( NULL, NULL );

    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );     // Get the rank of each process
    int auxSize;
    MPI_Comm_size( MPI_COMM_WORLD, &auxSize );  // Get the size of the calculation
    unsigned size = auxSize;

    // All ranks must have this data ---
    std::string fileAP = getString( argv[1], "fileAP" );
    std::string filePP = getString( argv[1], "filePP" );
    std::string fileAA = getString( argv[1], "fileAA" );

    unsigned rows;              // Number of rows contained in the files
    unsigned cols;              // Number of columns contained in the files
    unsigned timePoints;        // Number of time points contained in the simulation
    unsigned numBootstrap;      // Number of bootstrap iteration
    unsigned initSeed;          // Initial seed of the calculation
    unsigned COLRESCALE;        // Rescale or not the column selected
    unsigned selectCol;         // Select column to rescale
    double resFactor;           // Rescaling factor of the data

    if( rank == 0 ) {

        // Get values from input file
        unsigned* params = getVecUnsigned( argv[1], "params", 3 );

        rows = params[0];         // Number of rows in the PP channel
        cols = params[1];         // Number of cols in the AA channel
        timePoints = params[2];   // Number of points in the time direction

        // Number of bootstrap iterations in the calculation
        numBootstrap = getUnsigned( argv[1], "numBootstrap" );

        // Initial seed of the calculation
        initSeed = getUnsigned( argv[1], "initSeed" );

        // Rescale the data if needed
        COLRESCALE = getUnsigned( argv[1], "rescaleOrNot" );
        selectCol = getUnsigned( argv[1], "rescaledColumn" );

        double* valueRescale = getDouble( argv[1], "rescaledFactor", 1 );
        resFactor = valueRescale[0];
    }

    // Share the data from master rank
    MPI_Bcast( &rows, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    MPI_Bcast( &cols, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    MPI_Bcast( &timePoints, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    MPI_Bcast( &numBootstrap, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    MPI_Bcast( &COLRESCALE, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    MPI_Bcast( &selectCol, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    MPI_Bcast( &resFactor, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );

    // ----- ATTENTION -----
    MPI_Barrier( MPI_COMM_WORLD );

    // Each rank will calculate a portion of the bootstrapped data
    struct myMat dataFile, yData;

    numBootstrap = numBootstrap / size;

    // Then we calculate AA - g0g5 - g0g5 -------------------------------------------
    std::string saveAA = getString( argv[1], "saveAA" );

    dataFile = getLargeData( fileAA.c_str(), rows, cols,
                             COLRESCALE, selectCol, resFactor );

    yData = sliceCol( dataFile, 1 );

    obtBest( yData, numBootstrap, timePoints, initSeed, saveAA, rank, size );

    // ----- ATTENTION -----
    MPI_Barrier( MPI_COMM_WORLD );

    // First we calculate PP -  g5 - g5 ---------------------------------------------
    std::string savePP = getString( argv[1], "savePP" );

    dataFile = getLargeData( filePP.c_str(), rows, cols,
                             COLRESCALE, selectCol, resFactor );

    yData = sliceCol( dataFile, 1 );

    obtBest( yData, numBootstrap, timePoints, initSeed, savePP, rank, size );


    // ----- ATTENTION -----
    MPI_Barrier( MPI_COMM_WORLD );

    // Finally, we calculate the PCAC mass ------------------------------------------
    std::string savePCAC = getString( argv[1], "savePCAC" );

    struct myMat dataFileAP, yDataAP;

    dataFileAP = getLargeData( fileAP.c_str(), rows, cols,
                               COLRESCALE, selectCol, resFactor );

    yDataAP = sliceCol( dataFileAP, 1 );

    obtPCAC( yDataAP, yData, numBootstrap, timePoints, initSeed,
             savePCAC, rank, size );

    // Finalize the MPI communicator
    MPI_Finalize();

return 0;
}
