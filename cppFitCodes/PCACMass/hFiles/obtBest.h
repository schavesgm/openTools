#ifndef OBTBEST_H
#define OBTBEST_H

#include <string>
#include "mpi.h"

// Include header files
#include "./bootstrap.h"

void obtBest( const struct myMat yData, const unsigned numBootstrap,
              const unsigned timePoints, const unsigned initSeed,
              std::string saveName, int rank, int size ) {

    // Random engine to use in the calculation
    std::mt19937 rngEng;

    // Buffers needed
    double* gatherFile = new double[size * numBootstrap * timePoints];  // Whole data
    double* buffResamp = new double[numBootstrap * timePoints];         // Rank data

    int* seedsBoots =  new int[numBootstrap];                           // Boots seeds
    int* seedRank = new int[size];                                      // Rank seeds

    struct myMat sampleFile;                                            // Hold result

    // We need new seeds for each bootstrap iteration
    if ( rank == 0 ) {  // Generate in master, then scatter
        // uint32_t seed = initSeed;
        rngEng.seed( initSeed );
        std::uniform_int_distribution<int> dist( 0, 1e8 );

        for ( unsigned nB = 0; nB < numBootstrap; nB++ ) {
            seedsBoots[nB] = dist( rngEng );
        }
    }

    // Calculate the bootstrap resample
    for ( unsigned nB = 0; nB < numBootstrap; nB++ ) {

        // Generate random seed for each rank
        if ( rank == 0 ) {
            rngEng.seed( seedsBoots[nB] );
            std::uniform_int_distribution<int> dist( 0, 1e8 );

            for ( unsigned i = 0; i < size; i++ ) {
                seedRank[i] = dist( rngEng );
            }
        }

        int buffSeed;
        MPI_Barrier( MPI_COMM_WORLD );
        MPI_Scatter( seedRank, 1, MPI_INT, &buffSeed, 1, MPI_INT, 0, MPI_COMM_WORLD );
        MPI_Barrier( MPI_COMM_WORLD );

        // Calculate bootstrap resample of the whole data file
        sampleFile = bootstrap( yData, timePoints, 0, \
                                timePoints, buffSeed );

        // Hold the data into the buffer
        for ( unsigned t = 0; t < timePoints; t++ ) {
            buffResamp[t+nB*timePoints] = sampleFile.m_Matrix[t];
        }
    }

    // Now, we need to gather all the data into a single rank - Master
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Gather( buffResamp, numBootstrap * timePoints, MPI_DOUBLE,
                gatherFile, numBootstrap * timePoints, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Barrier( MPI_COMM_WORLD );

    // Now we calculate the average over the bootstrapped statistics in the MASTER rank
    if( rank == 0 ) {

        struct myMat meanFile = { gatherFile, size * numBootstrap, timePoints };
        struct myMat colData;

        double* avgRes = new double[timePoints];
        double* stdRes = new double[timePoints];

        for( unsigned i = 0; i < timePoints; i++ ) {
            colData = sliceCol( meanFile, i );
            avgRes[i] = meanValue( colData.m_Matrix, size * numBootstrap );
            stdRes[i] = stdeValue( colData.m_Matrix, size * numBootstrap );
        }

        // Write down results into a file
        std::ofstream saveStream;
        saveStream.open( saveName );
        saveStream << "# Time    <Corr>    sigma(Corr)" << std::endl;
        for( unsigned i = 0; i < timePoints; i++ )
            saveStream << i << "\t" << avgRes[i] << "\t" << stdRes[i] << std::endl;
        saveStream.close();
    }
}       // End of function

#endif
