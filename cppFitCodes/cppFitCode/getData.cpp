#include <iostream>
#include <iomanip>
#include <string>
#include <stdexcept>
#include <ctime>
#include "mpi.h"

// Include header files
#include "./hFiles/loadData.h"
#include "./hFiles/readInput.h"
#include "./hFiles/bootstrap.h"
#include "./hFiles/fitNM.h"
#include "./hFiles/linalg.h"

// Include headers used in the main file
#include <random>
#include <ctime>

int main( int argc, char* argv[] ) {

    // Initialize MPI engine
    // Initialize the MPI engine
    MPI_Init( NULL, NULL );

    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );     // Get the rank of each process
    int auxSize;
    MPI_Comm_size( MPI_COMM_WORLD, &auxSize );  // Get the size of the calculation
    unsigned size = auxSize;

    std::string fileName = getString( argv[1], "fileName" );    // FileName to fit

    unsigned rowSize;           // Number of columns contained in the dataFile
    unsigned colSize;           // Number of rows contained in the dataFile
    unsigned timePoints;        // Number of time points contained in the simulation
    unsigned numBootstrap;      // Number of bootstrap iteration
    unsigned initSeed;          // Initial seed of the calculation
    unsigned COLRESCALE;        // Rescale or not the column selected
    unsigned selectCol;         // Select column to rescale
    unsigned binSize;           // Binning size used
    double resFactor;           // Rescaling factor of the data

    if( rank == 0 ) {

        // Get values from input file
        unsigned* paramsFile = getVecUnsigned( argv[1], "paramsFile", 3 );
        rowSize = paramsFile[0];    // Number of rows in the file
        colSize = paramsFile[1];    // Number of cols in the file
        timePoints = paramsFile[2]; // Number of timePoints in each file

        // Number of bootstrap iterations in the calculation
        numBootstrap = getUnsigned( argv[1], "numBootstrap" );

        // Binning size of the calculation
        binSize = getUnsigned( argv[1], "binSize" );

        // Initial seed of the calculation
        initSeed = getUnsigned( argv[1], "initSeed" );

        // Rescale the data if needed
        COLRESCALE = getUnsigned( argv[1], "rescaleOrNot" );
        selectCol = getUnsigned( argv[1], "rescaledColumn" );

        double* valueRescale = getDouble( argv[1], "rescaledFactor", 1 );
        resFactor = valueRescale[0];
    }

    // Share the data from master rank
    MPI_Bcast( &rowSize, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    MPI_Bcast( &colSize, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    MPI_Bcast( &timePoints, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    MPI_Bcast( &numBootstrap, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    MPI_Bcast( &binSize, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    MPI_Bcast( &COLRESCALE, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    MPI_Bcast( &selectCol, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    MPI_Bcast( &resFactor, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Barrier( MPI_COMM_WORLD );

    struct myMat myData = getLargeData( fileName.c_str(), rowSize, colSize,
                                        COLRESCALE, selectCol, resFactor );

    struct myMat yData = sliceCol( myData, 1 );

    // Now we bin the data
    int numComp = yData.rowSize / ( binSize * timePoints );
    double* binnedPoints = new double[numComp*timePoints];
    double holdBin = 0;

    for( unsigned i = 0; i < numComp; i++ ) {
        for( unsigned j = 0; j < timePoints; j++ ) {
            for( unsigned k = 0; k < binSize; k++ ) {
                holdBin += yData.m_Matrix[i*binSize*timePoints+k*timePoints+j];
            }
            binnedPoints[i*timePoints+j] = holdBin / binSize;
            holdBin = 0;
        }
    }

    // Hold the new data inside the yData
    yData.m_Matrix = binnedPoints;
    yData.rowSize = numComp * timePoints;

    numBootstrap = numBootstrap / size;

    // GENERATE RANDOM SEEDS IN MASTER RANK
    int* bootsSeed = new int[numBootstrap];
    if( rank == 0 ) {
        typedef std::mt19937 MyRng;
        uint32_t seed_val = initSeed;
        MyRng rng;
        rng.seed(seed_val);
        std::uniform_int_distribution<int> dist(0,1e8);

        for( unsigned i = 0; i < numBootstrap; i++ )
            bootsSeed[i] = dist(rng);

    }

    // CALCULATE A BOOTSTRAP SAMPLE OF THE WHOLE DATA FILE
    double* holdFile = new double[numBootstrap*timePoints];
    for( unsigned i = 0; i < numBootstrap; i++ ) {

        // GENERATE RANDOM SEEDS IN MASTER RANK
        int* seedStorage = new int[size];
        if( rank == 0 ) {
            typedef std::mt19937 MyRng;
            uint32_t seed_val = bootsSeed[i];
            MyRng rng;
            rng.seed(seed_val);
            std::uniform_int_distribution<int> dist(0,1e8);

            for( unsigned i = 0; i < size; i++ )
                seedStorage[i] = dist(rng);
        }

        int bufferSeed;
        MPI_Barrier( MPI_COMM_WORLD );
        MPI_Scatter( seedStorage, 1, MPI_INT, &bufferSeed, 1,
                     MPI_INT, 0, MPI_COMM_WORLD );
        MPI_Barrier( MPI_COMM_WORLD );
        delete seedStorage;

        // Calculate bootstrap resample of the whole data file

        struct myMat sampleFile = pickAnalyzeBootstrap( yData, timePoints,
                                                        0, timePoints, bufferSeed );

        // Hold the data into a pointer
        for( unsigned j = 0; j < timePoints; j++ ) {
            holdFile[i*timePoints+j] = sampleFile.m_Matrix[j];
        }
        // delete sampleFile.m_Matrix;
    }
    // This pointer will contain all the data
    double* gatherFile = new double[size*numBootstrap*timePoints];

    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Gather( holdFile, numBootstrap * timePoints, MPI_DOUBLE,
                gatherFile, numBootstrap * timePoints, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Barrier( MPI_COMM_WORLD );
    // delete holdFile;

    // Now we calculate the average over the bootstrapped statistics in the MASTER rank
    if( rank == 0 ) {

        struct myMat meanFile = { gatherFile, size * numBootstrap, timePoints };
        struct myMat colData;

        double* fileMean = new double[timePoints];
        double* fileStde = new double[timePoints];

        for( unsigned i = 0; i < timePoints; i++ ) {
            colData = sliceCol( meanFile, i );
            fileMean[i] = meanValue( colData.m_Matrix, size * numBootstrap );
            fileStde[i] = stdeValue( colData.m_Matrix, size * numBootstrap );
        }

        std::string saveName = getString( argv[1], "saveName" );

        // Write down results into a file
        std::ofstream saveStream;
        saveStream.open( saveName );
        saveStream << "# Time    <Corr>    sigma(Corr)" << std::endl;
        for( unsigned i = 0; i < timePoints; i++ )
            saveStream << i << "\t" << fileMean[i] << "\t" << fileStde[i] << std::endl;

        saveStream.close();
    }

    delete gatherFile;
    MPI_Finalize();

return 0;
}
