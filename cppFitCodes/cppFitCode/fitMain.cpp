#include <iostream>
#include <iomanip>
#include <string>
#include <stdexcept>
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

#define ERR_FUNCTION 0
#define MAX_LIKELIHOOD 1
#define DO_NOT_CALCULATE 0
#define SIMPLE_PRINT 0
#define TIME_EXTENT 128
#define ANSATZ 0

// DEFINE HERE THE FUNCTION YOU WOULD LIKE TO FIT
#if ANSATZ == 0     // Cosh function
double f( double *params, double xData ) {
    return params[0] * ( std::exp( - params[1] * xData ) + \
                         std::exp( - params[1] * ( TIME_EXTENT - xData ) ) );

}
#elif ANSATZ == 1   // Exponential function
double f( double *params, double xData ) {
    return params[0] * ( std::exp( - params[1] * xData ) );
}
#endif

// -----------------------------------------------------------------------------------------------//

int main( int argc, char* argv[] ) {

    // Initialize MPI engine
    // Initialize the MPI engine
    MPI_Init( NULL, NULL );

    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );     // Get the rank of each process
    int auxSize;
    MPI_Comm_size( MPI_COMM_WORLD, &auxSize );  // Get the size of the calculation
    unsigned size = auxSize;

    // READ DATA FROM FILE INTO THE MASTER RANK --------------------//

    std::string fileName = getString( argv[1], "fileName" );    // Fit file
    std::string saveName = getString( argv[1], "saveName" );    // Save file

    unsigned rowSize;           // Number of columns contained in dataFile
    unsigned colSize;           // Number of rows contained in the dataFile
    unsigned timePoints;        // Number of time points
    unsigned inPosFit;          // Initital position of the fit in time
    unsigned ouPosFit;          // Final position of the fit in time
    unsigned dimParam;          // Dimension of the parameter space to fit
    double* initGuess;          // Initial guess for the fit
    double* explGuess;          // Initial exploring volume for the fit
    unsigned numBootstrap;      // Number of bootstrap iteration
    unsigned binSize;           // Binning size used in the calculation
    unsigned initSeed;          // Initial seed of the calculation
    unsigned calculateBest;     // Calculate best estimator or not for the data
    unsigned COLRESCALE;        // Rescale or not the column selected
    unsigned selectCol;         // Select which column to rescale
    double resFactor;           // Rescaling factor for the data
    double largestChiSquare;    // Largest chi-Square accepted
    unsigned defChiSquare;      // Definition of the chiSquare used.

    if( rank == 0 ) {

        // Get values from input file

        unsigned* paramsFile = getVecUnsigned( argv[1], "paramsFile", 3 );
        rowSize = paramsFile[0];    // Number of rows in the file
        colSize = paramsFile[1];    // Number of cols in the file
        timePoints = paramsFile[2]; // Number of timePoints in each file

        unsigned* winFit = getVecUnsigned( argv[1], "windowFit", 2 );
        inPosFit = winFit[0];   // Initial position of the fit in time direction
        ouPosFit = winFit[1];   // Final position of the fit in time direction

        if( ouPosFit <= inPosFit )  // Keep some consistency in the calculation
            ouPosFit = timePoints;

        // Get some properties of the fit
        dimParam = getUnsigned( argv[1], "dimParam" );
        initGuess = getDouble( argv[1], "initGuess", 2 );
        explGuess = getDouble( argv[1], "explGuess", 2 );

        // Number of bootstrap iterations in the calculation
        numBootstrap = getUnsigned( argv[1], "numBootstrap" );

        // Binning size of the calculation
        binSize = getUnsigned( argv[1], "binSize" );

        // Initial seed of the calculation
        initSeed = getUnsigned( argv[1], "initSeed" );

        // Calculate best estimator for the data or not
        calculateBest = getUnsigned( argv[1], "calculateBest" );

        // Rescale the data if needed
        COLRESCALE = getUnsigned( argv[1], "rescaleOrNot" );
        selectCol = getUnsigned( argv[1], "rescaledColumn" );
        double* valueRescale = getDouble( argv[1], "rescaledFactor", 1 );
        resFactor = valueRescale[0];

        // Get the largest chi-Square value accepted
        double* getLargest = getDouble( argv[1], "largestChiSquare", 1 );
        largestChiSquare = getLargest[0];

        // Get the definition of the chiSquare function used
        defChiSquare = getUnsigned( argv[1], "defChiSquare" );

    }

    // Share the data from master rank
    MPI_Bcast( &rowSize, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    MPI_Bcast( &colSize, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    MPI_Bcast( &timePoints, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    MPI_Bcast( &inPosFit, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    MPI_Bcast( &ouPosFit, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    MPI_Bcast( &dimParam, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    MPI_Bcast( &numBootstrap, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    MPI_Bcast( &binSize, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    MPI_Bcast( &calculateBest, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    MPI_Bcast( &COLRESCALE, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    MPI_Bcast( &selectCol, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    MPI_Bcast( &resFactor, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &defChiSquare, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    MPI_Barrier( MPI_COMM_WORLD );

    unsigned pointsUsed = ouPosFit - inPosFit;

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
    yData.colSize = dimParam;

    // Get the properties of the fit, initial guesses and exploring volume
    if( rank != 0 ) {
        initGuess = new double[dimParam];
        explGuess = new double[dimParam];
    }

    // Share the properties of the fit to each rank
    MPI_Bcast( initGuess, dimParam, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( explGuess, dimParam, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Barrier( MPI_COMM_WORLD );

    int* bootsSeed = new int[numBootstrap];     // Pointer that contains the data
    numBootstrap = numBootstrap / size;         // Each rank will calculate a portion

    // GENERATE THE FILE WITH THE BEST ESTIMATION OF THE OBSERVABLE USING BOOTSTRAP
    if( calculateBest != DO_NOT_CALCULATE ) {

        std::cout << "I am calculating the best estimator" << std::endl;

        // GENERATE RANDOM SEEDS IN MASTER RANK
        if( rank == 0 ) {
            typedef std::mt19937 MyRng;
            uint32_t seed_val = std::time(0);
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
            delete sampleFile.m_Matrix;
        }
        // This pointer will contain all the data
        double* gatherFile = new double[size*numBootstrap*timePoints];

        MPI_Barrier( MPI_COMM_WORLD );
        MPI_Gather( holdFile, numBootstrap * timePoints, MPI_DOUBLE,
                    gatherFile, numBootstrap * timePoints, MPI_DOUBLE,
                    0, MPI_COMM_WORLD );
        MPI_Barrier( MPI_COMM_WORLD );
        delete holdFile;

        // Calculate average of bootstrap resample
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


            // Write down results into a file
            std::ofstream saveStream;
            saveStream.open( saveName );
            saveStream << "# Time    <Corr>    sigma(Corr)" << std::endl;
            for( unsigned i = 0; i < timePoints; i++ )
                saveStream << i << "\t" << fileMean[i] << "\t" << fileStde[i] << std::endl;

            saveStream.close();
        }
        delete gatherFile;
    }

    // FIT THE DATA

    if( rank == 0 ) {
        typedef std::mt19937 MyRng;
        uint32_t seed_val = initSeed;
        MyRng rng;
        rng.seed(seed_val);
        std::uniform_int_distribution<int> dist(0,1e8);

        for( unsigned i = 0; i < numBootstrap; i++ )
            bootsSeed[i] = dist(rng);
    }

    // Generate xData, which will be the time points
    double* xPoints = new double[pointsUsed];
    // Fill the xData structure with the times used in the calculation
    unsigned auxPos = 0;
    for( int i = inPosFit; i < ouPosFit; i++ ) {
        xPoints[auxPos] = i;
        auxPos += 1;
    }
    auxPos = 0;
    struct myMat xData = { xPoints, pointsUsed, 1 };

    // Generate a error vector with the bootstrapped data
    struct myMat yError = loadData( saveName.c_str() );
    yError = sliceCol( yError, 2 );
    double* yErrorContainer = new double[pointsUsed];
    for( unsigned i = inPosFit; i < ouPosFit; i++ ) {
        yErrorContainer[auxPos] = yError.m_Matrix[i];
        auxPos += 1;
    }
    yError = { yErrorContainer, pointsUsed, 1 };    // Error contained in the file

    // Now we proceed with the fitting calculation
    double* holdValue = new double[numBootstrap*dimParam];
    double* holdChiSquare = new double[numBootstrap];

    for( unsigned i = 0; i < numBootstrap; i++ ) {

        // GENERATE RANDOM SEEDS IN MASTER RANK
        int* seedStorage = new int[size];
        if( rank == 0 ) {
            typedef std::mt19937 MyRng;
            uint32_t seed_val = bootsSeed[i];
            MyRng rng;
            rng.seed(seed_val);
            std::uniform_int_distribution<int> dist(0,1e8);

            for( unsigned i = 0; i < size; i++ )  {
                seedStorage[i] = dist(rng);
            }
        }

        int bufferSeed;
        MPI_Barrier( MPI_COMM_WORLD );
        MPI_Scatter( seedStorage, 1, MPI_INT, &bufferSeed, 1, MPI_INT, 0, MPI_COMM_WORLD );
        MPI_Barrier( MPI_COMM_WORLD );
        delete seedStorage;

        // Calculate the bootstrap resample of the data at each rank with different
        struct myMat sampleFit = pickAnalyzeBootstrap( yData, timePoints,
                                                inPosFit, ouPosFit, bufferSeed );

        // Fit the data using the resample at each rank
        struct GetVal getFit = fitNM( f, xData, sampleFit, yError,
                                      initGuess, explGuess, dimParam,
                                      bufferSeed * rank, defChiSquare );
        // Hold the data into a pointer
        for( unsigned j = 0; j < dimParam; j++ ) {
            holdValue[i*dimParam+j] = getFit.optmPoint[j];
        }
        holdChiSquare[i] = getFit.optmEval;

    }

    double* gatherFitValue = new double[size*numBootstrap*dimParam];
    double* gatherChiValue = new double[size*numBootstrap];

    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Gather( holdValue, numBootstrap * dimParam, MPI_DOUBLE,
                gatherFitValue, numBootstrap * dimParam, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Gather( holdChiSquare, numBootstrap, MPI_DOUBLE, gatherChiValue,
                numBootstrap, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Barrier( MPI_COMM_WORLD );

    // Now we calculate the average over the bootstrapped statistics in the MASTER rank
    if( rank == 0 ) {

        // Check for chi-square value and delete largest chi-square fits
        int* controlChi = new int[size*numBootstrap];     // One for each chi-Square
        unsigned countChecker = 0;
        for( unsigned i = 0; i < size * numBootstrap; i++ ) {
            if( gatherChiValue[i] / ( pointsUsed - dimParam )  > largestChiSquare ) {
                controlChi[i] = 1;
            }
            else {
                controlChi[i] = 0;
                countChecker += 1;
            }
        }

        // Get print selection
        unsigned printType = getUnsigned( argv[1], "printFancy" );

        if( printType != SIMPLE_PRINT )
            std::cout << "We got " << countChecker << " measurements from  " <<
                         size * numBootstrap <<
                         " initial calculations that obey reduced chiSquare < " <<
                         largestChiSquare << std::endl;

        // Eliminate those points from the final result to improve the calculation
        double* checkedValues = new double[countChecker*dimParam];
        double* checkedChi = new double[countChecker];
        int checkAux = 0;
        for( unsigned i = 0; i < size * numBootstrap; i++ ) {
            if( controlChi[i] == 0 ) {
                checkedChi[checkAux] = gatherChiValue[i] / ( pointsUsed - dimParam );
                for( unsigned j =0; j < dimParam; j++ ) {
                    checkedValues[checkAux*dimParam+j] = gatherFitValue[i*dimParam+j];
                }
                checkAux += 1;
            }
        }

        // Get the fitted parameters
        struct myMat meanData = { checkedValues, countChecker, dimParam };
        struct myMat colData;

        double* finalMean = new double[dimParam];
        double* finalStde = new double[dimParam];

        for( unsigned i = 0; i < dimParam; i++ ) {
            colData = sliceCol( meanData, i );
            finalMean[i] = meanValue( colData.m_Matrix, countChecker );
            finalStde[i] = stdeValue( colData.m_Matrix, countChecker );
        }

        // Get the chiSquare value
        double chiValue = 0.0;
        for( unsigned i = 0; i < countChecker; i++ )
            chiValue += checkedChi[i];
        chiValue = chiValue / ( countChecker );

        double chiError = 0.0;
        for( unsigned i = 0; i < countChecker; i++ )
            chiError += ( checkedChi[i] - chiValue ) * ( checkedChi[i] - chiValue );
        chiError = std::sqrt( chiError / ( countChecker - 1 ) );

        if( printType != SIMPLE_PRINT ) {
            for( unsigned i = 0; i < dimParam; i++ )
                std::cout << std::setprecision(12) << "Param " << i + 1 << ": " <<
                             finalMean[i] << " +- " << finalStde[i] << std::endl;

            std::cout << "Degrees of freedom: " << pointsUsed - dimParam << std::endl;
            std::cout << "Reduced chi-square value obtained is: " <<
                         chiValue <<
                         " +- " <<
                         chiError <<
                         std::endl;
        }
        else{
            for( unsigned i = 0; i < dimParam; i++ ) {
                std::cout << std::setprecision(12) << finalMean[i] << " " <<
                             finalStde[i] << " ";
            }
            std::cout << chiValue << " " << chiError << " " <<
                         pointsUsed - dimParam << std::endl;
        }

    }
    delete gatherFitValue;

    MPI_Finalize();

return 0;
}
