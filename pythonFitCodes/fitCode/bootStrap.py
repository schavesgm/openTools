import fitFunctionNM as fitf            ## Import my fitting functions

import numpy as np
import sys

from mpi4py import MPI                  ## Parallelize calculations using MPI

## Define communicator, size and ranks of the parallelization
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

def bootstrap( arrayBoot, numBootstrap, seed = 123412 ):
    '''
        Bootstrap function to resmple a number of 'numBootstrap' times a given array, it returns
        an array of meanValues calculated using the bootstrap. You can choose a seed in order to 
        randomize in parallel calculations 
    '''

    np.random.seed( seed )              ## Define the seed in each rank

    shape = arrayBoot.shape
    numRows = shape[0]    
    numCols = shape[1]

    bootstrapSamples = np.empty( [numRows, numCols] )
    meanBootstrap = np.empty( [numBootstrap,numCols] )

    for b in range(numBootstrap):
        for i in range( numRows ):
            for j in range( numCols ):
                randChoice = np.random.randint( 0, high = numRows )
                bootstrapSamples[i,j] = arrayBoot[randChoice,j]
        
        for k in range( numCols ):
            meanBootstrap[b,k] = np.mean( bootstrapSamples[:,k] )
    
    return meanBootstrap

if __name__ == '__main__':

    ''' DEFINE HERE YOUR FUNCTION TO FIT ''' 
    ''' ---------------------------------------------------------------------------------------'''
    def f(params, xData ):
        timeExtent = len(xData)
        return params[0] * ( np.exp( - params[1] * xData ) + \
                             np.exp( - params[1] * ( timeExtent - xData ) ) )
    
    ''' ---------------------------------------------------------------------------------------'''

    if rank == 0:        ## Master rank
        
        ''' DEFINE YOUR INITIAL PARAMETERS HERE: ----------------------------------------------'''
        initGuess = np.array([2e-16,0.01])      ## Initial guess of the parameter space
        nReps = 4                               ## Each rank will calculate nReps equal fits
        numBootstrap = 10000                    ## Number of bootstrap calculation in total


        ''' -----------------------------------------------------------------------------------'''
        
        ''' LOAD THE DATA FROM GIVEN FILE: ----------------------------------------------------'''
        fileName = sys.argv[1]                                          ## Data file to fit
        beginFit, endFit = int(sys.argv[2]), int(sys.argv[3])           ## Crop data if wanted

        dataLoad = np.loadtxt( fileName, skiprows = 1 )
        timeExtent  = len(dataLoad)
        if endFit <= beginFit:                                          ## Consistent limits
            endFit = timeExtent

        xData = dataLoad[beginFit:endFit,0]       ## Points in the x-axis
        yData = dataLoad[beginFit:endFit,1]       ## Points in the y-axis
        yError = dataLoad[beginFit:endFit,2]      ## Errors associated with y-axis
    
        paramSize = len(initGuess)                ## Dimension of parameter space


        ''' -----------------------------------------------------------------------------------'''
        
        ''' MPI OBJECTS NEEDED: ---------------------------------------------------------------'''
        
        initGuessMPI = []
        xDataMPI, yDataMPI, yErrorMPI = [], [], []
        seedMPI, nRepsMPI = [], []
        
        getFitParamsMPI = []

        for i in range( size ):
            initGuessMPI.append( initGuess * np.random.rand(paramSize) ) 
            xDataMPI.append( xData )
            yDataMPI.append( yData )
            yErrorMPI.append( yError )
            seedMPI.append( np.random.randint( 0, 10000 ) * i )
            nRepsMPI.append( nReps )

            getFitParamsMPI.append( 0 ) 

        ''' -----------------------------------------------------------------------------------'''
    else:

        ''' INITIALIZE SAME ARRAYS AT EACH RANK TO AVOID PROBLEMS: ----------------------------'''
        initGuessMPI = None
        xDataMPI, yDataMPI, yErrorMPI = None, None, None
        seedMPI, nRepsMPI = None, None
        
        getFitParamsMPI = None 

        ''' -----------------------------------------------------------------------------------'''

        ''' NEEDED FOR THE BOOTSTRAP CALCULATION AT THE END: ----------------------------------'''
        arrayBootMPI, numBootstrapMPI = None, None

        getBootMPI = None

        ''' -----------------------------------------------------------------------------------'''

    ''' SCATTER ALL THE DATA TO EACH RANK FROM MASTER RANK: -----------------------------------'''

    initGuessSCAT = comm.scatter( initGuessMPI, root = 0 )
    xDataSCAT = comm.scatter( xDataMPI, root = 0 )
    yDataSCAT = comm.scatter( yDataMPI, root = 0 )
    yErrorSCAT = comm.scatter( yErrorMPI, root = 0 )
    seedSCAT = comm.scatter( seedMPI, root = 0 )
    nRepsSCAT = comm.scatter( nRepsMPI, root = 0 )

    getFitParamsSCAT = comm.scatter( getFitParamsMPI, root = 0 )

    getFitParamsSCAT = fitf.fitNtimes( f, nRepsSCAT, initGuessSCAT, \
                                       xDataSCAT, yDataSCAT, yErrorSCAT, seedSCAT )
    
    resultFit = comm.gather( getFitParamsSCAT, root = 0 )
    
    ''' ---------------------------------------------------------------------------------------'''

    ''' COLLECT THE DATA INTO THE MASTER RANK TO PROCEED WITH BOOTSTRAP -----------------------'''
    if rank == 0:
        
        ## Collect the data into a matrix
        dataFit = []
        for i in range( size ):
            for j in range( nReps ):
                dataFit.append( resultFit[i][j] )
        
        dataFit = np.array( dataFit )           ## Convert the data into array for convenience

        
        ''' MPI OBJECTS NEEDED: ---------------------------------------------------------------'''
        arrayBootMPI, numBootstrapMPI = [], []                      ## We will use last seed array

        getBootMPI = []

        for i in range( size ):
            arrayBootMPI.append( dataFit )
            numBootstrapMPI.append( int( numBootstrap / size ) )    ## Sharing is caring

            getBootMPI.append( 0 )

    
    ''' ---------------------------------------------------------------------------------------'''

    ''' SCATTER ALL THE DATA TO EACH RANK FROM MASTER RANK: -----------------------------------'''

    arrayBootSCAT = comm.scatter( arrayBootMPI, root = 0 )
    numBootstrapSCAT = comm.scatter( numBootstrapMPI, root = 0 )

    getBootSCAT = comm.scatter( getBootMPI, root = 0 )

    getBootSCAT = bootstrap( arrayBootSCAT, numBootstrapSCAT, seedSCAT )

    getBootstrapped = comm.gather( getBootSCAT, root = 0 )

    ''' ---------------------------------------------------------------------------------------'''

    ''' COLLECT THE DATA INTO THE MASTER RANK TO FINALIZE THE CALCULATION: --------------------'''
    if rank == 0:
        
        # Collect the data into a matrix
        dataBoots = []
        for i in range( size ):
            for j in range( int(numBootstrap / size ) ) :
                dataBoots.append( getBootstrapped[i][j] )

        dataBoots = np.array( dataBoots ) 
        

        meanBootValue, stdeBootValue = [], []
        for i in range( paramSize ):
            meanBootValue.append( np.mean( dataBoots[:,i] ) )
            stdeBootValue.append( np.std( dataBoots[:,i] ) ) 
        
        
        for i in range( paramSize ):
            print( 'Value obtained for ', i+1, ' parameter is: ',\
                    meanBootValue[i], ' +- ', stdeBootValue[i] )
        
        ## PLOT THINGS TO CHECK THE RESULT
        import matplotlib.pyplot as plt

        plt.errorbar( xData, yData, yerr = yError )
        plt.plot( xData, f(meanBootValue, xData) )
        plt.show()

    ''' ---------------------------------------------------------------------------------------'''

