import numpy as np
import copy as cp
import sys

__all__ = ['fitOnceNM', 'fitNtimes']

def chiSquare(f, args, xData, yData, yError = 0 ):
    '''
        @arg f (function) : function that we would like to fit with the data.
                            Must return a double and be able to accept numpy
                            arrays as input
        @arg args (array) : values that we would like to fit
        @arg xData(array) : xData extracted from the experiment that we would
                            like to fit
        @arg yData(array) : yData extracted from the experiment that we would
                            like to fit
        @arg yError(array): error in the measured function used to calculate the
                            fit function
        return            : double with the chiSquare value at the current
                            status of the coeffs
    '''
    if yError.all() == 0:
        return np.sum( ( f(args, xData) - yData ) ** 2 ) / len(xData)
    else:
        return np.sum( ( (f(args, xData) - yData) / yError ) ** 2 ) / len(xData)


def fitOnceNM( f, initGuess, xData, yData, yError = 0,
               seed = 1231523, maxIter = 1e6, getChiSquare = 0 ):
    '''
        @arg f (function) : function that we would like to optimize. This
                            function must return a double and be able to get
                            numpy arrays as input

        @initGuess (array): array with the initial guess of the points
                            that minimize the function
        @maxIter (int)    : integer with the maximum number of iterations
                            that we accept
        return: tuple     : tuple( Points that minimize the function,
                            minimum value achieved)
    '''

    np.random.seed( seed )

    ## DEFINE THE CONSTANTS OF THE ALGORITHM
    alpha = 1.0
    gamma = 2.0
    rho = 0.5
    sigma = 0.5

    noImpThreshold = 1E-8
    noImprovementBreak = 1000

    # INIT THE ALGORITHM
    N = initGuess.shape[0];
    evalFunc = chiSquare(f,initGuess,xData,yData,yError)
    noImprov = 0
    imValues = [[initGuess, evalFunc]]

    for i in range( N ):
        newPoint = cp.copy(initGuess) + initGuess * np.random.random(N)
        imNewPoint= chiSquare(f,newPoint,xData,yData,yError)
        imValues.append([newPoint, imNewPoint])

    ### SIMPLEX ITERATION
    iterStep = 0
    while 1:

        # ORDER THE EVALUATED FUNCTION
        imValues.sort( key = lambda x: x[1] )   ## Sort using the evaluated function
        lowestValue = imValues[0][1]            ## Best value == lowest value

        # BREAK IF MAX ITERATION VALUE IS ACHIEVED
        if maxIter and iterStep >= maxIter:
            bestData = imValues[0][0]           ## Best minimum obtained
            evalChiSQ = imValues[0][1]          ## Value of the chiSquare achieved

            if getChiSquare == 0:                ## Control variable to print out
                                                 ## the chiSquare
                return bestData
            else:
                return bestData, evalChiSQ


        iterStep += 1

        # BREAK IN CASE OF NON-IMPORTANT IMPROVEMENT
        ## print( 'lowest value in iter ', iterStep, 'is ', lowestValue )

        if lowestValue <  evalFunc - noImpThreshold:
            noImprov = 0
            evalFunc = lowestValue
        else:
            noImprov += 1

        if noImprov >= noImprovementBreak:      ## Number of points without improvement
            bestData = imValues[0][0]           ## Best minimum obtained
            evalChiSQ = imValues[0][1]          ## Value of the chiSquare achieved

            if getChiSquare == 0:                ## Control variable to print
                                                 ## out the chiSquare
                return bestData
            else:
                return bestData, evalChiSQ

        # CALCULATE THE CENTROID OF THE POINTS
        centroidPoint = [0.] * N

        for tup in imValues[:-1]:               ## Get all the tuples in the 'matrix'
            for i, c in enumerate(tup[0]):      ## Disect both tuples, we just need
                                                ## the first one
                centroidPoint[i] += c / (len(imValues)-1)       ## Sum the centroid

        # REFLECTION PART OF THE ALGORITHM
        xRef = centroidPoint + alpha * ( centroidPoint - imValues[-1][0])
        refEval = chiSquare(f,xRef,xData,yData,yError)

        if imValues[0][1] <= refEval < imValues[-2][1]:
            del imValues[-1]                               ## Remove the last point
            imValues.append([xRef, refEval])               ## Append the reflection
                                                           ## as highest point
            continue

        # EXPANSION PART OF THE ALGORITHM
        if refEval < imValues[0][1]:

            xExp = centroidPoint + gamma * ( xRef - imValues[-1][0] )
            expValue = chiSquare(f,xExp,xData,yData,yError)
            if expValue < refEval:
                del imValues[-1]
                imValues.append([xExp, expValue])

                continue
            else:
                del imValues[-1]
                imValues.append([xRef, refEval])
                continue

        # CONTRACTION PART OF THE ALGORITHM
        xCont = centroidPoint + rho * ( centroidPoint - imValues[-1][0] )
        contEval = chiSquare(f,xCont,xData,yData,yError)

        if contEval < imValues[-1][1]:

            del imValues[-1]
            imValues.append([xCont, contEval])
            continue

        # REDUCTION PART
        nimValues = []

        for tup in imValues:
            xRed =  imValues[0][0] + sigma * ( tup[0] - imValues[0][0] )
            redEval = chiSquare(f,xRed,xData,yData,yError)
            nimValues.append( [xRed, redEval] )
        imValues = nimValues

def fitNtimes( f, nReps, initGuess, xData, yData, yError = 0, seed = 1231241 ):
    '''
        Routine that calls "fitOnceNM" nReps times in order to get some
        statistics. The value of initGuess is multiplied by some random value
        in order to randomize the nCalculations and explore the space better.
    '''

    getFit = []
    for i in range( nReps ):
        stepSpace = np.sqrt( np.dot( initGuess, initGuess ) )
        initGuess += stepSpace * np.random.random(len(initGuess))
        getFit.append( fitOnceNM( f, initGuess, xData, yData, yError, seed ) )

    return getFit

if __name__ == "__main__":
    pass
