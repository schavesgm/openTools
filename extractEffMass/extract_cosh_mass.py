import numpy as np
import scipy.optimize as opt
import re, sys

class CorrAnlz:

    def __init__( self, fileName ):

        # Load the data from the source and get the NT from it
        self.corrData = np.loadtxt( fileName )
        self.NT = int( re.findall( r'(\d+)x32', fileName )[0] )
        self.nConfigs = int( self.corrData.shape[0] / self.NT )

    def bootstrap_corr( self, numBoots ):
        """
            Generate the estimation of the correlation function using bootstrap
        """
        buffBoots = np.empty( self.NT * numBoots )
        for nB in range( numBoots ):
            buffBoots[ self.NT * nB : self.NT * ( nB + 1 ) ] = \
                    self.__res_estimate()

        resCorr = np.empty( [ self.NT, 2 ] )
        avg = self.__avg( buffBoots, self.NT, numBoots )
        std = self.__std( buffBoots, avg, self.NT, numBoots )

        resCorr[:,0] = avg  # Wrap the avg inside the result array
        resCorr[:,1] = std  # Wrap the std inside the result array
        return resCorr

    def bootstrap_effMass( self, numBoots ):
        """
            Generate the estimation of the effective mass using bootstrap
        """
        halfNT = int( self.NT / 2 - 1 )
        buffBoots = np.empty( halfNT * numBoots )
        for nB in range( numBoots ):
            buffBoots[ halfNT * nB  : halfNT * ( nB + 1 ) ] = \
                    self.__cosh_mass()

        resMass = np.empty( [ halfNT, 2 ] )
        avg = self.__avg( buffBoots, halfNT, numBoots )
        std = self.__std( buffBoots, avg, halfNT, numBoots )

        resMass[:,0] = avg  # Wrap the avg inside the result array
        resMass[:,1] = std  # Wrap the std inside the result array
        return resMass

        return self.__avg( buffBoots, halfNT, numBoots )

    def flush_data( self, data, nameFile ):
        """
            Flush a 2-dimensional array into a file with name nameFile
        """
        file = open( nameFile, "w" )
        file.write( '#\t\\tau\tData\tDelta(Data)\n' )
        for i in range( data.shape[0] ):
            file.write( str(i) + '\t' + str(data[i,0]) + '\t' + \
                        str(data[i,1]) + '\n' )

        file.close()

    def __res_estimate( self ):
        """
            Calculate a resampled estimation of the correlation function data
        """

        # Buffer to save the resampled data
        buffRes = np.empty( self.corrData.shape[0] )
        for nC in range( self.nConfigs ):
            # Get a random number from 0 to self.nConfigs - 1
            rC = np.random.randint( self.nConfigs )
            buffRes[ nC * self.NT : self.NT * ( 1 + nC ) ] = \
                self.corrData[ rC * self.NT : self.NT * ( 1 + rC ), 1 ]

        return self.__avg( buffRes, self.NT, self.nConfigs )

    def __cosh_mass( self ):
        """
            Calculate the effective mass for a given resample
        """

        def fCosh( meff, nt, corrData ):    # Function to find the roots to
            coshPart = np.cosh( meff * ( nt - self.NT / 2 ) ) / \
                       np.cosh( meff * ( nt + 1 - self.NT / 2 ) )
            dataPart = corrData[nt] / corrData[nt+1]
            return coshPart - dataPart

        # Array to store the data
        halfNT = int( self.NT / 2 - 1 )
        estMass = np.empty( halfNT )
        # Generate an estimation by resampling the data
        estCorr = self.__res_estimate()
        # Solve the equation for all nT
        for nt in range( halfNT ):
            try:
                logGuess = abs( np.log( estCorr[nt] / estCorr[nt+1] ) )
            except RuntimeWarning:
                logGuess = 0.5
            estMass[nt] = opt.fsolve( fCosh, x0 = logGuess, \
                                      args = ( nt, estCorr ) )

        return estMass


    def __avg( self, data, ntMax, numCopies ):
        """
            Function to calculate the average over an array with numCopies copies
            of ntMax data
        """
        buffAvg = np.empty( ntMax )
        auxAvg = 0
        for nt in range( ntMax ):
            for nC in range( numCopies ):
                auxAvg += data[ nC * ntMax + nt ]
            buffAvg[nt] = auxAvg / numCopies
            auxAvg = 0
        return buffAvg

    def __std( self, data, avgData, ntMax, numCopies ):
        """
            Function to calculate the standard deviation on a array known the average.
            The array contains numCopies copies of ntMax data
        """
        buffStd = np.empty( ntMax )
        auxStd = 0
        for nt in range( ntMax ):
            for nC in range( numCopies ):
                auxStd += ( data[ nC * ntMax + nt ] - avgData[nt] ) ** 2
            buffStd[nt] = np.sqrt( auxStd / ( numCopies - 1 ) )
            auxStd = 0
        return buffStd


if __name__ == '__main__':

    # Get the names of files to analyse and to save the data to
    fileName = sys.argv[1]
    fileSave = sys.argv[2]

    anlzFile = CorrAnlz( fileName )
    massData = anlzFile.bootstrap_effMass( numBoots = 1000 )
    anlzFile.flush_data( massData, fileSave )




