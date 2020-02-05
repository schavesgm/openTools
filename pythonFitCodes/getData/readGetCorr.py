import numpy as np
import sys
import os

'''
    Script to read raw correlation function data and obtain the mean value and standard deviation
    of the data provided. It acts on the real and imaginary part. Given the file ending in '.uu'
    it saves the results in a file finishing with 'AN.dat'.
'''

if __name__ == '__main__':

    fileName = sys.argv[1]                  ## Name of the file to analyze
    timeExt  = int(sys.argv[2])             ## Number of time points calculated 
    
    dataLoad = np.loadtxt(fileName, skiprows = 1)
    
    lenData = len(dataLoad)                      ## Number of points in the data file
    repData = int(lenData / timeExt)             ## Number of repetitions of the same data

    ## Transform file into columns
    dataReal = np.empty( [timeExt,repData] )
    dataImag = np.empty( [timeExt,repData] )

    for i in range(repData):
        dataReal[:,i] = dataLoad[i*timeExt:(i+1)*timeExt,1]
        dataImag[:,i] = dataLoad[i*timeExt:(i+1)*timeExt,2]
    
    ## Calculate the mean value, and its error
    meanValue, stdeValue = [], []
    for i in range(timeExt):
        meanReal = np.mean( dataReal[i,:] )
        meanImag = np.mean( dataImag[i,:] )
        stdeReal = np.std( dataReal[i,:] ) / np.sqrt(repData)
        stdeImag = np.std( dataImag[i,:] ) / np.sqrt(repData)
        meanValue.append([ meanReal, meanImag ])
        stdeValue.append([ stdeReal, stdeImag ])
    
    fileSave = fileName.replace( '.uu', '_AN.dat' )     ## Name of file to save data in
    subdir = 'resultsData'

    pathHere = os.path.dirname( os.path.relpath(__file__) )
    filePath = os.path.join( pathHere, subdir, fileSave ) 
        
    if not os.path.exists( subdir ):
        os.mkdir( os.path.join( pathHere, subdir ) )

    saveStream = open( filePath, 'w' )

    saveStream.write( '# Time' + '\t' + 'Re(Corr).mean' + '\t' + 'Re(Corr).stde' + '\t' \
                                      + 'Im(Corr).mean' + '\t' + 'Im(Corr).stde' + '\n' )

    for i in range(timeExt):
        saveStream.write( str(i) + '\t' + \
                          str(meanValue[i][0]) + '\t' + str(stdeValue[i][0]) + '\t' + \
                          str(meanValue[i][1]) + '\t' + str(stdeValue[i][1]) + '\n' )


    saveStream.close()
                                    



