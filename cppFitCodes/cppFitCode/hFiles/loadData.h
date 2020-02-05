//------------------------------------------------------------------------------------------------//
//    AUTHOR : SERGIO CHAVES GARCÍA-MASCARAQUE
//    E-MAIL : SERGIOZTESKATE@GMAIL.COM
//    DICIEMBRE DE 2018, SWANSEA, WALES - MADRID, SPAIN
//------------------------------------------------------------------------------------------------//


#ifndef LOAD_H
#define LOAD_H

#include <fstream>          // Used to load data from file
#include <string>           // Used to load data from file
#include <vector>           // Used to load data from file
#include <sstream>          // Used to load data from file
#include <iterator>         // Used to load large data from file
#include <stdexcept>        // Throw exceptions

#define NOTUSE 0

struct myMat {
    double* m_Matrix;
    unsigned rowSize;
    unsigned colSize;
};

struct myMat getLargeData( const char* fileName,
                           const unsigned rowSize,
                           const unsigned colSize,
                           const unsigned RESCALE = NOTUSE,
                           const unsigned colRescale = NOTUSE,
                           const double resFactor = NOTUSE,
                           std::string commChar = "#" ) {
    /*
     * Function to read large data files. It needs the number of rows and column
     * included in the file to work. It reads the file as a matrix with dimension
     * rowSize * colSize
     *
     * Args:
     * const char*      : Name of the file to extract
     * const int        : Number of rows contained in the file
     * const int        : Number of columns contained in the file
     * std::string      : Comment character. It must be at the beginning of the line
     *                    and just one character in each comment. Therefore '#' works,
     *                    '##' will not work.
     *
     * Return:
     * struct myMat     : Structure containing the data extracted from the file and its
     *                    properties.
     */

    if( RESCALE == NOTUSE ) {
        std::ifstream fileStream(fileName);
        std::string dummyLine;
        double* dataStore = new double[rowSize*colSize];
        double xTemp;

        int rowAux = 0;
        while( std::getline( fileStream, dummyLine ) ) {
            std::istringstream ss(dummyLine);
            std::istream_iterator<std::string> begin(ss), end;

            // Putting all tokens inside a vector
            std::vector<std::string> arrayTokens( begin, end );

            if( arrayTokens[0] != commChar ) {
                for( int i = 0; i < arrayTokens.size(); i++ ) {
                    dataStore[rowAux*colSize+i] = std::stod( arrayTokens[i] );
                }
            rowAux += 1;
            }
        }
        fileStream.close();
        struct myMat Return = { dataStore, rowSize, colSize };
        return Return;
    }
    else {
        std::ifstream fileStream(fileName);
        std::string dummyLine;
        double* dataStore = new double[rowSize*colSize];
        double xTemp;

        int rowAux = 0;
        while( std::getline( fileStream, dummyLine ) ) {
            std::istringstream ss(dummyLine);
            std::istream_iterator<std::string> begin(ss), end;

            // Putting all tokens inside a vector
            std::vector<std::string> arrayTokens( begin, end );

            if( arrayTokens[0] != commChar ) {
                for( int i = 0; i < arrayTokens.size(); i++ ) {
                    if( i == colRescale )
                        dataStore[rowAux*colSize+i] = resFactor * std::stod( arrayTokens[i] );
                    else
                        dataStore[rowAux*colSize+i] = std::stod( arrayTokens[i] );
                }
            rowAux += 1;
            }
        }
        fileStream.close();
        struct myMat Return = { dataStore, rowSize, colSize };
        return Return;
    }

}

struct myMat loadData( const char* fileName, const char* compChar = "#" ) {
    /*
     * Function to loadData automatically, it generates an structure containing the
     * data and the number of rows and columns. The function does not work for large
     * data files as it copies a vector iteratively.
     *
     * Args:
     * const char*      : Name of the file that we want to load
     * const char*      : Comment character in the file to be skipped. Note that the
     *                    character must be located at the beginning and just one char
     *
     * Return:
     * struct myMat     : Structure containing the data and its properties.
     */

    std::ifstream fileStream(fileName);

    // KEEPS TRACK OF THE COLUMN AND ROW SIZES
    unsigned colSize = 0;
    unsigned rowSize = 0;
    unsigned numComm = 0;

    // READ DATA USING STRINGS
    std::string dummyLine;
    unsigned numPoints = 0;
    double xTemp;
    double* vecStrorage = nullptr;

    if( fileStream.is_open() && fileStream.good() ) {
        while( std::getline(fileStream, dummyLine) ) {

            rowSize += 1;
            std::stringstream strStream(dummyLine);
            char* firstChar = &dummyLine[0];

            if( firstChar[0] != compChar[0] ) {

                colSize = 0;
                while( 1 ) {
                    strStream >> xTemp;
                    if( !strStream )
                        break;

                    double* tempArr = new double[numPoints+1];
                    colSize += 1;

                    // vecStorage grows from file
                    // Copy from vecStorage to vecStorage + numPoints the array tempArray
                    std::copy( vecStrorage, vecStrorage + numPoints, tempArr );
                    tempArr[numPoints] = xTemp;
                    vecStrorage = tempArr;
                    numPoints += 1;
                }
            }
            else
                numComm += 1;
        }
    }
    else
        std::invalid_argument("File name given is not correct. \n");

    unsigned realRow = rowSize - numComm;
    double* getData = new double[realRow*colSize];
    for( unsigned i = 0; i < realRow * colSize; i++ )
        getData[i] = vecStrorage[i];

    // Return values using Return structure
    myMat r = { getData, realRow, colSize };

    delete [] vecStrorage;
    fileStream.close();
    return r;
}

struct myMat sliceCol( myMat inMatrix,
                      const unsigned getCol,
                      const unsigned iniPos = 0,
                      unsigned finPos = 0 ) {
    /*
     * Function to slice a data matrix at one row. It needs data contained inside
     * a structure myMat.
     *
     * Args:
     * struct myMat     : Structure containing the data matrix to be sliced
     * const unsigned   : Column in the data to be extracted
     * const unsigned   : Initial row to be extracted
     * unsigned         : Final row to be extracted
     *
     * Return:
     * struct myMat     : Stucture containing a vector with the sliced data
     */

    unsigned rowSize = inMatrix.rowSize;
    unsigned colSize = inMatrix.colSize;

    if( finPos == 0 )
        finPos = rowSize;

    double* retCol = new double[finPos-iniPos];
    int auxRow = iniPos;
    for( unsigned i = 0; i < finPos - iniPos; i++ ) {
        retCol[i] = inMatrix.m_Matrix[auxRow*colSize+getCol];
        auxRow += 1;
    }

    myMat retMat = {retCol, finPos - iniPos, 1 };
    return retMat;

}

#endif
