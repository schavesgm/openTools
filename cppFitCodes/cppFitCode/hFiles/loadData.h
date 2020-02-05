//    AUTHOR : SERGIO CHAVES GARCÍA-MASCARAQUE
//    E-MAIL : SERGIOZTESKATE@GMAIL.COM
//    DICIEMBRE DE 2018, SWANSEA, WALES - MADRID, SPAIN


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
       Function to read large data files and avoid memory problems.

        Args:
            const char*:
                Name of the file to load.
            const unsigned:
                Number of rows in the file.
            const unsigned:
                Number of columns in the file.
            const unsigned:
                Control variable to rescale a column by a factor.
            const unsigned:
                Column to be rescaled.
            const double:
                Value to rescale a column with.
            std::string:
                String to skip lines beginning with the provided value.

        Return:
            struct myMat:
                Structure containing the content of the file.
                Value to multiply every value in a column.
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
                        dataStore[rowAux*colSize+i] = resFactor *
                                std::stod( arrayTokens[i] );
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
       Function to read a file automatically. It does not need the size of
       the data file to be read.

        Args:
            const char*:
                Name of the file to load.
            std::string:
                String to skip lines beginning with the provided value.

        Return:
            struct myMat:
                Structure containing the content of the file.
                Value to multiply every value in a column.
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
       Function to slice a column from a myMat structure.

        Args:
            myMat:
                Structure containg the data to be sliced.
            const unsiged:
                Column to be sliced.
            const unisgned:
                Initial row to be sliced.
            unsigned:
                Final positin to be sliced.

        Return:
            struct myMat:
                Structure containing the sliced column from the initial
                myMat structure.
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
