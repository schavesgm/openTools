//------------------------------------------------------------------------------------------------//
//    AUTHOR : SERGIO CHAVES GARCÍA-MASCARAQUE
//    E-MAIL : SERGIOZTESKATE@GMAIL.COM
//    ENERO DE 2019, SWANSEA, WALES - MADRID, SPAIN
//------------------------------------------------------------------------------------------------//

#ifndef READINPUT_H
#define READINPUT_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

double* getDouble( const char* fileName, std::string wordToFind, const int numValues ) {
    /*
     * Function to get double precision values from data file contained after some word
     * defined in wordToFind. We must provide the number of values that we look for. The
     * data must be separated by spaces.
     *
     * Args:
     * const char*      : Name of the file that we want to explore
     * std::string      : String containing the name of the word that will characterize the
     *                    regular expression that we are looking for
     * const int        : Number of values that will be located after the regular expression
     *
     * Return:
     * double*          : Pointer of array containing the data
     */

    std::ifstream inFile( fileName );
    std::string line;
    double* retValues = new double[numValues];

    // Extracting words from the file
    while( std::getline( inFile, line ) ) {
        std::string::size_type posWord = line.find(wordToFind);
        if( posWord == 0 ) {                                // Check if the word is the first one
            std::istringstream iss(line);
            std::string auxString;
            int auxValue = 0;
            while( std::getline( iss, auxString, ' ' ) ) {  // Split the line into words
                if( auxString == wordToFind )
                    continue;
                retValues[auxValue] = stod( auxString );
                auxValue += 1;
            }
        }
    }

    inFile.close();
    return retValues;
}

unsigned getUnsigned( const char* fileName, std::string wordToFind ) {
    /*
     * Function to get unsigned integer value from data file contained after some word
     * defined in wordToFind. The data must be separated by spaces.
     *
     * Args:
     * const char*      : Name of the file that we want to explore
     * std::string      : String containing the name of the word that will characterize the
     *                    regular expression that we are looking for
     *
     * Return:
     * unsigned         : Value containing the data wanted
     */

    std::ifstream inFile( fileName );
    std::string line;

    // Extracting words from the file
    while( std::getline( inFile, line ) ) {
        std::string::size_type posWord = line.find(wordToFind);
        if( posWord == 0 ) {
            std::istringstream iss(line);
            std::string auxString;
            while( std::getline( iss, auxString, ' ' ) ) {
                if( auxString == wordToFind )
                    continue;
                return stoul( auxString );
            }
        }
    }
    inFile.close();
    return 0;
}

unsigned* getVecUnsigned( const char* fileName, std::string wordToFind, unsigned numValues ) {
    /*
     * Function to get unsigned integer values from data file contained after some word
     * defined in wordToFind. We must provide the number of values that we look for. The
     * data must be separated by spaces.
     *
     * Args:
     * const char*      : Name of the file that we want to explore
     * std::string      : String containing the name of the word that will characterize the
     *                    regular expression that we are looking for
     * const int        : Number of values that will be located after the regular expression
     *
     * Return:
     * unsigned*          : Pointer of array containing the data
     */


    std::ifstream inFile( fileName );
    std::string line;
    unsigned* retValues = new unsigned[numValues];

    // Extracting words from the file
    while( std::getline( inFile, line ) ) {
        std::string::size_type posWord = line.find(wordToFind);
        if( posWord == 0 ) {                                // Check if the word is the first one
            std::istringstream iss(line);
            std::string auxString;
            int auxValue = 0;
            while( std::getline( iss, auxString, ' ' ) ) {  // Split the line into words
                if( auxString == wordToFind )
                    continue;
                retValues[auxValue] = stoul( auxString );
                auxValue += 1;
            }
        }
    }

    inFile.close();
    return retValues;
}

std::string getString( const char* fileName, std::string wordToFind ) {
    /*
     * Function to get unsigned string value from data file contained after some word
     * defined in wordToFind. The data must be separated by spaces.
     *
     * Args:
     * const char*      : Name of the file that we want to explore
     * std::string      : String containing the name of the word that will characterize the
     *                    regular expression that we are looking for
     *
     * Return:
     * std::string      : String contained after the regular expression
     */

    std::ifstream inFile( fileName );
    std::string line;

    // Extracting words from the file
    while( std::getline( inFile, line ) ) {
        std::string::size_type posWord = line.find(wordToFind);
        if( posWord == 0 ) {                                // Check if the word is the first one
            std::istringstream iss(line);
            std::string auxString;
            while( std::getline( iss, auxString, ' ' ) ) {  // Split the line into words
                if( auxString == wordToFind )
                    continue;
                return auxString;
            }
        }
    }
    inFile.close();
    return 0;
}

const char* getConstChar( const char* fileName, std::string wordToFind ) {
    /*
     * Function to get const char pointer from data file contained after some word
     * defined in wordToFind. The data must be separated by spaces.
     *
     * Args:
     * const char*      : Name of the file that we want to explore
     * std::string      : String containing the name of the word that will characterize the
     *                    regular expression that we are looking for
     *
     * Return:
     * const char*      : Const char pointer containing the data
     */

    std::ifstream inFile( fileName );
    std::string line;

    // Extracting words from the file
    while( std::getline( inFile, line ) ) {
        std::string::size_type posWord = line.find(wordToFind);
        if( posWord == 0 ) {                                // Check if the word is the first one
            std::istringstream iss(line);
            std::string auxString;
            while( std::getline( iss, auxString, ' ' ) ) {  // Split the line into words
                if( auxString == wordToFind )
                    continue;
                return auxString.c_str();
            }
        }
    }
    inFile.close();
    return 0;
}

#endif
