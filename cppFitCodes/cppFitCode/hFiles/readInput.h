//    AUTHOR : SERGIO CHAVES GARCÍA-MASCARAQUE
//    E-MAIL : SERGIOZTESKATE@GMAIL.COM
//    ENERO DE 2019, SWANSEA, WALES - MADRID, SPAIN

#ifndef READINPUT_H
#define READINPUT_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

double* getDouble( const char* fileName, std::string wordToFind, const int numValues ) {
    /*
       Function to get double precision values from a file contained after a
       word defined in wordToFind. We must provide the number of values that we
       look for. The data must be separated by spaces.

        Args:
            const char*:
                Name of the file that we want to explore.
            std::string:
                Name of the keyword that specifies the regular expression that
                we are looking for.
            const in:
                Number of values that will proceed the keyword wanted.

        Return:
            double*:
                Double precision values after keyword. Provided inside a file
                exists a keyword "my_nums 1 2 3", this function will return a
                pointer containing [1.0, 2.0, 3.0].
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

unsigned* getVecUnsigned( const char* fileName, std::string wordToFind,
                          unsigned numValues ) {
    /*
       Function to get unsigned values from a file contained after a
       word defined in wordToFind. We must provide the number of values that we
       look for. The data must be separated by spaces. Note that te values
       returned are going to be signless.

        Args:
            const char*:
                Name of the file that we want to explore.
            std::string:
                Name of the keyword that specifies the regular expression that
                we are looking for.
            const in:
                Number of values that will proceed the keyword wanted.

        Return:
            unsigned*:
                Unsigned values after keyword. Provided inside a file
                exists a keyword "my_nums 1 2 3", this function will return a
                pointer containing [1, 2, 3].
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

unsigned getUnsigned( const char* fileName, std::string wordToFind ) {
    /*
       Function to get unsigned value from a file contained after a
       word defined in wordToFind.  Note that the value returned is going to
       be signless.

        Args:
            const char*:
                Name of the file that we want to explore.
            std::string:
                Name of the keyword that specifies the regular expression that
                we are looking for.

        Return:
            unsigned:
                Unsigned precision values after keyword. Provided inside a file
                exists a keyword "my_nums 1", this function will return an
                unsigned value 1.
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


std::string getString( const char* fileName, std::string wordToFind ) {
    /*
       Function to parse strings from a file contained after a
       word defined in wordToFind.  Note that the value returned is going to
       be signless.

        Args:
            const char*:
                Name of the file that we want to explore.
            std::string:
                Name of the keyword that specifies the regular expression that
                we are looking for.

        Return:
            string:
                String containing the value after the keyword. Provided we
                can find inside the file a expression like "name MyFile", the
                function will return a string with the content "MyFile".
    */

    std::ifstream inFile( fileName );
    std::string line;

    // Extracting words from the file
    while( std::getline( inFile, line ) ) {
        std::string::size_type posWord = line.find(wordToFind);
        if( posWord == 0 ) {
            std::istringstream iss(line);
            std::string auxString;

            // Split the line into words
            while( std::getline( iss, auxString, ' ' ) ) {
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
       Function to parse char pointer from a file contained after a
       word defined in wordToFind.  Note that the value returned is going to
       be signless.

        Args:
            const char*:
                Name of the file that we want to explore.
            std::string:
                Name of the keyword that specifies the regular expression that
                we are looking for.

        Return:
            char*:
                Char pointer containing the value after the keyword. Provided we
                can find inside the file a expression like "name MyFile", the
                function will return a string with the content "MyFile".
    */

    std::ifstream inFile( fileName );
    std::string line;

    // Extracting words from the file
    while( std::getline( inFile, line ) ) {
        std::string::size_type posWord = line.find(wordToFind);
        // Check if the word is the first one
        if( posWord == 0 ) {
            std::istringstream iss(line);
            std::string auxString;
            // Split the line into words
            while( std::getline( iss, auxString, ' ' ) ) {
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
