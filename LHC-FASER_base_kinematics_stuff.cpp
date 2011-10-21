/*
 * LHC-FASER_base_kinematics_stuff.cpp
 *
 *  Created on: 02 Mar 2011
 *      Authors: Ben O'Leary (benjamin.oleary@gmail.com)
 *               Jonas Lindert (jonas.lindert@googlemail.com)
 *               Carsten Robens (carsten.robens@gmx.de)
 *      Copyright 2010 Ben O'Leary, Jonas Lindert, Carsten Robens
 *
 *      This file is part of LHC-FASER.
 *
 *      LHC-FASER is free software: you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation, either version 3 of the License, or
 *      (at your option) any later version.
 *
 *      LHC-FASER is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *
 *      You should have received a copy of the GNU General Public License
 *      along with LHC-FASER.  If not, see <http://www.gnu.org/licenses/>.
 *
 *      The GNU General Public License should be in GNU_public_license.txt
 *      the files of LHC-FASER are:
 *      LHC-FASER.hpp
 *      LHC-FASER.cpp
 *      LHC-FASER_electroweak_cascade_stuff.hpp
 *      LHC-FASER_electroweak_cascade_stuff.cpp
 *      LHC-FASER_full_cascade_stuff.hpp
 *      LHC-FASER_full_cascade_stuff.cpp
 *      LHC-FASER_global_stuff.hpp
 *      LHC-FASER_global_stuff.cpp
 *      LHC-FASER_input_handling_stuff.hpp
 *      LHC-FASER_input_handling_stuff.cpp
 *      LHC-FASER_kinematics_stuff.hpp
 *      LHC-FASER_kinematics_stuff.cpp
 *      LHC-FASER_lepton_distributions.hpp
 *      LHC-FASER_lepton_distributions.cpp
 *      LHC-FASER_signal_calculator_stuff.hpp
 *      LHC-FASER_signal_calculator_stuff.cpp
 *      LHC-FASER_signal_data_collection_stuff.hpp
 *      LHC-FASER_signal_data_collection_stuff.cpp
 *      LHC-FASER_sparticle_decay_stuff.hpp
 *      LHC-FASER_sparticle_decay_stuff.cpp
 *      and README.LHC-FASER.txt which describes the package.
 *
 *      LHC-FASER also requires CppSLHA. It should be found in a subdirectory
 *      included with this package.
 *
 *      LHC-FASER also requires grids of lookup acceptanceValues. These should also be
 *      found in a subdirectory included with this package.
 */

#include "LHC-FASER_base_kinematics_stuff.hpp"

namespace LHC_FASER
{
  acceptanceCutSet::acceptanceCutSet() :
    beamEnergy( -1 ),
    primaryLeptonCut( -1.0 ),
    secondaryLeptonCut( -1.0 ),
    jetCut( -1.0 ),
    excludedStandardModelProducts( NULL )
  // this constructor sets the acceptanceValues to -1.0, which is the default "unset"
  // value.
  {
    // just an initialization list.
  }

  acceptanceCutSet::acceptanceCutSet(
                                  acceptanceCutSet const* const copyPointer ) :
    beamEnergy( copyPointer->beamEnergy ),
    primaryLeptonCut( copyPointer->primaryLeptonCut ),
    secondaryLeptonCut( copyPointer->secondaryLeptonCut ),
    jetCut( copyPointer->jetCut ),
    excludedStandardModelProducts( copyPointer->excludedStandardModelProducts )
  // this constructor copies the acceptanceValues from a given acceptanceCutSet.
  {
    // just an initialization list.
  }

  acceptanceCutSet::~acceptanceCutSet()
  {
    // does nothing.
  }



  /* these are for which elements of the acceptanceValues vector of vectors
   * of vectors to use.
   * N.B.! this bit is HIGHLY format-dependent, since the acceptanceValues are
   * read in in the order:
   * 1: light - heavy
   * 2: light - medium
   * 3: light - light
   * 4: medium - heavy
   * 5: medium - medium
   * 6: heavy - heavy
   * though these are stored in the vector as elements
   * 0, 1, 2, 3, 4, & 5, which is why the comments below seem to not
   * match the code.
   */
  int const acceptanceGrid::indexForLightLightNeutralinoPair( 2 );
  int const acceptanceGrid::indexForLightMediumNeutralinoPair( 1 );
  int const acceptanceGrid::indexForLightHeavyNeutralinoPair( 0 );
  int const acceptanceGrid::indexForMediumMediumNeutralinoPair( 4 );
  int const acceptanceGrid::indexForMediumHeavyNeutralinoPair( 3 );
  int const acceptanceGrid::indexForHeavyHeavyNeutralinoPair( 5 );


  acceptanceGrid::acceptanceGrid( std::string const* const gridFileLocation,
                                  inputHandler const* const inputShortcut ) :
      inputShortcut( inputShortcut )
  {
    // debugging:
    /**std::cout
    << std::endl
    << "debugging: acceptanceGrid::acceptanceGrid( " << *gridFileLocation
    << " ) called";
    std::cout << std::endl;**/

    // read in the grid:
    if( NULL != gridFileLocation )
    {
      // open the file:
      int const maxLengthOfGridFile( 100000 );
      CppSLHA::hash_commented_file_reader inputFileReader( gridFileLocation,
                                                           maxLengthOfGridFile,
                                                           false );
      /* a hundred thousand lines should be long enough for any grid file. if
       * readIn() is reading in more lines than that, it's probably because
       * of some file corruption.
       * also, if it's not enough, change it here.
       */
      std::string inputLineAsString;
      // prepare a string to parse the line.
      std::stringstream inputLineAsStream;
      // prepare a stringstream to parse the line.
      // the masses from the last line read & the currently-read line:
      double lastSquarkMass( 0.0 );
      double lastGluinoMass( 0.0 );
      double currentSquarkMass( 0.0 );
      double currentGluinoMass( 0.0 );
      std::vector< double > neutralinoMasses;
      lowestSquarkMass = 0.0;
      lowestGluinoMass = 0.0;
      double currentValue;
      // the value currently being read in.
      bool gridSizeStillUnknown( true );
      // we start out not knowing the step size for the grid (but we do
      // assume that it is a square grid).
      std::vector< double >* currentAcceptanceVector( NULL );
      // this holds all the acceptance acceptanceValues for a grid point.
      std::vector< std::vector< double >* >* currentNeutralinoVector( NULL );
      // this holds all the acceptanceValues for a set of points with the same
      // squark & gluino masses, but varying neutralino masses.
      std::vector< std::vector< std::vector< double >* >* >*
      currentGluinoVector( NULL );
      // this holds all the acceptanceValues for a set of points with the same
      // squark masses, but varying gluino masses.

      while( inputFileReader.read_line( &inputLineAsString ) )
        // the evaluation of the conditional reads in the next line.
      {
        // read in the line for interpretting:
        inputLineAsStream.clear();
        inputLineAsStream.str( inputLineAsString );
        // update what the masses for the last line were:
        lastSquarkMass = currentSquarkMass;
        lastGluinoMass = currentGluinoMass;
        inputLineAsStream >> currentSquarkMass
                          >> currentGluinoMass;
        inputLineAsStream >> currentValue;
        if( gridSizeStillUnknown )
        {
          neutralinoMasses.push_back( currentValue );
        }
        inputLineAsStream >> currentValue;
        if( gridSizeStillUnknown )
        {
          neutralinoMasses.push_back( currentValue );
        }
        currentAcceptanceVector = new std::vector< double >;
        while( !(inputLineAsStream.eof()) )
        {
          inputLineAsStream >> currentValue;
          currentAcceptanceVector->push_back( currentValue );
        }
        if( currentSquarkMass > lastSquarkMass )
          // if the squark mass has changed, we need a new vector for the
          // acceptanceValues for varying gluino masses.
        {
          // 1st record the last set of points with the same squark mass:
          if( NULL != currentGluinoVector )
            // if we have a vector to add...
          {
            currentGluinoVector->push_back( currentNeutralinoVector );
            acceptanceValues.push_back( currentGluinoVector );

            // debugging:
            /**std::cout
            << std::endl
            << "pushed back a neutralino vector of size "
            << currentNeutralinoVector->size()
            << " and a gluino vector of size "
            << currentGluinoVector->size();**/
          }
          // prepare a new vector for the new squark mass:
          currentGluinoVector
          = new std::vector< std::vector< std::vector< double >* >* >;
          // prepare a new vector for the new gluino mass:
          currentNeutralinoVector = new std::vector< std::vector< double >* >;
        }
        else if( currentGluinoMass > lastGluinoMass )
          /* if the gluino mass has changed (but the squark mass has not
           * changed, we need a new vector for the acceptanceValues for varying
           * neutralino masses.
           */
        {
          // 1st record the last set of points with the same gluino mass:
          currentGluinoVector->push_back( currentNeutralinoVector );

          // debugging:
          /**std::cout
          << std::endl
          << "pushed back a neutralino vector of size "
          << currentNeutralinoVector->size();**/

          // prepare a new vector for the new gluino mass:
          currentNeutralinoVector = new std::vector< std::vector< double >* >;
        }

        /*if( !gridSizeStillUnknown
            &&
            ( acceptanceColumns > currentAcceptanceVector->size() ) )
        {
          acceptanceColumns = currentAcceptanceVector->size();
          std::cout
          << std::endl
          << "LHC-FASER warning! read in an acceptance grid line that"
          << " is shorter than the previous one, so the grid is"
          << " \"ragged\"; hence *all* the lines for this grid will be"
          << " truncated to this smaller length.";
          std::cout << std::endl;
        }*/
        currentNeutralinoVector->push_back( currentAcceptanceVector );

        if( gridSizeStillUnknown
            &&
            ( 0.0 < lastSquarkMass )
            &&
            ( 0.0 < lastGluinoMass )
            &&
            ( ( currentSquarkMass > lastSquarkMass )
              ||
              ( currentGluinoMass > lastGluinoMass ) ) )
          // if we still need to work out the grid step size & we have read
          // in enough to do so...
        {
          //acceptanceColumns = currentAcceptanceVector->size();
          // first we take note of the 1st point, which should have the
          // lowest co-ordinate acceptanceValues:
          lowestSquarkMass = lastSquarkMass;
          lowestGluinoMass = lastGluinoMass;
          // then we work out the grid step size:
          if( currentGluinoMass > lastGluinoMass )
            // if it was the gluino mass that incremented...
          {
            scoloredMassStepSize = ( currentGluinoMass - lastGluinoMass );
          }
          else
            // otherwise it was the squark mass that incremented...
          {
            scoloredMassStepSize = ( currentSquarkMass - lastSquarkMass );
          }
          // debugging:
          /**std::cout
          << std::endl
          << "scoloredMassStepSize = " << scoloredMassStepSize;**/

          // now we need to work out the neutralino mass ratios to the
          // lighter colored sparticle:
          double lowestNeutralinoMass( lowestSquarkMass );
          double middleNeutralinoMass( lowestSquarkMass );
          double highestNeutralinoMass( lowestSquarkMass );
          for( std::vector< double >::const_iterator
               massIterator( neutralinoMasses.begin() );
               neutralinoMasses.end() > massIterator;
               ++massIterator )
          {
            if( *massIterator < lowestNeutralinoMass )
            {
              highestNeutralinoMass = middleNeutralinoMass;
              middleNeutralinoMass = lowestNeutralinoMass;
              lowestNeutralinoMass = *massIterator;
            }
            if( ( *massIterator < middleNeutralinoMass )
                &&
                ( *massIterator > lowestNeutralinoMass ) )
            {
              highestNeutralinoMass = middleNeutralinoMass;
              middleNeutralinoMass = *massIterator;
            }
            if( ( *massIterator < highestNeutralinoMass )
                &&
                ( *massIterator > middleNeutralinoMass ) )
            {
              highestNeutralinoMass = *massIterator;
            }
            double lighterScoloredMass( lastSquarkMass );
            if( lastGluinoMass < lastSquarkMass )
            {
              lighterScoloredMass = lastGluinoMass;
            }
            lowNeutralinoMassRatio = ( lowestNeutralinoMass
                                       / lighterScoloredMass );
            middleNeutralinoMassRatio = ( middleNeutralinoMass
                                          / lighterScoloredMass );
            highNeutralinoMassRatio = ( highestNeutralinoMass
                                        / lighterScoloredMass );
          }
          gridSizeStillUnknown = false;
        }
      }  // end of while loop going over the lines of the file.
      // we still have to add the last vectors:
      //acceptanceColumns = currentAcceptanceVector->size();
      currentGluinoVector->push_back( currentNeutralinoVector );
      acceptanceValues.push_back( currentGluinoVector );

      // debugging:
      /**std::cout
      << std::endl
      << "pushed back last currentNeutralinoVector of size "
      << currentNeutralinoVector->size()
      << " & last currentGluinoVector of size "
      << currentGluinoVector->size();
      std::cout << std::endl;**/

      // the last masses should be the maximal acceptanceValues:
      highestSquarkMass = currentSquarkMass;
      highestGluinoMass = currentGluinoMass;
    }
    else
      // otherwise, a NULL pointer was given instead of a pointer to the
      // grid file name string.
    {
      std::cout
      << std::endl
      << "LHC-FASER::error! asked to open a file with a NULL pointer!";
      std::cout << std::endl;
    }
    // debugging:
    /**std::cout
    << std::endl
    << "debugging: acceptanceGrid::readIn() finished: "
    << std::endl << "lowestSquarkMass = " << lowestSquarkMass
    << std::endl << "highestSquarkMass = " << highestSquarkMass
    << std::endl << "lowestGluinoMass = " << lowestGluinoMass
    << std::endl << "highestGluinoMass = " << highestGluinoMass
    << std::endl << "lowNeutralinoMassRatio = " << lowNeutralinoMassRatio
    << std::endl << "middleNeutralinoMassRatio = "
    << middleNeutralinoMassRatio
    << std::endl << "highNeutralinoMassRatio = "
    << highNeutralinoMassRatio
    << std::endl << "value at ( " << lowestSquarkMass << ", "
    << lowestGluinoMass << ", "
    << ( lowNeutralinoMassRatio * lowestSquarkMass ) << ", "
    << ( highNeutralinoMassRatio * lowestSquarkMass ) << ", 0 ) is "
    << valueAt( lowestSquarkMass,
                lowestGluinoMass,
                ( lowNeutralinoMassRatio * lowestSquarkMass ),
                ( highNeutralinoMassRatio * lowestSquarkMass ),
                0 )
    << std::endl << "value at ( "
    << ( 0.5 * ( lowestSquarkMass + highestSquarkMass ) + 0.2 ) << ", "
    << ( 0.5 * ( lowestGluinoMass + highestGluinoMass ) + 0.3 ) << ", "
    << ( ( lowNeutralinoMassRatio + 0.1 )
         * ( 0.5 * ( lowestSquarkMass + highestSquarkMass ) + 0.2 ) )
    << ", "
    << ( ( middleNeutralinoMassRatio + 0.1 )
         * ( 0.5 * ( lowestSquarkMass + highestSquarkMass ) + 0.2 ) )
    << ", 0 ) is "
    << valueAt( ( 0.5 * ( lowestSquarkMass + highestSquarkMass ) + 0.2 ),
                ( 0.5 * ( lowestGluinoMass + highestGluinoMass ) + 0.3 ),
                ( ( lowNeutralinoMassRatio + 0.1 )
                  * ( 0.5 * ( lowestSquarkMass + highestSquarkMass ) + 0.2 ) ),
                ( ( middleNeutralinoMassRatio + 0.1 )
                  * ( 0.5 * ( lowestSquarkMass + highestSquarkMass ) + 0.2 ) ),
                0 )
    << std::endl << "printing out to test.dat";
    std::cout << std::endl;
    std::ofstream testOutput;
    testOutput.open( "test.dat" );
    double outputSquarkMass;
    double outputGluinoMass;
    double outputLighterScoloredMass;
    for( int squarkCounter( 0 );
         acceptanceValues.size() > squarkCounter;
         ++squarkCounter )
    {
      outputSquarkMass
      = ( lowestSquarkMass + ( scoloredMassStepSize * squarkCounter ) );
      for( int gluinoCounter( 0 );
           acceptanceValues.at( squarkCounter )->size() > gluinoCounter;
           ++gluinoCounter )
      {
        outputGluinoMass
        = ( lowestGluinoMass + ( scoloredMassStepSize * gluinoCounter ) );
        if( outputGluinoMass < outputSquarkMass )
        {
          outputLighterScoloredMass = outputGluinoMass;
        }
        else
        {
          outputLighterScoloredMass = outputSquarkMass;
        }
        for( int neutralinoCounter( 0 );
             6 > neutralinoCounter;
             ++neutralinoCounter )
        {
          testOutput
          << " " << outputSquarkMass
          << " " << outputGluinoMass;
          if( 0 == neutralinoCounter )
          {
            testOutput
            << " " << ( lowNeutralinoMassRatio * outputLighterScoloredMass )
            << " " << ( highNeutralinoMassRatio * outputLighterScoloredMass );
          }
          else if( 1 == neutralinoCounter )
          {
            testOutput
            << " " << ( lowNeutralinoMassRatio * outputLighterScoloredMass )
            << " "
            << ( middleNeutralinoMassRatio * outputLighterScoloredMass );
          }
          else if( 2 == neutralinoCounter )
          {
            testOutput
            << " " << ( lowNeutralinoMassRatio * outputLighterScoloredMass )
            << " " << ( lowNeutralinoMassRatio * outputLighterScoloredMass );
          }
          else if( 3 == neutralinoCounter )
          {
            testOutput
            << " " << ( middleNeutralinoMassRatio * outputLighterScoloredMass )
            << " " << ( highNeutralinoMassRatio * outputLighterScoloredMass );
          }
          else if( 4 == neutralinoCounter )
          {
            testOutput
            << " " << ( middleNeutralinoMassRatio * outputLighterScoloredMass )
            << " "
            << ( middleNeutralinoMassRatio * outputLighterScoloredMass );
          }
          else if( 5 == neutralinoCounter )
          {
            testOutput
            << " " << ( highNeutralinoMassRatio * outputLighterScoloredMass )
            << " " << ( highNeutralinoMassRatio * outputLighterScoloredMass );
          }
          for( std::vector< double >::const_iterator
               acceptanceIterator( acceptanceValues.at( squarkCounter )->at(
                                                           gluinoCounter )->at(
                                                neutralinoCounter )->begin() );
               acceptanceValues.at( squarkCounter )->at(
                                gluinoCounter )->at( neutralinoCounter )->end()
               > acceptanceIterator;
               ++acceptanceIterator )
          {
            testOutput << " " << *acceptanceIterator;
          }
          testOutput << std::endl;
        }
      }
    }
    testOutput.close();**/
  }

  acceptanceGrid::~acceptanceGrid()
  {
    // delete all the vectors held as pointers by the vector of vectors of
    // vectors:
    for( std::vector< std::vector< std::vector< std::vector< double
                                                           >* >* >* >::iterator
         gluinoVectorIterator( acceptanceValues.begin() );
         acceptanceValues.end() > gluinoVectorIterator;
         ++gluinoVectorIterator )
    {
      for( std::vector< std::vector< std::vector< double >* >* >::iterator
           neutralinoVectorIterator( (*gluinoVectorIterator)->begin() );
           (*gluinoVectorIterator)->end() > neutralinoVectorIterator;
           ++neutralinoVectorIterator )
      {
        for( std::vector< std::vector< double >* >::iterator
             deletionIterator( (*neutralinoVectorIterator)->begin() );
             (*neutralinoVectorIterator)->end() > deletionIterator;
             ++deletionIterator )
        {
          delete *deletionIterator;
        }
        delete *neutralinoVectorIterator;
      }
      delete *gluinoVectorIterator;
    }
  }


  bool
  acceptanceGrid::pointIsOnGrid( double const squarkMass,
                                 double const gluinoMass )
  /* this looks to see if the squark & gluino masses are on the grid, & if
   * so, sets lowerLeftVectorOfVectors, lowerRightVectorOfVectors,
   * upperRightVectorOfVectors, upperLeftVectorOfVectors, squarkMassFraction,
   * & gluinoMassFraction appropriately, returning true. otherwise, it
   * returns false.
   */
  {
    bool returnBool( false );
    if( ( 0.0 < scoloredMassStepSize )
        &&
        ( squarkMass >= lowestSquarkMass )
        &&
        ( gluinoMass >= lowestGluinoMass ) )
    {
      squarkFraction
      = ( ( squarkMass - lowestSquarkMass ) / scoloredMassStepSize );
      unsigned int squarkElement( (unsigned int)squarkFraction );
      if( ( acceptanceValues.size() == ( squarkElement + 1 ) )
          &&
          ( 0.0 == squarkFraction ) )
      {
        --squarkElement;
        squarkMassFraction = 1.0;
      }
      if( acceptanceValues.size() > ( squarkElement + 1 ) )
        // if the x co-ordinate is less than its maximal grid value...
      {
        gluinoFraction
        = ( ( gluinoMass - lowestGluinoMass ) / scoloredMassStepSize );
        unsigned int gluinoElement( (unsigned int)gluinoFraction );
        if( ( acceptanceValues.at( squarkElement ) == ( gluinoElement + 1 ) )
            &&
            ( 0.0 == gluinoFraction ) )
        {
          --gluinoElement;
          gluinoFraction = 1.0;
        }
        if( ( acceptanceValues.at( squarkElement )->size()
              > ( gluinoElement + 1 ) )
            &&
            ( acceptanceValues.at( squarkElement + 1 )->size()
              > ( gluinoElement + 1 ) ) )
          // if the y co-ordinate is less than its maximal grid value (on both
          // sides of the grid square)...
        {
          lowerLeftVectorOfVectors
          = acceptanceValues.at( squarkElement )->at( gluinoElement );
          lowerRightVectorOfVectors
          = acceptanceValues.at( squarkElement + 1 )->at( gluinoElement );
          upperRightVectorOfVectors
          = acceptanceValues.at( squarkElement + 1 )->at( gluinoElement + 1 );
          upperLeftVectorOfVectors
          = acceptanceValues.at( squarkElement )->at( gluinoElement + 1 );
          returnBool = true;
        }
        else
          // otherwise, it's off the grid:
        {
          if( inputShortcut->isVerbose() )
          {
            std::cout
            << std::endl
            << "LHC-FASER::warning! acceptance requested for a gluino mass"
            << " ( " << gluinoMass << " ) larger than allowed for ( "
            << highestGluinoMass << " ) in the lookup table!";
            std::cout << std::endl;
          }
        }

      }
      else
        // otherwise, it's off the grid:
      {
        if( inputShortcut->isVerbose() )
        {
          std::cout
          << std::endl
          << "LHC-FASER::warning! acceptance requested for a squark mass"
          << " ( " << squarkMass << " ) larger than allowed for ( "
          << highestSquarkMass << " ) in the lookup table!";
          std::cout << std::endl;
        }
      }
    }
    else
      // otherwise, it's off the grid:
    {
      if( inputShortcut->isVerbose() )
      {
        std::cout
        << std::endl
        << "LHC-FASER::warning! acceptance requested for a squark mass"
        << " ( " << squarkMass << " ) or gluino mass ( "
        << gluinoMass << " ) smaller than allowed for ( "
        << lowestSquarkMass << " or " << lowestGluinoMass << " respectively)"
        << " in the lookup table (or the table is malformed with a negative"
        << " grid size ( " << scoloredMassStepSize << " ))!";
        std::cout << std::endl;
      }
    }
    return returnBool;
  }

}  // end of LHC_FASER namespace.

