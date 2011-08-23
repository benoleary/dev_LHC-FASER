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
 *      LHC-FASER also requires grids of lookup values. These should also be
 *      found in a subdirectory included with this package.
 */

#include "LHC-FASER_base_kinematics_stuff.hpp"

namespace LHC_FASER
{
  acceptanceCutSet::acceptanceCutSet() :
    beamEnergy( -1.0 ),
    primaryLeptonCut( -1.0 ),
    secondaryLeptonCut( -1.0 ),
    jetCut( -1.0 ),
    excludedStandardModelProducts( NULL )
  // this constructor sets the values to -1.0, which is the default "unset"
  // value.
  {
    // just an initialization list.
  }

  acceptanceCutSet::acceptanceCutSet( acceptanceCutSet* const copyPointer ) :
    beamEnergy( copyPointer->beamEnergy ),
    primaryLeptonCut( copyPointer->primaryLeptonCut ),
    secondaryLeptonCut( copyPointer->secondaryLeptonCut ),
    jetCut( copyPointer->jetCut ),
    excludedStandardModelProducts( copyPointer->excludedStandardModelProducts )
  // this constructor copies the values from a given acceptanceCutSet.
  {
    // just an initialization list.
  }

  acceptanceCutSet::~acceptanceCutSet()
  {
    // does nothing.
  }



  acceptanceGrid::acceptanceGrid( std::string const* const gridFileLocation )
  {
    // read in the grid:

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: acceptanceGrid::acceptanceGrid( " << *grid_file_location
    << " ) called";
    std::cout << std::endl;**/

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
      lowestSquarkMass( 0.0 );
      lowestGluinoMass( 0.0 );
      double currentValue;
      // the value currently being read in.
      bool gridSizeStillUnknown( true );
      // we start out not knowing the step size for the grid (but we do
      // assume that it is a square grid).
      std::vector< double >* currentAcceptanceVector( NULL );
      // this holds all the acceptance values for a grid point.
      std::vector< std::vector< double >* >* currentNeutralinoVector( NULL );
      // this holds all the values for a set of points with the same
      // squark & gluino masses, but varying neutralino masses.
      std::vector< std::vector< std::vector< double >* >* >*
      currentGluinoVector( NULL );
      // this holds all the values for a set of points with the same
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
        currentAcceptanceVector = new std::vector< double >();
        while( !(inputLineAsStream.eof()) )
        {
          inputLineAsStream >> currentValue;
          currentAcceptanceVector->push_back( currentValue );
        }
        if( currentSquarkMass > lastSquarkMass )
          // if the squark mass has changed, we need a new vector for the
          // values for varying gluino masses.
        {
          // 1st record the last set of points with the same squark mass:
          if( NULL != currentGluinoVector )
            // if we have a vector to add...
          {
            currentGluinoVector->push_back( currentNeutralinoVector );
            values.push_back( currentGluinoVector );

            // debugging:
            /**std::cout
            << std::endl
            << "pushed back a neutralino vector of size "
            << current_neutralino_vector->size()
            << " and a gluino vector of size "
            << current_gluino_vector->size();**/
          }

          // prepare a new vector for the new squark mass:
          currentGluinoVector
          = new std::vector< std::vector< std::vector< double >* >* >;
          // prepare a new vector for the new gluino mass:
          currentNeutralinoVector = new std::vector< std::vector< double >* >;
        }
        else if( currentGluinoMass > lastGluinoMass )
          /* if the gluino mass has changed (but the squark mass has not
           * changed, we need a new vector for the values for varying
           * neutralino masses.
           */
        {
          // 1st record the last set of points with the same gluino mass:
          currentGluinoVector->push_back( currentNeutralinoVector );

          // debugging:
          /**std::cout
                << std::endl
                << "pushed back a neutralino vector of size "
                << current_neutralino_vector->size();**/

          // prepare a new vector for the new gluino mass:
          currentNeutralinoVector = new std::vector< std::vector< double >* >;
        }

        if( !gridSizeStillUnknown
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
        }
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
          acceptanceColumns = currentAcceptanceVector->size();
          // first we take note of the 1st point, which should have the
          // lowest co-ordinate values:
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
          double lowestNeutralinoMass = lowestSquarkMass;
          double middleNeutralinoMass = lowestSquarkMass;
          double highestNeutralinoMass = lowestSquarkMass;
          for( std::vector< double >::const_iterator
               massIterator = neutralinoMasses.begin();
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
            double lighterScoloredMass = lastSquarkMass;
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
      acceptanceColumns = currentAcceptanceVector->size();
      currentGluinoVector->push_back( currentNeutralinoVector );
      values.push_back( currentGluinoVector );

      // debugging:
      /**std::cout
        << std::endl
        << "pushed back last current_neutralino_vector of size "
        << current_neutralino_vector->size()
        << " & last current_gluino_vector of size "
        << current_gluino_vector->size();
        std::cout << std::endl;**/

      // the last masses should be the maximal values:
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
     * ( 0.5 * ( lowestSquarkMass + highestSquarkMass )
                       + 0.2 ) ),
                 ( ( middleNeutralinoMassRatio + 0.1 )
     * ( 0.5 * ( lowestSquarkMass + highestSquarkMass )
                       + 0.2 ) ),
                 0 )
    << std::endl << "printing out to test.dat";
    std::cout << std::endl;
    std::ofstream test_output;
    test_output.open( "test.dat" );
    double output_squark_mass;
    double output_gluino_mass;
    double output_lighter_scolored_mass;
    for( int squark_counter = 0;
         values.size() > squark_counter;
         ++squark_counter )
      {

        output_squark_mass
        = ( lowestSquarkMass
            + ( scoloredMassStepSize * squark_counter ) );

        for( int gluino_counter = 0;
            values.at( squark_counter )->size() > gluino_counter;
            ++gluino_counter )
          {

            output_gluino_mass
            = ( lowestGluinoMass
                + ( scoloredMassStepSize * gluino_counter ) );

            if( output_gluino_mass < output_squark_mass )
              {

                output_lighter_scolored_mass = output_gluino_mass;

              }
            else
              {

                output_lighter_scolored_mass = output_squark_mass;

              }

            for( int neutralino_counter = 0;
                 6 > neutralino_counter;
                 ++neutralino_counter )
              {

                test_output
                << " " << output_squark_mass
                << " " << output_gluino_mass;
                if( 0 == neutralino_counter )
                  {

                    test_output
                    << " " << ( lowNeutralinoMassRatio
     * output_lighter_scolored_mass )
                    << " " << ( highNeutralinoMassRatio
     * output_lighter_scolored_mass );

                  }
                else if( 1 == neutralino_counter )
                  {

                    test_output
                    << " " << ( lowNeutralinoMassRatio
     * output_lighter_scolored_mass )
                    << " " << ( middleNeutralinoMassRatio
     * output_lighter_scolored_mass );

                  }
                else if( 2 == neutralino_counter )
                  {

                    test_output
                    << " " << ( lowNeutralinoMassRatio
     * output_lighter_scolored_mass )
                    << " " << ( lowNeutralinoMassRatio
     * output_lighter_scolored_mass );
                  }
                else if( 3 == neutralino_counter )
                  {

                    test_output
                    << " " << ( middleNeutralinoMassRatio
     * output_lighter_scolored_mass )
                    << " " << ( highNeutralinoMassRatio
     * output_lighter_scolored_mass );

                  }
                else if( 4 == neutralino_counter )
                  {

                    test_output
                    << " " << ( middleNeutralinoMassRatio
     * output_lighter_scolored_mass )
                    << " " << ( middleNeutralinoMassRatio
     * output_lighter_scolored_mass );

                  }
                else if( 5 == neutralino_counter )
                  {

                    test_output
                    << " " << ( highNeutralinoMassRatio
     * output_lighter_scolored_mass )
                    << " " << ( highNeutralinoMassRatio
     * output_lighter_scolored_mass );

                  }

                for( std::vector< double >::const_iterator
                     acceptance_iterator
                     = values.at( squark_counter )->at( gluino_counter
                                        )->at( neutralino_counter )->begin();
                    values.at( squark_counter )->at( gluino_counter
                                             )->at( neutralino_counter )->end()
                    > acceptance_iterator;
                    ++acceptance_iterator )
                  {

                    test_output << " " << *acceptance_iterator;

                  }

                test_output << std::endl;

              }

          }

      }
    test_output.close();**/
  }

  acceptanceGrid::~acceptanceGrid()
  {
    // delete all the vectors held as pointers by the vector of vectors of
    // vectors:
    for( std::vector< std::vector< std::vector< std::vector< double
                                                           >* >* >* >::iterator
         gluinoVectorIterator( values.begin() );
         values.end() > gluinoVectorIterator;
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

  double
  acceptanceGrid::vectorElementAt( double const squarkMass,
                                   double const gluinoMass,
                                   int const neutralinoElement,
                                   int const acceptanceElement )
  const
  /* this takes the values for the element of the acceptance vector in the
   * specified element of the vector for neutralino masses requested for the
   * nearest points based on squark & gluino masses & interpolates to the given
   * masses.
   */
  {
    if( ( 0.0 < scoloredMassStepSize )
        &&
        ( squarkMass >= lowestSquarkMass )
        &&
        ( gluinoMass >= lowestGluinoMass ) )
    {
      double
      xSteps( ( ( squarkMass - lowestSquarkMass ) / scoloredMassStepSize ) );
      int lowerLeftX ( (int)xSteps );

      if( values.size() > ( lowerLeftX + 1 ) )
        // if the x co-ordinate is less than its maximal grid value...
      {
        double ySteps( ( ( gluinoMass - lowestGluinoMass )
                         / scoloredMassStepSize ) );
        int lowerLeftY( (int)ySteps );
        if( values.at( lowerLeftX )->size() > ( lowerLeftY + 1 ) )
          // if the y co-ordinate is less than its maximal grid value...
        {
          return lhcFaserGlobal::squareBilinearInterpolation(
                                                       ( xSteps - lowerLeftX ),
                                                       ( ySteps - lowerLeftY ),
                                        values.at( lowerLeftX )->at( lowerLeftY
                                                       )->at( neutralinoElement
                                                    )->at( acceptanceElement ),
                                values.at( ( lowerLeftX + 1 ) )->at( lowerLeftY
                                                       )->at( neutralinoElement
                                                    )->at( acceptanceElement ),
                        values.at( ( lowerLeftX + 1 ) )->at( ( lowerLeftY + 1 )
                                                       )->at( neutralinoElement
                                                    )->at( acceptanceElement ),
                                values.at( lowerLeftX )->at( ( lowerLeftY + 1 )
                                                       )->at( neutralinoElement
                                                  )->at( acceptanceElement ) );
        }
        else if( ( values.at( lowerLeftX )->size()
                   == ( lowerLeftY + 1 ) )
                 &&
                 ( (double)lowerLeftY == ySteps ) )
          // otherwise, if it's on the maximal y edge...
        {
          return lhcFaserGlobal::unitLinearInterpolation(
                                                       ( xSteps - lowerLeftX ),
                                        values.at( lowerLeftX )->at( lowerLeftY
                                                       )->at( neutralinoElement
                                                    )->at( acceptanceElement ),
                                values.at( ( lowerLeftX + 1 ) )->at( lowerLeftY
                                                       )->at( neutralinoElement
                                                  )->at( acceptanceElement ) );
        }
        else
          // otherwise, it's off the grid:
        {
          std::cout
          << std::endl
          << "LHC-FASER::warning! acceptance requested for a gluino mass"
          << " ( " << squarkMass << " ) larger than allowed for ( "
          << highestGluinoMass << " ) in the lookup table!";
          std::cout << std::endl;
          return CppSLHA::CppSLHA_global::really_wrong_value;
        }
      }
      else if( ( values.size() == ( lowerLeftX + 1 ) )
               &&
               ( (double)lowerLeftX == xSteps ) )
        // otherwise, if it's on the maximal x edge...
      {

        double ySteps( ( ( gluinoMass - lowestGluinoMass )
                         / scoloredMassStepSize ) );
        int lowerLeftY( (int)ySteps );

        if( values.at( lowerLeftX )->size() > ( lowerLeftY + 1 ) )
          // if the y co-ordinate is less than its maximal grid value...
        {
          return lhcFaserGlobal::unitLinearInterpolation(
                                                       ( ySteps - lowerLeftY ),
                                        values.at( lowerLeftX )->at( lowerLeftY
                                                       )->at( neutralinoElement
                                                    )->at( acceptanceElement ),
                                values.at( lowerLeftX )->at( ( lowerLeftY + 1 )
                                                       )->at( neutralinoElement
                                                  )->at( acceptanceElement ) );
        }
        else if( ( values.at( lowerLeftX )->size()
                   == ( lowerLeftY + 1 ) )
                 &&
                 ( (double)lowerLeftY == ySteps ) )
          // otherwise, if it's on the maximal x & y corner...
        {
          return values.at( lowerLeftX )->at( lowerLeftY )->at(
                                  neutralinoElement )->at( acceptanceElement );
        }
        else
          // otherwise, it's off the grid:
        {
          std::cout
          << std::endl
          << "LHC-FASER::warning! acceptance requested for a gluino mass"
          << " ( " << squarkMass << " ) larger than allowed for ( "
          << highestGluinoMass << " ) in the lookup table!";
          std::cout << std::endl;
          return CppSLHA::CppSLHA_global::really_wrong_value;
        }
      }
      else
        // otherwise, it's off the grid:
      {
        std::cout
        << std::endl
        << "LHC-FASER::warning! acceptance requested for a squark mass ( "
        << squarkMass << " ) larger than allowed for ( "
        << highestSquarkMass << " ) in the lookup table!";
        std::cout << std::endl;
        return CppSLHA::CppSLHA_global::really_wrong_value;
      }
    }
    else
      // otherwise, it's off the grid (or the grid step size was not positive,
      // which is still a problem...):
    {
      std::cout
      << std::endl
      << "LHC-FASER::warning! acceptance requested for a point outside the"
      << " lookup table!"
      << std::endl << "squark mass:"
      << "lower limit = " << lowestSquarkMass
      << ", requested = " << squarkMass
      << std::endl << "gluino mass:"
      << "lower limit = " << lowestGluinoMass
      << ", requested = " << gluinoMass
      << " (or possibly the lookup table itself is invalid: its step size,"
      << " which should be positive, is " << scoloredMassStepSize << ")";
      std::cout << std::endl;
      return CppSLHA::CppSLHA_global::really_wrong_value;
    }
  }

  double
  acceptanceGrid::valueAt( double const squarkMass,
                           double const gluinoMass,
                           double const firstNeutralinoMass,
                           double const secondNeutralinoMass,
                           int const acceptanceElement,
                           bool const heavyNeutralinoEdgeIsLighterScoloredMass,
                           bool const heavyNeutralinoAreaIsConstant )
  const
  /* this finds the grid square which the given point is in, & then uses
   * lhcFaserGlobal::squareBilinearInterpolation to get an interpolated
   * value, assuming that the heavy neutralino edge goes to 0.0 as the
   * heavier neutralino mass approaches the lighter scolored mass unless
   * heavy_neutralino_edge_is_lighter_scolored_mass is true, in which case
   * it interpolates to the lighter scolored mass, or unless
   * heavy_neutralino_area_is_constant is true, in which case it interpolates
   * to the same value as the grid points with the heaviest neutralino mass.
   */
  {
    // debugging:
    /**std::cout
    << std::endl
    << "debugging: acceptanceGrid::valueAt( " << squark_mass
    << ", " << gluino_mass << ", " << first_neutralino_mass << ", "
    << second_neutralino_mass << ", " << acceptance_element << ", "
    << heavy_neutralino_edge_is_lighter_scolored_mass
    << " ) called.";
    std::cout << std::endl;**/

    if( ( acceptanceElement < acceptanceColumns )
        &&
        ( 0.0 < scoloredMassStepSize )
        &&
        ( squarkMass >= lowestSquarkMass )
        &&
        ( gluinoMass >= lowestGluinoMass ) )
    {
      // 1st we find out which are the lighter sparticles:
      double lighterScolored( squarkMass );
      if( gluinoMass < squarkMass )
      {
        lighterScolored = gluinoMass;
      }
      double lighterNeutralino( firstNeutralinoMass );
      double heavierNeutralino( secondNeutralinoMass );
      if( secondNeutralinoMass < firstNeutralinoMass )
      {
        lighterNeutralino = secondNeutralinoMass;
        heavierNeutralino = firstNeutralinoMass;
      }
      lighterNeutralino *= ( 1.0 / lighterScolored );
      heavierNeutralino *= ( 1.0 / lighterScolored );
      // now lighter_neutralino & heavier_neutralino are mass ratios.

      // debugging:
      /**std::cout
        << std::endl
        << "lighter_scolored = " << lighter_scolored
        << std::endl
        << "lighter_neutralino = " << lighter_neutralino
        << std::endl
        << "heavier_neutralino = " << heavier_neutralino
        << std::endl
        << "lowNeutralinoMassRatio = " << lowNeutralinoMassRatio
        << std::endl
        << "middleNeutralinoMassRatio = " << middleNeutralinoMassRatio
        << std::endl
        << "highNeutralinoMassRatio = " << highNeutralinoMassRatio;
        std::cout << std::endl;**/

      if( ( 1.0 > lighterNeutralino )
          &&
          ( 1.0 > heavierNeutralino ) )
      {
        /* now we work out which elements of the values vector of vectors
         * of vectors to use.
         * N.B.! this bit is HIGHLY format-dependent, since the values are
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
        if( lowNeutralinoMassRatio > lighterNeutralino )
        {
          if( lowNeutralinoMassRatio > heavierNeutralino )
          {
            return vectorElementAt( squarkMass,
                                    gluinoMass,
                                    2,
                                    acceptanceElement );
            /*lower_left_point = 3;
            lower_right_point = 3;
            upper_right_point = 3;
            upper_left_point = 3;*/
          }
          else if( middleNeutralinoMassRatio > heavierNeutralino )
          {
            return lhcFaserGlobal::unitLinearInterpolation(
                               ( ( heavierNeutralino - lowNeutralinoMassRatio )
                    / ( middleNeutralinoMassRatio - lowNeutralinoMassRatio ) ),
                                                   vectorElementAt( squarkMass,
                                                                    gluinoMass,
                                                                    2,
                                                           acceptanceElement ),
                                                   vectorElementAt( squarkMass,
                                                                    gluinoMass,
                                                                    1,
                                                         acceptanceElement ) );
            /*lower_left_point = 3;
            lower_right_point = 3;
            upper_right_point = 2;
            upper_left_point = 2;*/
          }
          else if( highNeutralinoMassRatio > heavierNeutralino )
          {
            return lhcFaserGlobal::unitLinearInterpolation(
                            ( ( heavierNeutralino - middleNeutralinoMassRatio )
                   / ( highNeutralinoMassRatio - middleNeutralinoMassRatio ) ),
                                                   vectorElementAt( squarkMass,
                                                                    gluinoMass,
                                                                    1,
                                                           acceptanceElement ),
                                                   vectorElementAt( squarkMass,
                                                                    gluinoMass,
                                                                    0,
                                                         acceptanceElement ) );

            /*lower_left_point = 2;
            lower_right_point = 2;
            upper_right_point = 1;
            upper_left_point = 1;*/
          }
          else
          {
            if( heavyNeutralinoAreaIsConstant )
            {
              return vectorElementAt( squarkMass,
                                      gluinoMass,
                                      1,
                                      acceptanceElement );
            }
            else
            {
              double heavyEdgeValue( 0.0 );
              if( heavyNeutralinoEdgeIsLighterScoloredMass )
              {
                if( gluinoMass < squarkMass )
                {
                  heavyEdgeValue = gluinoMass;
                }
                else
                {
                  heavyEdgeValue = squarkMass;
                }
              }
            return lhcFaserGlobal::unitLinearInterpolation(
                              ( ( heavierNeutralino - highNeutralinoMassRatio )
                                         / ( 1.0 - highNeutralinoMassRatio ) ),
                                                   vectorElementAt( squarkMass,
                                                                    gluinoMass,
                                                                    1,
                                                           acceptanceElement ),
                                                              heavyEdgeValue );
            }
            /*lower_left_point = 1;
            lower_right_point = 1;
            upper_right_point = 0;
            upper_left_point = 0;*/
          }
        }
        else if( middleNeutralinoMassRatio > lighterNeutralino )
        {
          if( middleNeutralinoMassRatio > heavierNeutralino )
          {
            double lowerLeftValue( vectorElementAt( squarkMass,
                                                    gluinoMass,
                                                    2,
                                                    acceptanceElement ) );
            double upperRightValue( vectorElementAt( squarkMass,
                                                     gluinoMass,
                                                     4,
                                                     acceptanceElement ) );
            double upperLeftValue( vectorElementAt( squarkMass,
                                                    gluinoMass,
                                                    1,
                                                    acceptanceElement ) );
            double lowerRightValue( ( lowerLeftValue
                                      + upperRightValue
                                      - upperLeftValue ) );

            return lhcFaserGlobal::squareBilinearInterpolation(
                               ( ( lighterNeutralino - lowNeutralinoMassRatio )
                    / ( middleNeutralinoMassRatio - lowNeutralinoMassRatio ) ),
                        ( ( heavierNeutralino - lowNeutralinoMassRatio )
                    / ( middleNeutralinoMassRatio - lowNeutralinoMassRatio ) ),
                                                                lowerLeftValue,
                                                               lowerRightValue,
                                                               upperRightValue,
                                                              upperLeftValue );
            /*lower_left_point = 3;
            lower_right_point = -2;
            upper_right_point = 5;
            upper_left_point = 2;*/
          }
          else if( highNeutralinoMassRatio > heavierNeutralino )
          {
            return lhcFaserGlobal::squareBilinearInterpolation(
                               ( ( lighterNeutralino - lowNeutralinoMassRatio )
                    / ( middleNeutralinoMassRatio - lowNeutralinoMassRatio ) ),
                            ( ( heavierNeutralino - middleNeutralinoMassRatio )
                   / ( highNeutralinoMassRatio - middleNeutralinoMassRatio ) ),
                                                   vectorElementAt( squarkMass,
                                                                    gluinoMass,
                                                                    1,
                                                           acceptanceElement ),
                                                   vectorElementAt( squarkMass,
                                                                    gluinoMass,
                                                                    4,
                                                           acceptanceElement ),
                                                   vectorElementAt( squarkMass,
                                                                    gluinoMass,
                                                                    3,
                                                           acceptanceElement ),
                                                   vectorElementAt( squarkMass,
                                                                    gluinoMass,
                                                                    0,
                                                         acceptanceElement ) );
            /*lower_left_point = 2;
            lower_right_point = 5;
            upper_right_point = 4;
            upper_left_point = 1;*/
          }
          else
          {
            double heavyEdgeValue( 0.0 );
            if( heavyNeutralinoEdgeIsLighterScoloredMass )
            {
              if( gluinoMass < squarkMass )
              {
                heavyEdgeValue = gluinoMass;
              }
              else
              {
                heavyEdgeValue = squarkMass;
              }
            }
            else if( heavyNeutralinoAreaIsConstant )
            {
              heavyEdgeValue = vectorElementAt( squarkMass,
                                                gluinoMass,
                                                3,
                                                acceptanceElement );
            }
            return lhcFaserGlobal::squareBilinearInterpolation(
                               ( ( lighterNeutralino - lowNeutralinoMassRatio )
                    / ( middleNeutralinoMassRatio - lowNeutralinoMassRatio ) ),
                              ( ( heavierNeutralino - highNeutralinoMassRatio )
                                / ( 1.0 - highNeutralinoMassRatio ) ),
                                                   vectorElementAt( squarkMass,
                                                                    gluinoMass,
                                                                    0,
                                                           acceptanceElement ),
                                                   vectorElementAt( squarkMass,
                                                                    gluinoMass,
                                                                    3,
                                                           acceptanceElement ),
                                                                heavyEdgeValue,
                                                              heavyEdgeValue );
            /*lower_left_point = 1;
            lower_right_point = 4;
            upper_right_point = 0;
            upper_left_point = 0;*/
          }
        }
        else if( highNeutralinoMassRatio > lighterNeutralino )
        {
          if( highNeutralinoMassRatio > heavierNeutralino )
          {
            // debugging:
            /**std::cout
                    << std::endl
                    << "so, highNeutralinoMassRatio > both > middle"
                    << std::endl;
                    std::cout << std::endl;**/

            double lowerLeftValue( vectorElementAt( squarkMass,
                                                    gluinoMass,
                                                    4,
                                                    acceptanceElement ) );

            // debugging:
            /**std::cout
                    << std::endl
                    << "lower_left_value = " << lower_left_value;
                    std::cout << std::endl;**/

            double upperRightValue( vectorElementAt( squarkMass,
                                                     gluinoMass,
                                                     5,
                                                     acceptanceElement ) );

            // debugging:
            /**std::cout
                    << std::endl
                    << "upper_right_value = " << upper_right_value;
                    std::cout << std::endl;**/

            double upperLeftValue( vectorElementAt( squarkMass,
                                                    gluinoMass,
                                                    3,
                                                    acceptanceElement ) );

            // debugging:
            /**std::cout
                    << std::endl
                    << "upper_left_value = " << upper_left_value;
                    std::cout << std::endl;**/

            double lowerRightValue( ( lowerLeftValue
                                      + upperRightValue
                                      - upperLeftValue ) );

            // debugging:
            /**std::cout
                    << std::endl
                    << "lower_right_value = " << lower_right_value;
                    std::cout << std::endl;**/

            return lhcFaserGlobal::squareBilinearInterpolation(
                            ( ( lighterNeutralino - middleNeutralinoMassRatio )
                   / ( highNeutralinoMassRatio - middleNeutralinoMassRatio ) ),
                            ( ( heavierNeutralino - middleNeutralinoMassRatio )
                   / ( highNeutralinoMassRatio - middleNeutralinoMassRatio ) ),
                                                                lowerLeftValue,
                                                               lowerRightValue,
                                                               upperRightValue,
                                                              upperLeftValue );
            /*lower_left_point = 5;
           lower_right_point = -4;
           upper_right_point = 6;
           upper_left_point = 4;*/
          }
          else
          {
            double heavyEdgeValue( 0.0 );
            if( heavyNeutralinoEdgeIsLighterScoloredMass )
            {
              if( gluinoMass < squarkMass )
              {
                heavyEdgeValue = gluinoMass;
              }
              else
              {
                heavyEdgeValue = squarkMass;
              }
            }
            else if( heavyNeutralinoAreaIsConstant )
            {
              heavyEdgeValue = vectorElementAt( squarkMass,
                                                gluinoMass,
                                                3,
                                                acceptanceElement );
            }
            return lhcFaserGlobal::squareBilinearInterpolation(
                            ( ( lighterNeutralino - middleNeutralinoMassRatio )
                   / ( highNeutralinoMassRatio - middleNeutralinoMassRatio ) ),
                              ( ( heavierNeutralino - highNeutralinoMassRatio )
                                / ( 1.0 - highNeutralinoMassRatio ) ),
                                                   vectorElementAt( squarkMass,
                                                                    gluinoMass,
                                                                    3,
                                                           acceptanceElement ),
                                                   vectorElementAt( squarkMass,
                                                                    gluinoMass,
                                                                    5,
                                                           acceptanceElement ),
                                                                heavyEdgeValue,
                                                              heavyEdgeValue );
            /*lower_left_point = 4;
            lower_right_point = 6;
            upper_right_point = 0;
            upper_left_point = 0;*/
          }
        }
        else
        {
          double heavyEdgeValue( 0.0 );
          if( heavyNeutralinoEdgeIsLighterScoloredMass )
          {
            if( gluinoMass < squarkMass )
            {
              heavyEdgeValue = gluinoMass;
            }
            else
            {
              heavyEdgeValue = squarkMass;
            }
          }
          else if( heavyNeutralinoAreaIsConstant )
          {
            heavyEdgeValue = vectorElementAt( squarkMass,
                                                gluinoMass,
                                                5,
                                                acceptanceElement );
          }

          return lhcFaserGlobal::unitLinearInterpolation(
                              ( ( lighterNeutralino - highNeutralinoMassRatio )
                                / ( 1.0 - highNeutralinoMassRatio ) ),
                                                   vectorElementAt( squarkMass,
                                                                    gluinoMass,
                                                                    5,
                                                           acceptanceElement ),
                                                              heavyEdgeValue );
          /*lower_left_point = 6;
          lower_right_point = 6;
          upper_right_point = 0;
          upper_left_point = 0;*/
        }  // at this point, the points have been assigned.
      }
      else
      {
        return CppSLHA::CppSLHA_global::really_wrong_value;
      }
    }
    else
    {
      return CppSLHA::CppSLHA_global::really_wrong_value;
    }
  }

  double
  acceptanceGrid::valueAt( double const squarkMass,
                           double const gluinoMass,
                           double const degenerateNeutralinoMass,
                           int const acceptanceElement,
                           bool const heavyNeutralinoEdgeIsLighterScoloredMass,
                           bool const heavyNeutralinoAreaIsConstant )
  const
  /* this finds the grid square which the given point is in, & then uses
   * lhcFaserGlobal::squareBilinearInterpolation to get an interpolated
   * value, assuming that the heavy neutralino edge goes to 0.0 as the
   * heavier neutralino mass approaches the lighter scolored mass unless
   * heavy_neutralino_edge_is_lighter_scolored_mass is true, in which case
   * it interpolates to the lighter scolored mass, or unless
   * heavy_neutralino_area_is_constant is true, in which case it interpolates
   * to the same value as the grid points with the heaviest neutralino mass.
   * N.B.: this version is just to save a little calculation for the
   * lepton acceptance for a cascade because of the approximation that the
   * kinematics for the lepton acceptance of 1 cascade is independent of that
   * of the other cascade. the improvement is probably utterly negligible,
   * but I want to write this function anyway...
   */
  {
    // debugging:
    /**std::cout
    << std::endl
    << "debugging: acceptanceGrid::valueAt( " << squark_mass
    << ", " << gluino_mass << ", " << degenerate_neutralino_mass << ", "
    << acceptance_element << ", "
    << heavy_neutralino_edge_is_lighter_scolored_mass
    << " ) called.";
    std::cout << std::endl;**/

    if( ( acceptanceElement < acceptanceColumns )
        &&
        ( 0.0 < scoloredMassStepSize )
        &&
        ( squarkMass >= lowestSquarkMass )
        &&
        ( gluinoMass >= lowestGluinoMass ) )
    {
      // 1st we find out which are the lighter sparticles:
      double lighterScolored( squarkMass );
      if( gluinoMass < squarkMass )
      {

        lighterScolored = gluinoMass;

      }
      double neutralinoRatio = ( degenerateNeutralinoMass / lighterScolored );

      // debugging:
      /**std::cout
        << std::endl
        << "lighter_scolored = " << lighter_scolored
        << std::endl
        << "neutralino_ratio = " << neutralino_ratio
        << std::endl
        << "lowNeutralinoMassRatio = " << lowNeutralinoMassRatio
        << std::endl
        << "middleNeutralinoMassRatio = " << middleNeutralinoMassRatio
        << std::endl
        << "highNeutralinoMassRatio = " << highNeutralinoMassRatio;
        std::cout << std::endl;**/

      if( 1.0 > neutralinoRatio )
      {

        /* now we work out which elements of the values vector of vectors
         * of vectors to use.
         * N.B.! this bit is HIGHLY format-dependent, since the values are
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
        if( lowNeutralinoMassRatio > neutralinoRatio )
        {
          return vectorElementAt( squarkMass,
                                  gluinoMass,
                                  2,
                                  acceptanceElement );

          /*lower_left_point = 3;
          lower_right_point = 3;
          upper_right_point = 3;
          upper_left_point = 3;*/
        }
        else if( middleNeutralinoMassRatio > neutralinoRatio )
        {
          double lowerLeftValue( vectorElementAt( squarkMass,
                                                  gluinoMass,
                                                  2,
                                                  acceptanceElement ) );
          double upperRightValue( vectorElementAt( squarkMass,
                                                   gluinoMass,
                                                   4,
                                                   acceptanceElement ) );
          return lhcFaserGlobal::unitLinearInterpolation(
                                 ( ( neutralinoRatio - lowNeutralinoMassRatio )
                    / ( middleNeutralinoMassRatio - lowNeutralinoMassRatio ) ),
                                                              lowerLeftValue,
                                                             upperRightValue );
          /*lower_left_point = 3;
          lower_right_point = -2;
          upper_right_point = 5;
          upper_left_point = 2;*/
        }
        else if( highNeutralinoMassRatio > neutralinoRatio )
        {
          double lowerLeftValue( vectorElementAt( squarkMass,
                                                  gluinoMass,
                                                  4,
                                                  acceptanceElement ) );
          double upperRightValue( vectorElementAt( squarkMass,
                                                   gluinoMass,
                                                   5,
                                                   acceptanceElement ) );

          return lhcFaserGlobal::unitLinearInterpolation(
                              ( ( neutralinoRatio - middleNeutralinoMassRatio )
                   / ( highNeutralinoMassRatio - middleNeutralinoMassRatio ) ),
                                                              lowerLeftValue,
                                                             upperRightValue );
          /*lower_left_point = 5;
          lower_right_point = -4;
          upper_right_point = 6;
          upper_left_point = 4;*/
        }
        else
        {
          if( heavyNeutralinoAreaIsConstant )
          {
            return vectorElementAt( squarkMass,
                                    gluinoMass,
                                    5,
                                    acceptanceElement );
          }
          else
          {
            double heavyEdgeValue = 0.0;
            if( heavyNeutralinoEdgeIsLighterScoloredMass )
            {
              if( gluinoMass < squarkMass )
              {
                heavyEdgeValue = gluinoMass;
              }
              else
              {
                heavyEdgeValue = squarkMass;
              }
            }
            return lhcFaserGlobal::unitLinearInterpolation(
                             ( ( lighter_neutralino - highNeutralinoMassRatio )
                               / ( 1.0 - highNeutralinoMassRatio ) ),
                                                   vectorElementAt( squarkMass,
                                                                    gluinoMass,
                                                                    5,
                                                           acceptanceElement ),
                                                              heavyEdgeValue );
          }
          /*lower_left_point = 6;
          lower_right_point = 6;
          upper_right_point = 0;
          upper_left_point = 0;*/
        }  // at this point, the points have been assigned.
      }
      else
      {
        return CppSLHA::CppSLHA_global::really_wrong_value;
      }
    }
    else
    {
      return CppSLHA::CppSLHA_global::really_wrong_value;
    }
  }
}  // end of LHC_FASER namespace.

