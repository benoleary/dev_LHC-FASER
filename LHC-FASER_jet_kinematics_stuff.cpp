/*
 * LHC-FASER_kinematics_stuff.hpp
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

#include "LHC-FASER_kinematics_stuff.hpp"

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



  singleLeptonCut::singleLeptonCut() :
    acceptanceCutSet()
  {

    // just an initialization.

  }

  singleLeptonCut::~singleLeptonCut()
  {

    // does nothing.

  }



  leptonAcceptanceParameterSet::leptonAcceptanceCut( input_handler* const shortcut,
                                         leptonAcceptanceGrid const* acceptanceTable,
                          CppSLHA::particle_property_set const* const scolored,
                                       colored_cascade const* const cascade ) :
    readied_for_new_point( shortcut->get_readier() ),
    effectiveSquarkMassHolder(),
    shortcut( shortcut ),
    acceptanceTable( acceptanceTable ),
    scolored( scolored ),
    lookupPair( scolored,
                true,
                scolored,
                false ),
    cascade( cascade ),
    binSize( defaultBinSize ),
    transverseMomentumCut( defaultTransverseMomentumCut ),
    effectiveSquarkMass( 0.0 ),
    pseudorapidityAcceptance( 0.0 )
  {
    acceptanceBins.push_back( 0.0 );
    // this is to ensure that there is at least 1 entry for the transverse
    // momentum cut.
  }

  leptonAcceptanceParameterSet::leptonAcceptanceCut( input_handler* const shortcut,
                                         leptonAcceptanceGrid const* acceptanceTable,
                          CppSLHA::particle_property_set const* const scolored,
                                          colored_cascade const* const cascade,
                                            double const binSize,
                                         double const transverseMomentumCut ) :
    readied_for_new_point( shortcut->get_readier() ),
    effectiveSquarkMassHolder(),
    shortcut( shortcut ),
    acceptanceTable( acceptanceTable ),
    scolored( scolored ),
    cascade( cascade ),
    binSize( binSize ),
    transverseMomentumCut( transverseMomentumCut ),
    effectiveSquarkMass( 0.0 ),
    pseudorapidityAcceptance( 0.0 )
  {
    acceptanceBins.push_back( 0.0 );
    // this is to ensure that there is at least 1 entry for the transverse
    // momentum cut.
  }

  leptonAcceptanceParameterSet::~leptonAcceptanceParameterSet()
  {

    // does nothing.

  }


  void
  leptonAcceptanceParameterSet::resetValues()
  // this interpolates values from acceptanceTable to set up acceptanceBins
  // for the given colored sparticle.
  {
    acceptanceBins.clear();
    acceptanceCounter = 0;
    currentAcceptance = acceptanceTable->getValue( lookupPair,
                                                    cascade,
                                                    cascade,
                                                    acceptanceCounter,
                                                    true,
                                                    false );
    // the 1st element of the vector is the effective squark mass, so
    // should be extrapolated to the lighter scolored mass.
    effectiveSquarkMass = currentAcceptance;
    ++acceptanceCounter;
    currentAcceptance = acceptanceTable->getValue( lookupPair,
                                                    cascade,
                                                    cascade,
                                                    acceptanceCounter,
                                                    false,
                                                    true );
    // the 2nd element of the vector is the pseudorapidity cut acceptance, so
    // should be extrapolated to what it is for the heavier given values.
    pseudorapidityAcceptance = currentAcceptance;

    while( 0.0 <= currentAcceptance )
    {
      // debugging:
      /**std::cout
      << std::endl
      << "debugging: ( 0.0 <= current_acceptance ) is true:"
      << " current_acceptance = " << current_acceptance
      << " for bin energy " << ( 2.0 * ++bin_counter);
      std::cout << std::endl;**/

      ++acceptanceCounter;
      currentAcceptance = acceptanceTable->getValue( scoloredPair,
                                                      cascade,
                                                      cascade,
                                                      acceptanceCounter,
                                                      false,
                                                      false );
      if( 0.0 <= currentAcceptance )
      {
        acceptanceBins.push_back( currentAcceptance );
        // the remaining elements of the vector are acceptances, so should
        // be extrapolated to zero.
      }
      // debugging:
      /**std::cout
      << std::endl
      << "now current_acceptance = " << current_acceptance;
      std::cout << std::endl;**/
    }
  }


  double
  leptonAcceptanceParameterSet::calculateAcceptanceAt( double const givenEnergy,
                                              double const givenCut )
  const
  /* this interpolates the values in acceptanceBins to the requested value,
   * or returns pseudorapidityAcceptance if it's lower, scaled to the given
   * value for the transverse momentum cut.
   */
  {
    // debugging:
    /**std::cout
    << std::endl
    << "debugging: leptonAcceptanceParameterSet::acceptanceAt( "
    << energy << ", " << cut << " ) called";
    std::cout << std::endl;**/

    returnValue = 0.0;

    if( ( 0.0 < givenEnergy )
        &&
        ( 0.0 < givenCut ) )
    {
      binFraction = ( givenEnergy * ( transverseMomentumCut
                                      / ( givenCut * binSize ) ) );
      // now binFraction is givenEnergy scaled to the given momentum cut,
      // in units of binSize. importantly, binFraction > 0.0 still.
      lowerBin = (int)binFraction;
      /* now lowerBin is 1 more than the lower bin for givenEnergy because
       * "bin -1" is the acceptance for 0.0 GeV, which is 0.0 by assumption.
       * "bin 0" is for 2.0 GeV by default, "bin 1" for 4.0 GeV, & so on.
       */
      binFraction -= (double)lowerBin;
      // now binFraction is the fraction of the bin along from lowerBin
      // that givenEnergy corresponds to.
      --lowerBin;
      // now lowerBin corresponds to the correct bin for givenEnergy.

      // debugging:
      /**std::cout
      << std::endl
      << "binFraction = " << binFraction
      << ", lowerBin = " << lowerBin
      << ", acceptanceBins.size() = " << acceptanceBins.size();
      std::cout << std::endl;**/

      if( -1 == lowerBin )
        // if we're in the region 0.0 GeV to binSize GeV...
      {
        // debugging:
        /**std::cout
        << std::endl
        << "( -1 == lowerBin ) is true, binFraction = "
        << binFraction << ", acceptanceBins.front() = "
        << acceptanceBins.front();
        std::cout << std::endl;**/

        returnValue = ( givenEnergy * acceptanceBins.front() );
        // we assume that the acceptance at 0.0 GeV is 0.0, so the linear
        // interpolation is simple.

      }
      else if( ( 0 <= lowerBin )
               &&
               ( acceptanceBins.size() > ( lowerBin + 1 ) ) )
        // if we're in a region between 2 bin entries...
      {
        // debugging:
        /**std::cout
        << std::endl
        << "( ( ( 0 <= lowerBin ) &&"
        << " ( acceptanceBins.size() > ( lowerBin + 1 ) ) ) is true,"
        << " acceptanceBins.at( lowerBin ) = "
        << acceptanceBins.at( lowerBin )
        << ", binFraction = "
        << binFraction
        << ", acceptanceBins.at( ( lowerBin + 1 ) ) = "
        << acceptanceBins.at( ( lowerBin + 1 ) );
        std::cout << std::endl;**/

        returnValue
        = ( acceptanceBins.at( lowerBin )
            + binFraction * ( acceptanceBins.at( ( lowerBin + 1 ) )
                              - acceptanceBins.at( lowerBin ) ) );
        // we return a linear interpolation.

      }
      else if( acceptanceBins.size() <= ( lowerBin + 1 ) )
        // if we're in a region beyond the bins...
      {
        // debugging:
        /**std::cout
        << std::endl
        << "( acceptanceBins.size() <= ( lowerBin + 1 ) ) is true,"
        << " acceptanceBins.back() = " << acceptanceBins.back();
        std::cout << std::endl;**/

        return_value = acceptanceBins.back();
        // we assume the acceptance is dominated by the pseudorapidity cut
        // by this point, so is constant out to any higher energy.

      }

      // this is a sanity check that I have decided is unnecessary:
      /*if( pseudorapidityAcceptance < returnValue )
      {

        returnValue = pseudorapidityAcceptance;

      }*/

    }
    return returnValue;
  }



  acceptanceGrid::acceptanceGrid(
                                  std::string const* const grid_file_location )
  {

    // read in the grid:

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: acceptanceGrid::acceptanceGrid( " << *grid_file_location
    << " ) called";
    std::cout << std::endl;**/

    if( NULL != grid_file_location )
      {

        // open the file:
        int const max_length_of_grid_file = 100000;
        CppSLHA::hash_commented_file_reader
        input_file_reader( grid_file_location,
                           max_length_of_grid_file,
                           false );
        /* a hundred thousand lines should be long enough for any grid file. if
         * readIn() is reading in more lines than that, it's probably because
         * of some file corruption.
         * also, if it's not enough, change it here.
         */

        std::string input_line_as_string;
        // prepare a string to parse the line.
        std::stringstream input_line_as_stream;
        // prepare a stringstream to parse the line.

        // the masses from the last line read & the currently-read line:
        double last_squark_mass = 0.0;
        double last_gluino_mass = 0.0;
        double current_squark_mass = 0.0;
        double current_gluino_mass = 0.0;
        std::vector< double > neutralino_masses;
        lowestSquarkMass = 0.0;
        lowestGluinoMass = 0.0;

        double current_value;
        // the value currently being read in.

        bool grid_size_still_unknown = true;
        // we start out not knowing the step size for the grid (but we do
        // assume that it is a square grid).

        std::vector< double >* current_acceptance_vector = NULL;
        // this holds all the acceptance values for a grid point.
        std::vector< std::vector< double >* >*
        current_neutralino_vector = NULL;
        // this holds all the values for a set of points with the same
        // squark & gluino masses, but varying neutralino masses.
        std::vector< std::vector< std::vector< double >* >* >*
        current_gluino_vector = NULL;
        // this holds all the values for a set of points with the same
        // squark masses, but varying gluino masses.

        while( input_file_reader.read_line( &input_line_as_string ) )
          // the evaluation of the conditional reads in the next line.
          {

            // read in the line for interpretting:
            input_line_as_stream.clear();
            input_line_as_stream.str( input_line_as_string );

            // update what the masses for the last line were:
            last_squark_mass = current_squark_mass;
            last_gluino_mass = current_gluino_mass;

            input_line_as_stream >> current_squark_mass
                                 >> current_gluino_mass;
            input_line_as_stream >> current_value;
            if( grid_size_still_unknown )
              {
                neutralino_masses.push_back( current_value );
              }
            input_line_as_stream >> current_value;
            if( grid_size_still_unknown )
              {
                neutralino_masses.push_back( current_value );
              }

            current_acceptance_vector = new std::vector< double >();
            while( !(input_line_as_stream.eof()) )
              {

                input_line_as_stream >> current_value;
                current_acceptance_vector->push_back( current_value );

              }

            if( current_squark_mass > last_squark_mass )
              // if the squark mass has changed, we need a new vector for the
              // values for varying gluino masses.
              {

                // 1st record the last set of points with the same squark mass:
                if( NULL != current_gluino_vector )
                  // if we have a vector to add...
                  {

                    current_gluino_vector->push_back(
                                                   current_neutralino_vector );
                    values.push_back( current_gluino_vector );

                    // debugging:
                    /**std::cout
                    << std::endl
                    << "pushed back a neutralino vector of size "
                    << current_neutralino_vector->size()
                    << " and a gluino vector of size "
                    << current_gluino_vector->size();**/

                  }

                // prepare a new vector for the new squark mass:
                current_gluino_vector
                = new std::vector< std::vector< std::vector< double >* >* >();
                // prepare a new vector for the new gluino mass:
                current_neutralino_vector
                = new std::vector< std::vector< double >* >();

              }
            else if( current_gluino_mass > last_gluino_mass )
              /* if the gluino mass has changed (but the squark mass has not
               * changed, we need a new vector for the values for varying
               * neutralino masses.
               */
              {

                // 1st record the last set of points with the same gluino mass:
                current_gluino_vector->push_back( current_neutralino_vector );

                // debugging:
                /**std::cout
                << std::endl
                << "pushed back a neutralino vector of size "
                << current_neutralino_vector->size();**/

                // prepare a new vector for the new gluino mass:
                current_neutralino_vector
                = new std::vector< std::vector< double >* >();

              }

            if( !grid_size_still_unknown
                &&
                ( acceptanceColumns > current_acceptance_vector->size() ) )
              {

                acceptanceColumns = current_acceptance_vector->size();

                std::cout
                << std::endl
                << "LHC-FASER warning! read in an acceptance grid line that"
                << " is shorter than the previous one, so the grid is"
                << " \"ragged\"; hence *all* the lines for this grid will be"
                << " truncated to this smaller length.";
                std::cout << std::endl;

              }
            current_neutralino_vector->push_back( current_acceptance_vector );

            if( grid_size_still_unknown
                &&
                ( 0.0 < last_squark_mass )
                &&
                ( 0.0 < last_gluino_mass )
                &&
                ( ( current_squark_mass > last_squark_mass )
                  ||
                  ( current_gluino_mass > last_gluino_mass ) ) )
              // if we still need to work out the grid step size & we have read
              // in enough to do so...
              {

                acceptanceColumns = current_acceptance_vector->size();

                // first we take note of the 1st point, which should have the
                // lowest co-ordinate values:
                lowestSquarkMass = last_squark_mass;
                lowestGluinoMass = last_gluino_mass;

                // then we work out the grid step size:
                if( current_gluino_mass > last_gluino_mass )
                  // if it was the gluino mass that incremented...
                  {

                    scoloredMassStepSize
                    = ( current_gluino_mass - last_gluino_mass );

                  }
                else
                  // otherwise it was the squark mass that incremented...
                  {

                    scoloredMassStepSize
                    = ( current_squark_mass - last_squark_mass );

                  }

                // debugging:
                /**std::cout
                << std::endl
                << "scoloredMassStepSize = " << scoloredMassStepSize;**/

                // now we need to work out the neutralino mass ratios to the
                // lighter colored sparticle:
                double lowest_neutralino_mass = lowestSquarkMass;
                double middle_neutralino_mass = lowestSquarkMass;
                double highest_neutralino_mass = lowestSquarkMass;
                for( std::vector< double >::const_iterator
                     mass_iterator = neutralino_masses.begin();
                     neutralino_masses.end() > mass_iterator;
                     ++mass_iterator )
                  {

                    if( *mass_iterator < lowest_neutralino_mass )
                      {

                        highest_neutralino_mass = middle_neutralino_mass;
                        middle_neutralino_mass = lowest_neutralino_mass;
                        lowest_neutralino_mass = *mass_iterator;

                      }

                    if( ( *mass_iterator < middle_neutralino_mass )
                        &&
                        ( *mass_iterator > lowest_neutralino_mass ) )
                      {

                        highest_neutralino_mass = middle_neutralino_mass;
                        middle_neutralino_mass = *mass_iterator;

                      }

                    if( ( *mass_iterator < highest_neutralino_mass )
                        &&
                        ( *mass_iterator > middle_neutralino_mass ) )
                      {

                        highest_neutralino_mass = *mass_iterator;

                      }

                    double lighter_scolored_mass = last_squark_mass;
                    if( last_gluino_mass < last_squark_mass )
                      {

                        lighter_scolored_mass = last_gluino_mass;

                      }
                    lowNeutralinoMassRatio
                    = ( lowest_neutralino_mass / lighter_scolored_mass );
                    middleNeutralinoMassRatio
                    = ( middle_neutralino_mass / lighter_scolored_mass );
                    highNeutralinoMassRatio
                    = ( highest_neutralino_mass / lighter_scolored_mass );

                  }
                grid_size_still_unknown = false;

              }

          }  // end of while loop going over the lines of the file.

        // we still have to add the last vectors:
        acceptanceColumns = current_acceptance_vector->size();
        current_gluino_vector->push_back( current_neutralino_vector );
        values.push_back( current_gluino_vector );

        // debugging:
        /**std::cout
        << std::endl
        << "pushed back last current_neutralino_vector of size "
        << current_neutralino_vector->size()
        << " & last current_gluino_vector of size "
        << current_gluino_vector->size();
        std::cout << std::endl;**/

        // the last masses should be the maximal values:
        highestSquarkMass = current_squark_mass;
        highestGluinoMass = current_gluino_mass;

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
         gluino_vector_iterator = values.begin();
         values.end() > gluino_vector_iterator;
         ++gluino_vector_iterator )
      {

        for( std::vector< std::vector< std::vector< double >* >* >::iterator
               neutralino_vector_iterator = (*gluino_vector_iterator)->begin();
             (*gluino_vector_iterator)->end() > neutralino_vector_iterator;
             ++neutralino_vector_iterator )
          {

            for( std::vector< std::vector< double >* >::iterator
                   deletion_iterator = (*neutralino_vector_iterator)->begin();
                 (*neutralino_vector_iterator)->end() > deletion_iterator;
                 ++deletion_iterator )
              {

                delete *deletion_iterator;

              }

            delete *neutralino_vector_iterator;

          }

        delete *gluino_vector_iterator;

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
        x_steps
        = ( ( squarkMass - lowestSquarkMass ) / scoloredMassStepSize );
        int lower_left_x = (int)x_steps;

        if( values.size() > ( lower_left_x + 1 ) )
          // if the x co-ordinate is less than its maximal grid value...
          {

            double
            y_steps
            = ( ( gluinoMass - lowestGluinoMass )
                / scoloredMassStepSize );
            int lower_left_y = (int)y_steps;

            if( values.at( lower_left_x )->size() > ( lower_left_y + 1 ) )
              // if the y co-ordinate is less than its maximal grid value...
              {

                return LHC_FASER_global::square_bilinear_interpolation(
                                                    ( x_steps - lower_left_x ),
                                                    ( y_steps - lower_left_y ),
                                    values.at( lower_left_x )->at( lower_left_y
                                                      )->at( neutralinoElement
                                                   )->at( acceptanceElement ),
                            values.at( ( lower_left_x + 1 ) )->at( lower_left_y
                                                      )->at( neutralinoElement
                                                   )->at( acceptanceElement ),
                    values.at( ( lower_left_x + 1 ) )->at( ( lower_left_y + 1 )
                                                      )->at( neutralinoElement
                                                   )->at( acceptanceElement ),
                            values.at( lower_left_x )->at( ( lower_left_y + 1 )
                                                      )->at( neutralinoElement
                                                 )->at( acceptanceElement ) );

              }
            else if( ( values.at( lower_left_x )->size()
                       == ( lower_left_y + 1 ) )
                     &&
                     ( (double)lower_left_y == y_steps ) )
              // otherwise, if it's on the maximal y edge...
              {

                return LHC_FASER_global::unit_linear_interpolation(
                                                    ( x_steps - lower_left_x ),
                                    values.at( lower_left_x )->at( lower_left_y
                                                      )->at( neutralinoElement
                                                   )->at( acceptanceElement ),
                            values.at( ( lower_left_x + 1 ) )->at( lower_left_y
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
        else if( ( values.size() == ( lower_left_x + 1 ) )
                 &&
                 ( (double)lower_left_x == x_steps ) )
          // otherwise, if it's on the maximal x edge...
          {

            double
            y_steps
            = ( ( gluinoMass - lowestGluinoMass )
                / scoloredMassStepSize );
            int lower_left_y = (int)y_steps;

            if( values.at( lower_left_x )->size() > ( lower_left_y + 1 ) )
              // if the y co-ordinate is less than its maximal grid value...
              {

                return LHC_FASER_global::unit_linear_interpolation(
                                                    ( y_steps - lower_left_y ),
                                    values.at( lower_left_x )->at( lower_left_y
                                                      )->at( neutralinoElement
                                                   )->at( acceptanceElement ),
                            values.at( lower_left_x )->at( ( lower_left_y + 1 )
                                                      )->at( neutralinoElement
                                                 )->at( acceptanceElement ) );

              }
            else if( ( values.at( lower_left_x )->size()
                       == ( lower_left_y + 1 ) )
                     &&
                     ( (double)lower_left_y == y_steps ) )
              // otherwise, if it's on the maximal x & y corner...
              {

                return values.at( lower_left_x )->at( lower_left_y
                         )->at( neutralinoElement )->at( acceptanceElement );

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
   * LHC_FASER_global::square_bilinear_interpolation to get an interpolated
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
        double lighter_scolored = squarkMass;
        if( gluinoMass < squarkMass )
          {

            lighter_scolored = gluinoMass;

          }
        double lighter_neutralino = firstNeutralinoMass;
        double heavier_neutralino = secondNeutralinoMass;
        if( secondNeutralinoMass < firstNeutralinoMass )
          {

            lighter_neutralino = secondNeutralinoMass;
            heavier_neutralino = firstNeutralinoMass;

          }
        lighter_neutralino *= ( 1.0 / lighter_scolored );
        heavier_neutralino *= ( 1.0 / lighter_scolored );
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

        if( ( 1.0 > lighter_neutralino )
            &&
            ( 1.0 > heavier_neutralino ) )
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
            if( lowNeutralinoMassRatio > lighter_neutralino )
              {

                if( lowNeutralinoMassRatio > heavier_neutralino )
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
                else if( middleNeutralinoMassRatio > heavier_neutralino )
                  {

                    return LHC_FASER_global::unit_linear_interpolation(
                                                         ( ( heavier_neutralino
                                                  - lowNeutralinoMassRatio )
                                               / ( middleNeutralinoMassRatio
                                               - lowNeutralinoMassRatio ) ),
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
                else if( highNeutralinoMassRatio > heavier_neutralino )
                  {

                    return LHC_FASER_global::unit_linear_interpolation(
                                                         ( ( heavier_neutralino
                                               - middleNeutralinoMassRatio )
                                                 / ( highNeutralinoMassRatio
                                            - middleNeutralinoMassRatio ) ),
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

                    double heavy_edge_value = 0.0;

                    if( heavyNeutralinoEdgeIsLighterScoloredMass )
                      {

                        if( gluinoMass < squarkMass )
                          {

                            heavy_edge_value = gluinoMass;

                          }
                        else
                          {

                            heavy_edge_value = squarkMass;

                          }

                      }
                    else if( heavyNeutralinoAreaIsConstant )
                    {
                      return vectorElementAt( squarkMass,
                                                gluinoMass,
                                                1,
                                                acceptanceElement );
                    }

                    return LHC_FASER_global::unit_linear_interpolation(
                                                         ( ( heavier_neutralino
                                                 - highNeutralinoMassRatio )
                                                                        / ( 1.0
                                              - highNeutralinoMassRatio ) ),
                                                vectorElementAt( squarkMass,
                                                                   gluinoMass,
                                                                   1,
                                                          acceptanceElement ),
                                                            heavy_edge_value );

                    /*lower_left_point = 1;
                    lower_right_point = 1;
                    upper_right_point = 0;
                    upper_left_point = 0;*/

                  }

              }
            else if( middleNeutralinoMassRatio > lighter_neutralino )
              {

                if( middleNeutralinoMassRatio > heavier_neutralino )
                  {

                    double lower_left_value = vectorElementAt( squarkMass,
                                                                 gluinoMass,
                                                                 2,
                                                          acceptanceElement );
                    double upper_right_value = vectorElementAt( squarkMass,
                                                                  gluinoMass,
                                                                  4,
                                                          acceptanceElement );
                    double upper_left_value = vectorElementAt( squarkMass,
                                                                 gluinoMass,
                                                                 1,
                                                          acceptanceElement );
                    double lower_right_value = ( lower_left_value
                                                 + upper_right_value
                                                 - upper_left_value );

                    return LHC_FASER_global::square_bilinear_interpolation(
                                                         ( ( lighter_neutralino
                                                  - lowNeutralinoMassRatio )
                                               / ( middleNeutralinoMassRatio
                                               - lowNeutralinoMassRatio ) ),
                                                         ( ( heavier_neutralino
                                                  - lowNeutralinoMassRatio )
                                               / ( middleNeutralinoMassRatio
                                               - lowNeutralinoMassRatio ) ),
                                                              lower_left_value,
                                                             lower_right_value,
                                                             upper_right_value,
                                                            upper_left_value );

                    /*lower_left_point = 3;
                    lower_right_point = -2;
                    upper_right_point = 5;
                    upper_left_point = 2;*/

                  }
                else if( highNeutralinoMassRatio > heavier_neutralino )
                  {

                    return LHC_FASER_global::square_bilinear_interpolation(
                                                         ( ( lighter_neutralino
                                                  - lowNeutralinoMassRatio )
                                               / ( middleNeutralinoMassRatio
                                               - lowNeutralinoMassRatio ) ),
                                                         ( ( heavier_neutralino
                                               - middleNeutralinoMassRatio )
                                                 / ( highNeutralinoMassRatio
                                            - middleNeutralinoMassRatio ) ),
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

                    double heavy_edge_value = 0.0;

                    if( heavyNeutralinoEdgeIsLighterScoloredMass )
                      {

                        if( gluinoMass < squarkMass )
                          {

                            heavy_edge_value = gluinoMass;

                          }
                        else
                          {

                            heavy_edge_value = squarkMass;

                          }

                      }
                    else if( heavyNeutralinoAreaIsConstant )
                    {
                      heavy_edge_value = vectorElementAt( squarkMass,
                                                            gluinoMass,
                                                            3,
                                                          acceptanceElement );
                    }

                    return LHC_FASER_global::square_bilinear_interpolation(
                                                         ( ( lighter_neutralino
                                                  - lowNeutralinoMassRatio )
                                               / ( middleNeutralinoMassRatio
                                               - lowNeutralinoMassRatio ) ),
                                                         ( ( heavier_neutralino
                                                 - highNeutralinoMassRatio )
                                                                        / ( 1.0
                                              - highNeutralinoMassRatio ) ),
                                                vectorElementAt( squarkMass,
                                                                   gluinoMass,
                                                                   0,
                                                          acceptanceElement ),
                                                vectorElementAt( squarkMass,
                                                                   gluinoMass,
                                                                   3,
                                                          acceptanceElement ),
                                                              heavy_edge_value,
                                                            heavy_edge_value );

                    /*lower_left_point = 1;
                    lower_right_point = 4;
                    upper_right_point = 0;
                    upper_left_point = 0;*/

                  }

              }
            else if( highNeutralinoMassRatio > lighter_neutralino )
              {

                if( highNeutralinoMassRatio > heavier_neutralino )
                  {

                    // debugging:
                    /**std::cout
                    << std::endl
                    << "so, highNeutralinoMassRatio > both > middle"
                    << std::endl;
                    std::cout << std::endl;**/

                    double lower_left_value = vectorElementAt( squarkMass,
                                                                 gluinoMass,
                                                                 4,
                                                         acceptanceElement );

                    // debugging:
                    /**std::cout
                    << std::endl
                    << "lower_left_value = " << lower_left_value;
                    std::cout << std::endl;**/

                    double upper_right_value = vectorElementAt( squarkMass,
                                                                  gluinoMass,
                                                                  5,
                                                          acceptanceElement );

                    // debugging:
                    /**std::cout
                    << std::endl
                    << "upper_right_value = " << upper_right_value;
                    std::cout << std::endl;**/

                    double upper_left_value = vectorElementAt( squarkMass,
                                                                 gluinoMass,
                                                                 3,
                                                          acceptanceElement );

                    // debugging:
                    /**std::cout
                    << std::endl
                    << "upper_left_value = " << upper_left_value;
                    std::cout << std::endl;**/

                    double lower_right_value = ( lower_left_value
                                                 + upper_right_value
                                                 - upper_left_value );

                    // debugging:
                    /**std::cout
                    << std::endl
                    << "lower_right_value = " << lower_right_value;
                    std::cout << std::endl;**/

                    return LHC_FASER_global::square_bilinear_interpolation(
                                                         ( ( lighter_neutralino
                                               - middleNeutralinoMassRatio )
                                                 / ( highNeutralinoMassRatio
                                            - middleNeutralinoMassRatio ) ),
                                                         ( ( heavier_neutralino
                                               - middleNeutralinoMassRatio )
                                                 / ( highNeutralinoMassRatio
                                            - middleNeutralinoMassRatio ) ),
                                                              lower_left_value,
                                                             lower_right_value,
                                                             upper_right_value,
                                                            upper_left_value );

                    /*lower_left_point = 5;
                    lower_right_point = -4;
                    upper_right_point = 6;
                    upper_left_point = 4;*/

                  }
                else
                  {

                    double heavy_edge_value = 0.0;

                    if( heavyNeutralinoEdgeIsLighterScoloredMass )
                      {

                        if( gluinoMass < squarkMass )
                          {

                            heavy_edge_value = gluinoMass;

                          }
                        else
                          {

                            heavy_edge_value = squarkMass;

                          }

                      }
                    else if( heavyNeutralinoAreaIsConstant )
                    {
                      heavy_edge_value = vectorElementAt( squarkMass,
                                                            gluinoMass,
                                                            3,
                                                          acceptanceElement );
                    }

                    return LHC_FASER_global::square_bilinear_interpolation(
                                                         ( ( lighter_neutralino
                                               - middleNeutralinoMassRatio )
                                                 / ( highNeutralinoMassRatio
                                            - middleNeutralinoMassRatio ) ),
                                                         ( ( heavier_neutralino
                                                 - highNeutralinoMassRatio )
                                                                        / ( 1.0
                                              - highNeutralinoMassRatio ) ),
                                                vectorElementAt( squarkMass,
                                                                   gluinoMass,
                                                                   3,
                                                          acceptanceElement ),
                                                vectorElementAt( squarkMass,
                                                                   gluinoMass,
                                                                   5,
                                                          acceptanceElement ),
                                                              heavy_edge_value,
                                                            heavy_edge_value );

                    /*lower_left_point = 4;
                    lower_right_point = 6;
                    upper_right_point = 0;
                    upper_left_point = 0;*/

                  }

              }
            else
              {

                double heavy_edge_value = 0.0;

                if( heavyNeutralinoEdgeIsLighterScoloredMass )
                  {

                    if( gluinoMass < squarkMass )
                      {

                        heavy_edge_value = gluinoMass;

                      }
                    else
                      {

                        heavy_edge_value = squarkMass;

                      }

                  }
                else if( heavyNeutralinoAreaIsConstant )
                {
                  heavy_edge_value = vectorElementAt( squarkMass,
                                                        gluinoMass,
                                                        5,
                                                      acceptanceElement );
                }

                return LHC_FASER_global::unit_linear_interpolation(
                                                         ( ( lighter_neutralino
                                                 - highNeutralinoMassRatio )
                                                                    / ( 1.0
                                              - highNeutralinoMassRatio ) ),
                                                vectorElementAt( squarkMass,
                                                                   gluinoMass,
                                                                   5,
                                                          acceptanceElement ),
                                                            heavy_edge_value );

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
   * LHC_FASER_global::square_bilinear_interpolation to get an interpolated
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
        double lighter_scolored = squarkMass;
        if( gluinoMass < squarkMass )
          {

            lighter_scolored = gluinoMass;

          }
        double
        neutralino_ratio = ( degenerateNeutralinoMass / lighter_scolored );

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

        if( 1.0 > neutralino_ratio )
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
            if( lowNeutralinoMassRatio > neutralino_ratio )
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
            else if( middleNeutralinoMassRatio > neutralino_ratio )
              {

                double lower_left_value = vectorElementAt( squarkMass,
                                                             gluinoMass,
                                                             2,
                                                          acceptanceElement );
                double upper_right_value = vectorElementAt( squarkMass,
                                                              gluinoMass,
                                                              4,
                                                          acceptanceElement );

                return LHC_FASER_global::unit_linear_interpolation(
                                                         ( ( neutralino_ratio
                                                  - lowNeutralinoMassRatio )
                                               / ( middleNeutralinoMassRatio
                                               - lowNeutralinoMassRatio ) ),
                                                              lower_left_value,
                                                           upper_right_value );

                /*lower_left_point = 3;
                lower_right_point = -2;
                upper_right_point = 5;
                upper_left_point = 2;*/

              }
            else if( highNeutralinoMassRatio > neutralino_ratio )
              {

                double lower_left_value = vectorElementAt( squarkMass,
                                                             gluinoMass,
                                                             4,
                                                         acceptanceElement );

                double upper_right_value = vectorElementAt( squarkMass,
                                                              gluinoMass,
                                                              5,
                                                          acceptanceElement );

                return LHC_FASER_global::unit_linear_interpolation(
                                                           ( ( neutralino_ratio
                                               - middleNeutralinoMassRatio )
                                                 / ( highNeutralinoMassRatio
                                            - middleNeutralinoMassRatio ) ),
                                                              lower_left_value,
                                                           upper_right_value );

                /*lower_left_point = 5;
                lower_right_point = -4;
                upper_right_point = 6;
                upper_left_point = 4;*/

              }
            else
              {

                double heavy_edge_value = 0.0;

                if( heavyNeutralinoEdgeIsLighterScoloredMass )
                  {

                    if( gluinoMass < squarkMass )
                      {

                        heavy_edge_value = gluinoMass;

                      }
                    else
                      {

                        heavy_edge_value = squarkMass;

                      }

                  }
                else if( heavyNeutralinoAreaIsConstant )
                {
                  return vectorElementAt( squarkMass,
                                            gluinoMass,
                                            5,
                                            acceptanceElement );
                }

                return LHC_FASER_global::unit_linear_interpolation(
                                                         ( ( lighter_neutralino
                                                 - highNeutralinoMassRatio )
                                                                    / ( 1.0
                                              - highNeutralinoMassRatio ) ),
                                                vectorElementAt( squarkMass,
                                                                   gluinoMass,
                                                                   5,
                                                          acceptanceElement ),
                                                            heavy_edge_value );

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



  leptonAcceptanceGrid::leptonAcceptanceGrid( std::string const* const gridFileLocation,
                                input_handler const* const shortcut ) :
    shortcut( shortcut )
  {
    std::string gridFileName( *gridFileLocation );
    gridFileName.append( "/gluino+gluino_acceptance.dat");
    gluinoGluinoTable = new acceptanceGrid( &gridFileName );
    gridFileName.assign( *gridFileLocation );
    gridFileName.append( "/squark+gluino_acceptance.dat");
    squarkGluinoTable = new acceptanceGrid( &gridFileName );
    gridFileName.assign( *gridFileLocation );
    gridFileName.append( "/squark+antisquark_acceptance.dat");
    squarkAntisquarkTable = new acceptanceGrid( &gridFileName );
    gridFileName.assign( *gridFileLocation );
    gridFileName.append( "/squark+squark_acceptance.dat");
    squarkSquarkTable = new acceptanceGrid( &gridFileName );
  }

  leptonAcceptanceGrid::~leptonAcceptanceGrid()
  {
    delete gluinoGluinoTable;
    delete squarkGluinoTable;
    delete squarkAntisquarkTable;
    delete squarkSquarkTable;
  }


  double
  jetAcceptanceGrid::getValue(
                          signed_particle_shortcut_pair const* const scoloreds,
                             colored_cascade const* const first_sQCD_cascade,
                             colored_cascade const* const second_sQCD_cascade,
                             int const requested_column,
                     bool const heavy_neutralino_edge_is_lighter_scolored_mass,
                             bool const heavy_neutralino_area_is_constant )
  const
  /* this interpolates the requested column based on the squark, gluino, &
   * electroweakino masses. it fudges some cases that were not properly done in
   * the single-quark-flavor approximation.
   */
  {

    // debugging:
    /**std::cout
    << std::endl
    << "leptonAcceptanceGrid::getValue(...) called.";
    std::cout << std::endl;**/

    double appropriate_gluino_mass
    = shortcut->get_gluino()->get_absolute_mass();
    double appropriate_squark_mass;

    // gluino+gluino combinations:
    if( first_sQCD_cascade->is_gx()
        &&
        second_sQCD_cascade->is_gx() )
      {

        appropriate_squark_mass = shortcut->get_average_squarks4_mass();
        if( appropriate_gluino_mass > appropriate_squark_mass )
          {

            appropriate_squark_mass
            = ( appropriate_gluino_mass + ( 0.5 * appropriate_squark_mass ) );

          }

        return gluinoGluinoTable->valueAt( appropriate_squark_mass,
                                              appropriate_gluino_mass,
                          first_sQCD_cascade->get_EWino()->get_absolute_mass(),
                         second_sQCD_cascade->get_EWino()->get_absolute_mass(),
                                              requested_column,
                                heavy_neutralino_edge_is_lighter_scolored_mass,
                                           heavy_neutralino_area_is_constant );

      }
    else if( first_sQCD_cascade->is_gsx()
             &&
             second_sQCD_cascade->is_gx() )
      {

        appropriate_squark_mass
        = first_sQCD_cascade->get_product()->get_absolute_mass();

        return gluinoGluinoTable->valueAt( appropriate_squark_mass,
                                              appropriate_gluino_mass,
                          first_sQCD_cascade->get_EWino()->get_absolute_mass(),
                         second_sQCD_cascade->get_EWino()->get_absolute_mass(),
                                              requested_column,
                                heavy_neutralino_edge_is_lighter_scolored_mass,
                                           heavy_neutralino_area_is_constant );

      }
    else if( first_sQCD_cascade->is_gx()
             &&
             second_sQCD_cascade->is_gsx() )
      {

        appropriate_squark_mass
        = second_sQCD_cascade->get_product()->get_absolute_mass();

        return gluinoGluinoTable->valueAt( appropriate_squark_mass,
                                              appropriate_gluino_mass,
                          first_sQCD_cascade->get_EWino()->get_absolute_mass(),
                         second_sQCD_cascade->get_EWino()->get_absolute_mass(),
                                              requested_column,
                                heavy_neutralino_edge_is_lighter_scolored_mass,
                                           heavy_neutralino_area_is_constant );

      }
    else if( first_sQCD_cascade->is_gsx()
             &&
             second_sQCD_cascade->is_gsx() )
      {

        appropriate_squark_mass
        = ( 0.5 * ( first_sQCD_cascade->get_product()->get_absolute_mass()
                 + second_sQCD_cascade->get_product()->get_absolute_mass() ) );

        // debugging:
        /**std::cout
        << std::endl
        << "is gsx-gsx; appropriate_squark_mass = " << appropriate_squark_mass;
        std::cout << std::endl;**/

        return gluinoGluinoTable->valueAt( appropriate_squark_mass,
                                              appropriate_gluino_mass,
                          first_sQCD_cascade->get_EWino()->get_absolute_mass(),
                         second_sQCD_cascade->get_EWino()->get_absolute_mass(),
                                              requested_column,
                                heavy_neutralino_edge_is_lighter_scolored_mass,
                                           heavy_neutralino_area_is_constant );

      }
    // squark+gluino combinations:
    else if( first_sQCD_cascade->is_sx()
             &&
             second_sQCD_cascade->is_gx() )
      {

        appropriate_squark_mass
        = first_sQCD_cascade->get_decayer()->get_absolute_mass();

        return squarkGluinoTable->valueAt( appropriate_squark_mass,
                                              appropriate_gluino_mass,
                          first_sQCD_cascade->get_EWino()->get_absolute_mass(),
                         second_sQCD_cascade->get_EWino()->get_absolute_mass(),
                                              requested_column,
                                heavy_neutralino_edge_is_lighter_scolored_mass,
                                           heavy_neutralino_area_is_constant );

      }
    else if( first_sQCD_cascade->is_gx()
             &&
             second_sQCD_cascade->is_sx() )
      {

        appropriate_squark_mass
        = second_sQCD_cascade->get_decayer()->get_absolute_mass();

        return squarkGluinoTable->valueAt( appropriate_squark_mass,
                                              appropriate_gluino_mass,
                          first_sQCD_cascade->get_EWino()->get_absolute_mass(),
                         second_sQCD_cascade->get_EWino()->get_absolute_mass(),
                                              requested_column,
                                heavy_neutralino_edge_is_lighter_scolored_mass,
                                           heavy_neutralino_area_is_constant );

      }
    else if( first_sQCD_cascade->is_sx()
             &&
             second_sQCD_cascade->is_gsx() )
      {

        appropriate_squark_mass
        = ( 0.5 * ( first_sQCD_cascade->get_decayer()->get_absolute_mass()
                 + second_sQCD_cascade->get_product()->get_absolute_mass() ) );
        if( appropriate_squark_mass > appropriate_gluino_mass )
          {

            appropriate_gluino_mass
            = ( appropriate_squark_mass + ( 0.5 * appropriate_gluino_mass ) );

          }

        return squarkGluinoTable->valueAt( appropriate_squark_mass,
                                              appropriate_gluino_mass,
                          first_sQCD_cascade->get_EWino()->get_absolute_mass(),
                         second_sQCD_cascade->get_EWino()->get_absolute_mass(),
                                              requested_column,
                                heavy_neutralino_edge_is_lighter_scolored_mass,
                                           heavy_neutralino_area_is_constant );

      }
    else if( first_sQCD_cascade->is_gsx()
             &&
             second_sQCD_cascade->is_sx() )
      {

        appropriate_squark_mass
        = ( 0.5 * ( first_sQCD_cascade->get_product()->get_absolute_mass()
                 + second_sQCD_cascade->get_decayer()->get_absolute_mass() ) );
        if( appropriate_squark_mass > appropriate_gluino_mass )
          {

            appropriate_gluino_mass
            = ( appropriate_squark_mass + ( 0.5 * appropriate_gluino_mass ) );

          }

        return squarkGluinoTable->valueAt( appropriate_squark_mass,
                                              appropriate_gluino_mass,
                          first_sQCD_cascade->get_EWino()->get_absolute_mass(),
                         second_sQCD_cascade->get_EWino()->get_absolute_mass(),
                                              requested_column,
                                heavy_neutralino_edge_is_lighter_scolored_mass,
                                           heavy_neutralino_area_is_constant );

      }
    else if( first_sQCD_cascade->is_sgx()
             &&
             second_sQCD_cascade->is_gx() )
      {

        appropriate_squark_mass
        = first_sQCD_cascade->get_decayer()->get_absolute_mass();

        return squarkGluinoTable->valueAt( appropriate_squark_mass,
                                              appropriate_gluino_mass,
                          first_sQCD_cascade->get_EWino()->get_absolute_mass(),
                         second_sQCD_cascade->get_EWino()->get_absolute_mass(),
                                              requested_column,
                                heavy_neutralino_edge_is_lighter_scolored_mass,
                                           heavy_neutralino_area_is_constant );

      }
    else if( first_sQCD_cascade->is_gx()
             &&
             second_sQCD_cascade->is_sgx() )
      {

        appropriate_squark_mass
        = second_sQCD_cascade->get_decayer()->get_absolute_mass();

        return squarkGluinoTable->valueAt( appropriate_squark_mass,
                                              appropriate_gluino_mass,
                          first_sQCD_cascade->get_EWino()->get_absolute_mass(),
                         second_sQCD_cascade->get_EWino()->get_absolute_mass(),
                                              requested_column,
                                heavy_neutralino_edge_is_lighter_scolored_mass,
                                           heavy_neutralino_area_is_constant );

      }
    else if( first_sQCD_cascade->is_sgx()
             &&
             second_sQCD_cascade->is_gsx() )
      {

        appropriate_squark_mass
        = second_sQCD_cascade->get_product()->get_absolute_mass();

        return squarkGluinoTable->valueAt( appropriate_squark_mass,
                                              appropriate_gluino_mass,
                          first_sQCD_cascade->get_EWino()->get_absolute_mass(),
                         second_sQCD_cascade->get_EWino()->get_absolute_mass(),
                                              requested_column,
                                heavy_neutralino_edge_is_lighter_scolored_mass,
                                           heavy_neutralino_area_is_constant );

      }
    else if( first_sQCD_cascade->is_gsx()
             &&
             second_sQCD_cascade->is_sgx() )
      {

        appropriate_squark_mass
        = first_sQCD_cascade->get_product()->get_absolute_mass();

        return squarkGluinoTable->valueAt( appropriate_squark_mass,
                                              appropriate_gluino_mass,
                          first_sQCD_cascade->get_EWino()->get_absolute_mass(),
                         second_sQCD_cascade->get_EWino()->get_absolute_mass(),
                                              requested_column,
                                heavy_neutralino_edge_is_lighter_scolored_mass,
                                           heavy_neutralino_area_is_constant );

      }
    else if( first_sQCD_cascade->is_sgsx()
             &&
             second_sQCD_cascade->is_gx() )
      {

        appropriate_squark_mass
        = first_sQCD_cascade->get_product()->get_absolute_mass();

        return squarkGluinoTable->valueAt( appropriate_squark_mass,
                                              appropriate_gluino_mass,
                          first_sQCD_cascade->get_EWino()->get_absolute_mass(),
                         second_sQCD_cascade->get_EWino()->get_absolute_mass(),
                                              requested_column,
                                heavy_neutralino_edge_is_lighter_scolored_mass,
                                           heavy_neutralino_area_is_constant );

      }
    else if( first_sQCD_cascade->is_gx()
             &&
             second_sQCD_cascade->is_sgsx() )
      {

        appropriate_squark_mass
        = second_sQCD_cascade->get_product()->get_absolute_mass();

        return squarkGluinoTable->valueAt( appropriate_squark_mass,
                                              appropriate_gluino_mass,
                          first_sQCD_cascade->get_EWino()->get_absolute_mass(),
                         second_sQCD_cascade->get_EWino()->get_absolute_mass(),
                                              requested_column,
                                heavy_neutralino_edge_is_lighter_scolored_mass,
                                           heavy_neutralino_area_is_constant );

      }
    else if( first_sQCD_cascade->is_sgsx()
             &&
             second_sQCD_cascade->is_gsx() )
      {

        appropriate_squark_mass
        = ( 0.5 * ( first_sQCD_cascade->get_product()->get_absolute_mass()
                 + second_sQCD_cascade->get_product()->get_absolute_mass() ) );

        return gluinoGluinoTable->valueAt( appropriate_squark_mass,
                                              appropriate_gluino_mass,
                          first_sQCD_cascade->get_EWino()->get_absolute_mass(),
                         second_sQCD_cascade->get_EWino()->get_absolute_mass(),
                                              requested_column,
                                heavy_neutralino_edge_is_lighter_scolored_mass,
                                           heavy_neutralino_area_is_constant );

      }
    else if( first_sQCD_cascade->is_gsx()
             &&
             second_sQCD_cascade->is_sgsx() )
      {

        appropriate_squark_mass
        = ( 0.5 * ( first_sQCD_cascade->get_product()->get_absolute_mass()
                 + second_sQCD_cascade->get_product()->get_absolute_mass() ) );

        return gluinoGluinoTable->valueAt( appropriate_squark_mass,
                                              appropriate_gluino_mass,
                          first_sQCD_cascade->get_EWino()->get_absolute_mass(),
                         second_sQCD_cascade->get_EWino()->get_absolute_mass(),
                                              requested_column,
                                heavy_neutralino_edge_is_lighter_scolored_mass,
                                           heavy_neutralino_area_is_constant );

      }
    // squark+squark combinations:
    else if( first_sQCD_cascade->is_sx()
             &&
             second_sQCD_cascade->is_sx() )
      {

        appropriate_squark_mass
        = ( 0.5 * ( first_sQCD_cascade->get_decayer()->get_absolute_mass()
                 + second_sQCD_cascade->get_decayer()->get_absolute_mass() ) );
        if( appropriate_squark_mass > appropriate_gluino_mass )
          {

            appropriate_gluino_mass
            = ( appropriate_squark_mass + ( 0.5 * appropriate_gluino_mass ) );

          }

        if( scoloreds->first_is_not_antiparticle()
            &&
            scoloreds->second_is_not_antiparticle() )
          {

            return squarkSquarkTable->valueAt( appropriate_squark_mass,
                                                  appropriate_gluino_mass,
                          first_sQCD_cascade->get_EWino()->get_absolute_mass(),
                         second_sQCD_cascade->get_EWino()->get_absolute_mass(),
                                                  requested_column,
                                heavy_neutralino_edge_is_lighter_scolored_mass,
                                           heavy_neutralino_area_is_constant );

          }
        else if( ( scoloreds->first_is_not_antiparticle()
                   &&
                   !(scoloreds->second_is_not_antiparticle()) )
                 ||
                 ( !(scoloreds->first_is_not_antiparticle())
                   &&
                   scoloreds->second_is_not_antiparticle() ) )
          {

            return squarkAntisquarkTable->valueAt( appropriate_squark_mass,
                                                      appropriate_gluino_mass,
                          first_sQCD_cascade->get_EWino()->get_absolute_mass(),
                         second_sQCD_cascade->get_EWino()->get_absolute_mass(),
                                                      requested_column,
                                heavy_neutralino_edge_is_lighter_scolored_mass,
                                           heavy_neutralino_area_is_constant );

          }
        else
          {

            return 0.0;

          }

      }
    else if( first_sQCD_cascade->is_sgx()
             &&
             second_sQCD_cascade->is_sx() )
      {

        appropriate_squark_mass
        = second_sQCD_cascade->get_decayer()->get_absolute_mass();

        return squarkGluinoTable->valueAt( appropriate_squark_mass,
                                              appropriate_gluino_mass,
                          first_sQCD_cascade->get_EWino()->get_absolute_mass(),
                         second_sQCD_cascade->get_EWino()->get_absolute_mass(),
                                              requested_column,
                                heavy_neutralino_edge_is_lighter_scolored_mass,
                                           heavy_neutralino_area_is_constant );

      }
    else if( first_sQCD_cascade->is_sgsx()
             &&
             second_sQCD_cascade->is_sx() )
      {

        appropriate_squark_mass
        = ( 0.5 * ( first_sQCD_cascade->get_product()->get_absolute_mass()
                 + second_sQCD_cascade->get_decayer()->get_absolute_mass() ) );
        if( appropriate_squark_mass > appropriate_gluino_mass )
          {

            appropriate_gluino_mass
            = ( appropriate_squark_mass + ( 0.5 * appropriate_gluino_mass ) );

          }

        if( scoloreds->second_is_not_antiparticle() )
          {

            return squarkSquarkTable->valueAt( appropriate_squark_mass,
                                                  appropriate_gluino_mass,
                          first_sQCD_cascade->get_EWino()->get_absolute_mass(),
                         second_sQCD_cascade->get_EWino()->get_absolute_mass(),
                                                  requested_column,
                                heavy_neutralino_edge_is_lighter_scolored_mass,
                                           heavy_neutralino_area_is_constant );

          }
        else
          {

            return squarkAntisquarkTable->valueAt( appropriate_squark_mass,
                                                      appropriate_gluino_mass,
                          first_sQCD_cascade->get_EWino()->get_absolute_mass(),
                         second_sQCD_cascade->get_EWino()->get_absolute_mass(),
                                                      requested_column,
                                heavy_neutralino_edge_is_lighter_scolored_mass,
                                           heavy_neutralino_area_is_constant );


          }

      }
    else if( first_sQCD_cascade->is_sx()
             &&
             second_sQCD_cascade->is_sgx() )
      {

        appropriate_squark_mass
        = first_sQCD_cascade->get_decayer()->get_absolute_mass();

        return squarkGluinoTable->valueAt( appropriate_squark_mass,
                                              appropriate_gluino_mass,
                          first_sQCD_cascade->get_EWino()->get_absolute_mass(),
                         second_sQCD_cascade->get_EWino()->get_absolute_mass(),
                                              requested_column,
                                heavy_neutralino_edge_is_lighter_scolored_mass,
                                           heavy_neutralino_area_is_constant );

      }
    else if( first_sQCD_cascade->is_sgx()
             &&
             second_sQCD_cascade->is_sgx() )
      {

        appropriate_squark_mass
        = ( 0.5 * ( first_sQCD_cascade->get_decayer()->get_absolute_mass()
                 + second_sQCD_cascade->get_decayer()->get_absolute_mass() ) );

        if( scoloreds->first_is_not_antiparticle()
            &&
            scoloreds->second_is_not_antiparticle() )
          {

            return squarkSquarkTable->valueAt( appropriate_squark_mass,
                                                  appropriate_gluino_mass,
                          first_sQCD_cascade->get_EWino()->get_absolute_mass(),
                         second_sQCD_cascade->get_EWino()->get_absolute_mass(),
                                                  requested_column,
                                heavy_neutralino_edge_is_lighter_scolored_mass,
                                           heavy_neutralino_area_is_constant );

          }
        else if( ( scoloreds->first_is_not_antiparticle()
                   &&
                   !(scoloreds->second_is_not_antiparticle()) )
                 ||
                 ( !(scoloreds->first_is_not_antiparticle())
                   &&
                   scoloreds->second_is_not_antiparticle() ) )
          {

            return squarkAntisquarkTable->valueAt( appropriate_squark_mass,
                                                      appropriate_gluino_mass,
                          first_sQCD_cascade->get_EWino()->get_absolute_mass(),
                         second_sQCD_cascade->get_EWino()->get_absolute_mass(),
                                                      requested_column,
                                heavy_neutralino_edge_is_lighter_scolored_mass,
                                           heavy_neutralino_area_is_constant );

          }
        else
          {

            return 0.0;

          }

      }
    else if( first_sQCD_cascade->is_sgsx()
             &&
             second_sQCD_cascade->is_sgx() )
      {

        appropriate_squark_mass
        = first_sQCD_cascade->get_product()->get_absolute_mass();

        return squarkGluinoTable->valueAt( appropriate_squark_mass,
                                              appropriate_gluino_mass,
                          first_sQCD_cascade->get_EWino()->get_absolute_mass(),
                         second_sQCD_cascade->get_EWino()->get_absolute_mass(),
                                              requested_column,
                                heavy_neutralino_edge_is_lighter_scolored_mass,
                                           heavy_neutralino_area_is_constant );

      }
    else if( first_sQCD_cascade->is_sx()
             &&
             second_sQCD_cascade->is_sgsx() )
      {

        appropriate_squark_mass
        = ( 0.5 * ( first_sQCD_cascade->get_decayer()->get_absolute_mass()
                 + second_sQCD_cascade->get_product()->get_absolute_mass() ) );

        if( scoloreds->first_is_not_antiparticle()
            &&
            scoloreds->second_is_not_antiparticle() )
          {

            return squarkSquarkTable->valueAt( appropriate_squark_mass,
                                                  appropriate_gluino_mass,
                          first_sQCD_cascade->get_EWino()->get_absolute_mass(),
                         second_sQCD_cascade->get_EWino()->get_absolute_mass(),
                                                  requested_column,
                                heavy_neutralino_edge_is_lighter_scolored_mass,
                                           heavy_neutralino_area_is_constant );

          }
        else if( ( scoloreds->first_is_not_antiparticle()
                   &&
                   !(scoloreds->second_is_not_antiparticle()) )
                 ||
                 ( !(scoloreds->first_is_not_antiparticle())
                   &&
                   scoloreds->second_is_not_antiparticle() ) )
          {

            return squarkAntisquarkTable->valueAt( appropriate_squark_mass,
                                                      appropriate_gluino_mass,
                          first_sQCD_cascade->get_EWino()->get_absolute_mass(),
                         second_sQCD_cascade->get_EWino()->get_absolute_mass(),
                                                      requested_column,
                                heavy_neutralino_edge_is_lighter_scolored_mass,
                                           heavy_neutralino_area_is_constant );

          }
        else
          {

            return 0.0;

          }

      }
    else if( first_sQCD_cascade->is_sgx()
             &&
             second_sQCD_cascade->is_sgsx() )
      {

        appropriate_squark_mass
        = second_sQCD_cascade->get_product()->get_absolute_mass();

        return squarkGluinoTable->valueAt( appropriate_squark_mass,
                                              appropriate_gluino_mass,
                          first_sQCD_cascade->get_EWino()->get_absolute_mass(),
                         second_sQCD_cascade->get_EWino()->get_absolute_mass(),
                                              requested_column,
                                heavy_neutralino_edge_is_lighter_scolored_mass,
                                           heavy_neutralino_area_is_constant );

      }
    else if( first_sQCD_cascade->is_sgsx()
             &&
             second_sQCD_cascade->is_sgsx() )
      {

        appropriate_squark_mass
        = ( 0.5 * ( first_sQCD_cascade->get_product()->get_absolute_mass()
                 + second_sQCD_cascade->get_product()->get_absolute_mass() ) );

        return gluinoGluinoTable->valueAt( appropriate_squark_mass,
                                              appropriate_gluino_mass,
                          first_sQCD_cascade->get_EWino()->get_absolute_mass(),
                         second_sQCD_cascade->get_EWino()->get_absolute_mass(),
                                              requested_column,
                                heavy_neutralino_edge_is_lighter_scolored_mass,
                                           heavy_neutralino_area_is_constant );

      }
    else
      {

        return 0.0;

      }

  }



  jet_acceptance_value::jet_acceptance_value(
                                             leptonAcceptanceGrid const* given_table,
                                             int const given_acceptance_column,
                         signed_particle_shortcut_pair const* const given_pair,
                         colored_cascade const* const given_first_sQCD_cascade,
                     colored_cascade const* const given_second_sQCD_cascade ) :
    jet_table( given_table ),
    acceptance_column( given_acceptance_column ),
    scoloredPair( given_pair ),
    first_sQCD_cascade( given_first_sQCD_cascade ),
    second_sQCD_cascade( given_second_sQCD_cascade ),
    stored_acceptance( 0.0 ),
    not_already_calculated( true )
  {

    // just an initialization list.

  }

  jet_acceptance_value::~jet_acceptance_value()
  {

    // does nothing.

  }



  lepton_acceptance_value::lepton_acceptance_value(
                                             leptonAcceptanceGrid const* given_table,
                         signed_particle_shortcut_pair const* const given_pair,
                         colored_cascade const* const given_first_sQCD_cascade,
                     colored_cascade const* const given_second_sQCD_cascade ) :
    lepton_table( given_table ),
    scoloredPair( given_pair ),
    first_sQCD_cascade( given_first_sQCD_cascade ),
    second_sQCD_cascade( given_second_sQCD_cascade ),
    stored_acceptance(),
    not_already_calculated( true )
  {

    // just an initialization list.

  }

  lepton_acceptance_value::~lepton_acceptance_value()
  {

    // does nothing.

  }


  leptonAcceptanceParameterSet const*
  lepton_acceptance_value::get_acceptance()
  {

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: lepton_acceptance_value::get_acceptance() called."
    << std::endl << "not_already_calculated = " << not_already_calculated;
    std::cout << std::endl;
    int bin_counter = -2;**/

    if( not_already_calculated )
      {

        std::vector< double > acceptance_passer;
        int acceptance_counter = 0;
        double current_acceptance = lepton_table->getValue( scoloredPair,
                                                            first_sQCD_cascade,
                                                           second_sQCD_cascade,
                                                          acceptance_counter,
                                                             true );
        // the 1st element of the vector is the effective squark mass, so
        // should be extrapolated to the lighter scolored mass.

        // debugging:
        /**std::cout
        << std::endl
        << "starting current_acceptance = " << current_acceptance;
        std::cout << std::endl;**/

        while( 0.0 <= current_acceptance )
          {

            // debugging:
            /**std::cout
            << std::endl
            << "debugging: ( 0.0 <= current_acceptance ) is true:"
            << " current_acceptance = " << current_acceptance
            << " for bin energy " << ( 2.0 * ++bin_counter);
            std::cout << std::endl;**/

            acceptance_passer.push_back( current_acceptance );
            ++acceptance_counter;
            current_acceptance = lepton_table->getValue( scoloredPair,
                                                          first_sQCD_cascade,
                                                          second_sQCD_cascade,
                                                          acceptance_counter,
                                                          false );
            // the remaining elements of the vector are acceptances, so should
            // be extrapolated to zero.

            // debugging:
            /**std::cout
            << std::endl
            << "now current_acceptance = " << current_acceptance;
            std::cout << std::endl;**/

          }

        // debugging:
        /**std::cout
        << std::endl
        << "debugging: lepton_acceptance_value::get_acceptance() about to call"
        << " stored_acceptance.set_values( &acceptance_passer )";
        std::cout << std::endl;**/

        stored_acceptance.setValues( &acceptance_passer );

        not_already_calculated = false;

      }

    return &stored_acceptance;

  }



  jet_acceptance_table::jet_acceptance_table(
                                      std::string const* const given_directory,
                                           std::string const* const given_name,
                                 input_handler const* const given_shortcuts ) :
    readied_for_new_point( given_shortcuts->get_readier() ),
    table_name( *given_name ),
    acceptance_grids( given_directory,
                      given_shortcuts )
  {

    // just an initialization list.

  }

  jet_acceptance_table::~jet_acceptance_table()
  {

    for( std::vector< jet_acceptance_value* >::iterator
         deletion_iterator = acceptances.begin();
         acceptances.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

  }


  jet_acceptance_value*
  jet_acceptance_table::get_acceptance( int const given_acceptance_column,
                         signed_particle_shortcut_pair const* const given_pair,
                         colored_cascade const* const given_first_sQCD_cascade,
                       colored_cascade const* const given_second_sQCD_cascade )
  {

    if( needs_to_prepare_for_this_point() )
      {

        for( std::vector< jet_acceptance_value* >::iterator
             deletion_iterator = acceptances.begin();
             acceptances.end() > deletion_iterator;
             ++deletion_iterator )
          {

            delete *deletion_iterator;

          }

        acceptances.clear();
        finish_preparing_for_this_point();

      }

    jet_acceptance_value* return_pointer = NULL;

    for( std::vector< jet_acceptance_value* >::iterator
         value_iterator = acceptances.begin();
         acceptances.end() > value_iterator;
         ++value_iterator )
      {

        if( (*value_iterator)->is_requested( given_acceptance_column,
                                             given_pair,
                                             given_first_sQCD_cascade,
                                             given_second_sQCD_cascade ) )
          {

            return_pointer = *value_iterator;
            value_iterator = acceptances.end();

          }

      }

    if( NULL == return_pointer )
      {

        return_pointer = new jet_acceptance_value( &acceptance_grids,
                                                   given_acceptance_column,
                                                   given_pair,
                                                   given_first_sQCD_cascade,
                                                   given_second_sQCD_cascade );

        acceptances.push_back( return_pointer );

      }


    return return_pointer;

  }


  lepton_acceptance_table::lepton_acceptance_table(
                                      std::string const* const given_directory,
                                             int const given_LHC_energy_in_TeV,
                                 input_handler const* const given_shortcuts ) :
    readied_for_new_point( given_shortcuts->get_readier() ),
    LHC_energy_in_TeV( given_LHC_energy_in_TeV ),
    acceptance_grids( given_directory,
                      given_shortcuts )
  {

    // just an initialization list.

  }

  lepton_acceptance_table::~lepton_acceptance_table()
  {

    for( std::vector< lepton_acceptance_value* >::iterator
         deletion_iterator = acceptances.begin();
         acceptances.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

  }


  lepton_acceptance_value*
  lepton_acceptance_table::get_acceptance(
                         signed_particle_shortcut_pair const* const given_pair,
                         colored_cascade const* const given_first_sQCD_cascade,
                       colored_cascade const* const given_second_sQCD_cascade )
  {

    if( needs_to_prepare_for_this_point() )
      {

        for( std::vector< lepton_acceptance_value* >::iterator
             deletion_iterator = acceptances.begin();
             acceptances.end() > deletion_iterator;
             ++deletion_iterator )
          {

            delete *deletion_iterator;

          }

        acceptances.clear();
        finish_preparing_for_this_point();

      }

    lepton_acceptance_value* return_pointer = NULL;

    for( std::vector< lepton_acceptance_value* >::iterator
         value_iterator = acceptances.begin();
         acceptances.end() > value_iterator;
         ++value_iterator )
      {

        if( (*value_iterator)->is_requested( given_pair,
                                             given_first_sQCD_cascade,
                                             given_second_sQCD_cascade ) )
          {

            return_pointer = *value_iterator;
            value_iterator = acceptances.end();

          }

      }

    if( NULL == return_pointer )
      {

        return_pointer = new lepton_acceptance_value( &acceptance_grids,
                                                      given_pair,
                                                      given_first_sQCD_cascade,
                                                   given_second_sQCD_cascade );

        acceptances.push_back( return_pointer );

      }

    return return_pointer;

  }



  kinematics_table::kinematics_table(
                                  jet_acceptance_table* const given_jets_table,
                                      int const given_jet_acceptance_column,
                         lepton_acceptance_table* const given_leptons_table ) :
    jets_table( given_jets_table ),
    jet_acceptance_column( given_jet_acceptance_column ),
    leptons_table( given_leptons_table )
  {

    // just an initialization list.

  }

  kinematics_table::~kinematics_table()
  {

    // does nothing.

  }



  kinematics_table_set::kinematics_table_set(
                                             int const given_LHC_energy_in_TeV,
                                 std::string const* const given_grid_directory,
                               std::string const* const given_jet_subdirectory,
                             lepton_acceptance_table* const given_lepton_table,
                                 input_handler const* const given_shortcuts ) :
    LHC_energy_in_TeV( given_LHC_energy_in_TeV ),
    grid_directory( *given_grid_directory ),
    jet_grid_subdirectory( *given_grid_directory ),
    lepton_table( given_lepton_table ),
    shortcut( given_shortcuts )
  {

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: kinematics_table_set::kinematics_table_set( "
    << given_LHC_energy_in_TeV << ", " << *given_grid_directory
    << ", " << *given_jet_subdirectory << ", " << given_lepton_table << ", "
    << given_shortcuts << ") about to make a new jet_acceptance_table( "
    << energy_name << ", " << *given_grid_directory << ", "
    << given_shortcuts << " )";
    std::cout << std::endl;**/
    jet_grid_subdirectory.append( "/jets/" );
    jet_grid_subdirectory.append( *given_jet_subdirectory );
    jet_table = new jet_acceptance_table( &jet_grid_subdirectory,
                                          given_grid_directory,
                                          given_shortcuts );
    jet_grid_subdirectory.assign( *given_jet_subdirectory );
    // we set jet_grid_subdirectory back to *given_jet_subdirectory for later
    // comparisons.

  }

  kinematics_table_set::~kinematics_table_set()
  {

    for( std::vector< kinematics_table* >::iterator
         deletion_iterator = kinematics_tables.begin();
         kinematics_tables.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

    delete jet_table;

  }


  kinematics_table*
  kinematics_table_set::get_table( int const given_acceptance_column )
  // this returns the kinematics_table for the requested setup.
  {

    kinematics_table* return_pointer = NULL;
    // this starts as NULL so that we know if it wasn't found among the
    // existing instances.

    for( std::vector< kinematics_table* >::iterator
         table_iterator = kinematics_tables.begin();
         kinematics_tables.end() > table_iterator;
         ++table_iterator )
      {

        if( given_acceptance_column == (*table_iterator)->get_column() )
          // if we find the requested column...
          {

            return_pointer = *table_iterator;
            // note the table to return.
            table_iterator = kinematics_tables.end();
            // stop looking.

          }

      }

    if( NULL == return_pointer )
      // if we didn't find an existing match, we make a new instance, store
      // it & return a pointer to it:
      {

        return_pointer = new kinematics_table( jet_table,
                                               given_acceptance_column,
                                               lepton_table );

        kinematics_tables.push_back( return_pointer );

      }

    return return_pointer;

  }



  kinematics_handler::kinematics_handler(
                                 input_handler const* const given_shortcuts ) :
    shortcut( given_shortcuts )
  {

    // just an initialization list

  }

  kinematics_handler::~kinematics_handler()
  {

    for( std::vector< kinematics_table_set* >::iterator
         deletion_iterator = table_sets.begin();
         table_sets.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

    for( std::vector< lepton_acceptance_table* >::iterator
         deletion_iterator = lepton_tables.begin();
         lepton_tables.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

  }


  lepton_acceptance_table*
  kinematics_handler::get_lepton_acceptance_table(
                                                  int const LHC_energy_in_TeV )
  {

    // we have to check to see if the lepton_acceptance_table already exists:
    lepton_acceptance_table* return_pointer = NULL;
    for( std::vector< lepton_acceptance_table* >::iterator
         lepton_iterator = lepton_tables.begin();
         lepton_tables.end() > lepton_iterator;
         ++lepton_iterator )
      {

        if( LHC_energy_in_TeV == (*lepton_iterator)->get_energy() )
          {

            return_pointer = *lepton_iterator;
            lepton_iterator = lepton_tables.end();

          }

      }
    if( NULL == return_pointer )
      {

        std::string
        lepton_grids( *(shortcut->inspect_path_to_kinematics_grids()) );
        std::stringstream energy_stream( "" );
        energy_stream << "/" << LHC_energy_in_TeV << "TeV";
        lepton_grids.append( energy_stream.str() );
        lepton_grids.append( "/leptons/" );

        return_pointer = new lepton_acceptance_table( &lepton_grids,
                                                      LHC_energy_in_TeV,
                                                      shortcut );

        lepton_tables.push_back( return_pointer );

      }

    return return_pointer;

  }


  kinematics_table_set*
  kinematics_handler::get_table_set( int const LHC_energy_in_TeV,
                              std::string const* const given_jet_subdirectory )
  {

    kinematics_table_set* return_pointer = NULL;
    // this starts as NULL so that we know if it wasn't found among the
    // existing instances.

    for( std::vector< kinematics_table_set* >::iterator
         set_iterator = table_sets.begin();
         table_sets.end() > set_iterator;
         ++set_iterator )
      {

        if( (*set_iterator)->is_requested( LHC_energy_in_TeV,
                                           given_jet_subdirectory ) )
          {

            return_pointer = *set_iterator;
            set_iterator = table_sets.end();
            // stop looking.

          }

      }

    if( NULL == return_pointer )
      // if we didn't find an existing match, we make a new instance, store
      // it & return a pointer to it:
      {

        std::string
        energy_name( *(shortcut->inspect_path_to_kinematics_grids()) );
        std::stringstream energy_stream( "" );
        energy_stream << "/" << LHC_energy_in_TeV << "TeV";
        energy_name.append( energy_stream.str() );

        // we have to check to see if the lepton_acceptance_table already
        // exists:
        lepton_acceptance_table* lepton_pointer = NULL;
        for( std::vector< lepton_acceptance_table* >::iterator
             lepton_iterator = lepton_tables.begin();
             lepton_tables.end() > lepton_iterator;
             ++lepton_iterator )
          {

            if( LHC_energy_in_TeV == (*lepton_iterator)->get_energy() )
              {

                lepton_pointer = *lepton_iterator;
                lepton_iterator = lepton_tables.end();

              }

          }
        if( NULL == lepton_pointer )
          {

            std::string lepton_grids( energy_name );
            lepton_grids.append( "/leptons/" );

            lepton_pointer = new lepton_acceptance_table( &lepton_grids,
                                                          LHC_energy_in_TeV,
                                                          shortcut );

            lepton_tables.push_back( lepton_pointer );

          }

        return_pointer = new kinematics_table_set( LHC_energy_in_TeV,
                                                   &energy_name,
                                                   given_jet_subdirectory,
                                                   lepton_pointer,
                                                   shortcut );

        table_sets.push_back( return_pointer );

      }

    return return_pointer;

  }

}  // end of LHC_FASER namespace.

