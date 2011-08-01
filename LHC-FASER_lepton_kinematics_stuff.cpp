/*
 * LHC-FASER_lepton_kinematics_stuff.cpp
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

#include "LHC-FASER_lepton_kinematics_stuff.hpp"

namespace LHC_FASER
{
  singleLeptonCut::singleLeptonCut() :
    acceptanceCutSet()
  {
    // just an initialization list.
  }

  singleLeptonCut::~singleLeptonCut()
  {
    // does nothing.
  }



  leptonAcceptanceGrid::leptonAcceptanceGrid(
                                     std::string const* const gridFileLocation,
                                        input_handler const* const shortcut ) :
    shortcut( shortcut )
  {
    std::string gridFileName( *gridFileLocation );
    gridFileName.append( "/gluino+gluino_acceptance.dat");
    gluinoTable = new acceptanceGrid( &gridFileName );
    gridFileName.assign( *gridFileLocation );
    gridFileName.append( "/squark+antisquark_acceptance.dat");
    squarkTable = new acceptanceGrid( &gridFileName );
  }

  leptonAcceptanceGrid::~leptonAcceptanceGrid()
  {
    delete gluinoTable;
    delete squarkTable;
  }


  double
  leptonAcceptanceGrid::getSquarkValue(
                            CppSLHA::particle_property_set const* const squark,
                             CppSLHA::particle_property_set const* const ewino,
                                        int const requestedColumn,
                                        bool const lookingForPseudorapidity )
  const
  /* this interpolates the requested column based on the squark, gluino, &
   * electroweakino masses. it fudges some cases that were not properly done in
   * the single-quark-flavor approximation.
   */
  {
    // debugging:
    /**std::cout
    << std::endl
    << "leptonAcceptanceGrid::getSquarkValue(...) called.";
    std::cout << std::endl;**/

    double appropriateGluinoMass = shortcut->get_gluino()->get_absolute_mass();
    double appropriateSquarkMass = squark->get_absolute_mass();
    if( appropriateSquarkMass > appropriateGluinoMass )
    {
      appropriateGluinoMass = appropriateSquarkMass;
    }
    return squarkTable->valueAt( appropriateSquarkMass,
                                 appropriateGluinoMass,
                                 ewino->get_absolute_mass(),
                                 ewino->get_absolute_mass(),
                                 requestedColumn,
                                 false,
                                 lookingForPseudorapidity );
  }

  double
  leptonAcceptanceGrid::getGluinoValue(
                             CppSLHA::particle_property_set const* const ewino,
                                        int const requestedColumn,
                                      bool const lookingForEffectiveSquarkMass,
                                        bool const lookingForPseudorapidity )
  const
  /* this interpolates the requested column based on the squark, gluino, &
   * electroweakino masses. it fudges some cases that were not properly done in
   * the single-quark-flavor approximation.
   */
  {
    // debugging:
    /**std::cout
    << std::endl
    << "leptonAcceptanceGrid::getGluinoValue(...) called.";
    std::cout << std::endl;**/

    double appropriateGluinoMass = shortcut->get_gluino()->get_absolute_mass();
    double appropriateSquarkMass = squark->get_absolute_mass();
    if( appropriateGluinoMass > appropriateSquarkMass )
    {
      appropriateSquarkMass = appropriateGluinoMass;
    }
    return gluinoTable->valueAt( appropriateSquarkMass,
                                 appropriateGluinoMass,
                                 ewino->get_absolute_mass(),
                                 ewino->get_absolute_mass(),
                                 requestedColumn,
                                 lookingForEffectiveSquarkMass,
                                 lookingForPseudorapidity );
  }



  squarkBasedLeptonAcceptanceTable::squarkBasedLeptonAcceptanceTable(
                                  leptonAcceptanceGrid const* const lookupGrid,
                            CppSLHA::particle_property_set const* const squark,
                                        input_handler const* const shortcut ) :
    leptonAcceptanceTable( lookupGrid,
                           squark,
                           shortcut )
  {
    // just an initialization list.
  }

  squarkBasedLeptonAcceptanceTable::~squarkBasedLeptonAcceptanceTable()
  {
    // does nothing.
  }



  gluinoBasedLeptonAcceptanceTable::gluinoBasedLeptonAcceptanceTable(
                                  leptonAcceptanceGrid const* const lookupGrid,
                                        input_handler const* const shortcut ) :
    leptonAcceptanceTable( lookupGrid,
                           shortcut->get_gluino(),
                           shortcut )
  {
    // just an initialization list.
  }

  gluinoBasedLeptonAcceptanceTable::~gluinoBasedLeptonAcceptanceTable()
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
    effectiveSquarkMass = acceptanceTable->getValue( ewino,
                                                     0 );
    // the 1st element of the vector is the effective squark mass.
    pseudorapidityAcceptance = acceptanceTable->getValue( ewino,
                                                          1 );
    // the 2nd element of the vector is the pseudorapidity cut acceptance.
    acceptanceCounter = 1;
    currentAcceptance = 0.0;
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
      currentAcceptance = acceptanceTable->getValue( ewino,
                                                     acceptanceCounter );
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
  leptonAcceptanceParameterSet::calculateAcceptanceAt(
                                                      double const givenEnergy,
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

        returnValue = ( acceptanceBins.at( lowerBin )
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

