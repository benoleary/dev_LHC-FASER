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
                                        inputHandler const* const shortcut ) :
    shortcut( shortcut )
  {
    std::string gridFileName( *gridFileLocation );
    gridFileName.append( "/gluino+gluino_acceptance.dat" );
    gluinoTable = new acceptanceGrid( &gridFileName );
    gridFileName.assign( *gridFileLocation );
    gridFileName.append( "/squark+antisquark_acceptance.dat" );
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
                                        int const requestedColumn )
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

    double appropriateGluinoMass = shortcut->getGluinoMass();
    double appropriateSquarkMass = squark->get_absolute_mass();
    if( appropriateSquarkMass > appropriateGluinoMass )
    {
      appropriateGluinoMass = ( appropriateSquarkMass + 1.0 );
    }
    return squarkTable->neutralinoIndependentValueAt( appropriateSquarkMass,
                                                      appropriateGluinoMass,
                                                      requestedColumn );
    // lepton acceptance parameters should be independent of the neutralino
    // masses used.
  }

  double
  leptonAcceptanceGrid::getGluinoValue(
                             CppSLHA::particle_property_set const* const ewino,
                                        int const requestedColumn,
                                      bool const lookingForEffectiveSquarkMass,
                                  bool const lookingForLeptonAcceptanceNumber )
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

    double appropriateGluinoMass = shortcut->getGluinoMass();
    double appropriateSquarkMass = shortcut->getAverageSquarks4Mass();
    if( appropriateGluinoMass > appropriateSquarkMass )
    {
      appropriateSquarkMass = ( appropriateGluinoMass + 1.0 );
    }
    return gluinoTable->valueAt( appropriateSquarkMass,
                                 appropriateGluinoMass,
                                 ewino->get_absolute_mass(),
                                 ewino->get_absolute_mass(),
                                 requestedColumn,
                                 lookingForEffectiveSquarkMass,
                                 lookingForLeptonAcceptanceNumber );
  }



  squarkBasedLeptonAcceptanceTable::squarkBasedLeptonAcceptanceTable(
                                  leptonAcceptanceGrid const* const lookupGrid,
                            CppSLHA::particle_property_set const* const squark,
                                        inputHandler const* const shortcut ) :
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
                                        inputHandler const* const shortcut ) :
    leptonAcceptanceTable( lookupGrid,
                           shortcut->getGluino(),
                           shortcut )
  {
    // just an initialization list.
  }

  gluinoBasedLeptonAcceptanceTable::~gluinoBasedLeptonAcceptanceTable()
  {
    // does nothing.
  }




  leptonAcceptanceParameterSet::leptonAcceptanceParameterSet(
                                                 inputHandler* const shortcut,
                                  leptonAcceptanceTable const* acceptanceTable,
                          CppSLHA::particle_property_set const* const scolored,
                          CppSLHA::particle_property_set const* const ewino ) :
    getsReadiedForNewPoint( shortcut->getReadier() ),
    effectiveSquarkMassHolder(),
    shortcut( shortcut ),
    acceptanceTable( acceptanceTable ),
    scolored( scolored ),
    ewino( ewino ),
    binSize( defaultBinSize ),
    transverseMomentumCut( defaultTransverseMomentumCut ),
    effectiveSquarkMass( CppSLHA::CppSLHA_global::really_wrong_value ),
    pseudorapidityAcceptance( CppSLHA::CppSLHA_global::really_wrong_value )
  {
    acceptanceBins.push_back( CppSLHA::CppSLHA_global::really_wrong_value );
    // this is to ensure that there is at least 1 entry for the transverse
    // momentum cut.
  }

  leptonAcceptanceParameterSet::leptonAcceptanceParameterSet(
                                                 inputHandler* const shortcut,
                                  leptonAcceptanceTable const* acceptanceTable,
                          CppSLHA::particle_property_set const* const scolored,
                             CppSLHA::particle_property_set const* const ewino,
                                                          double const binSize,
                                         double const transverseMomentumCut ) :
    getsReadiedForNewPoint( shortcut->getReadier() ),
    effectiveSquarkMassHolder(),
    shortcut( shortcut ),
    acceptanceTable( acceptanceTable ),
    scolored( scolored ),
    ewino( ewino ),
    binSize( binSize ),
    transverseMomentumCut( transverseMomentumCut ),
    effectiveSquarkMass( CppSLHA::CppSLHA_global::really_wrong_value ),
    pseudorapidityAcceptance( CppSLHA::CppSLHA_global::really_wrong_value )
  {
    acceptanceBins.push_back( CppSLHA::CppSLHA_global::really_wrong_value );
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




  leptonAcceptancesForOneScolored::leptonAcceptancesForOneScolored(
                                           inputHandler const* const shortcut,
                          CppSLHA::particle_property_set const* const scolored,
                                    leptonAcceptanceGrid const* acceptanceGrid,
                                                          double const binSize,
                                         double const transverseMomentumCut ) :
    shortcut( shortcut ),
    scolored( scolored ),
    binSize( binSize ),
    transverseMomentumCut( transverseMomentumCut )
  {
    if( CppSLHA::PDG_code::gluino == scolored->get_PDG_code() )
    {
      acceptanceTable = new gluinoBasedLeptonAcceptanceTable( acceptanceGrid,
                                                              shortcut );
    }
    else
    {
      acceptanceTable = new squarkBasedLeptonAcceptanceTable( acceptanceGrid,
                                                              scolored,
                                                              shortcut );
    }
  }

  leptonAcceptancesForOneScolored::~leptonAcceptancesForOneScolored()
  {
    for( std::vector< leptonAcceptanceParameterSet* >::iterator
         deletionIterator = parameterSets->begin();
         parameterSets->end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
  }

  leptonAcceptanceParameterSet*
  leptonAcceptancesForOneScolored::getParameterSet(
                            CppSLHA::particle_property_set const* const ewino )
  // this returns the leptonAcceptanceParameterSet for the requested
  // electroweakino.
  {
    leptonAcceptanceParameterSet* returnPointer( NULL );
    // we look to see if we already have a leptonAcceptanceParameterSet for
    // this electroweakino:
    for( std::vector< leptonAcceptanceParameterSet* >::iterator
         searchIterator = parameterSets->begin();
         parameterSets->end() > searchIterator;
         ++searchIterator )
    {
      if( ewino == (*searchIterator)->getEwino() )
      {
        returnPointer = *searchIterator;
        searchIterator = parameterSets->end();
      }
    }
    if( NULL == returnPointer )
      // if we do not already have a leptonAcceptanceParameterSet for this
      // electroweakino, we make a new instance:
    {
      returnPointer = new leptonAcceptanceParameterSet( shortcut,
                                                        acceptanceTable,
                                                        scolored,
                                                        ewino,
                                                        binSize,
                                                       transverseMomentumCut );
      parameterSets.push_back( returnPointer );
    }
    return returnPointer;
  }


  leptonAcceptancesForOneBeamEnergy::leptonAcceptancesForOneBeamEnergy(
                                           inputHandler const* const shortcut,
                                                          int const beamEnergy,
                                  std::string const* const gridFileSetLocation,
                                                          double const binSize,
                                         double const transverseMomentumCut ) :
    shortcut( shortcut ),
    beamEnergy( beamEnergy ),
    binSize( binSize ),
    transverseMomentumCut( transverseMomentumCut )
  {
    std::string gridFileLocation( *gridFileSetLocation );
    std::stringstream energyStream( "" );
    energyStream << "/" << beamEnergy << "TeV/leptons/";
    gridFileLocation.append( energyStream.str() );
    acceptanceGrid = new leptonAcceptanceGrid( &gridFileLocation,
                                               shortcut );
  }

  leptonAcceptancesForOneBeamEnergy::~leptonAcceptancesForOneBeamEnergy()
  {
    for( std::vector< leptonAcceptancesForOneScolored* >::iterator
         deletionIterator = acceptanceSets->begin();
         acceptanceSets->end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
    delete acceptanceGrid;
  }

  leptonAcceptancesForOneScolored*
  leptonAcceptancesForOneBeamEnergy::getParameterSets(
                         CppSLHA::particle_property_set const* const scolored )
  // this returns the leptonAcceptancesForOneScolored for the requested
  // colored sparticle.
  {
    leptonAcceptancesForOneScolored* returnPointer( NULL );
    // we look to see if we already have a leptonAcceptancesForOneScolored for
    // this colored sparticle:
    for( std::vector< leptonAcceptancesForOneScolored* >::iterator
         searchIterator = acceptanceSets->begin();
         acceptanceSets->end() > searchIterator;
         ++searchIterator )
    {
      if( scolored == (*searchIterator)->getScolored() )
      {
        returnPointer = *searchIterator;
        searchIterator = acceptanceSets->end();
      }
    }
    if( NULL == returnPointer )
      // if we do not already have a leptonAcceptanceParameterSet for this
      // colored sparticle, we make a new instance:
    {
      returnPointer = new leptonAcceptancesForOneScolored( shortcut,
                                                           scolored,
                                                           acceptanceGrid,
                                                           binSize,
                                                       transverseMomentumCut );
      acceptanceSets.push_back( returnPointer );
    }
    return returnPointer;
  }



  leptonAcceptanceHandler::leptonAcceptanceHandler(
                                           inputHandler const* const shortcut,
                               std::string const* const gridFileSetLocation ) :
    shortcut( shortcut ),
    gridFileSetLocation( *gridFileSetLocation )
  {
    // just an initialization list.
  }

  leptonAcceptanceHandler::~leptonAcceptanceHandler()
  {
    for( std::vector< leptonAcceptancesForOneBeamEnergy* >::iterator
         deletionIterator = acceptanceTables->begin();
         acceptanceTables->end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
  }

  leptonAcceptancesForOneBeamEnergy*
  leptonAcceptanceHandler::getLeptonAcceptancesForOneBeamEnergy(
                                                          int const beamEnergy,
                                                          double const binSize,
                                           double const transverseMomentumCut )
  /* this looks to see if there is an existing
   * leptonAcceptancesForOneBeamEnergy with the requested values, & if not,
   * makes 1, & returns the pointer.
   */
  {
    leptonAcceptancesForOneBeamEnergy* returnPointer( NULL );
    // we look to see if we already have a leptonAcceptancesForOneBeamEnergy
    // for these values:
    for( std::vector< leptonAcceptancesForOneBeamEnergy* >::iterator
         searchIterator = acceptanceTables->begin();
         acceptanceTables->end() > searchIterator;
         ++searchIterator )
    {
      if( scolored == (*searchIterator)->isRequested( beamEnergy,
                                                      binSize,
                                                      transverseMomentumCut ) )
      {
        returnPointer = *searchIterator;
        searchIterator = acceptanceTables->end();
      }
    }
    if( NULL == returnPointer )
      // if we do not already have a leptonAcceptanceParameterSet for this
      // electroweakino, we make a new instance:
    {
      returnPointer = new leptonAcceptancesForOneBeamEnergy( shortcut,
                                                             beamEnergy,
                                                          &gridFileSetLocation,
                                                             binSize,
                                                       transverseMomentumCut );
      acceptanceTables.push_back( returnPointer );
    }
    return returnPointer;
  }
}  // end of LHC_FASER namespace.

