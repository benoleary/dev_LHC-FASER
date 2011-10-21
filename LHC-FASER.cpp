/*
 * LHC-FASER.cpp
 *
 *  Created on: 10 Nov 2010
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

#include "LHC-FASER.hpp"

namespace LHC_FASER
{
  /* constructor with specified CppSLHA & user-defined path to grids,
   * whether returned event rates should be in nb, pb or fb & whether to use
   * LO or NLO:
   */
  lhcFaser::lhcFaser( CppSLHA::CppSLHA0* const spectrumData,
                      std::string const pathToGrids,
                      std::string const crossSectionUnit,
                      bool const usingNlo ) :
    usingOwnCppSlha( false )
  {
    initialize( spectrumData,
                pathToGrids,
                crossSectionUnit,
                usingNlo );
  }

  /* constructor with specified SLHA file & user-defined path to grids,
   * whether returned event rates should be in nb, pb or fb & whether to use
   * LO or NLO:
   */
  lhcFaser::lhcFaser( std::string const slhaFileName,
                      std::string const pathToGrids,
                      std::string const crossSectionUnit,
                      bool const usingNlo ) :
    spectrumData( new CppSLHA::CppSLHA2( slhaFileName ) ),
    usingOwnCppSlha( true )
  {
    initialize( spectrumData,
                pathToGrids,
                crossSectionUnit,
                usingNlo );
  }

  lhcFaser::~lhcFaser()
  {
    if( usingOwnCppSlha )
    {
      delete spectrumData;
    }
    for( std::vector< signalHandler* >::iterator
         deletionIterator( signalSet.begin() );
         signalSet.end() > deletionIterator;
         deletionIterator++ )
    {
      delete *deletionIterator;
    }
    delete signalPreparationDefinitions;
    delete shortcut;
    delete inputSource;
    delete crossSectionSource;
    delete jetPlusMetAcceptanceSource;
    delete electroweakCascadeSource;
    delete readier;
  }


  void
  lhcFaser::initialize( CppSLHA::CppSLHA0* const spectrumData,
                        std::string const pathToGrids,
                        std::string const crossSectionUnit,
                        bool const usingNlo )
  // this is used by the constructors to do most of the construction.
  {
    this->spectrumData = spectrumData;
    this->pathToGrids.assign( pathToGrids );
    this->usingNlo = usingNlo;
    if( 0 == crossSectionUnit.compare( "pb" ) )
    {
      crossSectionUnitFactor = 1.0;
    }
    else if( 0 == crossSectionUnit.compare( "fb" ) )
    {
      crossSectionUnitFactor = 1000.0;
    }
    else if( 0 == crossSectionUnit.compare( "nb" ) )
    {
      crossSectionUnitFactor = 0.001;
    }
    else
    {
      std::cout
      << "LHC-FASER::error! given_unit has to be pb/fb/nb (all lowercase"
      << " characters). Carrying on, assuming fb.";
      std::cout << std::endl;
      crossSectionUnitFactor = 1000.0;
    }

    readier = new readierForNewPoint();
    inputSource = new inputHandler( spectrumData,
                                    spectrumData->get_particle_spectrum(),
                                    pathToGrids,
                                    readier );
    crossSectionSource = new crossSectionHandler( inputSource );
    jetPlusMetAcceptanceSource = new jetPlusMetAcceptanceHandler( inputSource,
                                                        &(this->pathToGrids) );
    electroweakCascadeSource = new electroweakCascadeHandler( inputSource,
                                                        &(this->pathToGrids) );

    shortcut = new signalShortcuts( inputSource,
                                    crossSectionSource,
                                    jetPlusMetAcceptanceSource,
                                    fullCascadeSetSource );

    signalPreparationDefinitions = new signalDefinitionSet( shortcut );

  }

}  // end of LHC_FASER namespace.


