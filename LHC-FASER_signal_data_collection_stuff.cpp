/*
 * LHC-FASER_signal_data_collection_stuff.cpp
 *
 *  Created on: 21 Mar 2011
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

#include "LHC-FASER_signal_data_collection_stuff.hpp"

namespace LHC_FASER
{
  signalShortcuts::signalShortcuts( inputHandler const* const inputShortcuts,
                                    crossSectionHandler* const crossSections,
                      jetPlusMetAcceptanceHandler* const jetPlusMetAcceptances,
                            //leptonAcceptanceHandler* const leptonAcceptances,
                                   fullCascadeSetFactory* const cascadeSets ) :
    inputShortcuts( inputShortcuts ),
    crossSections( given_cross_section_pointer ),
    jetPlusMetAcceptances( jetPlusMetAcceptances ),
    cascadeSets( cascadeSets )
  {
    // just an initialization list.
  }

  signalShortcuts::~signalShortcuts()
  {
    // does nothing.
  }


  double
  signalShortcuts::getUncertainty()
  /* this looks at the various mass differences in the hierarchy & then
   * guesses how uncertain the signal value is (as a multiplicative factor
   * > 1.0).
   */
  {
    // actually, right now it returns a blanket 1.3 factor...
    return 1.3;
  }



  signalDefinitionSet::signalDefinitionSet( signalShortcuts* const shortcut ) :
    acceptanceCutSet(),
    shortcut( shortcut ),
    crossSections( NULL ),
    leptonKinematics( NULL ),
    jetPlusMetAcceptance( NULL )
  {
    // just an initialization list.
  }

  signalDefinitionSet::signalDefinitionSet(
                                     signalDefinitionSet* const copySource  ) :
    acceptanceCutSet( (acceptanceCutSet)copySource ),
    shortcut( copySource->getShortcuts() ),
    crossSections( copySource->getCrossSections() ),
    leptonKinematics( copySource->getLeptonKinematics() ),
    jetPlusMetAcceptance( copySource->getJetPlusMetAcceptance() )
  {
    // just an initialization list.
  }

  signalDefinitionSet::~signalDefinitionSet()
  {
    // does nothing.
  }



  productionChannelPointerSet::productionChannelPointerSet(
                                  signalDefinitionSet* const signalDefinitions,
                       signedParticleShortcutPair const* const scoloredPair ) :
    crossSection( signalDefinitions->getCrossSections()->getTable(
                                                              scoloredPair ) ),
    firstCascadeSet( signalDefinitions->getShortcuts()->getCascadeSets(
                                                          )->getFullCascadeSet(
                                              scoloredPair->getFirstParticle(),
                                        signalDefinitions->getBeamEnergy() ) ),
    secondCascadeSet( signalDefinitions->getShortcuts()->getCascadeSets(
                                                          )->getFullCascadeSet(
                                             scoloredPair->getSecondParticle(),
                                         signalDefinitions->getBeamEnergy() ) )
  {
    // just an initialization list.
  }

  productionChannelPointerSet::~productionChannelPointerSet()
  {
    // does nothing.
  }
}  // end of LHC_FASER namespace
