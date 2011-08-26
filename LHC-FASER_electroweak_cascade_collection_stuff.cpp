/*
 * LHC-FASER_electroweak_cascade_collection_stuff.cpp
 *
 *  Created on: 26 Aug 2011
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

#include "LHC-FASER_electroweak_cascade_collection_stuff.hpp"

namespace LHC_FASER
{
  electroweakCascadeSet::electroweakCascadeSet(
                                leptonAcceptanceParameterSet* const kinematics,
                                          particlePointer const coloredDecayer,
                                      particlePointer const electroweakDecayer,
                                         inputHandler const* const shortcut ) :
    kinematics( kinematics ),
    coloredDecayer( coloredDecayer ),
    electroweakDecayer( electroweakDecayer ),
    lighterScolored( NULL ),
    shortcut( shortcut )
  // this is the version for electroweakinos.
  {
    effectiveSquarkMassHolder* effectiveSquarkMassPointer;
    if( shortcut->getGluino() == coloredDecayer )
    {
      effectiveSquarkMassPointer = kinematics;
    }
    else
    {
      effectiveSquarkMassPointer
      = shortcut->getOnShellEffectiveSquarkMass( coloredDecayer );
    }

    if( CppSLHA::PDG_code::neutralino_one
        == electroweakDecayer->get_PDG_code() )
    {
      currentCascade = new lightestNeutralinoCascade();
      cascades.push_back( currentCascade );
    }
    else if( shortcut->isIn( electroweakDecayer->get_PDG_code(),
                             shortcut->getUnstableNeutralinos() ) )
    {
      currentCascade = new neutralinoToSemuCascade( kinematics,
                                                    effectiveSquarkMassPointer,
                                                    coloredDecayer,
                                                    electroweakDecayer,
                                                    shortcut->getSelectronL(),
                                                    shortcut );
      cascades.push_back( currentCascade );
      currentCascade = new neutralinoToSemuCascade( kinematics,
                                                    effectiveSquarkMassPointer,
                                                    coloredDecayer,
                                                    electroweakDecayer,
                                                    shortcut->getSelectronR(),
                                                    shortcut );
      cascades.push_back( currentCascade );
      currentCascade = new neutralinoToSemuCascade( kinematics,
                                                    effectiveSquarkMassPointer,
                                                    coloredDecayer,
                                                    electroweakDecayer,
                                                    shortcut->getSmuonL(),
                                                    shortcut );
      cascades.push_back( currentCascade );
      currentCascade = new neutralinoToSemuCascade( kinematics,
                                                    effectiveSquarkMassPointer,
                                                    coloredDecayer,
                                                    electroweakDecayer,
                                                    shortcut->getSmuonR(),
                                                    shortcut );
      cascades.push_back( currentCascade );
      currentCascade = new neutralinoToStauCascade( kinematics,
                                                    effectiveSquarkMassPointer,
                                                    coloredDecayer,
                                                    electroweakDecayer,
                                                    shortcut->getStauOne(),
                                                    shortcut );
      cascades.push_back( currentCascade );
      currentCascade = new neutralinoToStauCascade( kinematics,
                                                    effectiveSquarkMassPointer,
                                                    coloredDecayer,
                                                    electroweakDecayer,
                                                    shortcut->getStauTwo(),
                                                    shortcut );
      cascades.push_back( currentCascade );
      currentCascade = new neutralinoToZCascade( kinematics,
                                                 effectiveSquarkMassPointer,
                                                 coloredDecayer,
                                                 electroweakDecayer,
                                                 shortcut );
      cascades.push_back( currentCascade );
      for( std::vector< particlePointer >::const_iterator
           bosonIterator(
                  shortcut->getNeutralEwsbScalarsAndPseudoscalars()->begin() );
           shortcut->getNeutralEwsbScalarsAndPseudoscalars()->end()
           > bosonIterator;
           ++bosonIterator )
      {
        currentCascade = new neutralinoToHiggsCascade( kinematics,
                                                    effectiveSquarkMassPointer,
                                                       coloredDecayer,
                                                       electroweakDecayer,
                                                       *bosonIterator,
                                                       shortcut );
        cascades.push_back( currentCascade );
      }
      currentCascade = new neutralinoVirtualCascade( kinematics,
                                                    effectiveSquarkMassPointer,
                                                     coloredDecayer,
                                                     electroweakDecayer,
                                                     shortcut );
      cascades.push_back( currentCascade );
    }
    else if( shortcut->isIn( electroweakDecayer->get_PDG_code(),
                             shortcut->getCharginos() ) )
    {
      currentCascade = new charginoToSemuOrEmuSnuCascade( kinematics,
                                                    effectiveSquarkMassPointer,
                                                          coloredDecayer,
                                                          electroweakDecayer,
                                                     shortcut->getSelectronL(),
                                                          shortcut,
                                                          false );
      cascades.push_back( currentCascade );
      currentCascade = new charginoToSemuOrEmuSnuCascade( kinematics,
                                                    effectiveSquarkMassPointer,
                                                          coloredDecayer,
                                                          electroweakDecayer,
                                                     shortcut->getSelectronR(),
                                                          shortcut,
                                                          false );
      cascades.push_back( currentCascade );
      currentCascade = new charginoToSemuOrEmuSnuCascade( kinematics,
                                                    effectiveSquarkMassPointer,
                                                          coloredDecayer,
                                                          electroweakDecayer,
                                             shortcut->getElectronSneutrinoL(),
                                                          shortcut,
                                                          true );
      cascades.push_back( currentCascade );
      currentCascade = new charginoToSemuOrEmuSnuCascade( kinematics,
                                                    effectiveSquarkMassPointer,
                                                          coloredDecayer,
                                                          electroweakDecayer,
                                             shortcut->getElectronSneutrinoR(),
                                                          shortcut,
                                                          true );
      cascades.push_back( currentCascade );
      currentCascade = new charginoToSemuOrEmuSnuCascade( kinematics,
                                                    effectiveSquarkMassPointer,
                                                          coloredDecayer,
                                                          electroweakDecayer,
                                                         shortcut->getSmuonL(),
                                                          shortcut,
                                                          false );
      cascades.push_back( currentCascade );
      currentCascade = new charginoToSemuOrEmuSnuCascade( kinematics,
                                                    effectiveSquarkMassPointer,
                                                          coloredDecayer,
                                                          electroweakDecayer,
                                                         shortcut->getSmuonR(),
                                                          shortcut,
                                                          false );
      cascades.push_back( currentCascade );
      currentCascade = new charginoToSemuOrEmuSnuCascade( kinematics,
                                                    effectiveSquarkMassPointer,
                                                          coloredDecayer,
                                                          electroweakDecayer,
                                                 shortcut->getMuonSneutrinoL(),
                                                          shortcut,
                                                          true );
      cascades.push_back( currentCascade );
      currentCascade = new charginoToSemuOrEmuSnuCascade( kinematics,
                                                    effectiveSquarkMassPointer,
                                                          coloredDecayer,
                                                          electroweakDecayer,
                                                 shortcut->getMuonSneutrinoR(),
                                                          shortcut,
                                                          true );
      cascades.push_back( currentCascade );
      currentCascade = new charginoToStauOrTauSnuCascade( kinematics,
                                          effectiveSquarkMassPointer,
                                                          coloredDecayer,
                                                          electroweakDecayer,
                                                        shortcut->getStauOne(),
                                                          shortcut,
                                                          false );
      cascades.push_back( currentCascade );
      currentCascade = new charginoToStauOrTauSnuCascade( kinematics,
                                          effectiveSquarkMassPointer,
                                                          coloredDecayer,
                                                          electroweakDecayer,
                                                        shortcut->getStauTwo(),
                                                          shortcut,
                                                          false );
      cascades.push_back( currentCascade );
      currentCascade = new charginoToStauOrTauSnuCascade( kinematics,
                                          effectiveSquarkMassPointer,
                                                          coloredDecayer,
                                                          electroweakDecayer,
                                                  shortcut->getTauSneutrinoL(),
                                                          shortcut,
                                                          true );
      cascades.push_back( currentCascade );
      currentCascade = new charginoToStauOrTauSnuCascade( kinematics,
                                          effectiveSquarkMassPointer,
                                                          coloredDecayer,
                                                          electroweakDecayer,
                                                  shortcut->getTauSneutrinoR(),
                                                          shortcut,
                                                          true );
      cascades.push_back( currentCascade );
      currentCascade = new charginoToWCascade( kinematics,
                                               effectiveSquarkMassPointer,
                                               coloredDecayer,
                                               electroweakDecayer,
                                               shortcut );
      cascades.push_back( currentCascade );
      currentCascade = new charginoToHiggsCascade( kinematics,
                                                   effectiveSquarkMassPointer,
                                                   coloredDecayer,
                                                   electroweakDecayer,
                                              shortcut->getChargedEwsbScalar(),
                                                   shortcut );
      cascades.push_back( currentCascade );
      currentCascade = new charginoVirtualCascade( kinematics,
                                                   effectiveSquarkMassPointer,
                                                   coloredDecayer,
                                                   electroweakDecayer,
                                                   shortcut );
      cascades.push_back( currentCascade );
    }
  }

  electroweakCascadeSet::electroweakCascadeSet(
                                leptonAcceptanceParameterSet* const kinematics,
                                          particlePointer const coloredDecayer,
                                      particlePointer const electroweakDecayer,
                                         particlePointer const lighterScolored,
                                         inputHandler const* const shortcut ) :
    kinematics( kinematics ),
    coloredDecayer( coloredDecayer ),
    electroweakDecayer( electroweakDecayer ),
    lighterScolored( lighterScolored ),
    shortcut( shortcut )
  // this is the version for vector bosons.
  {
    effectiveSquarkMassHolder* effectiveSquarkMassPointer;
    if( shortcut->getGluino() == coloredDecayer )
    {
      effectiveSquarkMassPointer = kinematics;
    }
    else
    {
      effectiveSquarkMassPointer
      = shortcut->getOnShellEffectiveSquarkMass( coloredDecayer );
    }

    if( CppSLHA::PDG_code::W_plus == electroweakDecayer->get_PDG_code() )
    {
      currentCascade = new scoloredToWPlusScoloredCascade( kinematics,
                                                    effectiveSquarkMassPointer,
                                                           coloredDecayer,
                                                           electroweakDecayer,
                                                           lighterScolored,
                                                           shortcut );
      cascades.push_back( currentCascade );
    }
    else if( CppSLHA::PDG_code::Z == electroweakDecayer->get_PDG_code() )
    {
      currentCascade = new scoloredToZPlusScoloredCascade( kinematics,
                                                    effectiveSquarkMassPointer,
                                                           coloredDecayer,
                                                           electroweakDecayer,
                                                           lighterScolored,
                                                           shortcut );
      cascades.push_back( currentCascade );
    }
  }

  electroweakCascadeSet::~electroweakCascadeSet()
  {
    for( std::vector< electroweakCascade* >::iterator
         deletionIterator = cascades.begin();
         cascades.end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
  }



  electroweakCascadeHandler::electroweakCascadeHandler(
                      leptonAcceptancesForOneBeamEnergy* const kinematicsTable,
                                                        int const beamEnergy,
                                         inputHandler const* const shortcut ) :
    kinematicsTable( kinematicsTable ),
    beamEnergy( beamEnergy ),
    shortcut( shortcut )
  {
    // just an initialization list.
  }

  electroweakCascadeHandler::~electroweakCascadeHandler()
  {
    for( std::vector< electroweakCascadeSet* >::iterator
         deletionIterator( electroweakinoCascadeSets.begin() );
         electroweakinoCascadeSets.end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
    for( std::vector< electroweakCascadeSet* >::iterator
         deletionIterator( bosonCascadeSets.begin() );
         bosonCascadeSets.end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
  }


  electroweakCascadeSet*
  electroweakCascadeHandler::getCascadeSet(
                                          particlePointer const coloredDecayer,
                                     particlePointer const electroweakDecayer )
  /* this looks to see if it already has an electroweakCascadeSet
   * corresponding to the requested pairing, & if it does, it returns a
   * pointer to it, & if it doesn't, it constructs a new
   * electroweakCascadeSet & returns a pointer to that.
   */
  {
    electroweakCascadeSet* returnPointer( NULL );
    for( std::vector< electroweakCascadeSet* >::iterator
         cascadeIterator( electroweakinoCascadeSets.begin() );
         electroweakinoCascadeSets.end() > cascadeIterator;
         ++cascadeIterator )
    {
      if( (*cascadeIterator)->isEquivalent( coloredDecayer,
                                            electroweakDecayer ) )
      {
        returnPointer = *cascadeIterator;
        cascadeIterator = electroweakinoCascadeSets.end();
      }
    }
    if( NULL == returnPointer )
    {
      returnPointer
      = new electroweakCascadeSet( kinematicsTable->getParameterSets(
                                             coloredDecayer )->getParameterSet(
                                                          electroweakDecayer ),
                                   coloredDecayer,
                                   electroweakDecayer,
                                   shortcut );
      electroweakinoCascadeSets.push_back( returnPointer );
    }
    return returnPointer;
  }

  electroweakCascadeSet*
  electroweakCascadeHandler::getCascadeSet(
                                          particlePointer const coloredDecayer,
                                        particlePointer const electroweakBoson,
                                        particlePointer const lighterScolored )
  /* this looks to see if it already has an electroweakCascadeSet
   * corresponding to the requested pairing, & if it does, it returns a
   * pointer to it, & if it doesn't, it constructs a new
   * electroweakCascadeSet & returns a pointer to that.
   */
  {
    electroweakCascadeSet* returnPointer( NULL );
    for( std::vector< electroweakCascadeSet* >::iterator
         cascadeIterator( bosonCascadeSets.begin() );
         bosonCascadeSets.end() > cascadeIterator;
         ++cascadeIterator )
    {
      if( (*cascadeIterator)->isEquivalent( coloredDecayer,
                                            electroweakBoson,
                                            lighterScolored ) )
      {
        returnPointer = *cascadeIterator;
        cascadeIterator = bosonCascadeSets.end();
      }
    }
    if( NULL == returnPointer )
    {
      returnPointer
      = new electroweakCascadeSet( kinematicsTable->getParameterSets(
                                             coloredDecayer )->getParameterSet(
                                                   shortcut->getNeutralinoOne()
                                   /* it shouldn't really matter which
                                    *  neutralino is used here, just as long as
                                    *   it's lighter than coloredDecayer */ ),
                                   coloredDecayer,
                                   electroweakDecayer,
                                   lighterScolored,
                                   shortcut );
      bosonCascadeSets.push_back( returnPointer );
    }
    return returnPointer;
  }

}  // end of LHC_FASER namespace.
