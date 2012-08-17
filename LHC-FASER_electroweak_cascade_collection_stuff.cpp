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
 *      LHC-FASER_base_electroweak_cascade_stuff.hpp
 *      LHC-FASER_base_electroweak_cascade_stuff.cpp
 *      LHC-FASER_base_kinematics_stuff.hpp
 *      LHC-FASER_base_kinematics_stuff.cpp
 *      LHC-FASER_base_lepton_distribution_stuff.hpp
 *      LHC-FASER_base_lepton_distribution_stuff.cpp
 *      LHC-FASER_charged_electroweak_cascade_stuff.hpp
 *      LHC-FASER_charged_electroweak_cascade_stuff.cpp
 *      LHC-FASER_cross-section_stuff.hpp
 *      LHC-FASER_cross-section_stuff.cpp
 *      LHC-FASER_derived_lepton_distributions.hpp
 *      LHC-FASER_derived_lepton_distributions.cpp
 *      LHC-FASER_electroweak_cascade_collection_stuff.hpp
 *      LHC-FASER_electroweak_cascade_collection_stuff.cpp
 *      LHC-FASER_full_cascade_stuff.hpp
 *      LHC-FASER_full_cascade_stuff.cpp
 *      LHC-FASER_global_stuff.hpp
 *      LHC-FASER_global_stuff.cpp
 *      LHC-FASER_input_handling_stuff.hpp
 *      LHC-FASER_input_handling_stuff.cpp
 *      LHC-FASER_jet_kinematics_stuff.hpp
 *      LHC-FASER_jet_kinematics_stuff.cpp
 *      LHC-FASER_lepton_kinematics_stuff.hpp
 *      LHC-FASER_lepton_kinematics_stuff.cpp
 *      LHC-FASER_neutral_electroweak_cascade_stuff.hpp
 *      LHC-FASER_neutral_electroweak_cascade_stuff.cpp
 *      LHC-FASER_signal_calculator_stuff.hpp
 *      LHC-FASER_signal_calculator_stuff.cpp
 *      LHC-FASER_signal_data_collection_stuff.hpp
 *      LHC-FASER_signal_data_collection_stuff.cpp
 *      LHC-FASER_sparticle_decay_stuff.hpp
 *      LHC-FASER_sparticle_decay_stuff.cpp
 *      LHC-FASER_template_classes.hpp
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
                           leptonAcceptanceParameterSet* const kinematicsTable,
                                          particlePointer const coloredDecayer,
                                      particlePointer const electroweakDecayer,
                                       inputHandler const* const inputShortcut,
                      effectiveSquarkMassHolder* effectiveSquarkMassPointer ) :
    kinematicsTable( kinematicsTable ),
    inputShortcut( inputShortcut ),
    coloredDecayer( coloredDecayer ),
    electroweakDecayer( electroweakDecayer ),
    lighterScolored( NULL )
  // this is the version for electroweakinos.
  {
    if( NULL == effectiveSquarkMassPointer )
    {
      if( coloredDecayer == inputShortcut->getGluino() )
      {
        effectiveSquarkMassPointer = kinematicsTable;
      }
      else
      {
        effectiveSquarkMassPointer
        = inputShortcut->getOnShellEffectiveSquarkMass( coloredDecayer );
      }
    }

    if( CppSLHA::PDG_code::neutralino_one
        == electroweakDecayer->get_PDG_code() )
    {
      currentCascade = new lightestNeutralinoCascade( inputShortcut );
      cascadeSet.push_back( currentCascade );
    }
    else if( inputShortcut->isIn( electroweakDecayer->get_PDG_code(),
                                  inputShortcut->getUnstableNeutralinos() ) )
    {
      currentCascade = new neutralinoToSemuCascade( kinematicsTable,
                                                    effectiveSquarkMassPointer,
                                                    coloredDecayer,
                                                    electroweakDecayer,
                                                inputShortcut->getSelectronL(),
                                                    inputShortcut );
      cascadeSet.push_back( currentCascade );
      currentCascade = new neutralinoToSemuCascade( kinematicsTable,
                                                    effectiveSquarkMassPointer,
                                                    coloredDecayer,
                                                    electroweakDecayer,
                                                inputShortcut->getSelectronR(),
                                                    inputShortcut );
      cascadeSet.push_back( currentCascade );
      currentCascade = new neutralinoToSemuCascade( kinematicsTable,
                                                    effectiveSquarkMassPointer,
                                                    coloredDecayer,
                                                    electroweakDecayer,
                                                    inputShortcut->getSmuonL(),
                                                    inputShortcut );
      cascadeSet.push_back( currentCascade );
      currentCascade = new neutralinoToSemuCascade( kinematicsTable,
                                                    effectiveSquarkMassPointer,
                                                    coloredDecayer,
                                                    electroweakDecayer,
                                                    inputShortcut->getSmuonR(),
                                                    inputShortcut );
      cascadeSet.push_back( currentCascade );
      currentCascade = new neutralinoToStauCascade( kinematicsTable,
                                                    effectiveSquarkMassPointer,
                                                    coloredDecayer,
                                                    electroweakDecayer,
                                                   inputShortcut->getStauOne(),
                                                    inputShortcut );
      cascadeSet.push_back( currentCascade );
      currentCascade = new neutralinoToStauCascade( kinematicsTable,
                                                    effectiveSquarkMassPointer,
                                                    coloredDecayer,
                                                    electroweakDecayer,
                                                   inputShortcut->getStauTwo(),
                                                    inputShortcut );
      cascadeSet.push_back( currentCascade );
      currentCascade = new neutralinoToZCascade( kinematicsTable,
                                                 effectiveSquarkMassPointer,
                                                 coloredDecayer,
                                                 electroweakDecayer,
                                                 inputShortcut );
      cascadeSet.push_back( currentCascade );
      for( std::vector< particlePointer >::const_iterator
           bosonIterator(
             inputShortcut->getNeutralEwsbScalarsAndPseudoscalars()->begin() );
           inputShortcut->getNeutralEwsbScalarsAndPseudoscalars()->end()
           > bosonIterator;
           ++bosonIterator )
      {
        currentCascade = new neutralinoToHiggsCascade( kinematicsTable,
                                                    effectiveSquarkMassPointer,
                                                       coloredDecayer,
                                                       electroweakDecayer,
                                                       *bosonIterator,
                                                       inputShortcut );
        cascadeSet.push_back( currentCascade );
      }
      currentCascade = new neutralinoVirtualCascade( kinematicsTable,
                                                    effectiveSquarkMassPointer,
                                                     coloredDecayer,
                                                     electroweakDecayer,
                                                     inputShortcut );
      cascadeSet.push_back( currentCascade );
    }
    else if( inputShortcut->isIn( electroweakDecayer->get_PDG_code(),
                                  inputShortcut->getCharginos() ) )
    {
      currentCascade = new charginoToSemuOrEmuSnuCascade( kinematicsTable,
                                                    effectiveSquarkMassPointer,
                                                          coloredDecayer,
                                                          electroweakDecayer,
                                                inputShortcut->getSelectronL(),
                                                          inputShortcut,
                                                          false );
      cascadeSet.push_back( currentCascade );
      currentCascade = new charginoToSemuOrEmuSnuCascade( kinematicsTable,
                                                    effectiveSquarkMassPointer,
                                                          coloredDecayer,
                                                          electroweakDecayer,
                                                inputShortcut->getSelectronR(),
                                                          inputShortcut,
                                                          false );
      cascadeSet.push_back( currentCascade );
      currentCascade = new charginoToSemuOrEmuSnuCascade( kinematicsTable,
                                                    effectiveSquarkMassPointer,
                                                          coloredDecayer,
                                                          electroweakDecayer,
                                        inputShortcut->getElectronSneutrinoL(),
                                                          inputShortcut,
                                                          true );
      cascadeSet.push_back( currentCascade );
      currentCascade = new charginoToSemuOrEmuSnuCascade( kinematicsTable,
                                                    effectiveSquarkMassPointer,
                                                          coloredDecayer,
                                                          electroweakDecayer,
                                        inputShortcut->getElectronSneutrinoR(),
                                                          inputShortcut,
                                                          true );
      cascadeSet.push_back( currentCascade );
      currentCascade = new charginoToSemuOrEmuSnuCascade( kinematicsTable,
                                                    effectiveSquarkMassPointer,
                                                          coloredDecayer,
                                                          electroweakDecayer,
                                                    inputShortcut->getSmuonL(),
                                                          inputShortcut,
                                                          false );
      cascadeSet.push_back( currentCascade );
      currentCascade = new charginoToSemuOrEmuSnuCascade( kinematicsTable,
                                                    effectiveSquarkMassPointer,
                                                          coloredDecayer,
                                                          electroweakDecayer,
                                                    inputShortcut->getSmuonR(),
                                                          inputShortcut,
                                                          false );
      cascadeSet.push_back( currentCascade );
      currentCascade = new charginoToSemuOrEmuSnuCascade( kinematicsTable,
                                                    effectiveSquarkMassPointer,
                                                          coloredDecayer,
                                                          electroweakDecayer,
                                            inputShortcut->getMuonSneutrinoL(),
                                                          inputShortcut,
                                                          true );
      cascadeSet.push_back( currentCascade );
      currentCascade = new charginoToSemuOrEmuSnuCascade( kinematicsTable,
                                                    effectiveSquarkMassPointer,
                                                          coloredDecayer,
                                                          electroweakDecayer,
                                            inputShortcut->getMuonSneutrinoR(),
                                                          inputShortcut,
                                                          true );
      cascadeSet.push_back( currentCascade );
      currentCascade = new charginoToStauOrTauSnuCascade( kinematicsTable,
                                          effectiveSquarkMassPointer,
                                                          coloredDecayer,
                                                          electroweakDecayer,
                                                   inputShortcut->getStauOne(),
                                                          inputShortcut,
                                                          false );
      cascadeSet.push_back( currentCascade );
      currentCascade = new charginoToStauOrTauSnuCascade( kinematicsTable,
                                          effectiveSquarkMassPointer,
                                                          coloredDecayer,
                                                          electroweakDecayer,
                                                   inputShortcut->getStauTwo(),
                                                          inputShortcut,
                                                          false );
      cascadeSet.push_back( currentCascade );
      currentCascade = new charginoToStauOrTauSnuCascade( kinematicsTable,
                                          effectiveSquarkMassPointer,
                                                          coloredDecayer,
                                                          electroweakDecayer,
                                             inputShortcut->getTauSneutrinoL(),
                                                          inputShortcut,
                                                          true );
      cascadeSet.push_back( currentCascade );
      currentCascade = new charginoToStauOrTauSnuCascade( kinematicsTable,
                                          effectiveSquarkMassPointer,
                                                          coloredDecayer,
                                                          electroweakDecayer,
                                             inputShortcut->getTauSneutrinoR(),
                                                          inputShortcut,
                                                          true );
      cascadeSet.push_back( currentCascade );
      currentCascade = new charginoToWCascade( kinematicsTable,
                                               effectiveSquarkMassPointer,
                                               coloredDecayer,
                                               electroweakDecayer,
                                               inputShortcut );
      cascadeSet.push_back( currentCascade );
      currentCascade = new charginoToHiggsCascade( kinematicsTable,
                                                   effectiveSquarkMassPointer,
                                                   coloredDecayer,
                                                   electroweakDecayer,
                                         inputShortcut->getChargedEwsbScalar(),
                                                   inputShortcut );
      cascadeSet.push_back( currentCascade );
      currentCascade = new charginoVirtualCascade( kinematicsTable,
                                                   effectiveSquarkMassPointer,
                                                   coloredDecayer,
                                                   electroweakDecayer,
                                                   inputShortcut );
      cascadeSet.push_back( currentCascade );
    }
  }

  electroweakCascadeSet::electroweakCascadeSet(
                           leptonAcceptanceParameterSet* const kinematicsTable,
                                          particlePointer const coloredDecayer,
                                      particlePointer const electroweakDecayer,
                                         particlePointer const lighterScolored,
                                       inputHandler const* const inputShortcut,
                      effectiveSquarkMassHolder* effectiveSquarkMassPointer ) :
    kinematicsTable( kinematicsTable ),
    inputShortcut( inputShortcut ),
    coloredDecayer( coloredDecayer ),
    electroweakDecayer( electroweakDecayer ),
    lighterScolored( lighterScolored )
  // this is the version for bosons.
  {
    if( NULL == effectiveSquarkMassPointer )
    {
      if( coloredDecayer == inputShortcut->getGluino() )
      {
        effectiveSquarkMassPointer = kinematicsTable;
      }
      else
      {
        effectiveSquarkMassPointer
        = inputShortcut->getOnShellEffectiveSquarkMass( coloredDecayer );
      }
    }

    if( CppSLHA::PDG_code::W_plus == electroweakDecayer->get_PDG_code() )
    {
      currentCascade = new scoloredToWPlusScoloredCascade( kinematicsTable,
                                                    effectiveSquarkMassPointer,
                                                           coloredDecayer,
                                                           electroweakDecayer,
                                                           lighterScolored,
                                                           inputShortcut );
      cascadeSet.push_back( currentCascade );
    }
    else if( CppSLHA::PDG_code::Z == electroweakDecayer->get_PDG_code() )
    {
      currentCascade = new scoloredToZPlusScoloredCascade( kinematicsTable,
                                                    effectiveSquarkMassPointer,
                                                           coloredDecayer,
                                                           lighterScolored,
                                                           inputShortcut );
      cascadeSet.push_back( currentCascade );
    }
    else if( CppSLHA::PDG_code::charged_EWSB_scalar
             == electroweakDecayer->get_PDG_code() )
    {
      currentCascade = new scoloredToHiggsPlusScoloredCascade( kinematicsTable,
                                                    effectiveSquarkMassPointer,
                                                               coloredDecayer,
                                                               lighterScolored,
                                                            electroweakDecayer,
                                                               inputShortcut );
      cascadeSet.push_back( currentCascade );
    }
    else if( inputShortcut->isIn( electroweakDecayer->get_PDG_code(),
                     inputShortcut->getNeutralEwsbScalarsAndPseudoscalars() ) )
    {
      currentCascade = new scoloredToHiggsPlusScoloredCascade( kinematicsTable,
                                                    effectiveSquarkMassPointer,
                                                               coloredDecayer,
                                                               lighterScolored,
                                                            electroweakDecayer,
                                                               inputShortcut );
      cascadeSet.push_back( currentCascade );
    }
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "cascadeSet.size() = " << cascadeSet.size();
    std::cout << std::endl;**/
  }

  electroweakCascadeSet::~electroweakCascadeSet()
  {
    for( std::vector< electroweakCascade* >::iterator
         deletionIterator( cascadeSet.begin() );
         cascadeSet.end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
  }



  electroweakCascadesForOneBeamEnergy::electroweakCascadesForOneBeamEnergy(
                               leptonAcceptanceHandler* const kinematicsSource,
                                                          int const beamEnergy,
                                    inputHandler const* const inputShortcut ) :
    kinematicsTable( kinematicsSource->getLeptonAcceptancesForOneBeamEnergy(
                                                                beamEnergy ) ),
    beamEnergy( beamEnergy ),
    inputShortcut( inputShortcut ),
    electroweakinoCascadeSets(),
    bosonCascadeSets(),
    virtualSdownCascadeSets(),
    virtualTopCascadeSets(),
    currentCascadeSet( NULL )
  {
    // just an initialization list.
  }

  electroweakCascadesForOneBeamEnergy::~electroweakCascadesForOneBeamEnergy()
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
    for( std::vector< electroweakCascadeSet* >::iterator
         deletionIterator( virtualSdownCascadeSets.begin() );
         virtualSdownCascadeSets.end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
    for( std::vector< electroweakCascadeSet* >::iterator
         deletionIterator( virtualTopCascadeSets.begin() );
         virtualTopCascadeSets.end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
  }


  electroweakCascadeSet*
  electroweakCascadesForOneBeamEnergy::getCascadeSet(
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
                                                                 coloredDecayer
                                                            )->getParameterSet(
                                                          electroweakDecayer ),
                                   coloredDecayer,
                                   electroweakDecayer,
                                   inputShortcut );
      electroweakinoCascadeSets.push_back( returnPointer );
    }
    return returnPointer;
  }

  electroweakCascadeSet*
  electroweakCascadesForOneBeamEnergy::getCascadeSet(
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
                                              inputShortcut->getNeutralinoOne()
                                   /* it shouldn't really matter which
                                    *  neutralino is used here, just as long as
                                    *   it's lighter than coloredDecayer */ ),
                                   coloredDecayer,
                                   electroweakBoson,
                                   lighterScolored,
                                   inputShortcut );
      bosonCascadeSets.push_back( returnPointer );
    }
    return returnPointer;
  }

  electroweakCascadeSet*
  electroweakCascadesForOneBeamEnergy::getCascadeSet(
                                          particlePointer const coloredDecayer,
                                      particlePointer const electroweakDecayer,
                         effectiveSquarkMassHolder* const effectiveSquarkMass )
  /* this looks to see if it already has an electroweakCascadeSet
   * corresponding to the requested pairing, & if it does, it returns a
   * pointer to it, & if it doesn't, it constructs a new
   * electroweakCascadeSet & returns a pointer to that.
   */
  {
    electroweakCascadeSet* returnPointer( NULL );
    for( std::vector< electroweakCascadeSet* >::iterator
         cascadeIterator( virtualSdownCascadeSets.begin() );
         virtualSdownCascadeSets.end() > cascadeIterator;
         ++cascadeIterator )
    {
      if( (*cascadeIterator)->isEquivalent( coloredDecayer,
                                            electroweakDecayer ) )
      {
        returnPointer = *cascadeIterator;
        cascadeIterator = virtualSdownCascadeSets.end();
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
                                   inputShortcut,
                                   effectiveSquarkMass );
      virtualSdownCascadeSets.push_back( returnPointer );
    }
    return returnPointer;
  }

  electroweakCascadeSet*
  electroweakCascadesForOneBeamEnergy::getCascadeSet(
                                          particlePointer const coloredDecayer,
                                        particlePointer const electroweakBoson,
                                         particlePointer const lighterScolored,
                         effectiveSquarkMassHolder* const effectiveSquarkMass )
  /* this looks to see if it already has an electroweakCascadeSet
   * corresponding to the requested pairing, & if it does, it returns a
   * pointer to it, & if it doesn't, it constructs a new
   * electroweakCascadeSet & returns a pointer to that.
   */
  {
    electroweakCascadeSet* returnPointer( NULL );
    for( std::vector< electroweakCascadeSet* >::iterator
         cascadeIterator( virtualTopCascadeSets.begin() );
         virtualTopCascadeSets.end() > cascadeIterator;
         ++cascadeIterator )
    {
      if( (*cascadeIterator)->isEquivalent( coloredDecayer,
                                            electroweakBoson,
                                            lighterScolored ) )
      {
        returnPointer = *cascadeIterator;
        cascadeIterator = virtualTopCascadeSets.end();
      }
    }
    if( NULL == returnPointer )
    {
      returnPointer
      = new electroweakCascadeSet( kinematicsTable->getParameterSets(
                                             coloredDecayer )->getParameterSet(
                                              inputShortcut->getNeutralinoOne()
                                   /* it shouldn't really matter which
                                    *  neutralino is used here, just as long as
                                    *   it's lighter than coloredDecayer */ ),
                                   coloredDecayer,
                                   electroweakBoson,
                                   lighterScolored,
                                   inputShortcut,
                                   effectiveSquarkMass );
      virtualTopCascadeSets.push_back( returnPointer );
    }
    return returnPointer;
  }


  electroweakCascadeHandler::electroweakCascadeHandler(
                                    inputHandler const* const inputShortcut ) :
    kinematicsSource( inputShortcut,
                      inputShortcut->inspectPathToKinematicsGrids() ),
    inputShortcut( inputShortcut )
  {
    // just an initialization list.
  }

  electroweakCascadeHandler::~electroweakCascadeHandler()
  {
    for( std::vector< electroweakCascadesForOneBeamEnergy* >::iterator
         deletionIterator( cascadeSetAtBeamEnergies.begin() );
         cascadeSetAtBeamEnergies.end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
  }


  electroweakCascadesForOneBeamEnergy*
  electroweakCascadeHandler::getElectroweakCascadesForOneBeamEnergy(
                                                         int const beamEnergy )
  /* this looks to see if there is an existing
   * electroweakCascadesForOneBeamEnergy with the requested beam energy, & if
   * not, makes 1, & returns the pointer.
   */
  {
    electroweakCascadesForOneBeamEnergy* returnPointer( NULL );
    // we look to see if we already have a leptonAcceptancesForOneBeamEnergy
    // for these acceptanceValues:
    for( std::vector< electroweakCascadesForOneBeamEnergy* >::iterator
         searchIterator( cascadeSetAtBeamEnergies.begin() );
         cascadeSetAtBeamEnergies.end() > searchIterator;
         ++searchIterator )
    {
      if( beamEnergy == (*searchIterator)->getBeamEnergy() )
      {
        returnPointer = *searchIterator;
        searchIterator = cascadeSetAtBeamEnergies.end();
      }
    }
    if( NULL == returnPointer )
      // if we do not already have a electroweakCascadeForOneBeamEnergy for
      // this beam energy, we make a new instance:
    {
      returnPointer
      = new electroweakCascadesForOneBeamEnergy( &kinematicsSource,
                                                 beamEnergy,
                                                 inputShortcut );
      cascadeSetAtBeamEnergies.push_back( returnPointer );
    }
    return returnPointer;
  }

}  // end of LHC_FASER namespace.
