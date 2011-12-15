/*
 * full_electroweak_cascade_testing.hpp
 *
 *  Created on: Oct 30, 2011
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

#ifndef FULL_ELECTROWEAK_CASCADE_TESTING_HPP_
#define FULL_ELECTROWEAK_CASCADE_TESTING_HPP_

#include "basic_testing_stuff.hpp"
#include "input_testing.hpp"
#include "../LHC-FASER_electroweak_cascade_collection_stuff.hpp"

namespace LHC_FASER
{
  class fullElectroweakCascadeTesting
  {
  public:
    fullElectroweakCascadeTesting( basicStuffTesting* basicStuff,
                                   inputTesting* inputTester ) :
        basicStuff( basicStuff ),
        testReadier( basicStuff->getReadier() ),
        inputTester( inputTester ),
        testInputHandler( inputTester->getInput() ),
        testAcceptanceCutSet(),
        noLeptonGapAcceptanceCutSet()
    {
      testAcceptanceCutSet.setJetCut( 40.0 );
      testAcceptanceCutSet.setPrimaryLeptonCut( 20.0 );
      testAcceptanceCutSet.setSecondaryLeptonCut( 10.0 );
      testAcceptanceCutSet.setExcludedStandardModelProducts(
                                           testInputHandler->getNotInJets4() );
      noLeptonGapAcceptanceCutSet.setJetCut( 40.0 );
      noLeptonGapAcceptanceCutSet.setPrimaryLeptonCut( 10.0 );
      noLeptonGapAcceptanceCutSet.setSecondaryLeptonCut( 10.0 );
      noLeptonGapAcceptanceCutSet.setExcludedStandardModelProducts(
                                           testInputHandler->getNotInJets4() );
    }
    ~fullElectroweakCascadeTesting()
    {
      // does nothing.
    }

    void
    performTest()
    {
      leptonAcceptanceHandler
      testLeptonAcceptanceHandler( testInputHandler,
                            testInputHandler->inspectPathToKinematicsGrids() );
      leptonAcceptanceParameterSet*
      cascadeKinematics(
                testLeptonAcceptanceHandler.getLeptonAcceptanceParameterSet( 7,
                                                testInputHandler->getScharmL(),
                                          testInputHandler->getNeutralinoTwo(),
                                                                           2.0,
                                                                      10.0 ) );
      electroweakCascadeSet*
      testEwinoCascadeSet( new electroweakCascadeSet( cascadeKinematics,
                                                 testInputHandler->getSdownL(),
                                          testInputHandler->getNeutralinoTwo(),
                                                      testInputHandler ) );
      electroweakCascadeSet*
      testBosonCascadeSet( new electroweakCascadeSet( cascadeKinematics,
                                             testInputHandler->getSbottomTwo(),
                                                  testInputHandler->getWPlus(),
                                                testInputHandler->getStopOne(),
                                                      testInputHandler ) );
      if( testEwinoCascadeSet->isEquivalent( testInputHandler->getSdownL(),
                                        testInputHandler->getNeutralinoTwo()) )
      {
        std::cout
        << std::endl
        << "testEwinoCascadeSet is equivalent for "
        << *(testInputHandler->getSdownL()->get_name()) << ", "
        << *(testInputHandler->getNeutralinoTwo()->get_name());
        std::cout << std::endl;
      }
      else
      {
        std::cout
        << std::endl
        << "testEwinoCascadeSet is not equivalent for "
        << *(testInputHandler->getSdownL()->get_name()) << ", "
        << *(testInputHandler->getNeutralinoTwo()->get_name()) << "?!";
        std::cout << std::endl;
      }
      if( testBosonCascadeSet->isEquivalent( testInputHandler->getSdownL(),
                                       testInputHandler->getNeutralinoTwo() ) )
      {
        std::cout
        << std::endl
        << "testBosonCascadeSet is equivalent for "
        << *(testInputHandler->getSdownL()->get_name()) << ", "
        << *(testInputHandler->getNeutralinoTwo()->get_name()) << "?!";
        std::cout << std::endl;
      }
      else
      {
        std::cout
        << std::endl
        << "testBosonCascadeSet is not equivalent for "
        << *(testInputHandler->getSdownL()->get_name()) << ", "
        << *(testInputHandler->getNeutralinoTwo()->get_name());
        std::cout << std::endl;
      }
      std::cout
      << std::endl
      << "testEwinoCascadeSet->getOssfMinusOsdf( ... ) = "
      << testEwinoCascadeSet->getOssfMinusOsdf( &testAcceptanceCutSet );
      std::cout << std::endl;
      printAcceptances( testEwinoCascadeSet,
                        &testAcceptanceCutSet );
      std::cout << std::endl;
      printAcceptances( testEwinoCascadeSet,
                        &noLeptonGapAcceptanceCutSet );
      std::cout
      << std::endl
      << "testEwinoCascadeSet->getOssfMinusOsdf( &testAcceptanceCutSet ) = "
      << testEwinoCascadeSet->getOssfMinusOsdf( &testAcceptanceCutSet );
      std::cout
      << std::endl
      << "testEwinoCascadeSet->getOssfMinusOsdf( &noLeptonGapAcceptanceCutSet"
      << " ) = "
      << testEwinoCascadeSet->getOssfMinusOsdf( &noLeptonGapAcceptanceCutSet );
      std::cout
      << std::endl
      << "testEwinoCascadeSet seems fine (remember that OSSF-OSDF requires"
      << " subtraction of OSDF :P )";
      std::cout << std::endl;

      std::cout << std::endl;
      std::cout
      << std::endl
      << "testBosonCascadeSet->getOssfMinusOsdf( ... ) = "
      << testBosonCascadeSet->getOssfMinusOsdf( &testAcceptanceCutSet );
      std::cout << std::endl;
      printAcceptances( testBosonCascadeSet,
                        &testAcceptanceCutSet );
      std::cout << std::endl;
      printAcceptances( testBosonCascadeSet,
                        &noLeptonGapAcceptanceCutSet );
      std::cout << std::endl;
      std::cout << std::endl;
      delete testBosonCascadeSet;
      delete testEwinoCascadeSet;

      std::cout
      << std::endl
      << "direct testEwinoCascadeSet & testBosonCascadeSet seem fine.";
      std::cout << std::endl;
      electroweakCascadeHandler
      testElectroweakCascadeHandler( testInputHandler );
      testBosonCascadeSet
      = testElectroweakCascadeHandler.getElectroweakCascadesForOneBeamEnergy( 7
                                                              )->getCascadeSet(
                                             testInputHandler->getSbottomTwo(),
                                                  testInputHandler->getWPlus(),
                                              testInputHandler->getStopOne() );
      testEwinoCascadeSet
      = testElectroweakCascadeHandler.getElectroweakCascadesForOneBeamEnergy( 7
                                                              )->getCascadeSet(
                                                 testInputHandler->getSdownL(),
                                        testInputHandler->getNeutralinoTwo() );
      std::cout << std::endl;
      printAcceptances( testEwinoCascadeSet,
                        &testAcceptanceCutSet );
      std::cout << std::endl;
      printAcceptances( testBosonCascadeSet,
                        &testAcceptanceCutSet );
      std::cout
      << std::endl
      << "checking caching:"
      << std::endl
      << "testEwinoCascadeSet->getAcceptance( &testAcceptanceCutSet,"
      << " 0, 1, 1, 0, 0 ) = "
      << testEwinoCascadeSet->getAcceptance( &testAcceptanceCutSet,
                                             0,
                                             1,
                                             1,
                                             0,
                                             0 )
      << std::endl
      << "testEwinoCascadeSet->getAcceptance( &testAcceptanceCutSet,"
      << " 0, 0, 0, 1, 1 ) = "
      << testEwinoCascadeSet->getAcceptance( &testAcceptanceCutSet,
                                             0,
                                             0,
                                             0,
                                             1,
                                             1 )
      << std::endl
      << "testEwinoCascadeSet->getAcceptance( &noLeptonGapAcceptanceCutSet,"
      << " 0, 0, 0, 1, 1 ) = "
      << testEwinoCascadeSet->getAcceptance( &noLeptonGapAcceptanceCutSet,
                                             0,
                                             0,
                                             0,
                                             1,
                                             1 );
      std::cout
      << std::endl
      << std::endl
      << "resetting readier.";
      std::cout << std::endl;
      testReadier->readyObserversForNewPoint();
      std::cout
      << std::endl
      << "testEwinoCascadeSet->getAcceptance( &testAcceptanceCutSet,"
      << " 0, 0, 0, 1, 1 ) = "
      << testEwinoCascadeSet->getAcceptance( &testAcceptanceCutSet,
                                             0,
                                             0,
                                             0,
                                             1,
                                             1 )
      << std::endl
      << "testEwinoCascadeSet->getAcceptance( &noLeptonGapAcceptanceCutSet,"
      << " 0, 1, 1, 0, 0 ) = "
      << testEwinoCascadeSet->getAcceptance( &noLeptonGapAcceptanceCutSet,
                                             0,
                                             1,
                                             1,
                                             0,
                                             0 );
      std::cout << std::endl
      << "testEwinoCascadeSet->getAcceptance( &testAcceptanceCutSet,"
      << " 0, 0, 0, 1, 1 ) = "
      << testEwinoCascadeSet->getAcceptance( &testAcceptanceCutSet,
                                             0,
                                             0,
                                             0,
                                             1,
                                             1 )
      << std::endl
      << "testEwinoCascadeSet->getAcceptance( &noLeptonGapAcceptanceCutSet,"
      << " 0, 1, 1, 0, 0 ) = "
      << testEwinoCascadeSet->getAcceptance( &noLeptonGapAcceptanceCutSet,
                                             0,
                                             1,
                                             1,
                                             0,
                                             0 );
      std::cout << std::endl;
      std::cout
      << std::endl
      << "testEwinoCascadeSet & testBosonCascadeSet from"
      << " testElectroweakCascadeHandler seem fine.";
      std::cout << std::endl;
    }


  protected:
    basicStuffTesting* basicStuff;
    readierForNewPoint* testReadier;
    inputTesting* inputTester;
    inputHandler const* const testInputHandler;
    acceptanceCutSet testAcceptanceCutSet;
    acceptanceCutSet noLeptonGapAcceptanceCutSet;

    void
    printAcceptances( electroweakCascadeSet* cascadeSet,
                      acceptanceCutSet const* const acceptanceCuts )
    {
      double acceptanceValue( 0.0 );
      std::cout
      << std::endl
      << "jet cut: " << acceptanceCuts->getJetCut()
      << ", primary lepton cut: " << acceptanceCuts->getPrimaryLeptonCut()
      << ", secondary lepton cut: " << acceptanceCuts->getSecondaryLeptonCut();

      for( int jetCounter( -1 );
           3 > jetCounter;
           ++jetCounter )
      {
        for( int eMinusCounter( -1 );
             3 > eMinusCounter;
             ++eMinusCounter )
        {
          for( int ePlusCounter( -1 );
               3 > ePlusCounter;
               ++ePlusCounter )
          {
            for( int muMinusCounter( -1 );
                 3 > muMinusCounter;
                 ++muMinusCounter )
            {
              for( int muPlusCounter( -1 );
                   3 > muPlusCounter;
                   ++muPlusCounter )
              {
                acceptanceValue = cascadeSet->getAcceptance( acceptanceCuts,
                                                             jetCounter,
                                                             eMinusCounter,
                                                             ePlusCounter,
                                                             muMinusCounter,
                                                             muPlusCounter );
                if( 0.0 != acceptanceValue )
                {
                  std::cout
                  << std::endl
                  << jetCounter << "j "
                  << eMinusCounter << "e- " << ePlusCounter << "e+ "
                  << muMinusCounter << "mu- " << muPlusCounter << "mu+ "
                  << "acceptance = " << acceptanceValue;
                }
              }
            }
          }
        }
      }
    }
  };

}  // end of LHC_FASER namespace

#endif /* FULL_ELECTROWEAK_CASCADE_TESTING_HPP_ */
