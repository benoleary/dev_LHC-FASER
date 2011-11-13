/*
 * signal_testing.hpp
 *
 *  Created on: Oct 30, 2011
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *              Jonas Lindert (jonas.lindert@googlemail.com)
 *              Carsten Robens (carsten.robens@gmx.de)
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

#ifndef SIGNAL_TESTING_HPP_
#define SIGNAL_TESTING_HPP_

#include "basic_testing_stuff.hpp"
#include "input_testing.hpp"
#include "../LHC-FASER_signal_data_collection_stuff.hpp"
#include "../LHC-FASER_signal_calculator_stuff.hpp"

namespace LHC_FASER
{
  class signalTesting
  {
  public:
    signalTesting( basicStuffTesting* basicStuff,
                                    inputTesting* inputTester ) :
        basicStuff( basicStuff ),
        inputTester( inputTester ),
        testInputHandler( inputTester->getInput() )
    {
      // just an initialization list.
    }
    ~signalTesting()
    {
      // does nothing.
    }


    void
    performTest()
    {
      crossSectionHandler crossSections( testInputHandler );
      jetPlusMetAcceptanceHandler
      jetPlusMetAcceptances( testInputHandler );
      electroweakCascadeHandler
      electroweakCascadeSource( testInputHandler );
      fullCascadeSetFactory cascadeSets( testInputHandler,
                                         &electroweakCascadeSource );
      signalShortcuts*
      testSignalShortcutsPointer( new signalShortcuts( testInputHandler,
                                                       &crossSections,
                                                       &jetPlusMetAcceptances,
                                                       &cascadeSets ) );
      std::cout
      << std::endl
      << "testSignalShortcutsPointer->getUncertainty() = "
      << testSignalShortcutsPointer->getUncertainty();
      std::cout << std::endl;
      delete testSignalShortcutsPointer;
      std::cout
      << std::endl
      << "signalShortcuts class seems to be OK.";
      std::cout << std::endl;

      signalShortcuts
      testSignalShortcuts( testInputHandler,
                           &crossSections,
                           &jetPlusMetAcceptances,
                           &cascadeSets );
      signalDefinitionSet*
      testSignalDefinitionSetPointer( new signalDefinitionSet(
                                                      &testSignalShortcuts ) );
      std::cout
      << std::endl
      << "testSignalDefinitionSetPointer->getBeamEnergy() = "
      << testSignalDefinitionSetPointer->getBeamEnergy();
      std::cout << std::endl;
      signalDefinitionSet
      testSignalDefinitionSet( testSignalDefinitionSetPointer );
      delete testSignalDefinitionSetPointer;
      std::cout
      << std::endl
      << "signalDefinitionSet class seems to be OK.";
      std::cout << std::endl;
      testSignalDefinitionSet.setBeamEnergy( 7 );

      particlePointer gluinoPointer( testInputHandler->getGluino() );
      particlePointer sdownLPointer( testInputHandler->getSdownL() );
      particlePointer sbottomOnePointer( testInputHandler->getSbottomOne() );
      particlePointer sbottomTwoPointer( testInputHandler->getSbottomTwo() );
      signedParticleShortcutPair gluinoWithSdownL( gluinoPointer,
                                                   true,
                                                   sdownLPointer,
                                                   true );
      signedParticleShortcutPair sbottomOneWithGluino( gluinoPointer,
                                                       true,
                                                       sbottomOnePointer,
                                                       true );
      signedParticleShortcutPair gluinoPair( gluinoPointer,
                                             true,
                                             gluinoPointer,
                                             true );
      signedParticleShortcutPair
      sbottomTwoWithAntisbottomTwo( sbottomTwoPointer,
                                    true,
                                    sbottomTwoPointer,
                                    false );
      signedParticleShortcutPair sdownLWithSdownL( sdownLPointer,
                                                   true,
                                                   sdownLPointer,
                                                   true );
      productionChannelPointerSet*
      testProductionChannelPointerSet( new productionChannelPointerSet(
                                                      &testSignalDefinitionSet,
                                                         &gluinoWithSdownL ) );
      std::cout
      << std::endl
      << "testProductionChannelPointerSet->getCrossSection() = "
      << testProductionChannelPointerSet->getCrossSection();
      std::cout << std::endl;
      delete testProductionChannelPointerSet;
      std::cout
      << std::endl
      << "deleted testProductionChannelPointerSet without problem.";
      std::cout << std::endl;
      productionChannelPointerSet
      gluinoWithSdownLPointerSet( &testSignalDefinitionSet,
                                  &gluinoWithSdownL );
      std::cout
      << std::endl
      << "gluinoWithSdownLPointerSet.getCrossSection() = "
      << gluinoWithSdownLPointerSet.getCrossSection();
      std::cout << std::endl;
      productionChannelPointerSet
      sbottomOneWithGluinoPointerSet( &testSignalDefinitionSet,
                                      &sbottomOneWithGluino );
      std::cout
      << std::endl
      << "sbottomOneWithGluinoPointerSet.getCrossSection() = "
      << sbottomOneWithGluinoPointerSet.getCrossSection();
      std::cout << std::endl;
      productionChannelPointerSet
      gluinoPairPointerSet( &testSignalDefinitionSet,
                            &gluinoPair );
      std::cout
      << std::endl
      << "gluinoPairPointerSet.getCrossSection() = "
      << gluinoPairPointerSet.getCrossSection();
      std::cout << std::endl;
      productionChannelPointerSet
      sbottomTwoWithAntisbottomTwoPointerSet( &testSignalDefinitionSet,
                                              &sbottomTwoWithAntisbottomTwo );
      std::cout
      << std::endl
      << "sbottomTwoWithAntisbottomTwoPointerSet.getCrossSection() = "
      << sbottomTwoWithAntisbottomTwoPointerSet.getCrossSection();
      std::cout << std::endl;
      productionChannelPointerSet
      sdownLWithSdownLPointerSet( &testSignalDefinitionSet,
                                  &sdownLWithSdownL );
      std::cout
      << std::endl
      << "sdownLWithSdownLPointerSet.getCrossSection() = "
      << sdownLWithSdownLPointerSet.getCrossSection();
      std::cout << std::endl;

      std::cout
      << std::endl
      << "productionChannelPointerSet class seems to be OK.";
      std::cout << std::endl;

      signalCalculatorClasses::reallyWrongCalculator*
      testReallyWrongCalculator(
                            new signalCalculatorClasses::reallyWrongCalculator(
                                                  &testSignalDefinitionSet ) );
      double signalValue;
      double uncertaintyFactor;
      std::cout
      << std::endl
      << "testReallyWrongCalculator->calculateValue( "
      << &signalValue << ", " << &uncertaintyFactor << " ) = "
      << testReallyWrongCalculator->calculateValue( &signalValue,
                                                    &uncertaintyFactor );
      std::cout
      << std::endl
      << "now signalValue = " << signalValue << ", uncertaintyFactor = "
      << uncertaintyFactor;
      std::cout << std::endl;

      std::cout << std::endl;
      delete testReallyWrongCalculator;
      std::cout
      << std::endl
      << "reallyWrongCalculator class seems to be OK.";
      std::cout << std::endl;

      testSignalDefinitionSet.setExcludedStandardModelProducts(
                                           testInputHandler->getNotInJets5() );
      testSignalDefinitionSet.setJetCut( 40.0 );
      testSignalDefinitionSet.setPrimaryLeptonCut( 20.0 );
      testSignalDefinitionSet.setSecondaryLeptonCut( 10.0 );
      signalHandler*
      testSignalHandlerPointer( new signalHandler( "sigmaBreakdownTest",
                                                   1000.0,
                                                  &testSignalDefinitionSet ) );
      std::cout
      << std::endl
      << "SPS1a: testSignalHandlerPointer->getValue() (in *fb*) = "
      << testSignalHandlerPointer->getValue();
      std::cout << std::endl;
      delete testSignalHandlerPointer;
      std::cout
      << std::endl
      << "sigmaBreakdownTest seems OK.";
      std::cout << std::endl;
      signalHandler
      testAtlasFourJetMetZeroLepton( "Atlas4jMET0l07TeV",
                                     1.0,
                                     &testSignalDefinitionSet );
      std::cout
      << std::endl
      << "SPS1a: testAtlasFourJetMetZeroLepton.getValue() (in *pb*) = "
      << testAtlasFourJetMetZeroLepton.getValue();
      std::cout << std::endl;
      std::cout
      << std::endl
      << "atlasFourJetMetZeroLepton seems OK.";
      std::cout << std::endl;
      signalHandler
      testAtlasFourJetMetZeroLeptonHighCut( "Atlas4jMET0l07TeVpTl50.0",
                                            1.0,
                                            &testSignalDefinitionSet );
      std::cout
      << std::endl
      << "SPS1a: testAtlasFourJetMetZeroLeptonHighCut.getValue() (in *pb*) = "
      << testAtlasFourJetMetZeroLeptonHighCut.getValue();
      std::cout << std::endl;
      std::cout
      << std::endl
      << "atlasFourJetMetZeroLepton seems OK.";
      std::cout << std::endl;
      signalHandler
      testAtlasThreeJetMetOneLepton( "Atlas3jMET1l07TeV",
                                     1.0,
                                     &testSignalDefinitionSet );
      std::cout
      << std::endl
      << "SPS1a: testAtlasThreeJetMetOneLepton.getValue() (in *pb*) = "
      << testAtlasThreeJetMetOneLepton.getValue();
      std::cout << std::endl;
      std::cout
      << std::endl
      << "atlasThreeJetMetOneLepton seems OK.";
      std::cout << std::endl;
      signalHandler
      testAtlasThreeJetMetOneLeptonHighCuts( "Atlas3jMET1l07TeVpTl60.8pTl18.5",
                                             1.0,
                                             &testSignalDefinitionSet );
      std::cout
      << std::endl
      << "SPS1a: testAtlasThreeJetMetOneLeptonHighCuts.getValue() (in *pb*) = "
      << testAtlasThreeJetMetOneLeptonHighCuts.getValue();
      std::cout << std::endl;
      std::cout
      << std::endl
      << "atlasThreeJetMetOneLeptonHighCuts seems OK.";
      std::cout << std::endl;

      std::cout
      << std::endl
      << "changing spectrum.";
      std::cout << std::endl;
      basicStuff->getSlha()->read_file( "SPS2_spectrum.out" );
      basicStuff->getReadier()->readyObserversForNewPoint();


      signalHandler testSigmaBreakdownTest( "sigmaBreakdownTest",
                                            1000.0,
                                            &testSignalDefinitionSet );
      std::cout
      << std::endl
      << "SPS2: testSigmaBreakdownTest.getValue() (in *fb*) = "
      << testSigmaBreakdownTest.getValue();
      std::cout << std::endl;
      std::cout
      << std::endl
      << "again sigmaBreakdownTest seems OK.";
      std::cout << std::endl;
      std::cout
      << std::endl
      << "SPS2: testAtlasFourJetMetZeroLepton.getValue() (in *pb*) = "
      << testAtlasFourJetMetZeroLepton.getValue();
      std::cout << std::endl;
      std::cout
      << std::endl
      << "again atlasFourJetMetZeroLepton seems OK.";
      std::cout << std::endl;
      std::cout
      << std::endl
      << "SPS2: testAtlasFourJetMetZeroLeptonHighCut.getValue() (in *pb*) = "
      << testAtlasFourJetMetZeroLeptonHighCut.getValue();
      std::cout << std::endl;
      std::cout
      << std::endl
      << "again atlasFourJetMetZeroLepton seems OK.";
      std::cout << std::endl;
      std::cout
      << std::endl
      << "SPS2: testAtlasThreeJetMetOneLepton.getValue() (in *pb*) = "
      << testAtlasThreeJetMetOneLepton.getValue();
      std::cout << std::endl;
      std::cout
      << std::endl
      << "again atlasThreeJetMetOneLepton seems OK.";
      std::cout << std::endl;
      std::cout
      << std::endl
      << "SPS2: testAtlasThreeJetMetOneLeptonHighCuts.getValue() (in *pb*) = "
      << testAtlasThreeJetMetOneLeptonHighCuts.getValue();
      std::cout << std::endl;
      std::cout
      << std::endl
      << "again atlasThreeJetMetOneLeptonHighCuts seems OK.";
      std::cout << std::endl;

      basicStuff->getSlha()->read_file( "SPS1a_spectrum.out" );
      basicStuff->getReadier()->readyObserversForNewPoint();
    }


  protected:
    basicStuffTesting* basicStuff;
    inputTesting* inputTester;
    inputHandler const* const testInputHandler;
  };

}  // end of LHC_FASER namespace

#endif /* SIGNAL_TESTING_HPP_ */
