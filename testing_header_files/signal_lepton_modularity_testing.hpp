/*
 * signal_lepton_modularity_testing.hpp
 *
 *  Created on: Mar 7, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SIGNAL_LEPTON_MODULARITY_TESTING_HPP_
#define SIGNAL_LEPTON_MODULARITY_TESTING_HPP_

#include "basic_testing_stuff.hpp"
#include "input_testing.hpp"
#include "../LHC-FASER_signal_data_collection_stuff.hpp"
#include "../LHC-FASER_signal_calculator_stuff.hpp"

namespace LHC_FASER
{
  class signalLeptonModularityTesting
  {
  public:
    signalLeptonModularityTesting( basicStuffTesting* basicStuff,
                                    inputTesting* inputTester ) :
        basicStuff( basicStuff ),
        inputTester( inputTester ),
        testInputHandler( inputTester->getInput() )
    {
      // just an initialization list.
    }
    ~signalLeptonModularityTesting()
    {
      // does nothing.
    }


    void
    performTest()
    {

      something();
      needToGetValidPairsOfCasades();
      crossSectionHandler crossSections( testInputHandler );
      jetPlusMetAcceptanceHandler
      jetPlusMetAcceptances( testInputHandler );
      electroweakCascadeHandler
      electroweakCascadeSource( testInputHandler );
      fullCascadeSetFactory cascadeSets( testInputHandler,
                                         &electroweakCascadeSource );
      /**/
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
      testAtlasFourJetMetZeroLepton( "Atlas4jMET_0l_07TeV",
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
      testAtlasFourJetMetZeroLeptonHighCut( "Atlas4jMET_0l_07TeV_pTl50.0GeV",
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
      testAtlasThreeJetMetOneLepton( "Atlas3jMET_1l_07TeV",
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
      testAtlasThreeJetMetOneLeptonHighCuts(
                                   "Atlas3jMET_1l_07TeV_pTl60.8GeV_pTl18.5GeV",
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
      signalHandler
      testSameSignDilepton( "sameSignDilepton_07TeV",
                            1.0,
                            &testSignalDefinitionSet );
      std::cout
      << std::endl
      << "SPS1a: testSameSignDilepton.getValue() (in *pb*) = "
      << testSameSignDilepton.getValue();
      std::cout << std::endl;
      std::cout
      << std::endl
      << "sameSignDilepton seems OK.";
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
      << "atlasFourJetMetZeroLepton seems OK.";
      std::cout << std::endl;
      std::cout
      << std::endl
      << "SPS2: testAtlasFourJetMetZeroLeptonHighCut.getValue() (in *pb*) = "
      << testAtlasFourJetMetZeroLeptonHighCut.getValue();
      std::cout << std::endl;
      std::cout
      << std::endl
      << "atlasFourJetMetZeroLepton seems OK.";
      std::cout << std::endl;
      std::cout
      << std::endl
      << "SPS2: testAtlasThreeJetMetOneLepton.getValue() (in *pb*) = "
      << testAtlasThreeJetMetOneLepton.getValue();
      std::cout << std::endl;
      std::cout
      << std::endl
      << "atlasThreeJetMetOneLepton seems OK.";
      std::cout << std::endl;
      std::cout
      << std::endl
      << "SPS2: testAtlasThreeJetMetOneLeptonHighCuts.getValue() (in *pb*) = "
      << testAtlasThreeJetMetOneLeptonHighCuts.getValue();
      std::cout << std::endl;
      std::cout
      << std::endl
      << "atlasThreeJetMetOneLeptonHighCuts seems OK.";
      std::cout << std::endl;
      std::cout
      << std::endl
      << "SPS2: testSameSignDilepton.getValue() (in *pb*) = "
      << testSameSignDilepton.getValue();
      std::cout << std::endl;
      std::cout
      << std::endl
      << "sameSignDilepton seems OK.";
      std::cout << std::endl;
      /**/
      signalShortcuts
      convolutedSignalShortcuts( testInputHandler,
                                 &crossSections,
                                 &jetPlusMetAcceptances,
                                 &cascadeSets );
      signalDefinitionSet
      convolutedSignalDefinitionSet( &convolutedSignalShortcuts );
      convolutedSignalDefinitionSet.setExcludedStandardModelProducts(
                                           testInputHandler->getNotInJets5() );
      convolutedSignalDefinitionSet.setJetCut( 40.0 );
      convolutedSignalDefinitionSet.setPrimaryLeptonCut( 20.0 );
      convolutedSignalDefinitionSet.setSecondaryLeptonCut( 10.0 );
      signalHandler convolutedSigmaBreakdownTest( "sigmaBreakdownTest",
                                                  1000.0,
                                              &convolutedSignalDefinitionSet );
      basicStuff->getSlha()->read_file( "convoluted_spectrum.out" );
      basicStuff->getReadier()->readyObserversForNewPoint();
      std::cout
      << std::endl
      << "convoluted: convolutedSigmaBreakdownTest.getValue() (in *fb*) = "
      << convolutedSigmaBreakdownTest.getValue();
      std::cout << std::endl;
      std::cout
      << std::endl
      << "don't worry about the BRs adding up to over 100% for the convoluted"
      << " point from sigmaBreakdownTest, since it overcounts channels where"
      << " the electroweakino can decay back into a colored cascade, for"
      << " example (so the assumption of sigmaBreakdownTest that the EW decays"
      << " of the charginos are 100% is wrong).";
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


#endif /* SIGNAL_LEPTON_MODULARITY_TESTING_HPP_ */
