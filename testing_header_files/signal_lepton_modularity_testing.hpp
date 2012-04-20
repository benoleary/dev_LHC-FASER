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
      leptonAcceptanceStyle::fullySpecified const testFullySpecified( 1,
                                                                      0,
                                                                      0,
                                                                      1 );
      leptonAcceptanceStyle::noLeptonCutNorExtraJetCut const
      testNoLeptonCutNorExtraJetCut;
      leptonAcceptanceStyle::chargeAndFlavorSummed const
      testChargeAndFlavorSummed( 2 );
      leptonAcceptanceStyle::chargeAndFlavorSummed const singleLepton( 1 );
      leptonAcceptanceStyle::chargeSummed const testChargeSummed( 0,
                                                                  1 );
      leptonAcceptanceStyle::chargeSummed const singleMuon( 0,
                                                            1 );
      leptonAcceptanceStyle::flavorSummed const testFlavorSummed( 0,
                                                                  1 );
      leptonAcceptanceStyle::flavorSummed const singleNegative( 1,
                                                                0 );
      leptonAcceptanceStyle::flavorSummed const singlePositive( 0,
                                                                1 );
      leptonAcceptanceStyle::ossfMinusOsdf const testOssfMinusOsdf;
      leptonAcceptanceStyle::sameSignSameFlavor const testSameSignSameFlavor;


      crossSectionHandler crossSections( testInputHandler );
      jetPlusMetAcceptanceHandler
      jetPlusMetAcceptances( testInputHandler );
      electroweakCascadeHandler
      electroweakCascadeSource( testInputHandler );
      fullCascadeSetFactory cascadeSets( testInputHandler,
                                         &electroweakCascadeSource );
      particlePointer sdownLPointer( testInputHandler->getSdownL() );
      std::vector< fullCascade* >*
      sdownLCascades( cascadeSets.getFullCascadeSet( sdownLPointer,
                                                     7 )->getOpenCascades() );
      particlePointer gluinoPointer( testInputHandler->getGluino() );
      std::vector< fullCascade* >*
      gluinoCascades( cascadeSets.getFullCascadeSet( gluinoPointer,
                                                     7 )->getOpenCascades() );
      signedParticleShortcutPair gluinoWithSdownL( gluinoPointer,
                                                   true,
                                                   sdownLPointer,
                                                   true );
      signalShortcuts testSignalShortcuts( testInputHandler,
                                           &crossSections,
                                           &jetPlusMetAcceptances,
                                           &cascadeSets );
      signalDefinitionSet testSignalDefinitionSet( &testSignalShortcuts );
      testSignalDefinitionSet.setBeamEnergy( 7 );
      testSignalDefinitionSet.setJetCut( 40.0 );
      testSignalDefinitionSet.setPrimaryLeptonCut( 20.0 );
      testSignalDefinitionSet.setSecondaryLeptonCut( 10.0 );

      basicStuff->getSlha()->read_file( "SPS1a_spectrum.out" );
      basicStuff->getReadier()->readyObserversForNewPoint();

      for( int sdownLIndex( sdownLCascades->size() - 1 );
           0 <= sdownLIndex;
           --sdownLIndex )
      {
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout
        << std::endl
        << "sdownLIndex = " << sdownLIndex;
        std::cout
        << std::endl
        << sdownLCascades->at( sdownLIndex )->getAsString();
        std::cout << std::endl;
        std::cout << std::endl;

        for( int gluinoIndex( gluinoCascades->size() - 1 );
             0 <= gluinoIndex;
             --gluinoIndex )
        {
          std::cout << std::endl;
          std::cout << std::endl;
          std::cout
          << std::endl
          << "gluinoIndex = " << gluinoIndex;
          std::cout
          << std::endl
          << gluinoCascades->at( gluinoIndex )->getAsString();
          std::cout << std::endl;

          for( int numberOfJets( 0 );
               5 > numberOfJets;
               ++numberOfJets )
          {
            std::cout << std::endl;
            std::cout
            << std::endl
            << "testFullySpecified.withExactlyNJets( ..., " << numberOfJets
            << ", ... ) = "
            << testFullySpecified.withExactlyNJets( &testSignalDefinitionSet,
                                                     numberOfJets,
                                             gluinoCascades->at( gluinoIndex ),
                                                     true,
                                             sdownLCascades->at( sdownLIndex ),
                                                     true );
            std::cout
            << std::endl
            << "testNoLeptonCutNorExtraJetCut.withExactlyNJets( ..., "
            << numberOfJets << ", ... ) = "
            << testNoLeptonCutNorExtraJetCut.withExactlyNJets(
                                                      &testSignalDefinitionSet,
                                                               numberOfJets,
                                             gluinoCascades->at( gluinoIndex ),
                                                               true,
                                             sdownLCascades->at( sdownLIndex ),
                                                               true );

            std::cout
            << std::endl
            << "testChargeAndFlavorSummed.withExactlyNJets( ..., "
            << numberOfJets << ", ... ) = "
            << testChargeAndFlavorSummed.withExactlyNJets(
                                                      &testSignalDefinitionSet,
                                                           numberOfJets,
                                             gluinoCascades->at( gluinoIndex ),
                                                           true,
                                             sdownLCascades->at( sdownLIndex ),
                                                           true );

            std::cout
            << std::endl
            << "singleLepton.withExactlyNJets( ..., "
            << numberOfJets << ", ... ) = "
            << singleLepton.withExactlyNJets( &testSignalDefinitionSet,
                                              numberOfJets,
                                             gluinoCascades->at( gluinoIndex ),
                                              true,
                                             sdownLCascades->at( sdownLIndex ),
                                              true );

            std::cout
            << std::endl
            << "testChargeSummed.withExactlyNJets( ..., " << numberOfJets
            << ", ... ) = "
            << testChargeSummed.withExactlyNJets( &testSignalDefinitionSet,
                                                  numberOfJets,
                                             gluinoCascades->at( gluinoIndex ),
                                                  true,
                                             sdownLCascades->at( sdownLIndex ),
                                                  true );

            std::cout
            << std::endl
            << "singleMuon.withExactlyNJets( ..., " << numberOfJets
            << ", ... ) = "
            << singleMuon.withExactlyNJets( &testSignalDefinitionSet,
                                            numberOfJets,
                                            gluinoCascades->at( gluinoIndex ),
                                            true,
                                            sdownLCascades->at( sdownLIndex ),
                                            true );

            std::cout
            << std::endl
            << "testFlavorSummed.withExactlyNJets( ..., " << numberOfJets
            << ", ... ) = "
            << testFlavorSummed.withExactlyNJets( &testSignalDefinitionSet,
                                                  numberOfJets,
                                             gluinoCascades->at( gluinoIndex ),
                                                  true,
                                             sdownLCascades->at( sdownLIndex ),
                                                  true );

            std::cout
            << std::endl
            << "singleNegative.withExactlyNJets( ..., " << numberOfJets
            << ", ... ) = "
            << singleNegative.withExactlyNJets( &testSignalDefinitionSet,
                                                numberOfJets,
                                             gluinoCascades->at( gluinoIndex ),
                                                true,
                                             sdownLCascades->at( sdownLIndex ),
                                                true );

            std::cout
            << std::endl
            << "singlePositive.withExactlyNJets( ..., " << numberOfJets
            << ", ... ) = "
            << singlePositive.withExactlyNJets( &testSignalDefinitionSet,
                                                numberOfJets,
                                             gluinoCascades->at( gluinoIndex ),
                                                true,
                                             sdownLCascades->at( sdownLIndex ),
                                                true );

            std::cout
            << std::endl
            << "testOssfMinusOsdf.withExactlyNJets( ..., " << numberOfJets
            << ", ... ) = "
            << testOssfMinusOsdf.withExactlyNJets( &testSignalDefinitionSet,
                                                   numberOfJets,
                                             gluinoCascades->at( gluinoIndex ),
                                                   true,
                                             sdownLCascades->at( sdownLIndex ),
                                                   true );

            std::cout
            << std::endl
            << "testSameSignSameFlavor.withExactlyNJets( ..., " << numberOfJets
            << ", ... ) = "
            << testSameSignSameFlavor.withExactlyNJets(
                                                      &testSignalDefinitionSet,
                                                        numberOfJets,
                                             gluinoCascades->at( gluinoIndex ),
                                                        true,
                                             sdownLCascades->at( sdownLIndex ),
                                                        true );
          }
        }
      }

      fullCascade* gluinoSbottomOneNeutralinoOne( NULL );
      fullCascade* gluinoSbottomOneNeutralinoTwo( NULL );
      fullCascade* gluinoSbottomOneNeutralinoThree( NULL );
      particlePointer sbottomOnePointer( testInputHandler->getSbottomOne() );
      particlePointer
      neutralinoOnePointer( testInputHandler->getNeutralinoOne() );
      particlePointer
      neutralinoTwoPointer( testInputHandler->getNeutralinoTwo() );
      particlePointer
      neutralinoThreePointer( testInputHandler->getNeutralinoThree() );
      std::vector< fullCascade::particleWithInt* > const* cascadeDefiner;

      std::cout
      << std::endl
      << "sbottomOnePointer = " << sbottomOnePointer;
      std::cout << std::endl;
      std::cout
      << "sbottomOnePointer->get_name() = " << sbottomOnePointer->get_name();
      std::cout << std::endl;
      std::cout
      << "*(sbottomOnePointer->get_name()) = "
      << *(sbottomOnePointer->get_name());
      std::cout << std::endl;
      std::cout
      << std::endl
      << "gluinoPointer = " << gluinoPointer;
      std::cout << std::endl;
      std::cout
      << "gluinoPointer->get_name() = " << gluinoPointer->get_name();
      std::cout << std::endl;
      std::cout
      << "*(gluinoPointer->get_name()) = "
      << *(gluinoPointer->get_name());
      std::cout << std::endl;

      for( int gluinoIndex( gluinoCascades->size() - 1 );
           0 <= gluinoIndex;
           --gluinoIndex )
      {
        cascadeDefiner
        = gluinoCascades->at( gluinoIndex )->getCascadeDefiner();
        if( 2 == cascadeDefiner->size() )
        {
          if( sbottomOnePointer == cascadeDefiner->back()->first )
          {
            if( neutralinoOnePointer == cascadeDefiner->front()->first )
            {
              gluinoSbottomOneNeutralinoOne
              = gluinoCascades->at( gluinoIndex );
            }
            else if( neutralinoTwoPointer == cascadeDefiner->front()->first )
            {
              gluinoSbottomOneNeutralinoTwo
              = gluinoCascades->at( gluinoIndex );
            }
            else if( neutralinoThreePointer == cascadeDefiner->front()->first )
            {
              gluinoSbottomOneNeutralinoThree
              = gluinoCascades->at( gluinoIndex );
            }
          }
        }
      }

      std::cout
      << std::endl
      << "gluinoSbottomOneNeutralinoOne = " << gluinoSbottomOneNeutralinoOne
      << ", gluinoSbottomOneNeutralinoTwo = " << gluinoSbottomOneNeutralinoTwo
      << ", gluinoSbottomOneNeutralinoThree = "
      << gluinoSbottomOneNeutralinoThree;
      std::cout << std::endl;
      std::cout
      << "gluinoSbottomOneNeutralinoOne = { "
      << gluinoSbottomOneNeutralinoOne->getAsString() << " }";
      std::cout << std::endl;
      std::cout
      << "gluinoSbottomOneNeutralinoTwo = { "
      << gluinoSbottomOneNeutralinoTwo->getAsString() << " }";
      std::cout << std::endl;
      std::cout
      << "goSb1N3 = { "
      << gluinoSbottomOneNeutralinoThree->getAsString() << " }";
      std::cout << std::endl;

      std::cout
      << std::endl
      << "gluinoSbottomOneNeutralinoTwo, gluinoSbottomOneNeutralinoThree, 1"
      << " jet: "
      << testChargeAndFlavorSummed.withExactlyNJets( &testSignalDefinitionSet,
                                                     1,
                                                 gluinoSbottomOneNeutralinoTwo,
                                                     true,
                                               gluinoSbottomOneNeutralinoThree,
                                                     true );
      std::cout
      << std::endl
      << "gluinoSbottomOneNeutralinoThree, gluinoSbottomOneNeutralinoTwo, 1"
      << " jet: "
      << testChargeAndFlavorSummed.withExactlyNJets( &testSignalDefinitionSet,
                                                     1,
                                               gluinoSbottomOneNeutralinoThree,
                                                     true,
                                                 gluinoSbottomOneNeutralinoTwo,
                                                     true );

      for( int numberOfJets( 0 );
           4 >= numberOfJets;
           ++numberOfJets )
      {
        std::cout
        << std::endl
        << "gluinoSbottomOneNeutralinoOne, gluinoSbottomOneNeutralinoOne: "
        << testChargeAndFlavorSummed.withExactlyNJets(
                                                      &testSignalDefinitionSet,
                                                       numberOfJets,
                                                 gluinoSbottomOneNeutralinoOne,
                                                       true,
                                                 gluinoSbottomOneNeutralinoOne,
                                                       true );
        std::cout
        << std::endl
        << "gluinoSbottomOneNeutralinoOne, gluinoSbottomOneNeutralinoTwo: "
        << testChargeAndFlavorSummed.withExactlyNJets(
                                                      &testSignalDefinitionSet,
                                                       numberOfJets,
                                                 gluinoSbottomOneNeutralinoOne,
                                                       true,
                                                 gluinoSbottomOneNeutralinoTwo,
                                                       true );
        std::cout
        << std::endl
        << "gluinoSbottomOneNeutralinoOne, gluinoSbottomOneNeutralinoThree: "
        << testChargeAndFlavorSummed.withExactlyNJets(
                                                      &testSignalDefinitionSet,
                                                       numberOfJets,
                                                 gluinoSbottomOneNeutralinoOne,
                                                       true,
                                               gluinoSbottomOneNeutralinoThree,
                                                       true );
        std::cout
        << std::endl
        << "gluinoSbottomOneNeutralinoTwo, gluinoSbottomOneNeutralinoOne: "
        << testChargeAndFlavorSummed.withExactlyNJets(
                                                      &testSignalDefinitionSet,
                                                       numberOfJets,
                                                 gluinoSbottomOneNeutralinoTwo,
                                                       true,
                                                 gluinoSbottomOneNeutralinoOne,
                                                       true );
        std::cout
        << std::endl
        << "gluinoSbottomOneNeutralinoTwo, gluinoSbottomOneNeutralinoTwo: "
        << testChargeAndFlavorSummed.withExactlyNJets(
                                                      &testSignalDefinitionSet,
                                                       numberOfJets,
                                                 gluinoSbottomOneNeutralinoTwo,
                                                       true,
                                                 gluinoSbottomOneNeutralinoTwo,
                                                       true );
        std::cout
        << std::endl
        << "gluinoSbottomOneNeutralinoTwo, gluinoSbottomOneNeutralinoThree: "
        << testChargeAndFlavorSummed.withExactlyNJets(
                                                      &testSignalDefinitionSet,
                                                       numberOfJets,
                                                 gluinoSbottomOneNeutralinoTwo,
                                                       true,
                                               gluinoSbottomOneNeutralinoThree,
                                                       true );
        std::cout
        << std::endl
        << "gluinoSbottomOneNeutralinoThree, gluinoSbottomOneNeutralinoOne: "
        << testChargeAndFlavorSummed.withExactlyNJets(
                                                      &testSignalDefinitionSet,
                                                       numberOfJets,
                                               gluinoSbottomOneNeutralinoThree,
                                                       true,
                                                 gluinoSbottomOneNeutralinoOne,
                                                       true );
        std::cout
        << std::endl
        << "gluinoSbottomOneNeutralinoThree, gluinoSbottomOneNeutralinoTwo: "
        << testChargeAndFlavorSummed.withExactlyNJets(
                                                      &testSignalDefinitionSet,
                                                       numberOfJets,
                                               gluinoSbottomOneNeutralinoThree,
                                                       true,
                                                 gluinoSbottomOneNeutralinoTwo,
                                                       true );
        std::cout
        << std::endl
        << "gluinoSbottomOneNeutralinoThree, gluinoSbottomOneNeutralinoThree: "
        << testChargeAndFlavorSummed.withExactlyNJets(
                                                      &testSignalDefinitionSet,
                                                       numberOfJets,
                                               gluinoSbottomOneNeutralinoThree,
                                                       true,
                                               gluinoSbottomOneNeutralinoThree,
                                                       true );
      }

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
