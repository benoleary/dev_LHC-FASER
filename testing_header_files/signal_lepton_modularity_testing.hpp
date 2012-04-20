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
      signalClasses::leptonAcceptanceForPairFactory testFactory;

      signalClasses::leptonAcceptanceForCascadePair const*
      testFullySpecified( testFactory.getFullySpecified( 1,
                                                         0,
                                                         0,
                                                         1 ) );
      signalClasses::leptonAcceptanceForCascadePair const*
      testNoLeptonCutNorExtraJetCut( testFactory.getNoLeptonCut() );
      signalClasses::leptonAcceptanceForCascadePair const*
      testChargeAndFlavorSummed( testFactory.getChargeAndFlavorSummed( 2 ) );
      signalClasses::leptonAcceptanceForCascadePair const*
      singleLepton( testFactory.getChargeAndFlavorSummed( 1 ) );
      signalClasses::leptonAcceptanceForCascadePair const*
      testChargeSummed( testFactory.getChargeSummed( 0,
                                                     1 ) );
      signalClasses::leptonAcceptanceForCascadePair const*
      singleMuon( testFactory.getChargeSummed( 0,
                                               1 ) );
      signalClasses::leptonAcceptanceForCascadePair const*
      testFlavorSummed( testFactory.getFlavorSummed( 0,
                                                     1 ) );
      signalClasses::leptonAcceptanceForCascadePair const*
      singleNegative( testFactory.getFlavorSummed( 1,
                                                   0 ) );
      signalClasses::leptonAcceptanceForCascadePair const*
      singlePositive( testFactory.getFlavorSummed( 0,
                                                   1 ) );
      signalClasses::leptonAcceptanceForCascadePair const*
      testOssfMinusOsdf( testFactory.getOssfMinusOsdf() );

      std::cout
      << std::endl
      << "testChargeAndFlavorSummed = " << testChargeAndFlavorSummed
      << ", singleLepton = " << singleLepton
      << ", testChargeSummed = " << testChargeSummed
      << ", singleMuon = " << singleMuon
      << ", testFlavorSummed = " << testFlavorSummed
      << ", singleNegative = " << singleNegative
      << ", singlePositive = " << singlePositive;
      std::cout << std::endl;


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
            << "testFullySpecified->withExactlyNJets( ..., " << numberOfJets
            << ", ... ) = "
            << testFullySpecified->withExactlyNJets( &testSignalDefinitionSet,
                                                     numberOfJets,
                                             gluinoCascades->at( gluinoIndex ),
                                                     true,
                                             sdownLCascades->at( sdownLIndex ),
                                                     true );
            std::cout
            << std::endl
            << "testNoLeptonCutNorExtraJetCut->withExactlyNJets( ..., "
            << numberOfJets << ", ... ) = "
            << testNoLeptonCutNorExtraJetCut->withExactlyNJets(
                                                      &testSignalDefinitionSet,
                                                                numberOfJets,
                                             gluinoCascades->at( gluinoIndex ),
                                                                true,
                                             sdownLCascades->at( sdownLIndex ),
                                                                true );

            std::cout
            << std::endl
            << "testChargeAndFlavorSummed->withExactlyNJets( ..., "
            << numberOfJets << ", ... ) = "
            << testChargeAndFlavorSummed->withExactlyNJets(
                                                      &testSignalDefinitionSet,
                                                            numberOfJets,
                                             gluinoCascades->at( gluinoIndex ),
                                                            true,
                                             sdownLCascades->at( sdownLIndex ),
                                                            true );

            std::cout
            << std::endl
            << "singleLepton->withExactlyNJets( ..., "
            << numberOfJets << ", ... ) = "
            << singleLepton->withExactlyNJets( &testSignalDefinitionSet,
                                               numberOfJets,
                                             gluinoCascades->at( gluinoIndex ),
                                               true,
                                             sdownLCascades->at( sdownLIndex ),
                                               true );

            std::cout
            << std::endl
            << "testChargeSummed->withExactlyNJets( ..., " << numberOfJets
            << ", ... ) = "
            << testChargeSummed->withExactlyNJets( &testSignalDefinitionSet,
                                                   numberOfJets,
                                             gluinoCascades->at( gluinoIndex ),
                                                   true,
                                             sdownLCascades->at( sdownLIndex ),
                                                   true );

            std::cout
            << std::endl
            << "singleMuon->withExactlyNJets( ..., " << numberOfJets
            << ", ... ) = "
            << singleMuon->withExactlyNJets( &testSignalDefinitionSet,
                                             numberOfJets,
                                             gluinoCascades->at( gluinoIndex ),
                                             true,
                                             sdownLCascades->at( sdownLIndex ),
                                             true );

            std::cout
            << std::endl
            << "testFlavorSummed->withExactlyNJets( ..., " << numberOfJets
            << ", ... ) = "
            << testFlavorSummed->withExactlyNJets( &testSignalDefinitionSet,
                                                   numberOfJets,
                                             gluinoCascades->at( gluinoIndex ),
                                                   true,
                                             sdownLCascades->at( sdownLIndex ),
                                                   true );

            std::cout
            << std::endl
            << "singleNegative->withExactlyNJets( ..., " << numberOfJets
            << ", ... ) = "
            << singleNegative->withExactlyNJets( &testSignalDefinitionSet,
                                                 numberOfJets,
                                             gluinoCascades->at( gluinoIndex ),
                                                 true,
                                             sdownLCascades->at( sdownLIndex ),
                                                 true );

            std::cout
            << std::endl
            << "singlePositive->withExactlyNJets( ..., " << numberOfJets
            << ", ... ) = "
            << singlePositive->withExactlyNJets( &testSignalDefinitionSet,
                                                 numberOfJets,
                                             gluinoCascades->at( gluinoIndex ),
                                                 true,
                                             sdownLCascades->at( sdownLIndex ),
                                                 true );

            std::cout
            << std::endl
            << "testOssfMinusOsdf->withExactlyNJets( ..., " << numberOfJets
            << ", ... ) = "
            << testOssfMinusOsdf->withExactlyNJets( &testSignalDefinitionSet,
                                                    numberOfJets,
                                             gluinoCascades->at( gluinoIndex ),
                                                    true,
                                             sdownLCascades->at( sdownLIndex ),
                                                    true );
          }
        }
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
