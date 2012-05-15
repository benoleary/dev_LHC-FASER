/*
 * LHC-FASER_signal_calculator_stuff.cpp
 *
 *  Created on: 17 Nov 2010
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

#include "LHC-FASER_signal_calculator_stuff.hpp"

namespace LHC_FASER
{
  int const leptonAcceptanceForCascadePair::maximumNumberOfJets( 2
                                     * fullCascade::maximumJetsFromEwCascade );

  leptonAcceptanceForCascadePair::leptonAcceptanceForCascadePair()
  {
    // does nothing.
  }

  leptonAcceptanceForCascadePair::~leptonAcceptanceForCascadePair()
  {
    // does nothing.
  }


  namespace leptonAcceptanceStyle
  {
    fullySpecified::fullySpecified( int const numberOfNegativeElectrons,
                                    int const numberOfPositiveElectrons,
                                    int const numberOfNegativeMuons,
                                    int const numberOfPositiveMuons ) :
      leptonAcceptanceForCascadePair(),
      numberOfNegativeElectrons( numberOfNegativeElectrons ),
      numberOfPositiveElectrons( numberOfPositiveElectrons ),
      numberOfNegativeMuons( numberOfNegativeMuons ),
      numberOfPositiveMuons( numberOfPositiveMuons )
    {
      // just an initialization list.
    }

    fullySpecified::~fullySpecified()
    {
      // does nothing.
    }


    double
    fullySpecified::withExactlyNJets(
                            signalDefinitionSet const* const signalDefinitions,
                                      int const exactNumberOfJets,
                                      fullCascade* const firstCascade,
                                      bool const firstIsNotAntiparticle,
                                      fullCascade* const secondCascade,
                                     bool const secondIsNotAntiparticle ) const
    {
      double returnValue( 0.0 );

      for( int jetsFromFirst( exactNumberOfJets );
           0 <= jetsFromFirst;
           --jetsFromFirst )
      {
        for( int negativeElectronsFromFirst( numberOfNegativeElectrons );
             0 <= negativeElectronsFromFirst;
             --negativeElectronsFromFirst )
        {
          for( int positiveElectronsFromFirst( numberOfPositiveElectrons );
               0 <= positiveElectronsFromFirst;
               --positiveElectronsFromFirst )
          {
            for( int negativeMuonsFromFirst( numberOfNegativeMuons );
                 0 <= negativeMuonsFromFirst;
                 --negativeMuonsFromFirst )
            {
              for( int positiveMuonsFromFirst( numberOfPositiveMuons );
                   0 <= positiveMuonsFromFirst;
                   --positiveMuonsFromFirst )
              {
                returnValue
                += ( firstCascade->getAcceptance( firstIsNotAntiparticle,
                                                  signalDefinitions,
                                                  jetsFromFirst,
                                                  negativeElectronsFromFirst,
                                                  positiveElectronsFromFirst,
                                                  negativeMuonsFromFirst,
                                                  positiveMuonsFromFirst )
                     * secondCascade->getAcceptance( secondIsNotAntiparticle,
                                                     signalDefinitions,
                                         ( exactNumberOfJets - jetsFromFirst ),
                    ( numberOfNegativeElectrons - negativeElectronsFromFirst ),
                    ( numberOfPositiveElectrons - positiveElectronsFromFirst ),
                            ( numberOfNegativeMuons - negativeMuonsFromFirst ),
                        ( numberOfPositiveMuons - positiveMuonsFromFirst ) ) );
              }
            }
          }
        }
      }
      return returnValue;
    }


    noLeptonCutNorExtraJetCut::noLeptonCutNorExtraJetCut() :
        leptonAcceptanceForCascadePair()
    {
      // just an initialization list.
    }

    noLeptonCutNorExtraJetCut::~noLeptonCutNorExtraJetCut()
    {
      // does nothing.
    }



    chargeAndFlavorSummed::chargeAndFlavorSummed( int const numberOfLeptons ) :
        leptonAcceptanceForCascadePair(),
        numberOfLeptons( numberOfLeptons )
    {
      // just an initialization list.
    }

    chargeAndFlavorSummed::~chargeAndFlavorSummed()
    {
      // does nothing.
    }


    double
    chargeAndFlavorSummed::withExactlyNJets(
                            signalDefinitionSet const* const signalDefinitions,
                                             int const exactNumberOfJets,
                                             fullCascade* const firstCascade,
                                             bool const firstIsNotAntiparticle,
                                             fullCascade* const secondCascade,
                                     bool const secondIsNotAntiparticle ) const
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "chargeAndFlavorSummed::withExactlyNJets( " << signalDefinitions
      << ", " << exactNumberOfJets << ", " << firstCascade << ", "
      << firstIsNotAntiparticle << ", " << secondCascade << ", "
      << secondIsNotAntiparticle << " ) called.";
      std::cout << std::endl;
      std::cout
      << "firstCascade = { "
      << firstCascade->getAsString() << " }, secondCascade = { "
      << secondCascade->getAsString() << " }";
      std::cout << std::endl;**/

      double returnValue( 0.0 );

      for( int jetsFromFirst( exactNumberOfJets );
           0 <= jetsFromFirst;
           --jetsFromFirst )
      {
        for( int leptonsFromFirst( numberOfLeptons );
             0 <= leptonsFromFirst;
             --leptonsFromFirst )
        {
          // debugging:
          /**std::cout << std::endl << "debugging:"
          << std::endl
          << "jetsFromFirst = " << jetsFromFirst
          << ", ( exactNumberOfJets - jetsFromFirst ) = "
          << ( exactNumberOfJets - jetsFromFirst )
          << ", leptonsFromFirst = " << leptonsFromFirst
          << ", ( numberOfLeptons - leptonsFromFirst ) = "
          << ( numberOfLeptons - leptonsFromFirst )
          << ", adding "
          << ( firstCascade->specifiedJetsSpecifiedChargeSummedLeptons(
                                                             signalDefinitions,
                                                             jetsFromFirst,
                                                             leptonsFromFirst )
               * secondCascade->specifiedJetsSpecifiedChargeSummedLeptons(
                                                             signalDefinitions,
                                         ( exactNumberOfJets - jetsFromFirst ),
                                    ( numberOfLeptons - leptonsFromFirst ) ) );
          std::cout << std::endl;**/

          returnValue
          += ( firstCascade->specifiedJetsSpecifiedChargeSummedLeptons(
                                                             signalDefinitions,
                                                                 jetsFromFirst,
                                                             leptonsFromFirst )
               * secondCascade->specifiedJetsSpecifiedChargeSummedLeptons(
                                                             signalDefinitions,
                                         ( exactNumberOfJets - jetsFromFirst ),
                                    ( numberOfLeptons - leptonsFromFirst ) ) );
        }
      }
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "chargeAndFlavorSummed::withExactlyNJets( " << signalDefinitions
      << ", " << exactNumberOfJets << ", { " << firstCascade->getAsString()
      << " }, " << firstIsNotAntiparticle << ", { "
      << secondCascade->getAsString() << " }, " << secondIsNotAntiparticle
      << " ) returning " << returnValue;
      std::cout << std::endl;**/

      return returnValue;
    }



    chargeSummed::chargeSummed( int const numberOfElectrons,
                                int const numberOfMuons ) :
        leptonAcceptanceForCascadePair(),
        numberOfElectrons( numberOfElectrons ),
        numberOfMuons( numberOfMuons )
    {
      // just an initialization list.
    }

    chargeSummed::~chargeSummed()
    {
      // does nothing.
    }


    double
    chargeSummed::withExactlyNJets(
                            signalDefinitionSet const* const signalDefinitions,
                                    int const exactNumberOfJets,
                                    fullCascade* const firstCascade,
                                    bool const firstIsNotAntiparticle,
                                    fullCascade* const secondCascade,
                                    bool const secondIsNotAntiparticle ) const
    {
      double returnValue( 0.0 );

      for( int jetsFromFirst( exactNumberOfJets );
           0 <= jetsFromFirst;
           --jetsFromFirst )
      {
        for( int electronsFromFirst( numberOfElectrons );
             0 <= electronsFromFirst;
             --electronsFromFirst )
        {
          for( int muonsFromFirst( numberOfMuons );
               0 <= muonsFromFirst;
               --muonsFromFirst )
          {
            returnValue
            += ( firstCascade->leptonChargeSummedWithSpecifiedJets(
                                                        firstIsNotAntiparticle,
                                                             signalDefinitions,
                                                                 jetsFromFirst,
                                                            electronsFromFirst,
                                                            muonsFromFirst )
               * secondCascade->leptonChargeSummedWithSpecifiedJets(
                                                       secondIsNotAntiparticle,
                                                             signalDefinitions,
                                         ( exactNumberOfJets - jetsFromFirst ),
                                    ( numberOfElectrons - electronsFromFirst ),
                                        ( numberOfMuons - muonsFromFirst ) ) );
          }
        }
      }
      return returnValue;
    }



    flavorSummed::flavorSummed( int const numberOfNegativeLeptons,
                                int const numberOfPositiveLeptons ) :
        leptonAcceptanceForCascadePair(),
        numberOfNegativeLeptons( numberOfNegativeLeptons ),
        numberOfPositiveLeptons( numberOfPositiveLeptons )
    {
      // just an initialization list.
    }

    flavorSummed::~flavorSummed()
    {
      // does nothing.
    }


    double
    flavorSummed::withExactlyNJets(
                            signalDefinitionSet const* const signalDefinitions,
                                    int const exactNumberOfJets,
                                    fullCascade* const firstCascade,
                                    bool const firstIsNotAntiparticle,
                                    fullCascade* const secondCascade,
                                    bool const secondIsNotAntiparticle ) const
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "flavorSummed::withExactlyNJets( " << signalDefinitions
      << ", " << exactNumberOfJets << ", " << firstCascade << ", "
      << firstIsNotAntiparticle << ", " << secondCascade << ", "
      << secondIsNotAntiparticle << " ) called.";
      std::cout << std::endl;
      std::cout
      << "firstCascade = { "
      << firstCascade->getAsString() << " }, secondCascade = { "
      << secondCascade->getAsString() << " }";
      std::cout << std::endl;**/

      double returnValue( 0.0 );

      for( int jetsFromFirst( exactNumberOfJets );
           0 <= jetsFromFirst;
           --jetsFromFirst )
      {
        for( int negativeLeptonsFromFirst( numberOfNegativeLeptons );
             0 <= negativeLeptonsFromFirst;
             --negativeLeptonsFromFirst )
        {
          for( int positiveLeptonsFromFirst( numberOfPositiveLeptons );
               0 <= positiveLeptonsFromFirst;
               --positiveLeptonsFromFirst )
          {
            // debugging:
            /**std::cout << std::endl << "debugging:"
            << std::endl
            << "jetsFromFirst = " << jetsFromFirst
            << ", ( exactNumberOfJets - jetsFromFirst ) = "
            << ( exactNumberOfJets - jetsFromFirst )
            << ", negativeLeptonsFromFirst = " << negativeLeptonsFromFirst
            << ", positiveLeptonsFromFirst = " << positiveLeptonsFromFirst
            << ", ( numberOfNegativeLeptons - negativeLeptonsFromFirst ) = "
            << ( numberOfNegativeLeptons - negativeLeptonsFromFirst )
            << ", ( numberOfPositiveLeptons - positiveLeptonsFromFirst ) = "
            << ( numberOfPositiveLeptons - positiveLeptonsFromFirst )
            << ", adding "
            << ( firstCascade->leptonFlavorSummedWithSpecifiedJets(
                                                        firstIsNotAntiparticle,
                                                             signalDefinitions,
                                                                 jetsFromFirst,
                                                      negativeLeptonsFromFirst,
                                                     positiveLeptonsFromFirst )
                 * secondCascade->leptonFlavorSummedWithSpecifiedJets(
                                                       secondIsNotAntiparticle,
                                                             signalDefinitions,
                                         ( exactNumberOfJets - jetsFromFirst ),
                        ( numberOfNegativeLeptons - negativeLeptonsFromFirst ),
                    ( numberOfPositiveLeptons - positiveLeptonsFromFirst ) ) );
            std::cout << std::endl;**/

            returnValue
            += ( firstCascade->leptonFlavorSummedWithSpecifiedJets(
                                                        firstIsNotAntiparticle,
                                                             signalDefinitions,
                                                                 jetsFromFirst,
                                                      negativeLeptonsFromFirst,
                                                     positiveLeptonsFromFirst )
                 * secondCascade->leptonFlavorSummedWithSpecifiedJets(
                                                       secondIsNotAntiparticle,
                                                             signalDefinitions,
                                         ( exactNumberOfJets - jetsFromFirst ),
                        ( numberOfNegativeLeptons - negativeLeptonsFromFirst ),
                    ( numberOfPositiveLeptons - positiveLeptonsFromFirst ) ) );
          }
        }
      }
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "flavorSummed::withExactlyNJets( " << signalDefinitions
      << ", " << exactNumberOfJets << ", { " << firstCascade->getAsString()
      << " }, " << firstIsNotAntiparticle << ", { "
      << secondCascade->getAsString() << " }, " << secondIsNotAntiparticle
      << " ) returning " << returnValue;
      std::cout << std::endl;**/

      return returnValue;
    }



    ossfMinusOsdf::ossfMinusOsdf() :
        leptonAcceptanceForCascadePair()
    {
      // just an initialization list.
    }

    ossfMinusOsdf::~ossfMinusOsdf()
    {
      // does nothing.
    }


    double
    ossfMinusOsdf::withExactlyNJets(
                            signalDefinitionSet const* const signalDefinitions,
                                     int const exactNumberOfJets,
                                     fullCascade* const firstCascade,
                                     bool const firstIsNotAntiparticle,
                                     fullCascade* const secondCascade,
                                     bool const secondIsNotAntiparticle ) const
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "ossfMinusOsdf::withExactlyNJets( " << signalDefinitions
      << ", " << exactNumberOfJets << ", " << firstCascade << ", "
      << firstIsNotAntiparticle << ", " << secondCascade << ", "
      << secondIsNotAntiparticle << " ) called.";
      std::cout << std::endl;
      std::cout
      << "firstCascade = { "
      << firstCascade->getAsString() << " }, secondCascade = { "
      << secondCascade->getAsString() << " }";
      std::cout << std::endl;**/

      double returnValue( 0.0 );

      for( int jetsFromFirst( exactNumberOfJets );
           0 <= jetsFromFirst;
           --jetsFromFirst )
      {
        // debugging:
        /**std::cout << std::endl << "debugging:"
        << std::endl
        << "firstCascade->specifiedJetsOneOssfMinusOsdfPair( "
        << signalDefinitions << ", " << jetsFromFirst << " ) = "
        << firstCascade->specifiedJetsOneOssfMinusOsdfPair( signalDefinitions,
                                                            jetsFromFirst )
        << ", secondCascade->getAcceptance( " << secondIsNotAntiparticle
        << ", " << signalDefinitions << ", "
        << ( exactNumberOfJets - jetsFromFirst ) << ", 0, 0, 0, 0 ) = "
        << secondCascade->getAcceptance( secondIsNotAntiparticle,
                                         signalDefinitions,
                                         ( exactNumberOfJets - jetsFromFirst ),
                                         0,
                                         0,
                                         0,
                                         0 )
        << ", firstCascade->getAcceptance( " << firstIsNotAntiparticle
        << ", " << signalDefinitions << ", "
        << jetsFromFirst << ", 0, 0, 0, 0 ) = "
        << firstCascade->getAcceptance( firstIsNotAntiparticle,
                                         signalDefinitions,
                                         jetsFromFirst,
                                         0,
                                         0,
                                         0,
                                         0 )
        << ", secondCascade->specifiedJetsOneOssfMinusOsdfPair( "
        << signalDefinitions << ", " << ( exactNumberOfJets - jetsFromFirst )
        << " ) = "
        << secondCascade->specifiedJetsOneOssfMinusOsdfPair( signalDefinitions,
                                       ( exactNumberOfJets - jetsFromFirst ) );
        std::cout << std::endl;**/

        returnValue
        += ( ( firstCascade->specifiedJetsOneOssfMinusOsdfPair(
                                                             signalDefinitions,
                                                                jetsFromFirst )
               * secondCascade->getAcceptance( secondIsNotAntiparticle,
                                               signalDefinitions,
                                         ( exactNumberOfJets - jetsFromFirst ),
                                               0,
                                               0,
                                               0,
                                               0 ) )
             + ( firstCascade->getAcceptance( firstIsNotAntiparticle,
                                              signalDefinitions,
                                              jetsFromFirst,
                                              0,
                                              0,
                                              0,
                                              0 )
                 * secondCascade->specifiedJetsOneOssfMinusOsdfPair(
                                                             signalDefinitions,
                                   ( exactNumberOfJets - jetsFromFirst ) ) ) );
      }
      return returnValue;
    }



    sameSignSameFlavor::sameSignSameFlavor() :
        leptonAcceptanceForCascadePair()
    {
      // just an initialization list.
    }

    sameSignSameFlavor::~sameSignSameFlavor()
    {
      // does nothing.
    }


    double
    sameSignSameFlavor::withExactlyNJets(
                            signalDefinitionSet const* const signalDefinitions,
                                          int const exactNumberOfJets,
                                          fullCascade* const firstCascade,
                                          bool const firstIsNotAntiparticle,
                                          fullCascade* const secondCascade,
                                     bool const secondIsNotAntiparticle ) const
    {
      double returnValue( 0.0 );

      for( int jetsFromFirst( exactNumberOfJets );
           0 <= jetsFromFirst;
           --jetsFromFirst )
      {
        for( int leptonsFromFirst( 2 );
             0 <= leptonsFromFirst;
             --leptonsFromFirst )
        {
          returnValue
          += ( ( firstCascade->getAcceptance( firstIsNotAntiparticle,
                                              signalDefinitions,
                                              jetsFromFirst,
                                              leptonsFromFirst,
                                              0,
                                              0,
                                              0 )
                 * secondCascade->getAcceptance( secondIsNotAntiparticle,
                                                 signalDefinitions,
                                         ( exactNumberOfJets - jetsFromFirst ),
                                                 ( 2 - leptonsFromFirst ),
                                                 0,
                                                 0,
                                                 0 ) )
                 + ( firstCascade->getAcceptance( firstIsNotAntiparticle,
                                                  signalDefinitions,
                                                  jetsFromFirst,
                                                  0,
                                                  leptonsFromFirst,
                                                  0,
                                                  0 )
                     * secondCascade->getAcceptance( secondIsNotAntiparticle,
                                                     signalDefinitions,
                                         ( exactNumberOfJets - jetsFromFirst ),
                                                     0,
                                                     ( 2 - leptonsFromFirst ),
                                                     0,
                                                     0 ) )
                 + ( firstCascade->getAcceptance( firstIsNotAntiparticle,
                                                  signalDefinitions,
                                                  jetsFromFirst,
                                                  0,
                                                  0,
                                                  leptonsFromFirst,
                                                  0 )
                     * secondCascade->getAcceptance( secondIsNotAntiparticle,
                                                     signalDefinitions,
                                         ( exactNumberOfJets - jetsFromFirst ),
                                                     0,
                                                     0,
                                                     ( 2 - leptonsFromFirst ),
                                                     0 ) )
                + ( firstCascade->getAcceptance( firstIsNotAntiparticle,
                                                 signalDefinitions,
                                                 jetsFromFirst,
                                                 0,
                                                 0,
                                                 0,
                                                 leptonsFromFirst )
                   * secondCascade->getAcceptance( secondIsNotAntiparticle,
                                                   signalDefinitions,
                                         ( exactNumberOfJets - jetsFromFirst ),
                                                   0,
                                                   0,
                                                   0,
                                                ( 2 - leptonsFromFirst ) ) ) );
        }
      }
      return returnValue;
    }

  }  // end of leptonAcceptanceStyle namespace



  signalCalculator::signalCalculator(
                         signalDefinitionSet const* const signalDefinitions ) :
    signalName( "noNameYet" ),
    signalDefinitions( signalDefinitions ),
    inputShortcut( signalDefinitions->getShortcuts()->getInputShortcuts() ),
    productionChannels(),
    firstCascades( NULL ),
    secondCascades( NULL ),
    excludedFinalStateParticles(),
    firstCascadeBrToEwino( CppSLHA::CppSLHA_global::really_wrong_value ),
    secondCascadeBrToEwino( CppSLHA::CppSLHA_global::really_wrong_value ),
    firstSparticleIsNotAntiparticle( true ),
    secondSparticleIsNotAntiparticle( true )
  {
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "signalCalculator::signalCalculator( " << signalDefinitions
    << " ) called.";
    std::cout << std::endl;
    std::cout
    << "signalDefinitions->getShortcuts()->getInputShortcuts("
    << ")->getScoloredProductionCombinations() = "
    << signalDefinitions->getShortcuts()->getInputShortcuts(
                                        )->getScoloredProductionCombinations();
    std::cout << std::endl;**/

    for( std::vector< signedParticleShortcutPair* >::const_iterator
         pairIterator( signalDefinitions->getShortcuts()->getInputShortcuts(
                             )->getScoloredProductionCombinations()->begin() );
         signalDefinitions->getShortcuts()->getInputShortcuts(
                                  )->getScoloredProductionCombinations()->end()
         > pairIterator;
         ++pairIterator )
    {
      // debugging:
      /**std::cout << std::endl << "debugging:" << std::endl;
      std::cout
      << "signalDefinitions->getCrossSections() = "
      << signalDefinitions->getCrossSections();
      std::cout << std::endl;
      std::cout
      << "signalDefinitions->getCrossSections()->getTable( *pairIterator ) = "
      << signalDefinitions->getCrossSections()->getTable( *pairIterator );
      std::cout << std::endl;
      std::cout
      << "signalDefinitions->getShortcuts() = "
      << signalDefinitions->getShortcuts();
      std::cout << std::endl;**/

      productionChannels.push_back( new productionChannelPointerSet(
                                                             signalDefinitions,
                                                             *pairIterator ) );
    }
  }

  // the version without any argument is used for reallyWrongCalculator:
  signalCalculator::signalCalculator() :
    signalName( "noNameYet" ),
    signalDefinitions( NULL ),
    inputShortcut( NULL ),
    productionChannels(),
    firstCascades( NULL ),
    secondCascades( NULL ),
    excludedFinalStateParticles(),
    firstCascadeBrToEwino( CppSLHA::CppSLHA_global::really_wrong_value ),
    secondCascadeBrToEwino( CppSLHA::CppSLHA_global::really_wrong_value ),
    firstSparticleIsNotAntiparticle( true ),
    secondSparticleIsNotAntiparticle( true )
  {
    // just an initialization list.
  }

  signalCalculator::~signalCalculator()
  {
    for( std::vector< productionChannelPointerSet* >::iterator
         deletionIterator( productionChannels.begin() );
         productionChannels.end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
  }


  bool
  signalCalculator::goThroughCascadesNormally( double* signalValue,
                                               double* uncertaintyFactor )
  // this calls valueForCurrentCascades(...) for each pair of cascades for
  // each production channel.
  {
    // debugging:
    /**std::cout
    << std::endl
    << "debugging: signalCalculator::goThroughCascadesNormally( "
    << signalValue << ", " << uncertaintyFactor << " ) called."
    << " signalDefinitions->getPrimaryLeptonCut() = "
    << signalDefinitions->getPrimaryLeptonCut();
    std::cout << std::endl;**/

    // start by setting the signal value to 0.0:
    *signalValue = 0.0;

    for( std::vector< productionChannelPointerSet* >::iterator
         channelIterator( productionChannels.begin() );
         productionChannels.end() > channelIterator;
         ++channelIterator )
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << *((*channelIterator)->getScoloredPair()->getFirstParticle(
                )->get_name_or_antiname( (*channelIterator)->getScoloredPair(
                                               )->firstIsNotAntiparticle() ))
      << " + "
      << *((*channelIterator)->getScoloredPair()->getSecondParticle(
                )->get_name_or_antiname( (*channelIterator)->getScoloredPair(
                                             )->secondIsNotAntiparticle() ))
      << "; cross-section = " << (*channelIterator)->getCrossSection();
      std::cout << std::endl;**/

      if( lhcFaserGlobal::negligibleSigma
          < (*channelIterator)->getCrossSection() )
        // if it's worth looking at this channel...
      {
        // we have to look at all the open cascade pairings:
        firstCascades
        = (*channelIterator)->getFirstCascadeSet()->getOpenCascades();
        secondCascades
        = (*channelIterator)->getSecondCascadeSet()->getOpenCascades();
        initialPair = (*channelIterator)->getScoloredPair();
        firstSparticleIsNotAntiparticle
        = initialPair->firstIsNotAntiparticle();
        secondSparticleIsNotAntiparticle
        = initialPair->secondIsNotAntiparticle();

        for( std::vector< fullCascade* >::iterator
             firstCascadeIterator( firstCascades->begin() );
             firstCascades->end() > firstCascadeIterator;
             ++firstCascadeIterator )
        {
          // debugging:
          /**std::cout << std::endl << "debugging:"
          << std::endl
          << "1st: " << (*firstCascadeIterator)->getAsString()
          << ", firstCascades->size() = " << firstCascades->size();
          std::cout << std::endl;**/

          firstCascadeBrToEwino = (*firstCascadeIterator)->getBrToEwino(
                                                &excludedFinalStateParticles );
          if( ( lhcFaserGlobal::negligibleBr
                * 2.0 * (double)(firstCascades->size()) )
              < firstCascadeBrToEwino )
            // this should mean we never throw away more than
            // 0.5 * lhcFaserGlobal::negligibleBr of acceptance.
          {
            for( std::vector< fullCascade* >::iterator
                 secondCascadeIterator( secondCascades->begin() );
                 secondCascades->end() > secondCascadeIterator;
                 ++secondCascadeIterator )
            {
              // debugging:
              /**std::cout << std::endl << "debugging:"
              << std::endl
              << "2nd: " << (*secondCascadeIterator)->getAsString();
              std::cout << std::endl;**/

              secondCascadeBrToEwino
              = (*secondCascadeIterator)->getBrToEwino(
                                                &excludedFinalStateParticles );
              if( ( lhcFaserGlobal::negligibleBr
                    * 2.0 * (double)(secondCascades->size()) )
                  < secondCascadeBrToEwino )
                // this should mean we never throw away more than
                // 0.5 * lhcFaserGlobal::negligibleBr of acceptance.
              {
                subchannelCrossSectionTimesBrToEwinos
                = ( (*channelIterator)->getCrossSection()
                    * firstCascadeBrToEwino * secondCascadeBrToEwino );

                // debugging:
                /**std::cout << std::endl << "debugging:"
                << std::endl
                << "(*channelIterator)->getCrossSection() = "
                << (*channelIterator)->getCrossSection()
                << ", firstCascadeBrToEwino = " << firstCascadeBrToEwino
                << ", secondCascadeBrToEwino = " << secondCascadeBrToEwino
                << ", subchannelCrossSectionTimesBrToEwinos = "
                << subchannelCrossSectionTimesBrToEwinos;
                std::cout << std::endl;**/

                *signalValue += valueForCurrentCascades( *firstCascadeIterator,
                                                      *secondCascadeIterator );
              }
              // end of if the BR of the 2nd cascade to its electroweakino was
              // not negligible.
            }  // end of loop over 2nd cascade.
          }
          // end of if the BR of the 1st cascade to its electroweakino was not
          // negligible.
        }  // end of loop over 1st cascade.
      }  // end of if the channel's cross-section was not negligible.
    }  // end of loop over productionChannels.
    // currently, we do the uncertainty very roughly:
    *uncertaintyFactor
    = ( *signalValue * signalDefinitions->getShortcuts()->getUncertainty() );
    return true;
  }

  bool
  signalCalculator::parseBeamEnergy( std::string& argumentString,
                                 signalDefinitionSet* const signalDefinitions )
  /* this looks for "_7TeV", "_07TeV", "_10TeV", or "_14TeV", removes it from
   * argumentString & sets the beam energy appropriately. returns false if it
   * could not find any (& does not modify argumentString in this case).
   */
  {
    if( findAndRemoveSubstring(  "_7TeV",
                                 argumentString ) )
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "signalCalculator::parseBeamEnergy setting to 7 TeV.";
      std::cout << std::endl;**/

      signalDefinitions->setBeamEnergy( 7 );
      return true;
    }
    else if( findAndRemoveSubstring(  "_07TeV",
                                      argumentString ) )
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "signalCalculator::parseBeamEnergy setting to 07 TeV.";
      std::cout << std::endl;**/

      signalDefinitions->setBeamEnergy( 7 );
      return true;
    }
    else if( findAndRemoveSubstring(  "_10TeV",
                                      argumentString ) )
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "signalCalculator::parseBeamEnergy setting to 10 TeV.";
      std::cout << std::endl;**/

      signalDefinitions->setBeamEnergy( 10 );
      return true;
    }
    else if( findAndRemoveSubstring(  "_14TeV",
                                      argumentString ) )
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "signalCalculator::parseBeamEnergy setting to 14 TeV.";
      std::cout << std::endl;**/

      signalDefinitions->setBeamEnergy( 14 );
      return true;
    }
    else
    {
      return false;
    }
  }

  leptonAcceptanceForCascadePair*
  signalCalculator::parseLeptonAcceptance( std::string& argumentString,
                                 signalDefinitionSet* const signalDefinitions )
  /* this looks for strings encoding the type of lepton cuts to use. the
   * strings are, where # stands for any string representing an integer, in the
   * order in which they are checked:
   * "_ossf" : ossfMinusOsdf
   * "_sssf" : sameSignSameFlavor
   * "_noExtraCut" : noLeptonCutNorExtraJetCut
   * "_#lm#lp" : flavorSummed
   * "_#l" : chargeAndFlavorSummed
   * "_#epm#mpm" : chargeSummed
   * "_#em#ep#mm#mp" : fullySpecified
   * parseLeptonTransverseMomentumCuts is then called on the remainder of
   * argumentString. finally, it creates a new leptonAcceptanceForCascadePair
   * & returns a pointer to it. NULL is returned if argumentString could not
   * be properly interpretted.
   */
  {
    leptonAcceptanceForCascadePair* returnPointer( NULL );

    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "signalCalculator::parseLeptonAcceptance( \"" << argumentString
    << "\", " << signalDefinitions << " ) called.";
    std::cout << std::endl;**/


    // 1st look for the easy cases:
    if( findAndRemoveSubstring( "_ossf",
                                argumentString ) )
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "signalCalculator::parseLeptonAcceptance found \"_ossf\".";
      std::cout << std::endl;**/
      returnPointer = new leptonAcceptanceStyle::ossfMinusOsdf();
    }
    else if( findAndRemoveSubstring( "_sssf",
                                     argumentString ) )
    {
      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "signalCalculator::parseLeptonAcceptance found \"_sssf\". remember"
      << " that the same-sign dilepton signal is probably flavor-summed; gotta"
      << " check this.";
      std::cout << std::endl;/**/
      returnPointer = new leptonAcceptanceStyle::sameSignSameFlavor();
      // the default behavior of sameSignSameFlavor is to set the secondary
      // lepton cut to be the same as the primary lepton cut:
      signalDefinitions->setSecondaryLeptonCut(
                                    signalDefinitions->getPrimaryLeptonCut() );
      // if the cuts need to be different, they can be over-written by
      // parseLeptonTransverseMomentumCuts later.
    }
    else if( findAndRemoveSubstring( "_noExtraCut",
                                     argumentString ) )
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "signalCalculator::parseLeptonAcceptance found \"_noExtraCut\".";
      std::cout << std::endl;**/
      returnPointer = new leptonAcceptanceStyle::noLeptonCutNorExtraJetCut();
    }
    else if( ( 0 == argumentString.compare( 0,
                                            1,
                                            "_" ) ) )
    {
      argumentString.erase( argumentString.begin() );

      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "argumentString now = \"" << argumentString
      << "\", about to look for an int...";
      std::cout << std::endl;**/

      int leptonNumber;
      if( separateOutIntegerSubstring( leptonNumber,
                                       argumentString ) )
        // otherwise look for a string of numeric characters...
      {

        // debugging:
        /**std::cout << std::endl << "debugging:"
        << std::endl
        << "argumentString now = \"" << argumentString
        << "\", int is " << leptonNumber;
        std::cout << std::endl;**/

        // ... followed by...
        if( 0 == argumentString.compare( 0,
                                         2,
                                         "lm" ) )
        {
          argumentString.erase( argumentString.begin(),
                                ( argumentString.begin() + 2 ) );
          int secondLeptonNumber;
          if( separateOutIntegerSubstring( secondLeptonNumber,
                                           argumentString )
              &&
              ( 0 == argumentString.compare( 0,
                                             2,
                                             "lp" ) ) )
          {
            // debugging:
            /**/std::cout << std::endl << "debugging:"
            << std::endl
            << "signalCalculator::parseLeptonAcceptance found \"_"
            << leptonNumber << "lm"
            << secondLeptonNumber
            << "lp\". remember that the same-sign dilepton signal is probably"
            << " flavor-summed; gotta check this.";
            std::cout << std::endl;/**/
            argumentString.erase( argumentString.begin(),
                                  ( argumentString.begin() + 2 ) );
            returnPointer
            = new leptonAcceptanceStyle::flavorSummed( leptonNumber,
                                                       secondLeptonNumber );
          }
        }
        else if( 0 == argumentString.compare( 0,
                                              1,
                                              "l" ) )
          // if leptonNumberSubstring ends with "_#l", it's
          // chargeAndFlavorSummed.
        {
          // debugging:
          /**std::cout << std::endl << "debugging:"
          << std::endl
          << "signalCalculator::parseLeptonAcceptance found \"_"
          << leptonNumber << "l\"";
          std::cout << std::endl;**/
          argumentString.erase( argumentString.begin() );
          returnPointer
          = new leptonAcceptanceStyle::chargeAndFlavorSummed( leptonNumber );

          // debugging:
          /**std::cout << std::endl << "debugging:"
          << std::endl
          << "returnPointer = " << returnPointer;
          std::cout << std::endl;**/
        }
        else if( 0 == argumentString.compare( 0,
                                              3,
                                              "epm" ) )
        {
          argumentString.erase( argumentString.begin(),
                                ( argumentString.begin() + 3 ) );
          int secondLeptonNumber;
          if( separateOutIntegerSubstring( secondLeptonNumber,
                                           argumentString )
              &&
              ( 0 == argumentString.compare( 0,
                                             3,
                                             "mpm" ) ) )
          {
            // debugging:
            /**std::cout << std::endl << "debugging:"
            << std::endl
            << "signalCalculator::parseLeptonAcceptance found \"_"
            << leptonNumber << "epm"
            << secondLeptonNumber << "mpm\"...";
            std::cout << std::endl;**/
            argumentString.erase( argumentString.begin(),
                                  ( argumentString.begin() + 3 ) );

            returnPointer
            = new leptonAcceptanceStyle::chargeSummed( leptonNumber,
                                                       secondLeptonNumber );
          }
        }
        else if( 0 == argumentString.compare( 0,
                                              2,
                                              "em" ) )
        {
          argumentString.erase( argumentString.begin(),
                                ( argumentString.begin() + 2 ) );
          int secondLeptonNumber;
          if( separateOutIntegerSubstring( secondLeptonNumber,
                                           argumentString )
              &&
              ( 0 == argumentString.compare( 0,
                                             2,
                                             "ep" ) ) )
          {
            argumentString.erase( argumentString.begin(),
                                  ( argumentString.begin() + 2 ) );
            int thirdLeptonNumber;
            if( separateOutIntegerSubstring( thirdLeptonNumber,
                                             argumentString )
                &&
                ( 0 == argumentString.compare( 0,
                                               2,
                                               "mm" ) ) )
            {
              argumentString.erase( argumentString.begin(),
                                    ( argumentString.begin() + 2 ) );
              int fourthLeptonNumber;
              if( separateOutIntegerSubstring( fourthLeptonNumber,
                                               argumentString )
                  &&
                  ( 0 == argumentString.compare( 0,
                                                 2,
                                                 "mp" ) ) )
              {
                // debugging:
                /**std::cout << std::endl << "debugging:"
                << std::endl
                << "signalCalculator::parseLeptonAcceptance found \"_"
                << leptonNumber << "em"
                << secondLeptonNumber << "ep"
                << thirdLeptonNumber << "mm"
                << fourthLeptonNumber << "mp\"...";
                std::cout << std::endl;**/

                argumentString.erase( argumentString.begin(),
                                      ( argumentString.begin() + 2 ) );
                returnPointer
                = new leptonAcceptanceStyle::fullySpecified( leptonNumber,
                                                            secondLeptonNumber,
                                                             thirdLeptonNumber,
                                                          fourthLeptonNumber );
              }
            }
          }
        }
      }
    }
    return returnPointer;
  }

  bool
  signalCalculator::separateOutIntegerSubstring( int& destinationInt,
                                                 std::string& argumentString )
  /* this checks if argumentString beings with a numeric character,
   * & if so destinationInt is set to the value of the int made from the
   * unbroken sequences of numeric characters from the start of
   * argumentString up to the first non-numeric character, & these characters
   * are removed from the start of argumentString, & true is returned.
   * otherwise argumentString is left as it is & false is returned.
   */
  {
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "signalCalculator::separateOutIntegerSubstring( " << destinationInt
    << ", \"" << argumentString << "\" ) called.";
    std::cout << std::endl;**/

    size_t charPosition( 0 );
    while( ( '0' <= argumentString[ charPosition ] )
           &&
           ( '9' >= argumentString[ charPosition ] ) )
    {
      ++charPosition;
    }

    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "charPosition = " << charPosition;
    std::cout << std::endl;**/

    if( 0 < charPosition )
    {
      std::stringstream numberParser( argumentString.substr( 0,
                                                             charPosition ) );
      numberParser >> destinationInt;
      argumentString.erase( argumentString.begin(),
                            ( argumentString.begin() + charPosition ) );

      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "argumentString now = \"" << argumentString << "\"";
      std::cout << std::endl;**/

      return true;
    }
    else
    {
      return false;
    }
  }

  bool
  signalCalculator::separateOutDecimalSubstring( double& destinationDouble,
                                                 std::string& argumentString )
  /* this uses separateOutIntegerSubstring to check if argumentString begins
   * with a numeric character or '.'. destinationDouble is then set to be
   * equal to the int that separateOutIntegerSubstring found (or 0 if
   * argumentString begins with '.'), & then if there is an int after the '.'
   * (if there is a '.'), that becomes the part of destinationDouble after
   * the decimal point.
   */
  {
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "signalCalculator::separateOutDecimalSubstring( " << destinationDouble
    << ", \"" << argumentString << "\" ) called.";
    std::cout << std::endl;**/

    int intPart( 0 );
    bool intPartIsNotZero( separateOutIntegerSubstring( intPart,
                                                        argumentString ) );

    if( 0 == argumentString.compare( 0,
                                     1,
                                     "." ) )
    {
      argumentString.erase( argumentString.begin() );
      int decimalPart;
      if( separateOutIntegerSubstring( decimalPart,
                                       argumentString ) )
      {
        destinationDouble = (double)decimalPart;
        while( 1.0 <= destinationDouble )
        {
          destinationDouble *= 0.1;
        }
      }
      else
      {
        destinationDouble = 0.0;
      }
      destinationDouble += (double)intPart;

      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "destinationDouble set to " << destinationDouble;
      std::cout << std::endl;**/

      return true;
    }
    else if( intPartIsNotZero )
    {
      destinationDouble = (double)intPart;

      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "destinationDouble set to " << destinationDouble;
      std::cout << std::endl;**/

      return true;
    }
    else
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "destinationDouble left as " << destinationDouble;
      std::cout << std::endl;**/

      return false;
    }
  }

  bool
  signalCalculator::parseLeptonTransverseMomentumCuts(
                                                   std::string& argumentString,
                                 signalDefinitionSet* const signalDefinitions )
  /* this looks for "_pTl" then a double then "GeV", then a subsequent double,
   * & a subsequent "GeV", interpretted as the primary & secondary lepton cuts
   * respectively. it returns true if there is nothing more in argumentString,
   * false otherwise.
   */
  {
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "signalCalculator::parseLeptonTransverseMomentumCuts( "
    << argumentString << ", " << signalDefinitions
    << " ) called. no further output implies that no cuts were found.";
    std::cout << std::endl;**/

    if( findAndRemoveSubstring( "_pTl",
                                argumentString ) )
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "signalCalculator::parseLeptonTransverseMomentumCuts found \"_pTl\".";
      std::cout << std::endl;**/
      double leptonCut;
      if( separateOutDecimalSubstring( leptonCut,
                                       argumentString )
          &&
          ( 0 == argumentString.compare( 0,
                                         3,
                                         "GeV" ) ) )
      {
        argumentString.erase( argumentString.begin(),
                              ( argumentString.begin() + 3 ) );
        signalDefinitions->setLeptonCut( leptonCut );

        // debugging:
        /**std::cout << std::endl << "debugging:"
        << std::endl
        << "signalCalculator::parseLeptonTransverseMomentumCuts calling"
        << " signalDefinitions->setLeptonCut( " << leptonCut << " ).";
        std::cout << std::endl;**/

        if( separateOutDecimalSubstring( leptonCut,
                                         argumentString )
            &&
            ( 0 == argumentString.compare( 0,
                                           3,
                                           "GeV" ) ) )
        {
          argumentString.erase( argumentString.begin(),
                                ( argumentString.begin() + 3 ) );
          signalDefinitions->setSecondaryLeptonCut( leptonCut );

          // debugging:
          /**std::cout << std::endl << "debugging:"
          << std::endl
          << "signalCalculator::parseLeptonTransverseMomentumCuts calling"
          << " signalDefinitions->setSecondaryLeptonCut( " << leptonCut
          << " ).";
          std::cout << std::endl;**/
        }
      }
    }
    return argumentString.empty();
  }

  bool
  signalCalculator::findAndRemoveSubstring( std::string const& soughtSubstring,
                                            std::string& searchedString )
  /* this removes soughtSubstring from searchedString if found as a
   * substring, & returns true. if it is not found, searchedString is left
   * untouched & false is returned.
   */
  {
    size_t foundPosition( searchedString.find( soughtSubstring ) );
    if( std::string::npos != foundPosition )
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "signalCalculator::findAndRemoveSubstring found \"" << soughtSubstring
      << "\" in \"" << soughtSubstring << "\"";
      std::cout << std::endl;**/

      std::string reducedString( searchedString.substr( 0,
                                                        foundPosition ) );
      reducedString.append( searchedString,
                            ( foundPosition + soughtSubstring.size() ),
                            searchedString.size() );
      searchedString.assign( reducedString );

      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "now searchedString = \"" << searchedString << "\"";
      std::cout << std::endl;**/

      return true;
    }
    else
    {
      return false;
    }
  }


  namespace signalClasses
  {
    reallyWrongCalculator::reallyWrongCalculator() :
        signalCalculator()
    {
      // just an initialization list.
    }

    reallyWrongCalculator::~reallyWrongCalculator()
    {
      // does nothing.
    }



    sigmaBreakdownTest::sigmaBreakdownTest(
                         signalDefinitionSet const* const signalDefinitions ) :
        signalCalculator( signalDefinitions )
    {
      // just an initialization list.
    }

    sigmaBreakdownTest::~sigmaBreakdownTest()
    {
      // does nothing.
    }


    bool
    sigmaBreakdownTest::calculateValue( double* const signalValue,
                                        double* const uncertaintyFactor )
    // this over-rides the base version because it does things a bit
    // differently in a few of the steps.
    {
      // start by setting the signal value to 0.0:
      *signalValue = 0.0;

      for( std::vector< productionChannelPointerSet* >::iterator
           channelIterator( productionChannels.begin() );
           productionChannels.end() > channelIterator;
           ++channelIterator )
      {
        if( lhcFaserGlobal::negligibleSigma
            < (*channelIterator)->getCrossSection() )
          // if it's worth looking at this channel...
        {
          std::cout
          << std::endl
          << "channel: "
          << *((*channelIterator)->getScoloredPair()->getFirstParticle(
                                                       )->get_name_or_antiname(
                                           (*channelIterator)->getScoloredPair(
                                                 )->firstIsNotAntiparticle() ))
          << " + "
          << *((*channelIterator)->getScoloredPair()->getSecondParticle(
                                                       )->get_name_or_antiname(
                                           (*channelIterator)->getScoloredPair(
                                               )->secondIsNotAntiparticle() ));
          std::cout << std::endl;
          channelBrTotal = 0.0;

          // we have to look at all the open cascade pairings:
          firstCascades
          = (*channelIterator)->getFirstCascadeSet()->getOpenCascades();
          // debugging:
          /**std::cout << std::endl << "debugging:"
          << std::endl
          << "got open cascades from "
          << *((*channelIterator)->getFirstCascadeSet()->getInitialSparticle(
                                                                )->get_name());
          std::cout << std::endl;**/
          secondCascades
          = (*channelIterator)->getSecondCascadeSet()->getOpenCascades();

          // debugging:
          /**std::cout << std::endl << "debugging:"
          << std::endl
          << "got both open cascades.";
          std::cout << std::endl;**/


          for( std::vector< fullCascade* >::iterator
               firstCascadeIterator( firstCascades->begin() );
               firstCascades->end() > firstCascadeIterator;
               ++firstCascadeIterator )
          {
            firstCascadeBrToEwino =  (*firstCascadeIterator)->getBrToEwino(
                                               inputShortcut->getEmptyList() );
            if( ( lhcFaserGlobal::negligibleBr
                  * 2.0 * (double)(firstCascades->size()) )
                < firstCascadeBrToEwino )
              // this should mean we never throw away more than
              // 0.5 * lhcFaserGlobal::negligibleBr of acceptance.
            {
              for( std::vector< fullCascade* >::iterator
                  secondCascadeIterator( secondCascades->begin() );
                  secondCascades->end() > secondCascadeIterator;
                  ++secondCascadeIterator )
              {
                secondCascadeBrToEwino
                = (*secondCascadeIterator)->getBrToEwino(
                                               inputShortcut->getEmptyList() );
                if( ( lhcFaserGlobal::negligibleBr
                      * 2.0 * (double)(secondCascades->size()) )
                    < secondCascadeBrToEwino )
                  // this should mean we never throw away more than
                  // 0.5 * lhcFaserGlobal::negligibleBr of acceptance.
                {
                  // debugging:
                  /**std::cout << std::endl << "debugging:"
                  << std::endl
                  << "adding" << std::endl
                  << (*firstCascadeIterator)->getAsString() << " BR = "
                  << firstCascadeBrToEwino << std::endl
                  << (*secondCascadeIterator)->getAsString() << " BR = "
                  << secondCascadeBrToEwino << std::endl << " => "
                  << ( firstCascadeBrToEwino * secondCascadeBrToEwino );
                  std::cout << std::endl;**/

                  channelBrTotal += ( firstCascadeBrToEwino
                                      * secondCascadeBrToEwino );
                }
                // end of if the BR of the 2nd cascade to its electroweakino
                // was not negligible.
              }  // end of loop over 2nd cascade.
            }
            // end of if the BR of the 1st cascade to its electroweakino was
            // not negligible.
          }  // end of loop over 1st cascade.
          std::cout
          << "BR total: " << channelBrTotal;
          std::cout << std::endl;
          *signalValue
          += ( channelBrTotal * (*channelIterator)->getCrossSection() );
        }  // end of if the channel's cross-section was not negligible.
      }  // end of loop over productionChannels.
      return true;
    }



    int const
    atlasFourJetMetPlusGivenLeptonCuts::jetAcceptanceGridTableColumn( 10 );
    // this is dependent on the format of the grids.
    double const
    atlasFourJetMetPlusGivenLeptonCuts::defaultExtraJetCut( 40.0 );
    // this is the standard cut for the jets beyond the hardest cut for this
    // signal as implemented in this code.
    double const
    atlasFourJetMetPlusGivenLeptonCuts::defaultPrimaryLeptonCut( 20.0 );
    double const
    atlasFourJetMetPlusGivenLeptonCuts::defaultSecondaryLeptonCut( 10.0 );
    // the default Atlas lepton transverse momentum cuts are 20.0 GeV
    // (for leptons to *pass*) & 10.0 GeV (for all others to *fail*).

    signalCalculator*
    atlasFourJetMetPlusGivenLeptonCuts::getCalculator(
                                             std::string const& argumentString,
                                 signalDefinitionSet* const signalDefinitions )
    /* this either returns a pointer to a new
     * atlasFourJetMetPlusGivenLeptonCuts with cuts decided by the given
     * string, or NULL if the string could not be parsed properly.
     */
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "atlasFourJetMetPlusGivenLeptonCuts::getCalculator( "
      << argumentString << ", ... ) called.";
      std::cout << std::endl;**/

      signalCalculator* returnPointer( NULL );

      if( 0 == argumentString.compare( 0,
                                       11,
                                       "Atlas4jMET_" ) )
      {
        std::string validArguments( argumentString.substr( 10 ) );
        // debugging:
        /**std::cout << std::endl << "debugging:"
        << std::endl
        << "validArguments = \"" << validArguments << "\"";
        std::cout << std::endl;**/

        // the lepton transverse momentum cuts are set to their defaults so
        // that parseLeptonAcceptance can over-write them if necessary.
        signalDefinitions->setPrimaryLeptonCut( defaultPrimaryLeptonCut );
        signalDefinitions->setSecondaryLeptonCut( defaultSecondaryLeptonCut );

        if( parseBeamEnergy( validArguments,
                             signalDefinitions ) )
        {
          leptonAcceptanceForCascadePair*
          leptonAcceptanceType( parseLeptonAcceptance( validArguments,
                                                       signalDefinitions ) );

          // debugging:
          /**std::cout << std::endl << "debugging:"
          << std::endl
          << "atlasFourJetMetPlusGivenLeptonCuts::getCalculator( "
          << argumentString << ", ... ): leptonAcceptanceType = "
          << leptonAcceptanceType;
          std::cout << std::endl;**/

          if( NULL != leptonAcceptanceType )
          {
            if( parseLeptonTransverseMomentumCuts( validArguments,
                                                   signalDefinitions ) )
            {
              returnPointer
              = new atlasFourJetMetPlusGivenLeptonCuts( signalDefinitions,
                                                        leptonAcceptanceType );
              /* this takes care of if there were lepton transverse momentum
               * cuts specified. it's OK to only set them after the new
               * leptonAcceptanceForCascadePair has been constructed because
               * signalDefinitions is still a pointer to the
               * signalDefinitionSet that returnPointer uses.
               */

              // debugging:
              /**returnPointer->setSignalName( argumentString );
              std::cout << std::endl << "debugging:"
              << std::endl
              << "set signal name to \"" << argumentString << "\".";
              std::cout << std::endl;**/
            }
          }
        }
      }
      return returnPointer;
    }

    atlasFourJetMetPlusGivenLeptonCuts::atlasFourJetMetPlusGivenLeptonCuts(
                            signalDefinitionSet const* const signalDefinitions,
     leptonAcceptanceForCascadePair const* const leptonAcceptanceCalculator ) :
        signalCalculator( signalDefinitions ),
        leptonAcceptanceCalculator( leptonAcceptanceCalculator )
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "atlasFourJetMetPlusGivenLeptonCuts constructing with primary lepton"
      << " cut = "
      << signalDefinitions->getPrimaryLeptonCut()
      << ", secondary lepton cut = "
      << signalDefinitions->getSecondaryLeptonCut();
      std::cout << std::endl;**/

      std::string jetGridName( "Atlas4jMET" );
      fourJetKinematics
      = signalDefinitions->getShortcuts()->getJetPlusMetAcceptances(
                                               )->getJetPlusMetAcceptanceTable(
                                            signalDefinitions->getBeamEnergy(),
                                                                  &jetGridName,
                                                jetAcceptanceGridTableColumn );
      jetGridName.assign( "Atlas3jMET" );
      threeJetKinematics
      = signalDefinitions->getShortcuts()->getJetPlusMetAcceptances(
                                               )->getJetPlusMetAcceptanceTable(
                                            signalDefinitions->getBeamEnergy(),
                                                                  &jetGridName,
                                                jetAcceptanceGridTableColumn );
      jetGridName.assign( "Atlas2jMET" );
      twoJetKinematics
      = signalDefinitions->getShortcuts()->getJetPlusMetAcceptances(
                                               )->getJetPlusMetAcceptanceTable(
                                            signalDefinitions->getBeamEnergy(),
                                                                  &jetGridName,
                                                jetAcceptanceGridTableColumn );
      excludedFinalStateParticles.push_back( CppSLHA::PDG_code::top );
      excludedFinalStateParticles.push_back( -(CppSLHA::PDG_code::top) );
    }

    atlasFourJetMetPlusGivenLeptonCuts::~atlasFourJetMetPlusGivenLeptonCuts()
    {
      // does nothing.
    }


    double
    atlasFourJetMetPlusGivenLeptonCuts::valueForCurrentCascades(
                                                     fullCascade* firstCascade,
                                                   fullCascade* secondCascade )
    // see base version's description.
    {
      subchannelZeroOrMoreJets
      = leptonAcceptanceCalculator->withAtLeastNJets( signalDefinitions,
                                                      0,
                                                      firstCascade,
                                               firstSparticleIsNotAntiparticle,
                                                      secondCascade,
                                            secondSparticleIsNotAntiparticle );
        // chance of the pair of cascades producing the correct number of
        // leptons between them.

      subchannelOneOrMoreJets = 0.0;
      subchannelTwoOrMoreJets = 0.0;
      // these default to zero if not calculated.

      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "opening logStream.";
      std::cout << std::endl;
      std::string logFilename( signalName );
      logFilename.append( ".log" );
      std::ofstream logStream( logFilename.c_str(),
                               std::ios::app );
      if( !logStream.good() )
      {
        std::cout
        << std::endl
        << "error! could not open \"" << logFilename << "\"!";
        std::cout << std::endl;
        exit( EXIT_FAILURE );
      }
      if( 0.0001 < subchannelCrossSectionTimesBrToEwinos ){
      logStream << std::endl << "debugging:"
      << std::endl
      << "first cascade: not-anti = " << firstSparticleIsNotAntiparticle
      << ", cascade = " << firstCascade->getAsString()
      << std::endl
      << "second cascade: not-anti = " << secondSparticleIsNotAntiparticle
      << ", cascade = " << secondCascade->getAsString()
      << std::endl
      << "subchannelCrossSectionTimesBrToEwinos = "
      << subchannelCrossSectionTimesBrToEwinos
      << ", subchannelZeroOrMoreJets = "
      << subchannelZeroOrMoreJets
      << ", subchannelOneOrMoreJets = "
      << subchannelOneOrMoreJets
      << ", subchannelTwoOrMoreJets = "
      << subchannelTwoOrMoreJets;
      std::cout << std::endl;}**/

      if( lhcFaserGlobal::negligibleBr < subchannelZeroOrMoreJets )
      {
        fourJetAcceptance = fourJetKinematics->getAcceptance( initialPair,
                                                              firstCascade,
                                                              secondCascade );
        subchannelValue = ( subchannelCrossSectionTimesBrToEwinos
                            * subchannelZeroOrMoreJets
                            * fourJetAcceptance );

        // debugging:
        /**
        if( 0.0001 < subchannelCrossSectionTimesBrToEwinos ){
        logStream
        << std::endl
        << "subchannelCrossSectionTimesBrToEwinos = "
        << subchannelCrossSectionTimesBrToEwinos
        << ", subchannelZeroOrMoreJets = " << subchannelZeroOrMoreJets
        << ", fourJetAcceptance = " << fourJetAcceptance
        << ", subchannelValue = " << subchannelValue
        << ", firstCascadeBrToEwino = " << firstCascadeBrToEwino
        << ", secondCascadeBrToEwino = " << secondCascadeBrToEwino;
        std::cout << std::endl;}**/


        subchannelOneOrMoreJets
        = ( subchannelZeroOrMoreJets
           - leptonAcceptanceCalculator->withExactlyNJets( signalDefinitions,
                                                           0,
                                                           firstCascade,
                                               firstSparticleIsNotAntiparticle,
                                                           secondCascade,
                                          secondSparticleIsNotAntiparticle ) );
        // subchannelOneOrMoreJets = subchannelZeroOrMoreJets
        // minus the amount with zero cascade jets.

        // it's not going to ever be the case where the acceptance for 1+ jets
        // is greater than 0+ jets...
        if( lhcFaserGlobal::negligibleBr < subchannelOneOrMoreJets )
        {
          threeJetAcceptance
          = threeJetKinematics->getAcceptance( initialPair,
                                               firstCascade,
                                               secondCascade );
          subchannelValue
          += ( subchannelCrossSectionTimesBrToEwinos
               * subchannelOneOrMoreJets
               * ( threeJetAcceptance
                   - fourJetAcceptance ) );
          // we assume that the 3-jet+MET acceptance from the grid is greater
          // than the 4-jet+MET acceptance.

          // debugging:
          /**
          if( 0.0001 < subchannelCrossSectionTimesBrToEwinos ){
          logStream
          << std::endl
          << "subchannelCrossSectionTimesBrToEwinos = "
          << subchannelCrossSectionTimesBrToEwinos
          << ", subchannelOneOrMoreJets = " << subchannelOneOrMoreJets
          << ", threeJetAcceptance = " << threeJetAcceptance
          << ", subchannelValue = " << subchannelValue;
          std::cout << std::endl;}**/

          subchannelTwoOrMoreJets
          = ( subchannelOneOrMoreJets
             - leptonAcceptanceCalculator->withExactlyNJets( signalDefinitions,
                                                             1,
                                                             firstCascade,
                                               firstSparticleIsNotAntiparticle,
                                                             secondCascade,
                                          secondSparticleIsNotAntiparticle ) );
          // subchannelTwoOrMoreJets = subchannelOneOrMoreJets
          // minus the amount with one cascade jets.

          if( lhcFaserGlobal::negligibleBr < subchannelTwoOrMoreJets )
          {
            subchannelValue
            += ( subchannelCrossSectionTimesBrToEwinos
                 * subchannelTwoOrMoreJets
                 * ( twoJetKinematics->getAcceptance( initialPair,
                                                      firstCascade,
                                                      secondCascade )
                     - threeJetAcceptance ) );
            // we assume that the 2-jet+MET acceptance from the grid is greater
            // than the 3-jet+MET acceptance.

            // debugging:
            /**
            if( 0.0001 < subchannelCrossSectionTimesBrToEwinos ){
            logStream
            << std::endl
            << "subchannelCrossSectionTimesBrToEwinos = "
            << subchannelCrossSectionTimesBrToEwinos
            << ", subchannelTwoOrMoreJets = " << subchannelTwoOrMoreJets
            << ", twoJetAcceptance would be = "
            << twoJetKinematics->getAcceptance( initialPair,
                                                firstCascade,
                                                secondCascade )
            << ", subchannelValue = " << subchannelValue;
            std::cout << std::endl;}**/
          }
          else
          {
            subchannelTwoOrMoreJets = 0.0;
            /* this stops false alarms due to slightly negative values instead
             * of zero due to rounding errors, but does prevent real alarms
             * about crazily negative values.
             */
          }
        }
        else
        {
          subchannelOneOrMoreJets = 0.0;
          threeJetAcceptance = 0.0;
          /* this stops false alarms due to slightly negative values instead
           * of zero due to rounding errors, but does prevent real alarms
           * about crazily negative values.
           */
        }

        // debugging:
        /**
        if( 0.0001 < subchannelCrossSectionTimesBrToEwinos ){
        logStream
        << std::endl
        << "=> adding subchannelValue = " << subchannelValue;
        std::cout << std::endl;}**/

        if( !( ( 1.0 >= subchannelZeroOrMoreJets )
               &&
               ( 0.0 <= subchannelZeroOrMoreJets )
               &&
               ( 1.0 >= subchannelOneOrMoreJets )
               &&
               ( 0.0 <= subchannelOneOrMoreJets )
               &&
               ( 1.0 >= subchannelTwoOrMoreJets )
               &&
               ( 0.0 <= subchannelTwoOrMoreJets )
               &&
               ( 1.0 >= fourJetAcceptance )
               &&
               ( 0.0 <= fourJetAcceptance )
               &&
               ( 1.0 >= threeJetAcceptance )
               &&
               ( 0.0 <= threeJetAcceptance )
               &&
               ( 0.0 <= subchannelValue ) ) )
        {
          std::cout
          << std::endl << "LHC-FASER::error! signal \""
          << signalName << "\":"
          << std::endl
          << *(initialPair->getFirstParticle()->get_name_or_antiname(
                                             firstSparticleIsNotAntiparticle ))
          << " + "
          << *(initialPair->getSecondParticle()->get_name_or_antiname(
                                           secondSparticleIsNotAntiparticle ));
          std::cout
          << std::endl
          << "1st cascade: " << firstCascade->getAsString()
          << std::endl
          << "2nd cascade: " << secondCascade->getAsString();
          std::cout << std::endl;
          std::cout
          << std::endl
          << "produced subchannelZeroOrMoreJetsZeroLeptons = "
          << subchannelZeroOrMoreJets
          << ", subchannelOneOrMoreJetsZeroLeptons = "
          << subchannelOneOrMoreJets
          << ", subchannelTwoOrMoreJets = "
          << subchannelTwoOrMoreJets
          << ", fourJetAcceptance = "
          << fourJetAcceptance
          << ", & threeJetAcceptance = "
          << threeJetAcceptance
          << ", all of which should be >= 0.0 & <= 1.0 for valid"
          << " input! subchannelValue = " << subchannelValue
          << ", which should be >= 0.0 for valid input!";
          std::cout << std::endl;
        }
        return subchannelValue;
      }
      else
      {
        return 0.0;
      }
    }



    int const
    atlasThreeJetMetPlusGivenLeptonCuts::jetAcceptanceGridTableColumn( 10 );
    // this is dependent on the format of the grids.
    double const
    atlasThreeJetMetPlusGivenLeptonCuts::defaultExtraJetCut( 40.0 );
    // this is the standard cut for the jets beyond the hardest cut for this
    // signal as implemented in this code.
    double const
    atlasThreeJetMetPlusGivenLeptonCuts::defaultPrimaryLeptonCut( 20.0 );
    double const
    atlasThreeJetMetPlusGivenLeptonCuts::defaultSecondaryLeptonCut( 10.0 );
    // the default Atlas lepton transverse momentum cuts are 20.0 GeV
    // (for leptons to *pass*) & 10.0 GeV (for all others to *fail*).

    signalCalculator*
    atlasThreeJetMetPlusGivenLeptonCuts::getCalculator(
                                             std::string const& argumentString,
                                 signalDefinitionSet* const signalDefinitions )
    /* this either returns a pointer to a new
     * atlasThreeJetMetPlusGivenLeptonCuts with cuts decided by the given
     * string, or NULL if the string could not be parsed properly.
     */
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "atlasThreeJetMetPlusGivenLeptonCuts::getCalculator( "
      << argumentString << ", ... ) called.";
      std::cout << std::endl;**/

      signalCalculator* returnPointer( NULL );

      if( 0 == argumentString.compare( 0,
                                       11,
                                       "Atlas3jMET_" ) )
      {
        std::string validArguments( argumentString.substr( 10 ) );
        // debugging:
        /**std::cout << std::endl << "debugging:"
        << std::endl
        << "validArguments = \"" << validArguments << "\"";
        std::cout << std::endl;**/

        // the lepton transverse momentum cuts are set to their defaults so
        // that parseLeptonAcceptance can over-write them if necessary.
        signalDefinitions->setPrimaryLeptonCut( defaultPrimaryLeptonCut );
        signalDefinitions->setSecondaryLeptonCut( defaultSecondaryLeptonCut );

        if( parseBeamEnergy( validArguments,
                             signalDefinitions ) )
        {
          leptonAcceptanceForCascadePair*
          leptonAcceptanceType( parseLeptonAcceptance( validArguments,
                                                       signalDefinitions ) );
          if( NULL != leptonAcceptanceType )
          {
            if( parseLeptonTransverseMomentumCuts( validArguments,
                                                   signalDefinitions ) )
            {
              returnPointer
              = new atlasThreeJetMetPlusGivenLeptonCuts( signalDefinitions,
                                                        leptonAcceptanceType );
              /* this takes care of if there were lepton transverse momentum
               * cuts specified. it's OK to only set them after the new
               * leptonAcceptanceForCascadePair has been constructed because
               * signalDefinitions is still a pointer to the
               * signalDefinitionSet that returnPointer uses.
               */

              // debugging:
              /**returnPointer->setSignalName( argumentString );
              std::cout << std::endl << "debugging:"
              << std::endl
              << "set signal name to \"" << argumentString << "\".";
              std::cout << std::endl;**/
            }
          }
        }
      }
      return returnPointer;
    }

    atlasThreeJetMetPlusGivenLeptonCuts::atlasThreeJetMetPlusGivenLeptonCuts(
                            signalDefinitionSet const* const signalDefinitions,
     leptonAcceptanceForCascadePair const* const leptonAcceptanceCalculator ) :
        signalCalculator( signalDefinitions ),
        leptonAcceptanceCalculator( leptonAcceptanceCalculator )
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "atlasThreeJetMetPlusGivenLeptonCuts constructing with primary lepton"
      << " cut = "
      << signalDefinitions->getPrimaryLeptonCut()
      << ", secondary lepton cut = "
      << signalDefinitions->getSecondaryLeptonCut();
      std::cout << std::endl;**/

      std::string jetGridName( "Atlas3jMET" );
      threeJetKinematics
      = signalDefinitions->getShortcuts()->getJetPlusMetAcceptances(
                                               )->getJetPlusMetAcceptanceTable(
                                            signalDefinitions->getBeamEnergy(),
                                                                  &jetGridName,
                                                jetAcceptanceGridTableColumn );
      jetGridName.assign( "Atlas2jMET" );
      twoJetKinematics
      = signalDefinitions->getShortcuts()->getJetPlusMetAcceptances(
                                               )->getJetPlusMetAcceptanceTable(
                                            signalDefinitions->getBeamEnergy(),
                                                                  &jetGridName,
                                                jetAcceptanceGridTableColumn );
      excludedFinalStateParticles.push_back( CppSLHA::PDG_code::top );
      excludedFinalStateParticles.push_back( -(CppSLHA::PDG_code::top) );
    }

    atlasThreeJetMetPlusGivenLeptonCuts::~atlasThreeJetMetPlusGivenLeptonCuts()
    {
      // does nothing.
    }


    double
    atlasThreeJetMetPlusGivenLeptonCuts::valueForCurrentCascades(
                                                     fullCascade* firstCascade,
                                                   fullCascade* secondCascade )
    // see base version's description.
    {
      subchannelZeroOrMoreJets
      = leptonAcceptanceCalculator->withAtLeastNJets( signalDefinitions,
                                                      0,
                                                      firstCascade,
                                               firstSparticleIsNotAntiparticle,
                                                      secondCascade,
                                            secondSparticleIsNotAntiparticle );
        // chance of the pair of cascades producing the correct number of
        // leptons between them.

      subchannelOneOrMoreJets
      = ( subchannelZeroOrMoreJets
         - leptonAcceptanceCalculator->withExactlyNJets( signalDefinitions,
                                                         0,
                                                         firstCascade,
                                               firstSparticleIsNotAntiparticle,
                                                         secondCascade,
                                          secondSparticleIsNotAntiparticle ) );
      // subchannelOneOrMoreJets = subchannelZeroOrMoreJets
      // minus the amount with zero cascade jets.

      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "subchannelZeroOrMoreJets = "
      << subchannelZeroOrMoreJets
      << ", subchannelOneOrMoreJets = "
      << subchannelOneOrMoreJets;
      std::cout << std::endl;**/

      if( lhcFaserGlobal::negligibleBr < subchannelZeroOrMoreJets )
      {
        threeJetAcceptance = threeJetKinematics->getAcceptance( initialPair,
                                                                firstCascade,
                                                               secondCascade );
        subchannelValue = ( subchannelCrossSectionTimesBrToEwinos
                            * subchannelZeroOrMoreJets
                            * threeJetAcceptance );

        // debugging:
        /**if( 0.05 < ( subchannelCrossSectionTimesBrToEwinos
                         * subchannelZeroOrMoreJets ) )
        {
          std::cout << std::endl << "debugging:"
          << std::endl
          << "adding" << std::endl;
          if( !firstSparticleIsNotAntiparticle )
          {
            std::cout << "anti-";
          }
          std::cout
          << firstCascade->getAsString() << " BR = "
          << firstCascadeBrToEwino << std::endl;
          if( !secondSparticleIsNotAntiparticle )
          {
            std::cout << "anti-";
          }
          std::cout
          << secondCascade->getAsString() << " BR = " << secondCascadeBrToEwino
          << std::endl
          << " => " << ( firstCascadeBrToEwino * secondCascadeBrToEwino )
          << std::endl
          << "subchannelCrossSectionTimesBrToEwinos = "
          << subchannelCrossSectionTimesBrToEwinos
          << std::endl
          << "subchannelZeroOrMoreJets = "
          << subchannelZeroOrMoreJets;
          std::cout << std::endl;
        }**/

        // it's not going to ever be the case where the acceptance for 1+ jets
        // is greater than 0+ jets...
        if( lhcFaserGlobal::negligibleBr < subchannelOneOrMoreJets )
        {
          subchannelValue
          += ( subchannelCrossSectionTimesBrToEwinos
               * subchannelOneOrMoreJets
               * ( twoJetKinematics->getAcceptance( initialPair,
                                                    firstCascade,
                                                    secondCascade )
                   - threeJetAcceptance ) );
          // we assume that the 2-jet+MET acceptance from the grid is greater
          // than the 3-jet+MET acceptance.
        }
        // debugging:
        /**std::cout << std::endl << "debugging:"
        << std::endl
        << "subchannelCrossSectionTimesBrToEwinos = "
        << subchannelCrossSectionTimesBrToEwinos
        << std::endl
        << "subchannelZeroOrMoreJets = "
        << subchannelZeroOrMoreJets
        << std::endl
        << "fourJetAcceptance = " << fourJetAcceptance
        << std::endl
        << "=> adding subchannelValue = " << subchannelValue;
        std::cout << std::endl;**/
        if( !( ( 1.0 >= subchannelZeroOrMoreJets )
               &&
               ( 0.0 <= subchannelZeroOrMoreJets )
               &&
               ( 1.0 >= subchannelOneOrMoreJets )
               &&
               ( 0.0 <= subchannelOneOrMoreJets )
               &&
               ( 1.0 >= threeJetAcceptance )
               &&
               ( 0.0 <= threeJetAcceptance )
               &&
               ( 0.0 <= subchannelValue ) ) )
        {
          std::cout
          << std::endl << "LHC-FASER::error! signal \""
          << signalName << "\":"
          << std::endl
          << *(initialPair->getFirstParticle()->get_name_or_antiname(
                                             firstSparticleIsNotAntiparticle ))
          << " + "
          << *(initialPair->getSecondParticle()->get_name_or_antiname(
                                           secondSparticleIsNotAntiparticle ));
          std::cout
          << std::endl
          << "1st cascade: " << firstCascade->getAsString()
          << std::endl
          << "2nd cascade: " << secondCascade->getAsString();
          std::cout << std::endl;
          std::cout
          << std::endl
          << "produced subchannelZeroOrMoreJetsZeroLeptons = "
          << subchannelZeroOrMoreJets
          << ", subchannelOneOrMoreJetsZeroLeptons = "
          << subchannelOneOrMoreJets
          << ", & threeJetAcceptance = "
          << threeJetAcceptance
          << ", all of which should be >= 0.0 & <= 1.0 for valid"
          << " input! subchannelValue = " << subchannelValue
          << ", which should be >= 0.0 for valid input!";
          std::cout << std::endl;
        }
        return subchannelValue;
      }
      else
      {
        return 0.0;
      }
    }



    int const
    cmsTwoJetAlphaTPlusGivenLeptonCuts::jetAcceptanceGridTableColumn( 10 );
    // this is dependent on the format of the grids.
    double const
    cmsTwoJetAlphaTPlusGivenLeptonCuts::defaultPrimaryLeptonCut( 20.0 );
    double const
    cmsTwoJetAlphaTPlusGivenLeptonCuts::defaultSecondaryLeptonCut( 10.0 );
    // the default CMS lepton transverse momentum cuts are 20.0 GeV
    // (for leptons to *pass*) & 10.0 GeV (for all others to *fail*).

    signalCalculator*
    cmsTwoJetAlphaTPlusGivenLeptonCuts::getCalculator(
                                             std::string const& argumentString,
                                 signalDefinitionSet* const signalDefinitions )
    /* this either returns a pointer to a new
     * cmsTwoJetAlphaTPlusGivenLeptonCuts with cuts decided by the given
     * string, or NULL if the string could not be parsed properly.
     */
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "cmsTwoJetAlphaTPlusGivenLeptonCuts::getCalculator( "
      << argumentString << ", ... ) called.";
      std::cout << std::endl;**/

      signalCalculator* returnPointer( NULL );

      if( 0 == argumentString.compare( 0,
                                       12,
                                       "CMS2jalphaT_" ) )
      {
        std::string validArguments( argumentString.substr( 11 ) );
        // debugging:
        /**std::cout << std::endl << "debugging:"
        << std::endl
        << "validArguments = \"" << validArguments << "\"";
        std::cout << std::endl;**/

        // the lepton transverse momentum cuts are set to their defaults so
        // that parseLeptonAcceptance can over-write them if necessary.
        signalDefinitions->setPrimaryLeptonCut( defaultPrimaryLeptonCut );
        signalDefinitions->setSecondaryLeptonCut( defaultSecondaryLeptonCut );

        if( parseBeamEnergy( validArguments,
                             signalDefinitions ) )
        {
          leptonAcceptanceForCascadePair*
          leptonAcceptanceType( parseLeptonAcceptance( validArguments,
                                                       signalDefinitions ) );
          if( NULL != leptonAcceptanceType )
          {
            if( parseLeptonTransverseMomentumCuts( validArguments,
                                                   signalDefinitions ) )
            {
              returnPointer
              = new cmsTwoJetAlphaTPlusGivenLeptonCuts( signalDefinitions,
                                                        leptonAcceptanceType );
              /* this takes care of if there were lepton transverse momentum
               * cuts specified. it's OK to only set them after the new
               * leptonAcceptanceForCascadePair has been constructed because
               * signalDefinitions is still a pointer to the
               * signalDefinitionSet that returnPointer uses.
               */

              // debugging:
              /**returnPointer->setSignalName( argumentString );
              std::cout << std::endl << "debugging:"
              << std::endl
              << "set signal name to \"" << argumentString << "\".";
              std::cout << std::endl;**/
            }
          }
        }
      }
      return returnPointer;
    }

    cmsTwoJetAlphaTPlusGivenLeptonCuts::cmsTwoJetAlphaTPlusGivenLeptonCuts(
                            signalDefinitionSet const* const signalDefinitions,
     leptonAcceptanceForCascadePair const* const leptonAcceptanceCalculator ) :
        signalCalculator( signalDefinitions ),
        leptonAcceptanceCalculator( leptonAcceptanceCalculator )
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "cmsTwoJetAlphaTPlusGivenLeptonCuts constructing with primary lepton"
      << " cut = "
      << signalDefinitions->getPrimaryLeptonCut()
      << ", secondary lepton cut = "
      << signalDefinitions->getSecondaryLeptonCut();
      std::cout << std::endl;**/

      std::string jetGridName( "CMS2jalphaT" );
      twoJetKinematics
      = signalDefinitions->getShortcuts()->getJetPlusMetAcceptances(
                                               )->getJetPlusMetAcceptanceTable(
                                            signalDefinitions->getBeamEnergy(),
                                                                  &jetGridName,
                                                jetAcceptanceGridTableColumn );
      excludedFinalStateParticles.push_back( CppSLHA::PDG_code::top );
      excludedFinalStateParticles.push_back( -(CppSLHA::PDG_code::top) );
    }

    cmsTwoJetAlphaTPlusGivenLeptonCuts::~cmsTwoJetAlphaTPlusGivenLeptonCuts()
    {
      // does nothing.
    }


    double
    cmsTwoJetAlphaTPlusGivenLeptonCuts::valueForCurrentCascades(
                                                     fullCascade* firstCascade,
                                                   fullCascade* secondCascade )
    // see base version's description.
    {
      subchannelZeroOrMoreJets
      = leptonAcceptanceCalculator->withAtLeastNJets( signalDefinitions,
                                                      0,
                                                      firstCascade,
                                               firstSparticleIsNotAntiparticle,
                                                      secondCascade,
                                            secondSparticleIsNotAntiparticle );
      // chance of the pair of cascades producing the correct number of
      // leptons between them.

      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "subchannelZeroOrMoreJets = "
      << subchannelZeroOrMoreJets;
      std::cout << std::endl;**/

      if( lhcFaserGlobal::negligibleBr < subchannelZeroOrMoreJets )
      {
        twoJetAcceptance = twoJetKinematics->getAcceptance( initialPair,
                                                            firstCascade,
                                                            secondCascade );
        subchannelValue = ( subchannelCrossSectionTimesBrToEwinos
                            * subchannelZeroOrMoreJets
                            * twoJetAcceptance );

        // debugging:
        /**if( 0.05 < ( subchannelCrossSectionTimesBrToEwinos
                         * subchannelZeroOrMoreJets ) )
        {
          std::cout << std::endl << "debugging:"
          << std::endl
          << "adding" << std::endl;
          if( !firstSparticleIsNotAntiparticle )
          {
            std::cout << "anti-";
          }
          std::cout
          << firstCascade->getAsString() << " BR = "
          << firstCascadeBrToEwino << std::endl;
          if( !secondSparticleIsNotAntiparticle )
          {
            std::cout << "anti-";
          }
          std::cout
          << secondCascade->getAsString() << " BR = " << secondCascadeBrToEwino
          << std::endl
          << " => " << ( firstCascadeBrToEwino * secondCascadeBrToEwino )
          << std::endl
          << "subchannelCrossSectionTimesBrToEwinos = "
          << subchannelCrossSectionTimesBrToEwinos
          << std::endl
          << "subchannelZeroOrMoreJets = "
          << subchannelZeroOrMoreJets;
          std::cout << std::endl;
        }**/

        if( !( ( 1.0 >= subchannelZeroOrMoreJets )
               &&
               ( 0.0 <= subchannelZeroOrMoreJets )
               &&
               ( 1.0 >= twoJetAcceptance )
               &&
               ( 0.0 <= twoJetAcceptance )
               &&
               ( 0.0 <= subchannelValue ) ) )
        {
          std::cout
          << std::endl << "LHC-FASER::error! signal \""
          << signalName << "\":"
          << std::endl
          << *(initialPair->getFirstParticle()->get_name_or_antiname(
                                             firstSparticleIsNotAntiparticle ))
          << " + "
          << *(initialPair->getSecondParticle()->get_name_or_antiname(
                                           secondSparticleIsNotAntiparticle ));
          std::cout
          << std::endl
          << "1st cascade: " << firstCascade->getAsString()
          << std::endl
          << "2nd cascade: " << secondCascade->getAsString();
          std::cout << std::endl;
          std::cout
          << std::endl
          << "produced subchannelZeroOrMoreJetsZeroLeptons = "
          << subchannelZeroOrMoreJets
          << ", & twoJetAcceptance = "
          << twoJetAcceptance
          << ", both of which should be >= 0.0 & <= 1.0 for valid"
          << " input! subchannelValue = " << subchannelValue
          << ", which should be >= 0.0 for valid input!";
          std::cout << std::endl;
        }
        return subchannelValue;
      }
      else
      {
        return 0.0;
      }
    }



    double const
    noJetMetButGivenLeptonCuts::defaultPrimaryLeptonCut( 20.0 );
    double const
    noJetMetButGivenLeptonCuts::defaultSecondaryLeptonCut( 10.0 );
    // the default Atlas lepton transverse momentum cuts are 20.0 GeV
    // (for N leptons to *pass*) & 10.0 GeV (for all others to *fail*).

    signalCalculator*
    noJetMetButGivenLeptonCuts::getCalculator(
                                             std::string const& argumentString,
                                 signalDefinitionSet* const signalDefinitions )
    /* this either returns a pointer to a new atlasFourJetMetNLeptons with
     * cuts decided by the given string, or a pointer to a
     * reallyWrongCalculator.
     */
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "noJetMetButGivenLeptonCuts::getCalculator( "
      << argumentString << ", ... ) called.";
      std::cout << std::endl;**/

      signalCalculator* returnPointer( NULL );

      if( 0 == argumentString.compare( 0,
                                       14,
                                       "noJetOrMETCut_" ) )
      {
        std::string validArguments( argumentString.substr( 13 ) );
        // debugging:
        /**std::cout << std::endl << "debugging:"
        << std::endl
        << "validArguments = \"" << validArguments << "\"";
        std::cout << std::endl;**/

        // the lepton transverse momentum cuts are set to their defaults so
        // that parseLeptonAcceptance can over-write them if necessary.
        signalDefinitions->setPrimaryLeptonCut( defaultPrimaryLeptonCut );
        signalDefinitions->setSecondaryLeptonCut( defaultSecondaryLeptonCut );

        std::string argumentRemainder;
        if( parseBeamEnergy( validArguments,
                             signalDefinitions ) )
        {
          leptonAcceptanceForCascadePair*
          leptonAcceptanceType( parseLeptonAcceptance( validArguments,
                                                       signalDefinitions ) );
          if( NULL != leptonAcceptanceType )
          {
            if( parseLeptonTransverseMomentumCuts( validArguments,
                                                   signalDefinitions ) )
            {
              returnPointer
              = new noJetMetButGivenLeptonCuts( signalDefinitions,
                                                leptonAcceptanceType );
              /* this takes care of if there were lepton transverse momentum
               * cuts specified. it's OK to only set them after the new
               * leptonAcceptanceForCascadePair has been constructed because
               * signalDefinitions is still a pointer to the
               * signalDefinitionSet that returnPointer uses.
               */

              // debugging:
              /**returnPointer->setSignalName( argumentString );
              std::cout << std::endl << "debugging:"
              << std::endl
              << "set signal name to \"" << argumentString << "\".";
              std::cout << std::endl;**/
            }
          }
        }
      }
      return returnPointer;
    }

    noJetMetButGivenLeptonCuts::noJetMetButGivenLeptonCuts(
                            signalDefinitionSet const* const signalDefinitions,
     leptonAcceptanceForCascadePair const* const leptonAcceptanceCalculator ) :
        signalCalculator( signalDefinitions ),
        leptonAcceptanceCalculator( leptonAcceptanceCalculator )
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "noJetMetButGivenLeptonCuts constructing with primary lepton"
      << " cut = "
      << signalDefinitions->getPrimaryLeptonCut()
      << ", secondary lepton cut = "
      << signalDefinitions->getSecondaryLeptonCut();
      std::cout << std::endl;**/

      excludedFinalStateParticles.push_back( CppSLHA::PDG_code::top );
      excludedFinalStateParticles.push_back( -(CppSLHA::PDG_code::top) );
    }

    noJetMetButGivenLeptonCuts::~noJetMetButGivenLeptonCuts()
    {
      // does nothing.
    }

  }  // end of signalClasses namespace.



  signalHandler::signalHandler( std::string const signalName,
                                double const crossSectionUnitFactor,
                         signalDefinitionSet const* const signalDefinitions ) :
    getsReadiedForNewPoint( signalDefinitions->getShortcuts(
                                        )->getInputShortcuts()->getReadier() ),
    signalName( signalName ),
    rateCalculator( NULL ),
    inputShortcut( signalDefinitions->getShortcuts() ),
    signalPreparationDefinitions( signalDefinitions ),
    signalValue( CppSLHA::CppSLHA_global::really_wrong_value ),
    uncertaintyFactor( CppSLHA::CppSLHA_global::really_wrong_value ),
    crossSectionUnitFactor( crossSectionUnitFactor )
  {
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "signalHandler::signalHandler( " << signalName << ", "
    << crossSectionUnitFactor << ", " << signalDefinitions << " ) called.";
    std::cout << std::endl;**/

    std::vector< signalCalculatorCreator > creationFunctions;
    creationFunctions.push_back(
         &(signalClasses::atlasFourJetMetPlusGivenLeptonCuts::getCalculator) );
    creationFunctions.push_back(
        &(signalClasses::atlasThreeJetMetPlusGivenLeptonCuts::getCalculator) );
    creationFunctions.push_back(
         &(signalClasses::cmsTwoJetAlphaTPlusGivenLeptonCuts::getCalculator) );
    creationFunctions.push_back(
                 &(signalClasses::noJetMetButGivenLeptonCuts::getCalculator) );
    creationFunctions.push_back(
                         &(signalClasses::sigmaBreakdownTest::getCalculator) );

    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "0 = Atlas4j, 1 = Atlas3j, 2 = noJetMet, 3 = sigma breakdown";
    std::cout << std::endl;**/

    for( int whichFunction( (int)(creationFunctions.size()) - 1 );
         ( ( 0 <= whichFunction )
           &&
         ( NULL == rateCalculator ) );
         --whichFunction )
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "calling for creationFunctions[ " << whichFunction << " ]";
      std::cout << std::endl;**/

      rateCalculator = (*(creationFunctions[ whichFunction ]))( signalName,
                                               &signalPreparationDefinitions );

      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "now rateCalculator = " << rateCalculator;
      std::cout << std::endl;**/
    }
    if( NULL == rateCalculator )
    {
      std::cout
      << std::endl
      << "LHC-FASER::error! signalHandler::signalHandler( "
      << signalName << ", " <<  inputShortcut << " ) was passed a"
      << " name it does not know ( \"" << signalName
      << "\" ). its calculator is being set to return only "
      << CppSLHA::CppSLHA_global::really_wrong_value_string << " for"
      << " every point.";
      std::cout << std::endl;

      rateCalculator = new signalClasses::reallyWrongCalculator();
    }
    rateCalculator->setSignalName( signalName );
  }

  signalHandler::~signalHandler()
  {
    delete rateCalculator;
  }

}  // end of LHC_FASER namespace.
