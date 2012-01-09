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
  signalCalculator::signalCalculator(
                               signalDefinitionSet* const signalDefinitions ) :
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
    for( std::vector< signedParticleShortcutPair* >::const_iterator
         pairIterator( signalDefinitions->getShortcuts()->getInputShortcuts(
                             )->getScoloredProductionCombinations()->begin() );
         signalDefinitions->getShortcuts()->getInputShortcuts(
                                  )->getScoloredProductionCombinations()->end()
         > pairIterator;
         ++pairIterator )
    {
      productionChannels.push_back( new productionChannelPointerSet(
                                                             signalDefinitions,
                                                             *pairIterator ) );
    }
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
    << " signalDefinitions.getPrimaryLeptonCut() = "
    << signalDefinitions.getPrimaryLeptonCut();
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
          << "1st: " << (*firstCascadeIterator)->getAsString();
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
    = ( *signalValue * signalDefinitions.getShortcuts()->getUncertainty() );
    return true;
  }

  void
  signalCalculator::parseLeptonTransverseMomentumCuts(
                                                  std::string const& cutString,
                                  signalDefinitionSet* const signalDefinitions,
                                          double const defaultPrimaryLeptonCut,
                                       double const defaultSecondaryLeptonCut )
  /* this looks for "pTl" then a double, then a subsequent "pTl" & a
   * subsequent double, interpretted as the primary & secondary lepton cuts
   * respectively.
   */
  {
    if( 0 == cutString.compare( 0,
                                3,
                                "pTl" ) )
      // if the 1st 3 characters denote a non-default primary lepton
      // transverse momentum cut...
    {
      /* the argument should be
       * "pTl[primary lepton cut]pTl[secondary lepton cut]
       * with the 2nd cut optional.
       */
      double leptonCut;
      std::string leptonCutSubstring( cutString,
                                      3 );
      // leptonCutSubstring is arguments with the
      // "pTl" bit removed.
      int failingLeptonCutPosition( leptonCutSubstring.find( "pTl" ) );
      // failingLeptonCutPosition is now the number of characters in
      // leptonCutSubstring before "pTl" (the 2nd "pTl" in arguments).
      if( std::string::npos != failingLeptonCutPosition )
      {
        std::stringstream
        leptonCutStream( leptonCutSubstring.substr( 0,
                                                failingLeptonCutPosition ) );
        // the substring provided to leptonCutStream is just all the
        // characters in leptonCutSubstring before "pTl".
        leptonCutStream >> leptonCut;
        signalDefinitions->setPrimaryLeptonCut( leptonCut );
        // these characters are now translated into a number which is passed
        // to signalDefinitions as the primary lepton cut.
        leptonCutStream.clear();
        leptonCutStream.str( leptonCutSubstring.substr(
                                          ( failingLeptonCutPosition + 3 ) ) );
        /* the substring now provided to leptonCutStream is all the
         * characters in leptonCutSubstring after "pTl", which should just be
         * the secondary lepton transverse momentum cut.
         */
        // debugging:
        /**std::cout << std::endl << "debugging:"
        << std::endl
        << "for secondary cut, leptonCutStream.str() = "
        << leptonCutStream.str();
        std::cout << std::endl;**/

        leptonCutStream >> leptonCut;
      }
      else
      {
        std::stringstream
        leptonCutStream( leptonCutSubstring );
        // if there is no subsequent "pTl", the remaining characters should
        // just define the lepton cut used for both.
        leptonCutStream >> leptonCut;
        signalDefinitions->setPrimaryLeptonCut( leptonCut );
        // these characters are now translated into a number which is passed
        // to signalDefinitions as the primary lepton cut.
      }
      /* if there was a secondary lepton cut specified, it is now in leptonCut,
       * otherwise leptonCut is still the same value as specified for the
       * primary cut.
       */
      signalDefinitions->setSecondaryLeptonCut( leptonCut );
    }
    else
    {
      signalDefinitions->setPrimaryLeptonCut( defaultPrimaryLeptonCut );
      signalDefinitions->setSecondaryLeptonCut( defaultSecondaryLeptonCut );
    }
  }



  signalHandler::signalHandler( std::string const signalName,
                                double const crossSectionUnitFactor,
                    signalDefinitionSet* const signalPreparationDefinitions ) :
    getsReadiedForNewPoint( signalPreparationDefinitions->getShortcuts(
                                        )->getInputShortcuts()->getReadier() ),
    signalName( signalName ),
    rateCalculator( NULL ),
    inputShortcut( signalPreparationDefinitions->getShortcuts() ),
    signalValue( CppSLHA::CppSLHA_global::really_wrong_value ),
    uncertaintyFactor( CppSLHA::CppSLHA_global::really_wrong_value ),
    crossSectionUnitFactor( crossSectionUnitFactor )
  {
    std::string calculatorArgument( "error" );
    if( 0 == signalName.compare( 0,
                                 10,
                                 "Atlas4jMET" ) )
      // if the signal is the Atlas 4-jet, missing transverse momentum,
      // 0 leptons event rate...
    {
      calculatorArgument.assign( ( signalName.begin() + 10 ),
                                 signalName.end() );
      rateCalculator
      = signalCalculatorClasses::atlasFourJetMetNLeptons::getCalculator(
                                                           &calculatorArgument,
                                                signalPreparationDefinitions );
    }
    else if( 0 == signalName.compare( 0,
                                      12,
                                      "Atlas3jMET1l" ) )
      // if the signal is the Atlas 3-jet, missing transverse momentum,
      // 1 lepton event rate...
    {
      calculatorArgument.assign( ( signalName.begin() + 12 ),
                                 signalName.end() );
      rateCalculator
      = signalCalculatorClasses::atlasThreeJetMetNLeptons::getCalculator(
                                                           &calculatorArgument,
                                                signalPreparationDefinitions );
    }
    else if( 0 == signalName.compare( "CMS2jMETanyl14TeV" ) )
      // if the signal is the CMS 2-jet, missing transverse momentum,
      // any amount of leptons event rate, for a beam energy of 14 TeV...
    {
      std::cout
      << std::endl
      << "LHC-FASER::error! signalHandler::signalHandler( "
      << signalName << ", "<<  inputShortcut << " ) was passed"
      << " CMS2jMETanyl14TeV, but unfortunately it has not yet been"
      << " implemented. its calculator is being set to return only "
      << CppSLHA::CppSLHA_global::really_wrong_value_string
      << " for every point instead.";

      rateCalculator
      = new signalCalculatorClasses::reallyWrongCalculator(
                                                signalPreparationDefinitions );
    }
    else if( 0 == signalName.compare( 0,
                                      16,
                                      "sameSignDilepton" ) )
      // if the signal is the same-sign dilepton event rate...
    {
      calculatorArgument.assign( ( signalName.begin() + 16 ),
                                 signalName.end() );
      rateCalculator
      = signalCalculatorClasses::sameSignDilepton::getCalculator(
                                                           &calculatorArgument,
                                                signalPreparationDefinitions );
    }
    else if( 0 == signalName.compare( "sigmaBreakdownTest" ) )
      // if the signal is the test to see if the cross-section breakdown works
      // as it should...
    {
      signalPreparationDefinitions->setBeamEnergy( 7 );
      rateCalculator
      = new signalCalculatorClasses::sigmaBreakdownTest(
                                                signalPreparationDefinitions );
    }
    else
      // otherwise, a malformed signal name was passed:
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

      rateCalculator
      = new signalCalculatorClasses::reallyWrongCalculator(
                                                signalPreparationDefinitions );
    }
  }

  signalHandler::~signalHandler()
  {
    delete rateCalculator;
  }



  namespace signalCalculatorClasses
  {
    reallyWrongCalculator::reallyWrongCalculator(
                               signalDefinitionSet* const signalDefinitions ) :
        signalCalculator( signalDefinitions )
    {
      // just an initialization list.
    }

    reallyWrongCalculator::~reallyWrongCalculator()
    {
      // does nothing.
    }



    sigmaBreakdownTest::sigmaBreakdownTest(
                               signalDefinitionSet* const signalDefinitions ) :
        signalCalculator( signalDefinitions )
    {
      this->signalDefinitions.setExcludedStandardModelProducts(
                          signalDefinitions->getShortcuts()->getInputShortcuts(
                                                          )->getNotInJets5() );
      //this->signalDefinitions.setBeamEnergy( 7 );
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
    atlasFourJetMetNLeptons::jetAcceptanceGridTableColumn( 10 );
    // this is dependent on the format of the grids.
    double const
    atlasFourJetMetNLeptons::defaultExtraJetCut( 40.0 );
    // this is the standard cut for the jets beyond the hardest cut for this
    // signal as implemented in this code.
    double const
    atlasFourJetMetNLeptons::defaultLeptonCut( 10.0 );
    // the default Atlas4jMET0l lepton transverse momentum cut is 10.0 GeV.

    signalCalculator*
    atlasFourJetMetNLeptons::getCalculator(
                                       std::string const* const argumentString,
                                 signalDefinitionSet* const signalDefinitions )
    /* this either returns a pointer to a new atlasFourJetMetNLeptons with
     * cuts decided by the given string, or a pointer to a
     * reallyWrongCalculator.
     */
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "atlasFourJetMetNLeptons::getCalculator( " << *argumentString
      << ", ... ) called.";
      std::cout << std::endl;**/

      signalDefinitions->setExcludedStandardModelProducts(
                          signalDefinitions->getShortcuts()->getInputShortcuts(
                                                          )->getNotInJets5() );
      int numberOfLeptons;
      signalCalculator* returnPointer( NULL );
      unsigned int characterSkipper;
      // this is used in interpreting the arguments string.
      if( 0 == argumentString->compare( 0,
                                        2,
                                        "0l" ) )
        // for zero leptons...
      {
        numberOfLeptons = 0;
        characterSkipper = 2;
      }
      else if( 0 == argumentString->compare( 0,
                                             2,
                                             "1l" ) )
        // for exactly one lepton...
      {
        numberOfLeptons = 1;
        characterSkipper = 2;
      }
      else if( 0 == argumentString->compare( 0,
                                             2,
                                             "2l" ) )
        // for exactly two leptons...
      {
        numberOfLeptons = 1;
        characterSkipper = 2;
      }
      else
      {
        std::cout
        << std::endl
        << "LHC-FASER::error! signalHandler::signalHandler( "
        << argumentString << ", " <<  signalDefinitions << " ) was passed"
        << " a name it does not know (started with \"Atlas4jMET\", but then"
        << " \"0l\", \"1l\", or \"2l\" are expected as the next characters,"
        << " before the energy). its calculator is being set to return "
        << CppSLHA::CppSLHA_global::really_wrong_value_string
        << " for every point.";
        returnPointer = new reallyWrongCalculator( signalDefinitions );
      }
      if( NULL == returnPointer )
        // if the signalCalculator hasn't been set to a reallyWrongCalculator...
      {
        // debugging:
        /**std::cout << std::endl << "debugging:"
        << std::endl
        << "characterSkipper = " << characterSkipper
        << ", substring = " << argumentString->substr( characterSkipper,
                                                       2 );
        std::cout << std::endl;**/

        if( 0 == argumentString->compare( characterSkipper,
                                          4,
                                          "7TeV" ) )
          // for a beam energy of 7 TeV...
        {
          signalDefinitions->setBeamEnergy( 7 );
          characterSkipper += 4;
        }
        else if( 0 == argumentString->compare( characterSkipper,
                                               5,
                                               "07TeV" ) )
          // for a beam energy of 7 TeV...
        {
          signalDefinitions->setBeamEnergy( 7 );
          characterSkipper += 5;
        }
        else if( 0 == argumentString->compare( characterSkipper,
                                               5,
                                               "10TeV" ) )
          // for a beam energy of 10 TeV...
        {
          signalDefinitions->setBeamEnergy( 10 );
          characterSkipper += 5;
        }
        else if( 0 == argumentString->compare( characterSkipper,
                                               5,
                                               "14TeV" ) )
          // for a beam energy of 14 TeV...
        {
          signalDefinitions->setBeamEnergy( 14 );
          characterSkipper += 5;
        }
        else
        {
          std::cout
          << std::endl
          << "LHC-FASER::error! signalHandler::signalHandler( "
          << argumentString << ", " <<  signalDefinitions << " ) was passed"
          << " a name it does not know (started with \"Atlas4jMET"
          << numberOfLeptons << "l\", but then \"7TeV\", \"07TeV\", \"10TeV\""
          << " or \"14TeV\" are expected as the next characters). its"
          << " calculator is being set to return "
          << CppSLHA::CppSLHA_global::really_wrong_value_string
          << " for every point.";
          returnPointer = new reallyWrongCalculator( signalDefinitions );
        }
      }
      if( NULL == returnPointer )
        // if the signalCalculator hasn't been set to a reallyWrongCalculator...
      {
        // debugging:
        /**std::cout << std::endl << "debugging:"
        << std::endl
        << "characterSkipper = " << characterSkipper
        << ", substring = " << argumentString->substr( characterSkipper,
                                                       3 );
        std::cout << std::endl;**/

        parseLeptonTransverseMomentumCuts(
                                    argumentString->substr( characterSkipper ),
                                           signalDefinitions,
                                           defaultPrimaryLeptonCut,
                                           defaultSecondaryLeptonCut );
        signalDefinitions->setJetCut( defaultExtraJetCut );
        returnPointer = new atlasFourJetMetNLeptons( signalDefinitions );
      }
      return returnPointer;
    }

    atlasFourJetMetNLeptons::atlasFourJetMetNLeptons(
                               signalDefinitionSet* const signalDefinitions,
                                                  int const numberOfLeptons ) :
        signalCalculator( signalDefinitions ),
        numberOfLeptons( numberOfLeptons )
    {
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
      excludedFinalStateParticles.push_back( CppSLHA::PDG_code::top );
      excludedFinalStateParticles.push_back( -(CppSLHA::PDG_code::top) );
    }

    atlasFourJetMetNLeptons::~atlasFourJetMetNLeptons()
    {
      // does nothing.
    }


    double
    atlasFourJetMetNLeptons::valueForCurrentCascades(
                                                     fullCascade* firstCascade,
                                                   fullCascade* secondCascade )
    // see base version's description.
    {
      subchannelZeroOrMoreJetsNLeptons = 0.0;
      subchannelOneOrMoreJetsNLeptons = 0.0;
      for( int leptonsInFirstCascade( 0 );
           numberOfLeptons > leptonsInFirstCascade;
           ++leptonsInFirstCascade )
      {
        subchannelZeroOrMoreJetsNLeptons
        += ( firstCascade->unspecifiedJetsSpecifiedChargeSummedLeptons(
                                                            &signalDefinitions,
                                                        leptonsInFirstCascade )
             * secondCascade->unspecifiedJetsSpecifiedChargeSummedLeptons(
                                                            &signalDefinitions,
                               ( numberOfLeptons - leptonsInFirstCascade ) ) );
        // chance of the pair of cascades producing the correct number of
        // leptons between them.

        subchannelOneOrMoreJetsNLeptons
        += ( subchannelZeroOrMoreJetsNLeptons
             - ( firstCascade->specifiedJetsSpecifiedChargeSummedLeptons(
                                                            &signalDefinitions,
                                                                          0,
                                                        leptonsInFirstCascade )
                 * secondCascade->specifiedJetsSpecifiedChargeSummedLeptons(
                                                            &signalDefinitions,
                                                                             0,
                             ( numberOfLeptons - leptonsInFirstCascade ) ) ) );
        // chance of no lepton rejection removing the portion without any extra
        // jet passing the cut.
      }

      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "subchannelZeroOrMoreJetsNLeptons = "
      << subchannelZeroOrMoreJetsNLeptons
      << ", subchannelOneOrMoreJetsNLeptons = "
      << subchannelOneOrMoreJetsNLeptons;
      std::cout << std::endl;**/

      if( lhcFaserGlobal::negligibleBr < subchannelZeroOrMoreJetsNLeptons )
      {
        fourJetAcceptance = fourJetKinematics->getAcceptance( initialPair,
                                                              firstCascade,
                                                              secondCascade );
        subchannelValue = ( subchannelCrossSectionTimesBrToEwinos
                            * subchannelZeroOrMoreJetsNLeptons
                            * fourJetAcceptance );

        // debugging:
        /**if( 0.05 < ( subchannelCrossSectionTimesBrToEwinos
                         * subchannelZeroOrMoreJetsNLeptons ) )
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
          << "subchannelZeroOrMoreJetsNLeptons = "
          << subchannelZeroOrMoreJetsNLeptons;
          std::cout << std::endl;
        }**/

        // it's not going to ever be the case where the acceptance for 1+ jets
        // is greater than 0+ jets...
        if( lhcFaserGlobal::negligibleBr < subchannelOneOrMoreJetsNLeptons )
        {
          subchannelValue
          += ( subchannelCrossSectionTimesBrToEwinos
               * subchannelOneOrMoreJetsNLeptons
               * ( threeJetKinematics->getAcceptance( initialPair,
                                                      firstCascade,
                                                      secondCascade )
                   - fourJetAcceptance ) );
          // we assume that the 3-jet+MET acceptance from the grid is greater
          // than the 4-jet+MET acceptance.
        }
        // debugging:
        /**std::cout << std::endl << "debugging:"
        << std::endl
        << "subchannelCrossSectionTimesBrToEwinos = "
        << subchannelCrossSectionTimesBrToEwinos
        << std::endl
        << "subchannelZeroOrMoreJetsNLeptons = "
        << subchannelZeroOrMoreJetsNLeptons
        << std::endl
        << "fourJetAcceptance = " << fourJetAcceptance
        << std::endl
        << "=> adding subchannelValue = " << subchannelValue;
        std::cout << std::endl;**/
        if( !( ( 1.0 >= subchannelZeroOrMoreJetsNLeptons )
               &&
               ( 0.0 <= subchannelZeroOrMoreJetsNLeptons )
               &&
               ( 1.0 >= subchannelOneOrMoreJetsNLeptons )
               &&
               ( 0.0 <= subchannelOneOrMoreJetsNLeptons )
               &&
               ( 1.0 >= fourJetAcceptance )
               &&
               ( 0.0 <= fourJetAcceptance )
               &&
               ( 0.0 <= subchannelValue ) ) )
        {
          std::cout << std::endl << "LHC-FASER::error!"
          << std::endl
          << *(initialPair->getFirstParticle()->get_name_or_antiname(
                                             firstSparticleIsNotAntiparticle ))
          << " + "
          << *(initialPair->getSecondParticle()->get_name_or_antiname(
                                           secondSparticleIsNotAntiparticle ));
          std::cout
          << std::endl
          << "1st: " << firstCascade->getAsString() << "; 0+j, Nl = "
          << firstCascade->unspecifiedJetsSpecifiedChargeSummedLeptons(
                                                            &signalDefinitions,
                                                                        0 );
          std::cout
          << std::endl
          << "2nd: " << secondCascade->getAsString() << "; 0+j, Nl = "
          << secondCascade->unspecifiedJetsSpecifiedChargeSummedLeptons(
                                                            &signalDefinitions,
                                                                         0 );
          std::cout << std::endl;
          std::cout
          << std::endl
          << "produced subchannelZeroOrMoreJetsZeroLeptons = "
          << subchannelZeroOrMoreJetsNLeptons
          << ", subchannelOneOrMoreJetsZeroLeptons = "
          << subchannelOneOrMoreJetsNLeptons
          << ", & fourJetAcceptance = "
          << fourJetAcceptance
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



    int const atlasThreeJetMetNLeptons::jetAcceptanceGridTableColumn( 10 );
    // this is dependent on the format of the grids.
    double const atlasThreeJetMetNLeptons::defaultExtraJetCut( 40.0 );
    // this is the standard cut for the jets beyond the hardest cut for this
    // signal as implemented in this code.
    double const atlasThreeJetMetNLeptons::defaultPrimaryLeptonCut( 20.0 );
    double const atlasThreeJetMetNLeptons::defaultSecondaryLeptonCut( 10.0 );
    // the default Atlas4jMET0l lepton transverse momentum cuts are 20.0 GeV
    // (for a single lepton to *pass*) & 10.0 GeV (for all others to *fail*).

    signalCalculator*
    atlasThreeJetMetNLeptons::getCalculator(
                                       std::string const* const argumentString,
                                 signalDefinitionSet* const signalDefinitions )
    {
      signalDefinitions->setExcludedStandardModelProducts(
                          signalDefinitions->getShortcuts()->getInputShortcuts(
                                                          )->getNotInJets5() );
      signalCalculator* returnPointer( NULL );
      int characterSkipper;
      // this is used in interpreting the arguments string.
      if( 0 == argumentString->compare( 0,
                                        4,
                                        "7TeV" ) )
        //for a beam energy of 7 TeV...
      {
        signalDefinitions->setBeamEnergy( 7 );
        characterSkipper = 4;
      }
      else if( 0 == argumentString->compare( 0,
                                             5,
                                             "07TeV" ) )
        //for a beam energy of 7 TeV...
      {
        signalDefinitions->setBeamEnergy( 7 );
        characterSkipper = 5;
      }
      else if( 0 == argumentString->compare( 0,
                                             5,
                                             "10TeV" ) )
        //for a beam energy of 10 TeV...
      {
        signalDefinitions->setBeamEnergy( 10 );
        characterSkipper = 5;
      }
      else if( 0 == argumentString->compare( 0,
                                             5,
                                             "14TeV" ) )
        //for a beam energy of 14 TeV...
      {
        signalDefinitions->setBeamEnergy( 14 );
        characterSkipper = 5;
      }
      else
      {
        std::cout
        << std::endl
        << "LHC-FASER::error! signalHandler::signalHandler( "
        << argumentString << ", " <<  signalDefinitions << " ) was passed"
        << " a name it does not know (started with \"Atlas3jMET\", but then"
        << " \"0l\", \"1l\", or \"2l\" are expected as the next characters,"
        << " before the energy). its calculator is being set to return "
        << CppSLHA::CppSLHA_global::really_wrong_value_string
        << " for every point.";
        returnPointer = new reallyWrongCalculator( signalDefinitions );
      }
      if( NULL == returnPointer )
        // if the signalCalculator hasn't been set to a reallyWrongCalculator...
      {
        parseLeptonTransverseMomentumCuts(
                                    argumentString->substr( characterSkipper ),
                                           signalDefinitions,
                                           defaultPrimaryLeptonCut,
                                           defaultSecondaryLeptonCut );
        signalDefinitions->setJetCut( defaultExtraJetCut );
        returnPointer = new atlasThreeJetMetNLeptons( signalDefinitions );
      }
      return returnPointer;
    }

    atlasThreeJetMetNLeptons::atlasThreeJetMetNLeptons(
                                  signalDefinitionSet* const signalDefinitions,
                                                  int const numberOfLeptons ) :
        signalCalculator( signalDefinitions ),
        numberOfLeptons( numberOfLeptons )
    {
      std::string jetGridName( "Atlas3jMET" );
      this->signalDefinitions.setJetPlusMetAcceptance(
                   signalDefinitions->getShortcuts()->getJetPlusMetAcceptances(
                                               )->getJetPlusMetAcceptanceTable(
                                            signalDefinitions->getBeamEnergy(),
                                                                  &jetGridName,
                                              jetAcceptanceGridTableColumn ) );
      excludedFinalStateParticles.push_back( CppSLHA::PDG_code::top );
      excludedFinalStateParticles.push_back( -(CppSLHA::PDG_code::top) );
    }

    atlasThreeJetMetNLeptons::~atlasThreeJetMetNLeptons()
    {
      // does nothing.
    }


    double
    atlasThreeJetMetNLeptons::valueForCurrentCascades(
                                                     fullCascade* firstCascade,
                                                   fullCascade* secondCascade )
    // see base version's description.
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "looking at " << std::endl;
      if( !firstSparticleIsNotAntiparticle )
      {
        std::cout << "anti-";
      }
      std::cout
      << firstCascade->getAsString() << " BR = " << firstCascadeBrToEwino
      << std::endl;
      if( !secondSparticleIsNotAntiparticle )
      {
        std::cout << "anti-";
      }
      std::cout
      << secondCascade->getAsString() << " BR = " << secondCascadeBrToEwino
      << std::endl << " => "
      << ( firstCascadeBrToEwino * secondCascadeBrToEwino )
      << std::endl
      << "subchannelCrossSectionTimesBrToEwinos = "
      << subchannelCrossSectionTimesBrToEwinos
      << std::endl
      << "subchannelZeroOrMoreJetsNLeptons = "
      << subchannelZeroOrMoreJetsNLeptons;
      std::cout << std::endl;**/

      subchannelZeroOrMoreJetsNLeptons
      = ( firstCascade->unspecifiedJetsSpecifiedChargeSummedLeptons(
                                                            &signalDefinitions,
                                                                     1 )
          * secondCascade->unspecifiedJetsSpecifiedChargeSummedLeptons(
                                                            &signalDefinitions,
                                                                        0 )
          + firstCascade->unspecifiedJetsSpecifiedChargeSummedLeptons(
                                                            &signalDefinitions,
                                                                       0 )
            * secondCascade->unspecifiedJetsSpecifiedChargeSummedLeptons(
                                                            &signalDefinitions,
                                                                         1 ) );
      /* chance for 1st cascade to produce the lepton * chance of 2nd to not
       * produce any further leptons (which would cause the event to be
       * rejected) + chance of lepton from 2nd cascade * chance of no lepton
       * from 1st cascade.
       */
      // debugging:
      /**if( 0.01 < ( subchannelCrossSectionTimesBrToEwinos
                       * subchannelZeroOrMoreJetsNLeptons ) )
      {
        std::cout << std::endl << "debugging:"
        << std::endl
        << "adding" << std::endl;
        if( !firstSparticleIsNotAntiparticle )
        {
          std::cout << "anti-";
        }
        std::cout
        << firstCascade->getAsString() << " BR = " << firstCascadeBrToEwino
        << std::endl;
        if( !secondSparticleIsNotAntiparticle )
        {
          std::cout << "anti-";
        }
        std::cout
        << secondCascade->getAsString() << " BR = " << secondCascadeBrToEwino
        << std::endl << " => "
        << ( firstCascadeBrToEwino * secondCascadeBrToEwino )
        << std::endl
        << "subchannelCrossSectionTimesBrToEwinos = "
        << subchannelCrossSectionTimesBrToEwinos
        << std::endl
        << "subchannelZeroOrMoreJetsNLeptons = "
        << subchannelZeroOrMoreJetsNLeptons;
        std::cout << std::endl;
      }**/

      if( lhcFaserGlobal::negligibleBr < subchannelZeroOrMoreJetsNLeptons )
      {
        return ( subchannelCrossSectionTimesBrToEwinos
                 * subchannelZeroOrMoreJetsNLeptons
                 * signalDefinitions.getJetPlusMetAcceptance()->getAcceptance(
                                                                   initialPair,
                                                                  firstCascade,
                                                             secondCascade ) );
      }
      else
      {
        return 0.0;
      }
    }



    double const sameSignDilepton::defaultLeptonCut( 20.0 );
    // the default same-sign dilepton lepton transverse momentum cut is
    // 20.0 GeV.

    signalCalculator*
    sameSignDilepton::getCalculator( std::string const* const argumentString,
                                 signalDefinitionSet* const signalDefinitions )
    {
      signalDefinitions->setExcludedStandardModelProducts(
                          signalDefinitions->getShortcuts()->getInputShortcuts(
                                                          )->getNotInJets5() );
      signalCalculator* returnPointer( NULL );
      int characterSkipper;
      // this is used in interpreting the arguments string.
      if( 0 == argumentString->compare( 0,
                                        4,
                                        "7TeV" ) )
        //for a beam energy of 7 TeV...
      {
        signalDefinitions->setBeamEnergy( 7 );
        characterSkipper = 4;
      }
      else if( 0 == argumentString->compare( 0,
                                             5,
                                             "07TeV" ) )
        //for a beam energy of 7 TeV...
      {
        signalDefinitions->setBeamEnergy( 7 );
        characterSkipper = 5;
      }
      else if( 0 == argumentString->compare( 0,
                                             5,
                                             "10TeV" ) )
        //for a beam energy of 10 TeV...
      {
        signalDefinitions->setBeamEnergy( 10 );
        characterSkipper = 5;
      }
      else if( 0 == argumentString->compare( 0,
                                             5,
                                             "14TeV" ) )
        //for a beam energy of 14 TeV...
      {
        signalDefinitions->setBeamEnergy( 14 );
        characterSkipper = 5;
      }
      else
      {
        std::cout
        << std::endl
        << "LHC-FASER::error! signalHandler::signalHandler( "
        << argumentString << ", " <<  signalDefinitions << " ) was passed"
        << " a name it does not know (started with \"sameSignDilepton\","
        << " but then \"7TeV\", \"07TeV\", \"10TeV\" or \"14TeV\" are expected"
        << " as the next characters). its calculator is being set to return "
        << CppSLHA::CppSLHA_global::really_wrong_value_string << " for"
        << " every point.";
        returnPointer = new reallyWrongCalculator( signalDefinitions );
      }
      if( NULL == returnPointer )
        // if the signalCalculator hasn't been set to a reallyWrongCalculator...
      {
        if( 0 == argumentString->compare( characterSkipper,
                                          3,
                                          "pTl" ) )
          // if the next 3 characters denote a non-default lepton
          // transverse momentum cut...
        {
          double leptonCut;
          characterSkipper += 3;
          std::stringstream leptonCutStream( *argumentString );
          leptonCutStream.ignore( characterSkipper );
          leptonCutStream >> leptonCut;
          signalDefinitions->setLeptonCut( leptonCut );
        }
        else
        {
          signalDefinitions->setLeptonCut( defaultLeptonCut );
        }
        returnPointer = new sameSignDilepton( signalDefinitions );
      }
      return returnPointer;
    }

    sameSignDilepton::sameSignDilepton(
                               signalDefinitionSet* const signalDefinitions ) :
        signalCalculator( signalDefinitions )
    {
      excludedFinalStateParticles.push_back( CppSLHA::PDG_code::top );
      excludedFinalStateParticles.push_back( -(CppSLHA::PDG_code::top) );
    }

    sameSignDilepton::~sameSignDilepton()
    {
      // does nothing.
    }


    double
    sameSignDilepton::getOneLeptonFromEachCascade( fullCascade* firstCascade,
                                             bool const firstIsNotAntiparticle,
                                                   fullCascade* secondCascade,
                                           bool const secondIsNotAntiparticle )
    {
      return ( ( firstCascade->leptonFlavorSummedWithMinimumJets(
                                                        firstIsNotAntiparticle,
                                                            &signalDefinitions,
                                                                  0,
                                                                  1,
                                                                  0 )
                 * secondCascade->leptonFlavorSummedWithMinimumJets(
                                                       secondIsNotAntiparticle,
                                                            &signalDefinitions,
                                                                     0,
                                                                     1,
                                                                     0 ) )
                 + ( firstCascade->leptonFlavorSummedWithMinimumJets(
                                                        firstIsNotAntiparticle,
                                                            &signalDefinitions,
                                                                      0,
                                                                      0,
                                                                      1 )
                     * secondCascade->leptonFlavorSummedWithMinimumJets(
                                                       secondIsNotAntiparticle,
                                                            &signalDefinitions,
                                                                         0,
                                                                         0,
                                                                       1 ) ) );
    }

    double
    sameSignDilepton::valueForCurrentCascades( fullCascade* firstCascade,
                                               fullCascade* secondCascade )
    // see base version's description.
    {
      subchannelDileptonAcceptance
      = ( getDileptonFromSingleCascade( firstCascade,
                                        firstSparticleIsNotAntiparticle,
                                        secondCascade,
                                        secondSparticleIsNotAntiparticle )
          + getDileptonFromSingleCascade( secondCascade,
                                          secondSparticleIsNotAntiparticle,
                                          firstCascade,
                                          firstSparticleIsNotAntiparticle )
          + getOneLeptonFromEachCascade( firstCascade,
                                         firstSparticleIsNotAntiparticle,
                                         secondCascade,
                                         secondSparticleIsNotAntiparticle ) );
      /* chance for 1st cascade to produce the same-sign leptons * chance of
       * 2nd to not produce any further leptons (which would cause the event to
       * be rejected?) + chance with cascades swapped + chance for one from
       * each.
       */
      // debugging:
      /**if( 0.0005 < ( subchannelCrossSectionTimesBrToEwinos
                         * subchannelDileptonAcceptance ) )
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
        << "subchannelDileptonAcceptance = "
        << subchannelDileptonAcceptance;
        std::cout << std::endl;
      }**/
      return ( subchannelCrossSectionTimesBrToEwinos
               * subchannelDileptonAcceptance );
    }

  }  // end of signalCalculatorClasses namespace.

}  // end of LHC_FASER namespace.
