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
    secondCascadeBrToEwino( CppSLHA::CppSLHA_global::really_wrong_value )
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
                                 12,
                                 "Atlas4jMET0l" ) )
      // if the signal is the Atlas 4-jet, missing transverse momentum,
      // 0 leptons event rate...
    {
      calculatorArgument.assign( ( signalName.begin()+12 ),
                                 signalName.end() );
      rateCalculator
      = signalCalculatorClasses::atlasFourJetMetZeroLepton::getCalculator(
                                                           &calculatorArgument,
                                                signalPreparationDefinitions );
    }
    else if( 0 == signalName.compare( 0,
                                      12,
                                      "Atlas3jMET1l" ) )
      // if the signal is the Atlas 3-jet, missing transverse momentum,
      // 1 lepton event rate...
    {
      calculatorArgument.assign( ( signalName.begin()+12 ),
                                 signalName.end() );
      rateCalculator
      = signalCalculatorClasses::atlasThreeJetMetOneLepton::getCalculator(
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
    /* this calculates the event rate for the signal & puts its value in
     * signalValue, & puts an estimate for the uncertainty into
     * uncertainty_factor, & returns true if it did all this successfully.
     */
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
                  /**/std::cout << std::endl << "debugging:"
                  << std::endl
                  << "adding" << std::endl
                  << (*firstCascadeIterator)->getAsString() << " BR = "
                  << firstCascadeBrToEwino << std::endl
                  << (*secondCascadeIterator)->getAsString() << " BR = "
                  << secondCascadeBrToEwino << std::endl << " => "
                  << ( firstCascadeBrToEwino * secondCascadeBrToEwino );
                  std::cout << std::endl;/**/

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
    atlasFourJetMetZeroLepton::jetAcceptanceGridTableColumn( 10 );
    // this is dependent on the format of the grids.
    double const
    atlasFourJetMetZeroLepton::defaultExtraJetCut( 40.0 );
    // this is the standard cut for the jets beyond the hardest cut for this
    // signal as implemented in this code.
    double const
    atlasFourJetMetZeroLepton::defaultLeptonCut( 10.0 );
    // the default Atlas4jMET0l lepton transverse momentum cut is 10.0 GeV.

    signalCalculator*
    atlasFourJetMetZeroLepton::getCalculator(
                                       std::string const* const argumentString,
                                 signalDefinitionSet* const signalDefinitions )
    /* this either returns a pointer to a new atlasFourJetMetZeroLepton with
     * cuts decided by the given string, or a pointer to a
     * reallyWrongCalculator.
     */
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "atlasFourJetMetZeroLepton::getCalculator( " << *argumentString
      << ", ... ) called.";
      std::cout << std::endl;**/

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
        << " a name it does not know (started with \"Atlas4jMET0l\","
        << " but then \"7TeV\", \"07TeV\", \"10TeV\" or \"14TeV\" are expected"
        << " as the next characters). its calculator is being set to return "
        << CppSLHA::CppSLHA_global::really_wrong_value_string << " for"
        << " every point.";
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
                                                       3 );
        std::cout << std::endl;**/

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
        signalDefinitions->setJetCut( defaultExtraJetCut );
        returnPointer = new atlasFourJetMetZeroLepton( signalDefinitions );
      }
      return returnPointer;
    }

    atlasFourJetMetZeroLepton::atlasFourJetMetZeroLepton(
        signalDefinitionSet* const signalDefinitions ) :
        signalCalculator( signalDefinitions )
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

    atlasFourJetMetZeroLepton::~atlasFourJetMetZeroLepton()
    {
      // does nothing.
    }


    bool
    atlasFourJetMetZeroLepton::calculateValue( double* const signalValue,
                                              double* const uncertaintyFactor )
    /* this calculates the event rate for the signal & puts its value in
     * signalValue, & puts an estimate for the uncertainty into
     * uncertaintyFactor, & returns true if it did all this successfully.
     */
    {
      // debugging:
      /**std::cout
      << std::endl
      << "debugging: atlasFourJetMetZeroLepton::calculateValue( "
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

          for( std::vector< fullCascade* >::iterator
               firstCascadeIterator( firstCascades->begin() );
               firstCascades->end() > firstCascadeIterator;
               ++firstCascadeIterator )
          {
            // debugging:
            /**std::cout << std::endl << "debugging:"
            << std::endl
            << "1st: "
            << (*firstCascadeIterator)->getAsString();
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
                << "2nd: "
                << (*secondCascadeIterator)->getAsString();
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
                  subchannelZeroOrMoreJetsZeroLeptons
                  = ( (*firstCascadeIterator
                                )->unspecifiedJetsSpecifiedChargeSummedLeptons(
                                                            &signalDefinitions,
                                                                            0 )
                       * (*secondCascadeIterator
                                )->unspecifiedJetsSpecifiedChargeSummedLeptons(
                                                            &signalDefinitions,
                                                                         0 ) );
                  // chance of both cascades not producing any leptons which
                  // pass the cut & thus cause the event to be rejected.

                  // debugging:
                  /**std::cout << std::endl << "debugging:"
                  << std::endl
                  << "subchannelZeroOrMoreJetsZeroLeptons = "
                  << subchannelZeroOrMoreJetsZeroLeptons;
                  std::cout << std::endl;**/

                  subchannelOneOrMoreJetsZeroLeptons
                  = ( subchannelZeroOrMoreJetsZeroLeptons
                      - ( (*firstCascadeIterator
                                  )->specifiedJetsSpecifiedChargeSummedLeptons(
                                                            &signalDefinitions,
                                                                             0,
                                                                            0 )
                          * (*secondCascadeIterator
                                  )->specifiedJetsSpecifiedChargeSummedLeptons(
                                                            &signalDefinitions,
                                                                             0,
                                                                       0 ) ) );
                  // chance of no lepton rejection removing the portion without
                  // any extra jet passing the cut.

                  // debugging:
                  /**std::cout << std::endl << "debugging:"
                  << std::endl
                  << "subchannelOneOrMoreJetsZeroLeptons = "
                  << subchannelOneOrMoreJetsZeroLeptons;
                  std::cout << std::endl;**/

                  if( lhcFaserGlobal::negligibleBr
                      < subchannelZeroOrMoreJetsZeroLeptons )
                  {
                    fourJetAcceptance = fourJetKinematics->getAcceptance(
                                         (*channelIterator)->getScoloredPair(),
                                                         *firstCascadeIterator,
                                                      *secondCascadeIterator );
                    subchannelValue
                    = ( subchannelCrossSectionTimesBrToEwinos
                        * subchannelZeroOrMoreJetsZeroLeptons
                        * fourJetAcceptance );

                    // it's not going to ever be the case where the acceptance
                    // for 1+ jets is greater than 0+ jets...
                    if( lhcFaserGlobal::negligibleBr
                        < subchannelOneOrMoreJetsZeroLeptons )
                    {
                      subchannelValue
                      += ( subchannelCrossSectionTimesBrToEwinos
                          * subchannelOneOrMoreJetsZeroLeptons
                          * ( threeJetKinematics->getAcceptance(
                                         (*channelIterator)->getScoloredPair(),
                                                         *firstCascadeIterator,
                                                       *secondCascadeIterator )
                              - fourJetAcceptance ) );
                      // we assume that the 3-jet+MET acceptance from the grid
                      // is greater than the 4-jet+MET acceptance.
                    }
                    // here is where we update *signalValue to account for
                    // this particular subchannel:
                    *signalValue += subchannelValue;
                    // debugging:
                    /**std::cout << std::endl << "debugging:"
                    << std::endl
                    << "subchannelCrossSectionTimesBrToEwinos = "
                    << subchannelCrossSectionTimesBrToEwinos
                    << std::endl
                    << "subchannelZeroOrMoreJetsZeroLeptons = "
                    << subchannelZeroOrMoreJetsZeroLeptons
                    << std::endl
                    << "fourJetAcceptance = " << fourJetAcceptance
                    << std::endl
                    << "=> adding subchannelValue = " << subchannelValue;
                    std::cout << std::endl;**/
                    if( !( ( 1.0 >= subchannelZeroOrMoreJetsZeroLeptons )
                           &&
                           ( 0.0 <= subchannelZeroOrMoreJetsZeroLeptons )
                           &&
                           ( 1.0 >= subchannelOneOrMoreJetsZeroLeptons )
                           &&
                           ( 0.0 <= subchannelOneOrMoreJetsZeroLeptons )
                           &&
                           ( 1.0 >= fourJetAcceptance )
                           &&
                           ( 0.0 <= fourJetAcceptance )
                           &&
                           ( 0.0 <= subchannelValue ) ) )
                    {
                      std::cout << std::endl << "LHC-FASER::error!"
                      << std::endl
                      << *((*channelIterator)->getScoloredPair(
                                   )->getFirstParticle()->get_name_or_antiname(
                                           (*channelIterator)->getScoloredPair(
                                                 )->firstIsNotAntiparticle() ))
                      << " + "
                      << *((*channelIterator)->getScoloredPair(
                                  )->getSecondParticle()->get_name_or_antiname(
                                           (*channelIterator)->getScoloredPair(
                                               )->secondIsNotAntiparticle() ));
                      std::cout
                      << std::endl
                      << "1st: "
                      << (*firstCascadeIterator)->getAsString()
                      << "; 0+j, 0l = "
                      << (*firstCascadeIterator
                                )->unspecifiedJetsSpecifiedChargeSummedLeptons(
                                                            &signalDefinitions,
                                                                           0 );
                      std::cout
                      << std::endl
                      << "2nd: "
                      << (*secondCascadeIterator)->getAsString()
                      << "; 0+j, 0l = "
                      << (*secondCascadeIterator
                                )->unspecifiedJetsSpecifiedChargeSummedLeptons(
                                                            &signalDefinitions,
                                                                           0 );
                      std::cout << std::endl;
                      std::cout
                      << std::endl
                      << "produced subchannelZeroOrMoreJetsZeroLeptons = "
                      << subchannelZeroOrMoreJetsZeroLeptons
                      << ", subchannelOneOrMoreJetsZeroLeptons = "
                      << subchannelOneOrMoreJetsZeroLeptons
                      << ", & fourJetAcceptance = "
                      << fourJetAcceptance
                      << ", all of which should be >= 0.0 & <= 1.0 for valid"
                      << " input! subchannelValue = " << subchannelValue
                      << ", which should be >= 0.0 for valid input!";
                      std::cout << std::endl;
                    }
                  }
                  // end of if the BR of the subchannel to zero leptons is not
                  // negligible.
                }
                // end of if the BR of the 2nd cascade to its electroweakino
                // was not negligible.
              }  // end of loop over 2nd cascade.
            }
            // end of if the BR of the 1st cascade to its electroweakino was
            // not negligible.
          }  // end of loop over 1st cascade.
        }  // end of if the channel's cross-section was not negligible.
      }  // end of loop over productionChannels.
      // currently, we do the uncertainty very roughly:
      *uncertaintyFactor
      = ( *signalValue * signalDefinitions.getShortcuts()->getUncertainty() );
      return true;
    }



    int const atlasThreeJetMetOneLepton::jetAcceptanceGridTableColumn( 10 );
    // this is dependent on the format of the grids.
    double const atlasThreeJetMetOneLepton::defaultExtraJetCut( 40.0 );
    // this is the standard cut for the jets beyond the hardest cut for this
    // signal as implemented in this code.
    double const atlasThreeJetMetOneLepton::defaultPrimaryLeptonCut( 20.0 );
    double const atlasThreeJetMetOneLepton::defaultSecondaryLeptonCut( 10.0 );
    // the default Atlas4jMET0l lepton transverse momentum cuts are 20.0 GeV
    // (for a single lepton to *pass*) & 10.0 GeV (for all others to *fail*).

    signalCalculator*
    atlasThreeJetMetOneLepton::getCalculator(
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
        << " a name it does not know (started with \"Atlas3jMET1l\","
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
          // if the next 3 characters denote a non-default primary lepton
          // transverse momentum cut...
        {
          // the argument should be
          // "[whatever]TeVpTl[primary lepton cut]pTl[secondary lepton cut]
          double leptonCut;
          characterSkipper += 3;
          std::string leptonCutSubstring( *argumentString,
                                          characterSkipper );
          // leptonCutSubstring is arguments with the
          // "Atlas3j1l[whatever]TeVpTl" bit removed.
          int failingLeptonCutPosition( leptonCutSubstring.find( "pTl" ) );
          // failingLeptonCutPosition is now the number of characters in
          // leptonCutSubstring before "pTl" (the 2nd "pTl" in arguments).
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
          signalDefinitions->setSecondaryLeptonCut( leptonCut );
        }
        else
        {
          signalDefinitions->setPrimaryLeptonCut( defaultPrimaryLeptonCut );
          signalDefinitions->setSecondaryLeptonCut( defaultSecondaryLeptonCut );
        }
        signalDefinitions->setJetCut( defaultExtraJetCut );
        returnPointer = new atlasThreeJetMetOneLepton( signalDefinitions );
      }
      return returnPointer;
    }

    atlasThreeJetMetOneLepton::atlasThreeJetMetOneLepton(
                               signalDefinitionSet* const signalDefinitions ) :
        signalCalculator( signalDefinitions )
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

    atlasThreeJetMetOneLepton::~atlasThreeJetMetOneLepton()
    {
      // does nothing.
    }


    bool
    atlasThreeJetMetOneLepton::calculateValue( double* const signalValue,
                                              double* const uncertaintyFactor )
    /* this calculates the event rate for the signal & puts its value in
     * signalValue, & puts an estimate for the uncertainty into
     * uncertaintyFactor, & returns true if it did all this successfully.
     */
    {
      // debugging:
      /**std::cout
      << std::endl
      << "debugging: atlasThreeJetMetOneLepton::calculateValue( "
      << signalValue << ", " << uncertaintyFactor << " ) called."
      << " signalDefinitions.getPrimaryLeptonCut() = "
      << signalDefinitions.getPrimaryLeptonCut()
      << ", signalDefinitions.getSecondaryLeptonCut() = "
      << signalDefinitions.getSecondaryLeptonCut();
      std::cout << std::endl;**/

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
          // we have to look at all the open cascade pairings:
          firstCascades
          = (*channelIterator)->getFirstCascadeSet()->getOpenCascades();
          secondCascades
          = (*channelIterator)->getSecondCascadeSet()->getOpenCascades();

          for( std::vector< fullCascade* >::iterator
               firstCascadeIterator( firstCascades->begin() );
               firstCascades->end() > firstCascadeIterator;
               ++firstCascadeIterator )
          {
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
                  subchannelZeroOrMoreJetsOneLepton
                  = ( (*firstCascadeIterator
                                )->unspecifiedJetsSpecifiedChargeSummedLeptons(
                                                            &signalDefinitions,
                                                                            1 )
                      * (*secondCascadeIterator
                                )->unspecifiedJetsSpecifiedChargeSummedLeptons(
                                                            &signalDefinitions,
                                                                            0 )
                      + (*firstCascadeIterator
                                )->unspecifiedJetsSpecifiedChargeSummedLeptons(
                                                            &signalDefinitions,
                                                                            0 )
                        * (*secondCascadeIterator
                                )->unspecifiedJetsSpecifiedChargeSummedLeptons(
                                                            &signalDefinitions,
                                                                         1 ) );
                  /* chance for 1st cascade to produce the lepton * chance
                   * of 2nd to not produce any further leptons (which would
                   * cause the event to be rejected) + chance of lepton from
                   * 2nd cascade * chance of no lepton from 1st cascade.
                   */

                  if( lhcFaserGlobal::negligibleBr
                      < subchannelZeroOrMoreJetsOneLepton )
                  {
                    // here is where we update *signalValue to account for
                    // this particular subchannel:
                    *signalValue
                    += ( subchannelCrossSectionTimesBrToEwinos
                        * subchannelZeroOrMoreJetsOneLepton
                        * signalDefinitions.getJetPlusMetAcceptance(
                                                              )->getAcceptance(
                                         (*channelIterator)->getScoloredPair(),
                                                         *firstCascadeIterator,
                                                    *secondCascadeIterator ) );
                  }
                  // end of if the BR of the subchannel to zero leptons is not
                  // negligible.
                }
                // end of if the BR of the 2nd cascade to its electroweakino
                // was not negligible.
              }  // end of loop over 2nd cascade.
            }
            // end of if the BR of the 1st cascade to its electroweakino was
            // not negligible.
          }  // end of loop over 1st cascade.
        }  // end of if the channel's cross-section was not negligible.
      }  // end of loop over productionChannels.
      // currently, we do the uncertainty very roughly:
      *uncertaintyFactor
      = ( *signalValue * signalDefinitions.getShortcuts()->getUncertainty() );
      return true;
    }
  }  // end of signalCalculatorClasses namespace.

}  // end of LHC_FASER namespace.
