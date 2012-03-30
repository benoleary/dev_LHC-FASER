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

  int
  signalCalculator::parseLeptonCutType( std::string const& argumentString )
  // this looks for "_0l", "_1l", or "_2l", & returns the appropriate int.
  // returns -1 if it could not find any.
  {
    if( std::string::npos != argumentString.find( "_0l" ) )
      // zero leptons...
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "signalCalculator::parseLeptonCutType returning 0 leptons.";
      std::cout << std::endl;**/

      return 0;
    }
    else if( std::string::npos != argumentString.find( "_1l" ) )
      // one lepton...
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "signalCalculator::parseLeptonCutType returning 1 lepton.";
      std::cout << std::endl;**/

      return 1;
    }
    else if( std::string::npos != argumentString.find( "_2l" ) )
      // two leptons
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "signalCalculator::parseLeptonCutType returning 2 leptons.";
      std::cout << std::endl;**/

      return 2;
    }
    else
    {
      std::cout
      << std::endl
      << "LHC-FASER::error! signalHandler::signalHandler( "
      << argumentString << " ) was passed a name it does not know (could not"
      << " find \"_0l\", \"_1l\", or \"_2l\", which are appropriate for this"
      << " type of signal). its calculator is being set to return "
      << CppSLHA::CppSLHA_global::really_wrong_value_string
      << " for every point.";
      return -1;
    }
  }

  bool
  signalCalculator::parseBeamEnergy( std::string const& argumentString,
                                 signalDefinitionSet* const signalDefinitions )
  // this looks for "_7TeV", "_07TeV", "_10TeV", or "_14TeV", & sets the beam
  // energy appropriately. returns false if it could not find any.
  {
    if( ( std::string::npos != argumentString.find( "_7TeV" ) )
        ||
        ( std::string::npos != argumentString.find( "_07TeV" ) ) )
      // for a beam energy of 7 TeV...
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "signalCalculator::parseBeamEnergy setting to 7 TeV.";
      std::cout << std::endl;**/

      signalDefinitions->setBeamEnergy( 7 );
      return true;
    }
    else if( std::string::npos != argumentString.find( "_10TeV" ) )
      // for a beam energy of 10 TeV...
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "signalCalculator::parseBeamEnergy setting to 10 TeV.";
      std::cout << std::endl;**/

      signalDefinitions->setBeamEnergy( 10 );
      return true;
    }
    else if( std::string::npos != argumentString.find( "_14TeV" ) )
      // for a beam energy of 14 TeV...
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
      std::cout
      << std::endl
      << "LHC-FASER::error! signalHandler::signalHandler( "
      << argumentString << ", " <<  signalDefinitions << " ) was passed"
      << " a name it does not know (could not find \"_7TeV\", \"_07TeV\","
      << " \"_10TeV\", or \"_14TeV\"). its calculator is being set to return "
      << CppSLHA::CppSLHA_global::really_wrong_value_string
      << " for every point.";
      return true;
    }
  }

  void
  signalCalculator::parseLeptonTransverseMomentumCuts(
                                             std::string const& argumentString,
                                  signalDefinitionSet* const signalDefinitions,
                                          double const defaultPrimaryLeptonCut,
                                       double const defaultSecondaryLeptonCut )
  /* this looks for "_pTl" then a double then "GeV", then a subsequent "_pTl",
   * a subsequent double, & a subsequent "GeV", interpretted as the primary &
   * secondary lepton cuts respectively.
   */
  {
    signalDefinitions->setPrimaryLeptonCut( defaultPrimaryLeptonCut );
    signalDefinitions->setSecondaryLeptonCut( defaultSecondaryLeptonCut );
    // set the cuts to the defaults in case we don't find any cuts in the
    // string.

    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "signalCalculator::parseLeptonTransverseMomentumCuts( "
    << argumentString << ", " << signalDefinitions << ", "
    << defaultPrimaryLeptonCut << ", " << defaultSecondaryLeptonCut
    << " ) called. no further output implies that no cuts were found.";
    std::cout << std::endl;**/

    size_t
    leptonCutStartPosition( argumentString.find( "_pTl" ) );
    if( std::string::npos != leptonCutStartPosition )
    {
      size_t leptonCutEndPosition( argumentString.find( "GeV" ) );
      if( std::string::npos != leptonCutEndPosition )
        // if there is a non-default primary lepton transverse momentum cut...
      {
        /* the argument should be
         * "_pTl[primary lepton cut]GeV_pTl[secondary lepton cut]GeV
         * with the 2nd cut optional.
         */
        leptonCutStartPosition += 4;
        // now leptonCutStartPosition has skipped over the "_pTl".
        double leptonCut;
        std::stringstream
        leptonCutStream( argumentString.substr( leptonCutStartPosition,
                         ( leptonCutEndPosition - leptonCutStartPosition ) ) );
        leptonCutStream >> leptonCut;
        signalDefinitions->setLeptonCut( leptonCut );
        // this also sets the secondary cut, in case we don't find another.

        // debugging:
        /**std::cout << std::endl << "debugging:"
        << std::endl
        << "primary leptonCut = " << leptonCut
        << ", argumentString.substr( " << leptonCutStartPosition
        << ", ( " << leptonCutEndPosition << " - " << leptonCutStartPosition
        << " ) ) = " << argumentString.substr( leptonCutStartPosition,
                           ( leptonCutEndPosition - leptonCutStartPosition ) );
        std::cout << std::endl;**/

        std::string
        possibleSecondaryCut( argumentString.substr( leptonCutEndPosition ) );
        leptonCutStartPosition = possibleSecondaryCut.find( "_pTl" );
        if( std::string::npos != leptonCutStartPosition )
        {
          leptonCutEndPosition = possibleSecondaryCut.find( "GeV" );
          if( std::string::npos != leptonCutEndPosition )
            // if there is a non-default secondary lepton transverse momentum
            // cut...
          {
            leptonCutStartPosition += 4;
            // now leptonCutStartPosition has skipped over the "_pTl".
            leptonCutStream.clear();
            leptonCutStream.str( possibleSecondaryCut.substr(
                                                        leptonCutStartPosition,
                         ( leptonCutEndPosition - leptonCutStartPosition ) ) );
            leptonCutStream >> leptonCut;
            signalDefinitions->setSecondaryLeptonCut( leptonCut );

            // debugging:
            /**std::cout << std::endl << "debugging:"
            << std::endl
            << "secondary leptonCut = " << leptonCut
            << ", possibleSecondaryCut.substr( " << leptonCutStartPosition
            << ", ( " << leptonCutEndPosition << " - "
            << leptonCutStartPosition << " ) ) = "
            << possibleSecondaryCut.substr( leptonCutStartPosition,
                           ( leptonCutEndPosition - leptonCutStartPosition ) );
            std::cout << std::endl;**/
          }
        }  // end of looking for a secondary cut.
      }
    }  // end of looking for a primary cut.
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
    std::vector< signalCalculatorCreator > creationFunctions;
    creationFunctions.push_back(
       &(signalClasses::atlasFourJetMetPlusGivenLeptonCuts::getCalculator()) );
    for( int whichFunction( (int)(creationFunctions.size()) - 1 );
         ( ( 0 <= whichFunction )
           &&
         ( NULL == rateCalculator ) );
         --whichFunction )
    {
      rateCalculator = (*(creationFunctions[ whichFunction ]))( signalName,
                                                signalPreparationDefinitions );
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

      rateCalculator = new signalClasses::reallyWrongCalculator(
                                                signalPreparationDefinitions );
    }

    below is just for reference while making getCalculator() for the derived
        classes:

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
      = signalClasses::atlasFourJetMetNLeptons::getCalculator(
                                                           &calculatorArgument,
                                                signalPreparationDefinitions );
    }
    else if( 0 == signalName.compare( 0,
                                      10,
                                      "Atlas3jMET" ) )
      // if the signal is the Atlas 3-jet, missing transverse momentum,
      // 1 lepton event rate...
    {
      calculatorArgument.assign( ( signalName.begin() + 10 ),
                                 signalName.end() );
      rateCalculator
      = signalClasses::atlasThreeJetMetNLeptons::getCalculator(
                                                           &calculatorArgument,
                                                signalPreparationDefinitions );
    }
    else if( 0 == signalName.compare( "CMS2jMET_anyl_14TeV" ) )
      // if the signal is the CMS 2-jet, missing transverse momentum,
      // any amount of leptons event rate, for a beam energy of 14 TeV...
    {
      std::cout
      << std::endl
      << "LHC-FASER::error! signalHandler::signalHandler( "
      << signalName << ", "<<  inputShortcut << " ) was passed"
      << " CMS2jMET_anyl_14TeV, but unfortunately it has not yet been"
      << " implemented. its calculator is being set to return only "
      << CppSLHA::CppSLHA_global::really_wrong_value_string
      << " for every point instead.";

      rateCalculator
      = new signalClasses::reallyWrongCalculator(
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
      = signalClasses::sameSignDilepton::getCalculator(
                                                           &calculatorArgument,
                                                signalPreparationDefinitions );
    }
    else if( 0 == signalName.compare( "sigmaBreakdownTest" ) )
      // if the signal is the test to see if the cross-section breakdown works
      // as it should...
    {
      signalPreparationDefinitions->setBeamEnergy( 7 );
      rateCalculator
      = new signalClasses::sigmaBreakdownTest(
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
      = new signalClasses::reallyWrongCalculator(
                                                signalPreparationDefinitions );
    }
  }

  signalHandler::~signalHandler()
  {
    delete rateCalculator;
  }



  namespace signalClasses
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
                                  signalDefinitionSet* const signalDefinitions,
                                        int const exactNumberOfJets,
                                        fullCascade* firstCascade,
                                        bool firstIsNotAntiparticle,
                                        fullCascade* secondCascade,
                                        bool secondIsNotAntiparticle ) const
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



      chargeAndFlavorSummed::chargeAndFlavorSummed(
                                                  int const numberOfLeptons ) :
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
                                  signalDefinitionSet* const signalDefinitions,
                                               int const exactNumberOfJets,
                                               fullCascade* firstCascade,
                                               bool firstIsNotAntiparticle,
                                               fullCascade* secondCascade,
                                           bool secondIsNotAntiparticle ) const
      {
        double returnValue( 0.0 );

        for( int jetsFromFirst( exactNumberOfJets );
             0 <= jetsFromFirst;
             --jetsFromFirst )
        {
          for( int leptonsFromFirst( numberOfLeptons );
               0 <= leptonsFromFirst;
               --leptonsFromFirst )
          {
            returnValue
            += ( firstCascade->specifiedJetsSpecifiedChargeSummedLeptons(
                                                             signalDefinitions,
                                                             exactNumberOfJets,
                                                             leptonsFromFirst )
                 * secondCascade->specifiedJetsSpecifiedChargeSummedLeptons(
                                                             signalDefinitions,
                                         ( exactNumberOfJets - jetsFromFirst ),
                                    ( numberOfLeptons - leptonsFromFirst ) ) );
          }
        }
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
                                  signalDefinitionSet* const signalDefinitions,
                                      int const exactNumberOfJets,
                                      fullCascade* firstCascade,
                                      bool firstIsNotAntiparticle,
                                      fullCascade* secondCascade,
                                      bool secondIsNotAntiparticle ) const
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
                                  signalDefinitionSet* const signalDefinitions,
                                      int const exactNumberOfJets,
                                      fullCascade* firstCascade,
                                      bool firstIsNotAntiparticle,
                                      fullCascade* secondCascade,
                                      bool secondIsNotAntiparticle ) const
      {
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
              returnValue
              += ( firstCascade->leptonChargeSummedWithSpecifiedJets(
                                                        firstIsNotAntiparticle,
                                                             signalDefinitions,
                                                             jetsFromFirst,
                                                      negativeLeptonsFromFirst,
                                                     positiveLeptonsFromFirst )
                   * secondCascade->leptonChargeSummedWithSpecifiedJets(
                                                       secondIsNotAntiparticle,
                                                             signalDefinitions,
                                         ( exactNumberOfJets - jetsFromFirst ),
                        ( numberOfNegativeLeptons - negativeLeptonsFromFirst ),
                    ( numberOfPositiveLeptons - positiveLeptonsFromFirst ) ) );
            }
          }
        }
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
                                  signalDefinitionSet* const signalDefinitions,
                                       int const exactNumberOfJets,
                                       fullCascade* firstCascade,
                                       bool firstIsNotAntiparticle,
                                       fullCascade* secondCascade,
                                       bool secondIsNotAntiparticle ) const
      {
        double returnValue( 0.0 );

        for( int jetsFromFirst( exactNumberOfJets );
             0 <= jetsFromFirst;
             --jetsFromFirst )
        {
          returnValue
          = ( ( firstCascade->specifiedJetsOneOssfMinusOsdfPair(
                                                             signalDefinitions,
                                                         exactNumberOfJets )
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

    }  // end of leptonAcceptanceStyle namespace



    leptonAcceptanceForPairFactory::leptonAcceptanceForPairFactory() :
        fullySpecifiedSet(),
        noLeptonOrExtraJetCutInstance(),
        chargeAndFlavorSummedSet(),
        chargeSummedSet(),
        flavorSummedSet(),
        ossfMinusOsdfInstance()
    {
      // just an initialization list.
    }

    leptonAcceptanceForPairFactory::~leptonAcceptanceForPairFactory()
    {
      for( int deletionIndex( fullySpecifiedSet.size() - 1 );
           0 <= deletionIndex;
           --deletionIndex )
      {
        delete fullySpecifiedSet[ deletionIndex ];
      }
      for( int deletionIndex( chargeAndFlavorSummedSet.size() - 1 );
           0 <= deletionIndex;
           --deletionIndex )
      {
        delete chargeAndFlavorSummedSet[ deletionIndex ];
      }
      for( int deletionIndex( chargeSummedSet.size() - 1 );
           0 <= deletionIndex;
           --deletionIndex )
      {
        delete chargeSummedSet[ deletionIndex ];
      }
      for( int deletionIndex( flavorSummedSet.size() - 1 );
           0 <= deletionIndex;
           --deletionIndex )
      {
        delete flavorSummedSet[ deletionIndex ];
      }
    }


    leptonAcceptanceForCascadePair const*
    leptonAcceptanceForPairFactory::getFullySpecified(
                                           int const numberOfNegativeElectrons,
                                           int const numberOfPositiveElectrons,
                                               int const numberOfNegativeMuons,
                                              int const numberOfPositiveMuons )
    {
      leptonAcceptanceForCascadePair* returnPointer( NULL );
      for( int searchIndex( fullySpecifiedSet.size() - 1 );
           0 <= searchIndex;
           --searchIndex )
      {
        if( fullySpecifiedSet[ searchIndex ]->isSameAs(
                                                     numberOfNegativeElectrons,
                                                     numberOfPositiveElectrons,
                                                        numberOfNegativeMuons,
                                                      numberOfPositiveMuons ) )
        {
          returnPointer = fullySpecifiedSet[ searchIndex ];
          searchIndex = -1;
          // stop looking.
        }
      }
      if( NULL == returnPointer )
      {
        fullySpecifiedSet.push_back( new leptonAcceptanceStyle::fullySpecified(
                                                     numberOfNegativeElectrons,
                                                     numberOfPositiveElectrons,
                                                         numberOfNegativeMuons,
                                                     numberOfPositiveMuons ) );
        returnPointer = fullySpecifiedSet.back();
      }
      return returnPointer;
    }

    leptonAcceptanceForCascadePair const*
    leptonAcceptanceForPairFactory::getChargeAndFlavorSummed(
                                                    int const numberOfLeptons )
    {
      leptonAcceptanceForCascadePair* returnPointer( NULL );
      for( int searchIndex( chargeAndFlavorSummedSet.size() - 1 );
           0 <= searchIndex;
           --searchIndex )
      {
        if( chargeAndFlavorSummedSet[ searchIndex ]->isSameAs(
                                                            numberOfLeptons ) )
        {
          returnPointer = chargeAndFlavorSummedSet[ searchIndex ];
          searchIndex = -1;
          // stop looking.
        }
      }
      if( NULL == returnPointer )
      {
        chargeAndFlavorSummedSet.push_back(
                              new leptonAcceptanceStyle::chargeAndFlavorSummed(
                                                           numberOfLeptons ) );
        returnPointer = chargeAndFlavorSummedSet.back();
      }
      return returnPointer;
    }

    leptonAcceptanceForCascadePair const*
    leptonAcceptanceForPairFactory::getChargeSummed(
                                                   int const numberOfElectrons,
                                                      int const numberOfMuons )
    {
      leptonAcceptanceForCascadePair* returnPointer( NULL );
      for( int searchIndex( chargeSummedSet.size() - 1 );
           0 <= searchIndex;
           --searchIndex )
      {
        if( chargeSummedSet[ searchIndex ]->isSameAs( numberOfElectrons,
                                                      numberOfMuons ) )
        {
          returnPointer = chargeSummedSet[ searchIndex ];
          searchIndex = -1;
          // stop looking.
        }
      }
      if( NULL == returnPointer )
      {
        chargeSummedSet.push_back( new leptonAcceptanceStyle::chargeSummed(
                                                             numberOfElectrons,
                                                             numberOfMuons ) );
        returnPointer = chargeSummedSet.back();
      }
      return returnPointer;
    }

    leptonAcceptanceForCascadePair const*
    leptonAcceptanceForPairFactory::getFlavorSummed(
                                             int const numberOfNegativeLeptons,
                                            int const numberOfPositiveLeptons )
    {
      leptonAcceptanceForCascadePair* returnPointer( NULL );
      for( int searchIndex( flavorSummedSet.size() - 1 );
           0 <= searchIndex;
           --searchIndex )
      {
        if( flavorSummedSet[ searchIndex ]->isSameAs( numberOfNegativeLeptons,
                                                    numberOfPositiveLeptons ) )
        {
          returnPointer = flavorSummedSet[ searchIndex ];
          searchIndex = -1;
          // stop looking.
        }
      }
      if( NULL == returnPointer )
      {
        flavorSummedSet.push_back( new leptonAcceptanceStyle::flavorSummed(
                                                       numberOfNegativeLeptons,
                                                   numberOfPositiveLeptons ) );
        returnPointer = flavorSummedSet.back();
      }
      return returnPointer;
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
    // (for N leptons to *pass*) & 10.0 GeV (for all others to *fail*).

    signalCalculator*
    atlasFourJetMetPlusGivenLeptonCuts::getCalculator(
                                       std::string const* const argumentString,
        leptonAcceptanceForCascadePair const* const leptonAcceptanceCalculator,
                                 signalDefinitionSet* const signalDefinitions )
    /* this either returns a pointer to a new atlasFourJetMetNLeptons with
     * cuts decided by the given string, or a pointer to a
     * reallyWrongCalculator.
     */
    {
      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "atlasFourJetMetPlusGivenLeptonCuts::getCalculator( "
      << *argumentString << ", ... ) called.";
      std::cout << std::endl;/**/

      signalDefinitions->setExcludedStandardModelProducts(
                          signalDefinitions->getShortcuts()->getInputShortcuts(
                                                          )->getNotInJets5() );
      signalCalculator* returnPointer( NULL );
      if( parseBeamEnergy( *argumentString,
                           signalDefinitions ) )
      {
        int const numberOfLeptons( 0 );
        double
        defaultPrimaryCutBasedOnLeptonNumber( defaultPrimaryLeptonCut );
        if( 0 == numberOfLeptons )
        {
          defaultPrimaryCutBasedOnLeptonNumber = defaultSecondaryLeptonCut;
        }
        parseLeptonTransverseMomentumCuts( *argumentString,
                                           signalDefinitions,
                                          defaultPrimaryCutBasedOnLeptonNumber,
                                           defaultSecondaryLeptonCut );
        signalDefinitions->setJetCut( defaultExtraJetCut );
        returnPointer
        = new atlasFourJetMetPlusGivenLeptonCuts( signalDefinitions,
                                                  leptonAcceptanceCalculator );
      }
      else
      {
        returnPointer = new reallyWrongCalculator( signalDefinitions );
      }
      return returnPointer;
    }

    atlasFourJetMetPlusGivenLeptonCuts::atlasFourJetMetPlusGivenLeptonCuts(
                               signalDefinitionSet* const signalDefinitions,
     leptonAcceptanceForCascadePair const* const leptonAcceptanceCalculator ) :
        signalCalculator( signalDefinitions ),
        leptonAcceptanceCalculator( leptonAcceptanceCalculator )
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "atlasFourJetMetNLeptons constructing with primary lepton cut = "
      << signalDefinitions->getPrimaryLeptonCut()
      << ", secondary lepton cut = "
      << signalDefinitions->getSecondaryLeptonCut() << ", numberOfLeptons = "
      << numberOfLeptons;
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
      subchannelZeroOrMoreJetsNLeptons = 0.0;
      subchannelOneOrMoreJetsNLeptons = 0.0;
      int const numberOfLeptons( 0 );
      for( int leptonsInFirstCascade( 0 );
           numberOfLeptons >= leptonsInFirstCascade;
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
        -= ( firstCascade->specifiedJetsSpecifiedChargeSummedLeptons(
                                                            &signalDefinitions,
                                                                      0,
                                                        leptonsInFirstCascade )
             * secondCascade->specifiedJetsSpecifiedChargeSummedLeptons(
                                                            &signalDefinitions,
                                                                         0,
                               ( numberOfLeptons - leptonsInFirstCascade ) ) );
        // 1st we sum up the overlap that needs to be taken away.
      }
      subchannelOneOrMoreJetsNLeptons += subchannelZeroOrMoreJetsNLeptons;
      // then we specify the value minus the overlap with the totals.

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



    int const
    atlasFourJetMetNLeptons::jetAcceptanceGridTableColumn( 10 );
    // this is dependent on the format of the grids.
    double const
    atlasFourJetMetNLeptons::defaultExtraJetCut( 40.0 );
    // this is the standard cut for the jets beyond the hardest cut for this
    // signal as implemented in this code.
    double const
    atlasFourJetMetNLeptons::defaultPrimaryLeptonCut( 20.0 );
    double const
    atlasFourJetMetNLeptons::defaultSecondaryLeptonCut( 10.0 );
    // the default Atlas lepton transverse momentum cuts are 20.0 GeV
    // (for N leptons to *pass*) & 10.0 GeV (for all others to *fail*).

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
      signalCalculator* returnPointer( NULL );
      int numberOfLeptons( parseLeptonCutType( *argumentString ) );
      if( -1 == numberOfLeptons )
      {
        returnPointer = new reallyWrongCalculator( signalDefinitions );
      }
      else
      {
        if( parseBeamEnergy( *argumentString,
                             signalDefinitions ) )
        {
          double
          defaultPrimaryCutBasedOnLeptonNumber( defaultPrimaryLeptonCut );
          if( 0 == numberOfLeptons )
          {
            defaultPrimaryCutBasedOnLeptonNumber = defaultSecondaryLeptonCut;
          }
          parseLeptonTransverseMomentumCuts( *argumentString,
                                             signalDefinitions,
                                          defaultPrimaryCutBasedOnLeptonNumber,
                                             defaultSecondaryLeptonCut );
          signalDefinitions->setJetCut( defaultExtraJetCut );
          returnPointer = new atlasFourJetMetNLeptons( signalDefinitions,
                                                       numberOfLeptons );
        }
        else
        {
          returnPointer = new reallyWrongCalculator( signalDefinitions );
        }
      }
      return returnPointer;
    }

    atlasFourJetMetNLeptons::atlasFourJetMetNLeptons(
                               signalDefinitionSet* const signalDefinitions,
                                                  int const numberOfLeptons ) :
        signalCalculator( signalDefinitions ),
        numberOfLeptons( numberOfLeptons )
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "atlasFourJetMetNLeptons constructing with primary lepton cut = "
      << signalDefinitions->getPrimaryLeptonCut()
      << ", secondary lepton cut = "
      << signalDefinitions->getSecondaryLeptonCut() << ", numberOfLeptons = "
      << numberOfLeptons;
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
           numberOfLeptons >= leptonsInFirstCascade;
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
        -= ( firstCascade->specifiedJetsSpecifiedChargeSummedLeptons(
                                                            &signalDefinitions,
                                                                      0,
                                                        leptonsInFirstCascade )
             * secondCascade->specifiedJetsSpecifiedChargeSummedLeptons(
                                                            &signalDefinitions,
                                                                         0,
                               ( numberOfLeptons - leptonsInFirstCascade ) ) );
        // 1st we sum up the overlap that needs to be taken away.
      }
      subchannelOneOrMoreJetsNLeptons += subchannelZeroOrMoreJetsNLeptons;
      // then we specify the value minus the overlap with the totals.

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
    // the default Atlas lepton transverse momentum cuts are 20.0 GeV
    // (for N leptons to *pass*) & 10.0 GeV (for all others to *fail*).

    signalCalculator*
    atlasThreeJetMetNLeptons::getCalculator(
                                       std::string const* const argumentString,
                                 signalDefinitionSet* const signalDefinitions )
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "atlasThreeJetMetNLeptons::getCalculator( " << *argumentString
      << ", ... ) called.";
      std::cout << std::endl;**/
      signalDefinitions->setExcludedStandardModelProducts(
                          signalDefinitions->getShortcuts()->getInputShortcuts(
                                                          )->getNotInJets5() );
      signalCalculator* returnPointer( NULL );
      int numberOfLeptons( parseLeptonCutType( *argumentString ) );
      if( -1 == numberOfLeptons )
      {
        returnPointer = new reallyWrongCalculator( signalDefinitions );
      }
      else
      {
        if( parseBeamEnergy( *argumentString,
                             signalDefinitions ) )
        {
          parseLeptonTransverseMomentumCuts( *argumentString,
                                             signalDefinitions,
                                             defaultPrimaryLeptonCut,
                                             defaultSecondaryLeptonCut );
          signalDefinitions->setJetCut( defaultExtraJetCut );
          returnPointer = new atlasThreeJetMetNLeptons( signalDefinitions,
                                                       numberOfLeptons );
        }
        else
        {
          returnPointer = new reallyWrongCalculator( signalDefinitions );
        }
      }
      return returnPointer;
    }

    atlasThreeJetMetNLeptons::atlasThreeJetMetNLeptons(
                                  signalDefinitionSet* const signalDefinitions,
                                                  int const numberOfLeptons ) :
        signalCalculator( signalDefinitions ),
        numberOfLeptons( numberOfLeptons )
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "atlasThreeJetMetNLeptons constructing with primary lepton cut = "
      << signalDefinitions->getPrimaryLeptonCut()
      << ", secondary lepton cut = "
      << signalDefinitions->getSecondaryLeptonCut() << ", numberOfLeptons = "
      << numberOfLeptons;
      std::cout << std::endl;**/

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

      subchannelZeroOrMoreJetsNLeptons = 0.0;
      for( int leptonsInFirstCascade( 0 );
           numberOfLeptons >= leptonsInFirstCascade;
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
      }
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
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "sameSignDilepton::getCalculator( " << *argumentString
      << ", ... ) called.";
      std::cout << std::endl;**/
      signalDefinitions->setExcludedStandardModelProducts(
                          signalDefinitions->getShortcuts()->getInputShortcuts(
                                                          )->getNotInJets5() );
      signalCalculator* returnPointer( NULL );
      if( parseBeamEnergy( *argumentString,
                           signalDefinitions ) )
      {
        parseLeptonTransverseMomentumCuts( *argumentString,
                                           signalDefinitions,
                                           defaultLeptonCut,
                                           defaultLeptonCut );
        returnPointer = new sameSignDilepton( signalDefinitions );
      }
      else
      {
        returnPointer = new reallyWrongCalculator( signalDefinitions );
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

  }  // end of signalClasses namespace.

}  // end of LHC_FASER namespace.
