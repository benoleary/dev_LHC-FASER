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
    firstCascades( NULL ),
    secondCascades( NULL )
  {
    for( std::vector< signed_particle_shortcut_pair* >::const_iterator
         pairIterator = signalDefinitions->getShortcuts()->get_input_shortcuts(
                           )->get_sparticle_production_combinations()->begin();
         signalDefinitions->getShortcuts()->get_input_shortcuts(
                              )->get_sparticle_production_combinations()->end()
         > pairIterator;
         ++pairIterator )
    {
      channels.push_back( new productionChannelPointerSet( signalDefinitions,
                                                           *pairIterator ) );
    }
  }

  signalCalculator::~signalCalculator()
  {
    for( std::vector< productionChannelPointerSet* >::iterator
         deletionIterator = channels.begin();
         channels.end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
  }



  signalHandler::signalHandler( std::string const signalName,
                                signal_shortcuts* const shortcut,
                                double const crossSectionUnitFactor ) :
    readied_for_new_point( shortcut->get_input_shortcuts()->get_readier() ),
    signalName( signalName ),
    shortcut( shortcut ),
    crossSectionUnitFactor( crossSectionUnitFactor )
  {
    // debugging
    /**std::cout
    << std::endl
    << "debugging: signal_handler::signal_handler( " << given_signal_name
    << ", " << given_shortcuts << " ) called.";
    std::cout << std::endl;**/

    rateCalculator = NULL;
    // we start by ensuring that the pointer is NULL.

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
      = Atlas4jMET0l_calculator::getCalculator( &calculatorArgument,
                                                shortcut );
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
      = Atlas3jMET1l_calculator::getCalculator( &calculatorArgument,
                                                shortcut );
    }
    else if( 0 == signalName.compare( "CMS2jMETanyl14TeV" ) )
      // if the signal is the CMS 2-jet, missing transverse momentum,
      // any amount of leptons event rate, for a beam energy of 14 TeV...
    {
      std::cout
      << std::endl
      << "LHC-FASER::error! signal_handler::signal_handler( "
      << signalName << ", "<<  shortcut << " ) was passed"
      << " CMS2jMETanyl14TeV, but unfortunately it has not yet been"
      << " implemented. its calculator is being set to return only "
      << CppSLHA::CppSLHA_global::really_wrong_value_string
      << " for every point instead.";

      rateCalculator = new reallyWrongCalculator( shortcut );

    }
    else if( 0 == signalName.compare( "sigmaBreakdownTest" ) )
      // if the signal is the test to see if the cross-section breakdown works
      // as it should...
    {

      rateCalculator
      = new sigmaBreakdownTestCalculator( shortcut,
                                          shortcut->get_cross_section_handler(
                                                       )->get_table_set( 7 ) );

    }
    else
      // otherwise, a malformed signal name was passed:
    {
      std::cout
      << std::endl
      << "LHC-FASER::error! signal_handler::signal_handler( "
      << signalName << ", " <<  shortcut << " ) was passed a"
      << " name it does not know ( \"" << signalName
      << "\" ). its calculator is being set to return only "
      << CppSLHA::CppSLHA_global::really_wrong_value_string << " for"
      << " every point.";
      std::cout << std::endl;

      rateCalculator = new reallyWrongCalculator( shortcut );
    }

  }

  signalHandler::~signalHandler()
  {

    delete rateCalculator;

  }



  reallyWrongCalculator::reallyWrongCalculator(
                                           signal_shortcuts* const shortcut ) :
    signalCalculator( shortcut )
  {
    // just an initialization list.
  }

  reallyWrongCalculator::~reallyWrongCalculator()
  {
    // does nothing.
  }



  sigmaBreakdownTestCalculator::sigmaBreakdownTestCalculator(
                               signalDefinitionSet* const signalDefinitions ) :
    signalDefinitions( signalDefinitions )
  {
    this->signalDefinitions->setExcludedStandardModelProducts(
                        signalDefinitions->getShortcuts()->get_input_shortcuts(
                                                       )->get_not_in_jets5() );
  }

  sigmaBreakdownTestCalculator::~sigmaBreakdownTestCalculator()
  {

    // does nothing.

  }

  bool
  sigmaBreakdownTestCalculator::calculate( double* const signalValue,
                                           double* const uncertaintyFactor )
  /* this calculates the event rate for the signal & puts its value in
   * signal_value, & puts an estimate for the uncertainty into
   * uncertainty_factor, & returns true if it did all this successfully.
   */
  {
    // debugging:
    /**std::cout
    << std::endl
    << "debugging: sigmaBreakdownTestCalculator::calculate( "
    << signalValue << ", " << uncertaintyFactor << " ) called.";
    std::cout << std::endl;**/

    // start by setting the signal value to 0.0:
    *signalValue = 0.0;

    for( std::vector< productionChannelPointerSet* >::iterator
         channelIterator = channels.begin();
         channels.end() > channelIterator;
         ++channelIterator )
    {
      if( LHC_FASER_global::negligible_sigma
          < (*channelIterator)->getCrossSection() )
        // if it's worth looking at this channel...
      {
        std::cout
        << std::endl
        << "channel: "
        << *((*channelIterator)->getScoloredPair()->get_first_pointer(
                                                       )->get_name_or_antiname(
                                           (*channelIterator)->getScoloredPair(
                                              )->first_is_not_antiparticle() ))
        << " + "
        << *((*channelIterator)->getScoloredPair()->get_second_pointer(
                                                       )->get_name_or_antiname(
                                           (*channelIterator)->getScoloredPair(
                                            )->second_is_not_antiparticle() ));
        std::cout << std::endl;
        channelBrTotal = 0.0;

        // we have to look at all the open cascade pairings:
        firstCascades
        = (*channelIterator)->getFirstCascadeSet()->getOpenCascades(
          (*channelIterator)->getScoloredPair()->first_is_not_antiparticle(),
                                                          &signalDefinitions );
        secondCascades
        = (*channelIterator)->getSecondCascadeSet()->getOpenCascades(
           (*channelIterator)->getScoloredPair()->second_is_not_antiparticle(),
                                                           &signalDefinitions);

        for( std::vector< fullCascade* >::iterator
             firstCascadeIterator = firstCascades->begin();
             firstCascades->end() > firstCascadeIterator;
             ++firstCascadeIterator )
        {
          if( ( LHC_FASER_global::negligible_BR
              * 2.0 * (double)(firstCascades->size()) )
              < (*firstCascadeIterator)->getBrToEwino() )
            // this should mean we never throw away more than
            // 0.5 * LHC_FASER_global::negligible_BR of acceptance.
          {
            for( std::vector< fullCascade* >::iterator
                 secondCascadeIterator = secondCascades->begin();
                 secondCascades->end() > secondCascadeIterator;
                 ++secondCascadeIterator )
            {
              if( ( LHC_FASER_global::negligible_BR
                  * 2.0 * (double)(secondCascades->size()) )
                  < (*secondCascadeIterator)->getBrToEwino() )
                // this should mean we never throw away more than
                // 0.5 * LHC_FASER_global::negligible_BR of acceptance.
              {
                channelBrTotal += ( (*firstCascadeIterator)->getBrToEwino()
                                  * (*secondCascadeIterator)->getBrToEwino() );
              }
              // end of if the BR of the 2nd cascade to its electroweakino was
              // not negligible.
            }  // end of loop over 2nd cascade.
          }
          // end of if the BR of the 1st cascade to its electroweakino was not
          // negligible.
        }  // end of loop over 1st cascade.
        std::cout
        << "BR total: " << channelBrTotal;
        std::cout << std::endl;
      }  // end of if the channel's cross-section was not negligible.
    }  // end of loop over channels.
    return true;
  }



  signalCalculator*
  Atlas4jMET0l_calculator::getCalculator( std::string const* const arguments,
                                 signalDefinitionSet* const signalDefinitions )
  // this either returns a pointer to a new Atlas4jMET0l_calculator with cuts
  // decided by the given string, or a pointer to a reallyWrongCalculator.
  {
    signalDefinitions->setExcludedStandardModelProducts(
                        signalDefinitions->getShortcuts()->get_input_shortcuts(
                                                       )->get_not_in_jets5() );
    signalCalculator* returnPointer = NULL;
    int characterSkipper;
    // this is used in interpreting the arguments string.
    if( 0 == arguments->compare( 0,
                                 4,
                                 "7TeV" ) )
      //for a beam energy of 7 TeV...
    {
      signalDefinitions->setBeamEnergy( 7 );
      characterSkipper = 4;
    }
    else if( 0 == arguments->compare( 0,
                                      5,
                                      "07TeV" ) )
      //for a beam energy of 7 TeV...
    {
      signalDefinitions->setBeamEnergy( 7 );
      characterSkipper = 5;
    }
    else if( 0 == arguments->compare( 0,
                                      5,
                                      "10TeV" ) )
      //for a beam energy of 10 TeV...
    {
      signalDefinitions->setBeamEnergy( 10 );
      characterSkipper = 5;
    }
    else if( 0 == arguments->compare( 0,
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
      << arguments << ", " <<  signalDefinitions << " ) was passed"
      << " a name it does not know (started with \"Atlas4jMET0l\","
      << " but then \"7TeV\", \"07TeV\", \"10TeV\" or \"14TeV\" are expected"
      << " as the next characters). its calculator is being set to return "
      << CppSLHA::CppSLHA_global::really_wrong_value_string << " for"
      << " every point.";
      returnPointer = new reallyWrongCalculator();
    }
    if( NULL == returnPointer )
      // if the signalCalculator hasn't been set to a reallyWrongCalculator...
    {
      if( 0 == arguments->compare( characterSkipper,
                                   ( characterSkipper + 3 ),
                                   "pTl" ) )
        // if the next 3 characters denote a non-default lepton
        // transverse momentum cut...
      {
        double leptonCut;
        characterSkipper += 3;
        std::stringstream leptonCutStream( signalName );
        leptonCutStream.ignore( characterSkipper );
        leptonCutStream >> leptonCut;
        signalDefinitions->setLeptonCut( leptonCut );
      }
      else
      {
        signalDefinitions->setLeptonCut( defaultLeptonCut );
      }
      signalDefinitions->setJetCut( defaultExtraJetCut );
      returnPointer = new Atlas4jMET0l_calculator( signalDefinitions );
    }
    return returnPointer;
  }


  Atlas4jMET0l_calculator::Atlas4jMET0l_calculator(
                               signalDefinitionSet* const signalDefinitions ) :
    signalDefinitions( signalDefinitions )
  {
    std::string jetGridName( "Atlas4jMET" );
    fourJetKinematics
    = signalDefinitions->getShortcuts()->get_kinematics_handler()->get_table(
                                            signalDefinitions->getBeamEnergy(),
                                                                  &jetGridName,
                                                jetAcceptanceGridTableColumn );

    jetGridName.assign( "Atlas3jMET" );
    threeJetKinematics
    = signalDefinitions->getShortcuts()->get_kinematics_handler()->get_table(
                                            signalDefinitions->getBeamEnergy(),
                                                                  &jetGridName,
                                                jetAcceptanceGridTableColumn );
    excludedFinalStateParticles.push_back( CppSLHA::PDG_code::top );
    excludedFinalStateParticles.push_back( -(CppSLHA::PDG_code::top) );
  }

  Atlas4jMET0l_calculator::~Atlas4jMET0l_calculator()
  {
    // does nothing.
  }


  bool
  Atlas4jMET0l_calculator::calculate( double* const signalValue,
                                      double* const uncertaintyFactor )
  /* this calculates the event rate for the signal & puts its value in
   * signal_value, & puts an estimate for the uncertainty into
   * uncertainty_factor, & returns true if it did all this successfully.
   */
  {
    // debugging:
    /**std::cout
    << std::endl
    << "debugging: Atlas4jMET0l_calculator::calculate( "
    << signalValue << ", " << uncertaintyFactor << " ) called.";
    std::cout << std::endl;**/

    // start by setting the signal value to 0.0:
    *signalValue = 0.0;

    for( std::vector< productionChannelPointerSet* >::iterator
         channelIterator = channels.begin();
         channels.end() > channelIterator;
         ++channelIterator )
    {
      if( LHC_FASER_global::negligible_sigma
          < (*channelIterator)->getCrossSection() )
        // if it's worth looking at this channel...
      {
        // we have to look at all the open cascade pairings:
        firstCascades
        = (*channelIterator)->getFirstCascadeSet()->getOpenCascades(
            (*channelIterator)->getScoloredPair()->first_is_not_antiparticle(),
                                                          &signalDefinitions );
        secondCascades
        = (*channelIterator)->getSecondCascadeSet()->getOpenCascades(
           (*channelIterator)->getScoloredPair()->second_is_not_antiparticle(),
                                                          &signalDefinitions );

        for( std::vector< fullCascade* >::iterator
             firstCascadeIterator = firstCascades->begin();
             firstCascades->end() > firstCascadeIterator;
             ++firstCascadeIterator )
        {
          firstCascadeBrToEwino = (*firstCascadeIterator)->getBrToEwino(
                                                &excludedFinalStateParticles );
          if( ( LHC_FASER_global::negligible_BR
              * 2.0 * (double)(firstCascades->size()) )
              < firstCascadeBrToEwino )
            // this should mean we never throw away more than
            // 0.5 * LHC_FASER_global::negligible_BR of acceptance.
          {
            for( std::vector< fullCascade* >::iterator
                 secondCascadeIterator = secondCascades->begin();
                 secondCascades->end() > secondCascadeIterator;
                 ++secondCascadeIterator )
            {
              secondCascadeBrToEwino = (*secondCascadeIterator)->getBrToEwino(
                                                &excludedFinalStateParticles );
              if( ( LHC_FASER_global::negligible_BR
                  * 2.0 * (double)(secondCascades->size()) )
                  < secondCascadeBrToEwino )
                // this should mean we never throw away more than
                // 0.5 * LHC_FASER_global::negligible_BR of acceptance.
              {
                subchannelCrossSectionTimesBrToEwinos
                = ( (*channelIterator)->getCrossSection()
                    * firstCascadeBrToEwino * secondCascadeBrToEwino );
                subchannelZeroOrMoreJetsZeroLeptons
                = ( (*firstCascadeIterator
                           )->unspecifiedJetsSpecifiedChargeSummedLeptons( 0 )
                    * (*secondCascadeIterator
                         )->unspecifiedJetsSpecifiedChargeSummedLeptons( 0 ) );
                // chance of both cascades not producing any leptons which pass
                // the cut & thus cause the event to be rejected.
                subchannelOneOrMoreJetsZeroLeptons
                = ( subchannelZeroOrMoreJetsZeroLeptons
                    - ( (*firstCascadeIterator
                               )->specifiedJetsSpecifiedChargeSummedLeptons( 0,
                                                                            0 )
                        * (*secondCascadeIterator
                               )->specifiedJetsSpecifiedChargeSummedLeptons( 0,
                                                                       0 ) ) );
                // chance of no lepton rejection removing the portion without
                // any extra jet passing the cut.

                if( LHC_FASER_global::negligible_BR
                    < subchannelZeroOrMoreJetsZeroLeptons )
                {
                  fourJetAcceptance
                  = fourJetKinematics->getJetPlusMetAcceptance(
                                         (*channelIterator)->getScoloredPair(),
                                                         *firstCascadeIterator,
                                                      *secondCascadeIterator );
                  subchannelValue
                  = ( subchannelCrossSectionTimesBrToEwinos
                      * subchannelZeroOrMoreJetsZeroLeptons
                      * fourJetAcceptance );

                  // it's not going to ever be the case where the acceptance
                  // for 1+ jets is greater than 0+ jets...
                  if( LHC_FASER_global::negligible_BR
                      < subchannelOneOrMoreJetsZeroLeptons )
                  {
                    subchannelValue
                    += ( subchannelCrossSectionTimesBrToEwinos
                         * subchannelOneOrMoreJetsZeroLeptons
                         * ( threeJetKinematics->getJetPlusMetAcceptance(
                                         (*channelIterator)->getScoloredPair(),
                                                         *firstCascadeIterator,
                                                       *secondCascadeIterator )
                             - fourJetAcceptance ) );
                    // we assume that the 3-jet+MET acceptance from the grid is
                    // greater than the 4-jet+MET acceptance.
                  }
                  // here is where we update *signal_value to account for this
                  // particular subchannel:
                  *signalValue += subchannelValue;
                }
                // end of if the BR of the subchannel to zero leptons is not
                // negligible.
              }
              // end of if the BR of the 2nd cascade to its electroweakino was
              // not negligible.
            }  // end of loop over 2nd cascade.
          }
          // end of if the BR of the 1st cascade to its electroweakino was not
          // negligible.
        }  // end of loop over 1st cascade.
      }  // end of if the channel's cross-section was not negligible.
    }  // end of loop over channels.
    // currently, we do the uncertainty very roughly:
    *uncertaintyFactor
    = ( *signalValue * signalDefinitions->get_uncertainty() );
    return true;
  }



  Atlas3jMET1l_calculator::getCalculator( std::string const* const arguments,
                                 signalDefinitionSet* const signalDefinitions )
  {
    signalDefinitions->setExcludedStandardModelProducts(
                        signalDefinitions->getShortcuts()->get_input_shortcuts(
                                                       )->get_not_in_jets5() );
    signalCalculator* returnPointer = NULL;
    int characterSkipper;
    // this is used in interpreting the arguments string.
    if( 0 == arguments->compare( 0,
                                 4,
                                 "7TeV" ) )
      //for a beam energy of 7 TeV...
    {
      signalDefinitions->setBeamEnergy( 7 );
      characterSkipper = 4;
    }
    else if( 0 == arguments->compare( 0,
                                      5,
                                      "07TeV" ) )
      //for a beam energy of 7 TeV...
    {
      signalDefinitions->setBeamEnergy( 7 );
      characterSkipper = 5;
    }
    else if( 0 == arguments->compare( 0,
                                      5,
                                      "10TeV" ) )
      //for a beam energy of 10 TeV...
    {
      signalDefinitions->setBeamEnergy( 10 );
      characterSkipper = 5;
    }
    else if( 0 == arguments->compare( 0,
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
      << "LHC-FASER::error! signal_handler::signal_handler( "
      << arguments << ", " <<  signalDefinitions << " ) was passed"
      << " a name it does not know (started with \"Atlas3jMET1l\","
      << " but then \"7TeV\", \"07TeV\", \"10TeV\" or \"14TeV\" are expected"
      << " as the next characters). its calculator is being set to return "
      << CppSLHA::CppSLHA_global::really_wrong_value_string << " for"
      << " every point.";
      returnPointer = new reallyWrongCalculator();
    }
    if( NULL == returnPointer )
      // if the signalCalculator hasn't been set to a reallyWrongCalculator...
    {
      if( 0 == arguments->compare( characterSkipper,
                                   ( characterSkipper + 3 ),
                                   "pTl" ) )
        // if the next 3 characters denote a non-default primary lepton
        // transverse momentum cut...
      {
        // the argument should be
        // "[whatever]TeVpTl[primary lepton cut]pTl[secondary lepton cut]
        double leptonCut;
        characterSkipper += 3;
        std::string leptonCutSubstring( arguments,
                                        characterSkipper );
        // leptonCutSubstring is arguments with the "Atlas3j1l[whatever]TeVpTl"
        // bit removed.
        int failingLeptonCutPosition = leptonCutSubstring.find( "pTl" );
        // failingLeptonCutPosition is now the number of characters in
        // leptonCutSubstring before "pTl" (the 2nd "pTl" in arguments).
        std::stringstream
        leptonCutStream( leptonCutSubstring.substr( 0,
                                                  failingLeptonCutPosition ) );
        // the substring provided to leptonCutStream is just all the characters
        // in leptonCutSubstring before "pTl".
        leptonCutStream >> leptonCut;
        signalDefinitions->setPrimaryLeptonCut( leptonCut );
        // these characters are now translated into a number which is passed to
        // signalDefinitions as the primary lepton cut.
        leptonCutStream.str( leptonCutSubstring.substr(
                                          ( failingLeptonCutPosition + 3 ) ) );
        /* the substring now provided to leptonCutStream is all the characters
         * in leptonCutSubstring after "pTl", which should just be the
         * secondary lepton transverse momentum cut.
         */
        leptonCutStream >> leptonCut;
        signalDefinitions->setSecondaryLeptonCut( leptonCut );
      }
      else
      {
        signalDefinitions->setPrimaryLeptonCut( defaultPrimaryLeptonCut );
        signalDefinitions->setSecondaryLeptonCut( defaultSecondaryLeptonCut );
      }
      signalDefinitions->setJetCut( defaultExtraJetCut );
      returnPointer = new Atlas3jMET1l_calculator( signalDefinitions );
    }
    return returnPointer;
  }


  Atlas3jMET1l_calculator::Atlas3jMET1l_calculator(
                               signalDefinitionSet* const signalDefinitions ) :
    signalDefinitions( signalDefinitions )
  {
    std::string jetGridName( "Atlas3jMET" );
    this->signalDefinitions.setJetPlusMetAcceptance(
         signalDefinitions.getShortcuts()->get_kinematics_handler()->get_table(
                                            signalDefinitions->getBeamEnergy(),
                                                                  &jetGridName,
                                              jetAcceptanceGridTableColumn ) );
    excludedFinalStateParticles.push_back( CppSLHA::PDG_code::top );
    excludedFinalStateParticles.push_back( -(CppSLHA::PDG_code::top) );
  }

  Atlas3jMET1l_calculator::~Atlas3jMET1l_calculator()
  {
    // does nothing.
  }


  bool
  Atlas3jMET1l_calculator::calculate( double* const signalValue,
                                      double* const uncertaintyFactor )
  /* this calculates the event rate for the signal & puts its value in
   * signal_value, & puts an estimate for the uncertainty into
   * uncertainty_factor, & returns true if it did all this successfully.
   */
  {
    // debugging:
    /**std::cout
    << std::endl
    << "debugging: Atlas3jMET1l_calculator::calculate( "
    << signalValue << ", " << uncertaintyFactor << " ) called.";
    std::cout << std::endl;**/

    // start by setting the signal value to 0.0:
    *signalValue = 0.0;

    for( std::vector< productionChannelPointerSet* >::iterator
         channelIterator = channels.begin();
         channels.end() > channelIterator;
         ++channelIterator )
    {
      if( LHC_FASER_global::negligible_sigma
          < (*channelIterator)->getCrossSection() )
        // if it's worth looking at this channel...
      {
        // we have to look at all the open cascade pairings:
        firstCascades
        = (*channelIterator)->getFirstCascadeSet()->getOpenCascades(
            (*channelIterator)->getScoloredPair()->first_is_not_antiparticle(),
                                                          &signalDefinitions );
        secondCascades
        = (*channelIterator)->getSecondCascadeSet()->getOpenCascades(
           (*channelIterator)->getScoloredPair()->second_is_not_antiparticle(),
                                                           &signalDefinitions);

        for( std::vector< fullCascade* >::iterator
             firstCascadeIterator = firstCascades->begin();
             firstCascades->end() > firstCascadeIterator;
             ++firstCascadeIterator )
        {
          firstCascadeBrToEwino = (*firstCascadeIterator)->getBrToEwino(
                                                &excludedFinalStateParticles );
          if( ( LHC_FASER_global::negligible_BR
              * 2.0 * (double)(firstCascades->size()) )
              < firstCascadeBrToEwino )
            // this should mean we never throw away more than
            // 0.5 * LHC_FASER_global::negligible_BR of acceptance.
          {
            for( std::vector< fullCascade* >::iterator
                 secondCascadeIterator = secondCascades->begin();
                 secondCascades->end() > secondCascadeIterator;
                 ++secondCascadeIterator )
            {
              secondCascadeBrToEwino = (*secondCascadeIterator)->getBrToEwino(
                                                &excludedFinalStateParticles );
              if( ( LHC_FASER_global::negligible_BR
                  * 2.0 * (double)(secondCascades->size()) )
                  < secondCascadeBrToEwino )
                // this should mean we never throw away more than
                // 0.5 * LHC_FASER_global::negligible_BR of acceptance.
              {
                subchannelCrossSectionTimesBrToEwinos
                = ( (*channelIterator)->getCrossSection()
                    * firstCascadeBrToEwino * secondCascadeBrToEwino );
                subchannelZeroOrMoreJetsOneLepton
                = ( (*firstCascadeIterator
                            )->unspecifiedJetsSpecifiedChargeSummedLeptons( 1 )
                    * (*secondCascadeIterator
                            )->unspecifiedJetsSpecifiedChargeSummedLeptons( 0 )
                    + (*firstCascadeIterator
                            )->unspecifiedJetsSpecifiedChargeSummedLeptons( 0 )
                      * (*secondCascadeIterator
                         )->unspecifiedJetsSpecifiedChargeSummedLeptons( 1 ) );
                /* chance for 1st cascade to produce the lepton * chance
                 * of 2nd to not produce any further leptons (which would
                 * cause the event to be rejected) + chance of lepton from
                 * 2nd cascade * chance of no lepton from 1st cascade.
                 */

                if( LHC_FASER_global::negligible_BR
                    < subchannelZeroOrMoreJetsOneLepton )
                {
                  // here is where we update *signal_value to account for this
                  // particular subchannel:
                  *signalValue += ( subchannelCrossSectionTimesBrToEwinos
                                    * subchannelZeroOrMoreJetsOneLepton
                                   * signalDefinitions.getJetPlusMetAcceptance(
                                         (*channelIterator)->getScoloredPair(),
                                                         *firstCascadeIterator,
                                                    *secondCascadeIterator ) );
                }
                // end of if the BR of the subchannel to zero leptons is not
                // negligible.
              }
              // end of if the BR of the 2nd cascade to its electroweakino was
              // not negligible.
            }  // end of loop over 2nd cascade.
          }
          // end of if the BR of the 1st cascade to its electroweakino was not
          // negligible.
        }  // end of loop over 1st cascade.
      }  // end of if the channel's cross-section was not negligible.
    }  // end of loop over channels.
    // currently, we do the uncertainty very roughly:
    *uncertaintyFactor
    = ( *signalValue * signalDefinitions->get_uncertainty() );
    return true;
  }



}  // end of LHC_FASER namespace.
