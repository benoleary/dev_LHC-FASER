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

  bool
  signalCalculator::parseBeamEnergy( std::string const& argumentString,
                                 signalDefinitionSet* const signalDefinitions,
                                     std::string& argumentRemainder )
  /* this looks for "_7TeV_", "_07TeV_", "_10TeV_", or "_14TeV_", & sets the
   * beam energy appropriately, putting argumentString from the 2nd '_' onwards
   * into argumentRemainder. returns false if it could not find any (& does not
   * modify argumentRemainder in this case).
   */
  {
    if( 0 == argumentString.compare( 0,
                                     6,
                                     "_7TeV_" ) )
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "signalCalculator::parseBeamEnergy setting to 7 TeV.";
      std::cout << std::endl;**/

      signalDefinitions->setBeamEnergy( 7 );
      argumentRemainder.assign( argumentString,
                                5 );
      return true;
    }
    else if( 0 == argumentString.compare( 0,
                                          7,
                                          "_07TeV_" ) )
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "signalCalculator::parseBeamEnergy setting to 07 TeV.";
      std::cout << std::endl;**/

      signalDefinitions->setBeamEnergy( 7 );
      argumentRemainder.assign( argumentString,
                                6 );
      return true;
    }
    else if( 0 == argumentString.compare( 0,
                                          7,
                                          "_10TeV_" ) )
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "signalCalculator::parseBeamEnergy setting to 10 TeV.";
      std::cout << std::endl;**/

      signalDefinitions->setBeamEnergy( 10 );
      argumentRemainder.assign( argumentString,
                                6 );
      return true;
    }
    else if( 0 == argumentString.compare( 0,
                                          7,
                                          "_14TeV_" ) )
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "signalCalculator::parseBeamEnergy setting to 14 TeV.";
      std::cout << std::endl;**/

      signalDefinitions->setBeamEnergy( 14 );
      argumentRemainder.assign( argumentString,
                                6 );
      return true;
    }
    else
    {
      return false;
    }
  }

  signalClasses::leptonAcceptanceForCascadePair*
  signalCalculator::parseLeptonAcceptance( std::string const& argumentString,
                                 signalDefinitionSet* const signalDefinitions )
  /* this looks for strings encoding the type of lepton cuts to use. the
   * strings are, where # stands for any string representing an integer, in the
   * order in which they are checked:
   * "_ossf" : ossfMinusOsdf
   * "_noExtraCut" : noLeptonCutNorExtraJetCut
   * "_#l" : chargeAndFlavorSummed
   * "_#lm#lp" : flavorSummed
   * "_#epm#mpm" : chargeSummed
   * "_#em#ep#mm#mp" : fullySpecified
   * parseLeptonTransverseMomentumCuts is then called on the remainder of
   * argumentString. finally, it creates a new leptonAcceptanceForCascadePair
   * & returns a pointer to it. NULL is returned if argumentString could not
   * be properly interpretted.
   */
  {
    signalClasses::leptonAcceptanceForCascadePair* returnPointer( NULL );

    // 1st look for the easy cases:
    if( 0 == argumentString.compare( 0,
                                     5,
                                     "_ossf" ) )
    {
      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "signalCalculator::parseLeptonAcceptance found \"_ossf\".";
      std::cout << std::endl;/**/
      returnPointer
      = new signalClasses::leptonAcceptanceStyle::ossfMinusOsdf();
    }
    else if( 0 == argumentString.compare( 0,
                                          11,
                                          "_noExtraCut" ) )
    {
      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "signalCalculator::parseLeptonAcceptance found \"_noExtraCut\".";
      std::cout << std::endl;/**/
      returnPointer
      = new signalClasses::leptonAcceptanceStyle::ossfMinusOsdf();
      argumentRemainder.assign( argumentString,
                                11 );
    }
    else if( 0 == argumentString.compare( 0,
                                          1,
                                          "_" ) )
      // otherwise look for a string of numeric characters...
    {
      size_t charPosition( 1 );
      std::string integerString( "" );
      int leptonNumber( parseOutLeptonNumber( argumentString,
                                              charPosition ) );
      // ... followed by...
      if( 'l' == argumentString[ charPosition ] )
      {
        if( ( argumentString.size() == charPosition )
            ||
            ( '_' == argumentString[ ( charPosition + 1 ) ] ) )
          // if argumentString ends with "_#l" or continues with "_", it's
          // chargeAndFlavorSummed.
        {
          // debugging:
          /**/std::cout << std::endl << "debugging:"
          << std::endl
          << "signalCalculator::parseLeptonAcceptance found \"_"
          << leptonNumber << "l\"...";
          std::cout << std::endl;/**/

          returnPointer
          = new signalClasses::leptonAcceptanceStyle::chargeAndFlavorSummed(
                                                                leptonNumber );
        }
        else if( 'm' == argumentString[ ( charPosition + 1 ) ] )
        {
          charPosition += 2;
          int secondLeptonNumber( parseOutLeptonNumber( argumentString,
                                                        charPosition ) );
          if( 0 == argumentString.compare( charPosition,
                                           2,
                                           "lp" ) )
          {
            // debugging:
            /**/std::cout << std::endl << "debugging:"
            << std::endl
            << "signalCalculator::parseLeptonAcceptance found \"_"
            << leptonNumber << "lm"
            << secondLeptonNumber << "lp\"...";
            std::cout << std::endl;/**/

            returnPointer
            = new signalClasses::leptonAcceptanceStyle::flavorSummed(
                                                                  leptonNumber,
                                                          secondLeptonNumber );
          }
        }
      }
      else if( 'e' == argumentString[ charPosition ] )
      {
        if( 0 == argumentString.compare( charPosition,
                                         3,
                                         "epm" ) )
        {
          charPosition += 3;
          int secondLeptonNumber( parseOutLeptonNumber( argumentString,
                                                        charPosition ) );
          if( 0 == argumentString.compare( charPosition,
                                           3,
                                           "mpm" ) )
          {
            // debugging:
            /**/std::cout << std::endl << "debugging:"
            << std::endl
            << "signalCalculator::parseLeptonAcceptance found \"_"
            << leptonNumber << "epm"
            << secondLeptonNumber << "mpm\"...";
            std::cout << std::endl;/**/

            returnPointer
            = new signalClasses::leptonAcceptanceStyle::chargeSummed(
                                                                  leptonNumber,
                                                          secondLeptonNumber );
          }
        }
        else if( 0 == argumentString.compare( charPosition,
                                              2,
                                              "em" ) )
        {
          charPosition += 2;
          int secondLeptonNumber( parseOutLeptonNumber( argumentString,
                                                        charPosition ) );
          if( 0 == argumentString.compare( charPosition,
                                           2,
                                           "ep" ) )
          {
            charPosition += 2;
            int thirdLeptonNumber( parseOutLeptonNumber( argumentString,
                                                         charPosition ) );
            if( 0 == argumentString.compare( charPosition,
                                             2,
                                             "mm" ) )
            {
              charPosition += 2;
              int fourthLeptonNumber( parseOutLeptonNumber( argumentString,
                                                            charPosition ) );
              if( 0 == argumentString.compare( charPosition,
                                               2,
                                               "mm" ) )
              {
                // debugging:
                /**/std::cout << std::endl << "debugging:"
                << std::endl
                << "signalCalculator::parseLeptonAcceptance found \"_"
                << leptonNumber << "em"
                << secondLeptonNumber << "ep"
                << thirdLeptonNumber << "mm"
                << fourthLeptonNumber << "mp\"...";
                std::cout << std::endl;/**/

                returnPointer
                = new signalClasses::leptonAcceptanceStyle::fullySpecified(
                                                                  leptonNumber,
                                                            secondLeptonNumber,
                                                             thirdLeptonNumber,
                                                          fourthLeptonNumber );
              }
            }
          }
        }
      }
    }

    if( NULL != returnPointer )
    {
      parseLeptonTransverseMomentumCuts( argumentString,
                                         signalDefinitions );
      /* this takes care of if there were lepton transverse momentum cuts
       * specified. it's OK to only set them after the new
       * leptonAcceptanceForCascadePair has been constructed because
       * signalDefinitions is still a pointer to the signalDefinitionSet that
       * returnPointer uses.
       */
    }
    return returnPointer;
  }

  void
  signalCalculator::parseLeptonTransverseMomentumCuts(
                                             std::string const& argumentString,
                                 signalDefinitionSet* const signalDefinitions )
  /* this looks for "_pTl" then a double then "GeV", then a subsequent "_pTl",
   * a subsequent double, & a subsequent "GeV", interpretted as the primary &
   * secondary lepton cuts respectively.
   */
  {
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "signalCalculator::parseLeptonTransverseMomentumCuts( "
    << argumentString << ", " << signalDefinitions
    << " ) called. no further output implies that no cuts were found.";
    std::cout << std::endl;/**/

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
        /**/std::cout << std::endl << "debugging:"
        << std::endl
        << "primary leptonCut = " << leptonCut
        << ", argumentString.substr( " << leptonCutStartPosition
        << ", ( " << leptonCutEndPosition << " - " << leptonCutStartPosition
        << " ) ) = " << argumentString.substr( leptonCutStartPosition,
                           ( leptonCutEndPosition - leptonCutStartPosition ) );
        std::cout << std::endl;/**/

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
            /**/std::cout << std::endl << "debugging:"
            << std::endl
            << "secondary leptonCut = " << leptonCut
            << ", possibleSecondaryCut.substr( " << leptonCutStartPosition
            << ", ( " << leptonCutEndPosition << " - "
            << leptonCutStartPosition << " ) ) = "
            << possibleSecondaryCut.substr( leptonCutStartPosition,
                           ( leptonCutEndPosition - leptonCutStartPosition ) );
            std::cout << std::endl;/**/
          }
        }  // end of looking for a secondary cut.
      }
    }  // end of looking for a primary cut.
  }



  signalHandler::signalHandler( std::string const signalName,
                                double const crossSectionUnitFactor,
                                signalShortcuts const* const inputShortcut ) :
    getsReadiedForNewPoint( signalPreparationDefinitions->getShortcuts(
                                        )->getInputShortcuts()->getReadier() ),
    signalName( signalName ),
    rateCalculator( NULL ),
    inputShortcut( inputShortcut ),
    signalPreparationDefinitions( inputShortcut ),
    signalValue( CppSLHA::CppSLHA_global::really_wrong_value ),
    uncertaintyFactor( CppSLHA::CppSLHA_global::really_wrong_value ),
    crossSectionUnitFactor( crossSectionUnitFactor )
  {
    std::vector< signalCalculatorCreator > creationFunctions;
    creationFunctions.push_back(
       &(signalClasses::atlasFourJetMetPlusGivenLeptonCuts::getCalculator()) );
    creationFunctions.push_back(
      &(signalClasses::atlasThreeJetMetPlusGivenLeptonCuts::getCalculator()) );
    creationFunctions.push_back(
               &(signalClasses::noJetMetButGivenLeptonCuts::getCalculator()) );
    creationFunctions.push_back(
                       &(signalClasses::sigmaBreakdownTest::getCalculator()) );
    for( int whichFunction( (int)(creationFunctions.size()) - 1 );
         ( ( 0 <= whichFunction )
           &&
           ( NULL == rateCalculator ) );
         --whichFunction )
    {
      rateCalculator = (*(creationFunctions[ whichFunction ]))( signalName,
                                               &signalPreparationDefinitions );
    }
    if( NULL == rateCalculator )
      // if none of the signal calculator creation functions returned a new
      // signalCalculator, a malformed signal name was passed:
    {
      rateCalculator = new signalClasses::reallyWrongCalculator( signalName,
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
                                                 std::string const& signalName,
                               signalDefinitionSet* const signalDefinitions ) :
        signalCalculator( signalDefinitions )
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
                                             std::string const& argumentString,
                                 signalDefinitionSet* const signalDefinitions )
    /* this either returns a pointer to a new
     * atlasFourJetMetPlusGivenLeptonCuts with cuts decided by the given
     * string, or NULL if the string could not be parsed properly.
     */
    {
      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "atlasFourJetMetPlusGivenLeptonCuts::getCalculator( "
      << argumentString << ", ... ) called.";
      std::cout << std::endl;/**/

      signalCalculator* returnPointer( NULL );

      if( 0 == argumentString.compare( 0,
                                       11,
                                       "Atlas4jMET_" ) )
      {
        std::string validArguments( argumentString.substr( 10 ) );
        // debugging:
        /**/std::cout << std::endl << "debugging:"
        << std::endl
        << "validArguments = \"" << validArguments << "\"";
        std::cout << std::endl;/**/

        // the lepton transverse momentum cuts are set to their defaults so
        // that parseLeptonAcceptance can over-write them if necessary.
        signalDefinitions->setPrimaryLeptonCut( defaultPrimaryLeptonCut );
        signalDefinitions->setSecondaryLeptonCut( defaultSecondaryLeptonCut );

        std::string argumentRemainder;
        if( parseBeamEnergy( validArguments,
                             signalDefinitions,
                             argumentRemainder ) )
        {
          leptonAcceptanceForCascadePair*
          leptonAcceptanceType( parseLeptonAcceptance( argumentRemainder,
                                                       signalDefinitions ) );
          if( NULL != leptonAcceptanceType )
          {
            returnPointer
            = new atlasFourJetMetPlusGivenLeptonCuts( signalDefinitions,
                                                      leptonAcceptanceType );
          }
        }
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
      << "atlasFourJetMetPlusGivenLeptonCuts constructing with primary lepton"
      << " cut = "
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
      subchannelZeroOrMoreJets
      = leptonAcceptanceCalculator->withAtLeastNJets( &signalDefinitions,
                                                      0,
                                                      firstCascade,
                                             firstDoesNotStartWithAntiparticle,
                                                      secondCascade,
                                          secondDoesNotStartWithAntiparticle );
        // chance of the pair of cascades producing the correct number of
        // leptons between them.

      subchannelOneOrMoreJets
      = ( subchannelZeroOrMoreJets
         - leptonAcceptanceCalculator->withExactlyNJets( &signalDefinitions,
                                                         0,
                                                         firstCascade,
                                             firstDoesNotStartWithAntiparticle,
                                                         secondCascade,
                                        secondDoesNotStartWithAntiparticle ) );
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
        fourJetAcceptance = fourJetKinematics->getAcceptance( initialPair,
                                                              firstCascade,
                                                              secondCascade );
        subchannelValue = ( subchannelCrossSectionTimesBrToEwinos
                            * subchannelZeroOrMoreJets
                            * fourJetAcceptance );

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
          << subchannelZeroOrMoreJets
          << ", subchannelOneOrMoreJetsZeroLeptons = "
          << subchannelOneOrMoreJets
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
    // (for N leptons to *pass*) & 10.0 GeV (for all others to *fail*).

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
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "atlasThreeJetMetPlusGivenLeptonCuts::getCalculator( "
      << argumentString << ", ... ) called.";
      std::cout << std::endl;/**/

      signalCalculator* returnPointer( NULL );

      if( 0 == argumentString.compare( 0,
                                       11,
                                       "Atlas3jMET_" ) )
      {
        std::string validArguments( argumentString.substr( 10 ) );
        // debugging:
        /**/std::cout << std::endl << "debugging:"
        << std::endl
        << "validArguments = \"" << validArguments << "\"";
        std::cout << std::endl;/**/

        // the lepton transverse momentum cuts are set to their defaults so
        // that parseLeptonAcceptance can over-write them if necessary.
        signalDefinitions->setPrimaryLeptonCut( defaultPrimaryLeptonCut );
        signalDefinitions->setSecondaryLeptonCut( defaultSecondaryLeptonCut );

        std::string argumentRemainder;
        if( parseBeamEnergy( validArguments,
                             signalDefinitions,
                             argumentRemainder ) )
        {
          leptonAcceptanceForCascadePair*
          leptonAcceptanceType( parseLeptonAcceptance( argumentRemainder,
                                                       signalDefinitions ) );
          if( NULL != leptonAcceptanceType )
          {
            returnPointer
            = new atlasFourJetMetPlusGivenLeptonCuts( signalDefinitions,
                                                      leptonAcceptanceType );
          }
        }
      }
      return returnPointer;
    }

    atlasThreeJetMetPlusGivenLeptonCuts::atlasThreeJetMetPlusGivenLeptonCuts(
                                  signalDefinitionSet* const signalDefinitions,
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
      << signalDefinitions->getSecondaryLeptonCut() << ", numberOfLeptons = "
      << numberOfLeptons;
      std::cout << std::endl;**/

      std::string jetGridName( "Atlas3jMET" );
      threeJetKinematics
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
      = leptonAcceptanceCalculator->withAtLeastNJets( &signalDefinitions,
                                                      0,
                                                      firstCascade,
                                             firstDoesNotStartWithAntiparticle,
                                                      secondCascade,
                                          secondDoesNotStartWithAntiparticle );
        // chance of the pair of cascades producing the correct number of
        // leptons between them.

      subchannelOneOrMoreJets
      = ( subchannelZeroOrMoreJets
         - leptonAcceptanceCalculator->withExactlyNJets( &signalDefinitions,
                                                         0,
                                                         firstCascade,
                                             firstDoesNotStartWithAntiparticle,
                                                         secondCascade,
                                        secondDoesNotStartWithAntiparticle ) );
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
        fourJetAcceptance = fourJetKinematics->getAcceptance( initialPair,
                                                              firstCascade,
                                                              secondCascade );
        subchannelValue = ( subchannelCrossSectionTimesBrToEwinos
                            * subchannelZeroOrMoreJets
                            * fourJetAcceptance );

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
          << subchannelZeroOrMoreJets
          << ", subchannelOneOrMoreJetsZeroLeptons = "
          << subchannelOneOrMoreJets
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
