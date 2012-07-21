/*
 * LHC-FASER_charged_electroweak_cascade_stuff.cpp
 *
 *  Created on: May 9, 2011
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

#include "LHC-FASER_charged_electroweak_cascade_stuff.hpp"

namespace LHC_FASER
{
  charginoToSemuOrEmuSnuCascade::charginoToSemuOrEmuSnuCascade(
                           leptonAcceptanceParameterSet* const kinematicsTable,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                                          particlePointer const coloredDecayer,
                                      particlePointer const electroweakDecayer,
                                     particlePointer const intermediateDecayer,
                                       inputHandler const* const inputShortcut,
                                                bool const sneutrinoVersion ) :
    electroweakCascade( kinematicsTable,
                        effectiveSquarkMass,
                        coloredDecayer,
                        electroweakDecayer,
                        intermediateDecayer,
                        false,
                        inputShortcut ),
    muonsNotElectrons( false ),
    sneutrinoVersion( sneutrinoVersion ),
    sameDistribution( NULL ),
    oppositeDistribution( NULL ),
    jetLeftHandedness( CppSLHA::CppSLHA_global::really_wrong_value ),
    nearLeptonLeftHandedness( CppSLHA::CppSLHA_global::really_wrong_value ),
    nearLeptonSameHandednessTimesBr(
                                 CppSLHA::CppSLHA_global::really_wrong_value ),
    samePass( CppSLHA::CppSLHA_global::really_wrong_value ),
    sameFail( CppSLHA::CppSLHA_global::really_wrong_value ),
    oppositePass( CppSLHA::CppSLHA_global::really_wrong_value ),
    oppositeFail( CppSLHA::CppSLHA_global::really_wrong_value )
  {
    if( ( CppSLHA::PDG_code::smuon_L == intermediateDecayer->get_PDG_code() )
        ||
        ( CppSLHA::PDG_code::smuon_R == intermediateDecayer->get_PDG_code() )
        ||
        ( CppSLHA::PDG_code::muon_sneutrino_L
          == intermediateDecayer->get_PDG_code() )
        ||
        ( CppSLHA::PDG_code::muon_sneutrino_L
          == intermediateDecayer->get_PDG_code() ) )
    {
      muonsNotElectrons = true;
    }
    // we set up the BR of the (only) channel:
    secondBr = inputShortcut->getExclusiveBrCalculator( intermediateDecayer,
                                             inputShortcut->getNeutralinoOne(),
                                                        true,
                                               inputShortcut->getEmptyList() );
    if( sneutrinoVersion )
    {
      firstBr = inputShortcut->getExclusiveBrCalculator( electroweakDecayer,
                                                         intermediateDecayer,
                                                         true,
                                               inputShortcut->getEmptyList() );
      sameDistribution
      = new sameChiralityNearMuon( inputShortcut->getReadier(),
                                   inputShortcut->getCppSlha(),
                                   coloredDecayer,
                                   effectiveSquarkMass,
                                   electroweakDecayer,
                                   intermediateDecayer );
      oppositeDistribution
      = new oppositeChiralityNearMuon( inputShortcut->getReadier(),
                                       inputShortcut->getCppSlha(),
                                       coloredDecayer,
                                       effectiveSquarkMass,
                                       electroweakDecayer,
                                       intermediateDecayer );
    }
    else
    {
      firstBr = inputShortcut->getExclusiveBrCalculator( electroweakDecayer,
                                                         intermediateDecayer,
                                                         false,
                                               inputShortcut->getEmptyList() );
      sameDistribution
      = new sameChiralityFarMuon( inputShortcut->getReadier(),
                                  inputShortcut->getCppSlha(),
                                  coloredDecayer,
                                  effectiveSquarkMass,
                                  electroweakDecayer,
                                  intermediateDecayer,
                                  inputShortcut->getNeutralinoOne() );
      oppositeDistribution
      = new oppositeChiralityFarMuon( inputShortcut->getReadier(),
                                      inputShortcut->getCppSlha(),
                                      coloredDecayer,
                                      effectiveSquarkMass,
                                      electroweakDecayer,
                                      intermediateDecayer,
                                      inputShortcut->getNeutralinoOne() );
    }
    activeDistributions.push_back( sameDistribution );
    activeDistributions.push_back( oppositeDistribution );
  }

  charginoToSemuOrEmuSnuCascade::~charginoToSemuOrEmuSnuCascade()
  {
    // does nothing.
  }


  bool
  charginoToSemuOrEmuSnuCascade::validSignal( int const numberOfJets,
                                           int const numberOfNegativeElectrons,
                                           int const numberOfPositiveElectrons,
                                              int const numberOfNegativeMuons,
                                              int const numberOfPositiveMuons )
  /* this returns true if numberOfJets == 0, & either
   * numberOfPositiveElectrons == 1 OR numberOfPositiveMuons == 1, the rest 0,
   * or if they are all 0.
   */
  {
    if( ( 0 == numberOfJets )
        &&
        ( 0 == numberOfNegativeElectrons )
        &&
        ( 0 == numberOfNegativeMuons )
        &&
        ( ( ( 0 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 1 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 0 == numberOfPositiveElectrons )
            &&
            ( 1 == numberOfPositiveMuons ) ) ) )
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  void
  charginoToSemuOrEmuSnuCascade::calculateAcceptance(
                                     acceptanceCutSet const* const currentCuts,
                                    acceptanceValues* const currentAcceptance )
  // this returns the appropriate acceptances multiplied by branching ratios
  // from the electroweakino through to the LSP and visible final particles.
  {
    cascadeBr = ( firstBr->getBr() * secondBr->getBr() );
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "charginoToSemuOrEmuSnuCascade::calculateAcceptance(...) called."
    << " cascadeBr = " << cascadeBr
    << ", this->intermediateDecayer->get_name() = "
    << *(this->intermediateDecayer->get_name());
    std::cout << std::endl;**/

    if( lhcFaserGlobal::negligibleBr < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      jetLeftHandedness = inputShortcut->quarkOrLeptonLeftHandedness(
                                                coloredDecayer->get_PDG_code(),
                                          electroweakDecayer->get_PDG_code() );
      if( inputShortcut->isIn( coloredDecayer->get_PDG_code(),
                               inputShortcut->getSdownTypes() ) )
        /* if we've calculated the wrong handedness, because isIn assumes the
         * correct charge of chargino for a positive-PDG-coded quark, which is
         * the wrong case for down-type quarks since we are assuming that the
         * chargino is positive...
         */
      {
        jetLeftHandedness = ( 1.0 - jetLeftHandedness );
      }
      nearLeptonLeftHandedness = inputShortcut->quarkOrLeptonLeftHandedness(
                                           intermediateDecayer->get_PDG_code(),
                                          electroweakDecayer->get_PDG_code() );
      if( sneutrinoVersion )
        // positive charginos decaying to sneutrinos produce antileptons, which
        // have the opposite helicity to their charge-conjugates.
      {
        nearLeptonLeftHandedness = ( 1.0 - nearLeptonLeftHandedness );
      }

      nearLeptonSameHandednessTimesBr
      = ( cascadeBr * ( ( jetLeftHandedness * nearLeptonLeftHandedness )
                        + ( ( 1.0 - jetLeftHandedness )
                            * ( 1.0 - nearLeptonLeftHandedness ) ) ) );

      samePass = integrateAcceptance( sameDistribution,
                                      currentCuts->getPrimaryLeptonCut() );
      oppositePass = integrateAcceptance( oppositeDistribution,
                                          currentCuts->getPrimaryLeptonCut() );
      if( currentCuts->leptonCutsAreEqual() )
      {
        sameFail = ( 1.0 - samePass );
        oppositeFail = ( 1.0 - oppositePass );
      }
      else
      {
        sameFail = ( 1.0 - integrateAcceptance( sameDistribution,
                                      currentCuts->getSecondaryLeptonCut() ) );
        oppositeFail = ( 1.0 - integrateAcceptance( oppositeDistribution,
                                      currentCuts->getSecondaryLeptonCut() ) );
      }

      currentAcceptance->setZeroJetsZeroLeptons(
                                              ( nearLeptonSameHandednessTimesBr
                                                * sameFail
                                                + ( cascadeBr
                                            - nearLeptonSameHandednessTimesBr )
                                                  * oppositeFail ) );
      if( muonsNotElectrons )
      {
        currentAcceptance->setZeroJetsOnePositiveMuon(
                                              ( nearLeptonSameHandednessTimesBr
                                                * samePass
                                                + ( cascadeBr
                                            - nearLeptonSameHandednessTimesBr )
                                                  * oppositePass ) );
      }
      else
      {
        currentAcceptance->setZeroJetsOnePositiveElectron(
                                              ( nearLeptonSameHandednessTimesBr
                                                * samePass
                                                + ( cascadeBr
                                            - nearLeptonSameHandednessTimesBr )
                                                  * oppositePass ) );
      }
    }  // end of if cascade branching ratio was large enough to care about.
    // debugging:
    /**std::cout
    << std::endl
    << "cascadeBr = " << cascadeBr;
    std::cout << std::endl;**/
  }



  charginoToStauOrTauSnuCascade::charginoToStauOrTauSnuCascade(
                           leptonAcceptanceParameterSet* const kinematicsTable,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                                          particlePointer const coloredDecayer,
                                      particlePointer const electroweakDecayer,
                                     particlePointer const intermediateDecayer,
                                       inputHandler const* const inputShortcut,
                                                bool const sneutrinoVersion ) :
    electroweakCascade( kinematicsTable,
                        effectiveSquarkMass,
                        coloredDecayer,
                        electroweakDecayer,
                        intermediateDecayer,
                        false,
                        inputShortcut ),
    sneutrinoVersion( sneutrinoVersion ),
    sameTauDistribution( NULL ),
    oppositeTauDistribution( NULL ),
    sameHardMuonDistribution( NULL ),
    sameSoftMuonDistribution( NULL ),
    sameHardPionDistribution( NULL ),
    sameSoftPionDistribution( NULL ),
    oppositeHardMuonDistribution( NULL ),
    oppositeSoftMuonDistribution( NULL ),
    oppositeHardPionDistribution( NULL ),
    oppositeSoftPionDistribution( NULL ),
    jetLeftHandedness( CppSLHA::CppSLHA_global::really_wrong_value ),
    nearLeptonLeftHandedness( CppSLHA::CppSLHA_global::really_wrong_value ),
    antitauRightHandedness( CppSLHA::CppSLHA_global::really_wrong_value ),
    muonPass( CppSLHA::CppSLHA_global::really_wrong_value ),
    muonFail( CppSLHA::CppSLHA_global::really_wrong_value ),
    pionPass( CppSLHA::CppSLHA_global::really_wrong_value ),
    pionFail( CppSLHA::CppSLHA_global::really_wrong_value )
  {
    // we set up the BR of the (only) channel:
    secondBr = inputShortcut->getExclusiveBrCalculator( intermediateDecayer,
                                             inputShortcut->getNeutralinoOne(),
                                                        true,
                                               inputShortcut->getEmptyList() );
    if( sneutrinoVersion )
    {
      firstBr = inputShortcut->getExclusiveBrCalculator( electroweakDecayer,
                                                         intermediateDecayer,
                                                         true,
                                               inputShortcut->getEmptyList() );
      sameTauDistribution
      = new sameChiralityNearMuon( inputShortcut->getReadier(),
                                   inputShortcut->getCppSlha(),
                                   coloredDecayer,
                                   effectiveSquarkMass,
                                   electroweakDecayer,
                                   intermediateDecayer );
      oppositeTauDistribution
      = new oppositeChiralityNearMuon( inputShortcut->getReadier(),
                                       inputShortcut->getCppSlha(),
                                       coloredDecayer,
                                       effectiveSquarkMass,
                                       electroweakDecayer,
                                       intermediateDecayer );
    }
    else
    {
      firstBr = inputShortcut->getExclusiveBrCalculator( electroweakDecayer,
                                                         intermediateDecayer,
                                                         false,
                                               inputShortcut->getEmptyList() );
      sameTauDistribution
      = new sameChiralityFarMuon( inputShortcut->getReadier(),
                                  inputShortcut->getCppSlha(),
                                  coloredDecayer,
                                  effectiveSquarkMass,
                                  electroweakDecayer,
                                  intermediateDecayer,
                                  inputShortcut->getNeutralinoOne() );
      oppositeTauDistribution
      = new oppositeChiralityFarMuon( inputShortcut->getReadier(),
                                      inputShortcut->getCppSlha(),
                                      coloredDecayer,
                                      effectiveSquarkMass,
                                      electroweakDecayer,
                                      intermediateDecayer,
                                      inputShortcut->getNeutralinoOne() );
    }
    activeDistributions.push_back( sameTauDistribution );
    activeDistributions.push_back( oppositeTauDistribution );

    sameHardMuonDistribution
    = new visibleTauDecayProduct( inputShortcut->getReadier(),
                                  sameTauDistribution,
                                  inputShortcut->getHardMuonFromTau() );
    sameSoftMuonDistribution
    = new visibleTauDecayProduct( inputShortcut->getReadier(),
                                  sameTauDistribution,
                                  inputShortcut->getSoftMuonFromTau() );
    sameHardPionDistribution
    = new visibleTauDecayProduct( inputShortcut->getReadier(),
                                  sameTauDistribution,
                                  inputShortcut->getHardPionFromTau() );
    sameSoftPionDistribution
    = new visibleTauDecayProduct( inputShortcut->getReadier(),
                                  sameTauDistribution,
                                  inputShortcut->getSoftPionFromTau() );
    oppositeHardMuonDistribution
    = new visibleTauDecayProduct( inputShortcut->getReadier(),
                                  oppositeTauDistribution,
                                  inputShortcut->getHardMuonFromTau() );
    oppositeSoftMuonDistribution
    = new visibleTauDecayProduct( inputShortcut->getReadier(),
                                  oppositeTauDistribution,
                                  inputShortcut->getSoftMuonFromTau() );
    oppositeHardPionDistribution
    = new visibleTauDecayProduct( inputShortcut->getReadier(),
                                  oppositeTauDistribution,
                                  inputShortcut->getHardPionFromTau() );
    oppositeSoftPionDistribution
    = new visibleTauDecayProduct( inputShortcut->getReadier(),
                                  oppositeTauDistribution,
                                  inputShortcut->getSoftPionFromTau() );
    activeDistributions.push_back( sameHardMuonDistribution );
    activeDistributions.push_back( sameSoftMuonDistribution );
    activeDistributions.push_back( sameHardPionDistribution );
    activeDistributions.push_back( sameSoftPionDistribution );
    activeDistributions.push_back( oppositeHardMuonDistribution );
    activeDistributions.push_back( oppositeSoftMuonDistribution );
    activeDistributions.push_back( oppositeHardPionDistribution );
    activeDistributions.push_back( oppositeSoftPionDistribution );
  }

  charginoToStauOrTauSnuCascade::~charginoToStauOrTauSnuCascade()
  {
    // does nothing.
  }


  bool
  charginoToStauOrTauSnuCascade::validSignal( int const numberOfJets,
                                           int const numberOfNegativeElectrons,
                                           int const numberOfPositiveElectrons,
                                              int const numberOfNegativeMuons,
                                              int const numberOfPositiveMuons )
  /* this returns true if all the arguments are 0, or all are 0 but for
   * numberOfJets OR numberOfPositiveElectrons OR numberOfPositiveMuons, which
   * is 1.
   */
  {
    if( ( 0 == numberOfNegativeElectrons )
        &&
        ( 0 == numberOfNegativeMuons )
        &&
        ( ( ( 0 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 1 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 0 == numberOfJets )
            &&
            ( 1 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 0 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 1 == numberOfPositiveMuons ) ) ) )
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  void
  charginoToStauOrTauSnuCascade::calculateAcceptance(
                                     acceptanceCutSet const* const currentCuts,
                                    acceptanceValues* const currentAcceptance )
  // this returns the appropriate acceptances multiplied by branching ratios
  // from the electroweakino through to the LSP and visible final particles.
  {
    cascadeBr = ( firstBr->getBr() * secondBr->getBr() );
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "charginoToStauOrTauSnuCascade::calculateAcceptance(...) called."
    << " cascadeBr = " << cascadeBr;
    std::cout << std::endl;**/
    if( lhcFaserGlobal::negligibleBr < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      jetLeftHandedness = inputShortcut->quarkOrLeptonLeftHandedness(
                                                coloredDecayer->get_PDG_code(),
                                          electroweakDecayer->get_PDG_code() );
      if( inputShortcut->isIn( coloredDecayer->get_PDG_code(),
                               inputShortcut->getSdownTypes() ) )
        /* if we've calculated the wrong handedness, because isIn assumes the
         * correct charge of chargino for a positive-PDG-coded quark, which is
         * the wrong case for down-type quarks since we are assuming that the
         * chargino is positive...
         */
      {
        jetLeftHandedness = ( 1.0 - jetLeftHandedness );
      }
      nearLeptonLeftHandedness = inputShortcut->quarkOrLeptonLeftHandedness(
                                           intermediateDecayer->get_PDG_code(),
                                          electroweakDecayer->get_PDG_code() );
      if( sneutrinoVersion )
        // positive charginos decaying to sneutrinos produce antileptons, which
        // have the opposite helicity to their charge-conjugates.
      {
        antitauRightHandedness = nearLeptonLeftHandedness;
        nearLeptonLeftHandedness = ( 1.0 - antitauRightHandedness );
      }
      else
      {
        antitauRightHandedness = inputShortcut->quarkOrLeptonLeftHandedness(
                                           intermediateDecayer->get_PDG_code(),
                           inputShortcut->getNeutralinoOne()->get_PDG_code() );
      }

      if( sneutrinoVersion )
      {
        muonPass
        = ( cascadeBr
            * ( jetLeftHandedness
                * ( nearLeptonLeftHandedness
                    * integrateAcceptance( sameSoftMuonDistribution,
                                           currentCuts->getPrimaryLeptonCut() )
                    + ( 1.0 - nearLeptonLeftHandedness )
                      * integrateAcceptance( oppositeHardMuonDistribution,
                                         currentCuts->getPrimaryLeptonCut() ) )
                + ( 1.0 - jetLeftHandedness )
                  * ( nearLeptonLeftHandedness
                      * integrateAcceptance( oppositeSoftMuonDistribution,
                                           currentCuts->getPrimaryLeptonCut() )
                      + ( 1.0 - nearLeptonLeftHandedness )
                        * integrateAcceptance( sameHardMuonDistribution,
                                    currentCuts->getPrimaryLeptonCut() ) ) ) );
        if( currentCuts->leptonCutsAreEqual() )
        {
          muonFail = ( cascadeBr * ( 1.0 - muonPass ) );
        }
        else
        {
          muonFail
          = ( cascadeBr
              * ( 1.0
                  - ( jetLeftHandedness
                      * ( nearLeptonLeftHandedness
                          * integrateAcceptance( sameSoftMuonDistribution,
                                         currentCuts->getSecondaryLeptonCut() )
                          + ( 1.0 - nearLeptonLeftHandedness )
                           * integrateAcceptance( oppositeHardMuonDistribution,
                                       currentCuts->getSecondaryLeptonCut() ) )
                      + ( 1.0 - jetLeftHandedness )
                        * ( nearLeptonLeftHandedness
                           * integrateAcceptance( oppositeSoftMuonDistribution,
                                         currentCuts->getSecondaryLeptonCut() )
                            + ( 1.0 - nearLeptonLeftHandedness )
                              * integrateAcceptance( sameHardMuonDistribution,
                                currentCuts->getSecondaryLeptonCut() ) ) ) ) );
        }
        pionPass
        = ( cascadeBr
            * ( jetLeftHandedness
                * ( nearLeptonLeftHandedness
                    * integrateAcceptance( sameHardPionDistribution,
                                           currentCuts->getJetCut() )
                    + ( 1.0 - nearLeptonLeftHandedness )
                      * integrateAcceptance( oppositeSoftPionDistribution,
                                             currentCuts->getJetCut() ) )
                + ( 1.0 - jetLeftHandedness )
                  * ( nearLeptonLeftHandedness
                      * integrateAcceptance( oppositeHardPionDistribution,
                                             currentCuts->getJetCut() )
                      + ( 1.0 - nearLeptonLeftHandedness )
                        * integrateAcceptance( sameSoftPionDistribution,
                                           currentCuts->getJetCut() ) ) ) );
      }
      else
      {
        muonPass
        = ( cascadeBr
            * ( jetLeftHandedness
                * ( antitauRightHandedness
                    * integrateAcceptance( sameHardMuonDistribution,
                                           currentCuts->getPrimaryLeptonCut() )
                    + ( 1.0 - antitauRightHandedness )
                      * integrateAcceptance( sameSoftMuonDistribution,
                                         currentCuts->getPrimaryLeptonCut() ) )
                + ( 1.0 - jetLeftHandedness )
                  * ( antitauRightHandedness
                      * integrateAcceptance( oppositeHardMuonDistribution,
                                           currentCuts->getPrimaryLeptonCut() )
                      + ( 1.0 - antitauRightHandedness )
                        * integrateAcceptance( oppositeSoftMuonDistribution,
                                    currentCuts->getPrimaryLeptonCut() ) ) ) );
        if( currentCuts->leptonCutsAreEqual() )
        {
          muonFail = ( cascadeBr * ( 1.0 - muonPass ) );
        }
        else
        {
          muonFail
                  = ( cascadeBr
                      * ( 1.0
                          - ( jetLeftHandedness
                              * ( antitauRightHandedness
                                  * integrateAcceptance(
                                                      sameHardMuonDistribution,
                                         currentCuts->getSecondaryLeptonCut() )
                                  + ( 1.0 - antitauRightHandedness )
                                    * integrateAcceptance(
                                                      sameSoftMuonDistribution,
                                       currentCuts->getSecondaryLeptonCut() ) )
                              + ( 1.0 - jetLeftHandedness )
                                * ( antitauRightHandedness
                                    * integrateAcceptance(
                                                  oppositeHardMuonDistribution,
                                         currentCuts->getSecondaryLeptonCut() )
                                    + ( 1.0 - antitauRightHandedness )
                                      * integrateAcceptance(
                                                  oppositeSoftMuonDistribution,
                                currentCuts->getSecondaryLeptonCut() ) ) ) ) );
        }
        pionPass
        = ( cascadeBr
            * ( jetLeftHandedness
                * ( antitauRightHandedness
                    * integrateAcceptance( sameSoftPionDistribution,
                                           currentCuts->getJetCut() )
                    + ( 1.0 - antitauRightHandedness )
                      * integrateAcceptance( sameHardPionDistribution,
                                             currentCuts->getJetCut() ) )
                + ( 1.0 - jetLeftHandedness )
                  * ( antitauRightHandedness
                      * integrateAcceptance( oppositeSoftPionDistribution,
                                             currentCuts->getJetCut() )
                      + ( 1.0 - antitauRightHandedness )
                        * integrateAcceptance( oppositeHardPionDistribution,
                                           currentCuts->getJetCut() ) ) ) );
      }


      currentAcceptance->setOneJetZeroLeptons( pionPass
                       * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR );
      currentAcceptance->setZeroJetsOnePositiveElectron( muonPass
                    * CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR );
      currentAcceptance->setZeroJetsOnePositiveMuon( muonPass
                        * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR );

      currentAcceptance->setZeroJetsZeroLeptons( ( cascadeBr - pionPass )
                          * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
                                                 + muonFail
                     * ( CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR
                      + CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR ) );
    }  // end of if cascade branching ratio was large enough to care about.
    // debugging:
    /**std::cout
    << std::endl
    << "cascadeBr = " << cascadeBr;
    std::cout << std::endl;**/
  }



  charginoToWCascade::charginoToWCascade(
                           leptonAcceptanceParameterSet* const kinematicsTable,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                                          particlePointer const coloredDecayer,
                                      particlePointer const electroweakDecayer,
                                    inputHandler const* const inputShortcut ) :
    electroweakCascade( kinematicsTable,
                        effectiveSquarkMass,
                        coloredDecayer,
                        electroweakDecayer,
                        inputShortcut->getWPlus(),
                        false,
                        inputShortcut ),
    downTypeQuark( false ),
    directJetDistribution( NULL ),
    leftHandedJetRightHandedAntimuonDistribution( NULL ),
    rightHandedJetRightHandedAntimuonDistribution( NULL ),
    leftHandedJetTauMuonDistribution( NULL ),
    rightHandedJetTauMuonDistribution( NULL ),
    leftHandedJetTauPionDistribution( NULL ),
    rightHandedJetTauPionDistribution( NULL ),
    currentMuonDistribution( NULL ),
    currentPionDistribution( NULL ),
    currentCuts( NULL ),
    currentAcceptance( NULL ),
    jetLeftHandedness( CppSLHA::CppSLHA_global::really_wrong_value ),
    directMuonPass( CppSLHA::CppSLHA_global::really_wrong_value ),
    directMuonFail( CppSLHA::CppSLHA_global::really_wrong_value ),
    directJetPass( CppSLHA::CppSLHA_global::really_wrong_value ),
    directJetFail( CppSLHA::CppSLHA_global::really_wrong_value ),
    antitauAntimuonPass( CppSLHA::CppSLHA_global::really_wrong_value ),
    antitauAntimuonFail( CppSLHA::CppSLHA_global::really_wrong_value ),
    antitauPionPass( CppSLHA::CppSLHA_global::really_wrong_value ),
    antitauPionFail( CppSLHA::CppSLHA_global::really_wrong_value ),
    currentPass( CppSLHA::CppSLHA_global::really_wrong_value ),
    configurationBr( CppSLHA::CppSLHA_global::really_wrong_value )
  {
    if( inputShortcut->isIn( coloredDecayer->get_PDG_code(),
                             inputShortcut->getSdownTypes() ) )
    {
      downTypeQuark = true;
    }

    // here firstBr & secondBr don't work as well as just getting the branching
    // ratios from the CppSLHA::particle_property_set pointers.

    directJetDistribution = new wMinusDirectJet( inputShortcut->getReadier(),
                                                 inputShortcut->getCppSlha(),
                                                 coloredDecayer,
                                                 effectiveSquarkMass,
                                                 electroweakDecayer,
                                                 intermediateDecayer,
                                           inputShortcut->getNeutralinoOne() );
    activeDistributions.push_back( directJetDistribution );

    /* in my infinite wisdom & boundless foresight, I wrote the distribution
     * for a negatively-charged chargino decaying to a W^- which decays to a
     * left-handed muon, given the handedness of the jet. however, here we want
     * the charge-conjugate process.
     */
    leftHandedJetRightHandedAntimuonDistribution
    = new wMinusHandedMuon( inputShortcut->getReadier(),
                            inputShortcut->getCppSlha(),
                            coloredDecayer,
                            effectiveSquarkMass,
                            electroweakDecayer,
                            intermediateDecayer,
                            inputShortcut->getNeutralinoOne(),
                            false /* see above */ );
    rightHandedJetRightHandedAntimuonDistribution
    = new wMinusHandedMuon( inputShortcut->getReadier(),
                            inputShortcut->getCppSlha(),
                            coloredDecayer,
                            effectiveSquarkMass,
                            electroweakDecayer,
                            intermediateDecayer,
                            inputShortcut->getNeutralinoOne(),
                            true /* see above */ );
    activeDistributions.push_back(
                                leftHandedJetRightHandedAntimuonDistribution );
    activeDistributions.push_back(
                               rightHandedJetRightHandedAntimuonDistribution );

    leftHandedJetTauMuonDistribution
    = new visibleTauDecayProduct( inputShortcut->getReadier(),
                                  leftHandedJetRightHandedAntimuonDistribution,
                                  inputShortcut->getHardMuonFromTau() );
    rightHandedJetTauMuonDistribution
    = new visibleTauDecayProduct( inputShortcut->getReadier(),
                                 rightHandedJetRightHandedAntimuonDistribution,
                                  inputShortcut->getHardMuonFromTau() );
    leftHandedJetTauPionDistribution
    = new visibleTauDecayProduct( inputShortcut->getReadier(),
                                  leftHandedJetRightHandedAntimuonDistribution,
                                  inputShortcut->getSoftPionFromTau() );
    rightHandedJetTauPionDistribution
    = new visibleTauDecayProduct( inputShortcut->getReadier(),
                                 rightHandedJetRightHandedAntimuonDistribution,
                                  inputShortcut->getSoftPionFromTau() );
    activeDistributions.push_back( leftHandedJetTauMuonDistribution );
    activeDistributions.push_back( rightHandedJetTauMuonDistribution );
    activeDistributions.push_back( leftHandedJetTauPionDistribution );
    activeDistributions.push_back( rightHandedJetTauPionDistribution );
  }

  charginoToWCascade::~charginoToWCascade()
  {
    // does nothing.
  }


  bool
  charginoToWCascade::validSignal( int const numberOfJets,
                                   int const numberOfNegativeElectrons,
                                   int const numberOfPositiveElectrons,
                                   int const numberOfNegativeMuons,
                                   int const numberOfPositiveMuons )
  /* this returns true if all the arguments are 0, or all are 0 but for
   * numberOfJets OR numberOfPositiveElectrons OR numberOfPositiveMuons, which
   * is 1, or all are jet but numberOfJets which is 2.
   */
  {
    if( ( 0 == numberOfNegativeElectrons )
        &&
        ( 0 == numberOfNegativeMuons )
        &&
        ( ( ( 0 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 1 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 2 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 0 == numberOfJets )
            &&
            ( 1 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 0 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 1 == numberOfPositiveMuons ) ) ) )
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  void
  charginoToWCascade::calculateAcceptance(
                                     acceptanceCutSet const* const currentCuts,
                                    acceptanceValues* const currentAcceptance )
  // this returns the appropriate acceptances multiplied by branching ratios
  // from the electroweakino through to the LSP and visible final particles.
  {
    currentAcceptance->setTwoJets( 0.0 );
    currentAcceptance->setOneJetZeroLeptons( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveElectron( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveMuon( 0.0 );
    currentAcceptance->setZeroJetsZeroLeptons( 0.0 );
    cascadeBr
    = ( electroweakDecayer->inspect_direct_decay_handler(
             )->get_branching_ratio_for_exact_match( CppSLHA::PDG_code::W_plus,
                                         CppSLHA::PDG_code::neutralino_one ) );
    // the branching ratios of the W into the various SM fermions are covered
    // by the following code.

    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "charginoToWCascade::calculateAcceptance(...) called."
    << " cascadeBr = " << cascadeBr;
    std::cout << std::endl;**/

    if( lhcFaserGlobal::negligibleBr < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      /* at some point, it'd be nice to break this down into transverse
       * contributions & longitudinal contributions, since actually the
       * Goldstone boson contribution should both be large (only needs the
       * Higgsino component of just 1 of the electroweakinos, rather than both)
       * & has a different (& easier to calculate) distribution.
       *
       * however, I think that for LHC-relevant energies, the Higgsino
       * components of all the neutralinos is still going to be larger than
       * m_tau / m_W, enough that I'm going to ignore the longitudinal
       * contributions.
       */

      jetLeftHandedness = inputShortcut->quarkOrLeptonLeftHandedness(
                                                coloredDecayer->get_PDG_code(),
                                          electroweakDecayer->get_PDG_code() );
      if( inputShortcut->isIn( coloredDecayer->get_PDG_code(),
                               inputShortcut->getSdownTypes() ) )
        /* if we've calculated the wrong handedness, because isIn assumes the
         * correct charge of chargino for a positive-PDG-coded quark, which is
         * the wrong case for down-type quarks since we are assuming that the
         * chargino is positive...
         */
      {
        jetLeftHandedness = ( 1.0 - jetLeftHandedness );
      }

      this->currentCuts = currentCuts;
      this->currentAcceptance = currentAcceptance;

      directMuonPass
      = ( cascadeBr
          * ( jetLeftHandedness
           * integrateAcceptance( leftHandedJetRightHandedAntimuonDistribution,
                                  currentCuts->getPrimaryLeptonCut() )
               + ( 1.0 - jetLeftHandedness )
          * integrateAcceptance( rightHandedJetRightHandedAntimuonDistribution,
                                 currentCuts->getPrimaryLeptonCut() ) ) );

      if( currentCuts->leptonCutsAreEqual() )
      {
        directMuonFail = ( cascadeBr * ( 1.0 - directMuonPass ) );
      }
      else
      {
        directMuonFail
        = ( cascadeBr
            * ( 1.0
                - ( jetLeftHandedness
                    * integrateAcceptance(
                                  leftHandedJetRightHandedAntimuonDistribution,
                                           currentCuts->getPrimaryLeptonCut() )
                    + ( 1.0 - jetLeftHandedness )
                      * integrateAcceptance(
                                 rightHandedJetRightHandedAntimuonDistribution,
                                  currentCuts->getSecondaryLeptonCut() ) ) ) );
      }

      directJetPass = integrateAcceptance( directJetDistribution,
                                           currentCuts->getJetCut() );
      directJetFail = ( 1.0 - directJetPass );

      configurationBr
      = ( cascadeBr * CppSLHA::PDG_data::W_plus_to_hadrons_BR
                    * directJetPass );
      currentAcceptance->setTwoJets( configurationBr * directJetPass );
      currentAcceptance->setOneJetZeroLeptons( 2.0 * configurationBr
                                                   * directJetFail );
      currentAcceptance->setZeroJetsOnePositiveElectron( directMuonPass
                     * CppSLHA::PDG_data::W_plus_to_neutrino_antielectron_BR );
      currentAcceptance->setZeroJetsOnePositiveMuon( directMuonPass
                         * CppSLHA::PDG_data::W_plus_to_neutrino_antimuon_BR );

      currentAcceptance->setZeroJetsZeroLeptons( cascadeBr
                                      * CppSLHA::PDG_data::W_plus_to_hadrons_BR
                                                 * directJetFail
                                                 * directJetFail
                      + ( CppSLHA::PDG_data::W_plus_to_neutrino_antielectron_BR
                          + CppSLHA::PDG_data::W_plus_to_neutrino_antimuon_BR )
                        * directMuonFail );


      // now W^+ decays to tau antileptons, followed by the decays of the tau
      // antileptons:
      cascadeBr *= CppSLHA::PDG_data::W_plus_to_neutrino_tau_antilepton_BR;
      // left-handed jet:
      configurationBr
      = ( cascadeBr * jetLeftHandedness );
      currentMuonDistribution = leftHandedJetTauMuonDistribution;
      currentPionDistribution = leftHandedJetTauPionDistribution;
      calculateForCurrentConfiguration();

      // right-handed jet:
      configurationBr
      = ( cascadeBr * ( 1.0 - jetLeftHandedness ) );
      currentMuonDistribution = rightHandedJetTauMuonDistribution;
      currentPionDistribution = rightHandedJetTauPionDistribution;
      calculateForCurrentConfiguration();
    }
    // debugging:
    /**std::cout
    << std::endl
    << "cascadeBr = " << cascadeBr;
    std::cout << std::endl;**/
  }

  void
  charginoToWCascade::calculateForCurrentConfiguration()
  {
    antitauAntimuonPass = integrateAcceptance( currentMuonDistribution,
                                          currentCuts->getPrimaryLeptonCut() );
    antitauPionPass = integrateAcceptance( currentPionDistribution,
                                       currentCuts->getJetCut() );
    if( currentCuts->leptonCutsAreEqual() )
    {
      antitauAntimuonFail = ( 1.0 - antitauAntimuonPass );
    }
    else
    {
      antitauAntimuonFail
      = ( 1.0 - integrateAcceptance( currentMuonDistribution,
                                     currentCuts->getSecondaryLeptonCut() ) );
    }
    antitauPionFail = ( 1.0 - antitauPionPass );

    currentAcceptance->addToOneJetZeroLeptons( configurationBr
                          * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
                                               * antitauPionPass );
    currentPass = ( configurationBr * antitauAntimuonPass );
    currentAcceptance->addToZeroJetsOnePositiveElectron( currentPass
                    * CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR );
    currentAcceptance->addToZeroJetsOnePositiveMuon( currentPass
                       * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR  );

    currentAcceptance->addToZeroJetsZeroLeptons( configurationBr
                        * ( CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
                                                 * antitauPionFail
                     + ( CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR
                         + CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR )
                       * antitauAntimuonFail ) );
  }



  charginoToHiggsCascade::charginoToHiggsCascade(
                           leptonAcceptanceParameterSet* const kinematicsTable,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                                          particlePointer const coloredDecayer,
                                      particlePointer const electroweakDecayer,
                                     particlePointer const intermediateDecayer,
                                    inputHandler const* const inputShortcut ) :
    electroweakCascade( kinematicsTable,
                        effectiveSquarkMass,
                        coloredDecayer,
                        electroweakDecayer,
                        intermediateDecayer,
                        false,
                        inputShortcut ),
    directAntimuonDistribution( NULL ),
    antitauAntimuonDistribution( NULL ),
    antitauPionDistribution( NULL ),
    currentCuts( NULL ),
    currentAcceptance( NULL ),
    jetLeftHandedness( CppSLHA::CppSLHA_global::really_wrong_value ),
    directMuonPass( CppSLHA::CppSLHA_global::really_wrong_value ),
    directMuonFail( CppSLHA::CppSLHA_global::really_wrong_value ),
    directJetPass( CppSLHA::CppSLHA_global::really_wrong_value ),
    directJetFail( CppSLHA::CppSLHA_global::really_wrong_value ),
    antitauAntimuonPass( CppSLHA::CppSLHA_global::really_wrong_value ),
    antitauAntimuonFail( CppSLHA::CppSLHA_global::really_wrong_value ),
    antitauPionPass( CppSLHA::CppSLHA_global::really_wrong_value ),
    antitauPionFail( CppSLHA::CppSLHA_global::really_wrong_value ),
    currentPass( CppSLHA::CppSLHA_global::really_wrong_value ),
    configurationBr( CppSLHA::CppSLHA_global::really_wrong_value ),
    intermediateDecayerDecays(
                         intermediateDecayer->inspect_direct_decay_handler() ),
    currentBrToHadrons( CppSLHA::CppSLHA_global::really_wrong_value ),
    currentBrToAntielectron( CppSLHA::CppSLHA_global::really_wrong_value ),
    currentBrToAntimuon( CppSLHA::CppSLHA_global::really_wrong_value ),
    currentBrToAntitau( CppSLHA::CppSLHA_global::really_wrong_value )
  {
    // here firstBr & secondBr don't work as well as just getting the branching
    // ratios from the CppSLHA::particle_property_set pointers.
    directAntimuonDistribution
    = new negativelyChargedHiggsMuon( inputShortcut->getReadier(),
                                      inputShortcut->getCppSlha(),
                                      coloredDecayer,
                                      effectiveSquarkMass,
                                      electroweakDecayer,
                                      intermediateDecayer,
                                      inputShortcut->getNeutralinoOne() );
    activeDistributions.push_back( directAntimuonDistribution );

    // the spin-0 nature of the boson combined with the left-handed nature of
    // the tau neutrino forces the tau antilepton to be left-handed.
    antitauAntimuonDistribution
    = new visibleTauDecayProduct( inputShortcut->getReadier(),
                                  directAntimuonDistribution,
                                  inputShortcut->getSoftMuonFromTau() );
    antitauPionDistribution
    = new visibleTauDecayProduct( inputShortcut->getReadier(),
                                  directAntimuonDistribution,
                                  inputShortcut->getHardPionFromTau() );
    activeDistributions.push_back( antitauAntimuonDistribution );
    activeDistributions.push_back( antitauPionDistribution );
  }

  charginoToHiggsCascade::~charginoToHiggsCascade()
  {
    // does nothing.
  }


  bool
  charginoToHiggsCascade::validSignal( int const numberOfJets,
                                       int const numberOfNegativeElectrons,
                                       int const numberOfPositiveElectrons,
                                       int const numberOfNegativeMuons,
                                       int const numberOfPositiveMuons )
  /* this returns true if all the arguments are 0, or all are 0 but for
   * numberOfJets OR numberOfPositiveElectrons OR numberOfPositiveMuons, which
   * is 1, or all are jet but numberOfJets which is 2.
   */
  {
    if( ( 0 == numberOfNegativeElectrons )
        &&
        ( 0 == numberOfNegativeMuons )
        &&
        ( ( ( 0 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 1 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 2 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 0 == numberOfJets )
            &&
            ( 1 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 0 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 1 == numberOfPositiveMuons ) ) ) )
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  void
  charginoToHiggsCascade::calculateAcceptance(
                                     acceptanceCutSet const* const currentCuts,
                                    acceptanceValues* const currentAcceptance )
  // this returns the appropriate acceptances multiplied by branching ratios
  // from the electroweakino through to the LSP and visible final particles.
  {
    currentAcceptance->setTwoJets( 0.0 );
    currentAcceptance->setOneJetZeroLeptons( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveElectron( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveMuon( 0.0 );
    currentAcceptance->setZeroJetsZeroLeptons( 0.0 );
    cascadeBr = ( electroweakDecayer->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                           intermediateDecayer->get_PDG_code(),
                                         CppSLHA::PDG_code::neutralino_one ) );
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "charginoToHiggsCascade::calculateAcceptance(...) called."
    << " cascadeBr = " << cascadeBr;
    std::cout << std::endl;**/

    // the branching ratios of the scalar into the various SM fermions are
    // covered by the following code.
    if( lhcFaserGlobal::negligibleBr < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      currentBrToHadrons
      = ( intermediateDecayerDecays->get_branching_ratio_for_exact_match(
                                                         CppSLHA::PDG_code::up,
                                                  -(CppSLHA::PDG_code::down ) )
          + intermediateDecayerDecays->get_branching_ratio_for_exact_match(
                                                         CppSLHA::PDG_code::up,
                                               -(CppSLHA::PDG_code::strange ) )
          + intermediateDecayerDecays->get_branching_ratio_for_exact_match(
                                                         CppSLHA::PDG_code::up,
                                                -(CppSLHA::PDG_code::bottom ) )
          + intermediateDecayerDecays->get_branching_ratio_for_exact_match(
                                                    CppSLHA::PDG_code::charm,
                                               -(CppSLHA::PDG_code::down ) )
          + intermediateDecayerDecays->get_branching_ratio_for_exact_match(
                                                      CppSLHA::PDG_code::charm,
                                               -(CppSLHA::PDG_code::strange ) )
          + intermediateDecayerDecays->get_branching_ratio_for_exact_match(
                                                      CppSLHA::PDG_code::charm,
                                             -(CppSLHA::PDG_code::bottom ) ) );
      currentBrToAntielectron
      = intermediateDecayerDecays->get_branching_ratio_for_exact_match(
                                          CppSLHA::PDG_code::electron_neutrino,
                                             -(CppSLHA::PDG_code::electron ) );
      currentBrToAntimuon
      = intermediateDecayerDecays->get_branching_ratio_for_exact_match(
                                              CppSLHA::PDG_code::muon_neutrino,
                                                 -(CppSLHA::PDG_code::muon ) );
      currentBrToAntitau
      = intermediateDecayerDecays->get_branching_ratio_for_exact_match(
                                               CppSLHA::PDG_code::tau_neutrino,
                                           -(CppSLHA::PDG_code::tau_lepton ) );

      directMuonPass
      = ( cascadeBr * integrateAcceptance( directAntimuonDistribution,
                                        currentCuts->getPrimaryLeptonCut() ) );
      directJetPass = integrateAcceptance( directAntimuonDistribution,
                                           currentCuts->getJetCut() );
      if( currentCuts->leptonCutsAreEqual() )
      {
        directMuonFail = ( 1.0 - directMuonPass );
      }
      else
      {
        directMuonFail
        = ( cascadeBr
            * ( 1.0 - integrateAcceptance( directAntimuonDistribution,
                                      currentCuts->getPrimaryLeptonCut() ) ) );
      }
      directJetFail = ( 1.0 - directJetPass );

      configurationBr = ( cascadeBr * currentBrToHadrons * directJetPass );
      currentAcceptance->setTwoJets( configurationBr * directJetPass );
      currentAcceptance->setOneJetZeroLeptons( 2.0 * configurationBr
                                                   * directJetFail );
      currentAcceptance->setZeroJetsOnePositiveElectron( directMuonPass
                                                   * currentBrToAntielectron );
      currentAcceptance->setZeroJetsOnePositiveMuon( directMuonPass
                                                     * currentBrToAntimuon );

      currentAcceptance->setZeroJetsZeroLeptons( cascadeBr
                                                 * currentBrToHadrons
                                                 * directJetFail
                                                 * directJetFail
                                                 + ( currentBrToAntielectron
                                                     + currentBrToAntimuon )
                                                   * directMuonFail );

      // now decays to tau antileptons, followed by the decays of the tau
      // antileptons:
      cascadeBr *= currentBrToAntitau;
      antitauAntimuonPass = integrateAcceptance( antitauAntimuonDistribution,
                                          currentCuts->getPrimaryLeptonCut() );
      antitauPionPass = integrateAcceptance( antitauPionDistribution,
                                             currentCuts->getJetCut() );
      if( currentCuts->leptonCutsAreEqual() )
      {
        antitauAntimuonFail = ( 1.0 - antitauAntimuonPass );
      }
      else
      {
        antitauAntimuonFail
        = ( 1.0 - integrateAcceptance( antitauAntimuonDistribution,
                                      currentCuts->getSecondaryLeptonCut() ) );
      }
      antitauPionFail = ( 1.0 - antitauPionPass );

      currentAcceptance->addToOneJetZeroLeptons( configurationBr
                          * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
                                                 * antitauPionPass );
      currentPass = ( configurationBr * antitauAntimuonPass );
      currentAcceptance->addToZeroJetsOnePositiveElectron( currentPass
                    * CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR );
      currentAcceptance->addToZeroJetsOnePositiveMuon( currentPass
                        * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR );

      currentAcceptance->addToZeroJetsZeroLeptons( configurationBr
                        * ( CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
                                                   * antitauPionFail
                     + ( CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR
                         + CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR )
                       * antitauAntimuonFail ) );
    }
    // debugging:
    /**std::cout
    << std::endl
    << "cascadeBr = " << cascadeBr;
    std::cout << std::endl;**/
  }



  charginoVirtualCascade::charginoVirtualCascade(
                           leptonAcceptanceParameterSet* const kinematicsTable,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                                          particlePointer const coloredDecayer,
                                      particlePointer const electroweakDecayer,
                                    inputHandler const* const inputShortcut ) :
    electroweakCascade( kinematicsTable,
                        effectiveSquarkMass,
                        coloredDecayer,
                        electroweakDecayer,
                        NULL,
                        false,
                        inputShortcut ),
    directDownUpDistribution( NULL ),
    directStrangeCharmDistribution( NULL ),
    directElectronDistribution( NULL ),
    directMuonDistribution( NULL ),
    directTauDistribution( NULL ),
    hardMuonDistribution( NULL ),
    softMuonDistribution( NULL ),
    hardPionDistribution( NULL ),
    softPionDistribution( NULL ),
    currentJetDistribution( NULL ),
    currentCuts( NULL ),
    currentAcceptance( NULL ),
    negativeTauLeftHandedness( CppSLHA::CppSLHA_global::really_wrong_value ),
    jetLeftHandedness( CppSLHA::CppSLHA_global::really_wrong_value ),
    directPass( CppSLHA::CppSLHA_global::really_wrong_value ),
    directFail( CppSLHA::CppSLHA_global::really_wrong_value ),
    hardMuonPass( CppSLHA::CppSLHA_global::really_wrong_value ),
    hardMuonFail( CppSLHA::CppSLHA_global::really_wrong_value ),
    softMuonPass( CppSLHA::CppSLHA_global::really_wrong_value ),
    softMuonFail( CppSLHA::CppSLHA_global::really_wrong_value ),
    averageMuonPass( CppSLHA::CppSLHA_global::really_wrong_value ),
    averageMuonFail( CppSLHA::CppSLHA_global::really_wrong_value ),
    hardPionPass( CppSLHA::CppSLHA_global::really_wrong_value ),
    hardPionFail( CppSLHA::CppSLHA_global::really_wrong_value ),
    softPionPass( CppSLHA::CppSLHA_global::really_wrong_value ),
    softPionFail( CppSLHA::CppSLHA_global::really_wrong_value ),
    averagePionPass( CppSLHA::CppSLHA_global::really_wrong_value ),
    averagePionFail( CppSLHA::CppSLHA_global::really_wrong_value ),
    currentPass( CppSLHA::CppSLHA_global::really_wrong_value ),
    electroweakDecayerDecays(
                           electroweakDecayer->inspect_direct_decay_handler() )
  {
    negativeTauLeftHandedness = 0.5;
    // this will probably have to be changed if we ever do the off-shell decays
    // properly.

    // here firstBr & secondBr don't work as well as just getting the branching
    // ratios from the CppSLHA::particle_property_set pointers.

    directDownUpDistribution
    = new charginoThreeBodyDecay( inputShortcut->getReadier(),
                                  inputShortcut->getCppSlha(),
                                  coloredDecayer,
                                  effectiveSquarkMass,
                                  electroweakDecayer,
                                  inputShortcut->getNeutralinoOne(),
                                  inputShortcut->getSupL(),
                                  inputShortcut->getSupR(),
                                  inputShortcut->getSdownL(),
                                  inputShortcut->getSdownR() );
    directStrangeCharmDistribution
    = new charginoThreeBodyDecay( inputShortcut->getReadier(),
                                  inputShortcut->getCppSlha(),
                                  coloredDecayer,
                                  effectiveSquarkMass,
                                  electroweakDecayer,
                                  inputShortcut->getNeutralinoOne(),
                                  inputShortcut->getScharmL(),
                                  inputShortcut->getScharmR(),
                                  inputShortcut->getSstrangeL(),
                                  inputShortcut->getSstrangeR() );
    directElectronDistribution
    = new charginoThreeBodyDecay( inputShortcut->getReadier(),
                                  inputShortcut->getCppSlha(),
                                  coloredDecayer,
                                  effectiveSquarkMass,
                                  electroweakDecayer,
                                  inputShortcut->getNeutralinoOne(),
                                  inputShortcut->getElectronSneutrinoL(),
                                  inputShortcut->getElectronSneutrinoR(),
                                  inputShortcut->getSelectronL(),
                                  inputShortcut->getSelectronR() );
    directMuonDistribution
    = new charginoThreeBodyDecay( inputShortcut->getReadier(),
                                  inputShortcut->getCppSlha(),
                                  coloredDecayer,
                                  effectiveSquarkMass,
                                  electroweakDecayer,
                                  inputShortcut->getNeutralinoOne(),
                                  inputShortcut->getMuonSneutrinoL(),
                                  inputShortcut->getMuonSneutrinoR(),
                                  inputShortcut->getSmuonL(),
                                  inputShortcut->getSmuonR() );
    directTauDistribution
    = new charginoThreeBodyDecay( inputShortcut->getReadier(),
                                  inputShortcut->getCppSlha(),
                                  coloredDecayer,
                                  effectiveSquarkMass,
                                  electroweakDecayer,
                                  inputShortcut->getNeutralinoOne(),
                                  inputShortcut->getTauSneutrinoL(),
                                  inputShortcut->getTauSneutrinoR(),
                                  inputShortcut->getStauOne(),
                                  inputShortcut->getStauTwo() );
    activeDistributions.push_back( directDownUpDistribution );
    activeDistributions.push_back( directStrangeCharmDistribution );
    activeDistributions.push_back( directElectronDistribution );
    activeDistributions.push_back( directMuonDistribution );
    activeDistributions.push_back( directTauDistribution );

    hardMuonDistribution
    = new visibleTauDecayProduct( inputShortcut->getReadier(),
                                  directTauDistribution,
                                  inputShortcut->getHardMuonFromTau() );
    softMuonDistribution
    = new visibleTauDecayProduct( inputShortcut->getReadier(),
                                  directTauDistribution,
                                  inputShortcut->getSoftMuonFromTau() );
    hardPionDistribution
    = new visibleTauDecayProduct( inputShortcut->getReadier(),
                                  directTauDistribution,
                                  inputShortcut->getHardPionFromTau() );
    softPionDistribution
    = new visibleTauDecayProduct( inputShortcut->getReadier(),
                                  directTauDistribution,
                                  inputShortcut->getSoftPionFromTau() );
    activeDistributions.push_back( hardMuonDistribution );
    activeDistributions.push_back( softMuonDistribution );
    activeDistributions.push_back( hardPionDistribution );
    activeDistributions.push_back( softPionDistribution );
  }

  charginoVirtualCascade::~charginoVirtualCascade()
  {
    // does nothing.
  }


  bool
  charginoVirtualCascade::validSignal( int const numberOfJets,
                                       int const numberOfNegativeElectrons,
                                       int const numberOfPositiveElectrons,
                                       int const numberOfNegativeMuons,
                                       int const numberOfPositiveMuons )
  /* this returns true if all the arguments are 0, or all are 0 but for
   * numberOfJets OR numberOfPositiveElectrons OR numberOfPositiveMuons, which
   * is 1, or all are jet but numberOfJets which is 2.
   */
  {
    if( ( 0 == numberOfNegativeElectrons )
        &&
        ( 0 == numberOfNegativeMuons )
        &&
        ( ( ( 0 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 1 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 2 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 0 == numberOfJets )
            &&
            ( 1 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 0 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 1 == numberOfPositiveMuons ) ) ) )
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  void
  charginoVirtualCascade::calculateAcceptance(
                                     acceptanceCutSet const* const currentCuts,
                                    acceptanceValues* const currentAcceptance )
  // this returns the appropriate acceptances multiplied by branching ratios
  // from the electroweakino through to the LSP and visible final particles.
  {
    currentAcceptance->setTwoJets( 0.0 );
    currentAcceptance->setOneJetZeroLeptons( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveElectron( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveMuon( 0.0 );
    currentAcceptance->setZeroJetsZeroLeptons( 0.0 );

    // we don't consider flavor non-diagonal decays: it's a lot of effort for
    // little gain.
    cascadeBr = electroweakDecayerDecays->get_branching_ratio_for_exact_match(
                                             CppSLHA::PDG_code::neutralino_one,
                                                         CppSLHA::PDG_code::up,
                                                  -(CppSLHA::PDG_code::down) );
    if( lhcFaserGlobal::negligibleBr < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      directPass = integrateAcceptance( directDownUpDistribution,
                                        currentCuts->getJetCut() );
      directFail = ( 1.0 - directPass );
      currentAcceptance->addToTwoJets( cascadeBr * directPass * directPass );
      currentAcceptance->addToOneJetZeroLeptons( 2.0 * cascadeBr * directPass
                                                                * directFail );
      currentAcceptance->addToZeroJetsZeroLeptons( cascadeBr * directFail );
    }
    cascadeBr = electroweakDecayerDecays->get_branching_ratio_for_exact_match(
                                             CppSLHA::PDG_code::neutralino_one,
                                                      CppSLHA::PDG_code::charm,
                                               -(CppSLHA::PDG_code::strange) );
    if( lhcFaserGlobal::negligibleBr < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      directPass = integrateAcceptance( directStrangeCharmDistribution,
                                        currentCuts->getJetCut() );
      directFail = ( 1.0 - directPass );
      currentAcceptance->addToTwoJets( cascadeBr * directPass * directPass );
      currentAcceptance->addToOneJetZeroLeptons( 2.0 * cascadeBr
                                                 * directPass * directFail );
      currentAcceptance->addToZeroJetsZeroLeptons( cascadeBr * directFail
                                                             * directFail );
    }
    cascadeBr = electroweakDecayerDecays->get_branching_ratio_for_exact_match(
                                             CppSLHA::PDG_code::neutralino_one,
                                          CppSLHA::PDG_code::electron_neutrino,
                                              -(CppSLHA::PDG_code::electron) );
    if( lhcFaserGlobal::negligibleBr < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      directPass = integrateAcceptance( directElectronDistribution,
                                        currentCuts->getPrimaryLeptonCut() );
      if( currentCuts->leptonCutsAreEqual() )
      {
        directFail = ( 1.0 - directPass );
      }
      else
      {
        directFail = ( 1.0 - integrateAcceptance( directElectronDistribution,
                                      currentCuts->getSecondaryLeptonCut() ) );
      }
      currentAcceptance->addToZeroJetsOnePositiveElectron( cascadeBr
                                                           * directPass );
      currentAcceptance->addToZeroJetsZeroLeptons( cascadeBr * directFail );
    }
    cascadeBr = electroweakDecayerDecays->get_branching_ratio_for_exact_match(
                                             CppSLHA::PDG_code::neutralino_one,
                                              CppSLHA::PDG_code::muon_neutrino,
                                                  -(CppSLHA::PDG_code::muon) );
    if( lhcFaserGlobal::negligibleBr < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      directPass = integrateAcceptance( directMuonDistribution,
                                        currentCuts->getPrimaryLeptonCut() );
      if( currentCuts->leptonCutsAreEqual() )
      {
        directFail = ( 1.0 - directPass );
      }
      else
      {
        directFail = ( 1.0 - integrateAcceptance( directMuonDistribution,
                                      currentCuts->getSecondaryLeptonCut() ) );
      }
      currentAcceptance->addToZeroJetsOnePositiveMuon( cascadeBr
                                                       * directPass );
      currentAcceptance->addToZeroJetsZeroLeptons( cascadeBr
                                                   * directFail );
    }
    cascadeBr = electroweakDecayerDecays->get_branching_ratio_for_exact_match(
                                             CppSLHA::PDG_code::neutralino_one,
                                               CppSLHA::PDG_code::tau_neutrino,
                                            -(CppSLHA::PDG_code::tau_lepton) );
    if( lhcFaserGlobal::negligibleBr < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      hardMuonPass = integrateAcceptance( hardMuonDistribution,
                                          currentCuts->getPrimaryLeptonCut() );
      softMuonPass = integrateAcceptance( hardMuonDistribution,
                                          currentCuts->getPrimaryLeptonCut() );
      averageMuonPass = ( 0.5 * ( hardMuonPass + softMuonPass ) );
      hardPionPass = integrateAcceptance( hardPionDistribution,
                                          currentCuts->getJetCut() );
      softPionPass = integrateAcceptance( softPionDistribution,
                                          currentCuts->getJetCut() );
      averagePionPass = ( 0.5 * ( hardPionPass + softPionPass ) );
      if( currentCuts->leptonCutsAreEqual() )
      {
        hardMuonFail = ( 1.0 - hardMuonPass );
        softMuonFail = ( 1.0 - softMuonPass );
      }
      else
      {
        hardMuonFail = ( 1.0 - integrateAcceptance( hardMuonDistribution,
                                      currentCuts->getSecondaryLeptonCut() ) );
        softMuonFail = ( 1.0 - integrateAcceptance( hardMuonDistribution,
                                      currentCuts->getSecondaryLeptonCut() ) );
      }
      averageMuonFail = ( 0.5 * ( hardMuonFail + softMuonFail ) );
      hardPionFail = ( 1.0 - hardPionPass );
      softPionFail = ( 1.0 - softPionPass );
      averagePionFail = ( 0.5 * ( hardPionFail + softPionFail ) );

      currentAcceptance->addToOneJetZeroLeptons( cascadeBr
                          * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
                                                 * averagePionPass );
      currentPass = ( cascadeBr * averageMuonPass );
      currentAcceptance->addToZeroJetsOnePositiveElectron( currentPass
                    * CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR );
      currentAcceptance->addToZeroJetsOnePositiveMuon( currentPass
                       * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR  );

      currentAcceptance->addToZeroJetsZeroLeptons( cascadeBr
                        * ( CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
                            * averagePionFail
                     + ( CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR
                         + CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR )
                         * averageMuonFail ) );
    }
    // debugging:
    /**std::cout
    << std::endl
    << "cascadeBr = " << cascadeBr;
    std::cout << std::endl;**/
  }



  scoloredToWPlusScoloredCascade::scoloredToWPlusScoloredCascade(
                           leptonAcceptanceParameterSet* const kinematicsTable,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                                          particlePointer const coloredDecayer,
                                      particlePointer const electroweakDecayer,
                                         particlePointer const lighterScolored,
                                    inputHandler const* const inputShortcut ) :
    electroweakCascade( kinematicsTable,
                        effectiveSquarkMass,
                        coloredDecayer,
                        lighterScolored,
                        inputShortcut->getWPlus(),
                        false,
                        inputShortcut ),
    directMuonDistribution( NULL ),
    antitauAntimuonDistribution( NULL ),
    antitauPionDistribution( NULL ),
    directMuonPass( CppSLHA::CppSLHA_global::really_wrong_value ),
    directMuonFail( CppSLHA::CppSLHA_global::really_wrong_value ),
    directJetPass( CppSLHA::CppSLHA_global::really_wrong_value ),
    directJetFail( CppSLHA::CppSLHA_global::really_wrong_value ),
    antitauAntimuonPass( CppSLHA::CppSLHA_global::really_wrong_value ),
    antitauAntimuonFail( CppSLHA::CppSLHA_global::really_wrong_value ),
    antitauPionPass( CppSLHA::CppSLHA_global::really_wrong_value ),
    antitauPionFail( CppSLHA::CppSLHA_global::really_wrong_value ),
    currentPass( CppSLHA::CppSLHA_global::really_wrong_value ),
    configurationBr( CppSLHA::CppSLHA_global::really_wrong_value ),
    wPdgCode( CppSLHA::PDG_code::W_plus )
  {
    // here firstBr & secondBr don't work as well as just getting the branching
    // ratios from the CppSLHA::particle_property_set pointers.
    if( inputShortcut->isIn( coloredDecayer->get_PDG_code(),
                             inputShortcut->getSdownTypes() ) )
    {
      wPdgCode = -(CppSLHA::PDG_code::W_plus);
    }

    directMuonDistribution
    = new vectorFromSquarkToMuon( inputShortcut->getReadier(),
                                  inputShortcut->getCppSlha(),
                                  coloredDecayer,
                                  effectiveSquarkMass,
                                  lighterScolored,
                                  intermediateDecayer );
    activeDistributions.push_back( directMuonDistribution );

    antitauAntimuonDistribution
    = new visibleTauDecayProduct( inputShortcut->getReadier(),
                                  directMuonDistribution,
                                  inputShortcut->getHardMuonFromTau() );
    antitauPionDistribution
    = new visibleTauDecayProduct( inputShortcut->getReadier(),
                                  directMuonDistribution,
                                  inputShortcut->getSoftPionFromTau() );
    activeDistributions.push_back( antitauAntimuonDistribution );
    activeDistributions.push_back( antitauPionDistribution );
  }

  scoloredToWPlusScoloredCascade::~scoloredToWPlusScoloredCascade()
  {
    // does nothing.
  }


  bool
  scoloredToWPlusScoloredCascade::validSignal( int const numberOfJets,
                                           int const numberOfNegativeElectrons,
                                           int const numberOfPositiveElectrons,
                                               int const numberOfNegativeMuons,
                                              int const numberOfPositiveMuons )
  /* this returns true if all the arguments are 0, or all are 0 but for
   * numberOfJets OR numberOfPositiveElectrons OR numberOfPositiveMuons, which
   * is 1, or all are jet but numberOfJets which is 2.
   */
  {
    if( ( 0 == numberOfNegativeElectrons )
        &&
        ( 0 == numberOfNegativeMuons )
        &&
        ( ( ( 0 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 1 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 2 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 0 == numberOfJets )
            &&
            ( 1 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 0 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 1 == numberOfPositiveMuons ) ) ) )
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "scoloredToWPlusScoloredCascade::validSignal() returning true";
      std::cout << std::endl;**/

      return true;
    }
    else
    {
      return false;
    }
  }

  void
  scoloredToWPlusScoloredCascade::calculateAcceptance(
                                     acceptanceCutSet const* const currentCuts,
                                    acceptanceValues* const currentAcceptance )
  // this returns the appropriate acceptances multiplied by branching ratios
  // from the electroweakino through to the LSP and visible final particles.
  {
    currentAcceptance->setTwoJets( 0.0 );
    currentAcceptance->setOneJetZeroLeptons( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveElectron( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveMuon( 0.0 );
    currentAcceptance->setZeroJetsZeroLeptons( 0.0 );
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "scoloredToWPlusScoloredCascade::calculateAcceptance( [ "
    << currentCuts->getJetCut() << "J, "
    << currentCuts->getPrimaryLeptonCut() << "L1, "
    << currentCuts->getSecondaryLeptonCut() << "L2 ] ) called."
    << " effectiveSquarkMass->getEffectiveSquarkMass() = "
    << effectiveSquarkMass->getEffectiveSquarkMass()
    << ", lighter squark mass = "
    << electroweakDecayer->get_absolute_mass();
    std::cout << std::endl;**/

    /* the branching ratios of the W into the various SM fermions are covered
     * by the following code.  the BR into the W itself should already be
     * covered.
     */
    /* at some point, it'd be nice to break this down into transverse
     * contributions & longitudinal contributions, since actually the
     * Goldstone boson contribution could both be large (for sbottom -> stop,
     * for example) & has a different (& easier to calculate) distribution.
     *
     * however, ( ( m_top * m_tau ) / m_W^2 ) is still pretty small, so it
     * should be fine to neglect the longitudinal contributions.
     */

    directMuonPass = integrateAcceptance( directMuonDistribution,
                                          currentCuts->getPrimaryLeptonCut() );
    directJetPass = integrateAcceptance( directMuonDistribution,
                                         currentCuts->getJetCut() );
    if( currentCuts->leptonCutsAreEqual() )
    {
      directMuonFail = ( 1.0 - directMuonPass );
    }
    else
    {
      directMuonFail
      = ( 1.0 - integrateAcceptance( directMuonDistribution,
                                     currentCuts->getSecondaryLeptonCut() ) );
    }
    directJetFail = ( 1.0 - directJetPass );

    configurationBr
    = ( CppSLHA::PDG_data::W_plus_to_hadrons_BR * directJetPass );
    currentAcceptance->setTwoJets( configurationBr * directJetPass );
    currentAcceptance->setOneJetZeroLeptons( 2.0 * configurationBr
                                                 * directJetFail );
    currentAcceptance->setZeroJetsOnePositiveElectron( directMuonPass
                     * CppSLHA::PDG_data::W_plus_to_neutrino_antielectron_BR );
    currentAcceptance->setZeroJetsOnePositiveMuon( directMuonPass
                         * CppSLHA::PDG_data::W_plus_to_neutrino_antimuon_BR );
    currentAcceptance->setZeroJetsZeroLeptons(
                                        CppSLHA::PDG_data::W_plus_to_hadrons_BR
                                        * directJetFail
                                        * directJetFail
                            + CppSLHA::PDG_data::W_plus_to_neutrino_antimuon_BR
                              * directMuonFail );

    // now W^+ decays to tau antileptons, followed by the decays of the tau
    // antileptons:
    antitauAntimuonPass = integrateAcceptance( antitauAntimuonDistribution,
                                          currentCuts->getPrimaryLeptonCut() );
    antitauPionPass = integrateAcceptance( antitauPionDistribution,
                                           currentCuts->getJetCut() );
    if( currentCuts->leptonCutsAreEqual() )
    {
      antitauAntimuonFail = ( 1.0 - antitauAntimuonPass );
    }
    else
    {
      antitauAntimuonFail
      = ( 1.0 - integrateAcceptance( antitauAntimuonDistribution,
                                     currentCuts->getSecondaryLeptonCut() ) );
    }
    antitauPionFail = ( 1.0 - antitauPionPass );

    currentAcceptance->addToOneJetZeroLeptons(
                        CppSLHA::PDG_data::W_plus_to_neutrino_tau_antilepton_BR
                        * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
                        * antitauPionPass );
    currentPass = ( CppSLHA::PDG_data::W_plus_to_neutrino_tau_antilepton_BR
                    * antitauAntimuonPass );
    currentAcceptance->addToZeroJetsOnePositiveElectron( currentPass
                    * CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR );
    currentAcceptance->addToZeroJetsOnePositiveMuon( currentPass
                       * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR  );

    currentAcceptance->addToZeroJetsZeroLeptons(
                        CppSLHA::PDG_data::W_plus_to_neutrino_tau_antilepton_BR
                        * ( CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
                            * antitauPionFail
                     + ( CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR
                         + CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR )
                       * antitauAntimuonFail ) );

    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "directMuonPass = " << directMuonPass
    << std::endl
    << "directMuonFail = " << directMuonFail
    << std::endl
    << "directJetPass = " << directJetPass
    << std::endl
    << "directJetFail = " << directJetFail
    << std::endl
    << "antitauAntimuonPass = " << antitauAntimuonPass
    << std::endl
    << "antitauAntimuonFail = " << antitauAntimuonFail
    << std::endl
    << "antitauPionPass = " << antitauPionPass
    << std::endl
    << "antitauPionFail = " << antitauPionFail
    << std::endl
    << "currentAcceptance->getTwoJets() = " << currentAcceptance->getTwoJets()
    << std::endl
    << "currentAcceptance->getOneJetZeroLeptons() = "
    << currentAcceptance->getOneJetZeroLeptons()
    << std::endl
    << "currentAcceptance->getZeroJetsOnePositiveElectron() = "
    << currentAcceptance->getZeroJetsOnePositiveElectron()
    << std::endl
    << "currentAcceptance->getZeroJetsOnePositiveMuon() = "
    << currentAcceptance->getZeroJetsOnePositiveMuon()
    << std::endl
    << "currentAcceptance->getZeroJetsZeroLeptons() = "
    << currentAcceptance->getZeroJetsZeroLeptons();
    std::cout << std::endl;**/
  }

}  // end of LHC_FASER namespace.
