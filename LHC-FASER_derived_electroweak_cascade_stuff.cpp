/*
 * LHC-FASER_derived_electroweak_cascade_stuff.cpp
 *
 *  Created on: May 9, 2011
 *      Author: bol
 */

#include "LHC-FASER_derived_electroweak_cascade_stuff.hpp"

namespace LHC_FASER
{

  lightestNeutralinoCascade::lightestNeutralinoCascade() :
    electroweakCascade( NULL,
                        NULL,
                        false,
                        NULL,
                        false,
                        NULL,
                        false,
                        NULL )
  {
    acceptances.setComparison( &acceptanceCutSet::justReturnTrue );
  }

  lightestNeutralinoCascade::~lightestNeutralinoCascade()
  {
    // does nothing.
  }



  neutralinoToSemuCascade::neutralinoToSemuCascade(
                                leptonAcceptanceParameterSet* const kinematics,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                                  //bool const coloredDecayerIsNotAntiparticle,
                CppSLHA::particle_property_set const* const electroweakDecayer,
                              //bool const electroweakDecayerIsNotAntiparticle,
               CppSLHA::particle_property_set const* const intermediateDecayer,
                                        input_handler const* const shortcut ) :
    electroweakCascade( kinematics,
                        coloredDecayer,
                        //coloredDecayerIsNotAntiparticle,
                        electroweakDecayer,
                        //electroweakDecayerIsNotAntiparticle,
                        intermediateDecayer,
                        true,
                        shortcut )
  {
    if( ( CppSLHA::PDG_code::smuon_L == intermediateDecayer->get_PDG_code() )
        ||
        ( CppSLHA::PDG_code::smuon_R == intermediateDecayer->get_PDG_code() ) )
    {
      muonsNotElectrons = true;
    }
    else
    {
      muonsNotElectrons = false;
    }
    // we set up the BR of the (only) channel:
    firstBr
    = shortcut->get_exclusive_BR_calculator( electroweakDecayer,
                                             intermediateDecayer,
                                             true,
                                             shortcut->get_empty_list() );
    secondBr
    = shortcut->get_exclusive_BR_calculator( intermediateDecayer,
                                             shortcut->get_neutralino_one(),
                                             true,
                                             shortcut->get_empty_list() );
    nearSameDistribution
    = new same_chirality_near_muon( shortcut->get_readier(),
                                    shortcut->get_CppSLHA(),
                                    coloredDecayer,
                                    kinematics,
                                    electroweakDecayer,
                                    intermediateDecayer );
    nearOppositeDistribution
    = new opposite_chirality_near_muon( shortcut->get_readier(),
                                        shortcut->get_CppSLHA(),
                                        coloredDecayer,
                                        kinematics,
                                        electroweakDecayer,
                                        intermediateDecayer );
    farSameDistribution
    = new same_chirality_far_muon( shortcut->get_readier(),
                                   shortcut->get_CppSLHA(),
                                   coloredDecayer,
                                   kinematics,
                                   electroweakDecayer,
                                   intermediateDecayer,
                                   shortcut->get_neutralino_one() );
    farOppositeDistribution
    = new opposite_chirality_far_muon( shortcut->get_readier(),
                                       shortcut->get_CppSLHA(),
                                       coloredDecayer,
                                       kinematics,
                                       electroweakDecayer,
                                       intermediateDecayer,
                                       shortcut->get_neutralino_one() );
    activeDistributions.push_back( nearSameDistribution );
    activeDistributions.push_back( nearOppositeDistribution );
    activeDistributions.push_back( farSameDistribution );
    activeDistributions.push_back( farOppositeDistribution );
  }

  neutralinoToSemuCascade::~neutralinoToSemuCascade()
  {
    // does nothing.
  }


  void
  neutralinoToSemuCascade::calculateAcceptance( acceptanceCutSet* const cuts,
                                    acceptanceValues* const currentAcceptance )
  // this returns the appropriate acceptances multiplied by branching ratios
  // from the electroweakino through the selectron or smuon to the LSP.
  {
    cascadeBr = ( first_BR->get_BR() * second_BR->get_BR() );
    if( LHC_FASER_global::negligible_BR < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      jetLeftHandedness
      = shortcut->quark_or_lepton_left_handedness(
                                                coloredDecayer->get_PDG_code(),
                                          electroweakDecayer->get_PDG_code() );
      nearNegativeLeptonLeftHandedness
      = shortcut->quark_or_lepton_left_handedness(
                                           intermediateDecayer->get_PDG_code(),
                                          electroweakDecayer->get_PDG_code() );
      nearNegativeLeptonSameHandednessTimesBr
      = ( cascadeBr * ( ( jetLeftHandedness
                          * nearNegativeLeptonLeftHandedness )
                        + ( ( 1.0 - jetLeftHandedness )
                            * ( 1.0 - nearNegativeLeptonLeftHandedness ) ) ) );

      nearSamePass = integrateAcceptance( nearSameDistribution,
                                          cuts->getPrimaryLeptonCut() );
      nearSameFail = ( 1.0 - integrateAcceptance( nearSameDistribution,
                                             cuts->getSecondaryLeptonCut() ) );
      farSamePass = integrateAcceptance( farSameDistribution,
                                         cuts->getPrimaryLeptonCut() );
      farSameFail = ( 1.0 - integrateAcceptance( farSameDistribution,
                                             cuts->getSecondaryLeptonCut() ) );
      nearOppositePass = integrateAcceptance( nearOppositeDistribution,
                                          cuts->getPrimaryLeptonCut() );
      nearOppositeFail = ( 1.0 - integrateAcceptance( nearOppositeDistribution,
                                             cuts->getSecondaryLeptonCut() ) );
      farOppositePass = integrateAcceptance( farOppositeDistribution,
                                         cuts->getPrimaryLeptonCut() );
      farOppositeFail = ( 1.0 - integrateAcceptance( farOppositeDistribution,
                                             cuts->getSecondaryLeptonCut() ) );

      currentAcceptance->setOssfMinusOsdf(
                    ( cascadeBr * ( nearSamePass * farSamePass
                                    + nearOppositePass * farOppositePass ) ) );
      currentAcceptance->setZeroJetsZeroLeptons(
                    ( cascadeBr * ( nearSameFail * farSameFail
                                    + nearOppositeFail * farOppositeFail ) ) );
      justOneNegativeLeptonPass
      = ( nearNegativeLeptonSameHandednessTimesBr
          * ( nearSamePass * farSameFail + farOppositePass * nearOppositeFail )
          + ( cascadeBr - nearNegativeLeptonSameHandednessTimesBr )
            * ( farSamePass * nearSameFail
                + nearOppositePass * farOppositeFail ) );
      justOnePositiveLeptonPass
      = ( nearNegativeLeptonSameHandednessTimesBr
          * ( farSamePass * nearSameFail + nearOppositePass * farOppositeFail )
          + ( cascadeBr - nearNegativeLeptonSameHandednessTimesBr )
            * ( nearSamePass * farSameFail
                + farOppositePass * nearOppositeFail ) );

      if( muonsNotElectrons )
      {
        currentAcceptance->setMuonPlusAntimuon(
                                       currentAcceptance->getOssfMinusOsdf() );
        currentAcceptance->setZeroJetsOneNegativeMuon(
                                                   justOneNegativeLeptonPass );
        currentAcceptance->setZeroJetsOnePositiveMuon(
                                                   justOnePositiveLeptonPass );
      }
      else
      {
        currentAcceptance->setElectronPlusAntielectron(
                                       currentAcceptance->getOssfMinusOsdf() );
        currentAcceptance->setZeroJetsOneNegativeElectron(
                                                   justOneNegativeLeptonPass );
        currentAcceptance->setZeroJetsOnePositiveElectron(
                                                   justOnePositiveLeptonPass );
      }
    }  // end of if cascade branching ratio was large enough to care about.
    // debugging:
    /**std::cout
    << std::endl
    << "cascadeBr = " << cascadeBr;
    std::cout << std::endl;**/
  }



  chargeSummedNeutralinoToSemuCascade::chargeSummedNeutralinoToSemuCascade(
                                         leptonAcceptanceParameterSet* const kinematics,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                                    bool const coloredDecayerIsNotAntiparticle,
                CppSLHA::particle_property_set const* const electroweakDecayer,
                                bool const electroweakDecayerIsNotAntiparticle,
               CppSLHA::particle_property_set const* const intermediateDecayer,
                                        input_handler const* const shortcut ) :
    electroweakCascade( kinematics,
                        coloredDecayer,
                        coloredDecayerIsNotAntiparticle,
                        electroweakDecayer,
                        electroweakDecayerIsNotAntiparticle,
                        intermediateDecayer,
                        true,
                        shortcut ),
    nearDistribution( NULL ),
    farDistribution( NULL )
  {
    if( ( CppSLHA::PDG_code::smuon_L == intermediateDecayer->get_PDG_code() )
        ||
        ( CppSLHA::PDG_code::smuon_R == intermediateDecayer->get_PDG_code() ) )
    {
      muonsNotElectrons = true;
    }
    else
    {
      muonsNotElectrons = false;
    }
    // we set up the BR of the (only) channel:
    firstBr
    = shortcut->get_exclusive_BR_calculator( electroweakDecayer,
                                             intermediateDecayer,
                                             true,
                                             shortcut->get_empty_list() );
    secondBr
    = shortcut->get_exclusive_BR_calculator( intermediateDecayer,
                                             shortcut->get_neutralino_one(),
                                             true,
                                             shortcut->get_empty_list() );
    nearDistribution
    = new same_chirality_near_muon( shortcut->get_readier(),
                                    shortcut->get_CppSLHA(),
                                    coloredDecayer,
                        kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                    electroweakDecayer,
                                    intermediateDecayer );
    farDistribution
    = new same_chirality_far_muon( shortcut->get_readier(),
                                   shortcut->get_CppSLHA(),
                                   coloredDecayer,
                        kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                   electroweakDecayer,
                                   intermediateDecayer,
                                   shortcut->get_neutralino_one() );
    activeDistributions.push_back( nearDistribution );
    activeDistributions.push_back( farDistribution );
  }

  chargeSummedNeutralinoToSemuCascade::~chargeSummedNeutralinoToSemuCascade()
  {
    // does nothing.
  }


  void
  chargeSummedNeutralinoToSemuCascade::calculateAcceptance(
                                                  acceptanceCutSet* const cuts,
                                    acceptanceValues* const currentAcceptance )
  // this returns the appropriate acceptances multiplied by branching ratios
  // from the electroweakino through the selectron or smuon to the LSP.
  {
    cascadeBr = ( first_BR->get_BR() * second_BR->get_BR() );
    // this is just for *1* of the possibilities (e.g. negatively-charged
    // slepton) - the charge-conjugate is accounted for with factors of 2.
    if( LHC_FASER_global::negligible_BR < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      nearPass = integrateAcceptance( nearDistribution,
                                      cuts->getPrimaryLeptonCut() );
      nearFail = ( 1.0 - integrateAcceptance( nearDistribution,
                                             cuts->getSecondaryLeptonCut() ) );
      farPass = integrateAcceptance( farDistribution,
                                     cuts->getPrimaryLeptonCut() );
      farFail = ( 1.0 - integrateAcceptance( farDistribution,
                                             cuts->getSecondaryLeptonCut() ) );
      bothPass = ( cascadeBr * ( 2.0 * nearPass * farPass ) );
      currentAcceptance->setOssfMinusOsdf( bothPass );
      // the factor of 2.0 accounts for the charge-conjugate cascade.
      currentAcceptance->setZeroJetsZeroLeptons(
                                ( cascadeBr * ( 2.0 * nearFail * farFail ) ) );
      // the factor of 2.0 accounts for the charge-conjugate cascade.
      onePassOneFail
      = ( cascadeBr * ( nearPass * farFail + nearFail * farPass ) );

      if( muonsNotElectrons )
      {
        currentAcceptance->setMuonPlusAntimuon( bothPass );
        currentAcceptance->setZeroJetsOneNegativeMuon( onePassOneFail );
        currentAcceptance->setZeroJetsOnePositiveMuon( onePassOneFail );
      }
      else
      {
        currentAcceptance->setElectronPlusAntielectron( bothPass );
        currentAcceptance->setZeroJetsOneNegativeElectron( onePassOneFail );
        currentAcceptance->setZeroJetsOnePositiveElectron( onePassOneFail );
      }
    }  // end of if cascade branching ratio was large enough to care about.
    // debugging:
    /**std::cout
    << std::endl
    << "cascade_BR = " << cascade_BR;
    std::cout << std::endl;**/
  }



  neutralinoToStauCascade::neutralinoToStauCascade(
                                leptonAcceptanceParameterSet* const kinematics,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                CppSLHA::particle_property_set const* const electroweakDecayer,
               CppSLHA::particle_property_set const* const intermediateDecayer,
                                        input_handler const* const shortcut ) :
    electroweakCascade( kinematics,
                        coloredDecayer,
                        electroweakDecayer,
                        intermediateDecayer,
                        false,
                        shortcut )
  {
    // we set up the BR of the (only-ish) channel:
    firstBr
    = shortcut->get_exclusive_BR_calculator( electroweakDecayer,
                                             intermediateDecayer,
                                             true,
                                             shortcut->get_empty_list() );
    secondBr
    = shortcut->get_exclusive_BR_calculator( intermediateDecayer,
                                             shortcut->get_neutralino_one(),
                                             true,
                                             shortcut->get_empty_list() );
    nearSameTauDistribution
    = new same_chirality_near_muon( shortcut->get_readier(),
                                    shortcut->get_CppSLHA(),
                                    coloredDecayer,
                                    kinematics,
                                    electroweakDecayer,
                                    intermediateDecayer );
    nearOppositeTauDistribution
    = new opposite_chirality_near_muon( shortcut->get_readier(),
                                        shortcut->get_CppSLHA(),
                                        coloredDecayer,
                                        kinematics,
                                        electroweakDecayer,
                                        intermediateDecayer );
    farSameTauDistribution
    = new same_chirality_far_muon( shortcut->get_readier(),
                                   shortcut->get_CppSLHA(),
                                   coloredDecayer,
                                   kinematics,
                                   electroweakDecayer,
                                   intermediateDecayer,
                                   shortcut->get_neutralino_one() );
    farOppositeTauDistribution
    = new opposite_chirality_far_muon( shortcut->get_readier(),
                                       shortcut->get_CppSLHA(),
                                       coloredDecayer,
                                       kinematics,
                                       electroweakDecayer,
                                       intermediateDecayer,
                                       shortcut->get_neutralino_one() );
    activeDistributions.push_back( nearSameTauDistribution );
    activeDistributions.push_back( nearOppositeTauDistribution );
    activeDistributions.push_back( farSameTauDistribution );
    activeDistributions.push_back( farOppositeTauDistribution );

    nearSameHardMuonDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     nearSameTauDistribution,
                                     shortcut->get_hard_muon_from_tau() );
    nearSameSoftMuonDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     nearSameTauDistribution,
                                     shortcut->get_soft_muon_from_tau() );
    nearSameHardPionDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     nearSameTauDistribution,
                                     shortcut->get_hard_pion_from_tau() );
    nearSameSoftPionDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     nearSameTauDistribution,
                                     shortcut->get_soft_pion_from_tau() );
    nearOppositeHardMuonDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     nearOppositeTauDistribution,
                                     shortcut->get_hard_muon_from_tau() );
    nearOppositeSoftMuonDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     nearOppositeTauDistribution,
                                     shortcut->get_soft_muon_from_tau() );
    nearOppositeHardPionDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     nearOppositeTauDistribution,
                                     shortcut->get_hard_pion_from_tau() );
    nearOppositeSoftPionDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     nearOppositeTauDistribution,
                                     shortcut->get_soft_pion_from_tau() );
    farSameHardMuonDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     farSameTauDistribution,
                                     shortcut->get_hard_muon_from_tau() );
    farSameSoftMuonDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     farSameTauDistribution,
                                     shortcut->get_soft_muon_from_tau() );
    farSameHardPionDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     farSameTauDistribution,
                                     shortcut->get_hard_pion_from_tau() );
    farSameSoftPionDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     farSameTauDistribution,
                                     shortcut->get_soft_pion_from_tau() );
    farOppositeHardMuonDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     farOppositeTauDistribution,
                                     shortcut->get_hard_muon_from_tau() );
    farOppositeSoftMuonDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     farOppositeTauDistribution,
                                     shortcut->get_soft_muon_from_tau() );
    farOppositeHardPionDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     farOppositeTauDistribution,
                                     shortcut->get_hard_pion_from_tau() );
    farOppositeSoftPionDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     farOppositeTauDistribution,
                                     shortcut->get_soft_pion_from_tau() );
    activeDistributions.push_back( nearSameHardMuonDistribution );
    activeDistributions.push_back( nearSameSoftMuonDistribution );
    activeDistributions.push_back( nearSameHardPionDistribution );
    activeDistributions.push_back( nearSameSoftPionDistribution );
    activeDistributions.push_back( nearOppositeHardMuonDistribution );
    activeDistributions.push_back( nearOppositeSoftMuonDistribution );
    activeDistributions.push_back( nearOppositeHardPionDistribution );
    activeDistributions.push_back( nearOppositeSoftPionDistribution );
    activeDistributions.push_back( farSameHardMuonDistribution );
    activeDistributions.push_back( farSameSoftMuonDistribution );
    activeDistributions.push_back( farSameHardPionDistribution );
    activeDistributions.push_back( farSameSoftPionDistribution );
    activeDistributions.push_back( farOppositeHardMuonDistribution );
    activeDistributions.push_back( farOppositeSoftMuonDistribution );
    activeDistributions.push_back( farOppositeHardPionDistribution );
    activeDistributions.push_back( farOppositeSoftPionDistribution );
  }

  neutralinoToStauCascade::~neutralinoToStauCascade()
  {
    // does nothing.
  }


  bool
  neutralinoToStauCascade::validSignal( int const numberOfJets,
                                        int const numberOfNegativeElectrons,
                                        int const numberOfPositiveElectrons,
                                        int const numberOfNegativeMuons,
                                        int const numberOfPositiveMuons )
  /* this returns true if a configuration where each of the signs of tau lepton
   * decayed either into a detected jet, detected lepton, or undetected
   * particle, & false otherwise.
   */
  {
    if( 2 == numberOfJets )
    {
      if( ( 0 == numberOfNegativeElectrons )
          &&
          ( 0 == numberOfPositiveElectrons )
          &&
          ( 0 == numberOfNegativeMuons )
          &&
          ( 0 == numberOfPositiveMuons ) )
      {
        return true;
      }
      else
      {
        return false;
      }
    }
    else if( ( 0 <= numberOfJets )
             &&
             ( 0 <= numberOfNegativeElectrons )
             &&
             ( 0 <= numberOfPositiveElectrons )
             &&
             ( 0 <= numberOfNegativeMuons )
             &&
             ( 0 <= numberOfPositiveMuons )
             &&
             ( 1 >= ( numberOfNegativeElectrons + numberOfNegativeMuons ) )
             &&
             ( 1 >= ( numberOfPositiveElectrons + numberOfPositiveMuons ) ) )
    {
      if( 2 >=
          ( numberOfJets
            + numberOfNegativeElectrons + numberOfPositiveElectrons
            + numberOfNegativeMuons + numberOfPositiveMuons ) )
      {
        return true;
      }
      else
      {
        return false;
      }
    }
    else
    {
      return false;
    }
  }

  void
  neutralinoToStauCascade::calculateAcceptance( acceptanceCutSet* const cuts,
                                    acceptanceValues* const currentAcceptance )
  // this returns the appropriate acceptances multiplied by branching ratios
  // from the electroweakino through the stau to the LSP.
  {
    currentAcceptance->setTwoJets( 0.0 );
    currentAcceptance->setOneJetOneNegativeElectron( 0.0 );
    currentAcceptance->setOneJetOnePositiveElectron( 0.0 );
    currentAcceptance->setOneJetOneNegativeMuon( 0.0 );
    currentAcceptance->setOneJetOnePositiveMuon( 0.0 );
    currentAcceptance->setElectronPlusAntielectron( 0.0 );
    currentAcceptance->setNegativeElectronPlusPositiveMuon( 0.0 );
    currentAcceptance->setNegativeMuonPlusPositiveElectron( 0.0 );
    currentAcceptance->setMuonPlusAntimuon( 0.0 );
    currentAcceptance->setOneJetZeroLeptons( 0.0 );
    currentAcceptance->setZeroJetsOneNegativeElectron( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveElectron( 0.0 );
    currentAcceptance->setZeroJetsOneNegativeMuon( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveMuon( 0.0 );
    currentAcceptance->setZeroJetsZeroLeptons( 0.0 );
    cascadeBr = ( first_BR->get_BR() * second_BR->get_BR() );
    if( LHC_FASER_global::negligible_BR < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      currentCuts = cuts;
      this->currentAcceptance = currentAcceptance;
      jetLeftHandedness
      = shortcut->quark_or_lepton_left_handedness(
                                                coloredDecayer->get_PDG_code(),
                                          electroweakDecayer->get_PDG_code() );
      jetRightHandedness = ( 1.0 - jetLeftHandedness );
      nearNegativeTauLeftHandedness
      = shortcut->quark_or_lepton_left_handedness(
                                            mediating_particle->get_PDG_code(),
                                              decaying_EWino->get_PDG_code() );
      nearNegativeTauRightHandedness = ( 1.0 - nearNegativeTauLeftHandedness );
      farNegativeTauLeftHandedness
      = shortcut->quark_or_lepton_left_handedness(
                                            mediating_particle->get_PDG_code(),
                              shortcut->get_neutralino_one()->get_PDG_code() );
      farNegativeTauRightHandedness = ( 1.0 - farNegativeTauLeftHandedness );


      /* in the following, XY(bar)Z(bar) is quark handedness, near tau
       * lepton handedness, far tau lepton handedness, YbarZ means that
       * the near tau lepton is positive while YZbar means that the far
       * tau lepton is positive.
       */
      // near lepton same handedness as quark:
      // hard near muons, hard far muons:
      currentNearMuonDistribution = nearSameHardMuonDistribution;
      currentNearPionDistribution = nearSameSoftPionDistribution;
      currentFarMuonDistribution = farSameHardMuonDistribution;
      currentFarPionDistribution = farSameSoftPionDistribution;
      // LLRbar
      negativeNearConfigurationBr = ( cascadeBr * jetLeftHandedness
                                                * nearNegativeTauLeftHandedness
                                              * farNegativeTauLeftHandedness );
      // RRbarL
      positiveNearConfigurationBr = ( cascadeBr * jetRightHandedness
                                                * nearNegativeTauLeftHandedness
                                              * farNegativeTauLeftHandedness );
      calculateForCurrentConfiguration();

      // hard near muons, soft far muons:
      currentNearMuonDistribution = nearSameHardMuonDistribution;
      currentNearPionDistribution = nearSameSoftPionDistribution;
      currentFarMuonDistribution = farSameSoftMuonDistribution;
      currentFarPionDistribution = farSameHardPionDistribution;
      // LLLbar
      negativeNearConfigurationBr = ( cascadeBr * jetLeftHandedness
                                                * nearNegativeTauLeftHandedness
                                             * farNegativeTauRightHandedness );
      // RRbarR
      positiveNearConfigurationBr = ( cascadeBr * jetRightHandedness
                                                * nearNegativeTauLeftHandedness
                                             * farNegativeTauRightHandedness );
      calculateForCurrentConfiguration();

      // soft near muons, hard far muons:
      currentNearMuonDistribution = nearSameSoftMuonDistribution;
      currentNearPionDistribution = nearSameHardPionDistribution;
      currentFarMuonDistribution = farSameHardMuonDistribution;
      currentFarPionDistribution = farSameSoftPionDistribution;
      // RRRbar
      negativeNearConfigurationBr = ( cascadeBr * jetRightHandedness
                                               * nearNegativeTauRightHandedness
                                              * farNegativeTauLeftHandedness );
      // LLbarL
      positiveNearConfigurationBr = ( cascadeBr * jetLeftHandedness
                                               * nearNegativeTauRightHandedness
                                              * farNegativeTauLeftHandedness );
      calculateForCurrentConfiguration();

      // soft near muons, soft far muons:
      currentNearMuonDistribution = nearSameSoftMuonDistribution;
      currentNearPionDistribution = nearSameHardPionDistribution;
      currentFarMuonDistribution = farSameSoftMuonDistribution;
      currentFarPionDistribution = farSameHardPionDistribution;
      // RRLbar
      negativeNearConfigurationBr = ( cascadeBr * jetRightHandedness
                                               * nearNegativeTauRightHandedness
                                             * farNegativeTauRightHandedness );
      // LLbarR
      positiveNearConfigurationBr = ( cascadeBr * jetLeftHandedness
                                               * nearNegativeTauRightHandedness
                                             * farNegativeTauRightHandedness );
      calculateForCurrentConfiguration();

      // near lepton opposite handedness to quark:
      // hard near muons, hard far muons:
      currentNearMuonDistribution = nearOppositeHardMuonDistribution;
      currentNearPionDistribution = nearOppositeSoftPionDistribution;
      currentFarMuonDistribution = farOppositeHardMuonDistribution;
      currentFarPionDistribution = farOppositeSoftPionDistribution;
      // RLRbar
      negativeNearConfigurationBr = ( cascadeBr * jetRightHandedness
                                                * nearNegativeTauLeftHandedness
                                              * farNegativeTauLeftHandedness );
      // LRbarL
      positiveNearConfigurationBr = ( cascadeBr * jetLeftHandedness
                                                * nearNegativeTauLeftHandedness
                                              * farNegativeTauLeftHandedness );
      calculateForCurrentConfiguration();

      // hard near muons, soft far muons:
      currentNearMuonDistribution = nearOppositeHardMuonDistribution;
      currentNearPionDistribution = nearOppositeSoftPionDistribution;
      currentFarMuonDistribution = farOppositeSoftMuonDistribution;
      currentFarPionDistribution = farOppositeHardPionDistribution;
      // RLLbar
      negativeNearConfigurationBr = ( cascadeBr * jetRightHandedness
                                                * nearNegativeTauLeftHandedness
                                             * farNegativeTauRightHandedness );
      // LRbarR
      positiveNearConfigurationBr = ( cascadeBr * jetLeftHandedness
                                                * nearNegativeTauLeftHandedness
                                             * farNegativeTauRightHandedness );
      calculateForCurrentConfiguration();

      // soft near muons, hard far muons:
      currentNearMuonDistribution = nearOppositeSoftMuonDistribution;
      currentNearPionDistribution = nearOppositeHardPionDistribution;
      currentFarMuonDistribution = farOppositeHardMuonDistribution;
      currentFarPionDistribution = farOppositeSoftPionDistribution;
      // LRRbar
      negativeNearConfigurationBr = ( cascadeBr * jetLeftHandedness
                                               * nearNegativeTauRightHandedness
                                              * farNegativeTauLeftHandedness );
      // RLbarL
      positiveNearConfigurationBr = ( cascadeBr * jetRightHandedness
                                               * nearNegativeTauRightHandedness
                                              * farNegativeTauLeftHandedness );
      calculateForCurrentConfiguration();

      // soft near muons, soft far muons:
      currentNearMuonDistribution = nearOppositeSoftMuonDistribution;
      currentNearPionDistribution = nearOppositeHardPionDistribution;
      currentFarMuonDistribution = farOppositeSoftMuonDistribution;
      currentFarPionDistribution = farOppositeHardPionDistribution;

      // LRLbar
      negativeNearConfigurationBr = ( cascadeBr * jetLeftHandedness
                                               * nearNegativeTauRightHandedness
                                             * farNegativeTauRightHandedness );
      // RLbarR
      positiveNearConfigurationBr = ( cascadeBr * jetRightHandedness
                                               * nearNegativeTauRightHandedness
                                             * farNegativeTauRightHandedness );
      calculateForCurrentConfiguration();
    }
    // debugging:
    /**std::cout
    << std::endl
    << "cascadeBr = " << cascadeBr;
    std::cout << std::endl;**/
  }

  void
  neutralinoToStauCascade::calculateForCurrentConfiguration()
  {
    nearMuonPass = integrateAcceptance( currentNearMuonDistribution,
                                        currentCuts->getPrimaryLeptonCut() );
    nearMuonFail = ( 1.0 - integrateAcceptance( currentNearMuonDistribution,
                                      currentCuts->getSecondaryLeptonCut() ) );
    nearPionPass = integrateAcceptance( currentNearPionDistribution,
                                        currentCuts->getJetCut() );
    nearPionFail = ( 1.0 - nearPionPass );
    farMuonPass = integrateAcceptance( currentFarMuonDistribution,
                                       currentCuts->getPrimaryLeptonCut() );
    farMuonFail = ( 1.0 - integrateAcceptance( currentFarMuonDistribution,
                                      currentCuts->getSecondaryLeptonCut() ) );
    farPionPass = integrateAcceptance( currentFarPionDistribution,
                                       currentCuts->getJetCut() );
    farPionFail = ( 1.0 - farPionPass );

    currentAcceptance->addToTwoJets( ( negativeNearConfigurationBr
                                       + positiveNearConfigurationBr )
                                     * tauPairToPionPairBr
                                     * nearPionPass * farPionPass );

    currentPass = ( negativeNearConfigurationBr * nearMuonPass * farPionPass
                  + positiveNearConfigurationBr * nearPionPass * farMuonPass );
    currentAcceptance->addToOneJetOneNegativeElectron(
                                 tauToPionTimesTauToElectronBr * currentPass );
    currentAcceptance->addToOneJetOneNegativeMuon( tauToPionTimesTauToMuonBr
                                                   * currentPass );

    currentPass = ( positiveNearConfigurationBr * nearMuonPass * farPionPass
                  + negativeNearConfigurationBr * nearPionPass * farMuonPass );
    currentAcceptance->addToOneJetOnePositiveElectron(
                                 tauToPionTimesTauToElectronBr * currentPass );
    currentAcceptance->addToOneJetOnePositiveMuon( tauToPionTimesTauToMuonBr
                                                   * currentPass );

    currentPass
    = ( ( negativeNearConfigurationBr + positiveNearConfigurationBr )
        * nearMuonPass * farMuonPass );
    currentAcceptance->addToElectronPlusAntielectron( tauPairToElectronPairBr
                                                      * currentPass );
    currentAcceptance->addToMuonPlusAntimuon( tauPairToMuonPairBr
                                              * currentPass );
    currentPass *= tauToElectronTimesTauToMuonBr;
    currentAcceptance->addToNegativeElectronPlusPositiveMuon( currentPass );
    currentAcceptance->addToNegativeMuonPlusPositiveElectron( currentPass );

    currentAcceptance->addToOneJetZeroLeptons( ( negativeNearConfigurationBr
                                                + positiveNearConfigurationBr )
                                               * ( tauPairToPionPairBr
                                                 * ( nearPionPass * farPionFail
                                                 + nearPionFail * farPionPass )
                                              + ( tauToPionTimesTauToElectronBr
                                                  + tauToPionTimesTauToMuonBr )
                                                 * ( nearPionPass * farMuonFail
                                            + nearMuonFail * farPionPass ) ) );

    currentNearFailSum
    = ( tauToPionTimesTauToElectronBr * nearPionFail
        + ( tauPairToElectronPairBr + tauToElectronTimesTauToMuonBr )
          * nearMuonFail );
    currentFarFailSum
    = ( tauToPionTimesTauToElectronBr * farPionFail
        + ( tauPairToElectronPairBr + tauToElectronTimesTauToMuonBr )
          * farMuonFail );
    currentAcceptance->addToZeroJetsOneNegativeElectron(
                                 ( negativeNearConfigurationBr * nearMuonPass )
                                 * currentFarFailSum
                                + ( positiveNearConfigurationBr * farMuonPass )
                                  * currentNearFailSum );
    currentAcceptance->addToZeroJetsOnePositiveElectron(
                                 ( positiveNearConfigurationBr * nearMuonPass )
                                 * currentFarFailSum
                                + ( negativeNearConfigurationBr * farMuonPass )
                                  * currentNearFailSum );

    currentNearFailSum
    = ( tauToPionTimesTauToMuonBr * nearPionFail
        + ( tauToElectronTimesTauToMuonBr + tauPairToMuonPairBr )
        * nearMuonFail );
    currentFarFailSum
    = ( tauToPionTimesTauToMuonBr * farPionFail
        + ( tauToElectronTimesTauToMuonBr + tauPairToMuonPairBr )
          * farMuonFail );
    currentAcceptance->addToZeroJetsOneNegativeMuon(
                                 ( negativeNearConfigurationBr * nearMuonPass )
                                 * currentFarFailSum
                                + ( positiveNearConfigurationBr * farMuonPass )
                                  * currentNearFailSum );
    currentAcceptance->addToZeroJetsOnePositiveMuon(
                                 ( positiveNearConfigurationBr * nearMuonPass )
                                 * currentFarFailSum
                                + ( negativeNearConfigurationBr * farMuonPass )
                                  * currentNearFailSum );

    currentAcceptance->addToZeroJetsZeroLeptons( ( negativeNearConfigurationBr
                                                + positiveNearConfigurationBr )
                                                 * ( tauPairToPionPairBr
                                                   * nearPionFail * farPionFail
                                              + ( tauToPionTimesTauToElectronBr
                                                  + tauToPionTimesTauToMuonBr )
                                                * ( nearPionFail * farMuonFail
                                                 + nearMuonFail * farPionFail )
                                                 + ( tauPairToElectronPairBr
                                          + 2.0 * tauToElectronTimesTauToMuonBr
                                                     + tauPairToMuonPairBr )
                                              * nearMuonFail * farMuonFail ) );
  }



  chargeSummedNeutralinoToStauCascade::chargeSummedNeutralinoToStauCascade(
                                leptonAcceptanceParameterSet* const kinematics,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                CppSLHA::particle_property_set const* const electroweakDecayer,
               CppSLHA::particle_property_set const* const intermediateDecayer,
                                        input_handler const* const shortcut ) :
    electroweakCascade( kinematics,
                        coloredDecayer,
                        electroweakDecayer,
                        intermediateDecayer,
                        false,
                        shortcut )
  {
    // we set up the BR of the (only-ish) channel:
    firstBr
    = shortcut->get_exclusive_BR_calculator( electroweakDecayer,
                                             intermediateDecayer,
                                             true,
                                             shortcut->get_empty_list() );
    secondBr
    = shortcut->get_exclusive_BR_calculator( intermediateDecayer,
                                             shortcut->get_neutralino_one(),
                                             true,
                                             shortcut->get_empty_list() );
    nearTauDistribution
    = new flat_near_muon_plus_antimuon( shortcut->get_readier(),
                                        shortcut->get_CppSLHA(),
                                        coloredDecayer,
                                        kinematics,
                                        electroweakDecayer,
                                        intermediateDecayer );
    farTauDistribution
    = new flat_far_muon_plus_antimuon( shortcut->get_readier(),
                                       shortcut->get_CppSLHA(),
                                       coloredDecayer,
                                       kinematics,
                                       electroweakDecayer,
                                       intermediateDecayer,
                                       shortcut->get_neutralino_one() );
    activeDistributions.push_back( nearTauDistribution );
    activeDistributions.push_back( farTauDistribution );

    nearHardMuonDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     nearTauDistribution,
                                     shortcut->get_hard_muon_from_tau() );
    nearSoftMuonDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     nearTauDistribution,
                                     shortcut->get_soft_muon_from_tau() );
    nearHardPionDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     nearTauDistribution,
                                     shortcut->get_hard_pion_from_tau() );
    nearSoftPionDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     nearTauDistribution,
                                     shortcut->get_soft_pion_from_tau() );
    farHardMuonDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     farTauDistribution,
                                     shortcut->get_hard_muon_from_tau() );
    farSoftMuonDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     farTauDistribution,
                                     shortcut->get_soft_muon_from_tau() );
    farHardPionDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     farTauDistribution,
                                     shortcut->get_hard_pion_from_tau() );
    farSoftPionDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     farTauDistribution,
                                     shortcut->get_soft_pion_from_tau() );
    activeDistributions.push_back( nearHardMuonDistribution );
    activeDistributions.push_back( nearSoftMuonDistribution );
    activeDistributions.push_back( nearHardPionDistribution );
    activeDistributions.push_back( nearSoftPionDistribution );
    activeDistributions.push_back( farHardMuonDistribution );
    activeDistributions.push_back( farSoftMuonDistribution );
    activeDistributions.push_back( farHardPionDistribution );
    activeDistributions.push_back( farSoftPionDistribution );
  }

  chargeSummedNeutralinoToStauCascade::~chargeSummedNeutralinoToStauCascade()
  {
    // does nothing.
  }


  void
  chargeSummedNeutralinoToStauCascade::calculateAcceptance(
                                                  acceptanceCutSet* const cuts,
                                    acceptanceValues* const currentAcceptance )
  // this returns the appropriate acceptances multiplied by branching ratios
  // from the electroweakino through the stau to the LSP.
  {
    currentAcceptance->setTwoJets( 0.0 );
    currentAcceptance->setOneJetOneNegativeElectron( 0.0 );
    currentAcceptance->setOneJetOnePositiveElectron( 0.0 );
    currentAcceptance->setOneJetOneNegativeMuon( 0.0 );
    currentAcceptance->setOneJetOnePositiveMuon( 0.0 );
    currentAcceptance->setElectronPlusAntielectron( 0.0 );
    currentAcceptance->setNegativeElectronPlusPositiveMuon( 0.0 );
    currentAcceptance->setNegativeMuonPlusPositiveElectron( 0.0 );
    currentAcceptance->setMuonPlusAntimuon( 0.0 );
    currentAcceptance->setOneJetZeroLeptons( 0.0 );
    currentAcceptance->setZeroJetsOneNegativeElectron( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveElectron( 0.0 );
    currentAcceptance->setZeroJetsOneNegativeMuon( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveMuon( 0.0 );
    currentAcceptance->setZeroJetsZeroLeptons( 0.0 );
    cascadeBr = ( first_BR->get_BR() * second_BR->get_BR() );
    if( LHC_FASER_global::negligible_BR < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      currentCuts = cuts;
      this->currentAcceptance = currentAcceptance;
      jetLeftHandedness
      = shortcut->quark_or_lepton_left_handedness(
                                                coloredDecayer->get_PDG_code(),
                                          electroweakDecayer->get_PDG_code() );
      jetRightHandedness = ( 1.0 - jetLeftHandedness );
      nearNegativeTauLeftHandedness
      = shortcut->quark_or_lepton_left_handedness(
                                            mediating_particle->get_PDG_code(),
                                              decaying_EWino->get_PDG_code() );
      nearNegativeTauRightHandedness = ( 1.0 - nearNegativeTauLeftHandedness );
      farNegativeTauLeftHandedness
      = shortcut->quark_or_lepton_left_handedness(
                                            mediating_particle->get_PDG_code(),
                              shortcut->get_neutralino_one()->get_PDG_code() );
      farNegativeTauRightHandedness = ( 1.0 - farNegativeTauLeftHandedness );

      /* 1st, left-handed quark, left-handed near negative tau lepton
       * with left-handed far positive tau lepton + right-handed near
       * positive tau lepton with right-handed far negative tau lepton,
       * plus the all-handedness-flipped versions:
       */
      currentNearMuonDistribution = nearHardMuonDistribution;
      currentNearPionDistribution = nearSoftPionDistribution;
      currentFarMuonDistribution = farSoftMuonDistribution;
      currentFarPionDistribution = farHardPionDistribution;
      configurationBr = ( cascadeBr * nearNegativeTauLeftHandedness
                                    * farNegativeTauRightHandedness );
      calculateForCurrentConfiguration();

      /* 2nd, left-handed quark, left-handed near negative tau lepton
       * with right-handed far positive tau lepton + right-handed near
       * positive tau lepton with left-handed far negative tau lepton,
       * plus the all-handedness-flipped versions:
       */
      currentNearMuonDistribution = nearHardMuonDistribution;
      currentNearPionDistribution = nearSoftPionDistribution;
      currentFarMuonDistribution = farHardMuonDistribution;
      currentFarPionDistribution = farSoftPionDistribution;
      configurationBr = ( cascadeBr * nearNegativeTauLeftHandedness
                                    * farNegativeTauLeftHandedness );
      calculateForCurrentConfiguration();

      /* 3rd, left-handed quark, right-handed near negative tau lepton
       * with left-handed far positive tau lepton + left-handed near
       * positive tau lepton with right-handed far negative tau lepton,
       * plus the all-handedness-flipped versions:
       */
      currentNearMuonDistribution = nearSoftMuonDistribution;
      currentNearPionDistribution = nearHardPionDistribution;
      currentFarMuonDistribution = farSoftMuonDistribution;
      currentFarPionDistribution = farHardPionDistribution;
      configurationBr = ( cascadeBr * nearNegativeTauRightHandedness
                                    * farNegativeTauRightHandedness );
      calculateForCurrentConfiguration();

      /* 4th, left-handed quark, right-handed near negative tau lepton
       * with right-handed far positive tau lepton + left-handed near
       * positive tau lepton with left-handed far negative tau lepton,
       * plus the all-handedness-flipped versions:
       */
      currentNearMuonDistribution = nearSoftMuonDistribution;
      currentNearPionDistribution = nearHardPionDistribution;
      currentFarMuonDistribution = farHardMuonDistribution;
      currentFarPionDistribution = farSoftPionDistribution;
      configurationBr = ( cascadeBr * nearNegativeTauRightHandedness
                                    * farNegativeTauLeftHandedness );
      calculateForCurrentConfiguration();
    }
    // debugging:
    /**std::cout
    << std::endl
    << "cascadeBr = " << cascadeBr;
    std::cout << std::endl;**/
  }

  void
  chargeSummedNeutralinoToStauCascade::calculateForCurrentConfiguration()
  {
    nearMuonPass = integrateAcceptance( currentNearMuonDistribution,
                                        currentCuts->getPrimaryLeptonCut() );
    nearMuonFail = ( 1.0 - integrateAcceptance( currentNearMuonDistribution,
                                      currentCuts->getSecondaryLeptonCut() ) );
    nearPionPass = integrateAcceptance( currentNearPionDistribution,
                                        currentCuts->getJetCut() );
    nearPionFail = ( 1.0 - nearPionPass );
    farMuonPass = integrateAcceptance( currentFarMuonDistribution,
                                       currentCuts->getPrimaryLeptonCut() );
    farMuonFail = ( 1.0 - integrateAcceptance( currentFarMuonDistribution,
                                      currentCuts->getSecondaryLeptonCut() ) );
    farPionPass = integrateAcceptance( currentFarPionDistribution,
                                       currentCuts->getJetCut() );
    farPionFail = ( 1.0 - farPionPass );

    currentAcceptance->addToTwoJets( 2.0 * configurationBr
                                         * tauPairToPionPairBr
                                         * nearPionPass * farPionPass );
    // the factor of 2.0 accounts for the charge-conjugate decay.

    currentPass = ( configurationBr * ( nearMuonPass * farPionPass
                                        + nearPionPass * farMuonPass ) );
    currentAcceptance->addToOneJetOneNegativeElectron(
                                 tauToPionTimesTauToElectronBr * currentPass );
    currentAcceptance->addToOneJetOnePositiveElectron(
                                 tauToPionTimesTauToElectronBr * currentPass );
    currentAcceptance->addToOneJetOnePositiveMuon(
                                     tauToPionTimesTauToMuonBr * currentPass );
    currentAcceptance->addToOneJetOneNegativeMuon(
                                     tauToPionTimesTauToMuonBr * currentPass );
    /* the charges are split evenly (chargeSummedNeutralinoToStauCascade should
     * only be used in situations where overall the charges are even, such as
     * in the case of a gluino decaying, which has equal branching ratio into
     * each charge version, so all the charge asymmetry cancels out, so we
     * might as well not bother calculating the asymmetry in the 1st place.
     */

    currentPass = ( 2.0 * configurationBr * nearMuonPass * farMuonPass );
    // the factor of 2.0 accounts for the charge-conjugate decay.
    currentAcceptance->addToElectronPlusAntielectron( tauPairToElectronPairBr
                                                      * currentPass );
    currentAcceptance->addToMuonPlusAntimuon( tauPairToMuonPairBr
                                              * currentPass );
    currentPass *= tauToElectronTimesTauToMuonBr;
    currentAcceptance->addToNegativeElectronPlusPositiveMuon( currentPass );
    currentAcceptance->addToNegativeMuonPlusPositiveElectron( currentPass );
    // as above, the charges are split evenly.

    currentAcceptance->addToOneJetZeroLeptons( 2.0 * configurationBr
                                                   * ( tauPairToPionPairBr
                                                 * ( nearPionPass * farPionFail
                                                 + nearPionFail * farPionPass )
                                              + ( tauToPionTimesTauToElectronBr
                                                  + tauToPionTimesTauToMuonBr )
                                                 * ( nearPionPass * farMuonFail
                                            + nearMuonFail * farPionPass ) ) );
    // the factor of 2.0 accounts for the charge-conjugate decay.

    currentPass
    = ( configurationBr
        * ( nearMuonPass
            * ( tauToPionTimesTauToElectronBr * nearPionFail
                + ( tauPairToElectronPairBr + tauToElectronTimesTauToMuonBr )
                  * nearMuonFail )
            + farMuonPass
              * ( tauToPionTimesTauToElectronBr * farPionFail
                  + ( tauPairToElectronPairBr + tauToElectronTimesTauToMuonBr )
                    * farMuonFail ) ) );
    currentAcceptance->addToZeroJetsOneNegativeElectron( currentPass );
    currentAcceptance->addToZeroJetsOnePositiveElectron( currentPass );
    currentPass
    = ( configurationBr
        * ( nearMuonPass
            * ( tauToPionTimesTauToMuonBr * nearPionFail
                + ( tauToElectronTimesTauToMuonBr + tauPairToMuonPairBr )
                  * nearMuonFail )
            + farMuonPass
              * ( tauToPionTimesTauToMuonBr * farPionFail
                  + ( tauToElectronTimesTauToMuonBr + tauPairToMuonPairBr )
                    * farMuonFail ) ) );
    currentAcceptance->addToZeroJetsOneNegativeMuon( currentPass );
    currentAcceptance->addToZeroJetsOnePositiveMuon( currentPass );

    currentAcceptance->addToZeroJetsZeroLeptons( configurationBr
                                                 * ( tauPairToPionPairBr
                                                   * nearPionFail * farPionFail
                                              + ( tauToPionTimesTauToElectronBr
                                                  + tauToPionTimesTauToMuonBr )
                                                 * ( nearPionFail * farMuonFail
                                                 + nearMuonFail * farPionFail )
                                                 + ( tauPairToElectronPairBr
                                          + 2.0 * tauToElectronTimesTauToMuonBr
                                                 + tauPairToMuonPairBr )
                                              * nearMuonFail * farMuonFail ) );
  }



  neutralinoToZCascade::neutralinoToZCascade(
                                leptonAcceptanceParameterSet* const kinematics,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                CppSLHA::particle_property_set const* const electroweakDecayer,
                                        input_handler const* const shortcut ) :
    electroweakCascade( kinematics,
                        coloredDecayer,
                        electroweakDecayer,
                        shortcut->get_Z(),
                        true,
                        shortcut )
  {
    double
    twiceSquareOfWeakCosine
    = ( 2.0 * shortcut->getWeakCosine() * shortcut->getWeakCosine() );
    negativeTauLeftHandedness
    = ( ( 1 - 2.0 * twiceSquareOfWeakCosine
          + twiceSquareOfWeakCosine * twiceSquareOfWeakCosine )
        / ( 5.0 - 2.0 * twiceSquareOfWeakCosine
            + 2.0 * twiceSquareOfWeakCosine * twiceSquareOfWeakCosine ) );
    // we set up the BR of the (only-ish) channel:
    firstBr
    = shortcut->get_exclusive_BR_calculator( electroweakDecayer,
                                             intermediateDecayer,
                                             true,
                                             shortcut->get_empty_list() );
    secondBr
    = shortcut->get_exclusive_BR_calculator( intermediateDecayer,
                                             shortcut->get_neutralino_one(),
                                             true,
                                             shortcut->get_empty_list() );
    directMuonDistribution = new Z_handed_muon( shortcut->get_readier(),
                                                shortcut->get_CppSLHA(),
                                                coloredDecayer,
                                                kinematics,
                                                electroweakDecayer,
                                                intermediateDecayer,
                                                shortcut->get_neutralino_one(),
                                                true,
                                                true );
    sameHandedTauDistribution = new Z_handed_muon( shortcut->get_readier(),
                                                   shortcut->get_CppSLHA(),
                                                   coloredDecayer,
                                                   kinematics,
                                                   electroweakDecayer,
                                                   intermediateDecayer,
                                                shortcut->get_neutralino_one(),
                                                   true,
                                                   false );
    oppositeHandedTauDistribution = new Z_handed_muon( shortcut->get_readier(),
                                                       shortcut->get_CppSLHA(),
                                                       coloredDecayer,
                                                       kinematics,
                                                       electroweakDecayer,
                                                       intermediateDecayer,
                                                shortcut->get_neutralino_one(),
                                                       false,
                                                       false );
    activeDistributions.push_back( directMuonDistribution );
    activeDistributions.push_back( sameHandedTauDistribution );
    activeDistributions.push_back( oppositeHandedTauDistribution );

    sameHardMuonDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     sameHandedTauDistribution,
                                     shortcut->get_hard_muon_from_tau() );
    sameSoftMuonDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     sameHandedTauDistribution,
                                     shortcut->get_soft_muon_from_tau() );
    sameHardPionDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     sameHandedTauDistribution,
                                     shortcut->get_hard_pion_from_tau() );
    sameSoftPionDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     sameHandedTauDistribution,
                                     shortcut->get_soft_pion_from_tau() );
    oppositeHardMuonDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     oppositeHandedTauDistribution,
                                     shortcut->get_hard_muon_from_tau() );
    oppositeSoftMuonDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     oppositeHandedTauDistribution,
                                     shortcut->get_soft_muon_from_tau() );
    oppositeHardPionDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     oppositeHandedTauDistribution,
                                     shortcut->get_hard_pion_from_tau() );
    oppositeSoftPionDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     oppositeHandedTauDistribution,
                                     shortcut->get_soft_pion_from_tau() );
    activeDistributions.push_back( sameHardMuonDistribution );
    activeDistributions.push_back( sameSoftMuonDistribution );
    activeDistributions.push_back( sameHardPionDistribution );
    activeDistributions.push_back( sameSoftPionDistribution );
    activeDistributions.push_back( oppositeHardMuonDistribution );
    activeDistributions.push_back( oppositeSoftMuonDistribution );
    activeDistributions.push_back( oppositeHardPionDistribution );
    activeDistributions.push_back( oppositeSoftPionDistribution );
  }

  neutralinoToZCascade::~neutralinoToZCascade()
  {
    // does nothing.
  }


  bool
  neutralinoToZCascade::validSignal( int const numberOfJets,
                                     int const numberOfNegativeElectrons,
                                     int const numberOfPositiveElectrons,
                                     int const numberOfNegativeMuons,
                                     int const numberOfPositiveMuons )
  /* this returns true if a configuration where each of the signs of tau lepton
   * decayed either into a detected jet, detected lepton, or undetected
   * particle, & false otherwise.
   */
  {
    if( 2 == numberOfJets )
    {
      if( ( 0 == numberOfNegativeElectrons )
          &&
          ( 0 == numberOfPositiveElectrons )
          &&
          ( 0 == numberOfNegativeMuons )
          &&
          ( 0 == numberOfPositiveMuons ) )
      {
        return true;
      }
      else
      {
        return false;
      }
    }
    else if( ( 0 <= numberOfJets )
             &&
             ( 0 <= numberOfNegativeElectrons )
             &&
             ( 0 <= numberOfPositiveElectrons )
             &&
             ( 0 <= numberOfNegativeMuons )
             &&
             ( 0 <= numberOfPositiveMuons )
             &&
             ( 1 >= ( numberOfNegativeElectrons + numberOfNegativeMuons ) )
             &&
             ( 1 >= ( numberOfPositiveElectrons + numberOfPositiveMuons ) ) )
    {
      if( 2 >=
          ( numberOfJets
            + numberOfNegativeElectrons + numberOfPositiveElectrons
            + numberOfNegativeMuons + numberOfPositiveMuons ) )
      {
        return true;
      }
      else
      {
        return false;
      }
    }
    else
    {
      return false;
    }
  }

  void
  neutralinoToZCascade::calculateAcceptance( acceptanceCutSet* const cuts,
                                    acceptanceValues* const currentAcceptance )
  // this returns the appropriate acceptances multiplied by branching ratios
  // from the electroweakino through the stau to the LSP.
  {
    currentAcceptance->setTwoJets( 0.0 );
    currentAcceptance->setOneJetOneNegativeElectron( 0.0 );
    currentAcceptance->setOneJetOnePositiveElectron( 0.0 );
    currentAcceptance->setOneJetOneNegativeMuon( 0.0 );
    currentAcceptance->setOneJetOnePositiveMuon( 0.0 );
    currentAcceptance->setElectronPlusAntielectron( 0.0 );
    currentAcceptance->setNegativeElectronPlusPositiveMuon( 0.0 );
    currentAcceptance->setNegativeMuonPlusPositiveElectron( 0.0 );
    currentAcceptance->setMuonPlusAntimuon( 0.0 );
    currentAcceptance->setOneJetZeroLeptons( 0.0 );
    currentAcceptance->setZeroJetsOneNegativeElectron( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveElectron( 0.0 );
    currentAcceptance->setZeroJetsOneNegativeMuon( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveMuon( 0.0 );
    currentAcceptance->setZeroJetsZeroLeptons( 0.0 );
    cascadeBr = ( first_BR->get_BR() * second_BR->get_BR() );
    if( LHC_FASER_global::negligible_BR < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      /* at some point, it'd be nice to break this down into transverse
       * contributions & longitudinal contributions, since actually the
       * Goldstone boson contribution should both be large (only needs the
       * Higgsino component of just 1 of the neutralinos, rather than both) &
       * has a different (& easier to calculate) distribution.
       *
       * however, I think that for LHC-relevant energies, the Higgsino
       * components of all the neutralinos is still going to be larger than
       * m_tau / m_W, enough that I'm going to ignore the longitudinal
       * contributions.
       */

      directMuonPass = integrateAcceptance( directMuonDistribution,
                                            cuts->getPrimaryLeptonCut() );
      directMuonFail = ( 1.0 - integrateAcceptance( directMuonDistribution,
                                             cuts->getSecondaryLeptonCut() ) );
      directJetPass = integrateAcceptance( directMuonDistribution,
                                           cuts->getJetCut() );
      directJetFail = ( 1.0 - directJetPass );

      configurationBr
      = ( cascadeBr * CppSLHA::PDG_data::Z_to_hadrons_BR * directJetPass );
      currentAcceptance->addToTwoJets( configurationBr * directJetPass );
      currentAcceptance->addToOneJetZeroLeptons( 2.0 * configurationBr
                                                     * directJetFail );
      currentAcceptance->addToOssfMinusOsdf( cascadeBr
                           * ( CppSLHA::PDG_data::Z_to_electron_antielectron_BR
                               + CppSLHA::PDG_data::Z_to_muon_antimuon_BR )
                             * directMuonPass * directMuonPass );
      configurationBr
      = ( cascadeBr * CppSLHA::PDG_data::Z_to_electron_antielectron_BR
          * directMuonPass );
      currentAcceptance->addToElectronPlusAntielectron( configurationBr
                                                        * directMuonPass );
      configurationBr *= directMuonFail;
      currentAcceptance->addToZeroJetsOneNegativeElectron( configurationBr );
      currentAcceptance->addToZeroJetsOnePositiveElectron( configurationBr );
      configurationBr
      = ( cascadeBr * CppSLHA::PDG_data::Z_to_muon_antimuon_BR
          * directMuonPass );
      currentAcceptance->addToMuonPlusAntimuon( configurationBr
                                                * directMuonPass );
      configurationBr *= directMuonFail;
      currentAcceptance->addToZeroJetsOneNegativeMuon( configurationBr );
      currentAcceptance->addToZeroJetsOnePositiveMuon( configurationBr );
      currentAcceptance->addToZeroJetsZeroLeptons( cascadeBr
                                         * ( CppSLHA::PDG_data::Z_to_hadrons_BR
                                             * directJetFail * directJetFail
                           + ( CppSLHA::PDG_data::Z_to_electron_antielectron_BR
                                   + CppSLHA::PDG_data::Z_to_muon_antimuon_BR )
                             * directMuonFail * directMuonFail ) );

      // now Z decays to tau-antitau pairs, followed by the decays of the taus:
      cascadeBr *= CppSLHA::PDG_data::Z_to_tau_lepton_antilepton_BR;
      jetLeftHandedness
      = shortcut->quark_or_lepton_left_handedness(
                                                coloredDecayer->get_PDG_code(),
                                          electroweakDecayer->get_PDG_code() );
      // left-handed quark, left-handed negative tau lepton:
      configurationBr
      = ( cascadeBr * jetLeftHandedness * negativeTauLeftHandedness );
      currentNegativeMuonDistribution = sameHardMuonDistribution;
      currentNegativePionDistribution = sameSoftPionDistribution;
      currentPositiveMuonDistribution = oppositeHardMuonDistribution;
      currentPositivePionDistribution = oppositeSoftPionDistribution;
      calculateForCurrentConfiguration();

      // left-handed quark, right-handed negative tau lepton:
      configurationBr
      = ( cascadeBr * jetLeftHandedness
                    * ( 1.0 - negativeTauLeftHandedness ) );
      currentNegativeMuonDistribution = oppositeSoftMuonDistribution;
      currentNegativePionDistribution = oppositeHardPionDistribution;
      currentPositiveMuonDistribution = sameSoftMuonDistribution;
      currentPositivePionDistribution = sameHardPionDistribution;
      calculateForCurrentConfiguration();

      // right-handed quark, left-handed negative tau lepton:
      configurationBr
      = ( cascadeBr * ( 1.0 - jetLeftHandedness )
          * negativeTauLeftHandedness );
      currentNegativeMuonDistribution = oppositeHardMuonDistribution;
      currentNegativePionDistribution = oppositeSoftPionDistribution;
      currentPositiveMuonDistribution = sameHardMuonDistribution;
      currentPositivePionDistribution = sameSoftPionDistribution;
      calculateForCurrentConfiguration();

      // right-handed quark, right-handed negative tau lepton:
      configurationBr
      = ( cascadeBr * ( 1.0 - jetLeftHandedness )
                    * ( 1.0 - negativeTauLeftHandedness ) );
      currentNegativeMuonDistribution = sameSoftMuonDistribution;
      currentNegativePionDistribution = sameHardPionDistribution;
      currentPositiveMuonDistribution = oppositeSoftMuonDistribution;
      currentPositivePionDistribution = oppositeHardPionDistribution;
      calculateForCurrentConfiguration();
    }
    // debugging:
    /**std::cout
    << std::endl
    << "cascadeBr = " << cascadeBr;
    std::cout << std::endl;**/
  }

  void
  neutralinoToZCascade::calculateForCurrentConfiguration()
  {
    negativeTauMuonPass = integrateAcceptance( currentNegativeMuonDistribution,
                                          currentCuts->getPrimaryLeptonCut() );
    negativeTauMuonFail
    = ( 1.0 - integrateAcceptance( currentNegativeMuonDistribution,
                                   currentCuts->getSecondaryLeptonCut() ) );
    negativeTauPionPass = integrateAcceptance( currentNegativePionDistribution,
                                               currentCuts->getJetCut() );
    negativeTauPionFail = ( 1.0 - negativeTauPionPass );
    positiveTauMuonPass = integrateAcceptance( currentPositiveMuonDistribution,
                                          currentCuts->getPrimaryLeptonCut() );
    positiveTauMuonFail
    = ( 1.0 - integrateAcceptance( currentPositiveMuonDistribution,
                                   currentCuts->getSecondaryLeptonCut() ) );
    positiveTauPionPass = integrateAcceptance( currentPositivePionDistribution,
                                               currentCuts->getJetCut() );
    positiveTauPionFail = ( 1.0 - positiveTauPionPass );

    currentAcceptance->addToTwoJets( configurationBr
                                     * tauPairToPionPairBr
                                     * negativeTauPionPass
                                     * positiveTauPionPass );
    currentPass = ( configurationBr * tauToPionTimesTauToElectronBr );
    currentAcceptance->addToOneJetOneNegativeElectron( currentPass
                                                       * positiveTauPionPass
                                                       * negativeTauMuonPass );
    currentAcceptance->addToOneJetOnePositiveElectron( currentPass
                                                       * negativeTauPionPass
                                                       * positiveTauMuonPass );
    currentPass = ( configurationBr * tauToPionTimesTauToMuonBr );
    currentAcceptance->addToOneJetOneNegativeMuon( currentPass
                                                   * positiveTauPionPass
                                                   * negativeTauMuonPass );
    currentAcceptance->addToOneJetOnePositiveMuon( currentPass
                                                   * negativeTauPionPass
                                                   * positiveTauMuonPass );
    currentPass
    = ( configurationBr * negativeTauMuonPass * positiveTauMuonPass );
    currentAcceptance->addToElectronPlusAntielectron( tauPairToElectronPairBr
                                                      * currentPass );
    currentAcceptance->addToMuonPlusAntimuon( tauPairToMuonPairBr
                                              * currentPass );
    currentPass *= tauToElectronTimesTauToMuonBr;
    currentAcceptance->addToNegativeElectronPlusPositiveMuon( currentPass );
    currentAcceptance->addToNegativeMuonPlusPositiveElectron( currentPass );

    currentAcceptance->addToOneJetZeroLeptons( configurationBr
                                               * ( tauPairToPionPairBr
                                                   * ( negativeTauPionPass
                                                       * positiveTauPionFail
                                                       + positiveTauPionPass
                                                        * negativeTauPionFail )
                                              + ( tauToPionTimesTauToElectronBr
                                                  + tauToPionTimesTauToMuonBr )
                                                * ( negativeTauPionPass
                                                    * positiveTauMuonFail
                                                    + positiveTauPionPass
                                                   * negativeTauMuonFail ) ) );

    currentAcceptance->addToZeroJetsOneNegativeElectron( configurationBr
                                                         * negativeTauMuonPass
                                              * ( tauToPionTimesTauToElectronBr
                                                  * positiveTauPionFail
                                                  + ( tauPairToElectronPairBr
                                              + tauToElectronTimesTauToMuonBr )
                                                     * positiveTauMuonFail ) );
    currentAcceptance->addToZeroJetsOnePositiveElectron( configurationBr
                                                         * positiveTauMuonPass
                                              * ( tauToPionTimesTauToElectronBr
                                                  * negativeTauPionFail
                                                  + ( tauPairToElectronPairBr
                                              + tauToElectronTimesTauToMuonBr )
                                                     * negativeTauMuonFail ) );
    currentAcceptance->addToZeroJetsOneNegativeMuon( configurationBr
                                                     * negativeTauMuonPass
                                                  * ( tauToPionTimesTauToMuonBr
                                                      * positiveTauPionFail
                                              + ( tauToElectronTimesTauToMuonBr
                                                  + tauPairToMuonPairBr )
                                                     * positiveTauMuonFail ) );
    currentAcceptance->addToZeroJetsOnePositiveMuon( configurationBr
                                                     * positiveTauMuonPass
                                                  * ( tauToPionTimesTauToMuonBr
                                                      * negativeTauPionFail
                                              + ( tauToElectronTimesTauToMuonBr
                                                  + tauPairToMuonPairBr )
                                                     * negativeTauMuonFail ) );

    currentAcceptance->addToZeroJetsZeroLeptons( configurationBr
                                                 * ( tauPairToPionPairBr
                                                     * negativeTauPionFail
                                                     * positiveTauPionFail
                                              + ( tauToPionTimesTauToElectronBr
                                                  + tauToPionTimesTauToMuonBr )
                                                  * ( negativeTauPionFail
                                                      * positiveTauMuonFail
                                                        + positiveTauPionFail
                                                        * negativeTauMuonFail )
                                                 + ( tauPairToElectronPairBr
                                          + 2.0 * tauToElectronTimesTauToMuonBr
                                                     + tauPairToMuonPairBr )
                                                   * negativeTauMuonFail
                                                   * negativeTauMuonFail ) );
  }



  chargeSummedNeutralinoToZCascade::chargeSummedNeutralinoToZCascade(
                                leptonAcceptanceParameterSet* const kinematics,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                CppSLHA::particle_property_set const* const electroweakDecayer,
                                        input_handler const* const shortcut ) :
    electroweakCascade( kinematics,
                        coloredDecayer,
                        electroweakDecayer,
                        shortcut->get_Z(),
                        true,
                        shortcut )
  {
    double
    twiceSquareOfWeakCosine
    = ( 2.0 * shortcut->getWeakCosine() * shortcut->getWeakCosine() );
    negativeTauLeftHandedness
    = ( ( 1 - 2.0 * twiceSquareOfWeakCosine
          + twiceSquareOfWeakCosine * twiceSquareOfWeakCosine )
        / ( 5.0 - 2.0 * twiceSquareOfWeakCosine
            + 2.0 * twiceSquareOfWeakCosine * twiceSquareOfWeakCosine ) );
    // we set up the BR of the (only-ish) channel:
    firstBr
    = shortcut->get_exclusive_BR_calculator( electroweakDecayer,
                                             intermediateDecayer,
                                             true,
                                             shortcut->get_empty_list() );
    secondBr
    = shortcut->get_exclusive_BR_calculator( intermediateDecayer,
                                             shortcut->get_neutralino_one(),
                                             true,
                                             shortcut->get_empty_list() );
    directMuonDistribution = new Z_handed_muon( shortcut->get_readier(),
                                                shortcut->get_CppSLHA(),
                                                coloredDecayer,
                                                kinematics,
                                                electroweakDecayer,
                                                intermediateDecayer,
                                                shortcut->get_neutralino_one(),
                                                true,
                                                true );
    activeDistributions.push_back( directMuonDistribution );

    hardMuonDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     directMuonDistribution,
                                     shortcut->get_hard_muon_from_tau() );
    softMuonDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     directMuonDistribution,
                                     shortcut->get_soft_muon_from_tau() );
    hardPionDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     directMuonDistribution,
                                     shortcut->get_hard_pion_from_tau() );
    softPionDistribution
    = new visible_tau_decay_product( shortcut->get_readier(),
                                     directMuonDistribution,
                                     shortcut->get_soft_pion_from_tau() );
    activeDistributions.push_back( hardMuonDistribution );
    activeDistributions.push_back( softMuonDistribution );
    activeDistributions.push_back( hardPionDistribution );
    activeDistributions.push_back( softPionDistribution );
  }

  chargeSummedNeutralinoToZCascade::~chargeSummedNeutralinoToZCascade()
  {
    // does nothing.
  }


  bool
  chargeSummedNeutralinoToZCascade::validSignal( int const numberOfJets,
                                           int const numberOfNegativeElectrons,
                                           int const numberOfPositiveElectrons,
                                               int const numberOfNegativeMuons,
                                              int const numberOfPositiveMuons )
  /* this returns true if a configuration where each of the signs of tau lepton
   * decayed either into a detected jet, detected lepton, or undetected
   * particle, & false otherwise.
   */
  {
    if( 2 == numberOfJets )
    {
      if( ( 0 == numberOfNegativeElectrons )
          &&
          ( 0 == numberOfPositiveElectrons )
          &&
          ( 0 == numberOfNegativeMuons )
          &&
          ( 0 == numberOfPositiveMuons ) )
      {
        return true;
      }
      else
      {
        return false;
      }
    }
    else if( ( 0 <= numberOfJets )
             &&
             ( 0 <= numberOfNegativeElectrons )
             &&
             ( 0 <= numberOfPositiveElectrons )
             &&
             ( 0 <= numberOfNegativeMuons )
             &&
             ( 0 <= numberOfPositiveMuons )
             &&
             ( 1 >= ( numberOfNegativeElectrons + numberOfNegativeMuons ) )
             &&
             ( 1 >= ( numberOfPositiveElectrons + numberOfPositiveMuons ) ) )
    {
      if( 2 >=
          ( numberOfJets
            + numberOfNegativeElectrons + numberOfPositiveElectrons
            + numberOfNegativeMuons + numberOfPositiveMuons ) )
      {
        return true;
      }
      else
      {
        return false;
      }
    }
    else
    {
      return false;
    }
  }

  void
  chargeSummedNeutralinoToZCascade::calculateAcceptance(
                                                  acceptanceCutSet* const cuts,
                                    acceptanceValues* const currentAcceptance )
  // this returns the appropriate acceptances multiplied by branching ratios
  // from the electroweakino through the stau to the LSP.
  {
    currentAcceptance->setTwoJets( 0.0 );
    currentAcceptance->setOneJetOneNegativeElectron( 0.0 );
    currentAcceptance->setOneJetOnePositiveElectron( 0.0 );
    currentAcceptance->setOneJetOneNegativeMuon( 0.0 );
    currentAcceptance->setOneJetOnePositiveMuon( 0.0 );
    currentAcceptance->setElectronPlusAntielectron( 0.0 );
    currentAcceptance->setNegativeElectronPlusPositiveMuon( 0.0 );
    currentAcceptance->setNegativeMuonPlusPositiveElectron( 0.0 );
    currentAcceptance->setMuonPlusAntimuon( 0.0 );
    currentAcceptance->setOneJetZeroLeptons( 0.0 );
    currentAcceptance->setZeroJetsOneNegativeElectron( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveElectron( 0.0 );
    currentAcceptance->setZeroJetsOneNegativeMuon( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveMuon( 0.0 );
    currentAcceptance->setZeroJetsZeroLeptons( 0.0 );
    cascadeBr = ( first_BR->get_BR() * second_BR->get_BR() );
    if( LHC_FASER_global::negligible_BR < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      /* at some point, it'd be nice to break this down into transverse
       * contributions & longitudinal contributions, since actually the
       * Goldstone boson contribution should both be large (only needs the
       * Higgsino component of just 1 of the neutralinos, rather than both) &
       * has a different (& easier to calculate) distribution.
       *
       * however, I think that for LHC-relevant energies, the Higgsino
       * components of all the neutralinos is still going to be larger than
       * m_tau / m_W, enough that I'm going to ignore the longitudinal
       * contributions.
       */

      directMuonPass = integrateAcceptance( directMuonDistribution,
                                            cuts->getPrimaryLeptonCut() );
      directMuonFail = ( 1.0 - integrateAcceptance( directMuonDistribution,
                                             cuts->getSecondaryLeptonCut() ) );
      directJetPass = integrateAcceptance( directMuonDistribution,
                                           cuts->getJetCut() );
      directJetFail = ( 1.0 - directJetPass );

      configurationBr
      = ( cascadeBr * CppSLHA::PDG_data::Z_to_hadrons_BR * directJetPass );
      currentAcceptance->addToTwoJets( configurationBr * directJetPass );
      currentAcceptance->addToOneJetZeroLeptons( 2.0 * configurationBr
                                                     * directJetFail );
      currentAcceptance->addToOssfMinusOsdf( cascadeBr
                           * ( CppSLHA::PDG_data::Z_to_electron_antielectron_BR
                               + CppSLHA::PDG_data::Z_to_muon_antimuon_BR )
                             * directMuonPass * directMuonPass );
      configurationBr
      = ( cascadeBr * CppSLHA::PDG_data::Z_to_electron_antielectron_BR
          * directMuonPass );
      currentAcceptance->addToElectronPlusAntielectron( configurationBr
                                                        * directMuonPass );
      configurationBr *= directMuonFail;
      currentAcceptance->addToZeroJetsOneNegativeElectron( configurationBr );
      currentAcceptance->addToZeroJetsOnePositiveElectron( configurationBr );
      configurationBr
      = ( cascadeBr * CppSLHA::PDG_data::Z_to_muon_antimuon_BR
          * directMuonPass );
      currentAcceptance->addToMuonPlusAntimuon( configurationBr
                                                * directMuonPass );
      configurationBr *= directMuonFail;
      currentAcceptance->addToZeroJetsOneNegativeMuon( configurationBr );
      currentAcceptance->addToZeroJetsOnePositiveMuon( configurationBr );
      currentAcceptance->addToZeroJetsZeroLeptons( cascadeBr
                                         * ( CppSLHA::PDG_data::Z_to_hadrons_BR
                                             * directJetFail * directJetFail
                           + ( CppSLHA::PDG_data::Z_to_electron_antielectron_BR
                                   + CppSLHA::PDG_data::Z_to_muon_antimuon_BR )
                             * directMuonFail * directMuonFail ) );

      // now Z decays to tau-antitau pairs, followed by the decays of the taus:
      cascadeBr *= CppSLHA::PDG_data::Z_to_tau_lepton_antilepton_BR;
      // left-handed negative tau lepton:
      configurationBr
      = ( cascadeBr * negativeTauLeftHandedness );
      currentMuonDistribution = hardMuonDistribution;
      currentPionDistribution = softPionDistribution;
      calculateForCurrentConfiguration();

      // right-handed negative tau lepton:
      configurationBr
      = ( cascadeBr * ( 1.0 - negativeTauLeftHandedness ) );
      currentMuonDistribution = softMuonDistribution;
      currentPionDistribution = hardPionDistribution;
      calculateForCurrentConfiguration();
    }
    // debugging:
    /**std::cout
    << std::endl
    << "cascadeBr = " << cascadeBr;
    std::cout << std::endl;**/
  }

  void
  chargeSummedNeutralinoToZCascade::calculateForCurrentConfiguration()
  {
    tauMuonPass = integrateAcceptance( currentMuonDistribution,
                                       currentCuts->getPrimaryLeptonCut() );
    tauMuonFail = ( 1.0 - integrateAcceptance( currentMuonDistribution,
                                      currentCuts->getSecondaryLeptonCut() ) );
    tauPionPass = integrateAcceptance( currentPionDistribution,
                                       currentCuts->getJetCut() );
    tauPionFail = ( 1.0 - tauPionPass );

    currentAcceptance->addToTwoJets( configurationBr
                                     * tauPairToPionPairBr
                                     * tauPionPass * tauPionPass );
    currentPass = ( configurationBr * tauToPionTimesTauToElectronBr
                    * tauPionPass * tauMuonPass );
    currentAcceptance->addToOneJetOneNegativeElectron( currentPass );
    currentAcceptance->addToOneJetOnePositiveElectron( currentPass );
    currentPass = ( configurationBr * tauToPionTimesTauToMuonBr
                    * tauPionPass * tauMuonPass );
    currentAcceptance->addToOneJetOneNegativeMuon( currentPass );
    currentAcceptance->addToOneJetOnePositiveMuon( currentPass );
    currentPass = ( configurationBr * tauMuonPass * tauMuonPass );
    currentAcceptance->addToElectronPlusAntielectron( tauPairToElectronPairBr
                                                      * currentPass );
    currentAcceptance->addToMuonPlusAntimuon( tauPairToMuonPairBr
                                              * currentPass );
    currentPass *= tauToElectronTimesTauToMuonBr;
    currentAcceptance->addToNegativeElectronPlusPositiveMuon( currentPass );
    currentAcceptance->addToNegativeMuonPlusPositiveElectron( currentPass );

    currentAcceptance->addToOneJetZeroLeptons( 2.0 * configurationBr
                                                   * tauPionPass
                                                   * ( tauPairToPionPairBr
                                                       * tauPionFail
                                              + ( tauToPionTimesTauToElectronBr
                                                  + tauToPionTimesTauToMuonBr )
                                                * tauMuonFail ) );
    currentPass = ( configurationBr * tauMuonPass
                    * ( tauToPionTimesTauToElectronBr * tauPionFail
                  + ( tauPairToElectronPairBr + tauToElectronTimesTauToMuonBr )
                    * tauMuonFail ) );
    currentAcceptance->addToZeroJetsOneNegativeElectron( currentPass );
    currentAcceptance->addToZeroJetsOnePositiveElectron( currentPass );
    currentPass = ( configurationBr * tauMuonPass
                    * ( tauToPionTimesTauToMuonBr * tauPionFail
                      + ( tauToElectronTimesTauToMuonBr + tauPairToMuonPairBr )
                        * tauMuonFail ) );
    currentAcceptance->addToZeroJetsOneNegativeMuon( currentPass );
    currentAcceptance->addToZeroJetsOnePositiveMuon( currentPass );

    currentAcceptance->addToZeroJetsZeroLeptons( configurationBr
                                                 * ( tauPairToPionPairBr
                                                   * nearPionFail * farPionFail
                                        + 2.0 * ( tauToPionTimesTauToElectronBr
                                                  + tauToPionTimesTauToMuonBr )
                                                  * tauPionFail * tauMuonFail
                                                 + ( tauPairToElectronPairBr
                                          + 2.0 * tauToElectronTimesTauToMuonBr
                                                     + tauPairToMuonPairBr )
                                              * nearMuonFail * farMuonFail ) );
  }
}  // end of LHC_FASER namespace.
