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
                                         leptonAcceptanceCut* const kinematics,
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
    nearSameDistribution( NULL ),
    farSameDistribution( NULL ),
    nearOppositeDistribution( NULL ),
    farOppositeDistribution( NULL )
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
                        kinematics->getAcceptance()->getEffectiveSquarkMass(),
                                    electroweakDecayer,
                                    intermediateDecayer );
    nearOppositeDistribution
    = new opposite_chirality_near_muon( shortcut->get_readier(),
                                        shortcut->get_CppSLHA(),
                                        coloredDecayer,
                        kinematics->getAcceptance()->getEffectiveSquarkMass(),
                                        electroweakDecayer,
                                        intermediateDecayer );
    farSameDistribution
    = new same_chirality_far_muon( shortcut->get_readier(),
                                   shortcut->get_CppSLHA(),
                                   coloredDecayer,
                        kinematics->getAcceptance()->getEffectiveSquarkMass(),
                                   electroweakDecayer,
                                   intermediateDecayer,
                                   shortcut->get_neutralino_one() );
    farOppositeDistribution
    = new opposite_chirality_far_muon( shortcut->get_readier(),
                                       shortcut->get_CppSLHA(),
                                       coloredDecayer,
                        kinematics->getAcceptance()->getEffectiveSquarkMass(),
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
    << "cascade_BR = " << cascade_BR;
    std::cout << std::endl;**/
  }



  chargeSummedNeutralinoToSemuCascade::chargeSummedNeutralinoToSemuCascade(
                                         leptonAcceptanceCut* const kinematics,
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



  neutralino_to_stau_cascade::neutralino_to_stau_cascade(
                                    readier_for_new_point* const given_readier,
                                                 bool const should_sum_charges,
                                                double const given_cuts->getPrimaryLeptonCut(),
                                              double const given_cuts->getSecondaryLeptonCut(),
                               lepton_acceptance_value* const given_kinematics,
           CppSLHA::particle_property_set const* const given_decaying_scolored,
                                 bool const given_scolored_is_not_antiparticle,
              CppSLHA::particle_property_set const* const given_decaying_EWino,
                                    bool const given_EWino_is_not_antiparticle,
          CppSLHA::particle_property_set const* const given_mediating_particle,
                                    input_handler const* const given_shortcuts,
                                       double* const given_two_jets_no_leptons,
                             double* const given_one_jet_one_negative_electron,
                             double* const given_one_jet_one_positive_electron,
                                 double* const given_one_jet_one_negative_muon,
                                 double* const given_one_jet_one_positive_muon,
                                     double* const given_no_jets_two_electrons,
                                         double* const given_no_jets_two_muons,
           double* const given_no_jets_one_positive_muon_one_negative_electron,
           double* const given_no_jets_one_negative_muon_one_positive_electron,
                                        double* const given_one_jet_no_leptons,
                             double* const given_no_jets_one_negative_electron,
                             double* const given_no_jets_one_positive_electron,
                                 double* const given_no_jets_one_negative_muon,
                                 double* const given_no_jets_one_positive_muon,
                                     double* const given_no_jets_no_leptons ) :
    channel_calculator( given_readier,
                        given_cuts->getPrimaryLeptonCut(),
                        given_cuts->getSecondaryLeptonCut(),
                        -1.0,
                        given_kinematics,
                        given_decaying_scolored,
                        given_scolored_is_not_antiparticle,
                        given_decaying_EWino,
                        given_EWino_is_not_antiparticle,
                        given_mediating_particle,
                        given_shortcuts ),
    should_sum_charges_flag( should_sum_charges ),
    near_same_distribution( NULL ),
    far_same_distribution( NULL ),
    near_opposite_distribution( NULL ),
    far_opposite_distribution( NULL ),
    near_same_hard_muon_distribution( NULL ),
    near_same_soft_muon_distribution( NULL ),
    near_same_hard_pion_distribution( NULL ),
    near_same_soft_pion_distribution( NULL ),
    near_opposite_hard_muon_distribution( NULL ),
    near_opposite_soft_muon_distribution( NULL ),
    near_opposite_hard_pion_distribution( NULL ),
    near_opposite_soft_pion_distribution( NULL ),
    far_same_hard_muon_distribution( NULL ),
    far_same_soft_muon_distribution( NULL ),
    far_same_hard_pion_distribution( NULL ),
    far_same_soft_pion_distribution( NULL ),
    far_opposite_hard_muon_distribution( NULL ),
    far_opposite_soft_muon_distribution( NULL ),
    far_opposite_hard_pion_distribution( NULL ),
    far_opposite_soft_pion_distribution( NULL ),
    two_jets_no_leptons( given_two_jets_no_leptons ),
    one_jet_one_negative_electron( given_one_jet_one_negative_electron ),
    one_jet_one_positive_electron( given_one_jet_one_positive_electron ),
    one_jet_one_negative_muon( given_one_jet_one_negative_muon ),
    one_jet_one_positive_muon( given_one_jet_one_positive_muon ),
    no_jets_two_electrons( given_no_jets_two_electrons ),
    no_jets_two_muons( given_no_jets_two_muons ),
    no_jets_one_positive_muon_one_negative_electron(
                       given_no_jets_one_positive_muon_one_negative_electron ),
    no_jets_one_negative_muon_one_positive_electron(
                       given_no_jets_one_negative_muon_one_positive_electron ),
    one_jet_no_leptons( given_one_jet_no_leptons ),
    no_jets_one_negative_electron( given_no_jets_one_negative_electron ),
    no_jets_one_positive_electron( given_no_jets_one_positive_electron ),
    no_jets_one_negative_muon( given_no_jets_one_negative_muon ),
    no_jets_one_positive_muon( given_no_jets_one_positive_muon ),
    no_jets_no_leptons( given_no_jets_no_leptons )
  {

    // we start by setting up the BR of the (only) channel:
    bool not_particle_antiparticle_flip = true;
    if( ( scolored_is_not_antiparticle
          &&
          !EWino_is_not_antiparticle )
        ||
        ( !scolored_is_not_antiparticle
          &&
          EWino_is_not_antiparticle ) )
      {

        not_particle_antiparticle_flip = false;

      }
    first_BR
    = shortcut->get_exclusive_BR_calculator( decaying_EWino,
                                             mediating_particle,
                                             not_particle_antiparticle_flip,
                                           given_shortcuts->get_empty_list() );
    second_BR
    = shortcut->get_exclusive_BR_calculator( mediating_particle,
                                             shortcut->get_neutralino_one(),
                                             true,
                                           given_shortcuts->get_empty_list() );

    if( should_sum_charges )
      {

        near_same_distribution
        = new same_chirality_near_muon( given_shortcuts->get_readier(),
                                        given_shortcuts->get_CppSLHA(),
                                        decaying_scolored,
                     leptonKinematics->getAcceptance()->getEffectiveSquarkMass(),
                                        decaying_EWino,
                                        mediating_particle );
        far_same_distribution
        = new same_chirality_far_muon( given_shortcuts->get_readier(),
                                       given_shortcuts->get_CppSLHA(),
                                       decaying_scolored,
                     leptonKinematics->getAcceptance()->getEffectiveSquarkMass(),
                                       decaying_EWino,
                                       mediating_particle,
                                       given_shortcuts->get_neutralino_one() );
        // near_same_distribution & far_same_distribution are used even though
        // the distributions sum over charges & helicities.
        active_distributions.push_back( near_same_distribution );
        active_distributions.push_back( far_same_distribution );

        near_same_hard_muon_distribution
        = new visible_tau_decay_product( given_shortcuts->get_readier(),
                                         near_same_distribution,
                                  given_shortcuts->get_hard_muon_from_tau() );
        near_same_soft_muon_distribution
        = new visible_tau_decay_product( given_shortcuts->get_readier(),
                                 near_same_distribution,
                          given_shortcuts->get_soft_muon_from_tau() );
        near_same_hard_pion_distribution
        = new visible_tau_decay_product( given_shortcuts->get_readier(),
                                         near_same_distribution,
                                  given_shortcuts->get_hard_pion_from_tau() );
        near_same_soft_pion_distribution
        = new visible_tau_decay_product( given_shortcuts->get_readier(),
                                         near_same_distribution,
                                  given_shortcuts->get_soft_pion_from_tau() );
        far_same_hard_muon_distribution
        = new visible_tau_decay_product( given_shortcuts->get_readier(),
                                         far_same_distribution,
                                  given_shortcuts->get_hard_muon_from_tau() );
        far_same_soft_muon_distribution
        = new visible_tau_decay_product( given_shortcuts->get_readier(),
                                         far_same_distribution,
                                  given_shortcuts->get_soft_muon_from_tau() );
        far_same_hard_pion_distribution
        = new visible_tau_decay_product( given_shortcuts->get_readier(),
                                         far_same_distribution,
                                  given_shortcuts->get_hard_pion_from_tau() );
        far_same_soft_pion_distribution
        = new visible_tau_decay_product( given_shortcuts->get_readier(),
                                         far_same_distribution,
                                  given_shortcuts->get_soft_pion_from_tau() );
        active_distributions.push_back( near_same_hard_muon_distribution );
        active_distributions.push_back( near_same_soft_muon_distribution );
        active_distributions.push_back( near_same_hard_pion_distribution );
        active_distributions.push_back( near_same_soft_pion_distribution );
        active_distributions.push_back( far_same_hard_muon_distribution );
        active_distributions.push_back( far_same_soft_muon_distribution );
        active_distributions.push_back( far_same_hard_pion_distribution );
        active_distributions.push_back( far_same_soft_pion_distribution );

      }
    else
      {

        near_same_distribution
        = new same_chirality_near_muon( shortcut->get_readier(),
                                        shortcut->get_CppSLHA(),
                                        decaying_scolored,
                     leptonKinematics->getAcceptance()->getEffectiveSquarkMass(),
                                        decaying_EWino,
                                        mediating_particle );
        near_opposite_distribution
        = new opposite_chirality_near_muon( shortcut->get_readier(),
                                            shortcut->get_CppSLHA(),
                                            decaying_scolored,
                     leptonKinematics->getAcceptance()->getEffectiveSquarkMass(),
                                            decaying_EWino,
                                            mediating_particle );
        far_same_distribution
        = new same_chirality_far_muon( shortcut->get_readier(),
                                       shortcut->get_CppSLHA(),
                                       decaying_scolored,
                     leptonKinematics->getAcceptance()->getEffectiveSquarkMass(),
                                       decaying_EWino,
                                       mediating_particle,
                                       given_shortcuts->get_neutralino_one() );
        far_opposite_distribution
        = new opposite_chirality_far_muon( shortcut->get_readier(),
                                           shortcut->get_CppSLHA(),
                                           decaying_scolored,
                     leptonKinematics->getAcceptance()->getEffectiveSquarkMass(),
                                           decaying_EWino,
                                           mediating_particle,
                                       given_shortcuts->get_neutralino_one() );
        active_distributions.push_back( near_same_distribution );
        active_distributions.push_back( near_opposite_distribution );
        active_distributions.push_back( far_same_distribution );
        active_distributions.push_back( far_opposite_distribution );

        near_same_hard_muon_distribution
        = new visible_tau_decay_product( given_shortcuts->get_readier(),
                                         near_same_distribution,
                                  given_shortcuts->get_hard_muon_from_tau() );
        near_same_soft_muon_distribution
        = new visible_tau_decay_product( given_shortcuts->get_readier(),
                                 near_same_distribution,
                          given_shortcuts->get_soft_muon_from_tau() );
        near_same_hard_pion_distribution
        = new visible_tau_decay_product( given_shortcuts->get_readier(),
                                         near_same_distribution,
                                  given_shortcuts->get_hard_pion_from_tau() );
        near_same_soft_pion_distribution
        = new visible_tau_decay_product( given_shortcuts->get_readier(),
                                         near_same_distribution,
                                  given_shortcuts->get_soft_pion_from_tau() );
        near_opposite_hard_muon_distribution
        = new visible_tau_decay_product( given_shortcuts->get_readier(),
                                         near_opposite_distribution,
                                  given_shortcuts->get_hard_muon_from_tau() );
        near_opposite_soft_muon_distribution
        = new visible_tau_decay_product( given_shortcuts->get_readier(),
                                         near_opposite_distribution,
                                   given_shortcuts->get_soft_muon_from_tau() );
        near_opposite_hard_pion_distribution
        = new visible_tau_decay_product( given_shortcuts->get_readier(),
                                         near_opposite_distribution,
                                  given_shortcuts->get_hard_pion_from_tau() );
        near_opposite_soft_pion_distribution
        = new visible_tau_decay_product( given_shortcuts->get_readier(),
                                         near_opposite_distribution,
                                  given_shortcuts->get_soft_pion_from_tau() );
        far_same_hard_muon_distribution
        = new visible_tau_decay_product( given_shortcuts->get_readier(),
                                         far_same_distribution,
                                  given_shortcuts->get_hard_muon_from_tau() );
        far_same_soft_muon_distribution
        = new visible_tau_decay_product( given_shortcuts->get_readier(),
                                         far_same_distribution,
                                  given_shortcuts->get_soft_muon_from_tau() );
        far_same_hard_pion_distribution
        = new visible_tau_decay_product( given_shortcuts->get_readier(),
                                         far_same_distribution,
                                  given_shortcuts->get_hard_pion_from_tau() );
        far_same_soft_pion_distribution
        = new visible_tau_decay_product( given_shortcuts->get_readier(),
                                         far_same_distribution,
                                  given_shortcuts->get_soft_pion_from_tau() );
        far_opposite_hard_muon_distribution
        = new visible_tau_decay_product( given_shortcuts->get_readier(),
                                         far_opposite_distribution,
                                  given_shortcuts->get_hard_muon_from_tau() );
        far_opposite_soft_muon_distribution
        = new visible_tau_decay_product( given_shortcuts->get_readier(),
                                         far_opposite_distribution,
                                  given_shortcuts->get_soft_muon_from_tau() );
        far_opposite_hard_pion_distribution
        = new visible_tau_decay_product( given_shortcuts->get_readier(),
                                         far_opposite_distribution,
                                  given_shortcuts->get_hard_pion_from_tau() );
        far_opposite_soft_pion_distribution
        = new visible_tau_decay_product( given_shortcuts->get_readier(),
                                         far_opposite_distribution,
                                  given_shortcuts->get_soft_pion_from_tau() );
        active_distributions.push_back( near_same_hard_muon_distribution );
        active_distributions.push_back( near_same_soft_muon_distribution );
        active_distributions.push_back( near_same_hard_pion_distribution );
        active_distributions.push_back( near_same_soft_pion_distribution );
        active_distributions.push_back( near_opposite_hard_muon_distribution );
        active_distributions.push_back( near_opposite_soft_muon_distribution );
        active_distributions.push_back( near_opposite_hard_pion_distribution );
        active_distributions.push_back( near_opposite_soft_pion_distribution );
        active_distributions.push_back( far_same_hard_muon_distribution );
        active_distributions.push_back( far_same_soft_muon_distribution );
        active_distributions.push_back( far_same_hard_pion_distribution );
        active_distributions.push_back( far_same_soft_pion_distribution );
        active_distributions.push_back( far_opposite_hard_muon_distribution );
        active_distributions.push_back( far_opposite_soft_muon_distribution );
        active_distributions.push_back( far_opposite_hard_pion_distribution );
        active_distributions.push_back( far_opposite_soft_pion_distribution );

      }

  }

  neutralino_to_stau_cascade::~neutralino_to_stau_cascade()
  {

    for( std::vector< leptonEnergyDistribution* >::iterator
         deletion_iterator = active_distributions.begin();
         active_distributions.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

  }


  void
  neutralino_to_stau_cascade::calculate()
  {

    cascadeBr = ( first_BR->get_BR() * second_BR->get_BR() );

    if( LHC_FASER_global::negligible_BR < cascadeBr )
      // if the branching ratio into this channel is not negligible...
      {

        near_negative_tau_left_handedness
        = shortcut->quark_or_lepton_left_handedness(
                                        mediating_particle->get_PDG_code(),
                                          decaying_EWino->get_PDG_code() );
        near_negative_tau_right_handedness
        = ( 1.0 - near_negative_tau_left_handedness );
        far_negative_tau_left_handedness
        = shortcut->quark_or_lepton_left_handedness(
                                        mediating_particle->get_PDG_code(),
                              shortcut->get_neutralino_one()->get_PDG_code() );
        far_negative_tau_right_handedness
        = ( 1.0 - far_negative_tau_left_handedness );

        if( should_sum_charges_flag )
          {

            /* 1st, left-handed quark, left-handed near negative tau lepton
             * with left-handed far positive tau lepton + right-handed near
             * positive tau lepton with right-handed far negative tau lepton,
             * plus the all-handedness-flipped versions:
             */
            current_near_muon_distribution = near_same_hard_muon_distribution;
            current_near_pion_distribution = near_same_soft_pion_distribution;
            current_far_muon_distribution = far_same_soft_muon_distribution;
            current_far_pion_distribution = far_same_hard_pion_distribution;
            negative_near_configuration_BR
            = ( cascadeBr * ( near_negative_tau_left_handedness
                               * far_negative_tau_right_handedness ) );
            charge_summed_calculate();

            /* 2nd, left-handed quark, left-handed near negative tau lepton
             * with right-handed far positive tau lepton + right-handed near
             * positive tau lepton with left-handed far negative tau lepton,
             * plus the all-handedness-flipped versions:
             */
            current_near_muon_distribution = near_same_hard_muon_distribution;
            current_near_pion_distribution = near_same_soft_pion_distribution;
            current_far_muon_distribution = far_same_hard_muon_distribution;
            current_far_pion_distribution = far_same_soft_pion_distribution;
            negative_near_configuration_BR
            = ( cascadeBr * ( near_negative_tau_left_handedness
                               * far_negative_tau_left_handedness ) );
            charge_summed_calculate();

            /* 3rd, left-handed quark, right-handed near negative tau lepton
             * with left-handed far positive tau lepton + left-handed near
             * positive tau lepton with right-handed far negative tau lepton,
             * plus the all-handedness-flipped versions:
             */
            current_near_muon_distribution
            = near_opposite_soft_muon_distribution;
            current_near_pion_distribution
            = near_opposite_hard_pion_distribution;
            current_far_muon_distribution
            = far_opposite_soft_muon_distribution;
            current_far_pion_distribution
            = far_opposite_hard_pion_distribution;
            negative_near_configuration_BR
            = ( cascadeBr * ( near_negative_tau_right_handedness
                               * far_negative_tau_right_handedness ) );
            charge_summed_calculate();

            /* 4th, left-handed quark, right-handed near negative tau lepton
             * with right-handed far positive tau lepton + left-handed near
             * positive tau lepton with left-handed far negative tau lepton,
             * plus the all-handedness-flipped versions:
             */
            current_near_muon_distribution
            = near_opposite_soft_muon_distribution;
            current_near_pion_distribution
            = near_opposite_hard_pion_distribution;
            current_far_muon_distribution
            = far_opposite_hard_muon_distribution;
            current_far_pion_distribution
            = far_opposite_soft_pion_distribution;
            negative_near_configuration_BR
            = ( cascadeBr * ( near_negative_tau_right_handedness
                               * far_negative_tau_left_handedness ) );
            charge_summed_calculate();

          }
        else
          {

            jet_left_handedness
            = shortcut->quark_or_lepton_left_handedness(
                                             decaying_scolored->get_PDG_code(),
                                              decaying_EWino->get_PDG_code() );
            jet_right_handedness = ( 1.0 - jet_left_handedness );

            /* in the following, XY(bar)Z(bar) is quark handedness, near tau
             * lepton handedness, far tau lepton handedness, YbarZ means that
             * the near tau lepton is positive while YZbar means that the far
             * tau lepton is positive.
             */

            // near lepton same handedness as quark:

            // hard near muons, hard far muons:
            current_near_muon_distribution = near_same_hard_muon_distribution;
            current_near_pion_distribution = near_same_soft_pion_distribution;
            current_far_muon_distribution = far_same_hard_muon_distribution;
            current_far_pion_distribution = far_same_soft_pion_distribution;

            // LLRbar
            negative_near_configuration_BR
            = ( cascadeBr * jet_left_handedness
                * near_negative_tau_left_handedness
                * far_negative_tau_left_handedness );
            // RRbarL
            positive_near_configuration_BR
            = ( cascadeBr * jet_right_handedness
                * near_negative_tau_left_handedness
                * far_negative_tau_left_handedness );
            charge_distinguished_calculate();


            // hard near muons, soft far muons:
            current_near_muon_distribution = near_same_hard_muon_distribution;
            current_near_pion_distribution = near_same_soft_pion_distribution;
            current_far_muon_distribution = far_same_soft_muon_distribution;
            current_far_pion_distribution = far_same_hard_pion_distribution;

            // LLLbar
            negative_near_configuration_BR
            = ( cascadeBr * jet_left_handedness
                * near_negative_tau_left_handedness
                * far_negative_tau_right_handedness );
            // RRbarR
            positive_near_configuration_BR
            = ( cascadeBr * jet_right_handedness
                * near_negative_tau_left_handedness
                * far_negative_tau_right_handedness );
            charge_distinguished_calculate();


            // soft near muons, hard far muons:
            current_near_muon_distribution = near_same_soft_muon_distribution;
            current_near_pion_distribution = near_same_hard_pion_distribution;
            current_far_muon_distribution = far_same_hard_muon_distribution;
            current_far_pion_distribution = far_same_soft_pion_distribution;

            // RRRbar
            negative_near_configuration_BR
            = ( cascadeBr * jet_right_handedness
                * near_negative_tau_right_handedness
                * far_negative_tau_left_handedness );
            // LLbarL
            positive_near_configuration_BR
            = ( cascadeBr * jet_left_handedness
                * near_negative_tau_right_handedness
                * far_negative_tau_left_handedness );
            charge_distinguished_calculate();


            // soft near muons, soft far muons:
            current_near_muon_distribution = near_same_soft_muon_distribution;
            current_near_pion_distribution = near_same_hard_pion_distribution;
            current_far_muon_distribution = far_same_soft_muon_distribution;
            current_far_pion_distribution = far_same_hard_pion_distribution;

            // RRLbar
            negative_near_configuration_BR
            = ( cascadeBr * jet_right_handedness
                * near_negative_tau_right_handedness
                * far_negative_tau_right_handedness );
            // LLbarR
            positive_near_configuration_BR
            = ( cascadeBr * jet_left_handedness
                * near_negative_tau_right_handedness
                * far_negative_tau_right_handedness );
            charge_distinguished_calculate();


            // near lepton opposite handedness to quark:

            // hard near muons, hard far muons:
            current_near_muon_distribution
            = near_opposite_hard_muon_distribution;
            current_near_pion_distribution
            = near_opposite_soft_pion_distribution;
            current_far_muon_distribution
            = far_opposite_hard_muon_distribution;
            current_far_pion_distribution
            = far_opposite_soft_pion_distribution;

            // RLRbar
            negative_near_configuration_BR
            = ( cascadeBr * jet_right_handedness
                * near_negative_tau_left_handedness
                * far_negative_tau_left_handedness );
            // LRbarL
            positive_near_configuration_BR
            = ( cascadeBr * jet_left_handedness
                * near_negative_tau_left_handedness
                * far_negative_tau_left_handedness );
            charge_distinguished_calculate();


            // hard near muons, soft far muons:
            current_near_muon_distribution
            = near_opposite_hard_muon_distribution;
            current_near_pion_distribution
            = near_opposite_soft_pion_distribution;
            current_far_muon_distribution
            = far_opposite_soft_muon_distribution;
            current_far_pion_distribution
            = far_opposite_hard_pion_distribution;

            // RLLbar
            negative_near_configuration_BR
            = ( cascadeBr * jet_right_handedness
                * near_negative_tau_left_handedness
                * far_negative_tau_right_handedness );
            // LRbarR
            positive_near_configuration_BR
            = ( cascadeBr * jet_left_handedness
                * near_negative_tau_left_handedness
                * far_negative_tau_right_handedness );
            charge_distinguished_calculate();


            // soft near muons, hard far muons:
            current_near_muon_distribution
            = near_opposite_soft_muon_distribution;
            current_near_pion_distribution
            = near_opposite_hard_pion_distribution;
            current_far_muon_distribution
            = far_opposite_hard_muon_distribution;
            current_far_pion_distribution
            = far_opposite_soft_pion_distribution;

            // LRRbar
            negative_near_configuration_BR
            = ( cascadeBr * jet_left_handedness
                * near_negative_tau_right_handedness
                * far_negative_tau_left_handedness );
            // RLbarL
            positive_near_configuration_BR
            = ( cascadeBr * jet_right_handedness
                * near_negative_tau_right_handedness
                * far_negative_tau_left_handedness );
            charge_distinguished_calculate();


            // soft near muons, soft far muons:
            current_near_muon_distribution
            = near_opposite_soft_muon_distribution;
            current_near_pion_distribution
            = near_opposite_hard_pion_distribution;
            current_far_muon_distribution
            = far_opposite_soft_muon_distribution;
            current_far_pion_distribution
            = far_opposite_hard_pion_distribution;

            // LRLbar
            negative_near_configuration_BR
            = ( cascadeBr * jet_left_handedness
                * near_negative_tau_right_handedness
                * far_negative_tau_right_handedness );
            // RLbarR
            positive_near_configuration_BR
            = ( cascadeBr * jet_right_handedness
                * near_negative_tau_right_handedness
                * far_negative_tau_right_handedness );
            charge_distinguished_calculate();


            // debugging:
            /**std::cout
            << std::endl
            << "cascade_BR = " << cascade_BR;
            std::cout << std::endl;**/

          }  // end of if summing over charges or not.

      }  // end of if cascade branching ratio was large enough to care about.

  }

  void
  neutralino_to_stau_cascade::charge_summed_calculate()
  /* this updates values based on the current near & far muon/electron & pion
   * distributions & configuration branching ratio, summing over tau lepton
   * charges.
   */
  {

    // at this point, ( 2 * negative_near_configuration_BR ) is the sum of the
    // BRs into this configuration.

    near_muon_pass = integrateAcceptance( current_near_muon_distribution, cuts->getPrimaryLeptonCut() );
    near_muon_fail = ( 1.0 - integrateAcceptance( current_near_muon_distribution, cuts->getSecondaryLeptonCut() ) );
    near_pion_pass = integrateAcceptance( current_near_pion_distribution, jet_cut );
    near_pion_fail = ( 1.0 - near_pion_pass );
    far_muon_pass = integrateAcceptance( current_far_muon_distribution, cuts->getPrimaryLeptonCut() );
    far_muon_fail = ( 1.0 - integrateAcceptance( current_far_muon_distribution, cuts->getSecondaryLeptonCut() ) );
    far_pion_pass = integrateAcceptance( current_far_pion_distribution, jet_cut );
    far_pion_fail = ( 1.0 - far_pion_pass );

    *two_jets_no_leptons
    += ( 2.0 * negative_near_configuration_BR
         * tau_pair_to_jet_pair_BR
         * near_pion_pass * far_pion_pass );
    // the factor of 2.0 accounts for the charge-conjugate.

    both_pass
    = ( 2.0 * ( near_pion_pass * far_muon_pass
                + near_muon_pass * far_pion_pass ) );
    // the factor of 2.0 accounts for the charge-conjugate.
    current_pass = ( negative_near_configuration_BR
                     * tau_pair_to_jet_electron_BR * both_pass );
    *one_jet_one_negative_electron += current_pass;
    *one_jet_one_positive_electron += current_pass;
    current_pass = ( negative_near_configuration_BR
                     * tau_pair_to_jet_muon_BR * both_pass );
    *one_jet_one_negative_muon += current_pass;
    *one_jet_one_positive_muon += current_pass;
    both_pass = ( 2.0 * near_muon_pass * far_muon_pass );
    // the factor of 2.0 accounts for the charge-conjugate.
    *no_jets_two_electrons += ( negative_near_configuration_BR
                                * tau_pair_to_electron_pair_BR * both_pass );
    *no_jets_two_muons += ( negative_near_configuration_BR
                            * tau_pair_to_muon_pair_BR * both_pass );
    current_pass = ( negative_near_configuration_BR
                     * tau_pair_to_muon_electron_BR * both_pass );
    *no_jets_one_positive_muon_one_negative_electron += current_pass;
    *no_jets_one_negative_muon_one_positive_electron += current_pass;

    *one_jet_no_leptons
    += ( 2.0 * negative_near_configuration_BR
         * ( tau_pair_to_jet_pair_BR
             * ( near_pion_pass * far_pion_fail
                 + near_pion_fail * far_pion_pass ) )
             + ( tau_pair_to_jet_muon_BR + tau_pair_to_jet_electron_BR )
               * ( near_pion_pass * far_muon_fail
                   + near_muon_fail * far_pion_pass ) );
    // the factor of 2.0 accounts for the charge-conjugate.
    current_pass
    = ( negative_near_configuration_BR
        * ( tau_pair_to_jet_electron_BR
            * ( near_muon_pass * far_pion_fail
                + near_pion_fail * far_muon_pass )
            + ( tau_pair_to_electron_pair_BR + tau_pair_to_muon_electron_BR )
              * ( near_muon_pass * far_muon_fail
                  + near_muon_fail * far_muon_pass ) ) );
    *no_jets_one_negative_electron += current_pass;
    *no_jets_one_positive_electron += current_pass;
    current_pass
    = ( negative_near_configuration_BR
        * ( tau_pair_to_jet_muon_BR
            * ( near_muon_pass * far_pion_fail
                + near_pion_fail * far_muon_pass )
            + ( tau_pair_to_muon_pair_BR + tau_pair_to_muon_electron_BR )
              * ( near_muon_pass * far_muon_fail
                  + near_muon_fail * far_muon_pass ) ) );
    *no_jets_one_negative_muon += current_pass;
    *no_jets_one_positive_muon += current_pass;

    *no_jets_no_leptons
    += ( 2.0 * negative_near_configuration_BR
         * ( tau_pair_to_jet_pair_BR * near_pion_fail * far_pion_fail
             + ( tau_pair_to_jet_muon_BR + tau_pair_to_jet_electron_BR )
               * ( near_pion_fail * far_muon_fail
                   + near_muon_fail * far_pion_fail )
             + ( tau_pair_to_muon_pair_BR + 2.0 * tau_pair_to_muon_electron_BR
                 + tau_pair_to_electron_pair_BR )
               * near_muon_fail * far_muon_fail ) );
    // the factor of 2.0 accounts for the charge-conjugate.

  }

  void
  neutralino_to_stau_cascade::charge_distinguished_calculate()
  /* this updates values based on the current near & far muon/electron & pion
   * distributions & configuration branching ratio, separately recording
   * values for differing charges.
   */
  {

    // at this point, negative_near_configuration_BR is the BR with a negative
    // near tau lepton.

    near_muon_pass = integrateAcceptance( current_near_muon_distribution, cuts->getPrimaryLeptonCut() );
    near_muon_fail = ( 1.0 - integrateAcceptance( current_near_muon_distribution, cuts->getSecondaryLeptonCut() ) );
    near_pion_pass = integrateAcceptance( current_near_pion_distribution, jet_cut );
    near_pion_fail = ( 1.0 - near_pion_pass );
    far_muon_pass = integrateAcceptance( current_far_muon_distribution, cuts->getPrimaryLeptonCut() );
    far_muon_fail = ( 1.0 - integrateAcceptance( current_far_muon_distribution, cuts->getSecondaryLeptonCut() ) );
    far_pion_pass = integrateAcceptance( current_far_pion_distribution, jet_cut );
    far_pion_fail = ( 1.0 - far_pion_pass );

    *two_jets_no_leptons
    += ( ( negative_near_configuration_BR + positive_near_configuration_BR )
         * tau_pair_to_jet_pair_BR
         * near_pion_pass * far_pion_pass );

    both_pass
    = ( negative_near_configuration_BR * near_muon_pass * far_pion_pass
        + positive_near_configuration_BR * near_pion_pass * far_muon_pass );
    *one_jet_one_negative_electron
    += ( tau_pair_to_jet_electron_BR * both_pass );
    *one_jet_one_negative_muon
    += ( tau_pair_to_jet_muon_BR * both_pass );
    both_pass
    = ( positive_near_configuration_BR * near_muon_pass * far_pion_pass
        + negative_near_configuration_BR * near_pion_pass * far_muon_pass );
    *one_jet_one_positive_electron
    += ( tau_pair_to_jet_electron_BR * both_pass );
    *one_jet_one_positive_muon
    += ( tau_pair_to_jet_muon_BR * both_pass );
    both_pass = ( near_muon_pass * far_muon_pass );
    current_pass
    = ( ( negative_near_configuration_BR + positive_near_configuration_BR )
        * both_pass );
    *no_jets_two_muons += ( tau_pair_to_muon_pair_BR * both_pass );
    *no_jets_two_electrons += ( tau_pair_to_electron_pair_BR * both_pass );
    current_pass *= tau_pair_to_muon_electron_BR;
    *no_jets_one_positive_muon_one_negative_electron += current_pass;
    *no_jets_one_negative_muon_one_positive_electron += current_pass;

    *one_jet_no_leptons
    += ( ( negative_near_configuration_BR + positive_near_configuration_BR )
         * ( tau_pair_to_jet_pair_BR
             * ( near_pion_pass * far_pion_fail
                 + near_pion_fail * far_pion_pass ) )
             + ( tau_pair_to_jet_muon_BR + tau_pair_to_jet_electron_BR )
               * ( near_pion_pass * far_muon_fail
                   + near_muon_fail * far_pion_pass ) );
    *no_jets_one_negative_muon
    += ( ( negative_near_configuration_BR * near_muon_pass )
         * ( tau_pair_to_jet_muon_BR * far_pion_fail
             + tau_pair_to_muon_pair_BR * far_muon_fail )
         + ( positive_near_configuration_BR * far_muon_pass )
           * ( tau_pair_to_jet_muon_BR * near_pion_fail
               + tau_pair_to_muon_pair_BR * near_muon_fail ) );
    *no_jets_one_negative_electron
    += ( ( negative_near_configuration_BR * near_muon_pass )
         * ( tau_pair_to_jet_electron_BR * far_pion_fail
             + tau_pair_to_electron_pair_BR * far_muon_fail )
         + ( positive_near_configuration_BR * far_muon_pass )
           * ( tau_pair_to_jet_electron_BR * near_pion_fail
               + tau_pair_to_electron_pair_BR * near_muon_fail ) );
    *no_jets_one_positive_muon
    += ( ( positive_near_configuration_BR * near_muon_pass )
         * ( tau_pair_to_jet_muon_BR * far_pion_fail
             + tau_pair_to_muon_pair_BR * far_muon_fail )
         + ( negative_near_configuration_BR * far_muon_pass )
           * ( tau_pair_to_jet_muon_BR * near_pion_fail
               + tau_pair_to_muon_pair_BR * near_muon_fail ) );
    *no_jets_one_positive_electron
    += ( ( positive_near_configuration_BR * near_muon_pass )
         * ( tau_pair_to_jet_electron_BR * far_pion_fail
             + tau_pair_to_electron_pair_BR * far_muon_fail )
         + ( negative_near_configuration_BR * far_muon_pass )
           * ( tau_pair_to_jet_electron_BR * near_pion_fail
               + tau_pair_to_electron_pair_BR * near_muon_fail ) );

    *no_jets_no_leptons
    += ( ( negative_near_configuration_BR + positive_near_configuration_BR )
         * ( tau_pair_to_jet_pair_BR * near_pion_fail * far_pion_fail
             + ( tau_pair_to_jet_muon_BR + tau_pair_to_jet_electron_BR )
               * ( near_pion_fail * far_muon_fail
                   + near_muon_fail * far_pion_fail )
             + ( tau_pair_to_muon_pair_BR + 2.0 * tau_pair_to_muon_electron_BR
                 + tau_pair_to_electron_pair_BR )
               * near_muon_fail * far_muon_fail ) );

  }


  neutralino_to_Z_cascade::neutralino_to_Z_cascade(
                                    readier_for_new_point* const given_readier,
                                                 bool const should_sum_charges,
                                                double const given_cuts->getPrimaryLeptonCut(),
                                              double const given_cuts->getSecondaryLeptonCut(),
                               lepton_acceptance_value* const given_kinematics,
           CppSLHA::particle_property_set const* const given_decaying_scolored,
                                 bool const given_scolored_is_not_antiparticle,
              CppSLHA::particle_property_set const* const given_decaying_EWino,
                                    bool const given_EWino_is_not_antiparticle,
          CppSLHA::particle_property_set const* const given_mediating_particle,
                                    input_handler const* const given_shortcuts,
                                       double* const given_two_jets_no_leptons,
                             double* const given_one_jet_one_negative_electron,
                             double* const given_one_jet_one_positive_electron,
                                 double* const given_one_jet_one_negative_muon,
                                 double* const given_one_jet_one_positive_muon,
                                   double* const given_OSSF_minus_OSDF_leptons,
                                     double* const given_no_jets_two_electrons,
                                         double* const given_no_jets_two_muons,
           double* const given_no_jets_one_positive_muon_one_negative_electron,
           double* const given_no_jets_one_negative_muon_one_positive_electron,
                                        double* const given_one_jet_no_leptons,
                             double* const given_no_jets_one_negative_electron,
                             double* const given_no_jets_one_positive_electron,
                                 double* const given_no_jets_one_negative_muon,
                                 double* const given_no_jets_one_positive_muon,
                                     double* const given_no_jets_no_leptons ) :
    channel_calculator( given_readier,
                        given_cuts->getPrimaryLeptonCut(),
                        given_cuts->getSecondaryLeptonCut(),
                        -1.0,
                        given_kinematics,
                        given_decaying_scolored,
                        given_scolored_is_not_antiparticle,
                        given_decaying_EWino,
                        given_EWino_is_not_antiparticle,
                        given_mediating_particle,
                        given_shortcuts ),
    should_sum_charges_flag( should_sum_charges ),
    direct_distribution( NULL ),
    hard_muon_distribution( NULL ),
    soft_muon_distribution( NULL ),
    hard_pion_distribution( NULL ),
    soft_pion_distribution( NULL ),
    two_jets_no_leptons( given_two_jets_no_leptons ),
    one_jet_one_negative_electron( given_one_jet_one_negative_electron ),
    one_jet_one_positive_electron( given_one_jet_one_positive_electron ),
    one_jet_one_negative_muon( given_one_jet_one_negative_muon ),
    one_jet_one_positive_muon( given_one_jet_one_positive_muon ),
    OSSF_minus_OSDF_leptons( given_OSSF_minus_OSDF_leptons ),
    no_jets_two_electrons( given_no_jets_two_electrons ),
    no_jets_two_muons( given_no_jets_two_muons ),
    no_jets_one_positive_muon_one_negative_electron(
                       given_no_jets_one_positive_muon_one_negative_electron ),
    no_jets_one_negative_muon_one_positive_electron(
                       given_no_jets_one_negative_muon_one_positive_electron ),
    one_jet_no_leptons( given_one_jet_no_leptons ),
    no_jets_one_negative_electron( given_no_jets_one_negative_electron ),
    no_jets_one_positive_electron( given_no_jets_one_positive_electron ),
    no_jets_one_negative_muon( given_no_jets_one_negative_muon ),
    no_jets_one_positive_muon( given_no_jets_one_positive_muon ),
    no_jets_no_leptons( given_no_jets_no_leptons )
  {

    double
    twice_square_of_weak_cosine = ( 2.0 * given_shortcuts->get_weak_cosine()
                                        * given_shortcuts->get_weak_cosine() );
    negative_tau_left_handedness
    = ( ( 1 - 2.0 * twice_square_of_weak_cosine
          + twice_square_of_weak_cosine * twice_square_of_weak_cosine )
        / ( 5.0 - 2.0 * twice_square_of_weak_cosine
            + 2.0 * twice_square_of_weak_cosine
                  * twice_square_of_weak_cosine ) );

    // setting up the BR of the (only) channel:
  }

  neutralino_to_Z_cascade::~neutralino_to_Z_cascade()
  {

    for( std::vector< leptonEnergyDistribution* >::iterator
         deletion_iterator = active_distributions.begin();
         active_distributions.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

  }


  void
  neutralino_to_Z_cascade::calculate()
  {

    cascadeBr = ( first_BR->get_BR() * second_BR->get_BR() );

    if( LHC_FASER_global::negligible_BR < cascadeBr )
      // if the branching ratio into this channel is not negligible...
      {

        near_negative_tau_left_handedness
        = shortcut->quark_or_lepton_left_handedness(
                                        mediating_particle->get_PDG_code(),
                                          decaying_EWino->get_PDG_code() );
        near_negative_tau_right_handedness
        = ( 1.0 - near_negative_tau_left_handedness );
        far_negative_tau_left_handedness
        = shortcut->quark_or_lepton_left_handedness(
                                        mediating_particle->get_PDG_code(),
                              shortcut->get_neutralino_one()->get_PDG_code() );
        far_negative_tau_right_handedness
        = ( 1.0 - far_negative_tau_left_handedness );

        if( should_sum_charges_flag )
          {

            /* 1st, left-handed quark, left-handed near negative tau lepton
             * with left-handed far positive tau lepton + right-handed near
             * positive tau lepton with right-handed far negative tau lepton,
             * plus the all-handedness-flipped versions:
             */
            current_near_muon_distribution = near_same_hard_muon_distribution;
            current_near_pion_distribution = near_same_soft_pion_distribution;
            current_far_muon_distribution = far_same_soft_muon_distribution;
            current_far_pion_distribution = far_same_hard_pion_distribution;
            negative_near_configuration_BR
            = ( cascadeBr * ( near_negative_tau_left_handedness
                               * far_negative_tau_right_handedness ) );
            charge_summed_calculate();

            /* 2nd, left-handed quark, left-handed near negative tau lepton
             * with right-handed far positive tau lepton + right-handed near
             * positive tau lepton with left-handed far negative tau lepton,
             * plus the all-handedness-flipped versions:
             */
            current_near_muon_distribution = near_same_hard_muon_distribution;
            current_near_pion_distribution = near_same_soft_pion_distribution;
            current_far_muon_distribution = far_same_hard_muon_distribution;
            current_far_pion_distribution = far_same_soft_pion_distribution;
            negative_near_configuration_BR
            = ( cascadeBr * ( near_negative_tau_left_handedness
                               * far_negative_tau_left_handedness ) );
            charge_summed_calculate();

            /* 3rd, left-handed quark, right-handed near negative tau lepton
             * with left-handed far positive tau lepton + left-handed near
             * positive tau lepton with right-handed far negative tau lepton,
             * plus the all-handedness-flipped versions:
             */
            current_near_muon_distribution
            = near_opposite_soft_muon_distribution;
            current_near_pion_distribution
            = near_opposite_hard_pion_distribution;
            current_far_muon_distribution
            = far_opposite_soft_muon_distribution;
            current_far_pion_distribution
            = far_opposite_hard_pion_distribution;
            negative_near_configuration_BR
            = ( cascadeBr * ( near_negative_tau_right_handedness
                               * far_negative_tau_right_handedness ) );
            charge_summed_calculate();

            /* 4th, left-handed quark, right-handed near negative tau lepton
             * with right-handed far positive tau lepton + left-handed near
             * positive tau lepton with left-handed far negative tau lepton,
             * plus the all-handedness-flipped versions:
             */
            current_near_muon_distribution
            = near_opposite_soft_muon_distribution;
            current_near_pion_distribution
            = near_opposite_hard_pion_distribution;
            current_far_muon_distribution
            = far_opposite_hard_muon_distribution;
            current_far_pion_distribution
            = far_opposite_soft_pion_distribution;
            negative_near_configuration_BR
            = ( cascadeBr * ( near_negative_tau_right_handedness
                               * far_negative_tau_left_handedness ) );
            charge_summed_calculate();

          }
        else
          {

            jetLeftHandedness
            = shortcut->quark_or_lepton_left_handedness(
                                             decaying_scolored->get_PDG_code(),
                                              decaying_EWino->get_PDG_code() );
            jet_right_handedness = ( 1.0 - jetLeftHandedness );

            /* in the following, XY(bar)Z(bar) is quark handedness, near tau
             * lepton handedness, far tau lepton handedness, YbarZ means that
             * the near tau lepton is positive while YZbar means that the far
             * tau lepton is positive.
             */

            // near lepton same handedness as quark:

            // hard near muons, hard far muons:
            current_near_muon_distribution = near_same_hard_muon_distribution;
            current_near_pion_distribution = near_same_soft_pion_distribution;
            current_far_muon_distribution = far_same_hard_muon_distribution;
            current_far_pion_distribution = far_same_soft_pion_distribution;

            // LLRbar
            negative_near_configuration_BR
            = ( cascadeBr * jetLeftHandedness
                * near_negative_tau_left_handedness
                * far_negative_tau_left_handedness );
            // RRbarL
            positive_near_configuration_BR
            = ( cascadeBr * jet_right_handedness
                * near_negative_tau_left_handedness
                * far_negative_tau_left_handedness );
            charge_distinguished_calculate();


            // hard near muons, soft far muons:
            current_near_muon_distribution = near_same_hard_muon_distribution;
            current_near_pion_distribution = near_same_soft_pion_distribution;
            current_far_muon_distribution = far_same_soft_muon_distribution;
            current_far_pion_distribution = far_same_hard_pion_distribution;

            // LLLbar
            negative_near_configuration_BR
            = ( cascadeBr * jetLeftHandedness
                * near_negative_tau_left_handedness
                * far_negative_tau_right_handedness );
            // RRbarR
            positive_near_configuration_BR
            = ( cascadeBr * jet_right_handedness
                * near_negative_tau_left_handedness
                * far_negative_tau_right_handedness );
            charge_distinguished_calculate();


            // soft near muons, hard far muons:
            current_near_muon_distribution = near_same_soft_muon_distribution;
            current_near_pion_distribution = near_same_hard_pion_distribution;
            current_far_muon_distribution = far_same_hard_muon_distribution;
            current_far_pion_distribution = far_same_soft_pion_distribution;

            // RRRbar
            negative_near_configuration_BR
            = ( cascadeBr * jet_right_handedness
                * near_negative_tau_right_handedness
                * far_negative_tau_left_handedness );
            // LLbarL
            positive_near_configuration_BR
            = ( cascadeBr * jetLeftHandedness
                * near_negative_tau_right_handedness
                * far_negative_tau_left_handedness );
            charge_distinguished_calculate();


            // soft near muons, soft far muons:
            current_near_muon_distribution = near_same_soft_muon_distribution;
            current_near_pion_distribution = near_same_hard_pion_distribution;
            current_far_muon_distribution = far_same_soft_muon_distribution;
            current_far_pion_distribution = far_same_hard_pion_distribution;

            // RRLbar
            negative_near_configuration_BR
            = ( cascadeBr * jet_right_handedness
                * near_negative_tau_right_handedness
                * far_negative_tau_right_handedness );
            // LLbarR
            positive_near_configuration_BR
            = ( cascadeBr * jetLeftHandedness
                * near_negative_tau_right_handedness
                * far_negative_tau_right_handedness );
            charge_distinguished_calculate();


            // near lepton opposite handedness to quark:

            // hard near muons, hard far muons:
            current_near_muon_distribution
            = near_opposite_hard_muon_distribution;
            current_near_pion_distribution
            = near_opposite_soft_pion_distribution;
            current_far_muon_distribution
            = far_opposite_hard_muon_distribution;
            current_far_pion_distribution
            = far_opposite_soft_pion_distribution;

            // RLRbar
            negative_near_configuration_BR
            = ( cascadeBr * jet_right_handedness
                * near_negative_tau_left_handedness
                * far_negative_tau_left_handedness );
            // LRbarL
            positive_near_configuration_BR
            = ( cascadeBr * jetLeftHandedness
                * near_negative_tau_left_handedness
                * far_negative_tau_left_handedness );
            charge_distinguished_calculate();


            // hard near muons, soft far muons:
            current_near_muon_distribution
            = near_opposite_hard_muon_distribution;
            current_near_pion_distribution
            = near_opposite_soft_pion_distribution;
            current_far_muon_distribution
            = far_opposite_soft_muon_distribution;
            current_far_pion_distribution
            = far_opposite_hard_pion_distribution;

            // RLLbar
            negative_near_configuration_BR
            = ( cascadeBr * jet_right_handedness
                * near_negative_tau_left_handedness
                * far_negative_tau_right_handedness );
            // LRbarR
            positive_near_configuration_BR
            = ( cascadeBr * jetLeftHandedness
                * near_negative_tau_left_handedness
                * far_negative_tau_right_handedness );
            charge_distinguished_calculate();


            // soft near muons, hard far muons:
            current_near_muon_distribution
            = near_opposite_soft_muon_distribution;
            current_near_pion_distribution
            = near_opposite_hard_pion_distribution;
            current_far_muon_distribution
            = far_opposite_hard_muon_distribution;
            current_far_pion_distribution
            = far_opposite_soft_pion_distribution;

            // LRRbar
            negative_near_configuration_BR
            = ( cascadeBr * jetLeftHandedness
                * near_negative_tau_right_handedness
                * far_negative_tau_left_handedness );
            // RLbarL
            positive_near_configuration_BR
            = ( cascadeBr * jet_right_handedness
                * near_negative_tau_right_handedness
                * far_negative_tau_left_handedness );
            charge_distinguished_calculate();


            // soft near muons, soft far muons:
            current_near_muon_distribution
            = near_opposite_soft_muon_distribution;
            current_near_pion_distribution
            = near_opposite_hard_pion_distribution;
            current_far_muon_distribution
            = far_opposite_soft_muon_distribution;
            current_far_pion_distribution
            = far_opposite_hard_pion_distribution;

            // LRLbar
            negative_near_configuration_BR
            = ( cascadeBr * jetLeftHandedness
                * near_negative_tau_right_handedness
                * far_negative_tau_right_handedness );
            // RLbarR
            positive_near_configuration_BR
            = ( cascadeBr * jet_right_handedness
                * near_negative_tau_right_handedness
                * far_negative_tau_right_handedness );
            charge_distinguished_calculate();


            // debugging:
            /**std::cout
            << std::endl
            << "cascade_BR = " << cascade_BR;
            std::cout << std::endl;**/

          }  // end of if summing over charges or not.

      }  // end of if cascade branching ratio was large enough to care about.

  }

  void
  neutralino_to_Z_cascade::charge_summed_calculate()
  /* this updates values based on the current near & far muon/electron & pion
   * distributions & configuration branching ratio, summing over tau lepton
   * charges.
   */
  {

    // at this point, ( 2 * negative_near_configuration_BR ) is the sum of the
    // BRs into this configuration.

    near_muon_pass = integrateAcceptance( current_near_muon_distribution, cuts->getPrimaryLeptonCut() );
    near_muon_fail = ( 1.0 - integrateAcceptance( current_near_muon_distribution, cuts->getSecondaryLeptonCut() ) );
    near_pion_pass = integrateAcceptance( current_near_pion_distribution, jet_cut );
    near_pion_fail = ( 1.0 - near_pion_pass );
    far_muon_pass = integrateAcceptance( current_far_muon_distribution, cuts->getPrimaryLeptonCut() );
    far_muon_fail = ( 1.0 - integrateAcceptance( current_far_muon_distribution, cuts->getSecondaryLeptonCut() ) );
    far_pion_pass = integrateAcceptance( current_far_pion_distribution, jet_cut );
    far_pion_fail = ( 1.0 - far_pion_pass );

    *two_jets_no_leptons
    += ( 2.0 * negative_near_configuration_BR
         * tau_pair_to_jet_pair_BR
         * near_pion_pass * far_pion_pass );
    // the factor of 2.0 accounts for the charge-conjugate.

    both_pass
    = ( 2.0 * ( near_pion_pass * far_muon_pass
                + near_muon_pass * far_pion_pass ) );
    // the factor of 2.0 accounts for the charge-conjugate.
    current_pass = ( negative_near_configuration_BR
                     * tau_pair_to_jet_electron_BR * both_pass );
    *one_jet_one_negative_electron += current_pass;
    *one_jet_one_positive_electron += current_pass;
    current_pass = ( negative_near_configuration_BR
                     * tau_pair_to_jet_muon_BR * both_pass );
    *one_jet_one_negative_muon += current_pass;
    *one_jet_one_positive_muon += current_pass;
    both_pass = ( 2.0 * near_muon_pass * far_muon_pass );
    // the factor of 2.0 accounts for the charge-conjugate.
    *no_jets_two_electrons += ( negative_near_configuration_BR
                                * tau_pair_to_electron_pair_BR * both_pass );
    *no_jets_two_muons += ( negative_near_configuration_BR
                            * tau_pair_to_muon_pair_BR * both_pass );
    current_pass = ( negative_near_configuration_BR
                     * tau_pair_to_muon_electron_BR * both_pass );
    *no_jets_one_positive_muon_one_negative_electron += current_pass;
    *no_jets_one_negative_muon_one_positive_electron += current_pass;

    *one_jet_no_leptons
    += ( 2.0 * negative_near_configuration_BR
         * ( tau_pair_to_jet_pair_BR
             * ( near_pion_pass * far_pion_fail
                 + near_pion_fail * far_pion_pass ) )
             + ( tau_pair_to_jet_muon_BR + tau_pair_to_jet_electron_BR )
               * ( near_pion_pass * far_muon_fail
                   + near_muon_fail * far_pion_pass ) );
    // the factor of 2.0 accounts for the charge-conjugate.
    current_pass
    = ( negative_near_configuration_BR
        * ( tau_pair_to_jet_electron_BR
            * ( near_muon_pass * far_pion_fail
                + near_pion_fail * far_muon_pass )
            + ( tau_pair_to_electron_pair_BR + tau_pair_to_muon_electron_BR )
              * ( near_muon_pass * far_muon_fail
                  + near_muon_fail * far_muon_pass ) ) );
    *no_jets_one_negative_electron += current_pass;
    *no_jets_one_positive_electron += current_pass;
    current_pass
    = ( negative_near_configuration_BR
        * ( tau_pair_to_jet_muon_BR
            * ( near_muon_pass * far_pion_fail
                + near_pion_fail * far_muon_pass )
            + ( tau_pair_to_muon_pair_BR + tau_pair_to_muon_electron_BR )
              * ( near_muon_pass * far_muon_fail
                  + near_muon_fail * far_muon_pass ) ) );
    *no_jets_one_negative_muon += current_pass;
    *no_jets_one_positive_muon += current_pass;

    *no_jets_no_leptons
    += ( 2.0 * negative_near_configuration_BR
         * ( tau_pair_to_jet_pair_BR * near_pion_fail * far_pion_fail
             + ( tau_pair_to_jet_muon_BR + tau_pair_to_jet_electron_BR )
               * ( near_pion_fail * far_muon_fail
                   + near_muon_fail * far_pion_fail )
             + ( tau_pair_to_muon_pair_BR + 2.0 * tau_pair_to_muon_electron_BR
                 + tau_pair_to_electron_pair_BR )
               * near_muon_fail * far_muon_fail ) );
    // the factor of 2.0 accounts for the charge-conjugate.

  }

  void
  neutralino_to_Z_cascade::charge_distinguished_calculate()
  /* this updates values based on the current near & far muon/electron & pion
   * distributions & configuration branching ratio, separately recording
   * values for differing charges.
   */
  {

    // at this point, negative_near_configuration_BR is the BR with a negative
    // near tau lepton.

    near_muon_pass = integrateAcceptance( current_near_muon_distribution, cuts->getPrimaryLeptonCut() );
    near_muon_fail = ( 1.0 - integrateAcceptance( current_near_muon_distribution, cuts->getSecondaryLeptonCut() ) );
    near_pion_pass = integrateAcceptance( current_near_pion_distribution, jet_cut );
    near_pion_fail = ( 1.0 - near_pion_pass );
    far_muon_pass = integrateAcceptance( current_far_muon_distribution, cuts->getPrimaryLeptonCut() );
    far_muon_fail = ( 1.0 - integrateAcceptance( current_far_muon_distribution, cuts->getSecondaryLeptonCut() ) );
    far_pion_pass = integrateAcceptance( current_far_pion_distribution, jet_cut );
    far_pion_fail = ( 1.0 - far_pion_pass );

    *two_jets_no_leptons
    += ( ( negative_near_configuration_BR + positive_near_configuration_BR )
         * tau_pair_to_jet_pair_BR
         * near_pion_pass * far_pion_pass );

    both_pass
    = ( negative_near_configuration_BR * near_muon_pass * far_pion_pass
        + positive_near_configuration_BR * near_pion_pass * far_muon_pass );
    *one_jet_one_negative_electron
    += ( tau_pair_to_jet_electron_BR * both_pass );
    *one_jet_one_negative_muon
    += ( tau_pair_to_jet_muon_BR * both_pass );
    both_pass
    = ( positive_near_configuration_BR * near_muon_pass * far_pion_pass
        + negative_near_configuration_BR * near_pion_pass * far_muon_pass );
    *one_jet_one_positive_electron
    += ( tau_pair_to_jet_electron_BR * both_pass );
    *one_jet_one_positive_muon
    += ( tau_pair_to_jet_muon_BR * both_pass );
    both_pass = ( near_muon_pass * far_muon_pass );
    current_pass
    = ( ( negative_near_configuration_BR + positive_near_configuration_BR )
        * both_pass );
    *no_jets_two_muons += ( tau_pair_to_muon_pair_BR * both_pass );
    *no_jets_two_electrons += ( tau_pair_to_electron_pair_BR * both_pass );
    current_pass *= tau_pair_to_muon_electron_BR;
    *no_jets_one_positive_muon_one_negative_electron += current_pass;
    *no_jets_one_negative_muon_one_positive_electron += current_pass;

    *one_jet_no_leptons
    += ( ( negative_near_configuration_BR + positive_near_configuration_BR )
         * ( tau_pair_to_jet_pair_BR
             * ( near_pion_pass * far_pion_fail
                 + near_pion_fail * far_pion_pass ) )
             + ( tau_pair_to_jet_muon_BR + tau_pair_to_jet_electron_BR )
               * ( near_pion_pass * far_muon_fail
                   + near_muon_fail * far_pion_pass ) );
    *no_jets_one_negative_muon
    += ( ( negative_near_configuration_BR * near_muon_pass )
         * ( tau_pair_to_jet_muon_BR * far_pion_fail
             + tau_pair_to_muon_pair_BR * far_muon_fail )
         + ( positive_near_configuration_BR * far_muon_pass )
           * ( tau_pair_to_jet_muon_BR * near_pion_fail
               + tau_pair_to_muon_pair_BR * near_muon_fail ) );
    *no_jets_one_negative_electron
    += ( ( negative_near_configuration_BR * near_muon_pass )
         * ( tau_pair_to_jet_electron_BR * far_pion_fail
             + tau_pair_to_electron_pair_BR * far_muon_fail )
         + ( positive_near_configuration_BR * far_muon_pass )
           * ( tau_pair_to_jet_electron_BR * near_pion_fail
               + tau_pair_to_electron_pair_BR * near_muon_fail ) );
    *no_jets_one_positive_muon
    += ( ( positive_near_configuration_BR * near_muon_pass )
         * ( tau_pair_to_jet_muon_BR * far_pion_fail
             + tau_pair_to_muon_pair_BR * far_muon_fail )
         + ( negative_near_configuration_BR * far_muon_pass )
           * ( tau_pair_to_jet_muon_BR * near_pion_fail
               + tau_pair_to_muon_pair_BR * near_muon_fail ) );
    *no_jets_one_positive_electron
    += ( ( positive_near_configuration_BR * near_muon_pass )
         * ( tau_pair_to_jet_electron_BR * far_pion_fail
             + tau_pair_to_electron_pair_BR * far_muon_fail )
         + ( negative_near_configuration_BR * far_muon_pass )
           * ( tau_pair_to_jet_electron_BR * near_pion_fail
               + tau_pair_to_electron_pair_BR * near_muon_fail ) );

    *no_jets_no_leptons
    += ( ( negative_near_configuration_BR + positive_near_configuration_BR )
         * ( tau_pair_to_jet_pair_BR * near_pion_fail * far_pion_fail
             + ( tau_pair_to_jet_muon_BR + tau_pair_to_jet_electron_BR )
               * ( near_pion_fail * far_muon_fail
                   + near_muon_fail * far_pion_fail )
             + ( tau_pair_to_muon_pair_BR + 2.0 * tau_pair_to_muon_electron_BR
                 + tau_pair_to_electron_pair_BR )
               * near_muon_fail * far_muon_fail ) );

  }

}  // end of LHC_FASER namespace.
