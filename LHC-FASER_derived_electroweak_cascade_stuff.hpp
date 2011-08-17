/*
 * LHC-FASER_derived_electroweak_cascade_stuff.hpp
 *
 *  Created on: May 9, 2011
 *      Author: bol
 */

#ifndef LHC_FASER_DERIVED_ELECTROWEAK_CASCADE_STUFF_HPP_
#define LHC_FASER_DERIVED_ELECTROWEAK_CASCADE_STUFF_HPP_

#include "LHC-FASER_base_electroweak_cascade_stuff.hpp"

namespace LHC_FASER
{

  /* this is a derived class that implements the channel where the
   * "decaying neutralino" is actually the lightest neutralino & hence assumed
   * to be stable
   */
  class lightestNeutralinoCascade : public electroweakCascade
  {
  public:
    lightestNeutralinoCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    virtual
    ~lightestNeutralinoCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    virtual bool
    validSignal( int const numberOfJets,
                 int const numberOfNegativeElectrons,
                 int const numberOfPositiveElectrons,
                 int const numberOfNegativeMuons,
                 int const numberOfPositiveMuons )
    // this returns false except for the case of 0 jets & 0 leptons, when it
    // returns true.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual void
    calculateAcceptance( acceptanceValues* const currentAcceptance )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };


  /* this is a derived class that implements the channel
   * decaying neutralino -> smuon or selectron -> lightest neutralino.
   * it also includes ... -> antismuon or antiselectron -> ... assuming that
   * the charge-conjugate process is identical for these purposes. it does
   * *not*, however, account for both chiralities of sfermion.
   */
  class neutralinoToSemuCascade : public electroweakCascade
  {
  public:
    neutralinoToSemuCascade( leptonAcceptanceParameterSet* const kinematics,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                             //bool const coloredDecayerIsNotAntiparticle,
                CppSLHA::particle_property_set const* const electroweakDecayer,
                             //bool const electroweakDecayerIsNotAntiparticle,
               CppSLHA::particle_property_set const* const intermediateDecayer,
                             input_handler const* const shortcut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    virtual
    ~neutralinoToSemuCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    bool muonsNotElectrons;
    leptonEnergyDistribution* nearSameDistribution;
    leptonEnergyDistribution* farSameDistribution;
    leptonEnergyDistribution* nearOppositeDistribution;
    leptonEnergyDistribution* farOppositeDistribution;
    double jetLeftHandedness;
    double nearNegativeLeptonLeftHandedness;
    double nearNegativeLeptonSameHandednessTimesBr;
    /* the above are for the case of a positively-charged up-type quark or
     * negatively-charged down-type quark as the jet.
     * nearNegativeLeptonSameHandednessTimesBr is the branching ratio for
     * the cascade times the proportion of decays where the diagram with a
     * negatively-charged near lepton has the same helicity as the jet. for the
     * diagram with a positively-charged near lepton, this number gives the
     * branching ratio times the proportion with opposite helicity.
     */

    double nearSamePass;
    double nearSameFail;
    double farSamePass;
    double farSameFail;
    double nearOppositePass;
    double nearOppositeFail;
    double farOppositePass;
    double farOppositeFail;
    double justOneNegativeLeptonPass;
    double justOnePositiveLeptonPass;
    // these are for holding the pass & fail rates for whichever configuration
    // is being calculated.

    virtual bool
    validSignal( int const numberOfJets,
                 int const numberOfNegativeElectrons,
                 int const numberOfPositiveElectrons,
                 int const numberOfNegativeMuons,
                 int const numberOfPositiveMuons )
    // this returns true if numberOfJets == 0, & either 1 of or both signs of
    // electron OR muon are asked for.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual void
    calculateAcceptance( acceptanceCutSet* const cuts,
                         acceptanceValues* const currentAcceptance )
    // this returns the appropriate acceptances multiplied by branching ratios
    // from the electroweakino through the selectron or smuon to the LSP.
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };


  /* this is a derived class that implements the channel
   * decaying neutralino -> smuon or selectron -> lightest neutralino.
   * it also includes ... -> antismuon or antiselectron -> ... assuming that
   * the charge-conjugate process is identical for these purposes. it does
   * *not*, however, account for both chiralities of sfermion.
   */
  class chargeSummedNeutralinoToSemuCascade : public electroweakCascade
  {

  public:
    chargeSummedNeutralinoToSemuCascade(
                                leptonAcceptanceParameterSet* const kinematics,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                CppSLHA::particle_property_set const* const electroweakDecayer,
               CppSLHA::particle_property_set const* const intermediateDecayer,
                                          input_handler const* const shortcut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    virtual
    ~chargeSummedNeutralinoToSemuCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    bool muonsNotElectrons;
    leptonEnergyDistribution* nearDistribution;
    leptonEnergyDistribution* farDistribution;

    double nearPass;
    double nearFail;
    double farPass;
    double farFail;
    double bothPass;
    double onePassOneFail;
    // these are for holding the pass & fail rates for whichever configuration
    // is being calculated.

    virtual bool
    validSignal( int const numberOfJets,
                 int const numberOfNegativeElectrons,
                 int const numberOfPositiveElectrons,
                 int const numberOfNegativeMuons,
                 int const numberOfPositiveMuons )
    // this returns true if numberOfJets == 0, & either 1 of or both signs of
    // electron OR muon are asked for.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual void
    calculateAcceptance( acceptanceCutSet* const cuts,
                         acceptanceValues* const currentAcceptance )
    // this returns the appropriate acceptances multiplied by branching ratios
    // from the electroweakino through the selectron or smuon to the LSP.
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };


  /* this is a derived class that implements the channel
   * decaying neutralino -> stau -> lightest neutralino.
   * it also includes ... -> antistau -> ... assuming that
   * the charge-conjugate process is identical for these purposes. it does
   * *not*, however, account for both stau mass eigenstates.
   */
  class neutralinoToStauCascade : public electroweakCascade
  {
  public:
    static double const tauPairToPionPairBr;
    static double const tauToPionTimesTauToElectronBr;
    static double const tauToPionTimesTauToMuonBr;
    static double const tauPairToElectronPairBr;
    static double const tauToElectronTimesTauToMuonBr;
    static double const tauPairToMuonPairBr;

    neutralinoToStauCascade( leptonAcceptanceParameterSet* const kinematics,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                CppSLHA::particle_property_set const* const electroweakDecayer,
               CppSLHA::particle_property_set const* const intermediateDecayer,
                             input_handler const* const shortcut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~neutralinoToStauCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    leptonEnergyDistribution* nearSameTauDistribution;
    leptonEnergyDistribution* farSameTauDistribution;
    leptonEnergyDistribution* nearOppositeTauDistribution;
    leptonEnergyDistribution* farOppositeTauDistribution;
    leptonEnergyDistribution* nearSameHardMuonDistribution;
    leptonEnergyDistribution* nearSameSoftMuonDistribution;
    leptonEnergyDistribution* nearSameHardPionDistribution;
    leptonEnergyDistribution* nearSameSoftPionDistribution;
    leptonEnergyDistribution* nearOppositeHardMuonDistribution;
    leptonEnergyDistribution* nearOppositeSoftMuonDistribution;
    leptonEnergyDistribution* nearOppositeHardPionDistribution;
    leptonEnergyDistribution* nearOppositeSoftPionDistribution;
    leptonEnergyDistribution* farSameHardMuonDistribution;
    leptonEnergyDistribution* farSameSoftMuonDistribution;
    leptonEnergyDistribution* farSameHardPionDistribution;
    leptonEnergyDistribution* farSameSoftPionDistribution;
    leptonEnergyDistribution* farOppositeHardMuonDistribution;
    leptonEnergyDistribution* farOppositeSoftMuonDistribution;
    leptonEnergyDistribution* farOppositeHardPionDistribution;
    leptonEnergyDistribution* farOppositeSoftPionDistribution;
    leptonEnergyDistribution* currentNearMuonDistribution;
    leptonEnergyDistribution* currentNearPionDistribution;
    leptonEnergyDistribution* currentFarMuonDistribution;
    leptonEnergyDistribution* currentFarPionDistribution;
    acceptanceCutSet* currentCuts;
    acceptanceValues* currentAcceptance;
    double jetLeftHandedness;
    double jetRightHandedness;
    double nearNegativeTauLeftHandedness;
    double farNegativeTauLeftHandedness;
    double nearNegativeTauRightHandedness;
    double farNegativeTauRightHandedness;
    // the above are for the case of a positively-charged up-type quark or
    // negatively-charged down-type quark as the jet.
    double negativeNearConfigurationBr;
    double positiveNearConfigurationBr;
    // these are the branching ratios for the particular polarization & charge
    // configurations of the tau leptons being calculated.

    double nearMuonPass;
    double nearMuonFail;
    double farMuonPass;
    double farMuonFail;
    double nearPionPass;
    double nearPionFail;
    double farPionPass;
    double farPionFail;
    double bothPass;
    double currentPass;
    double currentNearFailSum;
    double currentFarFailSum;
    // these are for holding the pass & fail rates for whichever configuration
    // is being calculated.

    virtual bool
    validSignal( int const numberOfJets,
                 int const numberOfNegativeElectrons,
                 int const numberOfPositiveElectrons,
                 int const numberOfNegativeMuons,
                 int const numberOfPositiveMuons )
    /* this returns true if a configuration where each of the signs of tau
     * lepton decayed either into a detected jet, detected lepton, or
     * undetected particle, & false otherwise.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual void
    calculateAcceptance( acceptanceCutSet* const cuts,
                         acceptanceValues* const currentAcceptance )
    // this returns the appropriate acceptances multiplied by branching ratios
    // from the electroweakino through the selectron or smuon to the LSP.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    calculateForCurrentConfiguration()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };

  /* this is a derived class that implements the channel
   * decaying neutralino -> stau -> lightest neutralino.
   * it also includes ... -> antistau -> ... assuming that
   * the charge-conjugate process is identical for these purposes. it does
   * *not*, however, account for both stau mass eigenstates.
   */
  class chargeSummedNeutralinoToStauCascade : public electroweakCascade
  {
  public:
    chargeSummedNeutralinoToStauCascade(
                                leptonAcceptanceParameterSet* const kinematics,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                CppSLHA::particle_property_set const* const electroweakDecayer,
               CppSLHA::particle_property_set const* const intermediateDecayer,
                                         input_handler const* const shortcut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~chargeSummedNeutralinoToStauCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    leptonEnergyDistribution* nearTauDistribution;
    leptonEnergyDistribution* farTauDistribution;
    leptonEnergyDistribution* nearHardMuonDistribution;
    leptonEnergyDistribution* nearSoftMuonDistribution;
    leptonEnergyDistribution* nearHardPionDistribution;
    leptonEnergyDistribution* nearSoftPionDistribution;
    leptonEnergyDistribution* farHardMuonDistribution;
    leptonEnergyDistribution* farSoftMuonDistribution;
    leptonEnergyDistribution* farHardPionDistribution;
    leptonEnergyDistribution* farSoftPionDistribution;
    leptonEnergyDistribution* currentNearMuonDistribution;
    leptonEnergyDistribution* currentNearPionDistribution;
    leptonEnergyDistribution* currentFarMuonDistribution;
    leptonEnergyDistribution* currentFarPionDistribution;
    acceptanceCutSet* currentCuts;
    acceptanceValues* currentAcceptance;
    double jetLeftHandedness;
    double jetRightHandedness;
    double nearNegativeTauLeftHandedness;
    double farNegativeTauLeftHandedness;
    double nearNegativeTauRightHandedness;
    double farNegativeTauRightHandedness;
    // the above are for the case of a positively-charged up-type quark or
    // negatively-charged down-type quark as the jet.
    double configurationBr;
    // this is the branching ratio for the particular polarization & charge
    // configuration of the tau leptons being calculated.

    double nearMuonPass;
    double nearMuonFail;
    double farMuonPass;
    double farMuonFail;
    double nearPionPass;
    double nearPionFail;
    double farPionPass;
    double farPionFail;
    double bothPass;
    double currentPass;
    double currentNearFailSum;
    double currentFarFailSum;
    // these are for holding the pass & fail rates for whichever configuration
    // is being calculated.

    virtual bool
    validSignal( int const numberOfJets,
                 int const numberOfNegativeElectrons,
                 int const numberOfPositiveElectrons,
                 int const numberOfNegativeMuons,
                 int const numberOfPositiveMuons )
    /* this returns true if a configuration where each of the signs of tau
     * lepton decayed either into a detected jet, detected lepton, or
     * undetected particle, & false otherwise.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual void
    calculateAcceptance( acceptanceCutSet* const cuts,
                         acceptanceValues* const currentAcceptance )
    // this returns the appropriate acceptances multiplied by branching ratios
    // from the electroweakino through the selectron or smuon to the LSP.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    calculateForCurrentConfiguration()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };

  // this is a derived class that implements the channel
  // decaying neutralino -> Z -> lightest neutralino
  class neutralino_to_Z_cascade : public channel_calculator
  {

  public:

    neutralino_to_Z_cascade( readier_for_new_point* const given_readier,
                             bool const should_sum_charges,
                             double const given_primary_cut,
                             double const given_secondary_cut,
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
                             double* const given_no_jets_no_leptons )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~neutralino_to_Z_cascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:

    bool const should_sum_charges_flag;

    leptonEnergyDistribution* direct_distribution;

    leptonEnergyDistribution* hard_muon_distribution;
    leptonEnergyDistribution* soft_muon_distribution;
    leptonEnergyDistribution* hard_pion_distribution;
    leptonEnergyDistribution* soft_pion_distribution;

    double const negative_tau_left_handedness;

    double negative_near_configuration_BR;
    double positive_near_configuration_BR;
    /* these are the branching ratios for the particular polarization & charge
     * configurations of the tau leptons being calculated. for the
     * charge-summed case, just negative_near_configuration_BR is used & it
     * covers both charges.
     */

    double near_muon_pass;
    double near_muon_fail;
    double far_muon_pass;
    double far_muon_fail;
    double near_pion_pass;
    double near_pion_fail;
    double far_pion_pass;
    double far_pion_fail;
    double both_pass;
    double current_pass;
    // these are for holding the pass & fail rates for whichever configuration
    // is being calculated.

    double* const two_jets_no_leptons;
    double* const one_jet_one_negative_electron;
    double* const one_jet_one_positive_electron;
    double* const one_jet_one_negative_muon;
    double* const one_jet_one_positive_muon;
    double* const OSSF_minus_OSDF_leptons;
    double* const no_jets_two_electrons;
    double* const no_jets_two_muons;
    double* const no_jets_one_positive_muon_one_negative_electron;
    double* const no_jets_one_negative_muon_one_positive_electron;

    double* const one_jet_no_leptons;
    double* const no_jets_one_negative_electron;
    double* const no_jets_one_positive_electron;
    double* const no_jets_one_negative_muon;
    double* const no_jets_one_positive_muon;

    double* const no_jets_no_leptons;


    void
    calculate()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    charge_summed_calculate()
    /* this updates values based on the current near & far muon/electron & pion
     * distributions & configuration branching ratio, summing over tau lepton
     * charges.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    charge_distinguished_calculate()
    /* this updates values based on the current near & far muon/electron & pion
     * distributions & configuration branching ratio, separately recording
     * values for differing charges.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };



  // this is a derived class that implements the channel
  // decaying neutralino -> Z -> lightest neutralino
  class charge_summed_neutralino_to_Z : public channel_calculator
  {

  public:

    charge_summed_neutralino_to_Z( double const given_primary_cut,
                                   double const given_secondary_cut,
                                   double const given_jet_cut,
                            lepton_acceptance_value* const given_kinematics,
           CppSLHA::particle_property_set const* const given_decaying_scolored,
                                 bool const given_scolored_is_not_antiparticle,
              CppSLHA::particle_property_set const* const given_decaying_EWino,
                                    bool const given_EWino_is_not_antiparticle,
          CppSLHA::particle_property_set const* const given_mediating_particle,
                                   input_handler const* const given_shortcuts,
                                   double* const given_two_jets_no_leptons,
                                   double* const given_one_jet_one_muon,
                                   double* const given_one_jet_one_electron,
                                   double* const given_OSSF_minus_OSDF_leptons,
                                   double* const given_no_jets_two_muons,
                                   double* const given_no_jets_two_electrons,
                             double* const given_no_jets_one_muon_one_electron,
                                   double* const given_one_jet_no_leptons,
                                   double* const given_no_jets_one_muon,
                                   double* const given_no_jets_one_electron,
                                   double* const given_no_jets_no_leptons )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~charge_summed_neutralino_to_Z()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:

    double* const two_jets_no_leptons;
    double* const one_jet_one_muon;
    double* const one_jet_one_electron;
    double* const OSSF_minus_OSDF_leptons;
    double* const no_jets_two_muons;
    double* const no_jets_two_electrons;
    double* const no_jets_one_muon_one_electron;

    double* const one_jet_no_leptons;
    double* const no_jets_one_muon;
    double* const no_jets_one_electron;

    double* const no_jets_no_leptons;


    void
    calculate()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };

  // this is a derived class that implements the channel
  // decaying neutralino -> neutral EWSB scalar -> lightest neutralino
  class neutralino_to_Higgs_cascade : public channel_calculator
  {

  public:

    neutralino_to_Higgs( double const given_primary_cut,
                         double const given_secondary_cut,
                         double const given_jet_cut,
                         lepton_acceptance_value* const given_kinematics,
           CppSLHA::particle_property_set const* const given_decaying_scolored,
                         bool const given_scolored_is_not_antiparticle,
              CppSLHA::particle_property_set const* const given_decaying_EWino,
                         bool const given_EWino_is_not_antiparticle,
          CppSLHA::particle_property_set const* const given_mediating_particle,
                         input_handler const* const given_shortcuts,
                         double* const given_two_jets_no_leptons,
                         double* const given_one_jet_one_negative_muon,
                         double* const given_one_jet_one_negative_electron,
                         double* const given_one_jet_one_positive_muon,
                         double* const given_one_jet_one_positive_electron,
                         double* const given_OSSF_minus_OSDF_leptons,
                         double* const given_no_jets_two_muons,
                         double* const given_no_jets_two_electrons,
           double* const given_no_jets_one_negative_muon_one_positive_electron,
           double* const given_no_jets_one_positive_muon_one_negative_electron,
                         double* const given_one_jet_no_leptons,
                         double* const given_no_jets_one_negative_muon,
                         double* const given_no_jets_one_negative_electron,
                         double* const given_no_jets_one_positive_muon,
                         double* const given_no_jets_one_positive_electron,
                         double* const given_no_jets_no_leptons )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~neutralino_to_Higgs_cascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    double* const two_jets_no_leptons;
    double* const one_jet_one_negative_muon;
    double* const one_jet_one_negative_electron;
    double* const one_jet_one_positive_muon;
    double* const one_jet_one_positive_electron;
    double* const OSSF_minus_OSDF_leptons;
    double* const no_jets_two_muons;
    double* const no_jets_two_electrons;
    double* const no_jets_one_negative_muon_one_positive_electron;
    double* const no_jets_one_positive_muon_one_negative_electron;

    double* const one_jet_no_leptons;
    double* const no_jets_one_negative_muon;
    double* const no_jets_one_negative_electron;
    double* const no_jets_one_positive_muon;
    double* const no_jets_one_positive_electron;

    double* const no_jets_no_leptons;


    void
    calculate()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };

  // this is a derived class that implements the channel
  // decaying neutralino -> lightest neutralino via a 3-body decay.
  class charge_summed_neutralino_three_body : public channel_calculator
  {

  public:

    charge_summed_neutralino_three_body( double const given_primary_cut,
                                         double const given_secondary_cut,
                                         double const given_jet_cut,
                               lepton_acceptance_value* const given_kinematics,
           CppSLHA::particle_property_set const* const given_decaying_scolored,
                                 bool const given_scolored_is_not_antiparticle,
              CppSLHA::particle_property_set const* const given_decaying_EWino,
                                    bool const given_EWino_is_not_antiparticle,
                                    input_handler const* const given_shortcuts,
                                       double* const given_two_jets_no_leptons,
                                         double* const given_one_jet_one_muon,
                                      double* const given_one_jet_one_electron,
                                   double* const given_OSSF_minus_OSDF_leptons,
                                         double* const given_no_jets_two_muons,
                                     double* const given_no_jets_two_electrons,
                             double* const given_no_jets_one_muon_one_electron,
                                        double* const given_one_jet_no_leptons,
                                         double* const given_no_jets_one_muon,
                                      double* const given_no_jets_one_electron,
                                     double* const given_no_jets_no_leptons )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~charge_summed_neutralino_three_body()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:

    double* const two_jets_no_leptons;
    double* const one_jet_one_muon;
    double* const one_jet_one_electron;
    double* const OSSF_minus_OSDF_leptons;
    double* const no_jets_two_muons;
    double* const no_jets_two_electrons;
    double* const no_jets_one_muon_one_electron;

    double* const one_jet_no_leptons;
    double* const no_jets_one_muon;
    double* const no_jets_one_electron;

    double* const no_jets_no_leptons;


    void
    calculate()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };


  /* this is a derived class that implements the channel
   * decaying neutralino -> lightest neutralino via a 3-body decay into tau
   * leptons which then decay into light leptons or jets.
   */
  class neutralino_three_body_to_taus : public channel_calculator
  {

  public:

    neutralino_three_body_to_taus( double const given_primary_cut,
                                         double const given_secondary_cut,
                                         double const given_jet_cut,
                               lepton_acceptance_value* const given_kinematics,
           CppSLHA::particle_property_set const* const given_decaying_scolored,
                                 bool const given_scolored_is_not_antiparticle,
              CppSLHA::particle_property_set const* const given_decaying_EWino,
                                    bool const given_EWino_is_not_antiparticle,
                                    input_handler const* const given_shortcuts,
                                       double* const given_two_jets_no_leptons,
                                         double* const given_one_jet_one_muon,
                                      double* const given_one_jet_one_electron,
                                   double* const given_OSSF_minus_OSDF_leptons,
                                         double* const given_no_jets_two_muons,
                                     double* const given_no_jets_two_electrons,
                             double* const given_no_jets_one_muon_one_electron,
                                        double* const given_one_jet_no_leptons,
                                         double* const given_no_jets_one_muon,
                                      double* const given_no_jets_one_electron,
                                     double* const given_no_jets_no_leptons )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~neutralino_three_body_to_taus()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    double* const two_jets_no_leptons;
    double* const one_jet_one_negative_electron;
    double* const one_jet_one_positive_electron;
    double* const one_jet_one_negative_muon;
    double* const one_jet_one_positive_muon;
    double* const no_jets_two_electrons;
    double* const no_jets_two_muons;
    double* const no_jets_one_positive_muon_one_negative_electron;
    double* const no_jets_one_negative_muon_one_positive_electron;

    double* const one_jet_no_leptons;
    double* const no_jets_one_negative_electron;
    double* const no_jets_one_positive_electron;
    double* const no_jets_one_negative_muon;
    double* const no_jets_one_positive_muon;

    double* const no_jets_no_leptons;


    void
    calculate()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };



  // this is a derived class that implements the channel
  // decaying chargino -> smuon or selectron -> lightest neutralino
  class chargino_to_semu : public channel_calculator
  {

  public:

    chargino_to_semu( double const given_primary_cut,
                      double const given_secondary_cut,
                      lepton_acceptance_value* const given_kinematics,
           CppSLHA::particle_property_set const* const given_decaying_scolored,
                      bool const given_scolored_is_not_antiparticle,
              CppSLHA::particle_property_set const* const given_decaying_EWino,
                      bool const given_EWino_is_not_antiparticle,
          CppSLHA::particle_property_set const* const given_mediating_particle,
                      input_handler const* const given_shortcuts,
                      double* const given_no_jets_one_lepton,
                      double* const given_no_jets_no_leptons )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~chargino_to_semu()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:

    double* const no_jets_one_lepton;

    double* const no_jets_no_leptons;

    void
    calculate()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };

  // this is a derived class that implements the channel
  // decaying chargino -> stau -> lightest neutralino
  class chargino_to_stau : public channel_calculator
  {

  public:

    chargino_to_stau( double const given_primary_cut,
                      double const given_secondary_cut,
                      double const given_jet_cut,
                      lepton_acceptance_value* const given_kinematics,
           CppSLHA::particle_property_set const* const given_decaying_scolored,
                      bool const given_scolored_is_not_antiparticle,
              CppSLHA::particle_property_set const* const given_decaying_EWino,
                      bool const given_EWino_is_not_antiparticle,
          CppSLHA::particle_property_set const* const given_mediating_particle,
                      input_handler const* const given_shortcuts,
                      double* const given_one_jet_no_leptons,
                      double* const given_no_jets_one_muon,
                      double* const given_no_jets_one_electron,
                      double* const given_no_jets_no_leptons )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~chargino_to_stau()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:

    double* const one_jet_no_leptons;
    double* const no_jets_one_muon;
    double* const no_jets_one_electron;

    double* const no_jets_no_leptons;


    void
    calculate()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };

  // this is a derived class that implements the channel
  // decaying chargino -> muon or electron sneutrino -> lightest neutralino
  class chargino_to_emu_sneutrino : public channel_calculator
  {

  public:

    chargino_to_emu_sneutrino( double const given_primary_cut,
                               double const given_secondary_cut,
                               lepton_acceptance_value* const given_kinematics,
           CppSLHA::particle_property_set const* const given_decaying_scolored,
                               bool const given_scolored_is_not_antiparticle,
              CppSLHA::particle_property_set const* const given_decaying_EWino,
                               bool const given_EWino_is_not_antiparticle,
          CppSLHA::particle_property_set const* const given_mediating_particle,
                               input_handler const* const given_shortcuts,
                               double* const given_no_jets_one_lepton,
                               double* const given_no_jets_no_leptons )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~chargino_to_emu_sneutrino()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:

    double* const no_jets_one_lepton;

    double* const no_jets_no_leptons;

    void
    calculate()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };

  // this is a derived class that implements the channel
  // decaying chargino -> tau sneutrino -> lightest neutralino
  class chargino_to_tau_sneutrino : public channel_calculator
  {

  public:

    chargino_to_tau_sneutrino( double const given_primary_cut,
                               double const given_secondary_cut,
                               double const given_jet_cut,
                               lepton_acceptance_value* const given_kinematics,
           CppSLHA::particle_property_set const* const given_decaying_scolored,
                               bool const given_scolored_is_not_antiparticle,
              CppSLHA::particle_property_set const* const given_decaying_EWino,
                               bool const given_EWino_is_not_antiparticle,
          CppSLHA::particle_property_set const* const given_mediating_particle,
                               input_handler const* const given_shortcuts,
                               double* const given_one_jet_no_leptons,
                               double* const given_no_jets_one_muon,
                               double* const given_no_jets_one_electron,
                               double* const given_no_jets_no_leptons )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~chargino_to_tau_sneutrino()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:

    double* const one_jet_no_leptons;
    double* const no_jets_one_muon;
    double* const no_jets_one_electron;

    double* const no_jets_no_leptons;


    void
    calculate()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };

  // this is a derived class that implements the channel
  // decaying chargino -> W -> lightest neutralino
  class chargino_to_W : public channel_calculator
  {

  public:

    chargino_to_W( double const given_primary_cut,
                   double const given_secondary_cut,
                   double const given_jet_cut,
                   lepton_acceptance_value* const given_kinematics,
           CppSLHA::particle_property_set const* const given_decaying_scolored,
                   bool const given_scolored_is_not_antiparticle,
              CppSLHA::particle_property_set const* const given_decaying_EWino,
                   bool const given_EWino_is_not_antiparticle,
          CppSLHA::particle_property_set const* const given_mediating_particle,
                   input_handler const* const given_shortcuts,
                   double* const given_two_jets_no_leptons,
                   double* const given_one_jet_no_leptons,
                   double* const given_no_jets_one_muon,
                   double* const given_no_jets_one_electron,
                   double* const given_no_jets_no_leptons )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~chargino_to_W()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:

    double* const two_jets_no_leptons;

    double* const one_jet_no_leptons;
    double* const no_jets_one_muon;
    double* const no_jets_one_electron;

    double* const no_jets_no_leptons;


    void
    calculate()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };

  // this is a derived class that implements the channel
  // decaying chargino -> charged EWSB scalar -> lightest neutralino
  class chargino_to_Higgs : public channel_calculator
  {

  public:

    chargino_to_Higgs( double const given_primary_cut,
                       double const given_secondary_cut,
                       double const given_jet_cut,
                       lepton_acceptance_value* const given_kinematics,
           CppSLHA::particle_property_set const* const given_decaying_scolored,
                       bool const given_scolored_is_not_antiparticle,
              CppSLHA::particle_property_set const* const given_decaying_EWino,
                       bool const given_EWino_is_not_antiparticle,
          CppSLHA::particle_property_set const* const given_mediating_particle,
                       input_handler const* const given_shortcuts,
                       double* const given_two_jets_no_leptons,
                       double* const given_one_jet_no_leptons,
                       double* const given_no_jets_one_muon,
                       double* const given_no_jets_one_electron,
                       double* const given_no_jets_no_leptons )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~chargino_to_Higgs()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:

    double* const two_jets_no_leptons;

    double* const one_jet_no_leptons;
    double* const no_jets_one_muon;
    double* const no_jets_one_electron;

    double* const no_jets_no_leptons;


    void
    calculate()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };

  // this is a derived class that implements the channel
  // decaying chargino -> lightest neutralino via a 3-body decay.
  class chargino_three_body : public channel_calculator
  {

  public:

    chargino_three_body( double const given_primary_cut,
                         double const given_secondary_cut,
                         double const given_jet_cut,
                         lepton_acceptance_value* const given_kinematics,
           CppSLHA::particle_property_set const* const given_decaying_scolored,
                         bool const given_scolored_is_not_antiparticle,
              CppSLHA::particle_property_set const* const given_decaying_EWino,
                         bool const given_EWino_is_not_antiparticle,
                         input_handler const* const given_shortcuts,
                         double* const given_two_jets_no_leptons,
                         double* const given_one_jet_no_leptons,
                         double* const given_no_jets_one_muon,
                         double* const given_no_jets_one_electron,
                         double* const given_no_jets_no_leptons )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~chargino_three_body()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:

    double* const two_jets_no_leptons;

    double* const one_jet_no_leptons;
    double* const no_jets_one_muon;
    double* const no_jets_one_electron;

    double* const no_jets_no_leptons;


    void
    calculate()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };






  // inline functions:

  inline bool
  lightestNeutralinoCascade::validSignal( int const numberOfJets,
                                          int const numberOfNegativeElectrons,
                                          int const numberOfPositiveElectrons,
                                          int const numberOfNegativeMuons,
                                          int const numberOfPositiveMuons )
  // this returns false except for the case of 0 jets & 0 leptons, when it
  // returns true.
  {
    if( ( 0 == numberOfJets )
        &&
        ( 0 == numberOfNegativeElectrons )
        &&
        ( 0 == numberOfPositiveElectrons )
        &&
        ( 0 == numberOfNegativeMuons )
        &&
        ( 0 == numberOfPositiveMuons ) )
    {
      return true;
      // once the cascade has reached the LSP, it stays as the LSP 100% of the
      // time, never triggering any jets or leptons.
    }
    else
    {
      return false;
    }
  }

  inline void
  lightestNeutralinoCascade::calculateAcceptance(
                                    acceptanceValues* const currentAcceptance )
  {
    currentAcceptance->setZeroJetsZeroLeptons( 1.0 );
    // once the cascade has reached the LSP, it stays as the LSP 100% of the
    // time, never triggering any jets or leptons.
  }


  inline bool
  neutralinoToSemuCascade::validSignal( int const numberOfJets,
                                        int const numberOfNegativeElectrons,
                                        int const numberOfPositiveElectrons,
                                        int const numberOfNegativeMuons,
                                        int const numberOfPositiveMuons )
  // this returns true if numberOfJets == 0, & either 1 of or both signs of
  // electron OR muon are asked for.
  {
    if( ( 0 == numberOfJets )
        // if we are asked for no jets...
        &&
        ( ( ( 0 == numberOfNegativeElectrons )
            &&
            ( 0 == numberOfPositiveElectrons ) )
          ||
          ( ( 0 == numberOfNegativeMuons )
              &&
              ( 0 == numberOfPositiveMuons ) ) ) )
        // AND if we are asked for EITHER electrons OR muons but NOT both...
    {
      return true;
    }
    else
    {
      return false;
    }
  }



  inline bool
  chargeSummedNeutralinoToSemuCascade::validSignal( int const numberOfJets,
                                           int const numberOfNegativeElectrons,
                                           int const numberOfPositiveElectrons,
                                               int const numberOfNegativeMuons,
                                              int const numberOfPositiveMuons )
  // this returns true if numberOfJets == 0, & either 1 of or both signs of
  // electron OR muon are asked for.
  {
    if( ( 0 == numberOfJets )
        // if we are asked for no jets...
        &&
        ( ( ( 0 == numberOfNegativeElectrons )
            &&
            ( 0 == numberOfPositiveElectrons ) )
          ||
          ( ( 0 == numberOfNegativeMuons )
              &&
              ( 0 == numberOfPositiveMuons ) ) ) )
        // AND if we are asked for EITHER electrons OR muons but NOT both...
    {
      return true;
    }
    else
    {
      return false;
    }
  }



  inline bool
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
}  // end of LHC_FASER namespace.

#endif /* LHC_FASER_DERIVED_ELECTROWEAK_CASCADE_STUFF_HPP_ */
