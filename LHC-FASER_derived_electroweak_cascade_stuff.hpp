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
  // decaying neutralino -> Z + lightest neutralino
  class neutralinoToZCascade : public electroweakCascade
  {
  public:
    neutralinoToZCascade( leptonAcceptanceParameterSet* const kinematics,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                CppSLHA::particle_property_set const* const electroweakDecayer,
                          input_handler const* const shortcut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~neutralinoToZCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    leptonEnergyDistribution* directMuonDistribution;
    leptonEnergyDistribution* sameHandedTauDistribution;
    leptonEnergyDistribution* oppositeHandedTauDistribution;
    leptonEnergyDistribution* sameHardMuonDistribution;
    leptonEnergyDistribution* sameSoftMuonDistribution;
    leptonEnergyDistribution* sameHardPionDistribution;
    leptonEnergyDistribution* sameSoftPionDistribution;
    leptonEnergyDistribution* oppositeHardMuonDistribution;
    leptonEnergyDistribution* oppositeSoftMuonDistribution;
    leptonEnergyDistribution* oppositeHardPionDistribution;
    leptonEnergyDistribution* oppositeSoftPionDistribution;
    leptonEnergyDistribution* currentNegativeMuonDistribution;
    leptonEnergyDistribution* currentPositivePionDistribution;
    leptonEnergyDistribution* currentPositiveMuonDistribution;
    leptonEnergyDistribution* currentNegativePionDistribution;
    double negativeTauLeftHandedness;
    double jetLeftHandedness;
    double directMuonPass;
    double directMuonFail;
    double directJetPass;
    double directJetFail;
    double negativeTauMuonPass;
    double negativeTauMuonFail;
    double negativeTauPionPass;
    double negativeTauPionFail;
    double positiveTauMuonPass;
    double positiveTauMuonFail;
    double positiveTauPionPass;
    double positiveTauPionFail;
    double currentPass;
    // these are for holding the pass & fail rates for whichever configuration
    // is being calculated.
    double configurationBr;
    // this is the branching ratio for the particular polarization & charge
    // configuration of the tau leptons being calculated.

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
  // decaying neutralino -> Z + lightest neutralino
  class chargeSummedNeutralinoToZCascade : public electroweakCascade
  {
  public:
    chargeSummedNeutralinoToZCascade(
                                leptonAcceptanceParameterSet* const kinematics,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                CppSLHA::particle_property_set const* const electroweakDecayer,
                          input_handler const* const shortcut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~chargeSummedNeutralinoToZCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    leptonEnergyDistribution* directMuonDistribution;
    leptonEnergyDistribution* hardMuonDistribution;
    leptonEnergyDistribution* softMuonDistribution;
    leptonEnergyDistribution* hardPionDistribution;
    leptonEnergyDistribution* softPionDistribution;
    leptonEnergyDistribution* currentMuonDistribution;
    leptonEnergyDistribution* currentPionDistribution;
    double negativeTauLeftHandedness;
    double directMuonPass;
    double directMuonFail;
    double directJetPass;
    double directJetFail;
    double tauMuonPass;
    double tauMuonFail;
    double tauPionPass;
    double tauPionFail;
    double currentPass;
    // these are for holding the pass & fail rates for whichever configuration
    // is being calculated.
    double configurationBr;
    // this is the branching ratio for the particular polarization & charge
    // configuration of the tau leptons being calculated.

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
  // decaying neutralino -> spin-0 boson + lightest neutralino
  class neutralinoToHiggsCascade : public electroweakCascade
  {
  public:
    neutralinoToHiggsCascade( leptonAcceptanceParameterSet* const kinematics,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                CppSLHA::particle_property_set const* const electroweakDecayer,
               CppSLHA::particle_property_set const* const intermediateDecayer,
                              input_handler const* const shortcut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~neutralinoToHiggsCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    leptonEnergyDistribution* directMuonDistribution;
    leptonEnergyDistribution* hardMuonDistribution;
    leptonEnergyDistribution* softMuonDistribution;
    leptonEnergyDistribution* hardPionDistribution;
    leptonEnergyDistribution* softPionDistribution;
    leptonEnergyDistribution* currentNegativeMuonDistribution;
    leptonEnergyDistribution* currentPositivePionDistribution;
    leptonEnergyDistribution* currentPositiveMuonDistribution;
    leptonEnergyDistribution* currentNegativePionDistribution;
    double negativeTauLeftHandedness;
    double directMuonPass;
    double directMuonFail;
    double directJetPass;
    double directJetFail;
    double negativeTauMuonPass;
    double negativeTauMuonFail;
    double negativeTauPionPass;
    double negativeTauPionFail;
    double positiveTauMuonPass;
    double positiveTauMuonFail;
    double positiveTauPionPass;
    double positiveTauPionFail;
    double currentPass;
    // these are for holding the pass & fail rates for whichever configuration
    // is being calculated.
    double configurationBr;
    // this is the branching ratio for the particular polarization & charge
    // configuration of the tau leptons being calculated.
    CppSLHA::particle_decay_set_handler const* intermediateDecayerDecays;
    double currentBrToHadrons;
    double currentBrToElectrons;
    double currentBrToMuons;
    double currentBrToTaus;

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
  // decaying neutralino -> SM fermion+antifermion pair + lightest neutralino
  class chargeSummedNeutralinoVirtualCascade : public electroweakCascade
  {
  public:
    chargeSummedNeutralinoVirtualCascade(
                                leptonAcceptanceParameterSet* const kinematics,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                CppSLHA::particle_property_set const* const electroweakDecayer,
                                          input_handler const* const shortcut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~chargeSummedNeutralinoVirtualCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    leptonEnergyDistribution* directDownDistribution;
    leptonEnergyDistribution* directUpDistribution;
    leptonEnergyDistribution* directStrangeDistribution;
    leptonEnergyDistribution* directCharmDistribution;
    leptonEnergyDistribution* directBottomDistribution;
    leptonEnergyDistribution* directElectronDistribution;
    leptonEnergyDistribution* directMuonDistribution;
    leptonEnergyDistribution* directTauDistribution;
    leptonEnergyDistribution* hardMuonDistribution;
    leptonEnergyDistribution* softMuonDistribution;
    leptonEnergyDistribution* hardPionDistribution;
    leptonEnergyDistribution* softPionDistribution;
    leptonEnergyDistribution* currentJetDistribution;
    double negativeTauLeftHandedness;
    double jetLeftHandedness;
    double directPass;
    double directFail;
    double hardMuonPass;
    double hardMuonFail;
    double softMuonPass;
    double softMuonFail;
    double averageMuonPass;
    double averageMuonFail;
    double hardPionPass;
    double hardPionFail;
    double softPionPass;
    double softPionFail;
    double averagePionPass;
    double averagePionFail;
    double currentPass;
    // these are for holding the pass & fail rates for whichever configuration
    // is being calculated.
    double configurationBr;
    // this is the branching ratio for the particular polarization & charge
    // configuration of the tau leptons being calculated.
    CppSLHA::particle_decay_set_handler const* electroweakDecayerDecays;
    double currentBrToHadrons;
    double currentBrToElectrons;
    double currentBrToMuons;
    double currentBrToTaus;
    acceptanceCutSet* currentCuts;
    acceptanceValues* currentAcceptance;

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
    calculateForCurrentJetConfiguration()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };

  // this is a derived class that implements the channel
  // decaying neutralino -> SM fermion+antifermion pair + lightest neutralino
  // currently this is a placeholder, just substituting in the charge-summed
  // version.
  class neutralinoThreeBodyCascade :
                                    public chargeSummedNeutralinoVirtualCascade
  {
  public:
    neutralinoThreeBodyCascade( leptonAcceptanceParameterSet* const kinematics,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                CppSLHA::particle_property_set const* const electroweakDecayer,
                                input_handler const* const shortcut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~neutralinoThreeBodyCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  //protected:
    // nothing.
  };


  /* this is a derived class that implements the channel
   * chargino -> smuon or selectron or electron snuetrino or muon sneutrino
   * -> lightest neutralino.
   */
  class charginoToSemuOrEmuSnuCascade : public electroweakCascade
  {
  public:
    charginoToSemuOrEmuSnuCascade(
                                leptonAcceptanceParameterSet* const kinematics,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                CppSLHA::particle_property_set const* const electroweakDecayer,
               CppSLHA::particle_property_set const* const intermediateDecayer,
                                   input_handler const* const shortcut,
                                   bool const sneutrinoVersion )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    virtual
    ~charginoToSemuOrEmuSnuCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    bool muonsNotElectrons;
    bool const sneutrinoVersion;
    leptonEnergyDistribution* sameDistribution;
    leptonEnergyDistribution* oppositeDistribution;
    double jetLeftHandedness;
    double nearLeptonLeftHandedness;
    double nearLeptonSameHandednessTimesBr;
    /* the above are for the case of a positively-charged up-type quark or
     * positively-charged down-type antiquark as the jet, since the chargino
     * assumed is that with the positive PDG code, hence in this case the
     * positive charge.
     */

    double samePass;
    double sameFail;
    double oppositePass;
    double oppositeFail;
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


  // this is a derived class that implements the channel
  // chargino -> stau or tau snuetrino -> lightest neutralino
  class charginoToStauOrTauSnuCascade : public electroweakCascade
  {
  public:
    charginoToStauOrTauSnuCascade(
                                leptonAcceptanceParameterSet* const kinematics,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                CppSLHA::particle_property_set const* const electroweakDecayer,
               CppSLHA::particle_property_set const* const intermediateDecayer,
                                   input_handler const* const shortcut,
                                   bool const sneutrinoVersion )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    virtual
    ~charginoToStauOrTauSnuCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    bool const sneutrinoVersion;
    leptonEnergyDistribution* sameTauDistribution;
    leptonEnergyDistribution* oppositeTauDistribution;
    leptonEnergyDistribution* sameHardMuonDistribution;
    leptonEnergyDistribution* sameSoftMuonDistribution;
    leptonEnergyDistribution* sameHardPionDistribution;
    leptonEnergyDistribution* sameSoftPionDistribution;
    leptonEnergyDistribution* oppositeHardMuonDistribution;
    leptonEnergyDistribution* oppositeSoftMuonDistribution;
    leptonEnergyDistribution* oppositeHardPionDistribution;
    leptonEnergyDistribution* oppositeSoftPionDistribution;
    double jetLeftHandedness;
    double nearLeptonLeftHandedness;
    double antitauRightHandedness;
    /* the above are for the case of a positively-charged up-type quark or
     * positively-charged down-type antiquark as the jet, since the chargino
     * assumed is that with the positive PDG code, hence in this case the
     * positive charge.
     */

    double muonPass;
    double muonFail;
    double pionPass;
    double pionFail;
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




  // this is a derived class that implements the channel
  // chargino -> W + lightest neutralino
  class charginoToWCascade : public electroweakCascade
  {
  public:
    charginoToWCascade( leptonAcceptanceParameterSet* const kinematics,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                CppSLHA::particle_property_set const* const electroweakDecayer,
                          input_handler const* const shortcut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~charginoToWCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    leptonEnergyDistribution* directJetDistribution;
    leptonEnergyDistribution* directMuonDistribution;
    leptonEnergyDistribution* sameHandedTauDistribution;
    leptonEnergyDistribution* oppositeHandedTauDistribution;
    leptonEnergyDistribution* sameHardMuonDistribution;
    leptonEnergyDistribution* sameSoftMuonDistribution;
    leptonEnergyDistribution* sameHardPionDistribution;
    leptonEnergyDistribution* sameSoftPionDistribution;
    leptonEnergyDistribution* oppositeHardMuonDistribution;
    leptonEnergyDistribution* oppositeSoftMuonDistribution;
    leptonEnergyDistribution* oppositeHardPionDistribution;
    leptonEnergyDistribution* oppositeSoftPionDistribution;
    leptonEnergyDistribution* currentNegativeMuonDistribution;
    leptonEnergyDistribution* currentPositivePionDistribution;
    leptonEnergyDistribution* currentPositiveMuonDistribution;
    leptonEnergyDistribution* currentNegativePionDistribution;
    double jetLeftHandedness;
    double directMuonPass;
    double directMuonFail;
    double directJetPass;
    double directJetFail;
    double antitauAntimuonPass;
    double antitauAntimuonFail;
    double antitauPionPass;
    double antitauPionFail;
    double currentPass;
    // these are for holding the pass & fail rates for whichever configuration
    // is being calculated.
    double configurationBr;
    // this is the branching ratio for the particular polarization & charge
    // configuration of the tau leptons being calculated.

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
}  // end of LHC_FASER namespace.

#endif /* LHC_FASER_DERIVED_ELECTROWEAK_CASCADE_STUFF_HPP_ */
