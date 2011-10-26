/*
 * LHC-FASER_neutral_electroweak_cascade_stuff.hpp
 *
 *  Created on: May 9, 2011
 *      Author: bol
 */

#ifndef LHC_FASER_NEUTRAL_ELECTROWEAK_CASCADE_STUFF_HPP_
#define LHC_FASER_NEUTRAL_ELECTROWEAK_CASCADE_STUFF_HPP_

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
    lightestNeutralinoCascade( inputHandler const* const inputShortcut );
    virtual
    ~lightestNeutralinoCascade();

  protected:
    virtual bool
    validSignal( int const numberOfJets,
                 int const numberOfNegativeElectrons,
                 int const numberOfPositiveElectrons,
                 int const numberOfNegativeMuons,
                 int const numberOfPositiveMuons );
    // this returns false except for the case of 0 jets & 0 leptons, when it
    // returns true.

    virtual void
    calculateAcceptance( acceptanceCutSet const* const currentCuts,
                         acceptanceValues* const currentAcceptance );
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
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                             particlePointer const coloredDecayer,
                             //bool const coloredDecayerIsNotAntiparticle,
                             particlePointer const electroweakDecayer,
                             //bool const electroweakDecayerIsNotAntiparticle,
                             particlePointer const intermediateDecayer,
                             inputHandler const* const shortcut )
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
    calculateAcceptance( acceptanceCutSet const* const currentCuts,
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
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                                         particlePointer const coloredDecayer,
                                      particlePointer const electroweakDecayer,
                                     particlePointer const intermediateDecayer,
                                         inputHandler const* const shortcut )
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
    calculateAcceptance( acceptanceCutSet const* const currentCuts,
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
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                             particlePointer const coloredDecayer,
                             particlePointer const electroweakDecayer,
                             particlePointer const intermediateDecayer,
                             inputHandler const* const shortcut )
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
    acceptanceCutSet const* currentCuts;
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
    calculateAcceptance( acceptanceCutSet const* const currentCuts,
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
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                                         particlePointer const coloredDecayer,
                                      particlePointer const electroweakDecayer,
                                     particlePointer const intermediateDecayer,
                                         inputHandler const* const shortcut )
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
    acceptanceCutSet const* currentCuts;
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
    calculateAcceptance( acceptanceCutSet const* const currentCuts,
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
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                          particlePointer const coloredDecayer,
                          particlePointer const electroweakDecayer,
                          inputHandler const* const shortcut )
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
    acceptanceCutSet const* currentCuts;
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
    calculateAcceptance( acceptanceCutSet const* const currentCuts,
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
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                                      particlePointer const coloredDecayer,
                                      particlePointer const electroweakDecayer,
                                      inputHandler const* const shortcut )
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
    acceptanceCutSet const* currentCuts;
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
    calculateAcceptance( acceptanceCutSet const* const currentCuts,
                         acceptanceValues* const currentAcceptance )
    // this returns the appropriate acceptances multiplied by branching ratios
    // from the electroweakino through the selectron or smuon to the LSP.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    calculateForCurrentConfiguration()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };


  // this is a derived class that implements the channel
  // decaying neutralino -> neutral spin-0 boson + lightest neutralino
  class neutralinoToHiggsCascade : public electroweakCascade
  {
  public:
    neutralinoToHiggsCascade( leptonAcceptanceParameterSet* const kinematics,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                              particlePointer const coloredDecayer,
                              particlePointer const electroweakDecayer,
                              particlePointer const intermediateDecayer,
                              inputHandler const* const shortcut )
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
    acceptanceCutSet const* currentCuts;
    acceptanceValues* currentAcceptance;
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
    calculateAcceptance( acceptanceCutSet const* const currentCuts,
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
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                                          particlePointer const coloredDecayer,
                                      particlePointer const electroweakDecayer,
                                          inputHandler const* const shortcut )
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
    acceptanceCutSet const* currentCuts;
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
    calculateAcceptance( acceptanceCutSet const* const currentCuts,
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
  class neutralinoVirtualCascade : public chargeSummedNeutralinoVirtualCascade
  {
  public:
    neutralinoVirtualCascade( leptonAcceptanceParameterSet* const kinematics,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                              particlePointer const coloredDecayer,
                              particlePointer const electroweakDecayer,
                              inputHandler const* const shortcut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~neutralinoVirtualCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  //protected:
    // nothing.
  };


  // this is a derived class that implements the channel
  // decaying neutralino -> Z + lightest neutralino
  class scoloredToZPlusScoloredCascade : public electroweakCascade
  {
  public:
    scoloredToZPlusScoloredCascade(
                                leptonAcceptanceParameterSet* const kinematics,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                                    particlePointer const coloredDecayer,
                                    particlePointer const lighterScolored,
                                    inputHandler const* const shortcut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~scoloredToZPlusScoloredCascade()
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
    acceptanceCutSet const* currentCuts;
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
    calculateAcceptance( acceptanceCutSet const* const currentCuts,
                         acceptanceValues* const currentAcceptance )
    // this returns the appropriate acceptances multiplied by branching ratios
    // from the electroweakino through the selectron or smuon to the LSP.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    calculateForCurrentConfiguration()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };


  // this is a derived class that implements the channel
  // decaying neutralino -> spin-0 neutral boson + lightest neutralino
  class scoloredToHiggsPlusScoloredCascade : public electroweakCascade
  {
  public:
    scoloredToHiggsPlusScoloredCascade(
                                leptonAcceptanceParameterSet* const kinematics,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                                        particlePointer const coloredDecayer,
                                        particlePointer const lighterScolored,
                                     particlePointer const intermediateDecayer,
                                        inputHandler const* const shortcut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~scoloredToHiggsPlusScoloredCascade()
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
    acceptanceCutSet const* currentCuts;
    acceptanceValues* currentAcceptance;
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
    calculateAcceptance( acceptanceCutSet const* const currentCuts,
                         acceptanceValues* const currentAcceptance )
    // this returns the appropriate acceptances multiplied by branching ratios
    // from the electroweakino through the selectron or smuon to the LSP.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    calculateForCurrentConfiguration()
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
                                     acceptanceCutSet const* const currentCuts,
                                    acceptanceValues* const currentAcceptance )
  {
    currentAcceptance->setZeroJetsZeroLeptons( 1.0 );
    // once the cascade has reached the LSP, it stays as the LSP 100% of the
    // time, never triggering any jets or leptons.
  }

}  // end of LHC_FASER namespace.

#endif /* LHC_FASER_NEUTRAL_ELECTROWEAK_CASCADE_STUFF_HPP_ */
