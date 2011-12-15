/*
 * LHC-FASER_charged_electroweak_cascade_stuff.hpp
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

#ifndef LHC_FASER_CHARGED_ELECTROWEAK_CASCADE_STUFF_HPP_
#define LHC_FASER_CHARGED_ELECTROWEAK_CASCADE_STUFF_HPP_

#include "LHC-FASER_base_electroweak_cascade_stuff.hpp"

namespace LHC_FASER
{
  /* this is a derived class that implements the channel
   * chargino -> smuon or selectron or electron sneutrino or muon sneutrino
   * -> lightest neutralino.
   */
  class charginoToSemuOrEmuSnuCascade : public electroweakCascade
  {
  public:
    charginoToSemuOrEmuSnuCascade(
                                leptonAcceptanceParameterSet* const kinematics,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                                   particlePointer const coloredDecayer,
                                   particlePointer const electroweakDecayer,
                                   particlePointer const intermediateDecayer,
                                   inputHandler const* const inputShortcut,
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
    calculateAcceptance( acceptanceCutSet const* const currentCuts,
                         acceptanceValues* const currentAcceptance )
    // this returns the appropriate acceptancesPerCutSet multiplied by branching ratios
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
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                                   particlePointer const coloredDecayer,
                                   particlePointer const electroweakDecayer,
                                   particlePointer const intermediateDecayer,
                                   inputHandler const* const inputShortcut,
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
    calculateAcceptance( acceptanceCutSet const* const currentCuts,
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
                        effectiveSquarkMassHolder* const effectiveSquarkMass,
                        particlePointer const coloredDecayer,
                        particlePointer const electroweakDecayer,
                        inputHandler const* const inputShortcut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~charginoToWCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    bool downTypeQuark;
    leptonEnergyDistribution* directJetDistribution;
    leptonEnergyDistribution* leftHandedJetRightHandedAntimuonDistribution;
    leptonEnergyDistribution* rightHandedJetRightHandedAntimuonDistribution;
    leptonEnergyDistribution* leftHandedJetTauMuonDistribution;
    leptonEnergyDistribution* rightHandedJetTauMuonDistribution;
    leptonEnergyDistribution* leftHandedJetTauPionDistribution;
    leptonEnergyDistribution* rightHandedJetTauPionDistribution;
    leptonEnergyDistribution* currentMuonDistribution;
    leptonEnergyDistribution* currentPionDistribution;
    acceptanceCutSet const* currentCuts;
    acceptanceValues* currentAcceptance;
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
    calculateAcceptance( acceptanceCutSet const* const currentCuts,
                         acceptanceValues* const currentAcceptance )
    // this returns the appropriate acceptancesPerCutSet multiplied by branching ratios
    // from the electroweakino through the selectron or smuon to the LSP.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    calculateForCurrentConfiguration()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };


  // this is a derived class that implements the channel
  // chargino -> neutral spin-0 boson + lightest neutralino
  class charginoToHiggsCascade : public electroweakCascade
  {
  public:
    charginoToHiggsCascade( leptonAcceptanceParameterSet* const kinematics,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                            particlePointer const coloredDecayer,
                            particlePointer const electroweakDecayer,
                            particlePointer const intermediateDecayer,
                            inputHandler const* const inputShortcut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~charginoToHiggsCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    leptonEnergyDistribution* directAntimuonDistribution;
    leptonEnergyDistribution* antitauAntimuonDistribution;
    leptonEnergyDistribution* antitauPionDistribution;
    acceptanceCutSet const* currentCuts;
    acceptanceValues* currentAcceptance;
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
    CppSLHA::particle_decay_set_handler const* intermediateDecayerDecays;
    double currentBrToHadrons;
    double currentBrToAntielectron;
    double currentBrToAntimuon;
    double currentBrToAntitau;

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
    // this returns the appropriate acceptancesPerCutSet multiplied by branching ratios
    // from the electroweakino through the selectron or smuon to the LSP.
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };


  // this is a derived class that implements the channel
  // chargino -> SM fermion+antifermion pair + lightest neutralino
  class charginoVirtualCascade : public electroweakCascade
  {
  public:
    charginoVirtualCascade( leptonAcceptanceParameterSet* const kinematics,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                            particlePointer const coloredDecayer,
                            particlePointer const electroweakDecayer,
                            inputHandler const* const inputShortcut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~charginoVirtualCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    leptonEnergyDistribution* directDownUpDistribution;
    leptonEnergyDistribution* directStrangeCharmDistribution;
    leptonEnergyDistribution* directElectronDistribution;
    leptonEnergyDistribution* directMuonDistribution;
    leptonEnergyDistribution* directTauDistribution;
    leptonEnergyDistribution* hardMuonDistribution;
    leptonEnergyDistribution* softMuonDistribution;
    leptonEnergyDistribution* hardPionDistribution;
    leptonEnergyDistribution* softPionDistribution;
    leptonEnergyDistribution* currentJetDistribution;
    acceptanceCutSet const* currentCuts;
    acceptanceValues* currentAcceptance;
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
    CppSLHA::particle_decay_set_handler const* electroweakDecayerDecays;

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
    // this returns the appropriate acceptancesPerCutSet multiplied by branching ratios
    // from the electroweakino through the selectron or smuon to the LSP.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    calculateForCurrentJetConfiguration()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };



  // this is a derived class that implements the channel
  // squark/gluino -> W + other squark/gluino
  class scoloredToWPlusScoloredCascade : public electroweakCascade
  {
  public:
    scoloredToWPlusScoloredCascade(
                                leptonAcceptanceParameterSet* const kinematics,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                                    particlePointer const coloredDecayer,
                                    particlePointer const electroweakDecayer,
                                    particlePointer const lighterScolored,
                                    inputHandler const* const inputShortcut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~scoloredToWPlusScoloredCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    leptonEnergyDistribution* directMuonDistribution;
    leptonEnergyDistribution* antitauAntimuonDistribution;
    leptonEnergyDistribution* antitauPionDistribution;
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
    int wPdgCode;

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
  };





  // inline functions:
  // none.

}  // end of LHC_FASER namespace.

#endif /* LHC_FASER_CHARGED_ELECTROWEAK_CASCADE_STUFF_HPP_ */
