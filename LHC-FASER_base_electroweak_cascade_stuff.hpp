/*
 * LHC-FASER_cascades_stuff.hpp
 *
 *  Created on: 17 Dec 2010
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

#ifndef LHC_FASER_BASE_ELECTROWEAK_CASCADE_STUFF_HPP_
#define LHC_FASER_BASE_ELECTROWEAK_CASCADE_STUFF_HPP_

#include "LHC-FASER_kinematics_stuff.hpp"
#include "LHC-FASER_derived_lepton_distributions.hpp"

namespace LHC_FASER
{

  /* this is a class to hold the various acceptances that can come from a
   * cascade decay (not necessarily ending with the LSP) which results in 2 SM
   * fermions.
   */
  class acceptanceValues
  {

  public:

    acceptanceValues( double const defaultUnsetValues )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~acceptanceValues()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    reset( double const unsetValues )
    // this sets all values to unsetValues & alreadyCalculatedFlag to false.
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    bool
    notAlreadyCalculated()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    flagAsAlreadyCalculated()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getTwoJets()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setTwoJets( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getOneJetOneNegativeElectron()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setOneJetOneNegativeElectron( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getOneJetOnePositiveElectron()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setOneJetOnePositiveElectron( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getOneJetOneNegativeMuon()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setOneJetOneNegativeMuon( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getOneJetOnePositiveMuon()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setOneJetOnePositiveMuon( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getOneJetZeroLeptons()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setOneJetZeroLeptons( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getOssfMinusOsdf()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setOssfMinusOsdf( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getElectronPlusAntielectron()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setElectronPlusAntielectron( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getNegativeElectronPlusPositiveMuon()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setNegativeElectronPlusPositiveMuon( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getNegativeMuonPlusPositiveElectron()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setNegativeMuonPlusPositiveElectron( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getMuonPlusAntimuon()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setMuonPlusAntimuon( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getZeroJetsOneNegativeElectron()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setZeroJetsOneNegativeElectron( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getZeroJetsOnePositiveElectron()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setZeroJetsOnePositiveElectron( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getZeroJetsOneNegativeMuon()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setZeroJetsOneNegativeMuon( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getZeroJetsOnePositiveMuon()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setZeroJetsOnePositiveMuon( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getZeroJetsZeroLeptons()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setZeroJetsZeroLeptons( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:
    bool alreadyCalculated;
    double twoJets;
    double oneJetOneNegativeElectron;
    double oneJetOnePositiveElectron;
    double oneJetOneNegativeMuon;
    double oneJetOnePositiveMuon;
    double oneJetZeroLeptons;
    double ossfMinusOsdf;
    double electronPlusAntielectron;
    double negativeElectronPlusPositiveMuon;
    double negativeMuonPlusPositiveElectron;
    double muonPlusAntimuon;
    double zeroJetsOneNegativeElectron;
    double zerJetsOnePositiveElectron;
    double zeroJetsOneNegativeMuon;
    double zeroJetsOnePositiveMuon;
    double zeroJetsZeroLeptons;
  };


  /* this is an abstract base class for objects that prepare lepton
   * distributions & cache values for their distribution combined with given
   * acceptance cuts.
   *
   * it isn't actually intended as a kind of interface to allow multiple
   * different derived classes to all appear the same to other objects, more
   * just as a convenient way to avoid reproducing common code.
   */
  class electroweakCascade : public readied_for_new_point
  {
  public:
    electroweakCascade( leptonAcceptanceCut* const kinematics,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                        bool const coloredDecayerIsNotAntiparticle,
                CppSLHA::particle_property_set const* const electroweakDecayer,
                        bool const electroweakDecayerIsNotAntiparticle,
              CppSLHA::particle_property_set const* const intermediateDecayer,
                        bool const canDoOssfMinusOsdf,
                        input_handler const* const shortcut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~electroweakCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    input_handler const*
    getShortcut()
    {
      return shortcut;
    }

    double
    getAcceptance( acceptanceCutSet const* const cuts,
                   int const numberOfAdditionalJets,
                   int const numberOfNegativeElectrons,
                   int const numberOfPositiveElectrons,
                   int const numberOfNegativeMuons,
                   int const numberOfPositiveMuons )
    /* this does the common job of checking to see if the point has been
     * updated before calling the relevant protected virtual function, which
     * returns the set of values for acceptance * branching ratio FOR THE
     * CASCADE FROM electroweakDecayer ONWARDS (does NOT include the BR for
     * coloredDecayer to electroweakDecayer - this is dealt with by the
     * cutSpecifiedFullCascade, because it depends on what flavors of quark jet
     * are being considered for the signal).
     * this function returns 0.0 if the cascade never results in the requested
     * set of jets + leptons.
     * the number of jets is the number produced from the decay of the
     * decaying *electroweak* particle (neutralino, chargino or vector boson),
     * NOT from the decay of the colored sparticle to the decaying electroweak
     * particle.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    getOssfMinusOsdf( acceptanceCutSet const* const cuts )
    /* this does the common job of checking to see if the point has been
     * updated before calling the relevant protected virtual function, which
     * returns the value for acceptance * branching ratio FOR THE
     * CASCADE FROM electroweakDecayer ONWARDS (does NOT include the BR for
     * coloredDecayer to electroweakDecayer - this is dealt with by the
     * cutSpecifiedFullCascade, because it depends on what flavors of quark jet
     * are being considered for the signal).
     * this function returns 0.0 if the cascade never results in an OSSF - OSDF
     * signal.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    static void
    cachePairConstruction( acceptanceValues** const cachedValuesAndFlag,
                           acceptanceCutSet** const cachedKey,
                           acceptanceCutSet* const constructionKey )
    /* this makes a new acceptanceValues instance to hold an acceptance value
     * along with a flag for whether it has been calculated or not, &
     * associates it with a key.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    static void
    cachePairReset( acceptanceValues* const cachedValuesAndFlag,
                    acceptanceCutSet* const cachedKey,
                    acceptanceCutSet* const constructionKey )
    /* this resets an acceptanceValues instance to hold an acceptance value
     * along with a flag for whether it has been calculated or not, &
     * associates it with a key.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    leptonAcceptanceCut* const kinematics;
    CppSLHA::particle_property_set const* const coloredDecayer;
    bool const coloredDecayerIsNotAntiparticle;
    CppSLHA::particle_property_set const* const electroweakDecayer;
    bool const electroweakDecayerIsNotAntiparticle;
    CppSLHA::particle_property_set const* const intermediateDecayer;
    bool const canDoOssfMinusOsdf;
    exclusive_BR_calculator* firstBr;
    exclusive_BR_calculator* secondBr;
    double cascadeBr;
    input_handler const* const shortcut;
    acceptanceValues* currentAcceptance;
    minimalAllocationKeyedVector< acceptanceValues,
                                  acceptanceCutSet > acceptances;
    std::vector< leptonEnergyDistribution* > activeDistributions;

    virtual bool
    validSignal( int const numberOfJets,
                 int const numberOfNegativeElectrons,
                 int const numberOfPositiveElectrons,
                 int const numberOfNegativeMuons,
                 int const numberOfPositiveMuons )
    = 0;
    // this should return true if the derived class can provide a non-zero
    // value for the requested signal, false otherwise.

    virtual void
    calculateAcceptance( acceptanceCutSet* const cuts,
                         acceptanceValues* const currentAcceptance )
    = 0;
    // this should set all the appropriate acceptance * branching ratio values
    // for the cascade.

    double
    integrateAcceptance( leptonEnergyDistribution* const leptonDistribution,
                          double const transverseMomentumCut )
    /* this numerically integrates the product of leptonDistribution with the
     * lepton acceptance histogram from kinematics, scaled to having
     * transverseMomentumCut as its transverse momentum cut.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };

  /* this is a class to hold a set of electroweakCascades associated with a
   * single decaying electroweakino or vector boson (from the decay of a single
   * colored sparticle). the decaying colored sparticle is assumed to be the
   * version with a PDG code which is POSITIVE, rather than the antiparticle
   * with negative code. hence for example an instance of this class would be
   * for a positively-charged up squark which decays into a negatively-charged
   * down quark & a positively-charged chargino, or for another example an
   * instance could be for a negatively-charged down squark, which decays into
   * a positively-charged up quark & a negatively-charged chargino (which has
   * negative PDG code!). the charge-conjugate process, which is assumed to
   * have the same branching ratios & acceptances, though for opposite-charge
   * SM fermions, must be taken care of by whatever objects are using this
   * class.
   */
  class electroweakCascadeSet
  {
  public:
    electroweakCascadeSet( leptonAcceptanceCut* const kinematics,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                           bool const coloredDecayerIsNotAntiparticle,
                CppSLHA::particle_property_set const* const electroweakDecayer,
                           bool const electroweakDecayerIsNotAntiparticle,
                           input_handler const* const shortcut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~electroweakCascadeSet()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    input_handler const*
    getShortcut()
    {
      return shortcut;
    }

    double
    getAcceptance( acceptanceCutSet const* const cuts,
                   int const numberOfAdditionalJets,
                   int const numberOfNegativeElectrons,
                   int const numberOfPositiveElectrons,
                   int const numberOfNegativeMuons,
                   int const numberOfPositiveMuons )
    /* this does the common job of checking to see if the point has been
     * updated before calling the relevant protected virtual function, which
     * returns the set of values for acceptance * branching ratio FOR THE
     * CASCADES FROM electroweakDecayer ONWARDS (does NOT include the BR for
     * coloredDecayer to electroweakDecayer - this is dealt with by the
     * cutSpecifiedFullCascade, because it depends on what flavors of quark jet
     * are being considered for the signal).
     * this function returns 0.0 if the cascades never result in the requested
     * set of jets + leptons.
     * the number of jets is the number produced from the decay of the
     * decaying *electroweak* particle (neutralino, chargino or vector boson),
     * NOT from the decay of the colored sparticle to the decaying electroweak
     * particle.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    getOssfMinusOsdf( acceptanceCutSet const* const cuts )
    /* this does the common job of checking to see if the point has been
     * updated before calling the relevant protected virtual function, which
     * returns the value for acceptance * branching ratio FOR THE
     * CASCADES FROM electroweakDecayer ONWARDS (does NOT include the BR for
     * coloredDecayer to electroweakDecayer - this is dealt with by the
     * cutSpecifiedFullCascade, because it depends on what flavors of quark jet
     * are being considered for the signal).
     * this function returns 0.0 if the cascades never result in an OSSF - OSDF
     * signal.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    leptonAcceptanceCut* const kinematics;
    CppSLHA::particle_property_set const* const coloredDecayer;
    bool const coloredDecayerIsNotAntiparticle;
    CppSLHA::particle_property_set const* const electroweakDecayer;
    bool const electroweakDecayerIsNotAntiparticle;
    input_handler const* const shortcut;
    std::vector< electroweakCascade* > cascades;
    electroweakCascade* currentCascade; // this is used for filling cascades.
  };


  /* this is an abstract base class for objects that calculate a set of
   * acceptances & update their respective values.
   */
  class channel_calculator : public readied_for_new_point
  {

  public:

    channel_calculator( double const given_primary_cut,
                        double const given_secondary_cut,
                        double const given_jet_cut,
                        lepton_acceptance_value* const given_kinematics,
           CppSLHA::particle_property_set const* const given_decaying_scolored,
                        bool const given_scolored_is_not_antiparticle,
              CppSLHA::particle_property_set const* const given_decaying_EWino,
                        bool const given_EWino_is_not_antiparticle,
          CppSLHA::particle_property_set const* const given_mediating_particle,
                        input_handler const* const given_shortcuts )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~channel_calculator()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    update_numbers()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    double const primary_cut;
    double const secondary_cut;
    double const jet_cut;
    lepton_acceptance_value* const kinematics;
    CppSLHA::particle_property_set const* const decaying_scolored;
    bool const scolored_is_not_antiparticle;
    CppSLHA::particle_property_set const* const decaying_EWino;
    bool const EWino_is_not_antiparticle;
    CppSLHA::particle_property_set const* const mediating_particle;
    input_handler const* const shortcut;

    exclusive_BR_calculator* first_BR;
    exclusive_BR_calculator* second_BR;
    double cascade_BR;
    // these refer to the branching ratios for the 1st & 2nd vertices of the
    // "main" (ignoring tau decay vertices) Feynmann diagrams of the cascade.

    std::vector< leptonEnergyDistribution* > active_distributions;
    // this is just to keep track of which lepton_energy_distribution pointers
    // have been allocated memory.


    /* this controls how fine the summation of the trapezia is for the
     * approximation of the integral of the product of the lepton energy
     * distribution with the acceptance:
     */
    static int const integration_bins = 10;

    // these are for ease of calculating branching ratios of distributions with
    // tau leptons:
    static double const tau_pair_to_muon_pair_BR
    = ( CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR
        * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR );
    static double const tau_pair_to_muon_electron_BR
    = ( CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR
        * CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR );
    static double const tau_pair_to_electron_pair_BR
    = ( CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR
        * CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR );
    static double const tau_pair_to_jet_pair_BR
    = ( CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
        * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR );
    static double const tau_pair_to_jet_muon_BR
    = ( CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
        * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR );
    static double const tau_pair_to_jet_electron_BR
    = ( CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
        * CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR );

    virtual void
    calculate()
    = 0
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    integrate_acceptance( double const given_cut,
                  leptonEnergyDistribution const* const lepton_distribution )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };


  // this holds pointers to all the channel_calculator instances associated
  // with a particular cascade acceptance value as well as the value itself.
  class cascade_acceptance_value : public readied_for_new_point
  {

  public:

    cascade_acceptance_value( double const given_starting_value,
                              readier_for_new_point* const given_readier )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~cascade_acceptance_value()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    double
    get_value()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double* const
    get_value_pointer()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    add_channel( channel_calculator* given_channel )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    double acceptance_value;

    std::vector< channel_calculator* > channels;

  };


  /* this sets up all the open channel_calculators & registers them with their
   * appropriate cascade_acceptance_values, uses the cascade_acceptance_values
   * to return values for various acceptances.
   * this does not keep track of different charges of lepton separately. I
   * intend to write another class to do that some day, when a relevant signal
   * needs it.
   */
  class charge_summed_cascade_calculator
  {

  public:

    charge_summed_cascade_calculator( double const given_primary_cut,
                                      double const given_secondary_cut,
                                      double const given_jet_cut,
                               lepton_acceptance_value* const given_kinematics,
           CppSLHA::particle_property_set const* const given_decaying_scolored,
                                 bool const given_scolored_is_not_antiparticle,
              CppSLHA::particle_property_set const* const given_decaying_EWino,
                                    bool const given_EWino_is_not_antiparticle,
                                   input_handler const* const given_shortcuts )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~charge_summed_cascade_calculator()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    bool
    is_requested( double const given_primary_cut,
                  double const given_secondary_cut,
                  double const given_jet_cut,
                  lepton_acceptance_value* const given_kinematics,
           CppSLHA::particle_property_set const* const given_decaying_scolored,
                  bool const given_scolored_is_not_antiparticle,
              CppSLHA::particle_property_set const* const given_decaying_EWino,
                  bool const given_EWino_is_not_antiparticle )
    /* this returns true if the relevant cuts, the kinematics, & the sparticles
     * match those held. relevant cuts are those over 0.0, any cuts of 0.0 or
     * less are assumed to be irrelevant to the signal requesting these values,
     * so this calculator can be given to the signal.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_OSSF_minus_OSDF_acceptance()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_acceptance( int const number_of_jets,
                    int const number_of_leptons )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_muon_acceptance( int const number_of_jets,
                         int const number_of_muons )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_electron_acceptance( int const number_of_jets,
                             int const number_of_electrons )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    std::vector< channel_calculator* > channels;

    cascade_acceptance_value OSSF_minus_OSDF_leptons;
    cascade_acceptance_value two_jets_no_leptons;
    cascade_acceptance_value one_jet_one_muon;
    cascade_acceptance_value one_jet_one_electron;
    cascade_acceptance_value no_jets_two_muons;
    cascade_acceptance_value no_jets_two_electrons;
    cascade_acceptance_value no_jets_one_muon_one_electron;

    cascade_acceptance_value one_jet_no_leptons;
    cascade_acceptance_value no_jets_one_muon;
    cascade_acceptance_value no_jets_one_electron;

    cascade_acceptance_value no_jets_no_leptons;

    // this stuff needs to be stored to identify the
    // charge_summed_cascade_calculator
    double const primary_cut;
    double const secondary_cut;
    double const jet_cut;
    lepton_acceptance_value* const kinematics;
    CppSLHA::particle_property_set const* const decaying_scolored;
    bool const scolored_is_not_antiparticle;
    CppSLHA::particle_property_set const* const decaying_EWino;
    bool const EWino_is_not_antiparticle;

  };


  // this class passes out pointers to cascade_calculators for given
  // kinematics & decaying sparticles.
  class cascade_handler : public readied_for_new_point
  {

  public:

    cascade_handler( input_handler const* const given_shortcuts )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~cascade_handler()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    charge_summed_cascade_calculator*
    get_charge_summed_calculator( double given_primary_cut,
                                  double given_secondary_cut,
                                  double given_jet_cut,
                               lepton_acceptance_value* const given_kinematics,
                 CppSLHA::particle_property_set const* const decaying_scolored,
                                  bool const scolored_is_not_antiparticle,
                    CppSLHA::particle_property_set const* const decaying_EWino,
                                  bool const EWino_is_not_antiparticle )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:

    input_handler const* const shortcut;

    std::vector< charge_summed_cascade_calculator* > calculators;

  };





  // inline functions:

  inline input_handler const*
  electroweakCascade::getShortcut()
  {
    return shortcut;
  }

  inline acceptanceValues const*
  electroweakCascade::getOssfMinusOsdf( acceptanceCutSet const* const cuts )
  /* this does the common job of checking to see if the point has been
   * updated before calling the relevant protected virtual function, which
   * returns the set of values for acceptance * branching ratio FOR THE
   * CASCADE FROM electroweakDecayer ONWARDS (does NOT include the BR for
   * coloredDecayer to electroweakDecayer - this is dealt with by the
   * cutSpecifiedFullCascade, because it depends on what flavors of quark jet
   * are being considered for the signal).
   * this function returns NULL if the cascade never results in the requested
   * set of jets + leptons.
   */
  {
    if( canDoOssfMinusOsdf )
    {
      if( needs_to_prepare_for_this_point() )
      {
        acceptances.clear();
      }
      currentAcceptance = acceptances.getPointer( cuts );
      if( currentAcceptance->notAlreadyCalculated() )
      {
        calculateAcceptance( currentAcceptance );
        currentAcceptance->flagAsAlreadyCalculated();
      }
      return currentAcceptance->getOssfMinusOsdf();
    }
    else
    {
      return 0.0;
    }
  }

  inline void
  electroweakCascade::cachePairConstruction(
                                  acceptanceValues** const cachedValuesAndFlag,
                                            acceptanceCutSet** const cachedKey,
                                      acceptanceCutSet* const constructionKey )
  // this makes a new std::pair to hold an acceptance value along with a flag
  // for whether it has been calculated or not, & associates it with a key.
  {
    *cachedValuesAndFlag
    = new acceptanceValues( CppSLHA::CppSLHA_global::really_wrong_value );
    // the new pair has the default "unset" values.
    *cachedKey = new acceptanceCutSet( constructionKey );
  }

  inline void
  electroweakCascade::cachePairReset(
                                   acceptanceValues* const cachedValuesAndFlag,
                                      acceptanceCutSet* const cachedKey,
                                      acceptanceCutSet* const constructionKey )
  // this resets a std::pair to hold an acceptance value along with a flag
  // for whether it has been calculated or not, & associates it with a key.
  {
    cachedValuesAndFlag->reset();
    // the new pair has the default "unset" values.
    cachedKey->becomeCopyOf( constructionKey );
  }

  inline acceptanceValues const*
  electroweakCascade::protectedGetAcceptance( int const numberOfJets,
                                           int const numberOfNegativeElectrons,
                                           int const numberOfPositiveElectrons,
                                              int const numberOfNegativeMuons,
                                              int const numberOfPositiveMuons )
  {
    return NULL;
  }

  inline acceptanceValues const*
  electroweakCascade::protectedGetOssfMinusOsdf()
  {
    return NULL;
  }




  inline void
  channel_calculator::update_numbers()
  {

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: channel_calculator::update_numbers() called."
    << std::endl << "this = " << this
    << std::endl << "needs_to_prepare_for_this_point() = "
    << needs_to_prepare_for_this_point();
    std::cout << std::endl;**/

    if( needs_to_prepare_for_this_point() )
      {

        calculate();
        finish_preparing_for_this_point();

      }

  }



  inline double
  cascade_acceptance_value::get_value()
  {

    if( needs_to_prepare_for_this_point() )
      {

        acceptance_value = 0.0;

        for( std::vector< channel_calculator* >::iterator
             channel_iterator = channels.begin();
             channels.end() > channel_iterator;
             ++channel_iterator )
          {

            // debugging:
            /**std::cout
            << std::endl
            << "debugging: " << *channel_iterator << "->update_numbers()";
            std::cout << std::endl;**/

            (*channel_iterator)->update_numbers();

          }

        finish_preparing_for_this_point();

      }

    return acceptance_value;

  }

  inline double* const
  cascade_acceptance_value::get_value_pointer()
  {

    return &acceptance_value;

  }

  inline void
  cascade_acceptance_value::add_channel( channel_calculator* given_channel )
  {

    channels.push_back( given_channel );

  }



  inline bool
  charge_summed_cascade_calculator::is_requested( double given_primary_cut,
                                                  double given_secondary_cut,
                                                  double given_jet_cut,
                               lepton_acceptance_value* const given_kinematics,
           CppSLHA::particle_property_set const* const given_decaying_scolored,
                                 bool const given_scolored_is_not_antiparticle,
              CppSLHA::particle_property_set const* const given_decaying_EWino,
                                   bool const given_EWino_is_not_antiparticle )
  /* this returns true if the relevant cuts, the kinematics, & the sparticles
   * match those held. relevant cuts are those over 0.0, any cuts of 0.0 or
   * less are assumed to be irrelevant to the signal requesting these values,
   * so this calculator can be given to the signal.
   */
  {

    if( ( ( given_primary_cut == primary_cut )
          ||
          ( given_primary_cut <= 0.0 ) )
        &&
        ( ( given_secondary_cut == secondary_cut )
          ||
          ( given_secondary_cut <= 0.0 ) )
        &&
        ( ( given_jet_cut == jet_cut )
          ||
          ( given_jet_cut <= 0.0 ) )
        &&
        ( given_kinematics == kinematics )
        &&
        ( given_decaying_scolored == decaying_scolored )
        &&
        ( given_scolored_is_not_antiparticle == scolored_is_not_antiparticle )
        &&
        ( given_decaying_EWino == decaying_EWino )
        &&
        ( given_EWino_is_not_antiparticle == EWino_is_not_antiparticle ) )
      {

        return true;

      }
    else
      {

        return false;

      }

  }


  inline double
  charge_summed_cascade_calculator::get_OSSF_minus_OSDF_acceptance()
  {

    return OSSF_minus_OSDF_leptons.get_value();

  }

}

#endif /* LHC_FASER_BASE_ELECTROWEAK_CASCADE_STUFF_HPP_ */
