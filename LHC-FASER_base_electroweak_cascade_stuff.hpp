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

#include "LHC-FASER_lepton_kinematics_stuff.hpp"
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
    void
    reset()
    // this sets all values to CppSLHA::CppSLHA_global::really_wrong_value &
    // alreadyCalculatedFlag to false.
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
    void
    addToTwoJets( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getOneJetOneNegativeElectron()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setOneJetOneNegativeElectron( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    addToOneJetOneNegativeElectron( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getOneJetOnePositiveElectron()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setOneJetOnePositiveElectron( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    addToOneJetOnePositiveElectron( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getOneJetOneNegativeMuon()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setOneJetOneNegativeMuon( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    addToOneJetOneNegativeMuon( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getOneJetOnePositiveMuon()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setOneJetOnePositiveMuon( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    addToOneJetOnePositiveMuon( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getOneJetZeroLeptons()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setOneJetZeroLeptons( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    addToOneJetZeroLeptons( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getOssfMinusOsdf()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setOssfMinusOsdf( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    addToOssfMinusOsdf( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getElectronPlusAntielectron()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setElectronPlusAntielectron( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    addToElectronPlusAntielectron( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getNegativeElectronPlusPositiveMuon()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setNegativeElectronPlusPositiveMuon( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    addToNegativeElectronPlusPositiveMuon( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getNegativeMuonPlusPositiveElectron()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setNegativeMuonPlusPositiveElectron( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    addToNegativeMuonPlusPositiveElectron( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getMuonPlusAntimuon()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setMuonPlusAntimuon( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    addToMuonPlusAntimuon( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getZeroJetsOneNegativeElectron()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setZeroJetsOneNegativeElectron( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    addToZeroJetsOneNegativeElectron( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getZeroJetsOnePositiveElectron()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setZeroJetsOnePositiveElectron( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    addToZeroJetsOnePositiveElectron( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getZeroJetsOneNegativeMuon()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setZeroJetsOneNegativeMuon( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    addToZeroJetsOneNegativeMuon( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getZeroJetsOnePositiveMuon()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setZeroJetsOnePositiveMuon( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    addToZeroJetsOnePositiveMuon( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getZeroJetsZeroLeptons()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setZeroJetsZeroLeptons( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    addToZeroJetsZeroLeptons( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:
    bool notAlreadyCalculatedFlag;
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
    double zeroJetsOnePositiveElectron;
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
    electroweakCascade( leptonAcceptanceParameterSet* const kinematics,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                CppSLHA::particle_property_set const* const electroweakDecayer,
               CppSLHA::particle_property_set const* const intermediateDecayer,
                        bool const canDoOssfMinusOsdf,
                        input_handler const* const shortcut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~electroweakCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    input_handler const*
    getShortcut()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

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
    getOssfMinusOsdf( acceptanceCutSet* const cuts )
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
    static double const tauPairToPionPairBr;
    static double const tauToPionTimesTauToElectronBr;
    static double const tauToPionTimesTauToMuonBr;
    static double const tauPairToElectronPairBr;
    static double const tauToElectronTimesTauToMuonBr;
    static double const tauPairToMuonPairBr;

    leptonAcceptanceParameterSet* const kinematics;
    CppSLHA::particle_property_set const* const coloredDecayer;
    //bool const coloredDecayerIsNotAntiparticle;
    CppSLHA::particle_property_set const* const electroweakDecayer;
    //bool const electroweakDecayerIsNotAntiparticle;
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
    electroweakCascadeSet( leptonAcceptanceParameterSet* const kinematics,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                CppSLHA::particle_property_set const* const electroweakDecayer,
                           input_handler const* const shortcut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~electroweakCascadeSet()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    bool
    isEquivalent( CppSLHA::particle_property_set const* const coloredDecayer,
               CppSLHA::particle_property_set const* const electroweakDecayer )
    // this returns true if the coloredDecayers & electroweakDecayers match.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

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
    getOssfMinusOsdf( acceptanceCutSet* const cuts )
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
    leptonAcceptanceParameterSet* const kinematics;
    input_handler const* const shortcut;
    CppSLHA::particle_property_set const* const coloredDecayer;
    CppSLHA::particle_property_set const* const electroweakDecayer;
    std::vector< electroweakCascade* > cascades;
    electroweakCascade* currentCascade; // this is used for filling cascades.
  };

  // this gives out pointers to electroweakCascadeSets at a fixed beam energy
  // based on the requested colored sparticle & electroweakino or vector boson.
  class electroweakCascadeHandler
  {
  public:
    electroweakCascadeHandler(
                      leptonAcceptancesForOneBeamEnergy* const kinematicsTable,
                               int const beamEnergy,
                               input_handler const* const shortcut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    virtual
    ~electroweakCascadeHandler()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    electroweakCascadeSet*
    getCascadeSet( CppSLHA::particle_property_set const* const coloredDecayer,
               CppSLHA::particle_property_set const* const electroweakDecayer )
    /* this looks to see if it already has an electroweakCascadeSet
     * corresponding to the requested pairing, & if it does, it returns a
     * pointer to it, & if it doesn't, it constructs a new
     * electroweakCascadeSet & returns a pointer to that.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    leptonAcceptancesForOneBeamEnergy* const kinematicsTable;
    int const beamEnergy;
    input_handler const* const shortcut;
    std::vector< electroweakCascadeSet* > cascadeSets;
    electroweakCascadeSet* currentCascadeSet;
    // this is used for filling cascadeSets.
  };





  // inline functions:

  inline void
  acceptanceValues::reset( double const unsetValues )
  // this sets all values to unsetValues & notAlreadyCalculatedFlag to true.
  {
    notAlreadyCalculatedFlag = true;
    twoJets = unsetValues;
    oneJetOneNegativeElectron = unsetValues;
    oneJetOnePositiveElectron = unsetValues;
    oneJetOneNegativeMuon = unsetValues;
    oneJetOnePositiveMuon = unsetValues;
    oneJetZeroLeptons = unsetValues;
    ossfMinusOsdf = unsetValues;
    electronPlusAntielectron = unsetValues;
    negativeElectronPlusPositiveMuon = unsetValues;
    negativeMuonPlusPositiveElectron = unsetValues;
    muonPlusAntimuon = unsetValues;
    zeroJetsOneNegativeElectron = unsetValues;
    zeroJetsOnePositiveElectron = unsetValues;
    zeroJetsOneNegativeMuon = unsetValues;
    zeroJetsOnePositiveMuon = unsetValues;
    zeroJetsZeroLeptons = unsetValues;
  }

  inline void
  acceptanceValues::reset()
  // this sets all values to CppSLHA::CppSLHA_global::really_wrong_value &
  // alreadyCalculatedFlag to false.
  {
    reset( CppSLHA::CppSLHA_global::really_wrong_value );
  }

  inline bool
  acceptanceValues::notAlreadyCalculated()
  {
    return notAlreadyCalculatedFlag;
  }

  inline void
  acceptanceValues::flagAsAlreadyCalculated()
  {
    notAlreadyCalculatedFlag = false;
  }

  inline double
  acceptanceValues::getTwoJets()
  const
  {
    return twoJets;
  }

  inline void
  acceptanceValues::setTwoJets( double const inputValue )
  {
    twoJets = inputValue;
  }

  inline void
  acceptanceValues::addToTwoJets( double const inputValue )
  {
    twoJets += inputValue;
  }

  inline double
  acceptanceValues::getOneJetOneNegativeElectron()
  const
  {
    return oneJetOneNegativeElectron;
  }

  inline void
  acceptanceValues::setOneJetOneNegativeElectron( double const inputValue )
  {
    oneJetOneNegativeElectron = inputValue;
  }

  inline void
  acceptanceValues::addToOneJetOneNegativeElectron( double const inputValue )
  {
    oneJetOneNegativeElectron += inputValue;
  }

  inline double
  acceptanceValues::getOneJetOnePositiveElectron()
  const
  {
    return oneJetOnePositiveElectron;
  }

  inline void
  acceptanceValues::setOneJetOnePositiveElectron( double const inputValue )
  {
    oneJetOnePositiveElectron = inputValue;
  }

  inline void
  acceptanceValues::addToOneJetOnePositiveElectron( double const inputValue )
  {
    oneJetOnePositiveElectron += inputValue;
  }

  inline double
  acceptanceValues::getOneJetOneNegativeMuon()
  const
  {
    return oneJetOneNegativeMuon;
  }

  inline void
  acceptanceValues::setOneJetOneNegativeMuon( double const inputValue )
  {
    oneJetOneNegativeMuon = inputValue;
  }

  inline void
  acceptanceValues::addToOneJetOneNegativeMuon( double const inputValue )
  {
    oneJetOneNegativeMuon += inputValue;
  }

  inline double
  acceptanceValues::getOneJetOnePositiveMuon()
  const
  {
    return oneJetOnePositiveMuon;
  }

  inline void
  acceptanceValues::setOneJetOnePositiveMuon( double const inputValue )
  {
    oneJetOnePositiveMuon = inputValue;
  }

  inline void
  acceptanceValues::addToOneJetOnePositiveMuon( double const inputValue )
  {
    oneJetOnePositiveMuon += inputValue;
  }

  inline double
  acceptanceValues::getOneJetZeroLeptons()
  const
  {
    return oneJetZeroLeptons;
  }

  inline void
  acceptanceValues::setOneJetZeroLeptons( double const inputValue )
  {
    oneJetZeroLeptons = inputValue;
  }

  inline void
  acceptanceValues::addToOneJetZeroLeptons( double const inputValue )
  {
    oneJetZeroLeptons += inputValue;
  }

  inline double
  acceptanceValues::getOssfMinusOsdf()
  const
  {
    return ossfMinusOsdf;
  }

  inline void
  acceptanceValues::setOssfMinusOsdf( double const inputValue )
  {
    ossfMinusOsdf = inputValue;
  }

  inline void
  acceptanceValues::addToOssfMinusOsdf( double const inputValue )
  {
    ossfMinusOsdf += inputValue;
  }

  inline double
  acceptanceValues::getElectronPlusAntielectron()
  const
  {
    return electronPlusAntielectron;
  }

  inline void
  acceptanceValues::setElectronPlusAntielectron( double const inputValue )
  {
    electronPlusAntielectron = inputValue;
  }

  inline void
  acceptanceValues::addToElectronPlusAntielectron( double const inputValue )
  {
    electronPlusAntielectron += inputValue;
  }

  inline double
  acceptanceValues::getNegativeElectronPlusPositiveMuon()
  const
  {
    return negativeElectronPlusPositiveMuon;
  }

  inline void
  acceptanceValues::setNegativeElectronPlusPositiveMuon(
                                                      double const inputValue )
  {
    negativeElectronPlusPositiveMuon = inputValue;
  }

  inline void
  acceptanceValues::addToNegativeElectronPlusPositiveMuon(
                                                      double const inputValue )
  {
    negativeElectronPlusPositiveMuon += inputValue;
  }

  inline double
  acceptanceValues::getNegativeMuonPlusPositiveElectron()
  const
  {
    return negativeMuonPlusPositiveElectron;
  }

  inline void
  acceptanceValues::setNegativeMuonPlusPositiveElectron(
                                                      double const inputValue )
  {
    negativeMuonPlusPositiveElectron = inputValue;
  }

  inline void
  acceptanceValues::addToNegativeMuonPlusPositiveElectron(
                                                      double const inputValue )
  {
    negativeMuonPlusPositiveElectron += inputValue;
  }

  inline double
  acceptanceValues::getMuonPlusAntimuon()
  const
  {
    return muonPlusAntimuon;
  }

  inline void
  acceptanceValues::setMuonPlusAntimuon( double const inputValue )
  {
    muonPlusAntimuon = inputValue;
  }

  inline void
  acceptanceValues::addToMuonPlusAntimuon( double const inputValue )
  {
    muonPlusAntimuon += inputValue;
  }

  inline double
  acceptanceValues::getZeroJetsOneNegativeElectron()
  const
  {
    return zeroJetsOneNegativeElectron;
  }

  inline void
  acceptanceValues::setZeroJetsOneNegativeElectron( double const inputValue )
  {
    zeroJetsOneNegativeElectron = inputValue;
  }

  inline void
  acceptanceValues::addToZeroJetsOneNegativeElectron( double const inputValue )
  {
    zeroJetsOneNegativeElectron += inputValue;
  }

  inline double
  acceptanceValues::getZeroJetsOnePositiveElectron()
  const
  {
    return zeroJetsOnePositiveElectron;
  }

  inline void
  acceptanceValues::setZeroJetsOnePositiveElectron( double const inputValue )
  {
    zeroJetsOnePositiveElectron = inputValue;
  }

  inline void
  acceptanceValues::addToZeroJetsOnePositiveElectron( double const inputValue )
  {
    zeroJetsOnePositiveElectron += inputValue;
  }

  inline double
  acceptanceValues::getZeroJetsOneNegativeMuon()
  const
  {
    return zeroJetsOneNegativeMuon;
  }

  inline void
  acceptanceValues::setZeroJetsOneNegativeMuon( double const inputValue )
  {
    zeroJetsOneNegativeMuon = inputValue;
  }

  inline void
  acceptanceValues::addToZeroJetsOneNegativeMuon( double const inputValue )
  {
    zeroJetsOneNegativeMuon += inputValue;
  }

  inline double
  acceptanceValues::getZeroJetsOnePositiveMuon()
  const
  {
    return zeroJetsOnePositiveMuon;
  }

  inline void
  acceptanceValues::setZeroJetsOnePositiveMuon( double const inputValue )
  {
    zeroJetsOnePositiveMuon = inputValue;
  }

  inline void
  acceptanceValues::addToZeroJetsOnePositiveMuon( double const inputValue )
  {
    zeroJetsOnePositiveMuon += inputValue;
  }

  inline double
  acceptanceValues::getZeroJetsZeroLeptons()
  const
  {
    return zeroJetsZeroLeptons;
  }

  inline void
  acceptanceValues::setZeroJetsZeroLeptons( double const inputValue )
  {
    zeroJetsZeroLeptons = inputValue;
  }

  inline void
  acceptanceValues::addToZeroJetsZeroLeptons( double const inputValue )
  {
    zeroJetsZeroLeptons += inputValue;
  }



  inline input_handler const*
  electroweakCascade::getShortcut()
  const
  {
    return shortcut;
  }

  inline double
  electroweakCascade::getOssfMinusOsdf( acceptanceCutSet* const cuts )
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



  inline bool
  electroweakCascadeSet::isEquivalent(
                    CppSLHA::particle_property_set const* const coloredDecayer,
               CppSLHA::particle_property_set const* const electroweakDecayer )
  // this returns true if the coloredDecayers & electroweakDecayers match.
  {
    if( ( coloredDecayer == this->coloredDecayer )
        &&
        ( electroweakDecayer == this->electroweakDecayer ) )
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  inline double
  electroweakCascadeSet::getAcceptance( acceptanceCutSet const* const cuts,
                                        int const numberOfAdditionalJets,
                                        int const numberOfNegativeElectrons,
                                        int const numberOfPositiveElectrons,
                                        int const numberOfNegativeMuons,
                                        int const numberOfPositiveMuons )
  // this adds up the acceptances from the stored electroweakCascades.
  {
    double returnValue( 0.0 );
    for( std::vector< electroweakCascade* >::iterator
         cascadeIterator = cascades.begin();
         cascades.end() > cascadeIterator;
         ++cascadeIterator )
    {
      returnValue += (*cascadeIterator)->getAcceptance( cuts,
                                                        numberOfAdditionalJets,
                                                     numberOfNegativeElectrons,
                                                     numberOfPositiveElectrons,
                                                        numberOfNegativeMuons,
                                                       numberOfPositiveMuons );
    }
    return returnValue;
  }
}

#endif /* LHC_FASER_BASE_ELECTROWEAK_CASCADE_STUFF_HPP_ */
