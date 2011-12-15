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
    acceptanceValues( double const defaultUnsetValues = 0.0 );
    ~acceptanceValues();

    void
    reset( double const unsetValues = 0.0 );
    // this sets all acceptances to unsetValues & alreadyCalculatedFlag to
    // false.
    bool
    needsToBeCalculated();
    double
    getTwoJets()
    const;
    void
    setTwoJets( double const inputValue );
    void
    addToTwoJets( double const inputValue );
    double
    getOneJetOneNegativeElectron()
    const;
    void
    setOneJetOneNegativeElectron( double const inputValue );
    void
    addToOneJetOneNegativeElectron( double const inputValue );
    double
    getOneJetOnePositiveElectron()
    const;
    void
    setOneJetOnePositiveElectron( double const inputValue );
    void
    addToOneJetOnePositiveElectron( double const inputValue );
    double
    getOneJetOneNegativeMuon()
    const;
    void
    setOneJetOneNegativeMuon( double const inputValue );
    void
    addToOneJetOneNegativeMuon( double const inputValue );
    double
    getOneJetOnePositiveMuon()
    const;
    void
    setOneJetOnePositiveMuon( double const inputValue );
    void
    addToOneJetOnePositiveMuon( double const inputValue );
    double
    getOneJetZeroLeptons()
    const;
    void
    setOneJetZeroLeptons( double const inputValue );
    void
    addToOneJetZeroLeptons( double const inputValue );
    double
    getOssfMinusOsdf()
    const;
    void
    setOssfMinusOsdf( double const inputValue );
    void
    addToOssfMinusOsdf( double const inputValue );
    double
    getElectronPlusAntielectron()
    const;
    void
    setElectronPlusAntielectron( double const inputValue );
    void
    addToElectronPlusAntielectron( double const inputValue );
    double
    getNegativeElectronPlusPositiveMuon()
    const;
    void
    setNegativeElectronPlusPositiveMuon( double const inputValue );
    void
    addToNegativeElectronPlusPositiveMuon( double const inputValue );
    double
    getNegativeMuonPlusPositiveElectron()
    const;
    void
    setNegativeMuonPlusPositiveElectron( double const inputValue );
    void
    addToNegativeMuonPlusPositiveElectron( double const inputValue );
    double
    getMuonPlusAntimuon()
    const;
    void
    setMuonPlusAntimuon( double const inputValue );
    void
    addToMuonPlusAntimuon( double const inputValue );
    double
    getZeroJetsOneNegativeElectron()
    const;
    void
    setZeroJetsOneNegativeElectron( double const inputValue );
    void
    addToZeroJetsOneNegativeElectron( double const inputValue );
    double
    getZeroJetsOnePositiveElectron()
    const;
    void
    setZeroJetsOnePositiveElectron( double const inputValue );
    void
    addToZeroJetsOnePositiveElectron( double const inputValue );
    double
    getZeroJetsOneNegativeMuon()
    const;
    void
    setZeroJetsOneNegativeMuon( double const inputValue );
    void
    addToZeroJetsOneNegativeMuon( double const inputValue );
    double
    getZeroJetsOnePositiveMuon()
    const;
    void
    setZeroJetsOnePositiveMuon( double const inputValue );
    void
    addToZeroJetsOnePositiveMuon( double const inputValue );
    double
    getZeroJetsZeroLeptons()
    const;
    void
    setZeroJetsZeroLeptons( double const inputValue );
    void
    addToZeroJetsZeroLeptons( double const inputValue );

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
  class electroweakCascade : public getsReadiedForNewPoint
  {
  public:
    electroweakCascade( leptonAcceptanceParameterSet* const kinematicsTable,
                        effectiveSquarkMassHolder* const effectiveSquarkMass,
                        particlePointer const coloredDecayer,
                        particlePointer const electroweakDecayer,
                        particlePointer const intermediateDecayer,
                        bool const canDoOssfMinusOsdf,
                        inputHandler const* const inputShortcut );
    virtual
    ~electroweakCascade();

    inputHandler const*
    getShortcut()
    const;
    double
    getAcceptance( acceptanceCutSet const* const currentCuts,
                   int const numberOfAdditionalJets,
                   int const numberOfNegativeElectrons,
                   int const numberOfPositiveElectrons,
                   int const numberOfNegativeMuons,
                   int const numberOfPositiveMuons );
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

    double
    getOssfMinusOsdf( acceptanceCutSet* const currentCuts );
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

    static void
    cachePairConstruction( acceptanceValues** const cachedValuesAndFlag,
                           acceptanceCutSet** const cachedKey,
                           acceptanceCutSet const* const constructionKey );
    /* this makes a new acceptanceValues instance to hold an acceptance value
     * along with a flag for whether it has been calculated or not, &
     * associates it with a key.
     */

    static void
    cachePairReset( acceptanceValues* const cachedValuesAndFlag,
                    acceptanceCutSet* const cachedKey,
                    acceptanceCutSet const* const constructionKey );
    /* this resets an acceptanceValues instance to hold an acceptance value
     * along with a flag for whether it has been calculated or not, &
     * associates it with a key.
     */

  protected:
    static double const tauPairToPionPairBr;
    static double const tauToPionTimesTauToElectronBr;
    static double const tauToPionTimesTauToMuonBr;
    static double const tauPairToElectronPairBr;
    static double const tauToElectronTimesTauToMuonBr;
    static double const tauPairToMuonPairBr;
    static int const numberOfIntegrationBins;
    static double const cachedAcceptanceResetValue;

    leptonAcceptanceParameterSet* const kinematicsTable;
    effectiveSquarkMassHolder* const effectiveSquarkMass;
    particlePointer const coloredDecayer;
    //bool const coloredDecayerIsNotAntiparticle;
    particlePointer const electroweakDecayer;
    //bool const electroweakDecayerIsNotAntiparticle;
    particlePointer const intermediateDecayer;
    bool const canDoOssfMinusOsdf;
    exclusiveBrCalculator* firstBr;
    exclusiveBrCalculator* secondBr;
    double cascadeBr;
    inputHandler const* const inputShortcut;
    acceptanceValues* currentAcceptance;
    minimalAllocationKeyedVector< acceptanceValues,
                                  acceptanceCutSet > acceptancesPerCutSet;
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
    calculateAcceptance( acceptanceCutSet const* const currentCuts,
                         acceptanceValues* const currentAcceptance )
    = 0;
    // this should set all the appropriate acceptance * branching ratio values
    // for the cascade.

    double
    integrateAcceptance( leptonEnergyDistribution* const leptonDistribution,
                         double const transverseMomentumCut );
    /* this numerically integrates leptonDistribution with the lepton
     * acceptance histogram from kinematics, scaled to having
     * transverseMomentumCut as its transverse momentum cut.
     */
  };





  // inline functions:


  inline void
  acceptanceValues::reset( double const unsetValues )
  // this sets all acceptances to unsetValues & notAlreadyCalculatedFlag to
  // true.
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

  inline bool
  acceptanceValues::needsToBeCalculated()
  {
    if( notAlreadyCalculatedFlag )
    {
      notAlreadyCalculatedFlag = false;
      return true;
    }
    else
    {
      return false;
    }
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



  inline inputHandler const*
  electroweakCascade::getShortcut()
  const
  {
    return inputShortcut;
  }

  inline double
  electroweakCascade::getOssfMinusOsdf( acceptanceCutSet* const currentCuts )
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
      if( needsToPrepareForThisPoint() )
      {
        acceptancesPerCutSet.clearEntries();
      }
      currentAcceptance = acceptancesPerCutSet.getPointer( currentCuts );
      if( currentAcceptance->needsToBeCalculated() )
      {
        calculateAcceptance( currentCuts,
                             currentAcceptance );
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
                                acceptanceCutSet const* const constructionKey )
  // this makes a new std::pair to hold an acceptance value along with a flag
  // for whether it has been calculated or not, & associates it with a key.
  {
    *cachedValuesAndFlag = new acceptanceValues( cachedAcceptanceResetValue );
    // the new pair has the default "unset" acceptanceValues.
    *cachedKey = new acceptanceCutSet( constructionKey );
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "electroweakCascade::cachePairConstruction( ..., "
    << constructionKey << " ) called.";
    std::cout << std::endl;**/

  }

  inline void
  electroweakCascade::cachePairReset(
                                   acceptanceValues* const cachedValuesAndFlag,
                                      acceptanceCutSet* const cachedKey,
                                acceptanceCutSet const* const constructionKey )
  // this resets a std::pair to hold an acceptance value along with a flag
  // for whether it has been calculated or not, & associates it with a key.
  {
    cachedValuesAndFlag->reset( cachedAcceptanceResetValue );
    // the new pair has the default "unset" acceptanceValues.
    cachedKey->becomeCopyOf( constructionKey );
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "electroweakCascade::cachePairReset( ..., "
    << constructionKey << " ) called.";
    std::cout << std::endl;**/
  }

}

#endif /* LHC_FASER_BASE_ELECTROWEAK_CASCADE_STUFF_HPP_ */
