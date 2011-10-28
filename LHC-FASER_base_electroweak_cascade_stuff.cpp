/*
 * LHC-FASER_cascade_stuff.cpp
 *
 *  Created on: 08 Mar 2011
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

#include "LHC-FASER_base_electroweak_cascade_stuff.hpp"

namespace LHC_FASER
{
  acceptanceValues::acceptanceValues( double const defaultUnsetValues ) :
      notAlreadyCalculatedFlag( true ),
      twoJets( defaultUnsetValues ),
      oneJetOneNegativeElectron( defaultUnsetValues ),
      oneJetOnePositiveElectron( defaultUnsetValues ),
      oneJetOneNegativeMuon( defaultUnsetValues ),
      oneJetOnePositiveMuon( defaultUnsetValues ),
      oneJetZeroLeptons( defaultUnsetValues ),
      ossfMinusOsdf( defaultUnsetValues ),
      electronPlusAntielectron( defaultUnsetValues ),
      negativeElectronPlusPositiveMuon( defaultUnsetValues ),
      negativeMuonPlusPositiveElectron( defaultUnsetValues ),
      muonPlusAntimuon( defaultUnsetValues ),
      zeroJetsOneNegativeElectron( defaultUnsetValues ),
      zeroJetsOnePositiveElectron( defaultUnsetValues ),
      zeroJetsOneNegativeMuon( defaultUnsetValues ),
      zeroJetsOnePositiveMuon( defaultUnsetValues ),
      zeroJetsZeroLeptons( defaultUnsetValues )
  {
    // just an initialization list.
  }

  acceptanceValues::~acceptanceValues()
  {
    // does nothing.
  }


  double const
  electroweakCascade::tauPairToPionPairBr(
                          ( CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
                     * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR ) );
  double const
  electroweakCascade::tauToPionTimesTauToElectronBr(
                          ( CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
                  * CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR ) );
  double const
  electroweakCascade::tauToPionTimesTauToMuonBr(
                          ( CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
                      * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR ) );
  double const
  electroweakCascade::tauPairToElectronPairBr(
                       ( CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR
                  * CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR ) );
  double const
  electroweakCascade::tauToElectronTimesTauToMuonBr(
                       ( CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR
                      * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR ) );
  double const
  electroweakCascade::tauPairToMuonPairBr(
                           ( CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR
                      * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR ) );
  int const
  electroweakCascade::numberOfIntegrationBins( 10 );
  double const
  electroweakCascade::cachedAcceptanceResetValue( 0.0 );

  electroweakCascade::electroweakCascade(
                                leptonAcceptanceParameterSet* const kinematics,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                                          particlePointer const coloredDecayer,
                                  //bool const coloredDecayerIsNotAntiparticle,
                                      particlePointer const electroweakDecayer,
                              //bool const electroweakDecayerIsNotAntiparticle,
                                     particlePointer const intermediateDecayer,
                                          bool const canDoOssfMinusOsdf,
                                    inputHandler const* const inputShortcut ) :
      getsReadiedForNewPoint( inputShortcut->getReadier() ),
      kinematics( kinematics ),
      effectiveSquarkMass( effectiveSquarkMass ),
      coloredDecayer( coloredDecayer ),
      //coloredDecayerIsNotAntiparticle( coloredDecayerIsNotAntiparticle ),
      electroweakDecayer( electroweakDecayer ),
      intermediateDecayer( intermediateDecayer ),
      canDoOssfMinusOsdf( canDoOssfMinusOsdf ),
      firstBr( NULL ),
      secondBr( NULL ),
      cascadeBr( CppSLHA::CppSLHA_global::really_wrong_value ),
      inputShortcut( inputShortcut ),
      currentAcceptance( NULL ),
      acceptancesPerCutSet( &acceptanceCutSet::compareJetAndBothLeptonCuts,
                            &electroweakCascade::cachePairConstruction,
                            &electroweakCascade::cachePairReset ),
      activeDistributions()
  {
    // just an initialization list.
  }

  electroweakCascade::~electroweakCascade()
  {
    for( std::vector< leptonEnergyDistribution* >::iterator
         deletionIterator( activeDistributions.begin() );
         activeDistributions.end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
  }


  double
  electroweakCascade::getAcceptance( acceptanceCutSet const* const currentCuts,
                                     int const numberOfAdditionalJets,
                                     int const numberOfNegativeElectrons,
                                     int const numberOfPositiveElectrons,
                                     int const numberOfNegativeMuons,
                                     int const numberOfPositiveMuons )
  /* this does the common job of checking to see if the point has been
   * updated before calling the relevant protected virtual function, which
   * returns the set of acceptanceValues for acceptance * branching ratio FOR THE
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
  {
    if( needsToPrepareForThisPoint() )
    {
      acceptancesPerCutSet.clearEntries();
    }
    if( validSignal( numberOfAdditionalJets,
                     numberOfNegativeElectrons,
                     numberOfPositiveElectrons,
                     numberOfNegativeMuons,
                     numberOfPositiveMuons ) )
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "electroweakCascade::getAcceptance( ... ) about to look for pointer"
      << " from acceptancesPerCutSet, [ "
      << *(coloredDecayer->get_name()) << ", "
      << *(electroweakDecayer->get_name()) << " ].";
      std::cout << std::endl;**/

      currentAcceptance = acceptancesPerCutSet.getPointer( currentCuts );
      // at this point, currentAcceptance is keyed to a copy of cuts.
      if( currentAcceptance->needsToBeCalculated() )
      {
        calculateAcceptance( currentCuts,
                             currentAcceptance );
        // now currentAcceptance has values corresponding to cuts as well as
        // being keyed to a copy of cuts.
      }
      /* now currentAcceptance has the values which are calculated by this
       * electroweakCascade in the valid entries. the entries in invalid
       * elements have leftover junk, but the condition of validSignal means
       * that they won't be looked at by this code.
       */
      if( 2 == numberOfAdditionalJets )
        // assuming that 2 jets + more than 0 leptons won't pass validSignal...
      {
        // we return the only 2-jet possibility:
        return currentAcceptance->getTwoJets();
      }
      else if( 1 == numberOfAdditionalJets )
        // assuming that 1 jet + more than 1 lepton won't pass validSignal...
      {
        // we return the appropriate 1-jet possibility:
        if( 1 == numberOfNegativeElectrons )
        {
          return currentAcceptance->getOneJetOneNegativeElectron();
        }
        else if( 1 == numberOfPositiveElectrons )
        {
          return currentAcceptance->getOneJetOnePositiveElectron();
        }
        else if( 1 == numberOfNegativeMuons )
        {
          return currentAcceptance->getOneJetOneNegativeMuon();
        }
        else if( 1 == numberOfPositiveMuons )
        {
          return currentAcceptance->getOneJetOnePositiveMuon();
        }
        else
        {
          return currentAcceptance->getOneJetZeroLeptons();
        }
      }
      // now we check the 2-lepton possibilities:
      else if( ( 1 == numberOfNegativeElectrons )
               &&
               ( 1 == numberOfPositiveElectrons ) )
      {
        return currentAcceptance->getElectronPlusAntielectron();
      }
      else if( ( 1 == numberOfNegativeElectrons )
               &&
               ( 1 == numberOfPositiveMuons ) )
      {
        return currentAcceptance->getNegativeElectronPlusPositiveMuon();
      }
      else if( ( 1 == numberOfNegativeMuons )
               &&
               ( 1 == numberOfPositiveElectrons ) )
      {
        return currentAcceptance->getNegativeMuonPlusPositiveElectron();
      }
      else if( ( 1 == numberOfNegativeMuons )
               &&
               ( 1 == numberOfPositiveMuons ) )
      {
        return currentAcceptance->getMuonPlusAntimuon();
      }
      /* now we check the remaining possibilities (assuming that validSignal
       * prevents numberOfAdditionalJets > 2, & any possibility with more than
       * 1 of any charge-&-flavor of lepton, & any possibility with more than 2
       * leptons in total, & any possibility with 2 leptons of the same sign):
       * hence now it is assumed that numberOfAdditionalJets == 0 & any number
       * of leptons besides that in the conditional is zero.
       */
      else if( 1 == numberOfNegativeElectrons )
      {
        return currentAcceptance->getZeroJetsOneNegativeElectron();
      }
      else if( 1 == numberOfPositiveElectrons )
      {
        return currentAcceptance->getZeroJetsOnePositiveElectron();
      }
      else if( 1 == numberOfNegativeMuons )
      {
        return currentAcceptance->getZeroJetsOneNegativeMuon();
      }
      else if( 1 == numberOfPositiveMuons )
      {
        return currentAcceptance->getZeroJetsOnePositiveMuon();
      }
      else
      {
        return currentAcceptance->getZeroJetsZeroLeptons();
      }
    }
    else
    {
      return 0.0;
    }
  }

  double
  electroweakCascade::integrateAcceptance(
                            leptonEnergyDistribution* const leptonDistribution,
                                           double const transverseMomentumCut )
  /* this numerically integrates leptonDistribution with the lepton acceptance
   * histogram from kinematics, scaled to having transverseMomentumCut as its
   * transverse momentum cut.
   */
  {
    if( 0.0 < transverseMomentumCut )
    {
      double binSize( ( ( leptonDistribution->getMaximumEnergy()
                          - leptonDistribution->getMinimumEnergy() )
                        / (double)numberOfIntegrationBins ) );
      double binEnergy( leptonDistribution->getMinimumEnergy() );

      /* I could do it so that it works out the area of each trapezium & sums
       * them, but it's equivalent to taking half a rectangle centered on the
       * starting value with height of the function at the starting value &
       * adding full rectangles for each bin value except the last, which
       * also is a half rectangle.
       */
      double
      returnValue( ( 0.5 * binSize
                         * kinematics->acceptanceAt( binEnergy,
                                                     transverseMomentumCut )
                         * leptonDistribution->valueAt( binEnergy ) ) );
      for( int binCounter( 1 );
           numberOfIntegrationBins > binCounter;
           ++binCounter )
      {
        binEnergy += binSize;
        returnValue += ( binSize
                         * kinematics->acceptanceAt( binEnergy,
                                                     transverseMomentumCut )
                         * leptonDistribution->valueAt( binEnergy ) );
      }
      binEnergy += binSize;
      returnValue += ( 0.5 * binSize
                           * kinematics->acceptanceAt( binEnergy,
                                                       transverseMomentumCut )
                           * leptonDistribution->valueAt( binEnergy ) );
      return returnValue;
    }
    else
      // a cut of 0.0 GeV or a negative value is interpretted as no cut at all,
      // so 1.0 is returned as the acceptance.
    {
      return 1.0;
    }
  }

}  // end of LHC_FASER namespace.
