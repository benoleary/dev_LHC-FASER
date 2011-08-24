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

  electroweakCascade::electroweakCascade(
                                leptonAcceptanceParameterSet* const kinematics,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                                  //bool const coloredDecayerIsNotAntiparticle,
                CppSLHA::particle_property_set const* const electroweakDecayer,
                              //bool const electroweakDecayerIsNotAntiparticle,
               CppSLHA::particle_property_set const* const intermediateDecayer,
                                          bool const canDoOssfMinusOsdf,
                                        inputHandler const* const shortcut ) :
    kinematics( kinematics ),
    effectiveSquarkMass( effectiveSquarkMass ),
    coloredDecayer( coloredDecayer ),
    //coloredDecayerIsNotAntiparticle( coloredDecayerIsNotAntiparticle ),
    electroweakDecayer( electroweakDecayer ),
    intermediateDecayer( intermediateDecayer ),
    firstBr( NULL ),
    secondBr( NULL ),
    cascadeBr( CppSLHA::CppSLHA_global::really_wrong_value ),
    shortcut( shortcut ),
    acceptances( &acceptanceCutSet::compareJetAndBothLeptonCuts,
                 &electroweakCascade::cachePairConstruction,
                 &electroweakCascade::cachePairReset )
  {
    // just an initialization list.
  }

  virtual
  electroweakCascade::~electroweakCascade()
  {
    for( std::vector< leptonEnergyDistribution* >::iterator
         deletionIterator = activeDistributions.begin();
         activeDistributions.end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
  }

  double
  electroweakCascade::getAcceptance( acceptanceCutSet* const cuts,
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
  {
    if( needsToPrepareForThisPoint() )
    {
      acceptances.clear();
    }
    if( validSignal( numberOfAdditionalJets,
                     numberOfNegativeElectrons,
                     numberOfPositiveElectrons,
                     numberOfNegativeMuons,
                     numberOfPositiveMuons ) )
    {
      currentAcceptance = acceptances.getPointer( cuts );
      // at this point, currentAcceptance is keyed to a copy of cuts.
      if( currentAcceptance->notAlreadyCalculated() )
      {
        calculateAcceptance( cuts,
                             currentAcceptance );
        // now currentAcceptance has values corresponding to cuts as well as
        // being keyed to a copy of cuts.
        currentAcceptance->flagAsAlreadyCalculated();
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
      double binEnergy = leptonDistribution->getMinimumEnergy();

      /* I could do it so that it works out the area of each trapezium & sums
       * them, but it's equivalent to taking half a rectangle centered on the
       * starting value with height of the function at the starting value &
       * adding full rectangles for each bin value except the last, which
       * also is a half rectangle.
       */
      double
      returnValue = ( 0.5 * binSize
                          * kinematics->acceptanceAt( binEnergy,
                                                      transverseMomentumCut )
                          * leptonDistribution->valueAt( binEnergy ) );
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
                           * lepton_distribution->value_at( binEnergy ) );
      return returnValue;
    }
    else
      // a cut of 0.0 GeV or a negative value is interpretted as no cut at all,
      // so 1.0 is returned as the acceptance.
    {
      return 1.0;
    }
  }



  electroweakCascadeSet::electroweakCascadeSet(
                                leptonAcceptanceParameterSet* const kinematics,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                CppSLHA::particle_property_set const* const electroweakDecayer,
                                        inputHandler const* const shortcut ) :
    kinematics( kinematics ),
    coloredDecayer( coloredDecayer ),
    electroweakDecayer( electroweakDecayer ),
    shortcut( shortcut )
  {
    effectiveSquarkMassHolder* effectiveSquarkMassPointer;
    if( shortcut->getGluino() == coloredDecayer )
    {
      effectiveSquarkMassPointer = kinematics;
    }
    else
    {
      effectiveSquarkMassPointer
      = shortcut->getOnShellEffectiveSquarkMass( coloredDecayer );
    }

    if( CppSLHA::PDG_code::neutralino_one
        == electroweakDecayer->get_PDG_code() )
    {
      currentCascade = new lightestNeutralinoCascade();
      cascades.push_back( currentCascade );
    }
    else if( shortcut->isIn( electroweakDecayer->get_PDG_code(),
                              shortcut->getUnstableNeutralinos() ) )
    {
      currentCascade = new neutralinoToSemuCascade( kinematics,
                                                    effectiveSquarkMassPointer,
                                                    coloredDecayer,
                                                    electroweakDecayer,
                                                   shortcut->getSelectronL(),
                                                    shortcut );
      cascades.push_back( currentCascade );
      currentCascade = new neutralinoToSemuCascade( kinematics,
                                                    effectiveSquarkMassPointer,
                                                    coloredDecayer,
                                                    electroweakDecayer,
                                                   shortcut->getSelectronR(),
                                                    shortcut );
      cascades.push_back( currentCascade );
      currentCascade = new neutralinoToSemuCascade( kinematics,
                                                    effectiveSquarkMassPointer,
                                                    coloredDecayer,
                                                    electroweakDecayer,
                                                    shortcut->getSmuonL(),
                                                    shortcut );
      cascades.push_back( currentCascade );
      currentCascade = new neutralinoToSemuCascade( kinematics,
                                                    effectiveSquarkMassPointer,
                                                    coloredDecayer,
                                                    electroweakDecayer,
                                                    shortcut->getSmuonR(),
                                                    shortcut );
      cascades.push_back( currentCascade );
      currentCascade = new neutralinoToStauCascade( kinematics,
                                                    effectiveSquarkMassPointer,
                                                    coloredDecayer,
                                                    electroweakDecayer,
                                                    shortcut->getStauOne(),
                                                    shortcut );
      cascades.push_back( currentCascade );
      currentCascade = new neutralinoToStauCascade( kinematics,
                                                    effectiveSquarkMassPointer,
                                                    coloredDecayer,
                                                    electroweakDecayer,
                                                    shortcut->getStauTwo(),
                                                    shortcut );
      cascades.push_back( currentCascade );
      currentCascade = new neutralinoToZCascade( kinematics,
                                                 effectiveSquarkMassPointer,
                                                 coloredDecayer,
                                                 electroweakDecayer,
                                                 shortcut );
      cascades.push_back( currentCascade );
      for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
           vevedBosonIterator
           = shortcut->getNeutralEwsbScalarsAndPseudoscalars()->begin();
           shortcut->getNeutralEwsbScalarsAndPseudoscalars()->end()
           > vevedBosonIterator;
           ++vevedBosonIterator )
      {
        currentCascade = new neutralinoToHiggsCascade( kinematics,
                                                    effectiveSquarkMassPointer,
                                                       coloredDecayer,
                                                       electroweakDecayer,
                                                       *vevedBosonIterator,
                                                       shortcut );
        cascades.push_back( currentCascade );
      }
      currentCascade = new neutralinoThreeBodyDecayCascade( kinematics,
                                                    effectiveSquarkMassPointer,
                                                            coloredDecayer,
                                                            electroweakDecayer,
                                                            shortcut );
      cascades.push_back( currentCascade );
      currentCascade = new neutralinoThreeBodyToTausCascade( kinematics,
                                                    effectiveSquarkMassPointer,
                                                             coloredDecayer,
                                                            electroweakDecayer,
                                                             shortcut );
      cascades.push_back( currentCascade );
    }
    else if( shortcut->isIn( electroweakDecayer->get_PDG_code(),
                              shortcut->getCharginos() ) )
    {
      currentCascade = new charginoToSemuCascade( kinematics,
                                                  effectiveSquarkMassPointer,
                                                  coloredDecayer,
                                               coloredDecayerIsNotAntiparticle,
                                                  electroweakDecayer,
                                           electroweakDecayerIsNotAntiparticle,
                                                  shortcut->getSelectronL(),
                                                  shortcut );
      cascades.push_back( currentCascade );
      currentCascade = new charginoToSemuCascade( kinematics,
                                                  effectiveSquarkMassPointer,
                                                  coloredDecayer,
                                                  electroweakDecayer,
                                                  shortcut->getSelectronR(),
                                                  shortcut );
      cascades.push_back( currentCascade );
      currentCascade = new charginoToEmuSneutrinoCascade( kinematics,
                                                    effectiveSquarkMassPointer,
                                                          coloredDecayer,
                                                          electroweakDecayer,
                                          shortcut->getElectronSneutrinoL(),
                                                          shortcut );
      cascades.push_back( currentCascade );
      cascades.push_back( currentCascade );
      currentCascade = new neutralinoToStauCascade( kinematics,
                                                    effectiveSquarkMassPointer,
                                                    coloredDecayer,
                                                    electroweakDecayer,
                                                    shortcut->getStauOne(),
                                                    shortcut );
      cascades.push_back( currentCascade );
      currentCascade = new neutralinoToStauCascade( kinematics,
                                                    effectiveSquarkMassPointer,
                                                    coloredDecayer,
                                                    electroweakDecayer,
                                                    shortcut->getStauTwo(),
                                                    shortcut );
      cascades.push_back( currentCascade );
      currentCascade = new neutralinoToZCascade( kinematics,
                                                 effectiveSquarkMassPointer,
                                                 coloredDecayer,
                                                 electroweakDecayer,
                                                 shortcut );
      cascades.push_back( currentCascade );
      for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
           vevedBosonIterator
           = shortcut->getNeutralEwsbScalarsAndPseudoscalars()->begin();
           shortcut->getNeutralEwsbScalarsAndPseudoscalars()->end()
           > vevedBosonIterator;
           ++vevedBosonIterator )
      {
        currentCascade = new neutralinoToHiggsCascade( kinematics,
                                                    effectiveSquarkMassPointer,
                                                       coloredDecayer,
                                                       electroweakDecayer,
                                                       *vevedBosonIterator,
                                                       shortcut );
        cascades.push_back( currentCascade );
      }
      currentCascade = new neutralinoThreeBodyDecayCascade( kinematics,
                                                    effectiveSquarkMassPointer,
                                                            coloredDecayer,
                                                            electroweakDecayer,
                                                            shortcut );
      cascades.push_back( currentCascade );
      currentCascade = new neutralinoThreeBodyToTausCascade( kinematics,
                                                    effectiveSquarkMassPointer,
                                                             coloredDecayer,
                                                            electroweakDecayer,
                                                             shortcut );
      cascades.push_back( currentCascade );
    }
    else if( CppSLHA::PDG_code::W_plus == electroweakDecayer->get_PDG_code() )
    {
      currentCascade = new wPlusCascade();
      cascades.push_back( currentCascade );
    }
    else if( CppSLHA::PDG_code::Z == electroweakDecayer->get_PDG_code() )
    {
      currentCascade = new zCascade();
      cascades.push_back( currentCascade );
    }
  }

  virtual
  electroweakCascadeSet::~electroweakCascadeSet()
  {
    for( std::vector< electroweakCascade* >::iterator
         deletionIterator = cascades.begin();
         cascades.end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
  }



  electroweakCascadeHandler::electroweakCascadeHandler(
                      leptonAcceptancesForOneBeamEnergy* const kinematicsTable,
                                                        int const beamEnergy,
                                        inputHandler const* const shortcut ) :
    kinematicsTable( kinematicsTable ),
    beamEnergy( beamEnergy ),
    shortcut( shortcut )
  {
    // just an initialization list.
  }

  electroweakCascadeHandler::~electroweakCascadeHandler()
  {
    for( std::vector< electroweakCascadeSet* >::iterator
         deletionIterator = cascadeSets.begin();
         cascadeSets.end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
  }

  electroweakCascadeSet*
  electroweakCascadeHandler::getCascadeSet(
                    CppSLHA::particle_property_set const* const coloredDecayer,
               CppSLHA::particle_property_set const* const electroweakDecayer )
  /* this looks to see if it already has an electroweakCascadeSet
   * corresponding to the requested pairing, & if it does, it returns a
   * pointer to it, & if it doesn't, it constructs a new
   * electroweakCascadeSet & returns a pointer to that.
   */
  {
    electroweakCascadeSet* returnPointer( NULL );
    for( std::vector< electroweakCascadeSet* >::iterator
         cascadeIterator = cascadeSets.begin();
         cascadeSets.end() > cascadeIterator;
         ++cascadeIterator )
    {
      if( (*cascadeIterator)->isEquivalent( coloredDecayer,
                                            electroweakDecayer ) )
      {
        returnPointer = *cascadeIterator;
        cascadeIterator = cascadeSets.end();
      }
    }
    if( NULL == returnPointer )
    {
      returnPointer
      = new electroweakCascadeSet( kinematicsTable->getParameterSets(
                                             coloredDecayer )->getParameterSet(
                                                          electroweakDecayer ),
                                   coloredDecayer,
                                   electroweakDecayer,
                                   shortcut );
      cascadeSets.push_back( returnPointer );
    }
    return returnPointer;
  }








  double
  channel_calculator::integrate_acceptance( double const given_cut,
                   leptonEnergyDistribution const* const lepton_distribution )
  {

    if( given_cut > 0.0 )
      {

        double bin_size = ( ( lepton_distribution->getMaximumEnergy()
                              - lepton_distribution->getMinimumEnergy() )
                            / integration_bins );

        double bin_energy = lepton_distribution->getMinimumEnergy();

        /* I could do it so that it works out the area of each trapezium & sums
         * them, but it's equivalent to taking half a rectangle centered on the
         * starting value with height of the function at the starting value &
         * adding full rectangles for each bin value except the last, which
         * also is a half rectangle.
         */

        double
        return_value
        = ( 0.5 * bin_size
                * kinematics->get_acceptance()->calculateAcceptanceAt( bin_energy,
                                                               given_cut )
                * lepton_distribution->valueAt( bin_energy ) );

        for( int bin_counter = 1;
             integration_bins > bin_counter;
             ++bin_counter )
          {

            bin_energy += bin_size;
            return_value
            += ( bin_size
                 * kinematics->get_acceptance()->calculateAcceptanceAt( bin_energy,
                                                                given_cut )
                 * lepton_distribution->valueAt( bin_energy ) );

          }

        bin_energy += bin_size;
        return_value
        += ( 0.5 * bin_size
                 * kinematics->get_acceptance()->calculateAcceptanceAt( bin_energy,
                                                                given_cut )
                 * lepton_distribution->valueAt( bin_energy ) );

        return return_value;

      }
    else
      // a cut of 0.0 GeV or a negative value is interpretted as no cut at all,
      // so 1.0 is returned as the acceptance.
      {

        return 1.0;

      }

  }


  void
  neutralino_to_Higgs::resetValues()
  {

    // we start by checking the BR of the channel to the EWSB scalar:
    double
    to_Higgs_BR
    = decaying_EWino->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                            mediating_particle->get_PDG_code(),
                                           CppSLHA::PDG_code::neutralino_one );

    if( lhcFaserGlobal::negligibleBr < to_Higgs_BR )
      // if the branching ratio into this channel is not negligible...
      {

        // we start with the direct distributions:
        Higgs_muon_plus_antimuon direct_lepton( shortcut->getCppSlha(),
                                                decaying_scolored,
                    leptonKinematics->getAcceptance()->getEffectiveSquarkMass(),
                                                decaying_EWino,
                                                mediating_particle,
                                              shortcut->getNeutralinoOne() );
        current_distribution = &direct_lepton;
        double direct_muon_pass = integrate_acceptance( primary_cut,
                                                        current_distribution );
        double
        direct_muon_fail = ( 1.0 - integrate_acceptance( secondary_cut,
                                                      current_distribution ) );
        double direct_jet_pass = integrate_acceptance( jetCut,
                                                       current_distribution );
        double direct_jet_fail = ( 1.0 - direct_jet_pass );


        double
        to_taus_BR
        = ( to_Higgs_BR
            * mediating_particle->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                                 CppSLHA::PDG_code::tau_lepton,
                                            -CppSLHA::PDG_code::tau_lepton ) );

        double tau_muon_pass = 0.0;
        double tau_muon_fail = 0.0;
        if( lhcFaserGlobal::negligibleBr
            < ( to_taus_BR
                * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR ) )
          {

            current_distribution = &direct_lepton;
            // for EWSB scalars, half the tau leptons are 1 chirality, the
            // other half are the other chirality.
            visible_tau_decay_product hard_tau_muon( current_distribution,
                                          shortcut->getHardMuonFromTau() );
            visible_tau_decay_product soft_tau_muon( current_distribution,
                                          shortcut->getSoftMuonFromTau() );
            current_distribution = &hard_tau_muon;
            tau_muon_pass
            += ( 0.5 * integrate_acceptance( primary_cut,
                                             current_distribution ) );
            tau_muon_fail
            += ( 0.5 * ( 1.0 - integrate_acceptance( secondary_cut,
                                                    current_distribution ) ) );
            current_distribution = &soft_tau_muon;
            tau_muon_pass
            += ( 0.5 * integrate_acceptance( primary_cut,
                                             current_distribution ) );
            tau_muon_fail
            += ( 0.5 * ( 1.0 - integrate_acceptance( secondary_cut,
                                                    current_distribution ) ) );

          }
        double tau_pion_pass = 0.0;
        if( lhcFaserGlobal::negligibleBr
            < ( to_taus_BR
                * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR ) )
          {

            current_distribution = &direct_lepton;
            visible_tau_decay_product soft_tau_pion( current_distribution,
                                          shortcut->getSoftPionFromTau() );
            visible_tau_decay_product hard_tau_pion( current_distribution,
                                          shortcut->getHardPionFromTau() );
            current_distribution = &soft_tau_pion;
            tau_pion_pass
            += ( 0.5 * integrate_acceptance( jetCut,
                                             current_distribution ) );
            current_distribution = &hard_tau_pion;
            tau_pion_pass
            += ( 0.5 * integrate_acceptance( jetCut,
                                             current_distribution ) );

          }
        double tau_pion_fail = ( 1.0 - tau_pion_pass );


        double
        Higgs_to_jets_BR
        = ( mediating_particle->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                                       CppSLHA::PDG_code::down,
                                                     -CppSLHA::PDG_code::down )
            + mediating_particle->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                                         CppSLHA::PDG_code::up,
                                                       -CppSLHA::PDG_code::up )
            + mediating_particle->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                                    CppSLHA::PDG_code::strange,
                                                  -CppSLHA::PDG_code::strange )
            + mediating_particle->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                                      CppSLHA::PDG_code::charm,
                                                    -CppSLHA::PDG_code::charm )
            + mediating_particle->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                                     CppSLHA::PDG_code::bottom,
                                                -CppSLHA::PDG_code::bottom ) );

        double
        BR_to_two_tau_jets
        = ( to_taus_BR
            * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
            * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR );
        *two_jets_no_leptons
        += ( to_Higgs_BR * Higgs_to_jets_BR
             * ( direct_jet_pass * direct_jet_pass )
             + BR_to_two_tau_jets * ( tau_pion_pass * tau_pion_pass ) );

        double
        BR_to_one_jet_one_muon
        = ( to_taus_BR
            * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
            * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR );
        double one_jet_one_negative_muon_value
        = ( BR_to_one_jet_one_muon
            * ( tau_pion_pass * tau_muon_pass ) );
        *oneJetOneNegativeMuon
        += one_jet_one_negative_muon_value;
        *oneJetOnePositiveMuon
        += one_jet_one_negative_muon_value;

        double
        BR_to_one_jet_one_electron
        = ( to_taus_BR
            * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
            * CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR );
        double one_jet_one_negative_electron_value
        = ( BR_to_one_jet_one_electron
            * ( tau_pion_pass * tau_muon_pass ) );
        *oneJetOneNegativeElectron
        += one_jet_one_negative_electron_value;
        *oneJetOnePositiveElectron
        += one_jet_one_negative_electron_value;

        double
        Higgs_to_muons_BR
        = mediating_particle->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                                       CppSLHA::PDG_code::muon,
                                                    -CppSLHA::PDG_code::muon );

        *OSSF_minus_OSDF_leptons
        += ( to_Higgs_BR * Higgs_to_muons_BR
             * ( direct_muon_pass * direct_muon_pass ) );

        double
        BR_to_two_tau_muons
        = ( to_taus_BR
            * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR
            * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR );
        *no_jets_two_muons
        += ( to_Higgs_BR * Higgs_to_muons_BR
             * ( direct_muon_pass * direct_muon_pass )
             + BR_to_two_tau_muons * ( tau_muon_pass * tau_muon_pass ) );

        double
        BR_to_two_tau_electrons
        = ( to_taus_BR
            * CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR
            * CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR );
        *no_jets_two_electrons
        += ( BR_to_two_tau_electrons * ( tau_muon_pass * tau_muon_pass ) );
        // the direct decay to electrons is too small to be worth recording.

        double
        BR_to_one_muon_one_electron
        = ( to_taus_BR
            * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR
            * CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR );
        double no_jets_one_muon_one_electron_value
        = ( BR_to_one_muon_one_electron
             * ( tau_muon_pass * tau_muon_pass ) );
        *no_jets_one_negative_muon_one_positive_electron
        += no_jets_one_muon_one_electron_value;
        *no_jets_one_positive_muon_one_negative_electron
        += no_jets_one_muon_one_electron_value;

        *oneJetZeroLeptons
        += ( 2.0 *
             ( to_Higgs_BR * Higgs_to_jets_BR
              * ( direct_jet_pass * direct_jet_fail )
              + BR_to_two_tau_jets
                * ( tau_pion_pass * tau_pion_fail ) ) );
        // factor of 2.0 accounts for either jet particle or antiparticle
        // passing & the other failing, or the other way around.

        double no_jets_one_negative_muon_value
        = ( ( to_Higgs_BR * Higgs_to_muons_BR
              * ( direct_muon_pass * direct_muon_fail )
                  + BR_to_one_jet_one_muon
                    * ( tau_muon_pass * tau_pion_fail )
                  + BR_to_two_tau_muons
                   * ( tau_muon_pass * tau_muon_fail ) ) );
        *zeroJetsOneNegativeMuon += no_jets_one_negative_muon_value;
        *zeroJetsOnePositiveMuon += no_jets_one_negative_muon_value;

        double no_jets_one_negative_electron_value
        = ( ( BR_to_one_jet_one_electron
              * ( tau_muon_pass * tau_pion_fail )
              + BR_to_two_tau_electrons
                * ( tau_muon_pass * tau_muon_fail ) ) );
        *zeroJetsOneNegativeElectron += no_jets_one_negative_electron_value;
        *zeroJetsOnePositiveElectron += no_jets_one_negative_electron_value;


        *zeroJetsZeroLeptons
        += ( to_Higgs_BR
            * ( Higgs_to_jets_BR
                * ( direct_jet_fail * direct_jet_fail )
                + Higgs_to_muons_BR
                  * ( direct_muon_fail * direct_muon_fail ) )
             + BR_to_two_tau_jets * ( tau_pion_fail * tau_pion_fail )
             + 2.0 * ( BR_to_one_jet_one_electron + BR_to_one_jet_one_muon )
                   * ( tau_muon_fail * tau_pion_fail )
             + ( BR_to_two_tau_electrons + 2.0 * BR_to_one_muon_one_electron
                 + BR_to_two_tau_muons ) * ( tau_muon_fail * tau_muon_fail ) );

      }

  }




  void
  charge_summed_neutralino_three_body::calculate()

  {

    // we start by checking the BR of the channel to the various SM particle
    // final states:
    double
    to_electrons_BR
    = decaying_EWino->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                             CppSLHA::PDG_code::neutralino_one,
                                                   CppSLHA::PDG_code::electron,
                                                -CppSLHA::PDG_code::electron );
    double
    to_muons_BR
    = decaying_EWino->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                             CppSLHA::PDG_code::neutralino_one,
                                                       CppSLHA::PDG_code::muon,
                                                    -CppSLHA::PDG_code::muon );
    double
    to_taus_BR
    = decaying_EWino->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                             CppSLHA::PDG_code::neutralino_one,
                                                 CppSLHA::PDG_code::tau_lepton,
                                              -CppSLHA::PDG_code::tau_lepton );
    double
    to_jets_BR
    = ( decaying_EWino->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                             CppSLHA::PDG_code::neutralino_one,
                                                       CppSLHA::PDG_code::down,
                                                     -CppSLHA::PDG_code::down )
        + decaying_EWino->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                             CppSLHA::PDG_code::neutralino_one,
                                             CppSLHA::PDG_code::up,
                                             -CppSLHA::PDG_code::up )
        + decaying_EWino->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_exact_match(
                                             CppSLHA::PDG_code::neutralino_one,
                                                    CppSLHA::PDG_code::strange,
                                                  -CppSLHA::PDG_code::strange )
        + decaying_EWino->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                             CppSLHA::PDG_code::neutralino_one,
                                                      CppSLHA::PDG_code::charm,
                                                    -CppSLHA::PDG_code::charm )
        + decaying_EWino->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                             CppSLHA::PDG_code::neutralino_one,
                                                     CppSLHA::PDG_code::bottom,
                                                -CppSLHA::PDG_code::bottom ) );

    double direct_electron_pass = 0.0;
    double direct_electron_fail = 0.0;
    if( lhcFaserGlobal::negligibleBr < to_electrons_BR )
      // if the branching ratio into this channel is not negligible...
      {

        neutralino_three_body_decay direct_lepton( shortcut->getCppSlha(),
                                                   decaying_scolored,
                    kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                                   decaying_EWino,
                                                shortcut->getNeutralinoOne(),
                                                   shortcut->getSelectronL(),
                                                 shortcut->getSelectronR() );
        current_distribution = &direct_lepton;
        direct_electron_pass = integrate_acceptance( primary_cut,
                                                     current_distribution );
        direct_electron_fail
        = ( 1.0 - integrate_acceptance( secondary_cut,
                                        current_distribution ) );

      }

    double direct_muon_pass = 0.0;
    double direct_muon_fail = 0.0;
    if( lhcFaserGlobal::negligibleBr < to_muons_BR )
      // if the branching ratio into this channel is not negligible...
      {

        neutralino_three_body_decay direct_lepton( shortcut->getCppSlha(),
                                             decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                             decaying_EWino,
                                             shortcut->getNeutralinoOne(),
                                             shortcut->getSmuonL(),
                                             shortcut->getSmuonR() );
        current_distribution = &direct_lepton;
        direct_muon_pass = integrate_acceptance( primary_cut,
                                                 current_distribution );
        direct_muon_fail = ( 1.0 - integrate_acceptance( secondary_cut,
                                                      current_distribution ) );

      }

    double direct_jet_pass = 0.0;
    if( lhcFaserGlobal::negligibleBr < to_jets_BR )
      // if the branching ratio into this channel is not negligible...
      {

        neutralino_three_body_decay direct_jet( shortcut->getCppSlha(),
                                                decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                                decaying_EWino,
                                                shortcut->getNeutralinoOne(),
                                                shortcut->getSdownL(),
                                                shortcut->getSdownR() );
        // actually, here I should use every nontop flavor of squark
        // separately, but it's unlikely to make any significant difference.
        current_distribution = &direct_jet;
        direct_jet_pass = integrate_acceptance( jet_cut,
                                                current_distribution );

      }
    double direct_jet_fail = ( 1.0 - direct_jet_pass );

    double tau_muon_pass = 0.0;
    double tau_muon_fail = 0.0;
    double tau_pion_pass = 0.0;
    if( lhcFaserGlobal::negligibleBr < to_taus_BR )
      // if the branching ratio into this channel is not negligible...
      {

        neutralino_three_body_decay direct_lepton( shortcut->getCppSlha(),
                                                   decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                                   decaying_EWino,
                                                shortcut->getNeutralinoOne(),
                                                   shortcut->getStauOne(),
                                                   shortcut->getStauTwo() );
        current_distribution = &direct_lepton;
        visible_tau_decay_product hard_tau_muon( current_distribution,
                                          shortcut->getHardMuonFromTau() );
        visible_tau_decay_product soft_tau_muon( current_distribution,
                                          shortcut->getSoftMuonFromTau() );
        visible_tau_decay_product hard_tau_pion( current_distribution,
                                          shortcut->getHardPionFromTau() );
        visible_tau_decay_product soft_tau_pion( current_distribution,
                                          shortcut->getSoftPionFromTau() );

        /* at the moment, I have no plans to provide a tau left-handedness
         * breakdown from a 3-body decay, so I assume that 50% of them are of
         * each handedness.
         */
        current_distribution = &hard_tau_muon;
        leptonEnergyDistribution* other_distribution = &soft_tau_muon;
        tau_muon_pass = ( 0.5 * integrate_acceptance( primary_cut,
                                                      current_distribution )
                          + 0.5 * integrate_acceptance( primary_cut,
                                                        other_distribution ) );
        tau_muon_fail = ( 0.5 * integrate_acceptance( secondary_cut,
                                                      current_distribution )
                          + 0.5 * integrate_acceptance( secondary_cut,
                                                        other_distribution ) );
        current_distribution = &hard_tau_pion;
        other_distribution = &soft_tau_pion;
        tau_pion_pass = ( 0.5 * integrate_acceptance( jet_cut,
                                                      current_distribution )
                          + 0.5 * integrate_acceptance( jet_cut,
                                                        other_distribution ) );

      }
    double tau_pion_fail = ( 1.0 - tau_pion_pass );

    double
    BR_to_two_tau_jets
    = ( to_taus_BR
        * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
        * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR );
    *two_jets_no_leptons
    += ( to_jets_BR
         * ( direct_jet_pass * direct_jet_pass )
         + BR_to_two_tau_jets * ( tau_pion_pass * tau_pion_pass ) );

    double
    BR_to_one_jet_one_muon
    = ( to_taus_BR
        * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
        * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR );
    *one_jet_one_muon
    += ( 2.0 * BR_to_one_jet_one_muon
         * ( tau_pion_pass * tau_muon_pass ) );
    // factor of 2.0 accounts for either tau lepton to muon & tau
    // antilepton to pion, or the other way around.

    double
    BR_to_one_jet_one_electron
    = ( to_taus_BR
        * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
        * CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR );
    *one_jet_one_electron
    += ( 2.0 * BR_to_one_jet_one_electron
         * ( tau_pion_pass * tau_muon_pass ) );
        // factor of 2.0 accounts for either tau lepton to electron & tau
        // antilepton to pion, or the other way around.

    *OSSF_minus_OSDF_leptons
    += ( to_muons_BR * ( direct_muon_pass * direct_muon_pass )
         + to_electrons_BR * ( direct_electron_pass * direct_electron_pass ) );

    double
    BR_to_two_tau_muons
    = ( to_taus_BR
        * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR
        * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR );
    *no_jets_two_muons
    += ( to_muons_BR * ( direct_muon_pass * direct_muon_pass )
         + BR_to_two_tau_muons * ( tau_muon_pass * tau_muon_pass ) );

    double
    BR_to_two_tau_electrons
    = ( to_taus_BR
        * CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR
        * CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR );
    *no_jets_two_electrons
    += ( to_electrons_BR * ( direct_electron_pass * direct_electron_pass )
         + BR_to_two_tau_electrons * ( tau_muon_pass * tau_muon_pass ) );

    double
    BR_to_one_muon_one_electron
    = ( to_taus_BR
        * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR
        * CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR );
    *no_jets_one_muon_one_electron
    += ( 2.0 * BR_to_one_muon_one_electron
         * ( tau_muon_pass * tau_muon_pass ) );
    // factor of 2.0 accounts for either tau lepton or antilepton to muon,
    // other to electron.

    *one_jet_no_leptons
    += ( 2.0 * ( to_jets_BR * ( direct_jet_pass * direct_jet_fail )
                 + BR_to_two_tau_jets * ( tau_pion_pass * tau_pion_fail ) ) );
        // factor of 2.0 accounts for either jet particle or antiparticle
        // passing & the other failing, or the other way around.

    *no_jets_one_muon
    += ( 2.0 * ( to_muons_BR * ( direct_muon_pass * direct_muon_fail )
                 + BR_to_one_jet_one_muon * ( tau_muon_pass * tau_pion_fail )
                 + BR_to_two_tau_muons * ( tau_muon_pass * tau_muon_fail ) ) );
    // factor of 2.0 accounts for either muon or antimuon passing, basically.

    *no_jets_one_electron
    += ( 2.0 * ( to_electrons_BR
                 * ( direct_electron_pass * direct_electron_fail )
                 + BR_to_one_jet_one_electron
                   * ( tau_muon_pass * tau_pion_fail )
                 + BR_to_two_tau_electrons
                   * ( tau_muon_pass * tau_muon_fail ) ) );
    // factor of 2.0 accounts for either electron or antielectron passing,
    // basically.

    *no_jets_no_leptons
    += ( to_jets_BR
         * ( direct_jet_fail * direct_jet_fail )
         + to_electrons_BR
           * ( direct_electron_fail * direct_electron_fail )
         + to_muons_BR
           * ( direct_muon_fail * direct_muon_fail )
         + BR_to_two_tau_jets * ( tau_pion_fail * tau_pion_fail )
         + 2.0 * ( BR_to_one_jet_one_electron + BR_to_one_jet_one_muon )
               * ( tau_muon_fail * tau_pion_fail )
         + ( BR_to_two_tau_electrons + 2.0 * BR_to_one_muon_one_electron
             + BR_to_two_tau_muons ) * ( tau_muon_fail * tau_muon_fail ) );

  }






  void
  chargino_to_semu::calculate()
  {

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: chargino_to_semu::calculate() called.";
    std::cout << std::endl;**/

    // we start by checking the BR of the (only) channel:
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
    exclusiveBrCalculator*
    first_BR
    = shortcut->getExclusiveBrCalculator( decaying_EWino,
                                             mediating_particle,
                                             not_particle_antiparticle_flip,
                                             shortcut->getEmptyList() );
    exclusiveBrCalculator*
    second_BR
    = shortcut->getExclusiveBrCalculator( mediating_particle,
                                             shortcut->getNeutralinoOne(),
                                             true,
                                             shortcut->getEmptyList() );
    double cascade_BR = ( first_BR->getBr() * second_BR->getBr() );

    if( lhcFaserGlobal::negligibleBr < cascade_BR )
      // if the branching ratio into this channel is not negligible...
      {

        // we just have a far lepton:
        same_chirality_far_muon same_lepton( shortcut->getCppSlha(),
                                             decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                             decaying_EWino,
                                             mediating_particle,
                                             shortcut->getNeutralinoOne() );
        opposite_chirality_far_muon opposite_lepton( shortcut->getCppSlha(),
                                                     decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                                     decaying_EWino,
                                                     mediating_particle,
                                              shortcut->getNeutralinoOne() );
        current_distribution = &same_lepton;

        double
        same_handedness
        = shortcut->quarkOrLeptonLeftHandedness(
                                             decaying_scolored->get_PDG_code(),
                                              decaying_EWino->get_PDG_code() );
        double
        same_pass = integrate_acceptance( primary_cut,
                                          current_distribution );
        double
        same_fail = ( 1.0 - integrate_acceptance( secondary_cut,
                                                  current_distribution ) );
        current_distribution = &opposite_lepton;
        double opposite_handedness = ( 1.0 - same_handedness );
        double
        opposite_pass = integrate_acceptance( primary_cut,
                                              current_distribution );
        double
        opposite_fail = ( 1.0 - integrate_acceptance( secondary_cut,
                                                      current_distribution ) );

        // debugging:
        /**std::cout
        << std::endl
        << "about to set values:";
        std::cout << std::endl;**/

        *no_jets_one_lepton
        += ( cascade_BR * ( same_handedness * same_pass
                            + opposite_handedness * opposite_pass ) );

        // debugging:
        /**std::cout
        << std::endl
        << "*no_jets_one_lepton set to" << *no_jets_one_lepton;
        std::cout << std::endl;**/

        *no_jets_no_leptons
        += ( cascade_BR * ( same_handedness * same_fail
                            + opposite_handedness * opposite_fail ) );

        // debugging:
        /**std::cout
        << std::endl
        << "*no_jets_no_leptons set to" << *no_jets_no_leptons;
        std::cout << std::endl;**/

      }

  }



  void
  chargino_to_stau::calculate()
  {

    // we start by checking the BR of the (only) channel:
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
    exclusiveBrCalculator*
    first_BR
    = shortcut->getExclusiveBrCalculator( decaying_EWino,
                                             mediating_particle,
                                             not_particle_antiparticle_flip,
                                             shortcut->getEmptyList() );
    exclusiveBrCalculator*
    second_BR
    = shortcut->getExclusiveBrCalculator( mediating_particle,
                                             shortcut->getNeutralinoOne(),
                                             true,
                                             shortcut->getEmptyList() );
    double cascade_BR = ( first_BR->getBr() * second_BR->getBr() );

    if( lhcFaserGlobal::negligibleBr < cascade_BR )
      // if the branching ratio into this channel is not negligible...
      {

        // we just have a far lepton:
        same_chirality_far_muon same_tau( shortcut->getCppSlha(),
                                          decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                          decaying_EWino,
                                          mediating_particle,
                                          shortcut->getNeutralinoOne() );
        opposite_chirality_far_muon opposite_tau( shortcut->getCppSlha(),
                                                  decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                                  decaying_EWino,
                                                  mediating_particle,
                                              shortcut->getNeutralinoOne() );

        double
        quark_left_handedness
        = shortcut->quarkOrLeptonLeftHandedness(
                                             decaying_scolored->get_PDG_code(),
                                              decaying_EWino->get_PDG_code() );
        double quark_right_handedness = ( 1.0 - quark_left_handedness );
        double
        tau_left_handedness
        = shortcut->quarkOrLeptonLeftHandedness(
                                            mediating_particle->get_PDG_code(),
                                              decaying_EWino->get_PDG_code() );
        double tau_right_handedness = ( 1.0 - tau_left_handedness );

        double muon_pass = 0.0;
        double muon_fail = 0.0;
        double pion_pass = 0.0;
        double pion_fail = 0.0;
        double pion_acceptance;
        if( lhcFaserGlobal::negligibleBr
            < ( cascade_BR
                * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR ) )
          {

            current_distribution = &same_tau;
            visible_tau_decay_product hard_same_muon( current_distribution,
                                          shortcut->getHardMuonFromTau() );
            leptonEnergyDistribution* visible_distribution = &hard_same_muon;
            muon_pass += ( quark_left_handedness * tau_left_handedness
                           * integrate_acceptance( primary_cut,
                                                   visible_distribution ) );
            muon_fail
            += ( quark_left_handedness * tau_left_handedness
                 * ( 1.0 - integrate_acceptance( secondary_cut,
                                                 visible_distribution ) ) );
            visible_tau_decay_product soft_same_muon( current_distribution,
                                          shortcut->getSoftMuonFromTau() );
            visible_distribution = &soft_same_muon;
            muon_pass += ( quark_right_handedness * tau_right_handedness
                           * integrate_acceptance( primary_cut,
                                                   visible_distribution ) );
            muon_fail
            += ( quark_right_handedness * tau_right_handedness
                 * ( 1.0 - integrate_acceptance( secondary_cut,
                                                 visible_distribution ) ) );
            visible_tau_decay_product soft_same_pion( current_distribution,
                                          shortcut->getSoftPionFromTau() );
            visible_distribution = &soft_same_pion;
            pion_acceptance = integrate_acceptance( jet_cut,
                                                    visible_distribution );
            pion_pass += ( quark_left_handedness * tau_left_handedness
                           * pion_acceptance );
            pion_fail
            += ( quark_left_handedness * tau_left_handedness
                 * ( 1.0 - pion_acceptance ) );
            visible_tau_decay_product hard_same_pion( current_distribution,
                                          shortcut->getHardPionFromTau() );
            visible_distribution = &hard_same_pion;
            pion_acceptance = integrate_acceptance( jet_cut,
                                                    visible_distribution );
            pion_pass += ( quark_right_handedness * tau_right_handedness
                           * pion_acceptance );
            pion_fail
            += ( quark_right_handedness * tau_right_handedness
                 * ( 1.0 - pion_acceptance ) );

            current_distribution = &opposite_tau;
            visible_tau_decay_product hard_opposite_muon( current_distribution,
                                          shortcut->getHardMuonFromTau() );
            visible_distribution = &hard_opposite_muon;
            muon_pass += ( quark_right_handedness * tau_left_handedness
                           * integrate_acceptance( primary_cut,
                                                   visible_distribution ) );
            muon_fail
            += ( quark_right_handedness * tau_left_handedness
                 * ( 1.0 - integrate_acceptance( secondary_cut,
                                                 visible_distribution ) ) );
            visible_tau_decay_product soft_opposite_muon( current_distribution,
                                          shortcut->getSoftMuonFromTau() );
            visible_distribution = &soft_opposite_muon;
            muon_pass += ( quark_left_handedness * tau_right_handedness
                           * integrate_acceptance( primary_cut,
                                                   visible_distribution ) );
            muon_fail
            += ( quark_left_handedness * tau_right_handedness
                 * ( 1.0 - integrate_acceptance( secondary_cut,
                                                 visible_distribution ) ) );
            visible_tau_decay_product soft_opposite_pion( current_distribution,
                                          shortcut->getSoftPionFromTau() );
            visible_distribution = &soft_same_pion;
            pion_acceptance = integrate_acceptance( jet_cut,
                                                    visible_distribution );
            pion_pass += ( quark_right_handedness * tau_left_handedness
                           * pion_acceptance );
            pion_fail
            += ( quark_right_handedness * tau_left_handedness
                 * ( 1.0 - pion_acceptance ) );
            visible_tau_decay_product hard_opposite_pion( current_distribution,
                                          shortcut->getHardPionFromTau() );
            visible_distribution = &hard_same_pion;
            pion_acceptance = integrate_acceptance( jet_cut,
                                                    visible_distribution );
            pion_pass += ( quark_left_handedness * tau_right_handedness
                           * pion_acceptance );
            pion_fail
            += ( quark_left_handedness * tau_right_handedness
                 * ( 1.0 - pion_acceptance ) );

          }

        *one_jet_no_leptons
        += ( cascade_BR * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
             * pion_pass );

        *no_jets_one_electron
        += ( cascade_BR
             * CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR
             * muon_pass );

        *no_jets_one_muon
        += ( cascade_BR * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR
             * muon_pass );

        *no_jets_no_leptons
        += ( cascade_BR
             * ( CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
                 * pion_fail
                 + ( CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR
                     + CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR )
                   * muon_fail ) );

      }

  }



  void
  chargino_to_emu_sneutrino::calculate()
  {

    // we start by checking the BR of the (only) channel:
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
    exclusiveBrCalculator*
    first_BR
    = shortcut->getExclusiveBrCalculator( decaying_EWino,
                                             mediating_particle,
                                             not_particle_antiparticle_flip,
                                             shortcut->getEmptyList() );
    exclusiveBrCalculator*
    second_BR
    = shortcut->getExclusiveBrCalculator( mediating_particle,
                                             shortcut->getNeutralinoOne(),
                                             true,
                                             shortcut->getEmptyList() );
    double cascade_BR = ( first_BR->getBr() * second_BR->getBr() );

    if( lhcFaserGlobal::negligibleBr < cascade_BR )
      // if the branching ratio into this channel is not negligible...
      {

        // we just have a near lepton:
        same_chirality_near_muon same_lepton( shortcut->getCppSlha(),
                                              decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                              decaying_EWino,
                                              mediating_particle );
        opposite_chirality_near_muon opposite_lepton( shortcut->getCppSlha(),
                                                      decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                                      decaying_EWino,
                                                      mediating_particle );
        current_distribution = &same_lepton;

        double
        quark_left_handedness = shortcut->quarkOrLeptonLeftHandedness(
                                             decaying_scolored->get_PDG_code(),
                                              decaying_EWino->get_PDG_code() );
        double
        lepton_left_handedness = shortcut->quarkOrLeptonLeftHandedness(
                                                decaying_EWino->get_PDG_code(),
                                          mediating_particle->get_PDG_code() );
        double
        same_handedness = ( quark_left_handedness * lepton_left_handedness
                            + ( 1.0 - quark_left_handedness )
                              *( 1.0 - lepton_left_handedness ) );
        double
        same_pass = integrate_acceptance( primary_cut,
                                          current_distribution );
        double
        same_fail = ( 1.0 - integrate_acceptance( secondary_cut,
                                                  current_distribution ) );
        current_distribution = &opposite_lepton;
        double opposite_handedness = ( 1.0 - same_handedness );
        double
        opposite_pass = integrate_acceptance( primary_cut,
                                              current_distribution );
        double
        opposite_fail = ( 1.0 - integrate_acceptance( secondary_cut,
                                                      current_distribution ) );

        *no_jets_one_lepton
        += ( cascade_BR * ( same_handedness * same_pass
                            + opposite_handedness * opposite_pass ) );

        *no_jets_no_leptons
        += ( cascade_BR * ( same_handedness * same_fail
                            + opposite_handedness * opposite_fail ) );

      }

  }




  void
  chargino_to_tau_sneutrino::calculate()
  {

    // we start by checking the BR of the (only) channel:
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
    exclusiveBrCalculator*
    first_BR
    = shortcut->getExclusiveBrCalculator( decaying_EWino,
                                             mediating_particle,
                                             not_particle_antiparticle_flip,
                                             shortcut->getEmptyList() );
    exclusiveBrCalculator*
    second_BR
    = shortcut->getExclusiveBrCalculator( mediating_particle,
                                             shortcut->getNeutralinoOne(),
                                             true,
                                             shortcut->getEmptyList() );
    double cascade_BR = ( first_BR->getBr() * second_BR->getBr() );

    if( lhcFaserGlobal::negligibleBr < cascade_BR )
      // if the branching ratio into this channel is not negligible...
      {

        // we just have a near lepton:
        same_chirality_near_muon same_tau( shortcut->getCppSlha(),
                                           decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                           decaying_EWino,
                                           mediating_particle );
        opposite_chirality_near_muon opposite_tau( shortcut->getCppSlha(),
                                                   decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                                   decaying_EWino,
                                                   mediating_particle );

        double
        quark_left_handedness
        = shortcut->quarkOrLeptonLeftHandedness(
                                             decaying_scolored->get_PDG_code(),
                                              decaying_EWino->get_PDG_code() );
        double quark_right_handedness = ( 1.0 - quark_left_handedness );
        double
        tau_left_handedness
        = shortcut->quarkOrLeptonLeftHandedness(
                                            mediating_particle->get_PDG_code(),
                                              decaying_EWino->get_PDG_code() );
        double tau_right_handedness = ( 1.0 - tau_left_handedness );

        double muon_pass = 0.0;
        double muon_fail = 0.0;
        double pion_pass = 0.0;
        double pion_fail = 0.0;
        double pion_acceptance;
        if( lhcFaserGlobal::negligibleBr
            < ( cascade_BR
                * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR ) )
          {

            current_distribution = &same_tau;
            visible_tau_decay_product hard_same_muon( current_distribution,
                                          shortcut->getHardMuonFromTau() );
            leptonEnergyDistribution* visible_distribution = &hard_same_muon;
            muon_pass += ( quark_left_handedness * tau_left_handedness
                           * integrate_acceptance( primary_cut,
                                                   visible_distribution ) );
            muon_fail
            += ( quark_left_handedness * tau_left_handedness
                 * ( 1.0 - integrate_acceptance( secondary_cut,
                                                 visible_distribution ) ) );
            visible_tau_decay_product soft_same_muon( current_distribution,
                                          shortcut->getSoftMuonFromTau() );
            visible_distribution = &soft_same_muon;
            muon_pass += ( quark_right_handedness * tau_right_handedness
                           * integrate_acceptance( primary_cut,
                                                   visible_distribution ) );
            muon_fail
            += ( quark_right_handedness * tau_right_handedness
                 * ( 1.0 - integrate_acceptance( secondary_cut,
                                                 visible_distribution ) ) );
            visible_tau_decay_product soft_same_pion( current_distribution,
                                          shortcut->getSoftPionFromTau() );
            visible_distribution = &soft_same_pion;
            pion_acceptance = integrate_acceptance( jet_cut,
                                                    visible_distribution );
            pion_pass += ( quark_left_handedness * tau_left_handedness
                           * pion_acceptance );
            pion_fail
            += ( quark_left_handedness * tau_left_handedness
                 * ( 1.0 - pion_acceptance ) );
            visible_tau_decay_product hard_same_pion( current_distribution,
                                          shortcut->getHardPionFromTau() );
            visible_distribution = &hard_same_pion;
            pion_acceptance = integrate_acceptance( jet_cut,
                                                    visible_distribution );
            pion_pass += ( quark_right_handedness * tau_right_handedness
                           * pion_acceptance );
            pion_fail
            += ( quark_right_handedness * tau_right_handedness
                 * ( 1.0 - pion_acceptance ) );

            current_distribution = &opposite_tau;
            visible_tau_decay_product hard_opposite_muon( current_distribution,
                                          shortcut->getHardMuonFromTau() );
            visible_distribution = &hard_opposite_muon;
            muon_pass += ( quark_right_handedness * tau_left_handedness
                           * integrate_acceptance( primary_cut,
                                                   visible_distribution ) );
            muon_fail
            += ( quark_right_handedness * tau_left_handedness
                 * ( 1.0 - integrate_acceptance( secondary_cut,
                                                 visible_distribution ) ) );
            visible_tau_decay_product soft_opposite_muon( current_distribution,
                                          shortcut->getSoftMuonFromTau() );
            visible_distribution = &soft_opposite_muon;
            muon_pass += ( quark_left_handedness * tau_right_handedness
                           * integrate_acceptance( primary_cut,
                                                   visible_distribution ) );
            muon_fail
            += ( quark_left_handedness * tau_right_handedness
                 * ( 1.0 - integrate_acceptance( secondary_cut,
                                                 visible_distribution ) ) );
            visible_tau_decay_product soft_opposite_pion( current_distribution,
                                          shortcut->getSoftPionFromTau() );
            visible_distribution = &soft_same_pion;
            pion_acceptance = integrate_acceptance( jet_cut,
                                                    visible_distribution );
            pion_pass += ( quark_right_handedness * tau_left_handedness
                           * pion_acceptance );
            pion_fail
            += ( quark_right_handedness * tau_left_handedness
                 * ( 1.0 - pion_acceptance ) );
            visible_tau_decay_product hard_opposite_pion( current_distribution,
                                          shortcut->getHardPionFromTau() );
            visible_distribution = &hard_same_pion;
            pion_acceptance = integrate_acceptance( jet_cut,
                                                    visible_distribution );
            pion_pass += ( quark_left_handedness * tau_right_handedness
                           * pion_acceptance );
            pion_fail
            += ( quark_left_handedness * tau_right_handedness
                 * ( 1.0 - pion_acceptance ) );

          }

        *one_jet_no_leptons
        += ( cascade_BR * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
             * pion_pass );

        *no_jets_one_electron
        += ( cascade_BR
             * CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR
             * muon_pass );

        *no_jets_one_muon
        += ( cascade_BR * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR
             * muon_pass );

        *no_jets_no_leptons
        += ( cascade_BR
             * ( CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
                 * pion_fail
                 + ( CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR
                     + CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR )
                   * muon_fail ) );

      }

  }



  void
  chargino_to_W::calculate()
  {

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: chargino_to_W::calculate() called.";
    std::cout << std::endl;**/

    // we start by checking the BR of the channel to the W:
    double
    to_W_BR
    = decaying_EWino->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                            mediating_particle->get_PDG_code(),
                                           CppSLHA::PDG_code::neutralino_one );
    // debugging:
    /**std::cout
    << std::endl
    << "to_W_BR = " << to_W_BR;
    std::cout << std::endl;**/

    if( lhcFaserGlobal::negligibleBr < to_W_BR )
      // if the branching ratio into this channel is not negligible...
      {

        // we start with the direct distributions:
        W_minus_direct_jet direct_jet( shortcut->getCppSlha(),
                                       decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                       decaying_EWino,
                                       mediating_particle,
                                       shortcut->getNeutralinoOne() );
        // debugging:
        /**std::cout
        << std::endl
        << "&direct_jet = " << &direct_jet;
        std::cout << std::endl;**/
        W_minus_plus_up_L_type_muon L_direct_lepton( shortcut->getCppSlha(),
                                                     decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                                     decaying_EWino,
                                                     mediating_particle,
                                              shortcut->getNeutralinoOne() );
        // debugging:
        /**std::cout
        << std::endl
        << "&L_direct_lepton = " << &L_direct_lepton;
        std::cout << std::endl;**/
        W_minus_plus_up_R_type_muon R_direct_lepton( shortcut->getCppSlha(),
                                                     decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                                     decaying_EWino,
                                                     mediating_particle,
                                              shortcut->getNeutralinoOne() );
        // debugging:
        /**std::cout
        << std::endl
        << "&R_direct_lepton = " << &R_direct_lepton;
        std::cout << std::endl;**/
        current_distribution = &direct_jet;
        // debugging:
        /**std::cout
        << std::endl
        << "current_distribution = (&direct_jet) " << current_distribution;
        std::cout << std::endl;**/
        double direct_jet_pass = integrate_acceptance( jet_cut,
                                                       current_distribution );
        double direct_jet_fail = ( 1.0 - direct_jet_pass );

        double
        quark_left_handedness
        = shortcut->quarkOrLeptonLeftHandedness(
                                             decaying_scolored->get_PDG_code(),
                                              decaying_EWino->get_PDG_code() );
        double quark_right_handedness = ( 1.0 - quark_left_handedness );
        current_distribution = &L_direct_lepton;
        // debugging:
        /**std::cout
        << std::endl
        << "current_distribution = (&L_direct_lepton) "
        << current_distribution;
        std::cout << std::endl;**/
        double
        direct_muon_pass
        = ( quark_left_handedness * integrate_acceptance( primary_cut,
                                                      current_distribution ) );
        double
        direct_muon_fail
        = ( quark_left_handedness
            * ( 1.0 - integrate_acceptance( secondary_cut,
                                            current_distribution ) ) );
        current_distribution = &R_direct_lepton;
        direct_muon_pass
        += ( quark_right_handedness * integrate_acceptance( primary_cut,
                                                      current_distribution ) );
        direct_muon_fail
        += ( quark_right_handedness
             * ( 1.0 - integrate_acceptance( secondary_cut,
                                             current_distribution ) ) );


        double
        to_tau_BR
        = ( to_W_BR
            * CppSLHA::PDG_data::W_plus_to_neutrino_tau_antilepton_BR );

        double tau_muon_pass = 0.0;
        double tau_muon_fail = 0.0;
        if( lhcFaserGlobal::negligibleBr
            < ( to_tau_BR
                * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR ) )
          {

            current_distribution = &L_direct_lepton;
            visible_tau_decay_product L_tau_muon( current_distribution,
                                          shortcut->getHardMuonFromTau() );
            current_distribution = &L_tau_muon;
            tau_muon_pass
            += ( quark_left_handedness * integrate_acceptance( primary_cut,
                                                      current_distribution ) );
            tau_muon_fail
            += ( quark_left_handedness
                * ( 1.0 - integrate_acceptance( secondary_cut,
                                                current_distribution ) ) );
            current_distribution = &R_direct_lepton;
            visible_tau_decay_product R_tau_muon( current_distribution,
                                          shortcut->getHardMuonFromTau() );
            current_distribution = &R_tau_muon;
            tau_muon_pass
            += ( quark_right_handedness * integrate_acceptance( primary_cut,
                                                      current_distribution ) );
            tau_muon_fail
            += ( quark_right_handedness
                * ( 1.0 - integrate_acceptance( secondary_cut,
                                                current_distribution ) ) );

          }
        double tau_pion_pass = 0.0;
        if( lhcFaserGlobal::negligibleBr
            < ( to_tau_BR
                * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR ) )
          {

            current_distribution = &L_direct_lepton;
            visible_tau_decay_product L_tau_pion( current_distribution,
                                          shortcut->getSoftPionFromTau() );
            current_distribution = &L_tau_pion;
            tau_pion_pass
            += ( quark_left_handedness * integrate_acceptance( jet_cut,
                                                      current_distribution ) );
            current_distribution = &R_direct_lepton;
            visible_tau_decay_product R_tau_pion( current_distribution,
                                          shortcut->getSoftPionFromTau() );
            current_distribution = &R_tau_pion;
            tau_pion_pass
            += ( quark_right_handedness * integrate_acceptance( jet_cut,
                                                      current_distribution ) );

          }
        double tau_pion_fail = ( 1.0 - tau_pion_pass );

        // debugging:
        /**std::cout
        << std::endl
        << "about to set values:";
        std::cout << std::endl;**/

        *two_jets_no_leptons
        += ( to_W_BR * CppSLHA::PDG_data::W_plus_to_hadrons_BR
             * ( direct_jet_pass * direct_jet_pass ) );

        // debugging:
        /**std::cout
        << std::endl
        << "*two_jets_no_leptons set to" << *two_jets_no_leptons;
        std::cout << std::endl;**/

        double
        BR_to_tau_jet
        = ( to_tau_BR
            * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR );

        *one_jet_no_leptons
        += ( 2.0 * to_W_BR * CppSLHA::PDG_data::W_plus_to_hadrons_BR
                   * ( direct_jet_pass * direct_jet_fail )
             + BR_to_tau_jet * tau_pion_pass );
        // factor of 2.0 accounts for either jet particle or antiparticle
        // passing & the other failing, or the other way around.

        // debugging:
        /**std::cout
        << std::endl
        << "*one_jet_no_leptons set to" << *one_jet_no_leptons;
        std::cout << std::endl;**/

        *no_jets_one_muon
        += ( to_W_BR * CppSLHA::PDG_data::W_plus_to_neutrino_antimuon_BR
             * direct_muon_pass
             + to_tau_BR * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR
               * tau_muon_pass );

        // debugging:
        /**std::cout
        << std::endl
        << "*no_jets_one_muon set to" << *no_jets_one_muon;
        std::cout << std::endl;**/

        *no_jets_one_electron
        += ( to_W_BR * CppSLHA::PDG_data::W_plus_to_neutrino_antielectron_BR
             * direct_muon_pass
             + to_tau_BR
               * CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR
               * tau_muon_pass );

        // debugging:
        /**std::cout
        << std::endl
        << "*no_jets_one_electron set to" << *no_jets_one_electron;
        std::cout << std::endl;**/

        *no_jets_no_leptons
        += ( to_W_BR * CppSLHA::PDG_data::W_plus_to_hadrons_BR
             * ( direct_jet_fail * direct_jet_fail )
             + BR_to_tau_jet * tau_pion_fail
             + to_W_BR
               * ( CppSLHA::PDG_data::W_plus_to_neutrino_antimuon_BR
                   + CppSLHA::PDG_data::W_plus_to_neutrino_antielectron_BR )
               * direct_muon_fail
             + to_tau_BR
               * ( CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR
                   + CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR )
               * tau_muon_fail );

        // debugging:
        /**std::cout
        << std::endl
        << "*no_jets_no_leptons set to" << *no_jets_no_leptons;
        std::cout << std::endl;**/

      }

  }



  void
  chargino_to_Higgs::calculate()
  {

    // we start by checking the BR of the channel to the W:
    double
    to_Higgs_BR
    = decaying_EWino->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                            mediating_particle->get_PDG_code(),
                                           CppSLHA::PDG_code::neutralino_one );

    if( lhcFaserGlobal::negligibleBr < to_Higgs_BR )
      // if the branching ratio into this channel is not negligible...
      {

        // we start with the direct distributions:
        negatively_charged_Higgs_jet direct_jet( shortcut->getCppSlha(),
                                                 decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                                 decaying_EWino,
                                                 mediating_particle,
                                              shortcut->getNeutralinoOne() );
        negatively_charged_Higgs_muon direct_lepton( shortcut->getCppSlha(),
                                                     decaying_scolored,
                    kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                                     decaying_EWino,
                                                     mediating_particle,
                                              shortcut->getNeutralinoOne() );
        current_distribution = &direct_jet;
        double direct_jet_pass = integrate_acceptance( jet_cut,
                                                       current_distribution );
        double direct_jet_fail = ( 1.0 - direct_jet_pass );

        current_distribution = &direct_lepton;
        double
        direct_muon_pass
        = integrate_acceptance( primary_cut,
                                current_distribution );
        double
        direct_muon_fail
        = ( 1.0 - integrate_acceptance( secondary_cut,
                                        current_distribution ) );

        double
        to_tau_BR
        = ( to_Higgs_BR
            * mediating_particle->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                               CppSLHA::PDG_code::tau_neutrino,
                                            -CppSLHA::PDG_code::tau_lepton ) );
        // since the positively charged EWSB scalar has the positive PDG
        // code...

        double tau_muon_pass = 0.0;
        double tau_muon_fail = 0.0;
        if( lhcFaserGlobal::negligibleBr
            < ( to_tau_BR
                * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR ) )
          {

            current_distribution = &direct_lepton;
            visible_tau_decay_product tau_muon( current_distribution,
                                          shortcut->getHardMuonFromTau() );
            current_distribution = &tau_muon;
            tau_muon_pass
            += integrate_acceptance( primary_cut,
                                     current_distribution );
            tau_muon_fail
            += ( 1.0 - integrate_acceptance( secondary_cut,
                                             current_distribution ) );

          }
        double tau_pion_pass = 0.0;
        if( lhcFaserGlobal::negligibleBr
            < ( to_tau_BR
                * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR ) )
          {

            current_distribution = &direct_lepton;
            visible_tau_decay_product tau_pion( current_distribution,
                                          shortcut->getSoftPionFromTau() );
            current_distribution = &tau_pion;
            tau_pion_pass
            += integrate_acceptance( jet_cut,
                                     current_distribution );

          }
        double tau_pion_fail = ( 1.0 - tau_pion_pass );

        double
        Higgs_to_jets_BR
        = ( mediating_particle->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                                         CppSLHA::PDG_code::up,
                                                     -CppSLHA::PDG_code::down )
            + mediating_particle->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                                      CppSLHA::PDG_code::charm,
                                               -CppSLHA::PDG_code::strange ) );

        *two_jets_no_leptons
        += ( Higgs_to_jets_BR
             * ( direct_jet_pass * direct_jet_pass ) );

        double
        BR_to_tau_jet
        = ( to_tau_BR
            * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR );

        *one_jet_no_leptons
        += ( 2.0 * Higgs_to_jets_BR
                   * ( direct_jet_pass * direct_jet_fail )
             + BR_to_tau_jet * tau_pion_pass );
        // factor of 2.0 accounts for either jet particle or antiparticle
        // passing & the other failing, or the other way around.


        double
        Higgs_to_muon_BR
        = mediating_particle->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                              CppSLHA::PDG_code::muon_neutrino,
                                                    -CppSLHA::PDG_code::muon );

        *no_jets_one_muon
        += ( Higgs_to_muon_BR
             * direct_muon_pass
             + to_tau_BR * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR
               * tau_muon_pass );

        *no_jets_one_electron
        += ( to_tau_BR
             * CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR
             * tau_muon_pass );

        *no_jets_no_leptons
        += ( Higgs_to_jets_BR
             * ( direct_jet_fail * direct_jet_fail )
             + BR_to_tau_jet * tau_pion_fail
             + Higgs_to_muon_BR
               * direct_muon_fail
             + to_tau_BR
               * ( CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR
                   + CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR )
               * tau_muon_fail );

      }

  }

  void
  chargino_three_body::calculate()
  {

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: chargino_three_body::calculate() called.";
    std::cout << std::endl;**/

    // we start by checking the BR of the channel to the various SM particle
    // final states:
    double
    to_electron_BR
    = decaying_EWino->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                             CppSLHA::PDG_code::neutralino_one,
                                          CppSLHA::PDG_code::electron_neutrino,
                                                -CppSLHA::PDG_code::electron );
    double
    to_muon_BR
    = decaying_EWino->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                             CppSLHA::PDG_code::neutralino_one,
                                              CppSLHA::PDG_code::muon_neutrino,
                                                    -CppSLHA::PDG_code::muon );
    double
    to_tau_BR
    = decaying_EWino->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                             CppSLHA::PDG_code::neutralino_one,
                                               CppSLHA::PDG_code::tau_neutrino,
                                              -CppSLHA::PDG_code::tau_lepton );
    double
    to_jets_BR
    = ( decaying_EWino->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                             CppSLHA::PDG_code::neutralino_one,
                                                         CppSLHA::PDG_code::up,
                                                     -CppSLHA::PDG_code::down )
        + decaying_EWino->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                             CppSLHA::PDG_code::neutralino_one,
                                                      CppSLHA::PDG_code::charm,
                                                -CppSLHA::PDG_code::strange ) );

    double direct_electron_pass = 0.0;
    double direct_electron_fail = 0.0;
    if( lhcFaserGlobal::negligibleBr < to_electron_BR )
      // if the branching ratio into this channel is not negligible...
      {

        chargino_three_body_decay direct_lepton( shortcut->getCppSlha(),
                                                 decaying_scolored,
                    kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                                 decaying_EWino,
                                                shortcut->getNeutralinoOne(),
                                                shortcut->getSelectronL(),
                                                shortcut->getSelectronR(),
                                          shortcut->getElectronSneutrinoL(),
                                        shortcut->getElectronSneutrinoR() );
        current_distribution = &direct_lepton;
        direct_electron_pass = integrate_acceptance( primary_cut,
                                                     current_distribution );
        direct_electron_fail
        = ( 1.0 - integrate_acceptance( secondary_cut,
                                        current_distribution ) );

      }

    double direct_muon_pass = 0.0;
    double direct_muon_fail = 0.0;
    if( lhcFaserGlobal::negligibleBr < to_muon_BR )
      // if the branching ratio into this channel is not negligible...
      {

        chargino_three_body_decay direct_lepton( shortcut->getCppSlha(),
                                                 decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                                 decaying_EWino,
                                                shortcut->getNeutralinoOne(),
                                                 shortcut->getSmuonL(),
                                                 shortcut->getSmuonR(),
                                              shortcut->getMuonSneutrinoL(),
                                            shortcut->getMuonSneutrinoR() );
        current_distribution = &direct_lepton;
        direct_muon_pass = integrate_acceptance( primary_cut,
                                                 current_distribution );
        direct_muon_fail = ( 1.0 - integrate_acceptance( secondary_cut,
                                                      current_distribution ) );

      }

    double direct_jet_pass = 0.0;
    if( lhcFaserGlobal::negligibleBr < to_jets_BR )
      // if the branching ratio into this channel is not negligible...
      {

        chargino_three_body_decay direct_jet( shortcut->getCppSlha(),
                                              decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                              decaying_EWino,
                                              shortcut->getNeutralinoOne(),
                                              shortcut->getSdownL(),
                                              shortcut->getSdownR(),
                                              shortcut->getSupL(),
                                              shortcut->getSupR() );
        // actually, here I should use every nontop flavor of squark
        // separately, but it's unlikely to make any significant difference.
        current_distribution = &direct_jet;
        direct_jet_pass = integrate_acceptance( jet_cut,
                                                current_distribution );

      }
    double direct_jet_fail = ( 1.0 - direct_jet_pass );

    double tau_muon_pass = 0.0;
    double tau_muon_fail = 0.0;
    double tau_pion_pass = 0.0;
    if( lhcFaserGlobal::negligibleBr < to_tau_BR )
      // if the branching ratio into this channel is not negligible...
      {

        chargino_three_body_decay direct_lepton( shortcut->getCppSlha(),
                                                 decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                                 decaying_EWino,
                                                shortcut->getNeutralinoOne(),
                                                 shortcut->getStauOne(),
                                                 shortcut->getStauTwo(),
                                               shortcut->getTauSneutrinoL(),
                                             shortcut->getTauSneutrinoR() );
        current_distribution = &direct_lepton;
        visible_tau_decay_product hard_tau_muon( current_distribution,
                                          shortcut->getHardMuonFromTau() );
        visible_tau_decay_product soft_tau_muon( current_distribution,
                                          shortcut->getSoftMuonFromTau() );
        visible_tau_decay_product hard_tau_pion( current_distribution,
                                          shortcut->getHardPionFromTau() );
        visible_tau_decay_product soft_tau_pion( current_distribution,
                                          shortcut->getSoftPionFromTau() );

        /* at the moment, I have no plans to provide a tau left-handedness
         * breakdown from a 3-body decay, so I assume that 50% of them are of
         * each handedness.
         */
        current_distribution = &hard_tau_muon;
        leptonEnergyDistribution* other_distribution = &soft_tau_muon;
        tau_muon_pass = ( 0.5 * integrate_acceptance( primary_cut,
                                                      current_distribution )
                          + 0.5 * integrate_acceptance( primary_cut,
                                                        other_distribution ) );
        tau_muon_fail = ( 0.5 * integrate_acceptance( secondary_cut,
                                                      current_distribution )
                          + 0.5 * integrate_acceptance( secondary_cut,
                                                        other_distribution ) );
        current_distribution = &hard_tau_pion;
        other_distribution = &soft_tau_pion;
        tau_pion_pass = ( 0.5 * integrate_acceptance( jet_cut,
                                                      current_distribution )
                          + 0.5 * integrate_acceptance( jet_cut,
                                                        other_distribution ) );

      }
    double tau_pion_fail = ( 1.0 - tau_pion_pass );

    *two_jets_no_leptons += ( to_jets_BR
                              * ( direct_jet_pass * direct_jet_pass ) );

    double
    BR_to_tau_jet
    = ( to_tau_BR * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR );

    *one_jet_no_leptons += ( 2.0 * to_jets_BR
                             * ( direct_jet_pass * direct_jet_fail )
                             + BR_to_tau_jet * tau_pion_pass );
    // factor of 2.0 accounts for either jet particle or antiparticle passing
    // & the other failing, or the other way around.

    *no_jets_one_muon
    += ( to_muon_BR * direct_muon_pass
         + to_tau_BR * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR
           * tau_muon_pass );

    *no_jets_one_electron
    += ( to_electron_BR * direct_electron_pass
        + to_tau_BR * CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR
          * tau_muon_pass );

    *no_jets_no_leptons
    += ( to_jets_BR * ( direct_jet_fail * direct_jet_fail )
         + BR_to_tau_jet * tau_pion_fail
         + to_muon_BR * direct_muon_fail
         + to_electron_BR * direct_electron_fail
         + to_tau_BR
           * ( CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR
               + CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR )
           * tau_muon_fail );

  }



}  // end of LHC_FASER namespace.
