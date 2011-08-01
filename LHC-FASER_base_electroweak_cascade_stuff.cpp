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
  electroweakCascade::electroweakCascade(
                                         leptonAcceptanceParameterSet* const kinematics,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                                  //bool const coloredDecayerIsNotAntiparticle,
                CppSLHA::particle_property_set const* const electroweakDecayer,
                              //bool const electroweakDecayerIsNotAntiparticle,
               CppSLHA::particle_property_set const* const intermediateDecayer,
                                          bool const canDoOssfMinusOsdf,
                                        input_handler const* const shortcut ) :
    kinematics( kinematics ),
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
  electroweakCascade::getAcceptance( acceptanceCutSet const* const cuts,
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
    if( needs_to_prepare_for_this_point() )
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



  electroweakCascadeSet::electroweakCascadeSet(
                                         leptonAcceptanceParameterSet* const kinematics,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                CppSLHA::particle_property_set const* const electroweakDecayer,
                                        input_handler const* const shortcut ) :
    kinematics( kinematics ),
    coloredDecayer( coloredDecayer ),
    electroweakDecayer( electroweakDecayer ),
    shortcut( shortcut )
  {
    if( CppSLHA::PDG_code::neutralino_one
        == electroweakDecayer->get_PDG_code() )
    {
      currentCascade = new lightestNeutralinoCascade();
      cascades.push_back( currentCascade );
    }
    else if( shortcut->is_in( electroweakDecayer->get_PDG_code(),
                              shortcut->get_unstable_neutralinos() ) )
    {
      currentCascade = new neutralinoToSemuCascade( kinematics,
                                                    coloredDecayer,
                                                    electroweakDecayer,
                                                   shortcut->get_selectron_L(),
                                                    shortcut );
      cascades.push_back( currentCascade );
      currentCascade = new neutralinoToSemuCascade( kinematics,
                                                    coloredDecayer,
                                                    electroweakDecayer,
                                                   shortcut->get_selectron_R(),
                                                    shortcut );
      cascades.push_back( currentCascade );
      currentCascade = new neutralinoToSemuCascade( kinematics,
                                                    coloredDecayer,
                                                    electroweakDecayer,
                                                    shortcut->get_smuon_L(),
                                                    shortcut );
      cascades.push_back( currentCascade );
      currentCascade = new neutralinoToSemuCascade( kinematics,
                                                    coloredDecayer,
                                                    electroweakDecayer,
                                                    shortcut->get_smuon_R(),
                                                    shortcut );
      cascades.push_back( currentCascade );
      currentCascade = new neutralinoToStauCascade( kinematics,
                                                    coloredDecayer,
                                                    electroweakDecayer,
                                                    shortcut->get_stau_one(),
                                                    shortcut );
      cascades.push_back( currentCascade );
      currentCascade = new neutralinoToStauCascade( kinematics,
                                                    coloredDecayer,
                                                    electroweakDecayer,
                                                    shortcut->get_stau_two(),
                                                    shortcut );
      cascades.push_back( currentCascade );
      currentCascade = new neutralinoToZCascade( kinematics,
                                                 coloredDecayer,
                                                 electroweakDecayer,
                                                 shortcut );
      cascades.push_back( currentCascade );
      for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
           vevedBosonIterator
           = shortcut->get_neutral_EWSB_scalars_and_pseudoscalars()->begin();
           shortcut->get_neutral_EWSB_scalars_and_pseudoscalars()->end()
           > vevedBosonIterator;
           ++vevedBosonIterator )
      {
        currentCascade = new neutralinoToHiggsCascade( kinematics,
                                                       coloredDecayer,
                                                       electroweakDecayer,
                                                       *vevedBosonIterator,
                                                       shortcut );
        cascades.push_back( currentCascade );
      }
      currentCascade = new neutralinoThreeBodyDecayCascade( kinematics,
                                                            coloredDecayer,
                                                            electroweakDecayer,
                                                            shortcut );
      cascades.push_back( currentCascade );
      currentCascade = new neutralinoThreeBodyToTausCascade( kinematics,
                                                             coloredDecayer,
                                                            electroweakDecayer,
                                                             shortcut );
      cascades.push_back( currentCascade );
    }
    else if( shortcut->is_in( electroweakDecayer->get_PDG_code(),
                              shortcut->get_charginos() ) )
    {
      currentCascade = new charginoToSemuCascade( kinematics,
                                                  coloredDecayer,
                                               coloredDecayerIsNotAntiparticle,
                                                  electroweakDecayer,
                                           electroweakDecayerIsNotAntiparticle,
                                                  shortcut->get_selectron_L(),
                                                  shortcut );
      cascades.push_back( currentCascade );
      currentCascade = new charginoToSemuCascade( kinematics,
                                                  coloredDecayer,
                                                  electroweakDecayer,
                                                  shortcut->get_selectron_R(),
                                                  shortcut );
      cascades.push_back( currentCascade );
      currentCascade = new charginoToEmuSneutrinoCascade( kinematics,
                                                          coloredDecayer,
                                                          electroweakDecayer,
                                          shortcut->get_electron_sneutrino_L(),
                                                          shortcut );
      cascades.push_back( currentCascade );
      cascades.push_back( currentCascade );
      currentCascade = new neutralinoToStauCascade( kinematics,
                                                    coloredDecayer,
                                                    electroweakDecayer,
                                                    shortcut->get_stau_one(),
                                                    shortcut );
      cascades.push_back( currentCascade );
      currentCascade = new neutralinoToStauCascade( kinematics,
                                                    coloredDecayer,
                                                    electroweakDecayer,
                                                    shortcut->get_stau_two(),
                                                    shortcut );
      cascades.push_back( currentCascade );
      currentCascade = new neutralinoToZCascade( kinematics,
                                                 coloredDecayer,
                                                 electroweakDecayer,
                                                 shortcut );
      cascades.push_back( currentCascade );
      for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
           vevedBosonIterator
           = shortcut->get_neutral_EWSB_scalars_and_pseudoscalars()->begin();
           shortcut->get_neutral_EWSB_scalars_and_pseudoscalars()->end()
           > vevedBosonIterator;
           ++vevedBosonIterator )
      {
        currentCascade = new neutralinoToHiggsCascade( kinematics,
                                                       coloredDecayer,
                                                       electroweakDecayer,
                                                       *vevedBosonIterator,
                                                       shortcut );
        cascades.push_back( currentCascade );
      }
      currentCascade = new neutralinoThreeBodyDecayCascade( kinematics,
                                                            coloredDecayer,
                                                            electroweakDecayer,
                                                            shortcut );
      cascades.push_back( currentCascade );
      currentCascade = new neutralinoThreeBodyToTausCascade( kinematics,
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
                                          kinematics_handler* const kinematics,
                                                        int const beamEnergy,
                                        input_handler const* const shortcut ) :
    kinematicsTable( kinematics->get_lepton_acceptance_table( beamEnergy ) ),
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
      = new electroweakCascadeSet( kinematicsTable->getAcceptanceValues(
                                                                coloredDecayer,
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

    if( LHC_FASER_global::negligible_BR < to_Higgs_BR )
      // if the branching ratio into this channel is not negligible...
      {

        // we start with the direct distributions:
        Higgs_muon_plus_antimuon direct_lepton( shortcut->get_CppSLHA(),
                                                decaying_scolored,
                    leptonKinematics->getAcceptance()->getEffectiveSquarkMass(),
                                                decaying_EWino,
                                                mediating_particle,
                                              shortcut->get_neutralino_one() );
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
        if( LHC_FASER_global::negligible_BR
            < ( to_taus_BR
                * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR ) )
          {

            current_distribution = &direct_lepton;
            // for EWSB scalars, half the tau leptons are 1 chirality, the
            // other half are the other chirality.
            visible_tau_decay_product hard_tau_muon( current_distribution,
                                          shortcut->get_hard_muon_from_tau() );
            visible_tau_decay_product soft_tau_muon( current_distribution,
                                          shortcut->get_soft_muon_from_tau() );
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
        if( LHC_FASER_global::negligible_BR
            < ( to_taus_BR
                * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR ) )
          {

            current_distribution = &direct_lepton;
            visible_tau_decay_product soft_tau_pion( current_distribution,
                                          shortcut->get_soft_pion_from_tau() );
            visible_tau_decay_product hard_tau_pion( current_distribution,
                                          shortcut->get_hard_pion_from_tau() );
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
        *zerJetsOnePositiveElectron += no_jets_one_negative_electron_value;


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
    if( LHC_FASER_global::negligible_BR < to_electrons_BR )
      // if the branching ratio into this channel is not negligible...
      {

        neutralino_three_body_decay direct_lepton( shortcut->get_CppSLHA(),
                                                   decaying_scolored,
                    kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                                   decaying_EWino,
                                                shortcut->get_neutralino_one(),
                                                   shortcut->get_selectron_L(),
                                                 shortcut->get_selectron_R() );
        current_distribution = &direct_lepton;
        direct_electron_pass = integrate_acceptance( primary_cut,
                                                     current_distribution );
        direct_electron_fail
        = ( 1.0 - integrate_acceptance( secondary_cut,
                                        current_distribution ) );

      }

    double direct_muon_pass = 0.0;
    double direct_muon_fail = 0.0;
    if( LHC_FASER_global::negligible_BR < to_muons_BR )
      // if the branching ratio into this channel is not negligible...
      {

        neutralino_three_body_decay direct_lepton( shortcut->get_CppSLHA(),
                                             decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                             decaying_EWino,
                                             shortcut->get_neutralino_one(),
                                             shortcut->get_smuon_L(),
                                             shortcut->get_smuon_R() );
        current_distribution = &direct_lepton;
        direct_muon_pass = integrate_acceptance( primary_cut,
                                                 current_distribution );
        direct_muon_fail = ( 1.0 - integrate_acceptance( secondary_cut,
                                                      current_distribution ) );

      }

    double direct_jet_pass = 0.0;
    if( LHC_FASER_global::negligible_BR < to_jets_BR )
      // if the branching ratio into this channel is not negligible...
      {

        neutralino_three_body_decay direct_jet( shortcut->get_CppSLHA(),
                                                decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                                decaying_EWino,
                                                shortcut->get_neutralino_one(),
                                                shortcut->get_sdown_L(),
                                                shortcut->get_sdown_R() );
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
    if( LHC_FASER_global::negligible_BR < to_taus_BR )
      // if the branching ratio into this channel is not negligible...
      {

        neutralino_three_body_decay direct_lepton( shortcut->get_CppSLHA(),
                                                   decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                                   decaying_EWino,
                                                shortcut->get_neutralino_one(),
                                                   shortcut->get_stau_one(),
                                                   shortcut->get_stau_two() );
        current_distribution = &direct_lepton;
        visible_tau_decay_product hard_tau_muon( current_distribution,
                                          shortcut->get_hard_muon_from_tau() );
        visible_tau_decay_product soft_tau_muon( current_distribution,
                                          shortcut->get_soft_muon_from_tau() );
        visible_tau_decay_product hard_tau_pion( current_distribution,
                                          shortcut->get_hard_pion_from_tau() );
        visible_tau_decay_product soft_tau_pion( current_distribution,
                                          shortcut->get_soft_pion_from_tau() );

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
    exclusive_BR_calculator*
    first_BR
    = shortcut->get_exclusive_BR_calculator( decaying_EWino,
                                             mediating_particle,
                                             not_particle_antiparticle_flip,
                                             shortcut->get_empty_list() );
    exclusive_BR_calculator*
    second_BR
    = shortcut->get_exclusive_BR_calculator( mediating_particle,
                                             shortcut->get_neutralino_one(),
                                             true,
                                             shortcut->get_empty_list() );
    double cascade_BR = ( first_BR->get_BR() * second_BR->get_BR() );

    if( LHC_FASER_global::negligible_BR < cascade_BR )
      // if the branching ratio into this channel is not negligible...
      {

        // we just have a far lepton:
        same_chirality_far_muon same_lepton( shortcut->get_CppSLHA(),
                                             decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                             decaying_EWino,
                                             mediating_particle,
                                             shortcut->get_neutralino_one() );
        opposite_chirality_far_muon opposite_lepton( shortcut->get_CppSLHA(),
                                                     decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                                     decaying_EWino,
                                                     mediating_particle,
                                              shortcut->get_neutralino_one() );
        current_distribution = &same_lepton;

        double
        same_handedness
        = shortcut->quark_or_lepton_left_handedness(
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
    exclusive_BR_calculator*
    first_BR
    = shortcut->get_exclusive_BR_calculator( decaying_EWino,
                                             mediating_particle,
                                             not_particle_antiparticle_flip,
                                             shortcut->get_empty_list() );
    exclusive_BR_calculator*
    second_BR
    = shortcut->get_exclusive_BR_calculator( mediating_particle,
                                             shortcut->get_neutralino_one(),
                                             true,
                                             shortcut->get_empty_list() );
    double cascade_BR = ( first_BR->get_BR() * second_BR->get_BR() );

    if( LHC_FASER_global::negligible_BR < cascade_BR )
      // if the branching ratio into this channel is not negligible...
      {

        // we just have a far lepton:
        same_chirality_far_muon same_tau( shortcut->get_CppSLHA(),
                                          decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                          decaying_EWino,
                                          mediating_particle,
                                          shortcut->get_neutralino_one() );
        opposite_chirality_far_muon opposite_tau( shortcut->get_CppSLHA(),
                                                  decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                                  decaying_EWino,
                                                  mediating_particle,
                                              shortcut->get_neutralino_one() );

        double
        quark_left_handedness
        = shortcut->quark_or_lepton_left_handedness(
                                             decaying_scolored->get_PDG_code(),
                                              decaying_EWino->get_PDG_code() );
        double quark_right_handedness = ( 1.0 - quark_left_handedness );
        double
        tau_left_handedness
        = shortcut->quark_or_lepton_left_handedness(
                                            mediating_particle->get_PDG_code(),
                                              decaying_EWino->get_PDG_code() );
        double tau_right_handedness = ( 1.0 - tau_left_handedness );

        double muon_pass = 0.0;
        double muon_fail = 0.0;
        double pion_pass = 0.0;
        double pion_fail = 0.0;
        double pion_acceptance;
        if( LHC_FASER_global::negligible_BR
            < ( cascade_BR
                * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR ) )
          {

            current_distribution = &same_tau;
            visible_tau_decay_product hard_same_muon( current_distribution,
                                          shortcut->get_hard_muon_from_tau() );
            leptonEnergyDistribution* visible_distribution = &hard_same_muon;
            muon_pass += ( quark_left_handedness * tau_left_handedness
                           * integrate_acceptance( primary_cut,
                                                   visible_distribution ) );
            muon_fail
            += ( quark_left_handedness * tau_left_handedness
                 * ( 1.0 - integrate_acceptance( secondary_cut,
                                                 visible_distribution ) ) );
            visible_tau_decay_product soft_same_muon( current_distribution,
                                          shortcut->get_soft_muon_from_tau() );
            visible_distribution = &soft_same_muon;
            muon_pass += ( quark_right_handedness * tau_right_handedness
                           * integrate_acceptance( primary_cut,
                                                   visible_distribution ) );
            muon_fail
            += ( quark_right_handedness * tau_right_handedness
                 * ( 1.0 - integrate_acceptance( secondary_cut,
                                                 visible_distribution ) ) );
            visible_tau_decay_product soft_same_pion( current_distribution,
                                          shortcut->get_soft_pion_from_tau() );
            visible_distribution = &soft_same_pion;
            pion_acceptance = integrate_acceptance( jet_cut,
                                                    visible_distribution );
            pion_pass += ( quark_left_handedness * tau_left_handedness
                           * pion_acceptance );
            pion_fail
            += ( quark_left_handedness * tau_left_handedness
                 * ( 1.0 - pion_acceptance ) );
            visible_tau_decay_product hard_same_pion( current_distribution,
                                          shortcut->get_hard_pion_from_tau() );
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
                                          shortcut->get_hard_muon_from_tau() );
            visible_distribution = &hard_opposite_muon;
            muon_pass += ( quark_right_handedness * tau_left_handedness
                           * integrate_acceptance( primary_cut,
                                                   visible_distribution ) );
            muon_fail
            += ( quark_right_handedness * tau_left_handedness
                 * ( 1.0 - integrate_acceptance( secondary_cut,
                                                 visible_distribution ) ) );
            visible_tau_decay_product soft_opposite_muon( current_distribution,
                                          shortcut->get_soft_muon_from_tau() );
            visible_distribution = &soft_opposite_muon;
            muon_pass += ( quark_left_handedness * tau_right_handedness
                           * integrate_acceptance( primary_cut,
                                                   visible_distribution ) );
            muon_fail
            += ( quark_left_handedness * tau_right_handedness
                 * ( 1.0 - integrate_acceptance( secondary_cut,
                                                 visible_distribution ) ) );
            visible_tau_decay_product soft_opposite_pion( current_distribution,
                                          shortcut->get_soft_pion_from_tau() );
            visible_distribution = &soft_same_pion;
            pion_acceptance = integrate_acceptance( jet_cut,
                                                    visible_distribution );
            pion_pass += ( quark_right_handedness * tau_left_handedness
                           * pion_acceptance );
            pion_fail
            += ( quark_right_handedness * tau_left_handedness
                 * ( 1.0 - pion_acceptance ) );
            visible_tau_decay_product hard_opposite_pion( current_distribution,
                                          shortcut->get_hard_pion_from_tau() );
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
    exclusive_BR_calculator*
    first_BR
    = shortcut->get_exclusive_BR_calculator( decaying_EWino,
                                             mediating_particle,
                                             not_particle_antiparticle_flip,
                                             shortcut->get_empty_list() );
    exclusive_BR_calculator*
    second_BR
    = shortcut->get_exclusive_BR_calculator( mediating_particle,
                                             shortcut->get_neutralino_one(),
                                             true,
                                             shortcut->get_empty_list() );
    double cascade_BR = ( first_BR->get_BR() * second_BR->get_BR() );

    if( LHC_FASER_global::negligible_BR < cascade_BR )
      // if the branching ratio into this channel is not negligible...
      {

        // we just have a near lepton:
        same_chirality_near_muon same_lepton( shortcut->get_CppSLHA(),
                                              decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                              decaying_EWino,
                                              mediating_particle );
        opposite_chirality_near_muon opposite_lepton( shortcut->get_CppSLHA(),
                                                      decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                                      decaying_EWino,
                                                      mediating_particle );
        current_distribution = &same_lepton;

        double
        quark_left_handedness = shortcut->quark_or_lepton_left_handedness(
                                             decaying_scolored->get_PDG_code(),
                                              decaying_EWino->get_PDG_code() );
        double
        lepton_left_handedness = shortcut->quark_or_lepton_left_handedness(
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
    exclusive_BR_calculator*
    first_BR
    = shortcut->get_exclusive_BR_calculator( decaying_EWino,
                                             mediating_particle,
                                             not_particle_antiparticle_flip,
                                             shortcut->get_empty_list() );
    exclusive_BR_calculator*
    second_BR
    = shortcut->get_exclusive_BR_calculator( mediating_particle,
                                             shortcut->get_neutralino_one(),
                                             true,
                                             shortcut->get_empty_list() );
    double cascade_BR = ( first_BR->get_BR() * second_BR->get_BR() );

    if( LHC_FASER_global::negligible_BR < cascade_BR )
      // if the branching ratio into this channel is not negligible...
      {

        // we just have a near lepton:
        same_chirality_near_muon same_tau( shortcut->get_CppSLHA(),
                                           decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                           decaying_EWino,
                                           mediating_particle );
        opposite_chirality_near_muon opposite_tau( shortcut->get_CppSLHA(),
                                                   decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                                   decaying_EWino,
                                                   mediating_particle );

        double
        quark_left_handedness
        = shortcut->quark_or_lepton_left_handedness(
                                             decaying_scolored->get_PDG_code(),
                                              decaying_EWino->get_PDG_code() );
        double quark_right_handedness = ( 1.0 - quark_left_handedness );
        double
        tau_left_handedness
        = shortcut->quark_or_lepton_left_handedness(
                                            mediating_particle->get_PDG_code(),
                                              decaying_EWino->get_PDG_code() );
        double tau_right_handedness = ( 1.0 - tau_left_handedness );

        double muon_pass = 0.0;
        double muon_fail = 0.0;
        double pion_pass = 0.0;
        double pion_fail = 0.0;
        double pion_acceptance;
        if( LHC_FASER_global::negligible_BR
            < ( cascade_BR
                * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR ) )
          {

            current_distribution = &same_tau;
            visible_tau_decay_product hard_same_muon( current_distribution,
                                          shortcut->get_hard_muon_from_tau() );
            leptonEnergyDistribution* visible_distribution = &hard_same_muon;
            muon_pass += ( quark_left_handedness * tau_left_handedness
                           * integrate_acceptance( primary_cut,
                                                   visible_distribution ) );
            muon_fail
            += ( quark_left_handedness * tau_left_handedness
                 * ( 1.0 - integrate_acceptance( secondary_cut,
                                                 visible_distribution ) ) );
            visible_tau_decay_product soft_same_muon( current_distribution,
                                          shortcut->get_soft_muon_from_tau() );
            visible_distribution = &soft_same_muon;
            muon_pass += ( quark_right_handedness * tau_right_handedness
                           * integrate_acceptance( primary_cut,
                                                   visible_distribution ) );
            muon_fail
            += ( quark_right_handedness * tau_right_handedness
                 * ( 1.0 - integrate_acceptance( secondary_cut,
                                                 visible_distribution ) ) );
            visible_tau_decay_product soft_same_pion( current_distribution,
                                          shortcut->get_soft_pion_from_tau() );
            visible_distribution = &soft_same_pion;
            pion_acceptance = integrate_acceptance( jet_cut,
                                                    visible_distribution );
            pion_pass += ( quark_left_handedness * tau_left_handedness
                           * pion_acceptance );
            pion_fail
            += ( quark_left_handedness * tau_left_handedness
                 * ( 1.0 - pion_acceptance ) );
            visible_tau_decay_product hard_same_pion( current_distribution,
                                          shortcut->get_hard_pion_from_tau() );
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
                                          shortcut->get_hard_muon_from_tau() );
            visible_distribution = &hard_opposite_muon;
            muon_pass += ( quark_right_handedness * tau_left_handedness
                           * integrate_acceptance( primary_cut,
                                                   visible_distribution ) );
            muon_fail
            += ( quark_right_handedness * tau_left_handedness
                 * ( 1.0 - integrate_acceptance( secondary_cut,
                                                 visible_distribution ) ) );
            visible_tau_decay_product soft_opposite_muon( current_distribution,
                                          shortcut->get_soft_muon_from_tau() );
            visible_distribution = &soft_opposite_muon;
            muon_pass += ( quark_left_handedness * tau_right_handedness
                           * integrate_acceptance( primary_cut,
                                                   visible_distribution ) );
            muon_fail
            += ( quark_left_handedness * tau_right_handedness
                 * ( 1.0 - integrate_acceptance( secondary_cut,
                                                 visible_distribution ) ) );
            visible_tau_decay_product soft_opposite_pion( current_distribution,
                                          shortcut->get_soft_pion_from_tau() );
            visible_distribution = &soft_same_pion;
            pion_acceptance = integrate_acceptance( jet_cut,
                                                    visible_distribution );
            pion_pass += ( quark_right_handedness * tau_left_handedness
                           * pion_acceptance );
            pion_fail
            += ( quark_right_handedness * tau_left_handedness
                 * ( 1.0 - pion_acceptance ) );
            visible_tau_decay_product hard_opposite_pion( current_distribution,
                                          shortcut->get_hard_pion_from_tau() );
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

    if( LHC_FASER_global::negligible_BR < to_W_BR )
      // if the branching ratio into this channel is not negligible...
      {

        // we start with the direct distributions:
        W_minus_direct_jet direct_jet( shortcut->get_CppSLHA(),
                                       decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                       decaying_EWino,
                                       mediating_particle,
                                       shortcut->get_neutralino_one() );
        // debugging:
        /**std::cout
        << std::endl
        << "&direct_jet = " << &direct_jet;
        std::cout << std::endl;**/
        W_minus_plus_up_L_type_muon L_direct_lepton( shortcut->get_CppSLHA(),
                                                     decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                                     decaying_EWino,
                                                     mediating_particle,
                                              shortcut->get_neutralino_one() );
        // debugging:
        /**std::cout
        << std::endl
        << "&L_direct_lepton = " << &L_direct_lepton;
        std::cout << std::endl;**/
        W_minus_plus_up_R_type_muon R_direct_lepton( shortcut->get_CppSLHA(),
                                                     decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                                     decaying_EWino,
                                                     mediating_particle,
                                              shortcut->get_neutralino_one() );
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
        = shortcut->quark_or_lepton_left_handedness(
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
        if( LHC_FASER_global::negligible_BR
            < ( to_tau_BR
                * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR ) )
          {

            current_distribution = &L_direct_lepton;
            visible_tau_decay_product L_tau_muon( current_distribution,
                                          shortcut->get_hard_muon_from_tau() );
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
                                          shortcut->get_hard_muon_from_tau() );
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
        if( LHC_FASER_global::negligible_BR
            < ( to_tau_BR
                * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR ) )
          {

            current_distribution = &L_direct_lepton;
            visible_tau_decay_product L_tau_pion( current_distribution,
                                          shortcut->get_soft_pion_from_tau() );
            current_distribution = &L_tau_pion;
            tau_pion_pass
            += ( quark_left_handedness * integrate_acceptance( jet_cut,
                                                      current_distribution ) );
            current_distribution = &R_direct_lepton;
            visible_tau_decay_product R_tau_pion( current_distribution,
                                          shortcut->get_soft_pion_from_tau() );
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

    if( LHC_FASER_global::negligible_BR < to_Higgs_BR )
      // if the branching ratio into this channel is not negligible...
      {

        // we start with the direct distributions:
        negatively_charged_Higgs_jet direct_jet( shortcut->get_CppSLHA(),
                                                 decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                                 decaying_EWino,
                                                 mediating_particle,
                                              shortcut->get_neutralino_one() );
        negatively_charged_Higgs_muon direct_lepton( shortcut->get_CppSLHA(),
                                                     decaying_scolored,
                    kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                                     decaying_EWino,
                                                     mediating_particle,
                                              shortcut->get_neutralino_one() );
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
        if( LHC_FASER_global::negligible_BR
            < ( to_tau_BR
                * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR ) )
          {

            current_distribution = &direct_lepton;
            visible_tau_decay_product tau_muon( current_distribution,
                                          shortcut->get_hard_muon_from_tau() );
            current_distribution = &tau_muon;
            tau_muon_pass
            += integrate_acceptance( primary_cut,
                                     current_distribution );
            tau_muon_fail
            += ( 1.0 - integrate_acceptance( secondary_cut,
                                             current_distribution ) );

          }
        double tau_pion_pass = 0.0;
        if( LHC_FASER_global::negligible_BR
            < ( to_tau_BR
                * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR ) )
          {

            current_distribution = &direct_lepton;
            visible_tau_decay_product tau_pion( current_distribution,
                                          shortcut->get_soft_pion_from_tau() );
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
    if( LHC_FASER_global::negligible_BR < to_electron_BR )
      // if the branching ratio into this channel is not negligible...
      {

        chargino_three_body_decay direct_lepton( shortcut->get_CppSLHA(),
                                                 decaying_scolored,
                    kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                                 decaying_EWino,
                                                shortcut->get_neutralino_one(),
                                                shortcut->get_selectron_L(),
                                                shortcut->get_selectron_R(),
                                          shortcut->get_electron_sneutrino_L(),
                                        shortcut->get_electron_sneutrino_R() );
        current_distribution = &direct_lepton;
        direct_electron_pass = integrate_acceptance( primary_cut,
                                                     current_distribution );
        direct_electron_fail
        = ( 1.0 - integrate_acceptance( secondary_cut,
                                        current_distribution ) );

      }

    double direct_muon_pass = 0.0;
    double direct_muon_fail = 0.0;
    if( LHC_FASER_global::negligible_BR < to_muon_BR )
      // if the branching ratio into this channel is not negligible...
      {

        chargino_three_body_decay direct_lepton( shortcut->get_CppSLHA(),
                                                 decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                                 decaying_EWino,
                                                shortcut->get_neutralino_one(),
                                                 shortcut->get_smuon_L(),
                                                 shortcut->get_smuon_R(),
                                              shortcut->get_muon_sneutrino_L(),
                                            shortcut->get_muon_sneutrino_R() );
        current_distribution = &direct_lepton;
        direct_muon_pass = integrate_acceptance( primary_cut,
                                                 current_distribution );
        direct_muon_fail = ( 1.0 - integrate_acceptance( secondary_cut,
                                                      current_distribution ) );

      }

    double direct_jet_pass = 0.0;
    if( LHC_FASER_global::negligible_BR < to_jets_BR )
      // if the branching ratio into this channel is not negligible...
      {

        chargino_three_body_decay direct_jet( shortcut->get_CppSLHA(),
                                              decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                              decaying_EWino,
                                              shortcut->get_neutralino_one(),
                                              shortcut->get_sdown_L(),
                                              shortcut->get_sdown_R(),
                                              shortcut->get_sup_L(),
                                              shortcut->get_sup_R() );
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
    if( LHC_FASER_global::negligible_BR < to_tau_BR )
      // if the branching ratio into this channel is not negligible...
      {

        chargino_three_body_decay direct_lepton( shortcut->get_CppSLHA(),
                                                 decaying_scolored,
                     kinematics->get_acceptance()->getEffectiveSquarkMass(),
                                                 decaying_EWino,
                                                shortcut->get_neutralino_one(),
                                                 shortcut->get_stau_one(),
                                                 shortcut->get_stau_two(),
                                               shortcut->get_tau_sneutrino_L(),
                                             shortcut->get_tau_sneutrino_R() );
        current_distribution = &direct_lepton;
        visible_tau_decay_product hard_tau_muon( current_distribution,
                                          shortcut->get_hard_muon_from_tau() );
        visible_tau_decay_product soft_tau_muon( current_distribution,
                                          shortcut->get_soft_muon_from_tau() );
        visible_tau_decay_product hard_tau_pion( current_distribution,
                                          shortcut->get_hard_pion_from_tau() );
        visible_tau_decay_product soft_tau_pion( current_distribution,
                                          shortcut->get_soft_pion_from_tau() );

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



  cascade_acceptance_value::cascade_acceptance_value(
                                 readier_for_new_point* const given_readier ) :
    readied_for_new_point( given_readier ),
    acceptance_value( CppSLHA::CppSLHA_global::really_wrong_value )
  {

    // just an initialization list.

  }

  cascade_acceptance_value::~cascade_acceptance_value()
  {

    // does nothing.

  }



  charge_summed_cascade_calculator::charge_summed_cascade_calculator(
                                                double const given_primary_cut,
                                              double const given_secondary_cut,
                                                    double const given_jet_cut,
                               lepton_acceptance_value* const given_kinematics,
           CppSLHA::particle_property_set const* const given_decaying_scolored,
                                 bool const given_scolored_is_not_antiparticle,
              CppSLHA::particle_property_set const* const given_decaying_EWino,
                                    bool const given_EWino_is_not_antiparticle,
                                 input_handler const* const given_shortcuts ) :
    primary_cut( given_primary_cut ),
    secondary_cut( given_secondary_cut ),
    jet_cut( given_jet_cut ),
    kinematics( given_kinematics ),
    decaying_scolored( given_decaying_scolored ),
    scolored_is_not_antiparticle( given_scolored_is_not_antiparticle ),
    decaying_EWino( given_decaying_EWino ),
    EWino_is_not_antiparticle( given_EWino_is_not_antiparticle ),
    OSSF_minus_OSDF_leptons( 0.0 ),
    two_jets_no_leptons( 0.0 ),
    one_jet_one_muon( 0.0 ),
    one_jet_one_electron( 0.0 ),
    no_jets_two_muons( 0.0 ),
    no_jets_two_electrons( 0.0 ),
    no_jets_one_muon_one_electron( 0.0 ),
    one_jet_no_leptons( 0.0 ),
    no_jets_one_muon( 0.0 ),
    no_jets_one_electron( 0.0 ),
    no_jets_no_leptons( 0.0 )
  {

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: charge_summed_cascade_calculator::"
    << "charge_summed_cascade_calculator( "
    << given_primary_cut << ", "
    << given_secondary_cut << ", "
    << given_jet_cut << ", ..."
    << *(given_decaying_scolored->get_name_or_antiname(
                                          given_scolored_is_not_antiparticle ))
    << ", "
    << *(given_decaying_EWino->get_name_or_antiname(
                                             given_EWino_is_not_antiparticle ))
    << ", ... ) called.";
    std::cout << std::endl;
    std::cout
    << std::endl
    << "given_kinematics->get_acceptance()->get_effective_squark_mass() = "
    << given_kinematics->get_acceptance()->get_effective_squark_mass();
    std::cout << std::endl;**/

    if( given_decaying_EWino == given_shortcuts->get_neutralino_one() )
      // if the electroweakino is the lightest neutralino, & hence assumed
      // stable, we set up the appropriate calculator...
      {

        // debugging:
        /**std::cout
        << std::endl
        << "adding neutralino_1 channel";
        std::cout << std::endl;**/

        channel_calculator*
        channel_adder
        = new lightestNeutralinoCascade( no_jets_no_leptons.get_value_pointer() );
        no_jets_no_leptons.add_channel( channel_adder );
        channels.push_back( channel_adder );

      }
    else if( given_kinematics->get_acceptance()->getEffectiveSquarkMass()
             > given_decaying_EWino->get_absolute_mass() )
      // if this decay can happen within our formalism...
      {

        // debugging:
        /**std::cout
        << std::endl
        << given_kinematics->get_acceptance()->get_effective_squark_mass()
        << " > " << given_decaying_EWino->get_absolute_mass();
        std::cout << std::endl;**/

        channel_calculator* channel_adder;

        if( given_shortcuts->is_in( given_decaying_EWino->get_PDG_code(),
                                given_shortcuts->get_unstable_neutralinos() ) )
          {

            if( given_decaying_EWino->get_absolute_mass()
                > given_shortcuts->get_selectron_L()->get_absolute_mass() )
              {

                // debugging:
                /**std::cout
                << std::endl
                << "adding selectron_L channel";
                std::cout << std::endl;**/

                channel_adder
                = new charge_summed_neutralino_to_semu( given_primary_cut,
                                                        given_secondary_cut,
                                                        given_kinematics,
                                                       given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                        given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                            given_shortcuts->get_selectron_L(),
                                                        given_shortcuts,
                                   OSSF_minus_OSDF_leptons.get_value_pointer(),
                                     no_jets_two_electrons.get_value_pointer(),
                                      no_jets_one_electron.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
                OSSF_minus_OSDF_leptons.add_channel( channel_adder );
                no_jets_two_electrons.add_channel( channel_adder );
                no_jets_one_electron.add_channel( channel_adder );
                no_jets_no_leptons.add_channel( channel_adder );
                channels.push_back( channel_adder );

              }
            if( given_decaying_EWino->get_absolute_mass()
                > given_shortcuts->get_selectron_R()->get_absolute_mass() )
              {

                // debugging:
                /**std::cout
                << std::endl
                << "adding selectron_R channel";
                std::cout << std::endl;**/

                channel_adder
                = new charge_summed_neutralino_to_semu( given_primary_cut,
                                                        given_secondary_cut,
                                                        given_kinematics,
                                                       given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                        given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                            given_shortcuts->get_selectron_R(),
                                                        given_shortcuts,
                                   OSSF_minus_OSDF_leptons.get_value_pointer(),
                                     no_jets_two_electrons.get_value_pointer(),
                                      no_jets_one_electron.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
                OSSF_minus_OSDF_leptons.add_channel( channel_adder );
                no_jets_two_electrons.add_channel( channel_adder );
                no_jets_one_electron.add_channel( channel_adder );
                no_jets_no_leptons.add_channel( channel_adder );
                channels.push_back( channel_adder );

              }
            if( given_decaying_EWino->get_absolute_mass()
                > given_shortcuts->get_smuon_L()->get_absolute_mass() )
              {

                // debugging:
                /**std::cout
                << std::endl
                << "adding smuon_L channel";
                std::cout << std::endl;**/

                channel_adder
                = new charge_summed_neutralino_to_semu( given_primary_cut,
                                                        given_secondary_cut,
                                                        given_kinematics,
                                                       given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                        given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                                given_shortcuts->get_smuon_L(),
                                                        given_shortcuts,
                                   OSSF_minus_OSDF_leptons.get_value_pointer(),
                                         no_jets_two_muons.get_value_pointer(),
                                          no_jets_one_muon.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
                OSSF_minus_OSDF_leptons.add_channel( channel_adder );
                no_jets_two_muons.add_channel( channel_adder );
                no_jets_one_muon.add_channel( channel_adder );
                no_jets_no_leptons.add_channel( channel_adder );
                channels.push_back( channel_adder );

              }
            if( given_decaying_EWino->get_absolute_mass()
                > given_shortcuts->get_smuon_R()->get_absolute_mass() )
              {

                // debugging:
                /**std::cout
                << std::endl
                << "adding smuon_R channel";
                std::cout << std::endl;**/

                channel_adder
                = new charge_summed_neutralino_to_semu( given_primary_cut,
                                                        given_secondary_cut,
                                                        given_kinematics,
                                                       given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                        given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                                given_shortcuts->get_smuon_R(),
                                                        given_shortcuts,
                                   OSSF_minus_OSDF_leptons.get_value_pointer(),
                                         no_jets_two_muons.get_value_pointer(),
                                          no_jets_one_muon.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
                OSSF_minus_OSDF_leptons.add_channel( channel_adder );
                no_jets_two_muons.add_channel( channel_adder );
                no_jets_one_muon.add_channel( channel_adder );
                no_jets_no_leptons.add_channel( channel_adder );
                channels.push_back( channel_adder );

              }
            if( given_decaying_EWino->get_absolute_mass()
                > given_shortcuts->get_stau_one()->get_absolute_mass() )
              {

                // debugging:
                /**std::cout
                << std::endl
                << "adding stau_1 channel";
                std::cout << std::endl;**/

                channel_adder
                = new charge_summed_neutralino_to_stau( given_primary_cut,
                                                        given_secondary_cut,
                                                        given_jet_cut,
                                                        given_kinematics,
                                                       given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                        given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                               given_shortcuts->get_stau_one(),
                                                        given_shortcuts,
                                       two_jets_no_leptons.get_value_pointer(),
                                          one_jet_one_muon.get_value_pointer(),
                                      one_jet_one_electron.get_value_pointer(),
                                         no_jets_two_muons.get_value_pointer(),
                                     no_jets_two_electrons.get_value_pointer(),
                             no_jets_one_muon_one_electron.get_value_pointer(),
                                        one_jet_no_leptons.get_value_pointer(),
                                          no_jets_one_muon.get_value_pointer(),
                                      no_jets_one_electron.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
                two_jets_no_leptons.add_channel( channel_adder );
                one_jet_one_muon.add_channel( channel_adder );
                one_jet_one_electron.add_channel( channel_adder );
                no_jets_two_muons.add_channel( channel_adder );
                no_jets_two_electrons.add_channel( channel_adder );
                no_jets_one_muon_one_electron.add_channel( channel_adder );
                one_jet_no_leptons.add_channel( channel_adder );
                no_jets_one_muon.add_channel( channel_adder );
                no_jets_one_electron.add_channel( channel_adder );
                no_jets_no_leptons.add_channel( channel_adder );
                channels.push_back( channel_adder );

              }
            if( given_decaying_EWino->get_absolute_mass()
                > given_shortcuts->get_stau_two()->get_absolute_mass() )
              {

                // debugging:
                /**std::cout
                << std::endl
                << "adding stau_2 channel";
                std::cout << std::endl;**/

                channel_adder
                = new charge_summed_neutralino_to_stau( given_primary_cut,
                                                        given_secondary_cut,
                                                        given_jet_cut,
                                                        given_kinematics,
                                                       given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                        given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                               given_shortcuts->get_stau_two(),
                                                        given_shortcuts,
                                       two_jets_no_leptons.get_value_pointer(),
                                          one_jet_one_muon.get_value_pointer(),
                                      one_jet_one_electron.get_value_pointer(),
                                         no_jets_two_muons.get_value_pointer(),
                                     no_jets_two_electrons.get_value_pointer(),
                             no_jets_one_muon_one_electron.get_value_pointer(),
                                        one_jet_no_leptons.get_value_pointer(),
                                          no_jets_one_muon.get_value_pointer(),
                                      no_jets_one_electron.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
                 two_jets_no_leptons.add_channel( channel_adder );
                 one_jet_one_muon.add_channel( channel_adder );
                 one_jet_one_electron.add_channel( channel_adder );
                 no_jets_two_muons.add_channel( channel_adder );
                 no_jets_two_electrons.add_channel( channel_adder );
                 no_jets_one_muon_one_electron.add_channel( channel_adder );
                 one_jet_no_leptons.add_channel( channel_adder );
                 no_jets_one_muon.add_channel( channel_adder );
                 no_jets_one_electron.add_channel( channel_adder );
                 no_jets_no_leptons.add_channel( channel_adder );
                 channels.push_back( channel_adder );

              }
            if( given_decaying_EWino->get_absolute_mass()
                > ( given_shortcuts->get_Z()->get_absolute_mass()
               + given_shortcuts->get_neutralino_one()->get_absolute_mass() ) )
              {

                // debugging:
                /**std::cout
                << std::endl
                << "adding Z channel";
                std::cout << std::endl;**/

                channel_adder
                = new charge_summed_neutralino_to_Z( given_primary_cut,
                                                     given_secondary_cut,
                                                     given_jet_cut,
                                                     given_kinematics,
                                                     given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                     given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                                     given_shortcuts->get_Z(),
                                                     given_shortcuts,
                                       two_jets_no_leptons.get_value_pointer(),
                                          one_jet_one_muon.get_value_pointer(),
                                      one_jet_one_electron.get_value_pointer(),
                                   OSSF_minus_OSDF_leptons.get_value_pointer(),
                                         no_jets_two_muons.get_value_pointer(),
                                     no_jets_two_electrons.get_value_pointer(),
                             no_jets_one_muon_one_electron.get_value_pointer(),
                                        one_jet_no_leptons.get_value_pointer(),
                                          no_jets_one_muon.get_value_pointer(),
                                      no_jets_one_electron.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
                two_jets_no_leptons.add_channel( channel_adder );
                one_jet_one_muon.add_channel( channel_adder );
                one_jet_one_electron.add_channel( channel_adder );
                OSSF_minus_OSDF_leptons.add_channel( channel_adder );
                no_jets_two_muons.add_channel( channel_adder );
                no_jets_two_electrons.add_channel( channel_adder );
                no_jets_one_muon_one_electron.add_channel( channel_adder );
                one_jet_no_leptons.add_channel( channel_adder );
                no_jets_one_muon.add_channel( channel_adder );
                no_jets_one_electron.add_channel( channel_adder );
                no_jets_no_leptons.add_channel( channel_adder );
                channels.push_back( channel_adder );

              }
            if( given_decaying_EWino->get_absolute_mass()
                > ( given_shortcuts->get_light_neutral_EWSB_scalar(
                                                         )->get_absolute_mass()
               + given_shortcuts->get_neutralino_one()->get_absolute_mass() ) )
              {

                // debugging:
                /**std::cout
                << std::endl
                << "adding h channel";
                std::cout << std::endl;**/

                channel_adder
                = new neutralino_to_Higgs( given_primary_cut,
                                           given_secondary_cut,
                                           given_jet_cut,
                                           given_kinematics,
                                           given_decaying_scolored,
                                           given_scolored_is_not_antiparticle,
                                           given_decaying_EWino,
                                           given_EWino_is_not_antiparticle,
                              given_shortcuts->get_light_neutral_EWSB_scalar(),
                                           given_shortcuts,
                                       two_jets_no_leptons.get_value_pointer(),
                                          one_jet_one_muon.get_value_pointer(),
                                      one_jet_one_electron.get_value_pointer(),
                                          one_jet_one_muon.get_value_pointer(),
                                      one_jet_one_electron.get_value_pointer(),
                                   OSSF_minus_OSDF_leptons.get_value_pointer(),
                                         no_jets_two_muons.get_value_pointer(),
                                     no_jets_two_electrons.get_value_pointer(),
                             no_jets_one_muon_one_electron.get_value_pointer(),
                             no_jets_one_muon_one_electron.get_value_pointer(),
                                        one_jet_no_leptons.get_value_pointer(),
                                          no_jets_one_muon.get_value_pointer(),
                                      no_jets_one_electron.get_value_pointer(),
                                          no_jets_one_muon.get_value_pointer(),
                                      no_jets_one_electron.get_value_pointer(),
                                       no_jets_no_leptons.get_value_pointer() );
                two_jets_no_leptons.add_channel( channel_adder );
                one_jet_one_muon.add_channel( channel_adder );
                one_jet_one_electron.add_channel( channel_adder );
                OSSF_minus_OSDF_leptons.add_channel( channel_adder );
                no_jets_two_muons.add_channel( channel_adder );
                no_jets_two_electrons.add_channel( channel_adder );
                no_jets_one_muon_one_electron.add_channel( channel_adder );
                one_jet_no_leptons.add_channel( channel_adder );
                no_jets_one_muon.add_channel( channel_adder );
                no_jets_one_electron.add_channel( channel_adder );
                no_jets_no_leptons.add_channel( channel_adder );
                channels.push_back( channel_adder );

              }
            if( given_decaying_EWino->get_absolute_mass()
                > ( given_shortcuts->get_heavy_neutral_EWSB_scalar(
                                                         )->get_absolute_mass()
               + given_shortcuts->get_neutralino_one()->get_absolute_mass() ) )
              {

                // debugging:
                /**std::cout
                << std::endl
                << "adding H channel";
                std::cout << std::endl;**/

                channel_adder
                = new neutralino_to_Higgs( given_primary_cut,
                                           given_secondary_cut,
                                           given_jet_cut,
                                           given_kinematics,
                                           given_decaying_scolored,
                                           given_scolored_is_not_antiparticle,
                                           given_decaying_EWino,
                                           given_EWino_is_not_antiparticle,
                              given_shortcuts->get_heavy_neutral_EWSB_scalar(),
                                           given_shortcuts,
                                       two_jets_no_leptons.get_value_pointer(),
                                          one_jet_one_muon.get_value_pointer(),
                                      one_jet_one_electron.get_value_pointer(),
                                          one_jet_one_muon.get_value_pointer(),
                                      one_jet_one_electron.get_value_pointer(),
                                   OSSF_minus_OSDF_leptons.get_value_pointer(),
                                         no_jets_two_muons.get_value_pointer(),
                                     no_jets_two_electrons.get_value_pointer(),
                             no_jets_one_muon_one_electron.get_value_pointer(),
                             no_jets_one_muon_one_electron.get_value_pointer(),
                                        one_jet_no_leptons.get_value_pointer(),
                                          no_jets_one_muon.get_value_pointer(),
                                      no_jets_one_electron.get_value_pointer(),
                                          no_jets_one_muon.get_value_pointer(),
                                      no_jets_one_electron.get_value_pointer(),
                                       no_jets_no_leptons.get_value_pointer() );
                two_jets_no_leptons.add_channel( channel_adder );
                one_jet_one_muon.add_channel( channel_adder );
                one_jet_one_electron.add_channel( channel_adder );
                OSSF_minus_OSDF_leptons.add_channel( channel_adder );
                no_jets_two_muons.add_channel( channel_adder );
                no_jets_two_electrons.add_channel( channel_adder );
                no_jets_one_muon_one_electron.add_channel( channel_adder );
                one_jet_no_leptons.add_channel( channel_adder );
                no_jets_one_muon.add_channel( channel_adder );
                no_jets_one_electron.add_channel( channel_adder );
                no_jets_no_leptons.add_channel( channel_adder );
                channels.push_back( channel_adder );

              }
            if( given_decaying_EWino->get_absolute_mass()
                > ( given_shortcuts->get_neutral_EWSB_pseudoscalar(
                                                         )->get_absolute_mass()
               + given_shortcuts->get_neutralino_one()->get_absolute_mass() ) )
              {

                // debugging:
                /**std::cout
                << std::endl
                << "adding A channel";
                std::cout << std::endl;**/

                channel_adder
                = new neutralino_to_Higgs( given_primary_cut,
                                           given_secondary_cut,
                                           given_jet_cut,
                                           given_kinematics,
                                           given_decaying_scolored,
                                           given_scolored_is_not_antiparticle,
                                           given_decaying_EWino,
                                           given_EWino_is_not_antiparticle,
                              given_shortcuts->get_neutral_EWSB_pseudoscalar(),
                                           given_shortcuts,
                                       two_jets_no_leptons.get_value_pointer(),
                                          one_jet_one_muon.get_value_pointer(),
                                      one_jet_one_electron.get_value_pointer(),
                                          one_jet_one_muon.get_value_pointer(),
                                      one_jet_one_electron.get_value_pointer(),
                                   OSSF_minus_OSDF_leptons.get_value_pointer(),
                                         no_jets_two_muons.get_value_pointer(),
                                     no_jets_two_electrons.get_value_pointer(),
                             no_jets_one_muon_one_electron.get_value_pointer(),
                             no_jets_one_muon_one_electron.get_value_pointer(),
                                        one_jet_no_leptons.get_value_pointer(),
                                          no_jets_one_muon.get_value_pointer(),
                                      no_jets_one_electron.get_value_pointer(),
                                          no_jets_one_muon.get_value_pointer(),
                                      no_jets_one_electron.get_value_pointer(),
                                       no_jets_no_leptons.get_value_pointer() );
                two_jets_no_leptons.add_channel( channel_adder );
                one_jet_one_muon.add_channel( channel_adder );
                one_jet_one_electron.add_channel( channel_adder );
                OSSF_minus_OSDF_leptons.add_channel( channel_adder );
                no_jets_two_muons.add_channel( channel_adder );
                no_jets_two_electrons.add_channel( channel_adder );
                no_jets_one_muon_one_electron.add_channel( channel_adder );
                one_jet_no_leptons.add_channel( channel_adder );
                no_jets_one_muon.add_channel( channel_adder );
                no_jets_one_electron.add_channel( channel_adder );
                no_jets_no_leptons.add_channel( channel_adder );
                channels.push_back( channel_adder );

              }
            if( given_decaying_EWino->get_absolute_mass()
                > given_shortcuts->get_neutralino_one()->get_absolute_mass() )
              {

                // debugging:
                /**std::cout
                << std::endl
                << "adding 3-body channel";
                std::cout << std::endl;**/

                channel_adder
                = new charge_summed_neutralino_three_body( given_primary_cut,
                                                           given_secondary_cut,
                                                           given_jet_cut,
                                                           given_kinematics,
                                                       given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                          given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                                           given_shortcuts,
                                       two_jets_no_leptons.get_value_pointer(),
                                          one_jet_one_muon.get_value_pointer(),
                                      one_jet_one_electron.get_value_pointer(),
                                   OSSF_minus_OSDF_leptons.get_value_pointer(),
                                         no_jets_two_muons.get_value_pointer(),
                                     no_jets_two_electrons.get_value_pointer(),
                             no_jets_one_muon_one_electron.get_value_pointer(),
                                        one_jet_no_leptons.get_value_pointer(),
                                          no_jets_one_muon.get_value_pointer(),
                                      no_jets_one_electron.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
                two_jets_no_leptons.add_channel( channel_adder );
                one_jet_one_muon.add_channel( channel_adder );
                one_jet_one_electron.add_channel( channel_adder );
                OSSF_minus_OSDF_leptons.add_channel( channel_adder );
                no_jets_two_muons.add_channel( channel_adder );
                no_jets_two_electrons.add_channel( channel_adder );
                no_jets_one_muon_one_electron.add_channel( channel_adder );
                one_jet_no_leptons.add_channel( channel_adder );
                no_jets_one_muon.add_channel( channel_adder );
                no_jets_one_electron.add_channel( channel_adder );
                no_jets_no_leptons.add_channel( channel_adder );
                channels.push_back( channel_adder );

              }

          }
        else if( given_shortcuts->is_in( given_decaying_EWino->get_PDG_code(),
                                           given_shortcuts->get_charginos() ) )
          {

            if( given_decaying_EWino->get_absolute_mass()
                > given_shortcuts->get_selectron_L()->get_absolute_mass() )
              {

                // debugging:
                /**std::cout
                << std::endl
                << "adding selectron_L channel";
                std::cout << std::endl;**/

                channel_adder
                = new chargino_to_semu( given_primary_cut,
                                        given_secondary_cut,
                                        given_kinematics,
                                        given_decaying_scolored,
                                        given_scolored_is_not_antiparticle,
                                        given_decaying_EWino,
                                        given_EWino_is_not_antiparticle,
                                        given_shortcuts->get_selectron_L(),
                                        given_shortcuts,
                                      no_jets_one_electron.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
                no_jets_one_electron.add_channel( channel_adder );
                no_jets_no_leptons.add_channel( channel_adder );
                channels.push_back( channel_adder );

              }
            if( given_decaying_EWino->get_absolute_mass()
                > given_shortcuts->get_selectron_R()->get_absolute_mass() )
              {

                // debugging:
                /**std::cout
                << std::endl
                << "adding selectron_R channel";
                std::cout << std::endl;**/

                channel_adder
                = new chargino_to_semu( given_primary_cut,
                                        given_secondary_cut,
                                        given_kinematics,
                                        given_decaying_scolored,
                                        given_scolored_is_not_antiparticle,
                                        given_decaying_EWino,
                                        given_EWino_is_not_antiparticle,
                                        given_shortcuts->get_selectron_R(),
                                        given_shortcuts,
                                      no_jets_one_electron.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
                no_jets_one_electron.add_channel( channel_adder );
                no_jets_no_leptons.add_channel( channel_adder );
                channels.push_back( channel_adder );

              }
            if( given_decaying_EWino->get_absolute_mass()
                > given_shortcuts->get_smuon_L()->get_absolute_mass() )
              {

                // debugging:
                /**std::cout
                << std::endl
                << "adding smuon_L channel";
                std::cout << std::endl;**/

                channel_adder
                = new chargino_to_semu( given_primary_cut,
                                        given_secondary_cut,
                                        given_kinematics,
                                        given_decaying_scolored,
                                        given_scolored_is_not_antiparticle,
                                        given_decaying_EWino,
                                        given_EWino_is_not_antiparticle,
                                        given_shortcuts->get_smuon_L(),
                                        given_shortcuts,
                                        no_jets_one_muon.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
                no_jets_one_muon.add_channel( channel_adder );
                no_jets_no_leptons.add_channel( channel_adder );
                channels.push_back( channel_adder );

              }
            if( given_decaying_EWino->get_absolute_mass()
                > given_shortcuts->get_smuon_R()->get_absolute_mass() )
              {

                // debugging:
                /**std::cout
                << std::endl
                << "adding smuon_R channel";
                std::cout << std::endl;**/

                channel_adder
                = new chargino_to_semu( given_primary_cut,
                                        given_secondary_cut,
                                        given_kinematics,
                                        given_decaying_scolored,
                                        given_scolored_is_not_antiparticle,
                                        given_decaying_EWino,
                                        given_EWino_is_not_antiparticle,
                                        given_shortcuts->get_smuon_R(),
                                        given_shortcuts,
                                        no_jets_one_muon.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
                no_jets_one_muon.add_channel( channel_adder );
                no_jets_no_leptons.add_channel( channel_adder );
                channels.push_back( channel_adder );

              }
            if( given_decaying_EWino->get_absolute_mass()
                > given_shortcuts->get_stau_one()->get_absolute_mass() )
              {

                // debugging:
                /**std::cout
                << std::endl
                << "adding stau_1 channel";
                std::cout << std::endl;**/

                channel_adder
                = new chargino_to_stau( given_primary_cut,
                                        given_secondary_cut,
                                        given_jet_cut,
                                        given_kinematics,
                                        given_decaying_scolored,
                                        given_scolored_is_not_antiparticle,
                                        given_decaying_EWino,
                                        given_EWino_is_not_antiparticle,
                                        given_shortcuts->get_stau_one(),
                                        given_shortcuts,
                                        one_jet_no_leptons.get_value_pointer(),
                                        no_jets_one_muon.get_value_pointer(),
                                      no_jets_one_electron.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
                one_jet_no_leptons.add_channel( channel_adder );
                no_jets_one_muon.add_channel( channel_adder );
                no_jets_one_electron.add_channel( channel_adder );
                no_jets_no_leptons.add_channel( channel_adder );
                channels.push_back( channel_adder );

              }
            if( given_decaying_EWino->get_absolute_mass()
                > given_shortcuts->get_stau_two()->get_absolute_mass() )
              {

                // debugging:
                /**std::cout
                << std::endl
                << "adding stau_2 channel";
                std::cout << std::endl;**/

                channel_adder
                = new chargino_to_stau( given_primary_cut,
                                        given_secondary_cut,
                                        given_jet_cut,
                                        given_kinematics,
                                        given_decaying_scolored,
                                        given_scolored_is_not_antiparticle,
                                        given_decaying_EWino,
                                        given_EWino_is_not_antiparticle,
                                        given_shortcuts->get_stau_two(),
                                        given_shortcuts,
                                        one_jet_no_leptons.get_value_pointer(),
                                        no_jets_one_muon.get_value_pointer(),
                                      no_jets_one_electron.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
                one_jet_no_leptons.add_channel( channel_adder );
                no_jets_one_muon.add_channel( channel_adder );
                no_jets_one_electron.add_channel( channel_adder );
                no_jets_no_leptons.add_channel( channel_adder );
                channels.push_back( channel_adder );

              }
            if( given_decaying_EWino->get_absolute_mass()
                > given_shortcuts->get_electron_sneutrino_L(
                                                       )->get_absolute_mass() )
              {

                // debugging:
                /**std::cout
                << std::endl
                << "adding e_sneutrino_L channel";
                std::cout << std::endl;**/

                channel_adder
                = new chargino_to_emu_sneutrino( given_primary_cut,
                                                 given_secondary_cut,
                                                 given_kinematics,
                                                 given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                 given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                   given_shortcuts->get_electron_sneutrino_L(),
                                                 given_shortcuts,
                                      no_jets_one_electron.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
                no_jets_one_electron.add_channel( channel_adder );
                no_jets_no_leptons.add_channel( channel_adder );
                channels.push_back( channel_adder );

              }
            if( given_shortcuts->get_electron_sneutrino_R(
                                                    )->has_mass_been_recorded()
                &&
                ( given_decaying_EWino->get_absolute_mass()
                  > given_shortcuts->get_electron_sneutrino_R(
                                                     )->get_absolute_mass() ) )
              {

                // debugging:
                /**std::cout
                << std::endl
                << "adding e_sneutrino_R channel";
                std::cout << std::endl;**/

                channel_adder
                = new chargino_to_emu_sneutrino( given_primary_cut,
                                                 given_secondary_cut,
                                                 given_kinematics,
                                                 given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                 given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                   given_shortcuts->get_electron_sneutrino_R(),
                                                 given_shortcuts,
                                      no_jets_one_electron.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
                no_jets_one_electron.add_channel( channel_adder );
                no_jets_no_leptons.add_channel( channel_adder );
                channels.push_back( channel_adder );

              }
            if( given_decaying_EWino->get_absolute_mass()
                > given_shortcuts->get_muon_sneutrino_L(
                                                       )->get_absolute_mass() )
              {

                // debugging:
                /**std::cout
                << std::endl
                << "adding mu_sneutrino_L channel";
                std::cout << std::endl;**/

                channel_adder
                = new chargino_to_emu_sneutrino( given_primary_cut,
                                                 given_secondary_cut,
                                                 given_kinematics,
                                                 given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                 given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                       given_shortcuts->get_muon_sneutrino_L(),
                                                 given_shortcuts,
                                          no_jets_one_muon.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
                no_jets_one_muon.add_channel( channel_adder );
                no_jets_no_leptons.add_channel( channel_adder );
                channels.push_back( channel_adder );

              }
            if( given_shortcuts->get_muon_sneutrino_R(
                                                    )->has_mass_been_recorded()
                &&
                ( given_decaying_EWino->get_absolute_mass()
                  > given_shortcuts->get_muon_sneutrino_R(
                                                     )->get_absolute_mass() ) )
              {

                // debugging:
                /**std::cout
                << std::endl
                << "adding mu_sneutrino_R channel";
                std::cout << std::endl;**/

                channel_adder
                = new chargino_to_emu_sneutrino( given_primary_cut,
                                                 given_secondary_cut,
                                                 given_kinematics,
                                                 given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                 given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                      given_shortcuts->get_muon_sneutrino_R(),
                                                 given_shortcuts,
                                          no_jets_one_muon.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
                no_jets_one_muon.add_channel( channel_adder );
                no_jets_no_leptons.add_channel( channel_adder );
                channels.push_back( channel_adder );

              }
            if( given_decaying_EWino->get_absolute_mass()
                > given_shortcuts->get_tau_sneutrino_L()->get_absolute_mass() )
              {

                // debugging:
                /**std::cout
                << std::endl
                << "adding tau_sneutrino_L channel";
                std::cout << std::endl;**/

                channel_adder
                = new chargino_to_tau_sneutrino( given_primary_cut,
                                                 given_secondary_cut,
                                                 given_jet_cut,
                                                 given_kinematics,
                                                 given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                 given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                        given_shortcuts->get_tau_sneutrino_L(),
                                                 given_shortcuts,
                                        one_jet_no_leptons.get_value_pointer(),
                                          no_jets_one_muon.get_value_pointer(),
                                      no_jets_one_electron.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
                one_jet_no_leptons.add_channel( channel_adder );
                no_jets_one_muon.add_channel( channel_adder );
                no_jets_one_electron.add_channel( channel_adder );
                no_jets_no_leptons.add_channel( channel_adder );
                channels.push_back( channel_adder );

              }
            if( given_shortcuts->get_tau_sneutrino_R(
                                                    )->has_mass_been_recorded()
                &&
                ( given_decaying_EWino->get_absolute_mass()
                  > given_shortcuts->get_tau_sneutrino_R(
                                                     )->get_absolute_mass() ) )
              {

                // debugging:
                /**std::cout
                << std::endl
                << "adding tau_sneutrino_R channel";
                std::cout << std::endl;**/

                channel_adder
                = new chargino_to_tau_sneutrino( given_primary_cut,
                                                 given_secondary_cut,
                                                 given_jet_cut,
                                                 given_kinematics,
                                                 given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                 given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                        given_shortcuts->get_tau_sneutrino_R(),
                                                 given_shortcuts,
                                        one_jet_no_leptons.get_value_pointer(),
                                          no_jets_one_muon.get_value_pointer(),
                                      no_jets_one_electron.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
                one_jet_no_leptons.add_channel( channel_adder );
                no_jets_one_muon.add_channel( channel_adder );
                no_jets_one_electron.add_channel( channel_adder );
                no_jets_no_leptons.add_channel( channel_adder );
                channels.push_back( channel_adder );

              }
            if( given_decaying_EWino->get_absolute_mass()
                > ( given_shortcuts->get_W_plus()->get_absolute_mass()
               + given_shortcuts->get_neutralino_one()->get_absolute_mass() ) )
              {

                // debugging:
                /**std::cout
                << std::endl
                << "adding W channel";
                std::cout << std::endl;**/

                channel_adder
                = new chargino_to_W( given_primary_cut,
                                     given_secondary_cut,
                                     given_jet_cut,
                                     given_kinematics,
                                     given_decaying_scolored,
                                     given_scolored_is_not_antiparticle,
                                     given_decaying_EWino,
                                     given_EWino_is_not_antiparticle,
                                     given_shortcuts->get_W_plus(),
                                     given_shortcuts,
                                     two_jets_no_leptons.get_value_pointer(),
                                     one_jet_no_leptons.get_value_pointer(),
                                     no_jets_one_muon.get_value_pointer(),
                                     no_jets_one_electron.get_value_pointer(),
                                     no_jets_no_leptons.get_value_pointer() );
                two_jets_no_leptons.add_channel( channel_adder );
                one_jet_no_leptons.add_channel( channel_adder );
                no_jets_one_muon.add_channel( channel_adder );
                no_jets_one_electron.add_channel( channel_adder );
                no_jets_no_leptons.add_channel( channel_adder );
                channels.push_back( channel_adder );

              }
            if( given_decaying_EWino->get_absolute_mass()
                > ( given_shortcuts->get_charged_EWSB_scalar(
                                                         )->get_absolute_mass()
                    + given_shortcuts->get_neutralino_one(
                                                     )->get_absolute_mass() ) )
              {

                // debugging:
                /**std::cout
                << std::endl
                << "adding H channel";
                std::cout << std::endl;**/

                channel_adder
                = new chargino_to_Higgs( given_primary_cut,
                                         given_secondary_cut,
                                         given_jet_cut,
                                         given_kinematics,
                                         given_decaying_scolored,
                                         given_scolored_is_not_antiparticle,
                                         given_decaying_EWino,
                                         given_EWino_is_not_antiparticle,
                                    given_shortcuts->get_charged_EWSB_scalar(),
                                         given_shortcuts,
                                       two_jets_no_leptons.get_value_pointer(),
                                        one_jet_no_leptons.get_value_pointer(),
                                         no_jets_one_muon.get_value_pointer(),
                                      no_jets_one_electron.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
                two_jets_no_leptons.add_channel( channel_adder );
                one_jet_no_leptons.add_channel( channel_adder );
                no_jets_one_muon.add_channel( channel_adder );
                no_jets_one_electron.add_channel( channel_adder );
                no_jets_no_leptons.add_channel( channel_adder );
                channels.push_back( channel_adder );

              }
            if( given_decaying_EWino->get_absolute_mass()
                > given_shortcuts->get_neutralino_one()->get_absolute_mass() )
              {

                // debugging:
                /**std::cout
                << std::endl
                << "adding 3-body channel";
                std::cout << std::endl;**/

                channel_adder
                = new chargino_three_body( given_primary_cut,
                                           given_secondary_cut,
                                           given_jet_cut,
                                           given_kinematics,
                                           given_decaying_scolored,
                                           given_scolored_is_not_antiparticle,
                                           given_decaying_EWino,
                                           given_EWino_is_not_antiparticle,
                                           given_shortcuts,
                                       two_jets_no_leptons.get_value_pointer(),
                                        one_jet_no_leptons.get_value_pointer(),
                                          no_jets_one_muon.get_value_pointer(),
                                      no_jets_one_electron.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
                two_jets_no_leptons.add_channel( channel_adder );
                one_jet_no_leptons.add_channel( channel_adder );
                no_jets_one_muon.add_channel( channel_adder );
                no_jets_one_electron.add_channel( channel_adder );
                no_jets_no_leptons.add_channel( channel_adder );
                channels.push_back( channel_adder );

              }

          }

      }

  }

  charge_summed_cascade_calculator::~charge_summed_cascade_calculator()
  {

    for( std::vector< channel_calculator* >::iterator
         deletion_iterator = channels.begin();
         channels.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

  }


  double
  charge_summed_cascade_calculator::get_acceptance( int const number_of_jets,
                                                  int const number_of_leptons )
  {

    // debugging:
    /**std::cout
    << std::endl
    << "charge_summed_cascade_calculator::get_acceptance( "
    << number_of_jets << ", " << number_of_leptons << ") called.";
    std::cout << std::endl;**/

    if( 2 == number_of_jets )
      {

        if( 0 == number_of_leptons )
          {

            // debugging:
            /**std::cout
            << "returning " << two_jets_no_leptons.getValue();
            std::cout << std::endl;**/

            return two_jets_no_leptons.get_value();

          }
        else
          {

            // debugging:
            /**std::cout
            << "returning 0.0";
            std::cout << std::endl;**/

            return 0.0;

          }

      }
    else if( 1 == number_of_jets )
      {

        if( 1 == number_of_leptons )
          {

            // debugging:
            /**std::cout
            << "returning " << ( one_jet_one_electron.getValue()
                                 + one_jet_one_muon.getValue() );
            std::cout << std::endl;**/

            return ( one_jet_one_electron.get_value()
                     + one_jet_one_muon.get_value() );

          }
        else if( 0 == number_of_leptons )
          {

            // debugging:
            /**std::cout
            << "returning " << one_jet_no_leptons.getValue();
            std::cout << std::endl;**/

            return one_jet_no_leptons.get_value();

          }
        else
          {

            // debugging:
            /**std::cout
            << "returning 0.0";
            std::cout << std::endl;**/

            return 0.0;

          }

      }
    else if( 0 == number_of_jets )
      {

        if( 2 == number_of_leptons )
          {

            // debugging:
            /**std::cout
            << "returning " << ( no_jets_two_electrons.getValue()
                                 + no_jets_two_muons.getValue()
                                 + no_jets_one_muon_one_electron.getValue() );
            std::cout << std::endl;**/

            return ( no_jets_two_electrons.get_value()
                     + no_jets_two_muons.get_value()
                     + no_jets_one_muon_one_electron.get_value() );

          }
        else if( 1 == number_of_leptons )
          {

            // debugging:
            /**std::cout
            << "returning " << ( no_jets_one_electron.getValue()
                                 + no_jets_one_muon.getValue() );
            std::cout << std::endl;**/

            return ( no_jets_one_electron.get_value()
                     + no_jets_one_muon.get_value() );

          }
        else if( 0 == number_of_leptons )
          {

            // debugging:
            /**std::cout
            << "returning " << no_jets_no_leptons.getValue();
            std::cout << std::endl;**/

            return no_jets_no_leptons.get_value();

          }
        else
          {

            // debugging:
            /**std::cout
            << "returning 0.0";
            std::cout << std::endl;**/

            return 0.0;

          }

      }
    else
      {

        // debugging:
        /**std::cout
        << "returning 0.0";
        std::cout << std::endl;**/

        return 0.0;

      }

  }

  double
  charge_summed_cascade_calculator::get_muon_acceptance(
                                                      int const number_of_jets,
                                                    int const number_of_muons )
  {

    if( 2 == number_of_jets )
      {

        if( 0 == number_of_muons )
          {

            return two_jets_no_leptons.get_value();

          }
        else
          {

            return 0.0;

          }

      }
    else if( 1 == number_of_jets )
      {

        if( 1 == number_of_muons )
          {

            return one_jet_one_muon.get_value();

          }
        else if( 0 == number_of_muons )
          {

            return one_jet_no_leptons.get_value();

          }
        else
          {

            return 0.0;

          }

      }
    else if( 0 == number_of_jets )
      {

        if( 2 == number_of_muons )
          {

            return no_jets_two_muons.get_value();

          }
        else if( 1 == number_of_muons )
          {

            return no_jets_one_muon.get_value();

          }
        else if( 0 == number_of_muons )
          {

            return no_jets_no_leptons.get_value();

          }
        else
          {

            return 0.0;

          }

      }
    else
      {

        return 0.0;

      }

  }

  double
  charge_summed_cascade_calculator::get_electron_acceptance(
                                                      int const number_of_jets,
                                                int const number_of_electrons )
  {

    if( 2 == number_of_jets )
      {

        if( 0 == number_of_electrons )
          {

            return two_jets_no_leptons.get_value();

          }
        else
          {

            return 0.0;

          }

      }
    else if( 1 == number_of_jets )
      {

        if( 1 == number_of_electrons )
          {

            return one_jet_one_electron.get_value();

          }
        else if( 0 == number_of_electrons )
          {

            return one_jet_no_leptons.get_value();

          }
        else
          {

            return 0.0;

          }

      }
    else if( 0 == number_of_jets )
      {

        if( 2 == number_of_electrons )
          {

            return no_jets_two_electrons.get_value();

          }
        else if( 1 == number_of_electrons )
          {

            return no_jets_one_electron.get_value();

          }
        else if( 0 == number_of_electrons )
          {

            return no_jets_no_leptons.get_value();

          }
        else
          {

            return 0.0;

          }

      }
    else
      {

        return 0.0;

      }

  }



  cascade_handler::cascade_handler(
                                 input_handler const* const given_shortcuts ) :
    readied_for_new_point( given_shortcuts->get_readier() ),
    shortcut( given_shortcuts )
  {

    // just an initialiation list.

  }

  cascade_handler::~cascade_handler()
  {

    for( std::vector< charge_summed_cascade_calculator* >::iterator
         deletion_iterator = calculators.begin();
         calculators.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

  }


  charge_summed_cascade_calculator*
  cascade_handler::get_charge_summed_calculator( double given_primary_cut,
                                                 double given_secondary_cut,
                                                 double given_jet_cut,
                               lepton_acceptance_value* const given_kinematics,
                 CppSLHA::particle_property_set const* const decaying_scolored,
                                       bool const scolored_is_not_antiparticle,
                    CppSLHA::particle_property_set const* const decaying_EWino,
                                         bool const EWino_is_not_antiparticle )
  {

    // debugging:
    /**std::cout
    << std::endl
    << "cascade_handler::get_charge_summed_calculator( "
    << given_primary_cut << ", "
    << given_secondary_cut << ", "
    << given_jet_cut << ", "
    << given_kinematics << ", "
    << *(decaying_scolored->get_name_or_antiname(
                                                scolored_is_not_antiparticle ))
    << ", "
    << *(decaying_EWino->get_name_or_antiname( EWino_is_not_antiparticle ))
    << " ) called.";
    std::cout << std::endl;**/

    charge_summed_cascade_calculator* return_pointer = NULL;

    if( needs_to_prepare_for_this_point() )
      // if we have a new point, we should drop all the held pointers, since
      // they are no longer relevant & cannot be re-used.
      {

        for( std::vector< charge_summed_cascade_calculator* >::iterator
             deletion_iterator = calculators.begin();
             calculators.end() > deletion_iterator;
             ++deletion_iterator )
          {

            delete *deletion_iterator;

          }

        calculators.clear();
        finish_preparing_for_this_point();

      }
    else
      // otherwise, we look through the existing pointers to see if they can
      // be re-used....
      {

        for( std::vector< charge_summed_cascade_calculator* >::iterator
             search_iterator = calculators.begin();
             calculators.end() > search_iterator;
             ++search_iterator )
          // look to see if this calculator already exists...
          {

            if( (*search_iterator)->is_requested( given_primary_cut,
                                                  given_secondary_cut,
                                                  given_jet_cut,
                                                  given_kinematics,
                                                  decaying_scolored,
                                                  scolored_is_not_antiparticle,
                                                  decaying_EWino,
                                                  EWino_is_not_antiparticle ) )
              {

                return_pointer = *search_iterator;
                // note the pre-existing calculator.

                search_iterator = calculators.end();
                // stop searching.

              }

          }

      }

    if( NULL == return_pointer )
      // if we did not find the calculator among those already constructed...
      {

        // all "no cut required" cuts are set to a standard value for
        // comparison:
        if( given_primary_cut <= 0.0 )
          {

            given_primary_cut = -1.0;

          }
        if( given_secondary_cut <= 0.0 )
          {

            given_secondary_cut = -1.0;

          }
        if( given_jet_cut <= 0.0 )
          {

            given_jet_cut = -1.0;

          }

        // make a new one:
        return_pointer
        = new charge_summed_cascade_calculator( given_primary_cut,
                                                given_secondary_cut,
                                                given_jet_cut,
                                                given_kinematics,
                                                decaying_scolored,
                                                scolored_is_not_antiparticle,
                                                decaying_EWino,
                                                EWino_is_not_antiparticle,
                                                shortcut );

        // keep hold of it:
        calculators.push_back( return_pointer );

      }


    // return the pointer to the appropriate object:
    return return_pointer;

  }

}  // end of LHC_FASER namespace.
