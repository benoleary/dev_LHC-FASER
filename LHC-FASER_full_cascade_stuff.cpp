/*
 * LHC-FASER_full_cascade_stuff.cpp
 *
 *  Created on: 30 Mar 2011
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

#include "LHC-FASER_full_cascade_stuff.hpp"

namespace LHC_FASER
{

  // this is here just for reference! it will be deleted later, when I am sure
  // that I have covered all the decays it covers.
  EW_decayer_values::EW_decayer_values( bool const given_should_sum_charges,
                                        double const given_primary_cut,
                                        double const given_secondary_cut,
                                        double const given_jet_cut,
                               lepton_acceptance_value* const given_kinematics,
           CppSLHA::particle_property_set const* const given_decaying_scolored,
                                 bool const given_scolored_is_not_antiparticle,
                  CppSLHA::particle_property_set const* const given_EW_decayer,
                                    bool const given_EWino_is_not_antiparticle,
                                 input_handler const* const given_shortcuts ) :
    should_sum_charges( given_should_sum_charges ),
    primary_cut( given_primary_cut ),
    secondary_cut( given_secondary_cut ),
    jetCut( given_jet_cut ),
    leptonKinematics( given_kinematics ),
    decaying_scolored( given_decaying_scolored ),
    EW_decayer( given_EW_decayer ),
    OSSF_minus_OSDF_leptons( CppSLHA::CppSLHA_global::really_wrong_value ),
    two_jets_no_leptons( CppSLHA::CppSLHA_global::really_wrong_value ),
    oneJetOneNegativeElectron(
                                 CppSLHA::CppSLHA_global::really_wrong_value ),
    oneJetOnePositiveElectron(
                                 CppSLHA::CppSLHA_global::really_wrong_value ),
    oneJetOneNegativeMuon( CppSLHA::CppSLHA_global::really_wrong_value ),
    oneJetOnePositiveMuon( CppSLHA::CppSLHA_global::really_wrong_value ),
    no_jets_two_electrons( CppSLHA::CppSLHA_global::really_wrong_value ),
    no_jets_one_positive_electron_one_negative_muon(
                                 CppSLHA::CppSLHA_global::really_wrong_value ),
    no_jets_one_negative_electron_one_positive_muon(
                                 CppSLHA::CppSLHA_global::really_wrong_value ),
    no_jets_two_muons( CppSLHA::CppSLHA_global::really_wrong_value ),
    oneJetZeroLeptons( CppSLHA::CppSLHA_global::really_wrong_value ),
    zeroJetsOneNegativeElectron(
                                 CppSLHA::CppSLHA_global::really_wrong_value ),
    zerJetsOnePositiveElectron(
                                 CppSLHA::CppSLHA_global::really_wrong_value ),
    zeroJetsOneNegativeMuon( CppSLHA::CppSLHA_global::really_wrong_value ),
    zeroJetsOnePositiveMuon( CppSLHA::CppSLHA_global::really_wrong_value ),
    zeroJetsZeroLeptons( CppSLHA::CppSLHA_global::really_wrong_value )
  {

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: EW_decayer_values::EW_decayer_values( "
    << given_primary_cut << ", "
    << given_secondary_cut << ", "
    << given_jet_cut << ", ..."
    << *(given_decaying_scolored->getName())
    << ", "
    << *(given_decaying_EWino->getName())
    << ", ... ) called.";
    std::cout << std::endl;
    std::cout << std::endl;**/


    if( given_should_sum_charges )
      {

        charge_summing_factor = 2.0;

      }
    else
      {

        charge_summing_factor = 1.0;

      }

    channel_calculator* channel_adder;

    if( given_decaying_EWino == given_shortcuts->get_neutralino_one() )
      // if the electroweakino is the lightest neutralino, & hence assumed
      // stable, we set up the appropriate calculator...
      {

        // debugging:
        /**std::cout
        << std::endl
        << "adding neutralino_1 channel";
        std::cout << std::endl;**/

        channel_adder
        = new lightestNeutralinoCascade( no_jets_no_leptons.get_value_pointer() );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );

      }
    else if( given_shortcuts->is_in( given_decaying_EWino->get_PDG_code(),
                                given_shortcuts->get_unstable_neutralinos() ) )
      {

        // debugging:
        /**std::cout
        << std::endl
        << "adding selectron_L channel";
        std::cout << std::endl;**/

        channel_adder = new neutralinoToSemuCascade( should_sum_charges,
                                                given_primary_cut,
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
                             no_jets_one_negative_electron.get_value_pointer(),
                             no_jets_one_positive_electron.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        OSSF_minus_OSDF_leptons.add_channel( channel_adder );
        no_jets_two_electrons.add_channel( channel_adder );
        no_jets_one_negative_electron.add_channel( channel_adder );
        no_jets_one_positive_electron.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );


        // debugging:
        /**std::cout
        << std::endl
        << "adding selectron_R channel";
        std::cout << std::endl;**/

        channel_adder = new neutralinoToSemuCascade( should_sum_charges,
                                                given_primary_cut,
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
                             no_jets_one_negative_electron.get_value_pointer(),
                             no_jets_one_positive_electron.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        OSSF_minus_OSDF_leptons.add_channel( channel_adder );
        no_jets_two_electrons.add_channel( channel_adder );
        no_jets_one_negative_electron.add_channel( channel_adder );
        no_jets_one_positive_electron.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );


        // debugging:
        /**std::cout
        << std::endl
        << "adding smuon_L channel";
        std::cout << std::endl;**/

        channel_adder = new neutralinoToSemuCascade( should_sum_charges,
                                                given_primary_cut,
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
                                 no_jets_one_negative_muon.get_value_pointer(),
                                 no_jets_one_positive_muon.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        OSSF_minus_OSDF_leptons.add_channel( channel_adder );
        no_jets_two_muons.add_channel( channel_adder );
        no_jets_one_negative_muon.add_channel( channel_adder );
        no_jets_one_positive_muon.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );


        // debugging:
        /**std::cout
        << std::endl
        << "adding smuon_R channel";
        std::cout << std::endl;**/

        channel_adder = new neutralinoToSemuCascade( should_sum_charges,
                                                given_primary_cut,
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
                                 no_jets_one_negative_muon.get_value_pointer(),
                                 no_jets_one_positive_muon.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        OSSF_minus_OSDF_leptons.add_channel( channel_adder );
        no_jets_two_muons.add_channel( channel_adder );
        no_jets_one_negative_muon.add_channel( channel_adder );
        no_jets_one_positive_muon.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );


        // debugging:
        /**std::cout
        << std::endl
        << "adding stau_one channels";
        std::cout << std::endl;**/

        channel_adder = new neutralino_to_stau_jj( should_sum_charges,
                                                   given_jet_cut,
                                                   given_kinematics,
                                                   given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                   given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                               given_shortcuts->get_stau_one(),
                                                   given_shortcuts,
                                       two_jets_no_leptons.get_value_pointer(),
                                        one_jet_no_leptons.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        two_jets_no_leptons.add_channel( channel_adder );
        one_jet_no_leptons.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );

        channel_adder = new neutralino_to_stau_jl( should_sum_charges,
                                                   given_primary_cut,
                                                   given_secondary_cut,
                                                   given_jet_cut,
                                                   given_kinematics,
                                                   given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                   given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                               given_shortcuts->get_stau_one(),
                                                   given_shortcuts,
                             one_jet_one_negative_electron.get_value_pointer(),
                             one_jet_one_positive_electron.get_value_pointer(),
                                 one_jet_one_negative_muon.get_value_pointer(),
                                 one_jet_one_positive_muon.get_value_pointer(),
                                        one_jet_no_leptons.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        one_jet_one_negative_electron.add_channel( channel_adder );
        one_jet_one_positive_electron.add_channel( channel_adder );
        one_jet_one_negative_muon.add_channel( channel_adder );
        one_jet_one_positive_muon.add_channel( channel_adder );
        one_jet_no_leptons.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );

        channel_adder = new neutralino_to_stau_ll( should_sum_charges,
                                                   given_primary_cut,
                                                   given_secondary_cut,
                                                   given_kinematics,
                                                   given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                   given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                               given_shortcuts->get_stau_one(),
                                                   given_shortcuts,
                                     no_jets_two_electrons.get_value_pointer(),
           no_jets_one_negative_electron_one_positive_muon.get_value_pointer(),
           no_jets_one_positive_electron_one_negative_muon.get_value_pointer(),
                                         no_jets_two_muons.get_value_pointer(),
                             no_jets_one_negative_electron.get_value_pointer(),
                             no_jets_one_positive_electron.get_value_pointer(),
                                 no_jets_one_negative_muon.get_value_pointer(),
                                 no_jets_one_positive_muon.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        no_jets_two_electrons.add_channel( channel_adder );
        no_jets_one_negative_electron_one_positive_muon.add_channel(
                                                               channel_adder );
        no_jets_one_positive_electron_one_negative_muon.add_channel(
                                                               channel_adder );
        no_jets_two_muons.add_channel( channel_adder );
        no_jets_one_negative_electron.add_channel( channel_adder );
        no_jets_one_positive_electron.add_channel( channel_adder );
        no_jets_one_negative_muon.add_channel( channel_adder );
        no_jets_one_positive_muon.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );


        // debugging:
        /**std::cout
        << std::endl
        << "adding stau_two channels";
        std::cout << std::endl;**/

        channel_adder = new neutralino_to_stau_jj( should_sum_charges,
                                                   given_jet_cut,
                                                   given_kinematics,
                                                   given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                   given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                               given_shortcuts->get_stau_two(),
                                                   given_shortcuts,
                                       two_jets_no_leptons.get_value_pointer(),
                                        one_jet_no_leptons.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        two_jets_no_leptons.add_channel( channel_adder );
        one_jet_no_leptons.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );

        channel_adder = new neutralino_to_stau_jl( should_sum_charges,
                                                   given_primary_cut,
                                                   given_secondary_cut,
                                                   given_jet_cut,
                                                   given_kinematics,
                                                   given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                   given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                               given_shortcuts->get_stau_two(),
                                                   given_shortcuts,
                             one_jet_one_negative_electron.get_value_pointer(),
                             one_jet_one_positive_electron.get_value_pointer(),
                                 one_jet_one_negative_muon.get_value_pointer(),
                                 one_jet_one_positive_muon.get_value_pointer(),
                                        one_jet_no_leptons.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        one_jet_one_negative_electron.add_channel( channel_adder );
        one_jet_one_positive_electron.add_channel( channel_adder );
        one_jet_one_negative_muon.add_channel( channel_adder );
        one_jet_one_positive_muon.add_channel( channel_adder );
        one_jet_no_leptons.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );

        channel_adder = new neutralino_to_stau_ll( should_sum_charges,
                                                   given_primary_cut,
                                                   given_secondary_cut,
                                                   given_kinematics,
                                                   given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                   given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                               given_shortcuts->get_stau_two(),
                                                   given_shortcuts,
                                     no_jets_two_electrons.get_value_pointer(),
           no_jets_one_negative_electron_one_positive_muon.get_value_pointer(),
           no_jets_one_positive_electron_one_negative_muon.get_value_pointer(),
                                         no_jets_two_muons.get_value_pointer(),
                             no_jets_one_negative_electron.get_value_pointer(),
                             no_jets_one_positive_electron.get_value_pointer(),
                                 no_jets_one_negative_muon.get_value_pointer(),
                                 no_jets_one_positive_muon.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        no_jets_two_electrons.add_channel( channel_adder );
        no_jets_one_negative_electron_one_positive_muon.add_channel(
                                                               channel_adder );
        no_jets_one_positive_electron_one_negative_muon.add_channel(
                                                               channel_adder );
        no_jets_two_muons.add_channel( channel_adder );
        no_jets_one_negative_electron.add_channel( channel_adder );
        no_jets_one_positive_electron.add_channel( channel_adder );
        no_jets_one_negative_muon.add_channel( channel_adder );
        no_jets_one_positive_muon.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );


        // debugging:
        /**std::cout
        << std::endl
        << "adding Z channels";
        std::cout << std::endl;**/

        channel_adder = new neutralino_to_Z_without_taus( given_primary_cut,
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
                                   OSSF_minus_OSDF_leptons.get_value_pointer(),
                                     no_jets_two_electrons.get_value_pointer(),
                                         no_jets_two_muons.get_value_pointer(),
                                        one_jet_no_leptons.get_value_pointer(),
                             no_jets_one_negative_electron.get_value_pointer(),
                             no_jets_one_positive_electron.get_value_pointer(),
                                 no_jets_one_negative_muon.get_value_pointer(),
                                 no_jets_one_positive_muon.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        two_jets_no_leptons.add_channel( channel_adder );
        OSSF_minus_OSDF_leptons.add_channel( channel_adder );
        no_jets_two_electrons.add_channel( channel_adder );
        no_jets_two_muons.add_channel( channel_adder );
        one_jet_no_leptons.add_channel( channel_adder );
        no_jets_one_negative_electron.add_channel( channel_adder );
        no_jets_one_positive_electron.add_channel( channel_adder );
        no_jets_one_negative_muon.add_channel( channel_adder );
        no_jets_one_positive_muon.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );

        channel_adder = new neutralino_to_Z_to_taus_jj( given_primary_cut,
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
                                        one_jet_no_leptons.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        two_jets_no_leptons.add_channel( channel_adder );
        one_jet_no_leptons.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );

        channel_adder = new neutralino_to_Z_to_taus_jl( given_primary_cut,
                                                        given_secondary_cut,
                                                        given_jet_cut,
                                                        given_kinematics,
                                                       given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                        given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                                      given_shortcuts->get_Z(),
                                                        given_shortcuts,
                             one_jet_one_negative_electron.get_value_pointer(),
                             one_jet_one_positive_electron.get_value_pointer(),
                                 one_jet_one_negative_muon.get_value_pointer(),
                                 one_jet_one_positive_muon.get_value_pointer(),
                                        one_jet_no_leptons.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        one_jet_one_negative_electron.add_channel( channel_adder );
        one_jet_one_positive_electron.add_channel( channel_adder );
        one_jet_one_negative_muon.add_channel( channel_adder );
        one_jet_one_positive_muon.add_channel( channel_adder );
        one_jet_no_leptons.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );

        channel_adder = new neutralino_to_Z_to_taus_ll( given_primary_cut,
                                                        given_secondary_cut,
                                                        given_jet_cut,
                                                        given_kinematics,
                                                       given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                        given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                                      given_shortcuts->get_Z(),
                                                        given_shortcuts,
                                     no_jets_two_electrons.get_value_pointer(),
           no_jets_one_negative_electron_one_positive_muon.get_value_pointer(),
           no_jets_one_positive_electron_one_negative_muon.get_value_pointer(),
                                         no_jets_two_muons.get_value_pointer(),
                             no_jets_one_negative_electron.get_value_pointer(),
                             no_jets_one_positive_electron.get_value_pointer(),
                                 no_jets_one_negative_muon.get_value_pointer(),
                                 no_jets_one_positive_muon.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        no_jets_two_electrons.add_channel( channel_adder );
        no_jets_one_negative_electron_one_positive_muon.add_channel(
                                                               channel_adder );
        no_jets_one_positive_electron_one_negative_muon.add_channel(
                                                               channel_adder );
        no_jets_two_muons.add_channel( channel_adder );
        no_jets_one_negative_electron.add_channel( channel_adder );
        no_jets_one_positive_electron.add_channel( channel_adder );
        no_jets_one_negative_muon.add_channel( channel_adder );
        no_jets_one_positive_muon.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );


        // debugging:
        /**std::cout
        << std::endl
        << "adding h channels";
        std::cout << std::endl;**/

        channel_adder
        = new neutralino_to_Higgs_without_taus( given_primary_cut,
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
                                   OSSF_minus_OSDF_leptons.get_value_pointer(),
                                         no_jets_two_muons.get_value_pointer(),
                                        one_jet_no_leptons.get_value_pointer(),
                                 no_jets_one_negative_muon.get_value_pointer(),
                                 no_jets_one_positive_muon.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        two_jets_no_leptons.add_channel( channel_adder );
        OSSF_minus_OSDF_leptons.add_channel( channel_adder );
        no_jets_two_muons.add_channel( channel_adder );
        one_jet_no_leptons.add_channel( channel_adder );
        no_jets_one_negative_muon.add_channel( channel_adder );
        no_jets_one_positive_muon.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );

        channel_adder
        = new neutralino_to_Higgs_to_taus_jj( given_primary_cut,
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
                                        one_jet_no_leptons.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        two_jets_no_leptons.add_channel( channel_adder );
        one_jet_no_leptons.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );

        channel_adder
        = new neutralino_to_Higgs_to_taus_jl( given_primary_cut,
                                              given_secondary_cut,
                                              given_jet_cut,
                                              given_kinematics,
                                              given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                              given_decaying_EWino,
                                              given_EWino_is_not_antiparticle,
                              given_shortcuts->get_light_neutral_EWSB_scalar(),
                                              given_shortcuts,
                             one_jet_one_negative_electron.get_value_pointer(),
                             one_jet_one_positive_electron.get_value_pointer(),
                                 one_jet_one_negative_muon.get_value_pointer(),
                                 one_jet_one_positive_muon.get_value_pointer(),
                                        one_jet_no_leptons.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        one_jet_one_negative_electron.add_channel( channel_adder );
        one_jet_one_positive_electron.add_channel( channel_adder );
        one_jet_one_negative_muon.add_channel( channel_adder );
        one_jet_one_positive_muon.add_channel( channel_adder );
        one_jet_no_leptons.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );

        channel_adder
        = new neutralino_to_Higgs_to_taus_ll( given_primary_cut,
                                              given_secondary_cut,
                                              given_jet_cut,
                                              given_kinematics,
                                              given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                              given_decaying_EWino,
                                              given_EWino_is_not_antiparticle,
                              given_shortcuts->get_light_neutral_EWSB_scalar(),
                                              given_shortcuts,
                                     no_jets_two_electrons.get_value_pointer(),
           no_jets_one_negative_electron_one_positive_muon.get_value_pointer(),
           no_jets_one_positive_electron_one_negative_muon.get_value_pointer(),
                                         no_jets_two_muons.get_value_pointer(),
                             no_jets_one_negative_electron.get_value_pointer(),
                             no_jets_one_positive_electron.get_value_pointer(),
                                 no_jets_one_negative_muon.get_value_pointer(),
                                 no_jets_one_positive_muon.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        no_jets_two_electrons.add_channel( channel_adder );
        no_jets_one_negative_electron_one_positive_muon.add_channel(
                                                               channel_adder );
        no_jets_one_positive_electron_one_negative_muon.add_channel(
                                                               channel_adder );
        no_jets_two_muons.add_channel( channel_adder );
        no_jets_one_negative_electron.add_channel( channel_adder );
        no_jets_one_positive_electron.add_channel( channel_adder );
        no_jets_one_negative_muon.add_channel( channel_adder );
        no_jets_one_positive_muon.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );


        // debugging:
        /**std::cout
        << std::endl
        << "adding H channels";
        std::cout << std::endl;**/

        channel_adder
        = new neutralino_to_Higgs_without_taus( given_primary_cut,
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
                                   OSSF_minus_OSDF_leptons.get_value_pointer(),
                                         no_jets_two_muons.get_value_pointer(),
                                        one_jet_no_leptons.get_value_pointer(),
                                 no_jets_one_negative_muon.get_value_pointer(),
                                 no_jets_one_positive_muon.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        two_jets_no_leptons.add_channel( channel_adder );
        OSSF_minus_OSDF_leptons.add_channel( channel_adder );
        no_jets_two_muons.add_channel( channel_adder );
        one_jet_no_leptons.add_channel( channel_adder );
        no_jets_one_negative_muon.add_channel( channel_adder );
        no_jets_one_positive_muon.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );

        channel_adder
        = new neutralino_to_Higgs_to_taus_jj( given_primary_cut,
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
                                        one_jet_no_leptons.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        two_jets_no_leptons.add_channel( channel_adder );
        one_jet_no_leptons.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );

        channel_adder
        = new neutralino_to_Higgs_to_taus_jl( given_primary_cut,
                                              given_secondary_cut,
                                              given_jet_cut,
                                              given_kinematics,
                                              given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                              given_decaying_EWino,
                                              given_EWino_is_not_antiparticle,
                              given_shortcuts->get_heavy_neutral_EWSB_scalar(),
                                              given_shortcuts,
                             one_jet_one_negative_electron.get_value_pointer(),
                             one_jet_one_positive_electron.get_value_pointer(),
                                 one_jet_one_negative_muon.get_value_pointer(),
                                 one_jet_one_positive_muon.get_value_pointer(),
                                        one_jet_no_leptons.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        one_jet_one_negative_electron.add_channel( channel_adder );
        one_jet_one_positive_electron.add_channel( channel_adder );
        one_jet_one_negative_muon.add_channel( channel_adder );
        one_jet_one_positive_muon.add_channel( channel_adder );
        one_jet_no_leptons.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );

        channel_adder
        = new neutralino_to_Higgs_to_taus_ll( given_primary_cut,
                                              given_secondary_cut,
                                              given_jet_cut,
                                              given_kinematics,
                                              given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                              given_decaying_EWino,
                                              given_EWino_is_not_antiparticle,
                              given_shortcuts->get_heavy_neutral_EWSB_scalar(),
                                              given_shortcuts,
                                     no_jets_two_electrons.get_value_pointer(),
           no_jets_one_negative_electron_one_positive_muon.get_value_pointer(),
           no_jets_one_positive_electron_one_negative_muon.get_value_pointer(),
                                         no_jets_two_muons.get_value_pointer(),
                             no_jets_one_negative_electron.get_value_pointer(),
                             no_jets_one_positive_electron.get_value_pointer(),
                                 no_jets_one_negative_muon.get_value_pointer(),
                                 no_jets_one_positive_muon.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        no_jets_two_electrons.add_channel( channel_adder );
        no_jets_one_negative_electron_one_positive_muon.add_channel(
                                                               channel_adder );
        no_jets_one_positive_electron_one_negative_muon.add_channel(
                                                               channel_adder );
        no_jets_two_muons.add_channel( channel_adder );
        no_jets_one_negative_electron.add_channel( channel_adder );
        no_jets_one_positive_electron.add_channel( channel_adder );
        no_jets_one_negative_muon.add_channel( channel_adder );
        no_jets_one_positive_muon.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );


        // debugging:
        /**std::cout
        << std::endl
        << "adding A channels";
        std::cout << std::endl;**/

        channel_adder
        = new neutralino_to_Higgs_without_taus( given_primary_cut,
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
                                   OSSF_minus_OSDF_leptons.get_value_pointer(),
                                         no_jets_two_muons.get_value_pointer(),
                                        one_jet_no_leptons.get_value_pointer(),
                                 no_jets_one_negative_muon.get_value_pointer(),
                                 no_jets_one_positive_muon.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        two_jets_no_leptons.add_channel( channel_adder );
        OSSF_minus_OSDF_leptons.add_channel( channel_adder );
        no_jets_two_muons.add_channel( channel_adder );
        one_jet_no_leptons.add_channel( channel_adder );
        no_jets_one_negative_muon.add_channel( channel_adder );
        no_jets_one_positive_muon.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );

        channel_adder
        = new neutralino_to_Higgs_to_taus_jj( given_primary_cut,
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
                                        one_jet_no_leptons.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        two_jets_no_leptons.add_channel( channel_adder );
        one_jet_no_leptons.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );

        channel_adder
        = new neutralino_to_Higgs_to_taus_jl( given_primary_cut,
                                              given_secondary_cut,
                                              given_jet_cut,
                                              given_kinematics,
                                              given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                              given_decaying_EWino,
                                              given_EWino_is_not_antiparticle,
                              given_shortcuts->get_neutral_EWSB_pseudoscalar(),
                                              given_shortcuts,
                             one_jet_one_negative_electron.get_value_pointer(),
                             one_jet_one_positive_electron.get_value_pointer(),
                                 one_jet_one_negative_muon.get_value_pointer(),
                                 one_jet_one_positive_muon.get_value_pointer(),
                                        one_jet_no_leptons.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        one_jet_one_negative_electron.add_channel( channel_adder );
        one_jet_one_positive_electron.add_channel( channel_adder );
        one_jet_one_negative_muon.add_channel( channel_adder );
        one_jet_one_positive_muon.add_channel( channel_adder );
        one_jet_no_leptons.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );

        channel_adder
        = new neutralino_to_Higgs_to_taus_ll( given_primary_cut,
                                              given_secondary_cut,
                                              given_jet_cut,
                                              given_kinematics,
                                              given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                              given_decaying_EWino,
                                              given_EWino_is_not_antiparticle,
                              given_shortcuts->get_neutral_EWSB_pseudoscalar(),
                                              given_shortcuts,
                                     no_jets_two_electrons.get_value_pointer(),
           no_jets_one_negative_electron_one_positive_muon.get_value_pointer(),
           no_jets_one_positive_electron_one_negative_muon.get_value_pointer(),
                                         no_jets_two_muons.get_value_pointer(),
                             no_jets_one_negative_electron.get_value_pointer(),
                             no_jets_one_positive_electron.get_value_pointer(),
                                 no_jets_one_negative_muon.get_value_pointer(),
                                 no_jets_one_positive_muon.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        no_jets_two_electrons.add_channel( channel_adder );
        no_jets_one_negative_electron_one_positive_muon.add_channel(
                                                               channel_adder );
        no_jets_one_positive_electron_one_negative_muon.add_channel(
                                                               channel_adder );
        no_jets_two_muons.add_channel( channel_adder );
        no_jets_one_negative_electron.add_channel( channel_adder );
        no_jets_one_positive_electron.add_channel( channel_adder );
        no_jets_one_negative_muon.add_channel( channel_adder );
        no_jets_one_positive_muon.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );


        // debugging:
        /**std::cout
        << std::endl
        << "adding 3-body channel";
        std::cout << std::endl;**/

        channel_adder
        = new neutralino_three_body_without_taus( given_primary_cut,
                                                  given_secondary_cut,
                                                  given_jet_cut,
                                                  given_kinematics,
                                                  given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                  given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                                  given_shortcuts,
                                       two_jets_no_leptons.get_value_pointer(),
                                   OSSF_minus_OSDF_leptons.get_value_pointer(),
                                     no_jets_two_electrons.get_value_pointer(),
                                         no_jets_two_muons.get_value_pointer(),
                                        one_jet_no_leptons.get_value_pointer(),
                             no_jets_one_negative_electron.get_value_pointer(),
                             no_jets_one_positive_electron.get_value_pointer(),
                                 no_jets_one_negative_muon.get_value_pointer(),
                                 no_jets_one_positive_muon.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        two_jets_no_leptons.add_channel( channel_adder );
        OSSF_minus_OSDF_leptons.add_channel( channel_adder );
        no_jets_two_electrons.add_channel( channel_adder );
        no_jets_two_muons.add_channel( channel_adder );
        one_jet_no_leptons.add_channel( channel_adder );
        no_jets_one_negative_electron.add_channel( channel_adder );
        no_jets_one_positive_electron.add_channel( channel_adder );
        no_jets_one_negative_muon.add_channel( channel_adder );
        no_jets_one_positive_muon.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );

        channel_adder
        = new neutralino_three_body_to_taus_jj( given_primary_cut,
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
                                      no_jets_no_leptons.get_value_pointer() );
        two_jets_no_leptons.add_channel( channel_adder );
        one_jet_no_leptons.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );

        channel_adder
        = new neutralino_three_body_to_taus_jl( given_primary_cut,
                                                given_secondary_cut,
                                                given_jet_cut,
                                                given_kinematics,
                                                given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                                given_shortcuts,
                             one_jet_one_negative_electron.get_value_pointer(),
                             one_jet_one_positive_electron.get_value_pointer(),
                                 one_jet_one_negative_muon.get_value_pointer(),
                                 one_jet_one_positive_muon.get_value_pointer(),
                                        one_jet_no_leptons.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        one_jet_one_negative_electron.add_channel( channel_adder );
        one_jet_one_positive_electron.add_channel( channel_adder );
        one_jet_one_negative_muon.add_channel( channel_adder );
        one_jet_one_positive_muon.add_channel( channel_adder );
        one_jet_no_leptons.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );

        channel_adder
        = new neutralino_three_body_to_taus_ll( given_primary_cut,
                                                given_secondary_cut,
                                                given_jet_cut,
                                                given_kinematics,
                                                given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                                given_shortcuts,
                                     no_jets_two_electrons.get_value_pointer(),
           no_jets_one_negative_electron_one_positive_muon.get_value_pointer(),
           no_jets_one_positive_electron_one_negative_muon.get_value_pointer(),
                                         no_jets_two_muons.get_value_pointer(),
                             no_jets_one_negative_electron.get_value_pointer(),
                             no_jets_one_positive_electron.get_value_pointer(),
                                 no_jets_one_negative_muon.get_value_pointer(),
                                 no_jets_one_positive_muon.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        no_jets_two_electrons.add_channel( channel_adder );
        no_jets_one_negative_electron_one_positive_muon.add_channel(
                                                               channel_adder );
        no_jets_one_positive_electron_one_negative_muon.add_channel(
                                                               channel_adder );
        no_jets_two_muons.add_channel( channel_adder );
        no_jets_one_negative_electron.add_channel( channel_adder );
        no_jets_one_positive_electron.add_channel( channel_adder );
        no_jets_one_negative_muon.add_channel( channel_adder );
        no_jets_one_positive_muon.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );

      }
    else if( given_shortcuts->is_in( given_decaying_EWino->get_PDG_code(),
                                     given_shortcuts->get_charginos() ) )
      {

        cascade_acceptance_value* no_jets_one_signed_electron;
        cascade_acceptance_value* no_jets_one_signed_muon;

        if( given_EWino_is_not_antiparticle )
          {

            no_jets_one_signed_electron = no_jets_one_positive_electron;
            no_jets_one_signed_muon = no_jets_one_positive_muon;

          }
        else
          {

            no_jets_one_signed_electron = no_jets_one_negative_electron;
            no_jets_one_signed_muon = no_jets_one_negative_muon;

          }


        // debugging:
        /**std::cout
        << std::endl
        << "adding selectron_L channel";
        std::cout << std::endl;**/

        channel_adder = new chargino_to_semu( given_primary_cut,
                                              given_secondary_cut,
                                              given_kinematics,
                                              given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                              given_decaying_EWino,
                                              given_EWino_is_not_antiparticle,
                                            given_shortcuts->get_selectron_L(),
                                              given_shortcuts,
                              no_jets_one_signed_electron->get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        no_jets_one_signed_electron->add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );


        // debugging:
        /**std::cout
        << std::endl
        << "adding selectron_R channel";
        std::cout << std::endl;**/

        channel_adder = new chargino_to_semu( given_primary_cut,
                                              given_secondary_cut,
                                              given_kinematics,
                                              given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                              given_decaying_EWino,
                                              given_EWino_is_not_antiparticle,
                                            given_shortcuts->get_selectron_R(),
                                              given_shortcuts,
                              no_jets_one_signed_electron->get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        no_jets_one_signed_electron->add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );


        // debugging:
        /**std::cout
        << std::endl
        << "adding smuon_L channel";
        std::cout << std::endl;**/

        channel_adder = new chargino_to_semu( given_primary_cut,
                                              given_secondary_cut,
                                              given_kinematics,
                                              given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                              given_decaying_EWino,
                                              given_EWino_is_not_antiparticle,
                                              given_shortcuts->get_smuon_L(),
                                              given_shortcuts,
                                  no_jets_one_signed_muon->get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        no_jets_one_signed_muon->add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );


        // debugging:
        /**std::cout
        << std::endl
        << "adding smuon_R channel";
        std::cout << std::endl;**/

        channel_adder = new chargino_to_semu( given_primary_cut,
                                              given_secondary_cut,
                                              given_kinematics,
                                              given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                              given_decaying_EWino,
                                              given_EWino_is_not_antiparticle,
                                              given_shortcuts->get_smuon_R(),
                                              given_shortcuts,
                                  no_jets_one_signed_muon->get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        no_jets_one_signed_muon->add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );


        // debugging:
        /**std::cout
        << std::endl
        << "adding stau_1 channels";
        std::cout << std::endl;**/

        channel_adder = new chargino_to_stau_j( given_jet_cut,
                                                given_kinematics,
                                                given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                               given_shortcuts->get_stau_one(),
                                                given_shortcuts,
                                        one_jet_no_leptons.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        one_jet_no_leptons.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );

        channel_adder = new chargino_to_stau_l( given_primary_cut,
                                                given_secondary_cut,
                                                given_kinematics,
                                                given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                               given_shortcuts->get_stau_one(),
                                                given_shortcuts,
                              no_jets_one_signed_electron->get_value_pointer(),
                                  no_jets_one_signed_muon->get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        no_jets_one_signed_electron->add_channel( channel_adder );
        no_jets_one_signed_muon->add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );


        // debugging:
        /**std::cout
        << std::endl
        << "adding stau_2 channels";
        std::cout << std::endl;**/

        channel_adder = new chargino_to_stau_j( given_jet_cut,
                                                given_kinematics,
                                                given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                               given_shortcuts->get_stau_two(),
                                                given_shortcuts,
                                        one_jet_no_leptons.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        one_jet_no_leptons.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );

        channel_adder = new chargino_to_stau_l( given_primary_cut,
                                                given_secondary_cut,
                                                given_kinematics,
                                                given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                               given_shortcuts->get_stau_two(),
                                                given_shortcuts,
                              no_jets_one_signed_electron->get_value_pointer(),
                                  no_jets_one_signed_muon->get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        no_jets_one_signed_electron->add_channel( channel_adder );
        no_jets_one_signed_muon->add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );


        // debugging:
        /**std::cout
        << std::endl
        << "adding e_sneutrino_L channel";
        std::cout << std::endl;**/

        channel_adder = new chargino_to_emu_sneutrino( given_primary_cut,
                                                       given_secondary_cut,
                                                       given_kinematics,
                                                       given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                       given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                   given_shortcuts->get_electron_sneutrino_L(),
                                                       given_shortcuts,
                              no_jets_one_signed_electron->get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        no_jets_one_signed_electron->add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );


        // debugging:
        /**std::cout
        << std::endl
        << "adding e_sneutrino_R channel";
        std::cout << std::endl;**/

        channel_adder = new chargino_to_emu_sneutrino( given_primary_cut,
                                                       given_secondary_cut,
                                                       given_kinematics,
                                                       given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                       given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                   given_shortcuts->get_electron_sneutrino_R(),
                                                       given_shortcuts,
                              no_jets_one_signed_electron->get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        no_jets_one_signed_electron->add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );


        // debugging:
        /**std::cout
        << std::endl
        << "adding mu_sneutrino_L channel";
        std::cout << std::endl;**/

        channel_adder = new chargino_to_emu_sneutrino( given_primary_cut,
                                                       given_secondary_cut,
                                                       given_kinematics,
                                                       given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                       given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                       given_shortcuts->get_muon_sneutrino_L(),
                                                       given_shortcuts,
                                  no_jets_one_signed_muon->get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        no_jets_one_signed_muon->add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );


        // debugging:
        /**std::cout
        << std::endl
        << "adding mu_sneutrino_R channel";
        std::cout << std::endl;**/

        channel_adder = new chargino_to_emu_sneutrino( given_primary_cut,
                                                       given_secondary_cut,
                                                       given_kinematics,
                                                       given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                       given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                       given_shortcuts->get_muon_sneutrino_R(),
                                                       given_shortcuts,
                                  no_jets_one_signed_muon->get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        no_jets_one_signed_muon->add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );


        // debugging:
        /**std::cout
        << std::endl
        << "adding tau_sneutrino_L channel";
        std::cout << std::endl;**/

        channel_adder = new chargino_to_tau_sneutrino_j( given_jet_cut,
                                                         given_kinematics,
                                                       given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                         given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                        given_shortcuts->get_tau_sneutrino_L(),
                                                         given_shortcuts,
                                        one_jet_no_leptons.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        one_jet_no_leptons.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );

        channel_adder = new chargino_to_tau_sneutrino_l( given_primary_cut,
                                                         given_secondary_cut,
                                                         given_kinematics,
                                                       given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                         given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                        given_shortcuts->get_tau_sneutrino_L(),
                                                         given_shortcuts,
                              no_jets_one_signed_electron->get_value_pointer(),
                                  no_jets_one_signed_muon->get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        no_jets_one_signed_electron->add_channel( channel_adder );
        no_jets_one_signed_muon->add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );


        // debugging:
        /**std::cout
        << std::endl
        << "adding tau_sneutrino_R channel";
        std::cout << std::endl;**/

        channel_adder = new chargino_to_tau_sneutrino_j( given_jet_cut,
                                                         given_kinematics,
                                                       given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                         given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                        given_shortcuts->get_tau_sneutrino_R(),
                                                         given_shortcuts,
                                        one_jet_no_leptons.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        one_jet_no_leptons.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );

        channel_adder = new chargino_to_tau_sneutrino_l( given_primary_cut,
                                                         given_secondary_cut,
                                                         given_kinematics,
                                                       given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                         given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                        given_shortcuts->get_tau_sneutrino_R(),
                                                         given_shortcuts,
                              no_jets_one_signed_electron->get_value_pointer(),
                                  no_jets_one_signed_muon->get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        no_jets_one_signed_electron->add_channel( channel_adder );
        no_jets_one_signed_muon->add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );


        // debugging:
        /**std::cout
        << std::endl
        << "adding W channels";
        std::cout << std::endl;**/

        channel_adder = new chargino_to_W_without_tau( given_primary_cut,
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
                              no_jets_one_signed_electron->get_value_pointer(),
                                  no_jets_one_signed_muon->get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        two_jets_no_leptons.add_channel( channel_adder );
        one_jet_no_leptons.add_channel( channel_adder );
        no_jets_one_signed_electron->add_channel( channel_adder );
        no_jets_one_signed_muon->add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );

        channel_adder = new chargino_to_W_to_tau_j( given_jet_cut,
                                                    given_kinematics,
                                                    given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                    given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                                 given_shortcuts->get_W_plus(),
                                                    given_shortcuts,
                                        one_jet_no_leptons.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        one_jet_no_leptons.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );

        channel_adder = new chargino_to_W_to_tau_l( given_primary_cut,
                                                    given_secondary_cut,
                                                    given_kinematics,
                                                    given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                    given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                                 given_shortcuts->get_W_plus(),
                                                    given_shortcuts,
                                        one_jet_no_leptons.get_value_pointer(),
                              no_jets_one_signed_electron->get_value_pointer(),
                                  no_jets_one_signed_muon->get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        no_jets_one_signed_electron->add_channel( channel_adder );
        no_jets_one_signed_muon->add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );


        // debugging:
        /**std::cout
        << std::endl
        << "adding H channel";
        std::cout << std::endl;**/

        channel_adder = new chargino_to_Higgs_without_tau( given_primary_cut,
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
                                  no_jets_one_signed_muon->get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        two_jets_no_leptons.add_channel( channel_adder );
        one_jet_no_leptons.add_channel( channel_adder );
        no_jets_one_signed_muon->add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );

        channel_adder = new chargino_to_Higgs_to_tau_j( given_jet_cut,
                                                        given_kinematics,
                                                       given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                        given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                    given_shortcuts->get_charged_EWSB_scalar(),
                                                        given_shortcuts,
                                        one_jet_no_leptons.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        one_jet_no_leptons.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );

        channel_adder = new chargino_to_Higgs_to_tau_l( given_primary_cut,
                                                        given_secondary_cut,
                                                        given_kinematics,
                                                       given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                        given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                    given_shortcuts->get_charged_EWSB_scalar(),
                                                        given_shortcuts,
                                        one_jet_no_leptons.get_value_pointer(),
                              no_jets_one_signed_electron->get_value_pointer(),
                                  no_jets_one_signed_muon->get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        no_jets_one_signed_electron->add_channel( channel_adder );
        no_jets_one_signed_muon->add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );


        // debugging:
        /**std::cout
        << std::endl
        << "adding 3-body channel";
        std::cout << std::endl;**/

        channel_adder = new chargino_three_body_without_tau( given_primary_cut,
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
                              no_jets_one_signed_electron->get_value_pointer(),
                                  no_jets_one_signed_muon->get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        two_jets_no_leptons.add_channel( channel_adder );
        one_jet_no_leptons.add_channel( channel_adder );
        no_jets_one_signed_electron->add_channel( channel_adder );
        no_jets_one_signed_muon->add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );

        channel_adder = new chargino_three_body_to_tau_j( given_jet_cut,
                                                          given_kinematics,
                                                       given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                          given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                                          given_shortcuts,
                                        one_jet_no_leptons.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        one_jet_no_leptons.add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );

        channel_adder = new chargino_three_body_to_tau_l( given_primary_cut,
                                                          given_secondary_cut,
                                                          given_kinematics,
                                                       given_decaying_scolored,
                                            given_scolored_is_not_antiparticle,
                                                          given_decaying_EWino,
                                               given_EWino_is_not_antiparticle,
                                                          given_shortcuts,
                                        one_jet_no_leptons.get_value_pointer(),
                              no_jets_one_signed_electron->get_value_pointer(),
                                  no_jets_one_signed_muon->get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
        no_jets_one_signed_electron->add_channel( channel_adder );
        no_jets_one_signed_muon->add_channel( channel_adder );
        no_jets_no_leptons.add_channel( channel_adder );
        EW_decay_channels.push_back( channel_adder );

      }

  }




  fullCascade::fullCascade( colorfulCascadeType const typeOfColorfulCascade,
                            int const firstDecayBodyNumber,
                            double const chargeConjugateSumFactor ) :
    shortcut( NULL ),
    initialScolored( NULL ),
    beamEnergy( CppSLHA::CppSLHA_global::really_wrong_value ),
    typeOfColorfulCascade( typeOfColorfulCascade ),
    firstDecayBodyNumber( firstDecayBodyNumber ),
    chargeConjugateSumFactor( chargeConjugateSumFactor ),
    subcascade( NULL ),
    cascadeDefiner(),
    ewinoCascades( NULL ),
    vectorCascades( NULL ),
    soughtDecayProductList( 1,
                            CppSLHA::CppSLHA_global::really_wrong_value )
  {
    // just an initialization list.
  }

  fullCascade::~fullCascade()
  {
    // does nothing.
  }

  double
  fullCascade::specifiedJetsOneOssfMinusOsdfPair( acceptanceCutSet* const cuts,
                                             int const numberOfAdditionalJets )
  /* this calls
   * ewinoCascade->getOssfMinusOsdf( scoloredIsNotAntiparticle,
   *                                 cuts ),
   * & if vectorCascade is not NULL, it calls
   * ewinoCascade->getAcceptance( [ { 0, 1, 2 } jets + 0 leptons for cuts ] )
   * & the same for vectorCascade, as well as ->getOssfMinusOsdf( cuts ).
   */
  {
    // we should really check that ewinoCascade is not NULL, but this function
    // should never be called on a fullCascade which has a NULL ewinoCascade.
    if( NULL == vectorCascades )
    {
      /* if there is no possibility of extra jets from a vector decay, we can
       * only get a single OSSF-OSDF pair from an electroweakino decay, which
       * cannot produce additional jets.
       */
      if( 0 == numberOfAdditionalJets )
      {
        return ewinoCascades->getOssfMinusOsdf( cuts );
      }
      else
      {
        return 0.0;
      }
    }
    else
    {
      return ( ewinoCascades->getOssfMinusOsdf( cuts )
               * vectorCascades->getAcceptance( true
                                             /* doesn't matter for 0 leptons*/,
                                               cuts,
                                               numberOfAdditionalJets,
                                               0,
                                               0,
                                               0,
                                               0 )
               + vectorCascades->getOssfMinusOsdf( cuts )
                 * ewinoCascades->getAcceptance( true
                                              /* doesn't matter for 0 leptons*/,
                                                cuts,
                                                numberOfAdditionalJets,
                                                0,
                                                0,
                                                0,
                                                0 ) );
    }
  }

  double
  fullCascade::unspecifiedJetsSpecifiedOssfMinusOsdfPairs(
                                                  acceptanceCutSet* const cuts,
                                                int const numberOfLeptonPairs )
  // this only calls ewinoCascade->getOssfMinusOsdf( cuts ),
  // & if vectorCascade is not NULL, vectorCascade->getOssfMinusOsdf( cuts ).
  {
    if( 0 == numberOfLeptonPairs )
    {
      return unspecifiedJetsSpecifiedChargeSummedLeptons( cuts,
                                                          0 );
    }
    else
    {
      // we should check that ewinoCascade is not NULL, but this function
      // should never be called on a fullCascade which has a NULL ewinoCascade.
      if( NULL == vectorCascades )
      {
        /* if there is no possibility of extra jets from a vector decay, we can
         * only get a single OSSF-OSDF pair from an electroweakino decay, which
         * cannot produce additional jets.
         */
        if( 1 == numberOfLeptonPairs )
        {
          return ewinoCascades->getOssfMinusOsdf( cuts );
        }
        else
          // if there's no vector cascade, the only possibilities for non-zero
          // results are 0 pairs or 1 pair.
        {
          return 0.0;
        }
      }
      else
      {
        if( 2 == numberOfLeptonPairs )
        {
          return ( ewinoCascades->getOssfMinusOsdf( cuts )
                   * vectorCascades->getOssfMinusOsdf( cuts ) );
        }
        else if( 1 == numberOfLeptonPairs )
        {
        return ( ewinoCascades->getOssfMinusOsdf( cuts )
                 * ( vectorCascades->getAcceptance( true
                                            /* doesn't matter for 0 leptons */,
                                                 cuts,
                                                 0,
                                                 0,
                                                 0,
                                                 0,
                                                 0 )
                     + vectorCascades->getAcceptance( true
                                            /* doesn't matter for 0 leptons */,
                                                     cuts,
                                                     1,
                                                     0,
                                                     0,
                                                     0,
                                                     0 )
                     + vectorCascades->getAcceptance( true
                                            /* doesn't matter for 0 leptons */,
                                                     cuts,
                                                     2,
                                                     0,
                                                     0,
                                                     0,
                                                     0 ) )
                 + vectorCascades->getOssfMinusOsdf( cuts )
                   * ( ewinoCascades->getAcceptance( true
                                            /* doesn't matter for 0 leptons */,
                                                    cuts,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0 )
                       + ewinoCascades->getAcceptance( true
                                            /* doesn't matter for 0 leptons */,
                                                      cuts,
                                                      1,
                                                      0,
                                                      0,
                                                      0,
                                                      0 )
                       + ewinoCascades->getAcceptance( true
                                            /* doesn't matter for 0 leptons */,
                                                      cuts,
                                                      2,
                                                      0,
                                                      0,
                                                      0,
                                                      0 ) ) );
        }
        else
          /* the case of 0 pairs has already been checked for, so by this
           * point, we've discovered that we were asked for a number of pairs
           * that is not 0, 1, or 2. thus we return 0.0 because we assume that
           * any cascades which could result in so many leptons are negligible.
           */
        {
          return 0.0;
        }
      }  // end of if vectorCascade is NULL or not.
    }  // end of if numberOfLeptonPairs == 0 or not.
  }

  double
  fullCascade::specifiedJetsSpecifiedChargeSummedLeptons(
                                                  acceptanceCutSet* const cuts,
                                              int const numberOfAdditionalJets,
                                                    int const numberOfLeptons )
  {
    if( ( 0 > ( numberOfLeptons + numberofAdditionalJets ) )
        ||
        ( numberOfSmFermionsFromElectroweakDecaysPerFullCascade
          < ( numberOfLeptons + numberofAdditionalJets ) ) )
    {
      return 0.0;
    }
    else
    {
      double returnDouble( 0.0 );
      for( int numberOfNegativeElectrons = numberOfLeptons;
           0 <= numberOfNegativeElectrons;
           --numberOfNegativeElectrons )
      {
        for( int numberOfPositiveElectrons
             = ( numberOfLeptons - numberOfNegativeElectrons );
             0 <= numberOfPositiveElectrons;
             --numberOfPositiveElectrons )
        {
          for( int numberOfNegativeMuons
               = ( numberOfLeptons - numberOfNegativeElectrons
                                   - numberOfPositiveElectrons );
               0 <= numberOfNegativeMuons;
               --numberOfNegativeMuons )
          {
            returnDouble += getAcceptance( true
                                /* doesn't matter for charge-summed leptons */,
                                           cuts,
                                           numberOfAdditionalJets,
                                           numberOfNegativeElectrons,
                                           numberOfPositiveElectrons,
                                           numberOfNegativeMuons,
                                           ( numberOfLeptons
                                             - numberOfNegativeElectrons
                                             - numberOfPositiveElectrons
                                             - numberOfNegativeMuons ) );
          }  // end of loop over negative muons.
        }  // end of loop over positive electrons.
      }  // end of loop over negative electrons.
      return returnDouble;
    }  // end of if numberOfLeptons was in the allowed range.
  }

  double
  fullCascade::unspecifiedJetsSpecifiedChargeSummedLeptons(
                                                  acceptanceCutSet* const cuts,
                                                    int const numberOfLeptons )
  {
    if( ( 0 < numberOfLeptons )
        ||
        ( numberOfSmFermionsFromElectroweakDecaysPerFullCascade
          > numberOfLeptons ) )
    {
      return 0.0;
    }
    else
    {
      double returnDouble( 0.0 );
      for( int numberOfJets
           = ( numberOfSmFermionsFromElectroweakDecaysPerFullCascade
               - numberOfLeptons );
           0 <= numberOfJets;
           --numberOfJets )
      {
        returnDouble += specifiedJetsSpecifiedChargeSummedLeptons( cuts,
                                                                  numberOfJets,
                                                             numberOfLeptons );
      }  // end of loop over jets.
      return returnDouble;
    }  // end of if numberOfLeptons was in the allowed range.
  }

  double
  fullCascade::getCombinedAcceptance(
                                    bool const vectorScoloredIsNotAntiparticle,
                                     bool const ewinoScoloredIsNotAntiparticle,
                                      acceptanceCutSet* const cuts,
                                      int const numberOfAdditionalJets,
                                      int const numberOfNegativeElectrons,
                                      int const numberOfPositiveElectrons,
                                      int const numberOfNegativeMuons,
                                      int const numberOfPositiveMuons )
  /* this combines the acceptances from ewinoCascade & vectorCascade. Whether
   * the relevant scoloreds are particles or antiparticles should be decided
   * by the derived class using this function.
   */
  {
    if( ( 0 > ( numberOfAdditionalJets
                + numberOfNegativeElectrons + numberOfPositiveElectrons
                + numberOfNegativeMuons + numberOfPositiveMuons ) )
        ||
        ( numberOfSmFermionsFromElectroweakDecaysPerFullCascade
          < ( numberOfAdditionalJets
              + numberOfNegativeElectrons + numberOfPositiveElectrons
              + numberOfNegativeMuons + numberOfPositiveMuons ) ) )
    {
      return 0.0;
    }
    else
    {
      double returnDouble( 0.0 );
      double ewinoAcceptance;
      double vectorAcceptance;
      for( int ewinoJets = numberOfAdditionalJets;
           0 <= ewinoJets;
           --ewinoJets )
      {
        for( int ewinoNegativeElectrons = numberOfNegativeElectrons;
             0 <= ewinoNegativeElectrons;
             --ewinoNegativeElectrons )
        {
          for( int ewinoPositiveElectrons = numberOfPositiveElectrons;
               0 <= ewinoPositiveElectrons;
               --ewinoPositiveElectrons )
          {
            for( int ewinoNegativeMuons = numberOfNegativeMuons;
                 0 <= ewinoNegativeMuons;
                 --ewinoNegativeMuons )
            {
              for( int ewinoPositiveMuons = numberOfPositiveMuons;
                   0 <= ewinoPositiveMuons;
                   --ewinoPositiveMuons )
              {
                if( ewinoScoloredIsNotAntiparticle )
                {
                  ewinoAcceptance = ewinoCascades->getAcceptance( cuts,
                                                                 ewinoJets,
                                                        ewinoNegativeElectrons,
                                                        ewinoPositiveElectrons,
                                                            ewinoNegativeMuons,
                                                          ewinoPositiveMuons );
                }
                else
                {
                  ewinoAcceptance = ewinoCascades->getAcceptance( cuts,
                                                                 ewinoJets,
                                                        ewinoPositiveElectrons,
                                                        ewinoNegativeElectrons,
                                                            ewinoPositiveMuons,
                                                          ewinoNegativeMuons );
                }
                if( vectorScoloredIsNotAntiparticle )
                {
                  vectorAcceptance = vectorCascades->getAcceptance( cuts,
                                        ( numberOfAdditionalJets - ewinoJets ),
                        ( numberOfNegativeElectrons - ewinoNegativeElectrons ),
                        ( numberOfPositiveElectrons - ewinoPositiveElectrons ),
                                ( numberOfNegativeMuons - ewinoNegativeMuons ),
                              ( numberOfPositiveMuons - ewinoPositiveMuons ) );
                }
                else
                {
                  vectorAcceptance = vectorCascades->getAcceptance( cuts,
                                        ( numberOfAdditionalJets - ewinoJets ),
                        ( numberOfPositiveElectrons - ewinoPositiveElectrons ),
                        ( numberOfNegativeElectrons - ewinoNegativeElectrons ),
                                ( numberOfPositiveMuons - ewinoPositiveMuons ),
                              ( numberOfNegativeMuons - ewinoNegativeMuons ) );
                }
              returnDouble += ( ewinoAcceptance * vectorAcceptance );
              }  // end of loop over positive electrons.
            }  // end of loop over negative muons.
          }  // end of loop over positive electrons.
        }  // end of loop over negative electrons.
      }  // end of loop over jets.
      return returnDouble;
    }  // end of if numberOfLeptons was in the allowed range.
  }



  sxFullCascade::sxFullCascade() :
    fullCascade( sx,
                 2,
                 1.0 )
  {
    // just an initialization list.
  }

  sxFullCascade::~sxFullCascade()
  {
    // does nothing.
  }



  gxFullCascade::gxFullCascade() :
    fullCascade( gx,
                 3,
                 2.0 )
  {
    // just an initialization list.
  }

  gxFullCascade::~gxFullCascade()
  {
    // does nothing.
  }

  double
  gxFullCascade::getAcceptance( bool const scoloredIsNotAntiparticle,
                                acceptanceCutSet* const cuts,
                                int const numberOfAdditionalJets,
                                int const numberOfNegativeElectrons,
                                int const numberOfPositiveElectrons,
                                int const numberOfNegativeMuons,
                                int const numberOfPositiveMuons )
  /* this calls the appropriate functions on ewinoCascade, summing over the 2
   * charge versions if the electroweakino is a chargino. it ignores the value
   * of scoloredIsNotAntiparticle.
   */
  {
    if( ewinoCascades->getElectroweakDecayer()->counts_as_self_conjugate() )
      // if the gluino has only 1 charge version of this decay...
    {
      return ewinoCascades->getAcceptance( cuts,
                                          numberOfAdditionalJets,
                                          numberOfNegativeElectrons,
                                          numberOfPositiveElectrons,
                                          numberOfNegativeMuons,
                                          numberOfPositiveMuons );
      // the summing over charge-conjugate decays is already included in the
      // BR to the electroweakino.
    }
    else
      // otherwise if the gluino has 2 charge versions of this decay...
    {
      return ( 0.5 * ( ewinoCascades->getAcceptance( cuts,
                                                    numberOfAdditionalJets,
                                                    numberOfNegativeElectrons,
                                                    numberOfPositiveElectrons,
                                                    numberOfNegativeMuons,
                                                    numberOfPositiveMuons )
                       + ewinoCascades->getAcceptance( cuts,
                                                      numberOfAdditionalJets,
                                                     numberOfPositiveElectrons,
                                                     numberOfNegativeElectrons,
                                                      numberOfPositiveMuons,
                                                   numberOfNegativeMuons ) ) );
      // the summing over charge-conjugate decays is already included in the BR
      // to the electroweakino.
    }
  }



  sjgxFullCascade::sjgxFullCascade() :
    fullCascade( sjgx,
                 2,
                 1.0 )
  // sjgx also means that the initial decay is 2-body.
  {
    // just an initialization list.
  }

  sjgxFullCascade::~sjgxFullCascade()
  {
    // does nothing.
  }



  gjsxFullCascade::gjsxFullCascade() :
    fullCascade( gjsx,
                 2,
                 2.0 )
  // gjsx also means that the initial decay is 2-body.
  {
    // just an initialization list.
  }

  gjsxFullCascade::~gjsxFullCascade()
  {
    // does nothing.
  }



  sjgjsxFullCascade::sjgjsxFullCascade() :
    fullCascade( sjgjsx,
                 2,
                 1.0 )
  // sjgjsx also means that the initial decay is 2-body.
  {
    // just an initialization list.
  }

  sjgjsxFullCascade::~sjgjsxFullCascade()
  {
    // does nothing.
  }



  svsxFullCascade::svsxFullCascade() :
    fullCascade( svsx,
                 2,
                 1.0 )
  // svsx also means that the initial decay is 2-body.
  {
    soughtDecayProductList.push_back(
                                 CppSLHA::CppSLHA_global::really_wrong_value );
    // a cascade beginning with sv needs to look for 2 specific decay products.
  }

  svsxFullCascade::~svsxFullCascade()
  {
    // does nothing.
  }



  gvsxFullCascade::gvsxFullCascade() :
    fullCascade( gvsx,
                 3,
                 2.0 )
  // gvsx also means that the initial decay is 3-body.
  {
    soughtDecayProductList.push_back(
                                 CppSLHA::CppSLHA_global::really_wrong_value );
    // a cascade beginning with gv needs to look for 2 specific decay products.
  }

  gvsxFullCascade::~gvsxFullCascade()
  {
    // does nothing.
  }



  gjsvsxFullCascade::gjsvsxFullCascade() :
    fullCascade( gjsvsx,
                 2,
                 2.0 )
  // gjsvsx also means that the initial decay is 2-body.
  {
    // just an initialization list.
  }

  gjsvsxFullCascade::~gjsvsxFullCascade()
  {
    // does nothing.
  }



  svgxFullCascade::svgxFullCascade() :
    fullCascade( svgx,
                 3,
                 1.0 )
  // svgx also means that the initial decay is 3-body.
  {
    // just an initialization list.
  }

  svgxFullCascade::~svgxFullCascade()
  {
    // does nothing.
  }



  svsjgxFullCascade::svsjgxFullCascade() :
    fullCascade( svsjgx,
                 2,
                 1.0 )
  // svsjgx also means that the initial decay is 2-body.
  {
    // just an initialization list.
  }

  svsjgxFullCascade::~svsjgxFullCascade()
  {
    // does nothing.
  }



  svsjgjsxFullCascade::svsjgjsxFullCascade() :
    fullCascade( svsjgjsx,
                 2,
                 1.0 )
  // svsjgjsx also means that the initial decay is 2-body.
  {
    // just an initialization list.
  }

  svsjgjsxFullCascade::~svsjgjsxFullCascade()
  {
    // does nothing.
  }



  svgjsxFullCascade::svgjsxFullCascade() :
    fullCascade( svgjsx,
                 3,
                 1.0 )
  // svgjsx also means that the initial decay is 3-body.
  {
    // just an initialization list.
  }

  svgjsxFullCascade::~svgjsxFullCascade()
  {
    // does nothing.
  }



  sjgvsxFullCascade::sjgvsxFullCascade() :
    fullCascade( sjgvsx,
                 2,
                 1.0 )
  // sjgvsx also means that the initial decay is 2-body.
  {
    // just an initialization list.
  }

  sjgvsxFullCascade::~sjgvsxFullCascade()
  {
    // does nothing.
  }



  sjgjsvsxFullCascade::sjgjsvsxFullCascade() :
    fullCascade( sjgjsvsx,
                 2,
                 1.0 )
  // sjgjsvsx also means that the initial decay is 2-body.
  {
    // just an initialization list.
  }

  sjgjsvsxFullCascade::~sjgjsvsxFullCascade()
  {
    // does nothing.
  }



  fullCascadeSet::fullCascadeSet( input_handler const* const shortcut,
                   CppSLHA::particle_property_set const* const initialScolored,
                   electroweakCascadeHandler* const electroweakCascadeHandling,
                      std::list< fullCascadeSet* >* const squarkCascadeSetList,
                                  double const beamEnergy ) :
    readied_for_new_point( shortcut->get_readier() ),
    initialScolored( initialScolored ),
    sxCascades( &fullCascade::getNewFullCascade() ),
    gxCascades( &fullCascade::getNewFullCascade() ),
    sjgxCascades( &fullCascade::getNewFullCascade() ),
    gjsxCascades( &fullCascade::getNewFullCascade() ),
    sjgjsxCascades( &fullCascade::getNewFullCascade() ),
    svsxCascades( &fullCascade::getNewFullCascade() ),
    gvsxCascades( &fullCascade::getNewFullCascade() ),
    gjsvsxCascades( &fullCascade::getNewFullCascade() ),
    svgxCascades( &fullCascade::getNewFullCascade() ),
    svsjgxCascades( &fullCascade::getNewFullCascade() ),
    svsjgjsxCascades( &fullCascade::getNewFullCascade() ),
    svgjsxCascades( &fullCascade::getNewFullCascade() ),
    sjgvsxCascades( &fullCascade::getNewFullCascade() ),
    sjgjsvsxCascades( &fullCascade::getNewFullCascade() ),
    shortcut( shortcut ),
    squarkCascadeSetList( squarkCascadeSetList ),
    //squarkCascadeSetListNotYetOrdered( true ),
    beamEnergy( beamEnergy )
  {
    // just an initialization list.
  }

  fullCascadeSet::~fullCascadeSet()
  {
    // does nothing.
  }

  void
  fullCascadeSet::setUpCascades()
  // this clears openCascades, then sets it to be filled with all open
  // fullCascades that initialScolored has for this point.
  {
    openCascades.clear();
    for( std::vector< sxFullCascade* >::iterator
         sxIterator = sxCascades->getVector()->begin();
         sxCascades->getVector()->end() > sxIterator;
         ++sxIterator )
    {
      if( (*sxIterator)->isOpen() )
      {
        openCascades.push_back( *sxIterator );
      }
    }
    for( std::vector< gxFullCascade* >::iterator
         gxIterator = gxCascades->getVector()->begin();
         gxCascades->getVector()->end() > gxIterator;
         ++gxIterator )
    {
      if( (*gxIterator)->isOpen() )
      {
        openCascades.push_back( *gxIterator );
      }
    }
    buildLongerCascades();
  }



  squarkFullCascadeSet::squarkFullCascadeSet(
                                           input_handler const* const shortcut,
                   electroweakCascadeHandler* const electroweakCascadeHandling,
                   CppSLHA::particle_property_set const* const initialScolored,
                      std::list< fullCascadeSet* >* const squarkCascadeSetList,
                                        fullCascadeSet* const gluinoCascadeSet,
                                              double const beamEnergy ) :
    fullCascadeSet( shortcut,
                    initialScolored,
                    electroweakCascadeHandling,
                    squarkCascadeSetList,
                    beamEnergy ),
    gluinoCascadeSet( gluinoCascadeSet )
  {
    // we have to set up sxCascades now:
    for( std::vector< CppSLHA::particle_property_set const* >::iterator
         ewinoIterator = shortcut->get_electroweakinos()->begin();
         shortcut->get_unstable_electroweakinos()->end() > ewinoIterator;
         ++ewinoIterator )
    {
      sxCascades.addNewAtEnd()->setProperties( shortcut,
                                               initialScolored,
                                               beamEnergy,
                                     electroweakCascadeHandling->getCascadeSet(
                                                               initialScolored,
                                                             ewinoIterator ) );
    }
  }

  squarkFullCascadeSet::~squarkFullCascadeSet()
  {
    // does nothing.
  }

  void
  squarkFullCascadeSet::buildLongerCascades()
  {
    squarkCascadeSetList->sort( &massOrdered );
    // now squarkCascadeSetList has been ordered such that the cascades of the
    // lighter squarks come before those of the heavier squarks.

    sjgxCascades.clear();
    sjgjsxCascades.clear();
    sjgvsxCascades.clear();
    sjgjsvsxCascades.clear();
    svgxCascades.clear();
    svgjsxCascades.clear();
    if( initialScolored->get_absolute_mass()
        > shortcut->get_gluino()->get_absolute_mass() )
    {
      for( std::vector< gxFullCascade* >::iterator
           gxIterator
           = gluinoCascadeSet->getGxCascades()->getVector()->begin();
           gluinoCascadeSet->getGxCascades()->getVector()->end() > gxIterator;
           ++gxIterator )
      {
        sjgxCascades.addNewAtEnd()->setProperties( initialScolored,
                                                   *gxIterator );
        openCascades.push_back( sjgxCascades.getBack() );
      }
      for( std::vector< gjsxFullCascade* >::iterator
           gjsxIterator
           = gluinoCascadeSet->getGjsxCascades()->getVector()->begin();
           gluinoCascadeSet->getGjsxCascades()->getVector()->end()
           > gjsxIterator;
           ++gjsxIterator )
      {
        sjgjsxCascades.addNewAtEnd()->setProperties( initialScolored,
                                                     *gjsxIterator );
        openCascades.push_back( sjgjsxCascades.getBack() );
      }
      for( std::vector< gvsxFullCascade* >::iterator
           gvsxIterator
           = gluinoCascadeSet->getGvsxCascades()->getVector()->begin();
           gluinoCascadeSet->getGvsxCascades()->getVector()->end()
           > gvsxIterator;
           ++gvsxIterator )
      {
        sjgvsxCascades.addNewAtEnd()->setProperties( initialScolored,
                                                     *gvsxIterator );
        openCascades.push_back( sjgvsxCascades.getBack() );
      }
      for( std::vector< gjsvsxFullCascade* >::iterator
           gjsvsxIterator
           = gluinoCascadeSet->getGjsvsxCascades()->getVector()->begin();
           gluinoCascadeSet->getGjsvsxCascades()->getVector()->end()
           > gjsvsxIterator;
           ++gjsvsxIterator )
      {
        sjgjsvsxCascades.addNewAtEnd()->setProperties( initialScolored,
                                                       *gjsvsxIterator );
        openCascades.push_back( sjgjsvsxCascades.getBack() );
      }
      for( std::vector< CppSLHA::particle_property_set const* >::iterator
           ewIterator = shortcut->get_EW_veved_and_vector_bosons()->begin();
           shortcut->get_EW_veved_and_vector_bosons()->end() > ewIterator;
           ++ewIterator )
      {
        soughtPositivePdgCodeList.clear();
        soughtPositivePdgCodeList.push_back( (*ewIterator)->get_PDG_code() );
        soughtPositivePdgCodeList.push_back( CppSLHA::PDG_code::gluino );
        soughtNegativePdgCodeList.clear();
        soughtNegativePdgCodeList.push_back(
                                            -((*ewIterator)->get_PDG_code()) );
        soughtNegativePdgCodeList.push_back( CppSLHA::PDG_code::gluino );
        if( ( initialScolored->get_absolute_mass()
              > ( shortcut->get_gluino()->get_absolute_mass()
                  + (*ewIterator)->get_absolute_mass() ) )
            &&
            ( ( LHC_FASER_global::negligible_BR
                < initialScolored->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                                 &soughtPositivePdgCodeList ) )
              || ( LHC_FASER_global::negligible_BR
                               < initialScolored->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                             &soughtNegativePdgCodeList ) ) ) )
        {
          for( std::vector< gxFullCascade* >::iterator
               gxIterator
               = gluinoCascadeSet->getGxCascades()->getVector()->begin();
               gluinoCascadeSet->getGxCascades()->getVector()->end()
               > gxIterator;
               ++gxIterator )
          {
            svgxCascades.addNewAtEnd()->setProperties( initialScolored,
                                     electroweakCascadeHandling->getCascadeSet(
                                                               initialScolored,
                                                                  *ewIterator),
                                                       *gxIterator );
            openCascades.push_back( svgxCascades.getBack() );
          }
          for( std::vector< gxFullCascade* >::iterator
               gjsxIterator
               = gluinoCascadeSet->getGjsxCascades()->getVector()->begin();
               gluinoCascadeSet->getGjsxCascades()->getVector()->end()
               > gjsxIterator;
               ++gjsxIterator )
          {
            svgjsxCascades.addNewAtEnd()->setProperties( initialScolored,
                                     electroweakCascadeHandling->getCascadeSet(
                                                               initialScolored,
                                                                  *ewIterator),
                                                         *gjsxIterator );
            openCascades.push_back( svgjsxCascades.getBack() );
          }
        }
      }  // end of loop over if svg is open
    }  // end of if squark mass is greater than gluino mass

    // now we look for decays into lighter squarks:
    svsxCascades.clear();
    svsjgxCascades.clear();
    svsjgjsxCascades.clear();
    for( std::list< fullCascadeSet* >::iterator
         cascadeIterator = squarkCascadeSetList->begin();
         (*cascadeIterator)->initialScolored->get_absolute_mass()
         < initialScolored->get_absolute_mass();
         // this conditional should prevent the iterator falling off the end
         // of the list provided that the list is as it should be.
         ++cascadeIterator )
    {
      for( std::vector< CppSLHA::particle_property_set const* >::iterator
           ewIterator = shortcut->get_EW_veved_and_vector_bosons()->begin();
           shortcut->get_EW_veved_and_vector_bosons()->end() > ewIterator;
           ++ewIterator )
      {
        soughtPositivePdgCodeList.clear();
        soughtPositivePdgCodeList.push_back(
                         (*cascadeIterator)->initialScolored->get_PDG_code() );
        soughtPositivePdgCodeList.push_back( (*ewIterator)->get_PDG_code() );
        soughtNegativePdgCodeList.clear();
        soughtNegativePdgCodeList.push_back(
                         (*cascadeIterator)->initialScolored->get_PDG_code() );
        soughtNegativePdgCodeList.push_back(
                                            -((*ewIterator)->get_PDG_code()) );
        if( ( initialScolored->get_absolute_mass()
              > ( (*cascadeIterator)->initialScolored->get_absolute_mass()
                  + (*ewIterator)->get_absolute_mass() ) )
            &&
            ( ( LHC_FASER_global::negligible_BR
                < initialScolored->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                                 &soughtPositivePdgCodeList ) )
              || ( LHC_FASER_global::negligible_BR
                               < initialScolored->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                             &soughtNegativePdgCodeList ) ) ) )
        {
          for( std::vector< sxFullCascade* >::iterator
               sxIterator
               = (*cascadeIterator)->getSxCascades()->getVector()->begin();
              (*cascadeIterator)->getSxCascades()->getVector()->end()
               > sxIterator;
               ++sxIterator )
          {
            svsxCascades.addNewAtEnd()->setProperties( initialScolored,
                                     electroweakCascadeHandling->getCascadeSet(
                                                               initialScolored,
                                                                  *ewIterator),
                                                       *sxIterator );
            openCascades.push_back( svsxCascades.getBack() );
          }
          for( std::vector< sjgxFullCascade* >::iterator
               sjgxIterator
               = (*cascadeIterator)->getSjgxCascades()->getVector()->begin();
              (*cascadeIterator)->getSjgxCascades()->getVector()->end()
               > sjgxIterator;
               ++sjgxIterator )
          {
            svsjgxCascades.addNewAtEnd()->setProperties( initialScolored,
                                     electroweakCascadeHandling->getCascadeSet(
                                                               initialScolored,
                                                                  *ewIterator),
                                                         *sjgxIterator );
            openCascades.push_back( svsjgxCascades.getBack() );
          }
          for( std::vector< sjgxFullCascade* >::iterator
               sjgjsxIterator
               = (*cascadeIterator)->getSjgjsxCascades()->getVector()->begin();
              (*cascadeIterator)->getSjgjsxCascades()->getVector()->end()
               > sjgjsxIterator;
               ++sjgjsxIterator )
          {
            svsjgjsxCascades.addNewAtEnd()->setProperties( initialScolored,
                                     electroweakCascadeHandling->getCascadeSet(
                                                               initialScolored,
                                                                  *ewIterator),
                                                           *sjgjsxIterator );
            openCascades.push_back( svsjgjsxCascades.getBack() );
          }
        }
      }  // end of loop over if svs is open
    }
  }



  gluinoFullCascadeSet::gluinoFullCascadeSet(
                                           input_handler const* const shortcut,
                   electroweakCascadeHandler* const electroweakCascadeHandling,
                      std::list< fullCascadeSet* >* const squarkCascadeSetList,
                                               double const beamEnergy ) :
    fullCascadeSet( shortcut,
                    shortcut->get_gluino(),
                    electroweakCascadeHandling,
                    squarkCascadeSetList,
                    beamEnergy )
  {
    // we have to set up gxCascades now:
    for( std::vector< CppSLHA::particle_property_set const* >::iterator
         ewinoIterator = shortcut->get_electroweakinos()->begin();
         shortcut->get_unstable_electroweakinos()->end() > ewinoIterator;
         ++ewinoIterator )
    {
      gxCascades.addNewAtEnd()->setProperties( shortcut,
                                               beamEnergy,
                                     electroweakCascadeHandling->getCascadeSet(
                                                               initialScolored,
                                                             ewinoIterator ) );
    }
  }

  gluinoFullCascadeSet::~gluinoFullCascadeSet()
  {
    // does nothing.
  }

  void
  gluinoFullCascadeSet::buildLongerCascades()
  {
    squarkCascadeSetList->sort( &massOrdered );
    // now squarkCascadeSetList has been ordered such that the cascades of the
    // lighter squarks come before those of the heavier squarks.

    gjsxCascades.clear();
    gvsxCascades.clear();
    gjsvsxCascades.clear();
    for( std::list< fullCascadeSet* >::iterator
         cascadeIterator = squarkCascadeSetList->begin();
         ( ( squarkCascadeSetList->end() != cascadeIterator )
           &&
           (*cascadeIterator)->initialScolored->get_absolute_mass()
           < initialScolored->get_absolute_mass() );
         ++cascadeIterator )
    {
      // each of the squarks looked at in this loop are lighter than the
      // gluino, so these cascades should be built:
      for( std::vector< sxFullCascade* >::iterator
           sxIterator
           = (*cascadeIterator)->getSxCascades()->getVector()->begin();
           (*cascadeIterator)->getSxCascades()->getVector()->end()
           > sxIterator;
           ++sxIterator )
      {
        gjsxCascades.addNewAtEnd()->setProperties( *sxIterator );
        openCascades.push_back( gjsxCascades.getBack() );
      }
      for( std::vector< svsxFullCascade* >::iterator
           svsxIterator
           = (*cascadeIterator)->getSvsxCascades()->getVector()->begin();
           (*cascadeIterator)->getSvsxCascades()->getVector()->end()
           > svsxIterator;
           ++svsxIterator )
      {
        gjsvsxCascades.addNewAtEnd()->setProperties( *svsxIterator );
        openCascades.push_back( gjsvsxCascades.getBack() );
      }

      // now we check for gvs:
      for( std::vector< CppSLHA::particle_property_set const* >::iterator
           ewIterator = shortcut->get_EW_veved_and_vector_bosons()->begin();
           shortcut->get_EW_veved_and_vector_bosons()->end() > ewIterator;
           ++ewIterator )
      {
        soughtPositivePdgCodeList.clear();
        soughtPositivePdgCodeList.push_back(
                         (*cascadeIterator)->initialScolored->get_PDG_code() );
        soughtPositivePdgCodeList.push_back( (*ewIterator)->get_PDG_code() );
        soughtNegativePdgCodeList.clear();
        soughtNegativePdgCodeList.push_back(
                         (*cascadeIterator)->initialScolored->get_PDG_code() );
        soughtNegativePdgCodeList.push_back(
                                            -((*ewIterator)->get_PDG_code()) );
        if( ( initialScolored->get_absolute_mass()
              > ( (*cascadeIterator)->initialScolored->get_absolute_mass()
                  + (*ewIterator)->get_absolute_mass() ) )
            &&
            ( ( LHC_FASER_global::negligible_BR
                < initialScolored->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                                 &soughtPositivePdgCodeList ) )
              || ( LHC_FASER_global::negligible_BR
                               < initialScolored->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                             &soughtNegativePdgCodeList ) ) ) )
        {
          for( std::vector< sxFullCascade* >::iterator
               sxIterator
               = (*cascadeIterator)->getSxCascades()->getVector()->begin();
              (*cascadeIterator)->getSxCascades()->getVector()->end()
               > sxIterator;
               ++sxIterator )
          {
            gvsxCascades.addNewAtEnd()->setProperties(
                                     electroweakCascadeHandling->getCascadeSet(
                                                               initialScolored,
                                                                  *ewIterator),
                                                       *sxIterator );
            openCascades.push_back( gvsxCascades.getBack() );
          }
        }
      }  // end of loop over if gvs is open
    }  // end of loop over squarkFullCascadeSets
  }
}  // end of LHC_FASER namespace.
