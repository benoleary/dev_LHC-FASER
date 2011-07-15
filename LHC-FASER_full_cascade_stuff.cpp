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
    << *(given_decaying_scolored->get_name())
    << ", "
    << *(given_decaying_EWino->get_name())
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

  EW_decayer_values::~EW_decayer_values()
  {

    for( std::vector< channel_calculator* >::iterator
         deletion_iterator = EW_decay_channels.begin();
         EW_decay_channels.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

  }



  full_cascade::full_cascade(
            CppSLHA::particle_property_set const* const given_initial_scolored,
                  CppSLHA::particle_property_set const* const given_EW_decayer,
                              kinematics_handler* const given_kinematics,
                              input_handler const* const given_shortcuts ) :
    readied_for_new_point( given_shortcuts->get_readier() ),
    initial_scolored( given_initial_scolored ),
    decaying_scolored( NULL ),
    EW_decayer( given_EW_decayer ),
    leptonKinematics( given_kinematics ),
    sQCD_BR( CppSLHA::CppSLHA_global::really_wrong_value ),
    shortcut( given_shortcuts )
  {

    // just an initialization list.

  }

  full_cascade::~full_cascade()
  {

    // does nothing.

  }


  EW_decayer_values*
  full_cascade::get_EW_decay( double const given_primary_cut,
                              double const given_secondary_cut,
                              double const given_jet_cut,
                              lepton_acceptance_value* const given_kinematics )
  /* this checks to see if it already has a EW_decayer_values instance
   * matching the arguments, & if it does, it returns a pointer to it,
   * otherwise, it makes a new one, stores it, & returns a pointer to it.
   */
  {

    EW_decayer_values* return_pointer = NULL;
    // this starts as NULL so that we know if it wasn't found among the
    // existing instances.

    for( std::vector< EW_decayer_values* >::iterator
         set_iterator = EW_decays.begin();
         EW_decays.end() > set_iterator;
         ++set_iterator )
      {

        if( (*set_iterator)->is_requested( given_primary_cut,
                                           given_secondary_cut,
                                           given_jet_cut,
                                           given_kinematics ) )
          {

            return_pointer = *set_iterator;
            set_iterator = EW_decays.end();
            // stop looking.

          }

      }

    if( NULL == return_pointer )
      // if we didn't find an existing match, we make a new instance, store
      // it & return a pointer to it:
      {

        return_pointer = new EW_decayer_values( given_primary_cut,
                                                given_secondary_cut,
                                                given_jet_cut,
                                                given_kinematics,
                                                decaying_scolored,
                                                decaying_EWino,
                                                shortcut );
        EW_decays.push_back( return_pointer );

      }

    return return_pointer;

  }


  squark_directly_to_EWino_cascade::squark_directly_to_EWino_cascade(
            CppSLHA::particle_property_set const* const given_initial_scolored,
                  CppSLHA::particle_property_set const* const given_EW_decayer,
                                    kinematics_handler* const given_kinematics,
                                  std::list< int > const* const given_excluded,
                                 input_handler const* const given_shortcuts ) :
    full_cascade( given_initial_scolored,
                  given_EW_decayer,
                  given_kinematics,
                  given_shortcuts )
  {

    cascade_type = sx;
    decaying_scolored = given_initial_scolored;

    bool decay_to_EWino_not_antiEWino = true;
    // assume that we have either a squark-to-neutralino or a
    // sup-type-to-chargino decay.

    bool decay_to_neutralino;

    if( given_shortcuts->is_in( given_EW_decayer,
                                given_shortcuts->get_charginos() ) )
      // if given_EW_decayer is a chargino...
      {

        decay_to_neutralino = false;

        if( given_shortcuts->is_in( given_initial_scolored->get_PDG_code(),
                                    given_shortcuts->get_sdown_types() ) )
          // if given_initial_scolored is sdown-type...
          {

            decay_to_EWino_not_antiEWino = false;
            // then we need to look for decays with a negative PDG code for
            // given_EW_decayer amongst given_initial_scolored's decays.

          }

      }
    /*else if( CppSLHA::PDG_code::W_plus == given_EW_decayer->get_PDG_code() )
      // if given_EW_decayer is a W boson...
      {

        cascade_type = sw;

        decay_to_neutralino = false;

        if( given_shortcuts->is_in( given_initial_scolored->get_PDG_code(),
                                    given_shortcuts->get_sdown_types() ) )
          // if given_initial_scolored is sdown-type...
          {

            decay_to_EWino_not_antiEWino = false;
            // then we need to look for decays with a negative PDG code for
            // given_EW_decayer amongst given_initial_scolored's decays.

          }

      }*/
    else
      // otherwise, we should have a neutralino...
      {

        decay_to_neutralino = true;

      }

    squark_to_EWino_BR
    = given_shortcuts->get_exclusive_BR_calculator( given_initial_scolored,
                                                    given_EW_decayer,
                                                  decay_to_EWino_not_antiEWino,
                                                    given_excluded );

    lepton_acceptance
    = given_kinematics->get_lepton_acceptance( given_initial_scolored,
                                               given_initial_scolored,
                                               given_EW_decayer,
                                               cascade_type );

  }

  squark_directly_to_EWino_cascade::~squark_directly_to_EWino_cascade()
  {

    // does nothing.

  }


  void
  squark_directly_to_EWino_cascade::update_self()
  // this overrides the base full_cascade::update_self() appropriately.
  {

    sQCD_BR = squark_to_EWino_BR->get_BR();

  }


  gluino_directly_to_EWino_cascade::gluino_directly_to_EWino_cascade(
            CppSLHA::particle_property_set const* const given_initial_scolored,
                  CppSLHA::particle_property_set const* const given_EW_decayer,
                                    kinematics_handler* const given_kinematics,
                                  std::list< int > const* const given_excluded,
                                 input_handler const* const given_shortcuts ) :
    full_cascade( given_initial_scolored,
                  given_EW_decayer,
                  given_kinematics,
                  given_shortcuts )
  {

    cascade_type = gx;
    decaying_scolored = given_initial_scolored;

    gluino_to_EWino_BR
    = given_shortcuts->get_exclusive_BR_calculator( given_initial_scolored,
                                                    given_EW_decayer,
                                                    true,
                                                    given_excluded );

    lepton_acceptance
    = given_kinematics->get_lepton_acceptance( given_initial_scolored,
                                               given_initial_scolored,
                                               given_EW_decayer,
                                               cascade_type );

  }

  gluino_directly_to_EWino_cascade::~gluino_directly_to_EWino_cascade()
  {

    // does nothing.

  }


  void
  gluino_directly_to_EWino_cascade::update_self()
  // this overrides the base full_cascade::update_self() appropriately.
  {

    sQCD_BR = gluino_to_EWino_BR->get_BR();

  }


  fullCascade::fullCascade( colorfulCascadeType const typeOfColorfulCascade ) :
    shortcut( NULL ),
    initialScolored( NULL ),
    beamEnergy( CppSLHA::CppSLHA_global::really_wrong_value ),
    typeOfColorfulCascade( typeOfColorfulCascade ),
    cascadeDefiner(),
    ewinoCascade( NULL ),
    vectorCascade( NULL ),
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
    if( NULL == vectorCascade )
    {
      /* if there is no possibility of extra jets from a vector decay, we can
       * only get a single OSSF-OSDF pair from an electroweakino decay, which
       * cannot produce additional jets.
       */
      if( 0 == numberOfAdditionalJets )
      {
        return ewinoCascade->getOssfMinusOsdf( cuts );
      }
      else
      {
        return 0.0;
      }
    }
    else
    {
      return ( ewinoCascade->getOssfMinusOsdf( cuts )
               * vectorCascade->getAcceptance( true
                                             /* doesn't matter for 0 leptons*/,
                                               cuts,
                                               numberOfAdditionalJets,
                                               0,
                                               0,
                                               0,
                                               0 )
               + vectorCascade->getOssfMinusOsdf( cuts )
                 * ewinoCascade->getAcceptance( true
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
      if( NULL == vectorCascade )
      {
        /* if there is no possibility of extra jets from a vector decay, we can
         * only get a single OSSF-OSDF pair from an electroweakino decay, which
         * cannot produce additional jets.
         */
        if( 1 == numberOfLeptonPairs )
        {
          return ewinoCascade->getOssfMinusOsdf( cuts );
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
          return ( ewinoCascade->getOssfMinusOsdf( cuts )
                   * vectorCascade->getOssfMinusOsdf( cuts ) );
        }
        else if( 1 == numberOfLeptonPairs )
        {
        return ( ewinoCascade->getOssfMinusOsdf( cuts )
                 * ( vectorCascade->getAcceptance( true
                                            /* doesn't matter for 0 leptons */,
                                                 cuts,
                                                 0,
                                                 0,
                                                 0,
                                                 0,
                                                 0 )
                     + vectorCascade->getAcceptance( true
                                            /* doesn't matter for 0 leptons */,
                                                     cuts,
                                                     1,
                                                     0,
                                                     0,
                                                     0,
                                                     0 )
                     + vectorCascade->getAcceptance( true
                                            /* doesn't matter for 0 leptons */,
                                                     cuts,
                                                     2,
                                                     0,
                                                     0,
                                                     0,
                                                     0 ) )
                 + vectorCascade->getOssfMinusOsdf( cuts )
                   * ( ewinoCascade->getAcceptance( true
                                            /* doesn't matter for 0 leptons */,
                                                    cuts,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0 )
                       + ewinoCascade->getAcceptance( true
                                            /* doesn't matter for 0 leptons */,
                                                      cuts,
                                                      1,
                                                      0,
                                                      0,
                                                      0,
                                                      0 )
                       + ewinoCascade->getAcceptance( true
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
                  ewinoAcceptance = ewinoCascade->getAcceptance( cuts,
                                                                 ewinoJets,
                                                        ewinoNegativeElectrons,
                                                        ewinoPositiveElectrons,
                                                            ewinoNegativeMuons,
                                                          ewinoPositiveMuons );
                }
                else
                {
                  ewinoAcceptance = ewinoCascade->getAcceptance( cuts,
                                                                 ewinoJets,
                                                        ewinoPositiveElectrons,
                                                        ewinoNegativeElectrons,
                                                            ewinoPositiveMuons,
                                                          ewinoNegativeMuons );
                }
                if( vectorScoloredIsNotAntiparticle )
                {
                  vectorAcceptance = vectorCascade->getAcceptance( cuts,
                                        ( numberOfAdditionalJets - ewinoJets ),
                        ( numberOfNegativeElectrons - ewinoNegativeElectrons ),
                        ( numberOfPositiveElectrons - ewinoPositiveElectrons ),
                                ( numberOfNegativeMuons - ewinoNegativeMuons ),
                              ( numberOfPositiveMuons - ewinoPositiveMuons ) );
                }
                else
                {
                  vectorAcceptance = vectorCascade->getAcceptance( cuts,
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
    fullCascade( sx )
  {
    // just an initialization list.
  }

  sxFullCascade::~sxFullCascade()
  {
    // does nothing.
  }

  double
  sxFullCascade::getAcceptance( bool const scoloredIsNotAntiparticle,
                                acceptanceCutSet* const cuts,
                                int const numberOfAdditionalJets,
                                int const numberOfNegativeElectrons,
                                int const numberOfPositiveElectrons,
                                int const numberOfNegativeMuons,
                                int const numberOfPositiveMuons )
  /* this calls the appropriate functions on ewinoCascade &, if not NULL,
   * vectorCascade, to build the required acceptance, taking into account
   * whether the charges should be swapped if scoloredIsNotAntiparticle is
   * false.
   */
  {
    if( scoloredIsNotAntiparticle )
      // if we have a particle, we return the acceptance for the particle.
    {
      return ewinoCascade->getAcceptance( cuts,
                                          numberOfAdditionalJets,
                                          numberOfNegativeElectrons,
                                          numberOfPositiveElectrons,
                                          numberOfNegativeMuons,
                                          numberOfPositiveMuons );
    }
    else
      // if we have an antiparticle, we swap the charges.
    {
      return ewinoCascade->getAcceptance( cuts,
                                          numberOfAdditionalJets,
                                          numberOfPositiveElectrons,
                                          numberOfNegativeElectrons,
                                          numberOfPositiveMuons,
                                          numberOfNegativeMuons );
    }
  }



  gxFullCascade::gxFullCascade() :
    fullCascade( gx )
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
    if( ewinoCascade->getElectroweakDecayer()->counts_as_self_conjugate() )
      // if the gluino has only 1 charge version of this decay...
    {
      return ewinoCascade->getAcceptance( cuts,
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
      return ( 0.5 * ( ewinoCascade->getAcceptance( cuts,
                                                    numberOfAdditionalJets,
                                                    numberOfNegativeElectrons,
                                                    numberOfPositiveElectrons,
                                                    numberOfNegativeMuons,
                                                    numberOfPositiveMuons )
                       + ewinoCascade->getAcceptance( cuts,
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
    fullCascade( sjgx )
  {
    // just an initialization list.
  }

  sjgxFullCascade::~sjgxFullCascade()
  {
    // does nothing.
  }


  gjsxFullCascade::gjsxFullCascade() :
    fullCascade( gjsx )
  {
    // just an initialization list.
  }

  gjsxFullCascade::~gjsxFullCascade()
  {
    // does nothing.
  }

  inline void
  gjsxFullCascade::setProperties( sxFullCascade* const sxCascade )
  {
    this->sxCascade = sxCascade;
    this->shortcut = sxCascade->shortcut;
    this->initialScolored = shortcut->get_gluino();
    this->beamEnergy = sxCascade->beamEnergy;
    this->ewinoCascade = sxCascade->ewinoCascade;
    soughtDecayProductList.front()
    = sxCascade->getInitialScolored()->get_PDG_code();
    copyCascadeDefiner( sxCascade->getCascadeDefiner() );
    // reset cascadeDefiner & fill it with a copy of sxCascade's version.
    cascadeSegment = cascadeDefiner.addNewAtEnd();
    cascadeSegment->first = sxCascade->getInitialScolored();
    // sjgx means 1 more decay to be recorded on top of those of the given
    // cascade.
    cascadeSegment->second = 2;
    // sjgx also means that the initial decay is 2-body.
  }

  inline double
  gjsxFullCascade::getBrToEwino( std::list< int > const* excludedSmParticles )
  /* this works out the branching ratio for the decays of colored sparticles
   * down to the 1st electroweakino in the cascade (including any decays to
   * vector bosons + squarks), using the given list of excluded SM particles
   * to exclude unwanted parts of the BR.
   */
  {
    // we return the BR to a neutralino, or, if the electroweakino is a
    // chargino, we return the sum of BRs to both charges of chargino.
    return ( initialScolored->inspect_direct_decay_handler(
                    )->get_branching_ratio_for_subset( &soughtDecayProductList,
                                                       excludedSmParticles )
             * 2.0 * sxCascade->getBrToEwino( excludedSmParticles ) );
  }

  double
  gjsxFullCascade::getAcceptance( bool const scoloredIsNotAntiparticle,
                                  acceptanceCutSet* const cuts,
                                  int const numberOfAdditionalJets,
                                  int const numberOfNegativeElectrons,
                                  int const numberOfPositiveElectrons,
                                  int const numberOfNegativeMuons,
                                  int const numberOfPositiveMuons )
  /* the assumptions made by this code mean that the branching ratios to
   * electroweakinos multiplied by acceptances through squark & antisquark
   * factorize. if the electroweakino is a chargino, it is assumed that only
   * either the squark or the antisquark has a non-zero BR to the chargino,
   * hence the compound BR from the gluino to the positively-charged chargino
   * is given by getBrToEwino, which will be the BR through the appropriate
   * charge of squark.
   */
  {
    return ( 0.5 * ( sxCascade->getAcceptance( true,
                                               cuts,
                                               numberOfAdditionalJets,
                                               numberOfNegativeElectrons,
                                               numberOfPositiveElectrons,
                                               numberOfNegativeMuons,
                                               numberOfPositiveMuons )
                     + sxCascade->getAcceptance( false,
                                                 cuts,
                                                 numberOfAdditionalJets,
                                                 numberOfNegativeElectrons,
                                                 numberOfPositiveElectrons,
                                                 numberOfNegativeMuons,
                                                 numberOfPositiveMuons ) ) );
    // return the weighted acceptances of the squark & the antisquark versions.
  }

  sjgjsxFullCascade::sjgjsxFullCascade() :
    fullCascade( sjgjsx )
  {
    // just an initialization list.
  }

  sjgjsxFullCascade::~sjgjsxFullCascade()
  {
    // does nothing.
  }

  inline void
  sjgjsxFullCascade::setProperties(
                     CppSLHA::particle_property_set const* const initialSquark,
                                    gjsxFullCascade* const gjsxCascade )
  {
    this->gjsxCascade = gjsxCascade;
    this->shortcut = gjsxCascade->shortcut;
    this->initialScolored = initialSquark;
    this->beamEnergy = gjsxCascade->beamEnergy;
    this->ewinoCascade = gjsxCascade->ewinoCascade;
    soughtDecayProductList.front() = CppSLHA::PDG_code::gluino;
    copyCascadeDefiner( gjsxCascade->getCascadeDefiner() );
    // reset cascadeDefiner & fill it with a copy of gjsxCascade's version.
    cascadeSegment = cascadeDefiner.addNewAtEnd();
    cascadeSegment->first = gjsxCascade->getInitialScolored();
    // sjgjsx means 1 more decay to be recorded on top of those of the given
    // cascade.
    cascadeSegment->second = 2;
    // sjgjsx also means that the initial decay is 2-body.
  }

  inline double
  sjgjsxFullCascade::getBrToEwino(
                                  std::list< int > const* excludedSmParticles )
  /* this works out the branching ratio for the decays of colored sparticles
   * down to the 1st electroweakino in the cascade (including any decays to
   * vector bosons + squarks), using the given list of excluded SM particles
   * to exclude unwanted parts of the BR.
   */
  {
    return ( initialScolored->inspect_direct_decay_handler(
                    )->get_branching_ratio_for_subset( &soughtDecayProductList,
                                                       excludedSmParticles )
             * gjsxCascade->getBrToEwino( excludedSmParticles ) );
  }

  inline double
  sjgjsxFullCascade::getAcceptance( bool const scoloredIsNotAntiparticle,
                                    acceptanceCutSet* const cuts,
                                    int const numberOfAdditionalJets,
                                    int const numberOfNegativeElectrons,
                                    int const numberOfPositiveElectrons,
                                    int const numberOfNegativeMuons,
                                    int const numberOfPositiveMuons )
  {
    return gjsxCascade->getAcceptance( scoloredIsNotAntiparticle,
                                       cuts,
                                       numberOfAdditionalJets,
                                       numberOfNegativeElectrons,
                                       numberOfPositiveElectrons,
                                       numberOfNegativeMuons,
                                       numberOfPositiveMuons );
    // gjsxCascade takes care of counting over different charginos if needed.
  }

  svsxFullCascade::svsxFullCascade() :
    fullCascade( svsx )
  {
    soughtDecayProductList.push_back(
                                 CppSLHA::CppSLHA_global::really_wrong_value );
    // a cascade beginning with sv needs to look for 2 specific decay products.
  }

  svsxFullCascade::~svsxFullCascade()
  {
    // does nothing.
  }

  void
  svsxFullCascade::setProperties(
                     CppSLHA::particle_property_set const* const initialSquark,
                                  electroweakCascade* const vectorCascade,
                                  sxFullCascade* const sxCascade )
  {
    this->sxCascade = sxCascade;
    this->shortcut = sxCascade->shortcut;
    this->initialScolored = initialSquark;
    this->beamEnergy = sxCascade->beamEnergy;
    this->ewinoCascade = sxCascade->ewinoCascade;
    this->vectorCascade = vectorCascade;
    soughtDecayProductList.front()
    = sxCascade->getInitialScolored()->get_PDG_code();
    if( ( CppSLHA::PDG_code::W_plus
          == vectorCascade->getElectroweakDecayer()->get_PDG_code() )
        &&
        ( shortcut->is_in( initialSquark->get_PDG_code(),
                           shortcut->get_sdown_types() ) ) )
      // if we have to worry about which sign of PDG code to use...
    {
      soughtDecayProductList.back() = -CppSLHA::PDG_code::W_plus;
    }
    else
    {
    soughtDecayProductList.back()
    = vectorCascade->getElectroweakDecayer()->get_PDG_code();
    }
    copyCascadeDefiner( sxCascade->getCascadeDefiner() );
    // reset cascadeDefiner & fill it with a copy of gjsxCascade's version.
    cascadeSegment = cascadeDefiner.addNewAtEnd();
    cascadeSegment->first = sxCascade->getInitialScolored();
    // svsx means 1 more decay to be recorded on top of those of the given
    // cascade.
    cascadeSegment->second = 2;
    // svsx also means that the initial decay is 2-body.
  }

  double
  svsxFullCascade::getBrToEwino( std::list< int > const* excludedSmParticles )
  /* this works out the branching ratio for the decays of colored sparticles
   * down to the 1st electroweakino in the cascade (including any decays to
   * vector bosons + squarks), using the given list of excluded SM particles
   * to exclude unwanted parts of the BR.
   */
  {
    return ( initialScolored->inspect_direct_decay_handler(
                    )->get_branching_ratio_for_subset( &soughtDecayProductList,
                                                       excludedSmParticles )
             * sxCascade->getBrToEwino( excludedSmParticles ) );
  }

  inline double
  svsxFullCascade::getAcceptance( bool const scoloredIsNotAntiparticle,
                                  acceptanceCutSet* const cuts,
                                  int const numberOfAdditionalJets,
                                  int const numberOfNegativeElectrons,
                                  int const numberOfPositiveElectrons,
                                  int const numberOfNegativeMuons,
                                  int const numberOfPositiveMuons )
  /* this calls the appropriate functions on ewinoCascade &, if not NULL,
   * vectorCascade, to build the required acceptance, taking into account
   * whether the charges should be swapped if scoloredIsNotAntiparticle is
   * false.
   */
  {
    return getCombinedAcceptance( scoloredIsNotAntiparticle,
                                  scoloredIsNotAntiparticle,
                                  cuts,
                                  numberOfAdditionalJets,
                                  numberOfNegativeElectrons,
                                  numberOfPositiveElectrons,
                                  numberOfNegativeMuons,
                                  numberOfPositiveMuons );
    // in the used conventions, an incoming squark decays to a vector boson
    // plus a squark, hence the same bool is used for both.
  }


  gvsxFullCascade::gvsxFullCascade() :
    fullCascade( svsx )
  {
    soughtDecayProductList.push_back(
                                 CppSLHA::CppSLHA_global::really_wrong_value );
    // a cascade beginning with gv needs to look for 2 specific decay products.
  }

  gvsxFullCascade::~gvsxFullCascade()
  {
    // does nothing.
  }

  void
  gvsxFullCascade::setProperties( electroweakCascade* const vectorCascade,
                                  sxFullCascade* const sxCascade )
  {
    this->sxCascade = sxCascade;
    this->shortcut = sxCascade->shortcut;
    this->initialScolored = shortcut->get_gluino();
    this->beamEnergy = sxCascade->beamEnergy;
    this->ewinoCascade = sxCascade->ewinoCascade;
    this->vectorCascade = vectorCascade;
    soughtDecayProductList.front()
    = sxCascade->getInitialScolored()->get_PDG_code();
    if( ( CppSLHA::PDG_code::W_plus
          == vectorCascade->getElectroweakDecayer()->get_PDG_code() )
        &&
        ( shortcut->is_in( soughtDecayProductList.front(),
                           shortcut->get_sup_types() ) ) )
      // if we have to worry about which sign of PDG code to use...
    {
      soughtDecayProductList.back() = -CppSLHA::PDG_code::W_plus;
    }
    else
    {
    soughtDecayProductList.back()
    = vectorCascade->getElectroweakDecayer()->get_PDG_code();
    }
    copyCascadeDefiner( sxCascade->getCascadeDefiner() );
    // reset cascadeDefiner & fill it with a copy of gjsxCascade's version.
    cascadeSegment = cascadeDefiner.addNewAtEnd();
    cascadeSegment->first = sxCascade->getInitialScolored();
    // gvsx means 1 more decay to be recorded on top of those of the given
    // cascade.
    cascadeSegment->second = 3;
    // gvsx also means that the initial decay is 3-body.
  }

  double
  gvsxFullCascade::getBrToEwino( std::list< int > const* excludedSmParticles )
  /* this works out the branching ratio for the decays of colored sparticles
   * down to the 1st electroweakino in the cascade (including any decays to
   * vector bosons + squarks), using the given list of excluded SM particles
   * to exclude unwanted parts of the BR.
   */
  {
    return ( initialScolored->inspect_direct_decay_handler(
                    )->get_branching_ratio_for_subset( &soughtDecayProductList,
                                                       excludedSmParticles )
             * 2.0 * sxCascade->getBrToEwino( excludedSmParticles ) );
    // we sum over squark and antisquark.
  }

  inline double
  gvsxFullCascade::getAcceptance( bool const scoloredIsNotAntiparticle,
                                  acceptanceCutSet* const cuts,
                                  int const numberOfAdditionalJets,
                                  int const numberOfNegativeElectrons,
                                  int const numberOfPositiveElectrons,
                                  int const numberOfNegativeMuons,
                                  int const numberOfPositiveMuons )
  /* this calls the appropriate functions on ewinoCascade &, if not NULL,
   * vectorCascade, to build the required acceptance, taking into account
   * whether the charges should be swapped if scoloredIsNotAntiparticle is
   * false.
   */
  {
    return ( 0.5 * ( getCombinedAcceptance( true,
                                            true,
                                            cuts,
                                            numberOfAdditionalJets,
                                            numberOfNegativeElectrons,
                                            numberOfPositiveElectrons,
                                            numberOfNegativeMuons,
                                            numberOfPositiveMuons )
                     + getCombinedAcceptance( false,
                                              false,
                                              cuts,
                                              numberOfAdditionalJets,
                                              numberOfNegativeElectrons,
                                              numberOfPositiveElectrons,
                                              numberOfNegativeMuons,
                                              numberOfPositiveMuons ) ) );
    // getBrToEwino gives the sum of the decays via squark & antisquark, so the
    // acceptances factorize this way.
  }


  gjsvsxFullCascade::gjsvsxFullCascade() :
    fullCascade( gjsvsx )
  {
    // just an initialization list.
  }

  gjsvsxFullCascade::~gjsvsxFullCascade()
  {
    // does nothing.
  }

  void
  gjsvsxFullCascade::setProperties( svsxFullCascade* const svsxCascade )
  {
    this->svsxCascade = svsxCascade;
    this->shortcut = svsxCascade->shortcut;
    this->initialScolored = shortcut->get_gluino();
    this->beamEnergy = svsxCascade->beamEnergy;
    this->ewinoCascade = svsxCascade->ewinoCascade;
    soughtDecayProductList.front()
    = svsxCascade->getInitialScolored()->get_PDG_code();
    copyCascadeDefiner( svsxCascade->getCascadeDefiner() );
    // reset cascadeDefiner & fill it with a copy of sxCascade's version.
    cascadeSegment = cascadeDefiner.addNewAtEnd();
    cascadeSegment->first = svsxCascade->getInitialScolored();
    // gjsvsx means 1 more decay to be recorded on top of those of the given
    // cascade.
    cascadeSegment->second = 2;
    // gjsvsx also means that the initial decay is 2-body.
  }

  double
  gjsvsxFullCascade::getBrToEwino(
                                  std::list< int > const* excludedSmParticles )
  /* this works out the branching ratio for the decays of colored sparticles
   * down to the 1st electroweakino in the cascade (including any decays to
   * vector bosons + squarks), using the given list of excluded SM particles
   * to exclude unwanted parts of the BR.
   */
  {
    return ( initialScolored->inspect_direct_decay_handler(
                    )->get_branching_ratio_for_subset( &soughtDecayProductList,
                                                       excludedSmParticles )
             * 2.0 * svsxCascade->getBrToEwino( excludedSmParticles ) );
    // we sum over squark and antisquark.
  }

  virtual double
  gjsvsxFullCascade::getAcceptance( bool const scoloredIsNotAntiparticle,
                                    acceptanceCutSet* const cuts,
                                    int const numberOfAdditionalJets,
                                    int const numberOfNegativeElectrons,
                                    int const numberOfPositiveElectrons,
                                    int const numberOfNegativeMuons,
                                    int const numberOfPositiveMuons )
  /* the assumptions made by this code mean that the branching ratios to
   * electroweakinos multiplied by acceptances through squark & antisquark
   * factorize. if the electroweakino is a chargino, it is assumed that only
   * either the squark or the antisquark has a non-zero BR to the chargino,
   * hence the compound BR from the gluino to the positively-charged chargino
   * is given by getBrToEwino, which will be the BR through the appropriate
   * charge of squark.
   */
  {
    return ( 0.5 * ( svsxCascade->getAcceptance( true,
                                                 cuts,
                                                 numberOfAdditionalJets,
                                                 numberOfNegativeElectrons,
                                                 numberOfPositiveElectrons,
                                                 numberOfNegativeMuons,
                                                 numberOfPositiveMuons )
                     + svsxCascade->getAcceptance( false,
                                                   cuts,
                                                   numberOfAdditionalJets,
                                                   numberOfNegativeElectrons,
                                                   numberOfPositiveElectrons,
                                                   numberOfNegativeMuons,
                                                   numberOfPositiveMuons ) ) );
    // return the weighted acceptances of the squark & the antisquark versions.
  }



}  // end of LHC_FASER namespace.
