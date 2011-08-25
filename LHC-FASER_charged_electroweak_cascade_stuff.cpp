/*
 * LHC-FASER_derived_electroweak_cascade_stuff.cpp
 *
 *  Created on: May 9, 2011
 *      Author: bol
 */

#include "LHC-FASER_derived_electroweak_cascade_stuff.hpp"

namespace LHC_FASER
{
  // this is here just for reference! it will be deleted later, when I am sure
  // that I have covered all the decays it covers.
  EW_decayer_values::EW_decayer_values( bool const given_should_sum_charges,
                                        double const given_primary_cut,
                                        double const given_secondary_cut,
                                        double const given_jet_cut,
                               lepton_acceptance_value* const given_kinematics,
           particlePointer const given_decaying_scolored,
                                 bool const given_scolored_is_not_antiparticle,
                  particlePointer const given_EW_decayer,
                                    bool const given_EWino_is_not_antiparticle,
                                 inputHandler const* const given_shortcuts ) :
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
    zeroJetsOnePositiveElectron(
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

    if( given_decaying_EWino == given_shortcuts->getNeutralinoOne() )
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
    else if( given_shortcuts->isIn( given_decaying_EWino->get_PDG_code(),
                                given_shortcuts->getUnstableNeutralinos() ) )
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
                                            given_shortcuts->getSelectronL(),
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
                                            given_shortcuts->getSelectronR(),
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
                                                given_shortcuts->getSmuonL(),
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
                                                given_shortcuts->getSmuonR(),
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
                                               given_shortcuts->getStauOne(),
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
                                               given_shortcuts->getStauOne(),
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
                                               given_shortcuts->getStauOne(),
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
                                               given_shortcuts->getStauTwo(),
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
                                               given_shortcuts->getStauTwo(),
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
                                               given_shortcuts->getStauTwo(),
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
                                                      given_shortcuts->getZ(),
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
                                                      given_shortcuts->getZ(),
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
                                                      given_shortcuts->getZ(),
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
                                                      given_shortcuts->getZ(),
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
                              given_shortcuts->getLightNeutralEwsbScalar(),
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
                              given_shortcuts->getLightNeutralEwsbScalar(),
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
                              given_shortcuts->getLightNeutralEwsbScalar(),
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
                              given_shortcuts->getLightNeutralEwsbScalar(),
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
                              given_shortcuts->getHeavyNeutralEwsbScalar(),
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
                              given_shortcuts->getHeavyNeutralEwsbScalar(),
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
                              given_shortcuts->getHeavyNeutralEwsbScalar(),
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
                              given_shortcuts->getHeavyNeutralEwsbScalar(),
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
                              given_shortcuts->getNeutralEwsbPseudoscalar(),
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
                              given_shortcuts->getNeutralEwsbPseudoscalar(),
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
                              given_shortcuts->getNeutralEwsbPseudoscalar(),
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
                              given_shortcuts->getNeutralEwsbPseudoscalar(),
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
    else if( given_shortcuts->isIn( given_decaying_EWino->get_PDG_code(),
                                     given_shortcuts->getCharginos() ) )
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
                                            given_shortcuts->getSelectronL(),
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
                                            given_shortcuts->getSelectronR(),
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
                                              given_shortcuts->getSmuonL(),
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
                                              given_shortcuts->getSmuonR(),
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
                                               given_shortcuts->getStauOne(),
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
                                               given_shortcuts->getStauOne(),
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
                                               given_shortcuts->getStauTwo(),
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
                                               given_shortcuts->getStauTwo(),
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
                                   given_shortcuts->getElectronSneutrinoL(),
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
                                   given_shortcuts->getElectronSneutrinoR(),
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
                                       given_shortcuts->getMuonSneutrinoL(),
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
                                       given_shortcuts->getMuonSneutrinoR(),
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
                                        given_shortcuts->getTauSneutrinoL(),
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
                                        given_shortcuts->getTauSneutrinoL(),
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
                                        given_shortcuts->getTauSneutrinoR(),
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
                                        given_shortcuts->getTauSneutrinoR(),
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
                                                 given_shortcuts->getWPlus(),
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
                                                 given_shortcuts->getWPlus(),
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
                                                 given_shortcuts->getWPlus(),
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
                                    given_shortcuts->getChargedEwsbScalar(),
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
                                    given_shortcuts->getChargedEwsbScalar(),
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
                                    given_shortcuts->getChargedEwsbScalar(),
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




  charginoToSemuOrEmuSnuCascade::charginoToSemuOrEmuSnuCascade(
                                leptonAcceptanceParameterSet* const kinematics,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                CppSLHA::particle_property_set const* const electroweakDecayer,
               CppSLHA::particle_property_set const* const intermediateDecayer,
                                           inputHandler const* const shortcut,
                                                bool const sneutrinoVersion ) :
    electroweakCascade( kinematics,
                        effectiveSquarkMass,
                        coloredDecayer,
                        electroweakDecayer,
                        intermediateDecayer,
                        false,
                        shortcut ),
    sneutrinoVersion( sneutrinoVersion )
  {
    if( ( CppSLHA::PDG_code::smuon_L == intermediateDecayer->get_PDG_code() )
        ||
        ( CppSLHA::PDG_code::smuon_R == intermediateDecayer->get_PDG_code() )
        ||
        ( CppSLHA::PDG_code::muon_sneutrino_L
          == intermediateDecayer->get_PDG_code() )
        ||
        ( CppSLHA::PDG_code::muon_sneutrino_L
          == intermediateDecayer->get_PDG_code() ) )
    {
      muonsNotElectrons = true;
    }
    else
    {
      muonsNotElectrons = false;
    }
    // we set up the BR of the (only) channel:
    firstBr
    = shortcut->getExclusiveBrCalculator( electroweakDecayer,
                                             intermediateDecayer,
                                             true,
                                             shortcut->getEmptyList() );
    secondBr
    = shortcut->getExclusiveBrCalculator( intermediateDecayer,
                                             shortcut->getNeutralinoOne(),
                                             true,
                                             shortcut->getEmptyList() );
    if( sneutrinoVersion )
    {
      sameDistribution = new same_chirality_near_muon( shortcut->getReadier(),
                                                       shortcut->getCppSlha(),
                                                       coloredDecayer,
                                                       effectiveSquarkMass,
                                                       electroweakDecayer,
                                                       intermediateDecayer,
                                              shortcut->getNeutralinoOne() );
      oppositeDistribution
      = new opposite_chirality_near_muon( shortcut->getReadier(),
                                          shortcut->getCppSlha(),
                                          coloredDecayer,
                                          effectiveSquarkMass,
                                          electroweakDecayer,
                                          intermediateDecayer,
                                          shortcut->getNeutralinoOne() );
    }
    else
    {
      sameDistribution = new same_chirality_far_muon( shortcut->getReadier(),
                                                      shortcut->getCppSlha(),
                                                      coloredDecayer,
                                                      effectiveSquarkMass,
                                                      electroweakDecayer,
                                                      intermediateDecayer,
                                              shortcut->getNeutralinoOne() );
      oppositeDistribution
      = new opposite_chirality_far_muon( shortcut->getReadier(),
                                         shortcut->getCppSlha(),
                                         coloredDecayer,
                                         effectiveSquarkMass,
                                         electroweakDecayer,
                                         intermediateDecayer,
                                         shortcut->getNeutralinoOne() );
    }
    activeDistributions.push_back( sameDistribution );
    activeDistributions.push_back( oppositeDistribution );
  }

  charginoToSemuOrEmuSnuCascade::~charginoToSemuOrEmuSnuCascade()
  {
    // does nothing.
  }


  bool
  charginoToSemuOrEmuSnuCascade::validSignal( int const numberOfJets,
                                           int const numberOfNegativeElectrons,
                                           int const numberOfPositiveElectrons,
                                              int const numberOfNegativeMuons,
                                              int const numberOfPositiveMuons )
  /* this returns true if numberOfJets == 0, & either
   * numberOfPositiveElectrons == 1 OR numberOfPositiveMuons == 1, the rest 0,
   * or if they are all 0.
   */
  {
    if( ( 0 == numberOfJets )
        &&
        ( 0 == numberOfNegativeElectrons )
        &&
        ( 0 == numberOfNegativeMuons )
        &&
        ( ( ( 0 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 1 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 0 == numberOfPositiveElectrons )
            &&
            ( 1 == numberOfPositiveMuons ) ) ) )
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  void
  charginoToSemuOrEmuSnuCascade::calculateAcceptance(
                                                  acceptanceCutSet* const cuts,
                                    acceptanceValues* const currentAcceptance )
  // this returns the appropriate acceptances multiplied by branching ratios
  // from the electroweakino through the slepton to the LSP.
  {
    cascadeBr = ( first_BR->getBr() * second_BR->getBr() );
    if( lhcFaserGlobal::negligibleBr < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      jetLeftHandedness
      = shortcut->quarkOrLeptonLeftHandedness(
                                                coloredDecayer->get_PDG_code(),
                                          electroweakDecayer->get_PDG_code() );
      if( shortcut->isIn( coloredDecayer->get_PDG_code(),
                           shortcut->getSdownTypes() ) )
        /* if we've calculated the wrong handedness, because isIn assumes the
         * correct charge of chargino for a positive-PDG-coded quark, which is
         * the wrong case for down-type quarks since we are assuming that the
         * chargino is positive...
         */
      {
        jetLeftHandedness = ( 1.0 - jetLeftHandedness );
      }
      nearLeptonLeftHandedness
      = shortcut->quarkOrLeptonLeftHandedness(
                                           intermediateDecayer->get_PDG_code(),
                                          electroweakDecayer->get_PDG_code() );
      if( sneutrinoVersion )
        // positive charginos decaying to sneutrinos produce antileptons, which
        // have the opposite helicity to their charge-conjugates.
      {
        nearLeptonLeftHandedness = ( 1.0 - nearLeptonLeftHandedness );
      }

      nearLeptonSameHandednessTimesBr
      = ( cascadeBr * ( ( jetLeftHandedness * nearLeptonLeftHandedness )
                        + ( ( 1.0 - jetLeftHandedness )
                            * ( 1.0 - nearLeptonLeftHandedness ) ) ) );

      samePass = integrateAcceptance( sameDistribution,
                                      cuts->getPrimaryLeptonCut() );
      sameFail = ( 1.0 - integrateAcceptance( sameDistribution,
                                             cuts->getSecondaryLeptonCut() ) );
      oppositePass = integrateAcceptance( oppositeDistribution,
                                          cuts->getPrimaryLeptonCut() );
      oppositeFail = ( 1.0 - integrateAcceptance( oppositeDistribution,
                                             cuts->getSecondaryLeptonCut() ) );

      currentAcceptance->setZeroJetsZeroLeptons(
                                              ( nearLeptonSameHandednessTimesBr
                                                * sameFail
                                                   + ( cascadeBr
                                            - nearLeptonSameHandednessTimesBr )
                                                     * oppositeFail ) );
      if( muonsNotElectrons )
      {
        currentAcceptance->setZeroJetsOnePositiveMuon(
                                              ( nearLeptonSameHandednessTimesBr
                                                * samePass
                                                + ( cascadeBr
                                            - nearLeptonSameHandednessTimesBr )
                                                  * oppositePass ) );
      }
      else
      {
        currentAcceptance->setZeroJetsOnePositiveElectron(
                                              ( nearLeptonSameHandednessTimesBr
                                                * samePass
                                                + ( cascadeBr
                                            - nearLeptonSameHandednessTimesBr )
                                                  * oppositePass ) );
      }
    }  // end of if cascade branching ratio was large enough to care about.
    // debugging:
    /**std::cout
    << std::endl
    << "cascadeBr = " << cascadeBr;
    std::cout << std::endl;**/
  }



  charginoToStauOrTauSnuCascade::charginoToStauOrTauSnuCascade(
                                leptonAcceptanceParameterSet* const kinematics,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                CppSLHA::particle_property_set const* const electroweakDecayer,
               CppSLHA::particle_property_set const* const intermediateDecayer,
                                           inputHandler const* const shortcut,
                                                bool const sneutrinoVersion ) :
    electroweakCascade( kinematics,
                        effectiveSquarkMass,
                        coloredDecayer,
                        electroweakDecayer,
                        intermediateDecayer,
                        false,
                        shortcut ),
    sneutrinoVersion( sneutrinoVersion )
  {
    // we set up the BR of the (only) channel:
    firstBr
    = shortcut->getExclusiveBrCalculator( electroweakDecayer,
                                             intermediateDecayer,
                                             true,
                                             shortcut->getEmptyList() );
    secondBr
    = shortcut->getExclusiveBrCalculator( intermediateDecayer,
                                             shortcut->getNeutralinoOne(),
                                             true,
                                             shortcut->getEmptyList() );
    if( sneutrinoVersion )
    {
      sameTauDistribution
      = new same_chirality_near_muon( shortcut->getReadier(),
                                      shortcut->getCppSlha(),
                                      coloredDecayer,
                                      effectiveSquarkMass,
                                      electroweakDecayer,
                                      intermediateDecayer,
                                      shortcut->getNeutralinoOne() );
      oppositeTauDistribution
      = new opposite_chirality_near_muon( shortcut->getReadier(),
                                          shortcut->getCppSlha(),
                                          coloredDecayer,
                                          effectiveSquarkMass,
                                          electroweakDecayer,
                                          intermediateDecayer,
                                          shortcut->getNeutralinoOne() );
    }
    else
    {
      sameTauDistribution
      = new same_chirality_far_muon( shortcut->getReadier(),
                                     shortcut->getCppSlha(),
                                     coloredDecayer,
                                     effectiveSquarkMass,
                                     electroweakDecayer,
                                     intermediateDecayer,
                                     shortcut->getNeutralinoOne() );
      oppositeTauDistribution
      = new opposite_chirality_far_muon( shortcut->getReadier(),
                                         shortcut->getCppSlha(),
                                         coloredDecayer,
                                         effectiveSquarkMass,
                                         electroweakDecayer,
                                         intermediateDecayer,
                                         shortcut->getNeutralinoOne() );
    }
    activeDistributions.push_back( sameTauDistribution );
    activeDistributions.push_back( oppositeTauDistribution );

    sameHardMuonDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     sameTauDistribution,
                                     shortcut->getHardMuonFromTau() );
    sameSoftMuonDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     sameTauDistribution,
                                     shortcut->getSoftMuonFromTau() );
    sameHardPionDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     sameTauDistribution,
                                     shortcut->getHardPionFromTau() );
    sameSoftPionDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     sameTauDistribution,
                                     shortcut->getSoftPionFromTau() );
    oppositeHardMuonDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     oppositeTauDistribution,
                                     shortcut->getHardMuonFromTau() );
    oppositeSoftMuonDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     oppositeTauDistribution,
                                     shortcut->getSoftMuonFromTau() );
    oppositeHardPionDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     oppositeTauDistribution,
                                     shortcut->getHardPionFromTau() );
    oppositeSoftPionDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     oppositeTauDistribution,
                                     shortcut->getSoftPionFromTau() );
    activeDistributions.push_back( sameHardMuonDistribution );
    activeDistributions.push_back( sameSoftMuonDistribution );
    activeDistributions.push_back( sameHardPionDistribution );
    activeDistributions.push_back( sameSoftPionDistribution );
    activeDistributions.push_back( oppositeHardMuonDistribution );
    activeDistributions.push_back( oppositeSoftMuonDistribution );
    activeDistributions.push_back( oppositeHardPionDistribution );
    activeDistributions.push_back( oppositeSoftPionDistribution );
  }

  charginoToStauOrTauSnuCascade::~charginoToStauOrTauSnuCascade()
  {
    // does nothing.
  }


  bool
  charginoToStauOrTauSnuCascade::validSignal( int const numberOfJets,
                                           int const numberOfNegativeElectrons,
                                           int const numberOfPositiveElectrons,
                                              int const numberOfNegativeMuons,
                                              int const numberOfPositiveMuons )
  /* this returns true if all the arguments are 0, or all are 0 but for
   * numberOfJets OR numberOfPositiveElectrons OR numberOfPositiveMuons, which
   * is 1.
   */
  {
    if( ( 0 == numberOfNegativeElectrons )
        &&
        ( 0 == numberOfNegativeMuons )
        &&
        ( ( ( 0 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 1 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 0 == numberOfJets )
            &&
            ( 1 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 0 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 1 == numberOfPositiveMuons ) ) ) )
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  void
  charginoToStauOrTauSnuCascade::calculateAcceptance(
                                                  acceptanceCutSet* const cuts,
                                    acceptanceValues* const currentAcceptance )
  // this returns the appropriate acceptances multiplied by branching ratios
  // from the electroweakino through the slepton to the LSP.
  {
    cascadeBr = ( first_BR->getBr() * second_BR->getBr() );
    if( lhcFaserGlobal::negligibleBr < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      jetLeftHandedness
      = shortcut->quarkOrLeptonLeftHandedness(
                                                coloredDecayer->get_PDG_code(),
                                          electroweakDecayer->get_PDG_code() );
      if( shortcut->isIn( coloredDecayer->get_PDG_code(),
                           shortcut->getSdownTypes() ) )
        /* if we've calculated the wrong handedness, because isIn assumes the
         * correct charge of chargino for a positive-PDG-coded quark, which is
         * the wrong case for down-type quarks since we are assuming that the
         * chargino is positive...
         */
      {
        jetLeftHandedness = ( 1.0 - jetLeftHandedness );
      }
      nearLeptonLeftHandedness
      = shortcut->quarkOrLeptonLeftHandedness(
                                           intermediateDecayer->get_PDG_code(),
                                          electroweakDecayer->get_PDG_code() );
      if( sneutrinoVersion )
        // positive charginos decaying to sneutrinos produce antileptons, which
        // have the opposite helicity to their charge-conjugates.
      {
        antitauRightHandedness = nearLeptonLeftHandedness;
        nearLeptonLeftHandedness = ( 1.0 - antitauLeftHandedness );
      }
      else
      {
        antitauRightHandedness = shortcut->quarkOrLeptonLeftHandedness(
                                            mediating_particle->get_PDG_code(),
                              shortcut->getNeutralinoOne()->get_PDG_code() );
      }

      if( sneutrinoVersion )
      {
        muonPass
        = ( cascadeBr
            * ( jetLeftHandedness
                * ( nearLeptonLeftHandedness
                    * integrateAcceptance( sameSoftMuonDistribution,
                                           cuts->getPrimaryLeptonCut() )
                    + ( 1.0 - nearLeptonLeftHandedness )
                      * integrateAcceptance( oppositeHardMuonDistribution,
                                             cuts->getPrimaryLeptonCut() ) )
                + ( 1.0 - jetLeftHandedness )
                  * ( nearLeptonLeftHandedness
                      * integrateAcceptance( oppositeSoftMuonDistribution,
                                             cuts->getPrimaryLeptonCut() )
                      + ( 1.0 - nearLeptonLeftHandedness )
                        * integrateAcceptance( sameHardMuonDistribution,
                                           cuts->getPrimaryLeptonCut() ) ) ) );
        muonFail
        = ( cascadeBr
            * ( 1.0
                - ( jetLeftHandedness
                    * ( nearLeptonLeftHandedness
                        * integrateAcceptance( sameSoftMuonDistribution,
                                               cuts->getSecondaryLeptonCut() )
                        + ( 1.0 - nearLeptonLeftHandedness )
                          * integrateAcceptance( oppositeHardMuonDistribution,
                                              cuts->getSecondaryLeptonCut() ) )
                    + ( 1.0 - jetLeftHandedness )
                      * ( nearLeptonLeftHandedness
                          * integrateAcceptance( oppositeSoftMuonDistribution,
                                                cuts->getSecondaryLeptonCut() )
                          + ( 1.0 - nearLeptonLeftHandedness )
                            * integrateAcceptance( sameHardMuonDistribution,
                                       cuts->getSecondaryLeptonCut() ) ) ) ) );
        pionPass
        = ( cascadeBr
            * ( jetLeftHandedness
                * ( nearLeptonLeftHandedness
                    * integrateAcceptance( sameHardPionDistribution,
                                           cuts->getJetCut() )
                    + ( 1.0 - nearLeptonLeftHandedness )
                      * integrateAcceptance( oppositeSoftPionDistribution,
                                             cuts->getJetCut() ) )
                + ( 1.0 - jetLeftHandedness )
                  * ( nearLeptonLeftHandedness
                      * integrateAcceptance( oppositeHardPionDistribution,
                                             cuts->getJetCut() )
                      + ( 1.0 - nearLeptonLeftHandedness )
                        * integrateAcceptance( sameSoftPionDistribution,
                                           cuts->getJetCut() ) ) ) );
      }
      else
      {
        muonPass
        = ( cascadeBr
            * ( jetLeftHandedness
                * ( antitauRightHandedness
                    * integrateAcceptance( sameHardMuonDistribution,
                                           cuts->getPrimaryLeptonCut() )
                    + ( 1.0 - antitauRightHandedness )
                      * integrateAcceptance( sameSoftMuonDistribution,
                                             cuts->getPrimaryLeptonCut() ) )
                + ( 1.0 - jetLeftHandedness )
                  * ( antitauRightHandedness
                      * integrateAcceptance( oppositeHardMuonDistribution,
                                             cuts->getPrimaryLeptonCut() )
                      + ( 1.0 - antitauRightHandedness )
                        * integrateAcceptance( oppositeSoftMuonDistribution,
                                           cuts->getPrimaryLeptonCut() ) ) ) );
        muonFail
        = ( cascadeBr
            * ( 1.0
                - ( jetLeftHandedness
                    * ( antitauRightHandedness
                        * integrateAcceptance( sameHardMuonDistribution,
                                               cuts->getSecondaryLeptonCut() )
                        + ( 1.0 - antitauRightHandedness )
                          * integrateAcceptance( sameSoftMuonDistribution,
                                              cuts->getSecondaryLeptonCut() ) )
                    + ( 1.0 - jetLeftHandedness )
                      * ( antitauRightHandedness
                          * integrateAcceptance( oppositeHardMuonDistribution,
                                                cuts->getSecondaryLeptonCut() )
                          + ( 1.0 - antitauRightHandedness )
                           * integrateAcceptance( oppositeSoftMuonDistribution,
                                       cuts->getSecondaryLeptonCut() ) ) ) ) );
        pionPass
        = ( cascadeBr
            * ( jetLeftHandedness
                * ( antitauRightHandedness
                    * integrateAcceptance( sameSoftPionDistribution,
                                           cuts->getJetCut() )
                    + ( 1.0 - antitauRightHandedness )
                      * integrateAcceptance( sameHardPionDistribution,
                                             cuts->getJetCut() ) )
                + ( 1.0 - jetLeftHandedness )
                  * ( antitauRightHandedness
                      * integrateAcceptance( oppositeSoftPionDistribution,
                                             cuts->getJetCut() )
                      + ( 1.0 - antitauRightHandedness )
                        * integrateAcceptance( oppositeHardPionDistribution,
                                           cuts->getJetCut() ) ) ) );
      }


      currentAcceptance->setOneJetZeroLeptons( pionPass
                       * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR );
      currentAcceptance->setZeroJetsOnePositiveElectron( muonPass
                    * CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR );
      currentAcceptance->setZeroJetsOnePositiveMuon( muonPass
                        * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR );

      currentAcceptance->setZeroJetsZeroLeptons( ( 1.0 - pionPass )
                          * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
                                                 + muonFail
                     * ( CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR
                      + CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR ) );
    }  // end of if cascade branching ratio was large enough to care about.
    // debugging:
    /**std::cout
    << std::endl
    << "cascadeBr = " << cascadeBr;
    std::cout << std::endl;**/
  }



  charginoToWCascade::charginoToWCascade(
                                leptonAcceptanceParameterSet* const kinematics,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                CppSLHA::particle_property_set const* const electroweakDecayer,
                                        inputHandler const* const shortcut ) :
    electroweakCascade( kinematics,
                        effectiveSquarkMass,
                        coloredDecayer,
                        electroweakDecayer,
                        shortcut->getWPlus(),
                        false,
                        shortcut )
  {
    if( shortcut->isIn( coloredDecayer->get_PDG_code(),
                         shortcut->getSdownTypes() ) )
    {
      downTypeQuark = true;
    }
    else
    {
      downTypeQuark = false;
    }

    // here firstBr & secondBr don't work as well as just getting the branching
    // ratios from the CppSLHA::particle_property_set pointers.

    directJetDistribution = new W_minus_direct_jet( shortcut->getReadier(),
                                                    shortcut->getCppSlha(),
                                                    coloredDecayer,
                                                    effectiveSquarkMass,
                                                    electroweakDecayer,
                                                    intermediateDecayer,
                                              shortcut->getNeutralinoOne() );
    activeDistributions.push_back( directJetDistribution );

    /* in my infinite wisdom & boundless foresight, I wrote the distribution
     * for a negatively-charged chargino decaying to a W^- which decays to a
     * left-handed muon, given the handedness of the jet. however, here we want
     * the charge-conjugate process.
     */
    leftHandedJetRightHandedAntimuonDistribution
    = new W_minus_handed_muon( shortcut->getReadier(),
                               shortcut->getCppSlha(),
                               coloredDecayer,
                               effectiveSquarkMass,
                               electroweakDecayer,
                               intermediateDecayer,
                               shortcut->getNeutralinoOne(),
                               false /* see above */ );
    rightHandedJetRightHandedAntimuonDistribution
    = new W_minus_handed_muon( shortcut->getReadier(),
                               shortcut->getCppSlha(),
                               coloredDecayer,
                               effectiveSquarkMass,
                               electroweakDecayer,
                               intermediateDecayer,
                               shortcut->getNeutralinoOne(),
                               true /* see above */ );
    activeDistributions.push_back(
                                leftHandedJetRightHandedAntimuonDistribution );
    activeDistributions.push_back(
                               rightHandedJetRightHandedAntimuonDistribution );

    leftHandedJetTauMuonDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                  leftHandedJetRightHandedAntimuonDistribution,
                                     shortcut->getHardMuonFromTau() );
    rightHandedJetTauMuonDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                 rightHandedJetRightHandedAntimuonDistribution,
                                     shortcut->getHardMuonFromTau() );
    leftHandedJetTauPionDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                  leftHandedJetRightHandedAntimuonDistribution,
                                     shortcut->getSoftPionFromTau() );
    rightHandedJetTauPionDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                 rightHandedJetRightHandedAntimuonDistribution,
                                     shortcut->getSoftPionFromTau() );
    activeDistributions.push_back( leftHandedJetTauMuonDistribution );
    activeDistributions.push_back( rightHandedJetTauMuonDistribution );
    activeDistributions.push_back( leftHandedJetTauPionDistribution );
    activeDistributions.push_back( rightHandedJetTauPionDistribution );
  }

  charginoToWCascade::~charginoToWCascade()
  {
    // does nothing.
  }


  bool
  charginoToWCascade::validSignal( int const numberOfJets,
                                   int const numberOfNegativeElectrons,
                                   int const numberOfPositiveElectrons,
                                   int const numberOfNegativeMuons,
                                   int const numberOfPositiveMuons )
  /* this returns true if all the arguments are 0, or all are 0 but for
   * numberOfJets OR numberOfPositiveElectrons OR numberOfPositiveMuons, which
   * is 1, or all are jet but numberOfJets which is 2.
   */
  {
    if( ( 0 == numberOfNegativeElectrons )
        &&
        ( 0 == numberOfNegativeMuons )
        &&
        ( ( ( 0 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 1 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 2 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 0 == numberOfJets )
            &&
            ( 1 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 0 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 1 == numberOfPositiveMuons ) ) ) )
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  void
  charginoToWCascade::calculateAcceptance( acceptanceCutSet* const cuts,
                                    acceptanceValues* const currentAcceptance )
  // this returns the appropriate acceptances multiplied by branching ratios.
  {
    currentAcceptance->setTwoJets( 0.0 );
    currentAcceptance->setOneJetZeroLeptons( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveElectron( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveMuon( 0.0 );
    currentAcceptance->setZeroJetsZeroLeptons( 0.0 );
    cascadeBr
    = ( electroweakDecayer->inspect_direct_decay_handler(
             )->get_branching_ratio_for_exact_match( CppSLHA::PDG_code::W_plus,
                                         CppSLHA::PDG_code::neutralino_one ) );
    // the branching ratios of the W into the various SM fermions are covered
    // by the following code.
    if( lhcFaserGlobal::negligibleBr < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      /* at some point, it'd be nice to break this down into transverse
       * contributions & longitudinal contributions, since actually the
       * Goldstone boson contribution should both be large (only needs the
       * Higgsino component of just 1 of the electroweakinos, rather than both)
       * & has a different (& easier to calculate) distribution.
       *
       * however, I think that for LHC-relevant energies, the Higgsino
       * components of all the neutralinos is still going to be larger than
       * m_tau / m_W, enough that I'm going to ignore the longitudinal
       * contributions.
       */

      directMuonPass
      = ( cascadeBr
          * ( jetLeftHandedness
           * integrateAcceptance( leftHandedJetRightHandedAntimuonDistribution,
                                  cuts->getPrimaryLeptonCut() )
               + ( 1.0 - jetLeftHandedness )
          * integrateAcceptance( rightHandedJetRightHandedAntimuonDistribution,
                                 cuts->getPrimaryLeptonCut() ) ) );
      directMuonFail
      = ( cascadeBr
          * ( 1.0
              - ( jetLeftHandedness
           * integrateAcceptance( leftHandedJetRightHandedAntimuonDistribution,
                                  cuts->getPrimaryLeptonCut() )
                  + ( 1.0 - jetLeftHandedness )
          * integrateAcceptance( rightHandedJetRightHandedAntimuonDistribution,
                                 cuts->getSecondaryLeptonCut() ) ) ) );
      directJetPass = integrateAcceptance( directJetDistribution,
                                           cuts->getJetCut() );
      directJetFail = ( 1.0 - directJetPass );

      configurationBr
      = ( cascadeBr * CppSLHA::PDG_data::W_plus_to_hadrons_BR
                    * directJetPass );
      currentAcceptance->setTwoJets( configurationBr * directJetPass );
      currentAcceptance->setOneJetZeroLeptons( 2.0 * configurationBr
                                                   * directJetFail );
      currentAcceptance->setZeroJetsOnePositiveElectron( directMuonPass
                     * CppSLHA::PDG_data::W_plus_to_neutrino_antielectron_BR );
      currentAcceptance->setZeroJetsOnePositiveMuon( directMuonPass
                         * CppSLHA::PDG_data::W_plus_to_neutrino_antimuon_BR );

      currentAcceptance->setZeroJetsZeroLeptons( cascadeBr
                                      * CppSLHA::PDG_data::W_plus_to_hadrons_BR
                                                 * directJetFail
                                                 * directJetFail
                      + ( CppSLHA::PDG_data::W_plus_to_neutrino_antielectron_BR
                          + CppSLHA::PDG_data::W_plus_to_neutrino_antimuon_BR )
                        * directMuonFail );


      // now W^+ decays to tau antileptons, followed by the decays of the tau
      // antileptons:
      cascadeBr *= CppSLHA::PDG_data::W_plus_to_neutrino_tau_antilepton_BR;
      // left-handed jet:
      configurationBr
      = ( cascadeBr * jetLeftHandedness );
      currentMuonDistribution = leftHandedJetTauMuonDistribution;
      currentPionDistribution = leftHandedJetTauPionDistribution;
      calculateForCurrentConfiguration();

      // right-handed jet:
      configurationBr
      = ( cascadeBr * ( 1.0 - jetLeftHandedness ) );
      currentMuonDistribution = rightHandedJetTauMuonDistribution;
      currentPionDistribution = rightHandedJetTauPionDistribution;
      calculateForCurrentConfiguration();
    }
    // debugging:
    /**std::cout
    << std::endl
    << "cascadeBr = " << cascadeBr;
    std::cout << std::endl;**/
  }

  void
  charginoToWCascade::calculateForCurrentConfiguration()
  {
    antitauAntimuonPass = integrateAcceptance( currentMuonDistribution,
                                          currentCuts->getPrimaryLeptonCut() );
    antitauAntimuonFail = ( 1.0 - integrateAcceptance( currentMuonDistribution,
                                      currentCuts->getSecondaryLeptonCut() ) );
    antitauPionPass = integrateAcceptance( currentPionDistribution,
                                       currentCuts->getJetCut() );
    antitauPionFail = ( 1.0 - antitauPionPass );

    currentAcceptance->addToOneJetZeroLeptons( configurationBr
                          * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
                                               * antitauPionPass );
    currentPass = ( configurationBr * antitauAntimuonPass );
    currentAcceptance->addToZeroJetsOnePositiveElectron( currentPass
                    * CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR );
    currentAcceptance->addToZeroJetsOnePositiveMuon( currentPass
                       * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR  );

    currentAcceptance->addToZeroJetsZeroLeptons( configurationBr
                        * ( CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
                                                 * antitauPionFail
                     + ( CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR
                         + CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR )
                       * antitauAntimuonFail ) );
  }



  charginoToHiggsCascade::charginoToHiggsCascade(
                                leptonAcceptanceParameterSet* const kinematics,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                CppSLHA::particle_property_set const* const electroweakDecayer,
               CppSLHA::particle_property_set const* const intermediateDecayer,
                                        inputHandler const* const shortcut ) :
    electroweakCascade( kinematics,
                        effectiveSquarkMass,
                        coloredDecayer,
                        electroweakDecayer,
                        intermediateDecayer,
                        false,
                        shortcut ),
  intermediateDecayerDecays(
                          intermediateDecayer->inspect_direct_decay_handler() )
  {
    // here firstBr & secondBr don't work as well as just getting the branching
    // ratios from the CppSLHA::particle_property_set pointers.

    directAntimuonDistribution
    = new negatively_charged_Higgs_muon( shortcut->getReadier(),
                                         shortcut->getCppSlha(),
                                         coloredDecayer,
                                         effectiveSquarkMass,
                                         electroweakDecayer,
                                         intermediateDecayer,
                                         shortcut->getNeutralinoOne() );
    activeDistributions.push_back( directAntimuonDistribution );

    // the spin-0 nature of the boson combined with the left-handed nature of
    // the tau neutrino forces the tau antilepton to be left-handed.
    antitauAntimuonDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     directAntimuonDistribution,
                                     shortcut->getSoftMuonFromTau() );
    antitauPionDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     directAntimuonDistribution,
                                     shortcut->getHardPionFromTau() );
    activeDistributions.push_back( antitauAntimuonDistribution );
    activeDistributions.push_back( antitauPionDistribution );
  }

  charginoToHiggsCascade::~charginoToHiggsCascade()
  {
    // does nothing.
  }


  bool
  charginoToHiggsCascade::validSignal( int const numberOfJets,
                                       int const numberOfNegativeElectrons,
                                       int const numberOfPositiveElectrons,
                                       int const numberOfNegativeMuons,
                                       int const numberOfPositiveMuons )
  /* this returns true if all the arguments are 0, or all are 0 but for
   * numberOfJets OR numberOfPositiveElectrons OR numberOfPositiveMuons, which
   * is 1, or all are jet but numberOfJets which is 2.
   */
  {
    if( ( 0 == numberOfNegativeElectrons )
        &&
        ( 0 == numberOfNegativeMuons )
        &&
        ( ( ( 0 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 1 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 2 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 0 == numberOfJets )
            &&
            ( 1 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 0 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 1 == numberOfPositiveMuons ) ) ) )
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  void
  charginoToHiggsCascade::calculateAcceptance( acceptanceCutSet* const cuts,
                                    acceptanceValues* const currentAcceptance )
  // this returns the appropriate acceptances multiplied by branching ratios.
  {
    currentAcceptance->setTwoJets( 0.0 );
    currentAcceptance->setOneJetZeroLeptons( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveElectron( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveMuon( 0.0 );
    currentAcceptance->setZeroJetsZeroLeptons( 0.0 );
    cascadeBr
    = ( electroweakDecayer->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                           intermediateDecayer->get_PDG_code(),
                                         CppSLHA::PDG_code::neutralino_one ) );
    // the branching ratios of the W into the various SM fermions are covered
    // by the following code.
    if( lhcFaserGlobal::negligibleBr < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      currentBrToHadrons
      = ( intermediateDecayerDecays->get_branching_ratio_for_exact_match(
                                                         CppSLHA::PDG_code::up,
                                                  -(CppSLHA::PDG_code::down ) )
          + intermediateDecayerDecays->get_branching_ratio_for_exact_match(
                                                         CppSLHA::PDG_code::up,
                                               -(CppSLHA::PDG_code::strange ) )
          + intermediateDecayerDecays->get_branching_ratio_for_exact_match(
                                                         CppSLHA::PDG_code::up,
                                                -(CppSLHA::PDG_code::bottom ) )
          + intermediateDecayerDecays->get_branching_ratio_for_exact_match(
                                                    CppSLHA::PDG_code::charm,
                                               -(CppSLHA::PDG_code::down ) )
          + intermediateDecayerDecays->get_branching_ratio_for_exact_match(
                                                      CppSLHA::PDG_code::charm,
                                               -(CppSLHA::PDG_code::strange ) )
          + intermediateDecayerDecays->get_branching_ratio_for_exact_match(
                                                      CppSLHA::PDG_code::charm,
                                             -(CppSLHA::PDG_code::bottom ) ) );
      currentBrToAntielectron
      = intermediateDecayerDecays->get_branching_ratio_for_exact_match(
                                          CppSLHA::PDG_code::electron_neutrino,
                                             -(CppSLHA::PDG_code::electron ) );
      currentBrToAntimuon
      = intermediateDecayerDecays->get_branching_ratio_for_exact_match(
                                              CppSLHA::PDG_code::muon_neutrino,
                                                 -(CppSLHA::PDG_code::muon ) );
      currentBrToAntitau
      = intermediateDecayerDecays->get_branching_ratio_for_exact_match(
                                               CppSLHA::PDG_code::tau_neutrino,
                                           -(CppSLHA::PDG_code::tau_lepton ) );

      directMuonPass
      = ( cascadeBr * integrateAcceptance( directAntimuonDistribution,
                                           cuts->getPrimaryLeptonCut() ) );
      directMuonFail
      = ( cascadeBr * ( 1.0 - integrateAcceptance( directAntimuonDistribution,
                                             cuts->getPrimaryLeptonCut() ) ) );
      directJetPass = integrateAcceptance( directAntimuonDistribution,
                                           cuts->getJetCut() );
      directJetFail = ( 1.0 - directJetPass );

      configurationBr
      = ( cascadeBr * currentBrToHadrons * directJetPass );
      currentAcceptance->setTwoJets( configurationBr * directJetPass );
      currentAcceptance->setOneJetZeroLeptons( 2.0 * configurationBr
                                                   * directJetFail );
      currentAcceptance->setZeroJetsOnePositiveElectron( directMuonPass
                     * currentBrToAntielectron );
      currentAcceptance->setZeroJetsOnePositiveMuon( directMuonPass
                         * currentBrToAntimuon );

      currentAcceptance->setZeroJetsZeroLeptons( cascadeBr
                                                 * currentBrToHadrons
                                                 * directJetFail
                                                 * directJetFail
                                                 + ( currentBrToAntielectron
                                                     + currentBrToAntimuon )
                                                   * directMuonFail );

      // now decays to tau antileptons, followed by the decays of the tau
      // antileptons:
      cascadeBr *= currentBrToAntitau;
      antitauAntimuonPass = integrateAcceptance( antitauAntimuonDistribution,
                                          currentCuts->getPrimaryLeptonCut() );
      antitauAntimuonFail
      = ( 1.0 - integrateAcceptance( antitauAntimuonDistribution,
                                     currentCuts->getSecondaryLeptonCut() ) );
      antitauPionPass = integrateAcceptance( antitauPionDistribution,
                                             currentCuts->getJetCut() );
      antitauPionFail = ( 1.0 - antitauPionPass );

      currentAcceptance->addToOneJetZeroLeptons( configurationBr
                          * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
                                                 * antitauPionPass );
      currentPass = ( configurationBr * antitauAntimuonPass );
      currentAcceptance->addToZeroJetsOnePositiveElectron( currentPass
                    * CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR );
      currentAcceptance->addToZeroJetsOnePositiveMuon( currentPass
                        * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR );

      currentAcceptance->addToZeroJetsZeroLeptons( configurationBr
                        * ( CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
                                                   * antitauPionFail
                     + ( CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR
                         + CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR )
                       * antitauAntimuonFail ) );
    }
    // debugging:
    /**std::cout
    << std::endl
    << "cascadeBr = " << cascadeBr;
    std::cout << std::endl;**/
  }



  charginoVirtualCascade::charginoVirtualCascade(
                                leptonAcceptanceParameterSet* const kinematics,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                CppSLHA::particle_property_set const* const electroweakDecayer,
                                        inputHandler const* const shortcut ) :
    electroweakCascade( kinematics,
                        effectiveSquarkMass,
                        coloredDecayer,
                        electroweakDecayer,
                        NULL,
                        false,
                        shortcut ),
    electroweakDecayerDecays(
                           electroweakDecayer->inspect_direct_decay_handler() )
  {
    negativeTauLeftHandedness = 0.5;
    // this will probably have to be changed if we ever do the off-shell decays
    // properly.

    // here firstBr & secondBr don't work as well as just getting the branching
    // ratios from the CppSLHA::particle_property_set pointers.

    directDownUpDistribution
    = new chargino_three_body_decay( shortcut->getReadier(),
                                     shortcut->getCppSlha(),
                                     coloredDecayer,
                                     effectiveSquarkMass,
                                     electroweakDecayer,
                                     shortcut->getNeutralinoOne(),
                                     shortcut->getSupL(),
                                     shortcut->getSupR(),
                                     shortcut->getSdownL(),
                                     shortcut->getSdownR() );
    directStrangeCharmDistribution
    = new chargino_three_body_decay( shortcut->getReadier(),
                                     shortcut->getCppSlha(),
                                     coloredDecayer,
                                     effectiveSquarkMass,
                                     electroweakDecayer,
                                     shortcut->getNeutralinoOne(),
                                     shortcut->getScharmL(),
                                     shortcut->getScharmR(),
                                     shortcut->getSstrangeL(),
                                     shortcut->getSstrangeR() );
    directElectronDistribution
    = new neutralino_three_body_decay( shortcut->getReadier(),
                                       shortcut->getCppSlha(),
                                       coloredDecayer,
                                       effectiveSquarkMass,
                                       electroweakDecayer,
                                       shortcut->getNeutralinoOne(),
                                       shortcut->getElectronSneutrinoL(),
                                       shortcut->getElectronSneutrinoR(),
                                       shortcut->getSelectronL(),
                                       shortcut->getSelectronR() );
    directMuonDistribution
    = new neutralino_three_body_decay( shortcut->getReadier(),
                                       shortcut->getCppSlha(),
                                       coloredDecayer,
                                       effectiveSquarkMass,
                                       electroweakDecayer,
                                       shortcut->getNeutralinoOne(),
                                       shortcut->getMuonSneutrinoL(),
                                       shortcut->getMuonSneutrinoR(),
                                       shortcut->getSmuonL(),
                                       shortcut->getSmuonR() );
    directTauDistribution
    = new neutralino_three_body_decay( shortcut->getReadier(),
                                       shortcut->getCppSlha(),
                                       coloredDecayer,
                                       effectiveSquarkMass,
                                       electroweakDecayer,
                                       shortcut->getNeutralinoOne(),
                                       shortcut->getTauSneutrinoL(),
                                       shortcut->getTauSneutrinoR(),
                                       shortcut->getStauOne(),
                                       shortcut->getStauTwo() );
    activeDistributions.push_back( directDownUpDistribution );
    activeDistributions.push_back( directStrangeCharmDistribution );
    activeDistributions.push_back( directElectronDistribution );
    activeDistributions.push_back( directMuonDistribution );
    activeDistributions.push_back( directTauDistribution );

    hardMuonDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     directTauDistribution,
                                     shortcut->getHardMuonFromTau() );
    softMuonDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     directTauDistribution,
                                     shortcut->getSoftMuonFromTau() );
    hardPionDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     directTauDistribution,
                                     shortcut->getHardPionFromTau() );
    softPionDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     directTauDistribution,
                                     shortcut->getSoftPionFromTau() );
    activeDistributions.push_back( hardMuonDistribution );
    activeDistributions.push_back( softMuonDistribution );
    activeDistributions.push_back( hardPionDistribution );
    activeDistributions.push_back( softPionDistribution );
  }

  charginoVirtualCascade::~charginoVirtualCascade()
  {
    // does nothing.
  }


  bool
  charginoVirtualCascade::validSignal( int const numberOfJets,
                                       int const numberOfNegativeElectrons,
                                       int const numberOfPositiveElectrons,
                                       int const numberOfNegativeMuons,
                                       int const numberOfPositiveMuons )
  /* this returns true if all the arguments are 0, or all are 0 but for
   * numberOfJets OR numberOfPositiveElectrons OR numberOfPositiveMuons, which
   * is 1, or all are jet but numberOfJets which is 2.
   */
  {
    if( ( 0 == numberOfNegativeElectrons )
        &&
        ( 0 == numberOfNegativeMuons )
        &&
        ( ( ( 0 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 1 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 2 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 0 == numberOfJets )
            &&
            ( 1 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 0 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 1 == numberOfPositiveMuons ) ) ) )
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  void
  charginoVirtualCascade::calculateAcceptance( acceptanceCutSet* const cuts,
                                    acceptanceValues* const currentAcceptance )
  // this returns the appropriate acceptances multiplied by branching ratios.
  {
    currentAcceptance->setTwoJets( 0.0 );
    currentAcceptance->setOneJetZeroLeptons( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveElectron( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveMuon( 0.0 );
    currentAcceptance->setZeroJetsZeroLeptons( 0.0 );

    // we don't consider flavor non-diagonal decays: it's a lot of effort for
    // little gain.
    cascadeBr
    =  electroweakDecayerDecays->get_branching_ratio_for_exact_match(
                                             CppSLHA::PDG_code::neutralino_one,
                                                         CppSLHA::PDG_code::up,
                                                  -(CppSLHA::PDG_code::down) );
    if( lhcFaserGlobal::negligibleBr < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      directPass = integrateAcceptance( directDownUpDistribution,
                                        cuts->getJetCut() );
      directFail = ( 1.0 - directPass );
      currentAcceptance->addToTwoJets( cascadeBr * directPass * directPass );
      currentAcceptance->addToOneJetZeroLeptons( 2.0 * cascadeBr * directPass
                                                                 * directFail );
      currentAcceptance->addToZeroJetsZeroLeptons( cascadeBr * directFail );
    }
    cascadeBr
    =  electroweakDecayerDecays->get_branching_ratio_for_exact_match(
                                             CppSLHA::PDG_code::neutralino_one,
                                                      CppSLHA::PDG_code::charm,
                                               -(CppSLHA::PDG_code::strange) );
    if( lhcFaserGlobal::negligibleBr < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      directPass = integrateAcceptance( directStrangeCharmDistribution,
                                        cuts->getJetCut() );
      directFail = ( 1.0 - directPass );
      currentAcceptance->addToTwoJets( cascadeBr * directPass * directPass );
      currentAcceptance->addToOneJetZeroLeptons( 2.0 * cascadeBr * directPass
                                                                 * directFail );
      currentAcceptance->addToZeroJetsZeroLeptons( cascadeBr * directFail
                                                             * directFail );
    }
    cascadeBr
    =  electroweakDecayerDecays->get_branching_ratio_for_exact_match(
                                             CppSLHA::PDG_code::neutralino_one,
                                          CppSLHA::PDG_code::electron_neutrino,
                                              -(CppSLHA::PDG_code::electron) );
    if( lhcFaserGlobal::negligibleBr < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      directPass = integrateAcceptance( directElectronDistribution,
                                        cuts->getPrimaryLeptonCut() );
      directFail = ( 1.0 - integrateAcceptance( directElectronDistribution,
                                             cuts->getSecondaryLeptonCut() ) );
      currentAcceptance->addToZeroJetsOnePositiveElectron( cascadeBr
                                                           * directPass );
      currentAcceptance->addToZeroJetsZeroLeptons( cascadeBr * directFail );
    }
    cascadeBr
    =  electroweakDecayerDecays->get_branching_ratio_for_exact_match(
                                             CppSLHA::PDG_code::neutralino_one,
                                              CppSLHA::PDG_code::muon_neutrino,
                                                  -(CppSLHA::PDG_code::muon) );
    if( lhcFaserGlobal::negligibleBr < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      directPass = integrateAcceptance( directMuonDistribution,
                                        cuts->getPrimaryLeptonCut() );
      directFail = ( 1.0 - integrateAcceptance( directMuonDistribution,
                                             cuts->getSecondaryLeptonCut() ) );
      currentAcceptance->addToZeroJetsOnePositiveMuon( cascadeBr
                                                       * directPass );
      currentAcceptance->addToZeroJetsZeroLeptons( cascadeBr
                                                   * directFail );
    }
    cascadeBr
    =  electroweakDecayerDecays->get_branching_ratio_for_exact_match(
                                             CppSLHA::PDG_code::neutralino_one,
                                               CppSLHA::PDG_code::tau_neutrino,
                                            -(CppSLHA::PDG_code::tau_lepton) );
    if( lhcFaserGlobal::negligibleBr < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      hardMuonPass = integrateAcceptance( hardMuonDistribution,
                                          cuts->getPrimaryLeptonCut() );
      hardMuonFail = ( 1.0 - integrateAcceptance( hardMuonDistribution,
                                             cuts->getSecondaryLeptonCut() ) );
      softMuonPass = integrateAcceptance( hardMuonDistribution,
                                          cuts->getPrimaryLeptonCut() );
      softMuonFail = ( 1.0 - integrateAcceptance( hardMuonDistribution,
                                             cuts->getSecondaryLeptonCut() ) );
      averageMuonPass = ( 0.5 * ( hardMuonPass + softMuonPass ) );
      averageMuonFail = ( 0.5 * ( hardMuonFail + softMuonFail ) );
      hardPionPass = integrateAcceptance( hardPionDistribution,
                                          cuts->getJetCut() );
      hardPionFail = ( 1.0 - hardPionPass );
      softPionPass = integrateAcceptance( softPionDistribution,
                                          cuts->getJetCut() );
      softPionFail = ( 1.0 - softPionPass );
      averagePionPass = ( 0.5 * ( hardPionPass + softPionPass ) );
      averagePionFail = ( 0.5 * ( hardPionFail + softPionFail ) );

      currentAcceptance->addToOneJetZeroLeptons( configurationBr
                            * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
                                                 * averagePionPass );
      currentPass = ( configurationBr * averageMuonPass );
      currentAcceptance->addToZeroJetsOnePositiveElectron( currentPass
                      * CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR );
      currentAcceptance->addToZeroJetsOnePositiveMuon( currentPass
                         * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR  );

      currentAcceptance->addToZeroJetsZeroLeptons( configurationBr
                        * ( CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
                                                   * averagePionFail
                     + ( CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR
                         + CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR )
                         * averageMuonFail ) );
    }
    // debugging:
    /**std::cout
    << std::endl
    << "cascadeBr = " << cascadeBr;
    std::cout << std::endl;**/
  }



  scoloredToWPlusScoloredCascade::scoloredToWPlusScoloredCascade(
                                leptonAcceptanceParameterSet* const kinematics,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                   CppSLHA::particle_property_set const* const lighterScolored,
                                        inputHandler const* const shortcut ) :
    electroweakCascade( kinematics,
                        effectiveSquarkMass,
                        coloredDecayer,
                        lighterScolored,
                        shortcut->getWPlus(),
                        false,
                        shortcut )
  {
    // here firstBr & secondBr don't work as well as just getting the branching
    // ratios from the CppSLHA::particle_property_set pointers.

    directMuonDistribution
    = new vector_from_squark_to_muon( shortcut->getReadier(),
                                       shortcut->getCppSlha(),
                                       coloredDecayer,
                                       effectiveSquarkMass,
                                       lighterScolored,
                                       intermediateDecayer );
    activeDistributions.push_back( directMuonDistribution );

    antitauAntimuonDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     directMuonDistribution,
                                     shortcut->getHardMuonFromTau() );
    antitauPionDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     directMuonDistribution,
                                     shortcut->getSoftPionFromTau() );
    activeDistributions.push_back( antitauAntimuonDistribution );
    activeDistributions.push_back( antitauPionDistribution );
  }

  scoloredToWPlusScoloredCascade::~scoloredToWPlusScoloredCascade()
  {
    // does nothing.
  }


  bool
  scoloredToWPlusScoloredCascade::validSignal( int const numberOfJets,
                                           int const numberOfNegativeElectrons,
                                           int const numberOfPositiveElectrons,
                                               int const numberOfNegativeMuons,
                                              int const numberOfPositiveMuons )
  /* this returns true if all the arguments are 0, or all are 0 but for
   * numberOfJets OR numberOfPositiveElectrons OR numberOfPositiveMuons, which
   * is 1, or all are jet but numberOfJets which is 2.
   */
  {
    if( ( 0 == numberOfNegativeElectrons )
        &&
        ( 0 == numberOfNegativeMuons )
        &&
        ( ( ( 0 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 1 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 2 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 0 == numberOfJets )
            &&
            ( 1 == numberOfPositiveElectrons )
            &&
            ( 0 == numberOfPositiveMuons ) )
          ||
          ( ( 0 == numberOfJets )
            &&
            ( 0 == numberOfPositiveElectrons )
            &&
            ( 1 == numberOfPositiveMuons ) ) ) )
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  void
  scoloredToWPlusScoloredCascade::calculateAcceptance(
                                                  acceptanceCutSet* const cuts,
                                    acceptanceValues* const currentAcceptance )
  // this returns the appropriate acceptances multiplied by branching ratios.
  {
    currentAcceptance->setTwoJets( 0.0 );
    currentAcceptance->setOneJetZeroLeptons( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveElectron( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveMuon( 0.0 );
    currentAcceptance->setZeroJetsZeroLeptons( 0.0 );
    cascadeBr
    = ( electroweakDecayer->inspect_direct_decay_handler(
             )->get_branching_ratio_for_exact_match( CppSLHA::PDG_code::W_plus,
                                             electroweakDecayer->get_PDG_code()
      /* the lighter scolored got put into electroweakDecayer, while the W
       * boson got put into intermediateDecayer. */ ) );
    // the branching ratios of the W into the various SM fermions are covered
    // by the following code.
    if( lhcFaserGlobal::negligibleBr < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      /* at some point, it'd be nice to break this down into transverse
       * contributions & longitudinal contributions, since actually the
       * Goldstone boson contribution could both be large (for sbottom -> stop,
       * for example) & has a different (& easier to calculate) distribution.
       *
       * however, ( ( m_top * m_tau ) / m_W^2 ) is still pretty small, so it
       * should be fine to neglect the longitudinal contributions.
       */

      directMuonPass
      = ( cascadeBr * integrateAcceptance( directMuonDistribution,
                                           cuts->getPrimaryLeptonCut() ) );
      directMuonFail
      = ( cascadeBr * ( 1.0 - integrateAcceptance( directMuonDistribution,
                                           cuts->getSecondaryLeptonCut() ) ) );
      directJetPass = integrateAcceptance( directMuonDistribution,
                                           cuts->getJetCut() );
      directJetFail = ( 1.0 - directJetPass );

      configurationBr
      = ( cascadeBr * CppSLHA::PDG_data::W_plus_to_hadrons_BR
                    * directJetPass );
      currentAcceptance->setTwoJets( configurationBr * directJetPass );
      currentAcceptance->setOneJetZeroLeptons( 2.0 * configurationBr
                                                   * directJetFail );
      currentAcceptance->setZeroJetsOnePositiveElectron( directMuonPass
                     * CppSLHA::PDG_data::W_plus_to_neutrino_antielectron_BR );
      currentAcceptance->setZeroJetsOnePositiveMuon( directMuonPass
                         * CppSLHA::PDG_data::W_plus_to_neutrino_antimuon_BR );

      currentAcceptance->setZeroJetsZeroLeptons( cascadeBr
                                      * CppSLHA::PDG_data::W_plus_to_hadrons_BR
                                                 * directJetFail
                                                 * directJetFail
                            + CppSLHA::PDG_data::W_plus_to_neutrino_antimuon_BR
                              * directMuonFail );

      // now W^+ decays to tau antileptons, followed by the decays of the tau
      // antileptons:
      cascadeBr *= CppSLHA::PDG_data::W_plus_to_neutrino_tau_antilepton_BR;
      antitauAntimuonPass = integrateAcceptance( antitauAntimuonDistribution,
                                          currentCuts->getPrimaryLeptonCut() );
      antitauAntimuonFail
      = ( 1.0 - integrateAcceptance( antitauAntimuonDistribution,
                                     currentCuts->getSecondaryLeptonCut() ) );
      antitauPionPass = integrateAcceptance( antitauPionDistribution,
                                             currentCuts->getJetCut() );
      antitauPionFail = ( 1.0 - antitauPionPass );

      currentAcceptance->addToOneJetZeroLeptons( configurationBr
                          * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
                                                 * antitauPionPass );
      currentPass = ( configurationBr * antitauAntimuonPass );
      currentAcceptance->addToZeroJetsOnePositiveElectron( currentPass
                      * CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR );
      currentAcceptance->addToZeroJetsOnePositiveMuon( currentPass
                         * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR  );

      currentAcceptance->addToZeroJetsZeroLeptons( configurationBr
                        * ( CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
                                                   * antitauPionFail
                     + ( CppSLHA::PDG_data::tau_lepton_to_neutrinos_electron_BR
                         + CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR )
                         * antitauAntimuonFail ) );
    }
    // debugging:
    /**std::cout
    << std::endl
    << "cascadeBr = " << cascadeBr;
    std::cout << std::endl;**/
  }
}  // end of LHC_FASER namespace.
