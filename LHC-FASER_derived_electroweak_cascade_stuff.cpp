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



  lightestNeutralinoCascade::lightestNeutralinoCascade() :
    electroweakCascade( NULL,
                        NULL,
                        false,
                        NULL,
                        false,
                        NULL,
                        false,
                        NULL )
  {
    acceptances.setComparison( &acceptanceCutSet::justReturnTrue );
  }

  lightestNeutralinoCascade::~lightestNeutralinoCascade()
  {
    // does nothing.
  }



  neutralinoToSemuCascade::neutralinoToSemuCascade(
                                leptonAcceptanceParameterSet* const kinematics,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                                  //bool const coloredDecayerIsNotAntiparticle,
                CppSLHA::particle_property_set const* const electroweakDecayer,
                              //bool const electroweakDecayerIsNotAntiparticle,
               CppSLHA::particle_property_set const* const intermediateDecayer,
                                        inputHandler const* const shortcut ) :
    electroweakCascade( kinematics,
                        effectiveSquarkMass,
                        coloredDecayer,
                        //coloredDecayerIsNotAntiparticle,
                        electroweakDecayer,
                        //electroweakDecayerIsNotAntiparticle,
                        intermediateDecayer,
                        true,
                        shortcut )
  {
    if( ( CppSLHA::PDG_code::smuon_L == intermediateDecayer->get_PDG_code() )
        ||
        ( CppSLHA::PDG_code::smuon_R == intermediateDecayer->get_PDG_code() ) )
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
    nearSameDistribution
    = new same_chirality_near_muon( shortcut->getReadier(),
                                    shortcut->getCppSlha(),
                                    coloredDecayer,
                                    effectiveSquarkMass,
                                    electroweakDecayer,
                                    intermediateDecayer );
    nearOppositeDistribution
    = new opposite_chirality_near_muon( shortcut->getReadier(),
                                        shortcut->getCppSlha(),
                                        coloredDecayer,
                                        effectiveSquarkMass,
                                        electroweakDecayer,
                                        intermediateDecayer );
    farSameDistribution
    = new same_chirality_far_muon( shortcut->getReadier(),
                                   shortcut->getCppSlha(),
                                   coloredDecayer,
                                   effectiveSquarkMass,
                                   electroweakDecayer,
                                   intermediateDecayer,
                                   shortcut->getNeutralinoOne() );
    farOppositeDistribution
    = new opposite_chirality_far_muon( shortcut->getReadier(),
                                       shortcut->getCppSlha(),
                                       coloredDecayer,
                                       effectiveSquarkMass,
                                       electroweakDecayer,
                                       intermediateDecayer,
                                       shortcut->getNeutralinoOne() );
    activeDistributions.push_back( nearSameDistribution );
    activeDistributions.push_back( nearOppositeDistribution );
    activeDistributions.push_back( farSameDistribution );
    activeDistributions.push_back( farOppositeDistribution );
  }

  neutralinoToSemuCascade::~neutralinoToSemuCascade()
  {
    // does nothing.
  }


  bool
  neutralinoToSemuCascade::validSignal( int const numberOfJets,
                                        int const numberOfNegativeElectrons,
                                        int const numberOfPositiveElectrons,
                                        int const numberOfNegativeMuons,
                                        int const numberOfPositiveMuons )
  // this returns true if numberOfJets == 0, & either 1 of or both signs of
  // electron OR muon are asked for, or zero leptons are asked for.
  {
    if( ( 0 == numberOfJets )
        // if we are asked for 0 jets...
        &&
        ( 0 <= numberOfNegativeElectrons )
        &&
        ( 0 <= numberOfPositiveElectrons )
        &&
        ( 0 <= numberOfNegativeMuons )
        &&
        ( 0 <= numberOfPositiveMuons )
        &&
        ( 1 >= ( numberOfNegativeElectrons + numberOfNegativeMuons ) )
        &&
        ( 1 >= ( numberOfNegativeElectrons + numberOfPositiveMuons ) )
        &&
        ( 1 >= ( numberOfPositiveElectrons + numberOfNegativeMuons ) )
        &&
        ( 1 >= ( numberOfPositiveElectrons + numberOfPositiveMuons ) ) )
        // AND if we are asked for EITHER electrons OR muons but NOT both...
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  void
  neutralinoToSemuCascade::calculateAcceptance( acceptanceCutSet* const cuts,
                                    acceptanceValues* const currentAcceptance )
  // this returns the appropriate acceptances multiplied by branching ratios
  // from the electroweakino through the selectron or smuon to the LSP.
  {
    cascadeBr = ( first_BR->getBr() * second_BR->getBr() );
    if( lhcFaserGlobal::negligibleBr < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      jetLeftHandedness
      = shortcut->quarkOrLeptonLeftHandedness( coloredDecayer->get_PDG_code(),
                                          electroweakDecayer->get_PDG_code() );
      nearNegativeLeptonLeftHandedness
      = shortcut->quarkOrLeptonLeftHandedness(
                                           intermediateDecayer->get_PDG_code(),
                                          electroweakDecayer->get_PDG_code() );
      nearNegativeLeptonSameHandednessTimesBr
      = ( cascadeBr * ( ( jetLeftHandedness
                          * nearNegativeLeptonLeftHandedness )
                        + ( ( 1.0 - jetLeftHandedness )
                            * ( 1.0 - nearNegativeLeptonLeftHandedness ) ) ) );

      nearSamePass = integrateAcceptance( nearSameDistribution,
                                          cuts->getPrimaryLeptonCut() );
      nearSameFail = ( 1.0 - integrateAcceptance( nearSameDistribution,
                                             cuts->getSecondaryLeptonCut() ) );
      farSamePass = integrateAcceptance( farSameDistribution,
                                         cuts->getPrimaryLeptonCut() );
      farSameFail = ( 1.0 - integrateAcceptance( farSameDistribution,
                                             cuts->getSecondaryLeptonCut() ) );
      nearOppositePass = integrateAcceptance( nearOppositeDistribution,
                                          cuts->getPrimaryLeptonCut() );
      nearOppositeFail = ( 1.0 - integrateAcceptance( nearOppositeDistribution,
                                             cuts->getSecondaryLeptonCut() ) );
      farOppositePass = integrateAcceptance( farOppositeDistribution,
                                         cuts->getPrimaryLeptonCut() );
      farOppositeFail = ( 1.0 - integrateAcceptance( farOppositeDistribution,
                                             cuts->getSecondaryLeptonCut() ) );

      currentAcceptance->setOssfMinusOsdf( ( cascadeBr * nearSamePass
                                                       * farSamePass
                                             + ( cascadeBr
                                    - nearNegativeLeptonSameHandednessTimesBr )
                                               * nearOppositePass
                                               * farOppositePass ) );
      currentAcceptance->setZeroJetsZeroLeptons( ( cascadeBr * nearSameFail
                                                             * farSameFail
                                                 + ( cascadeBr
                                    - nearNegativeLeptonSameHandednessTimesBr )
                                                   * nearOppositeFail
                                                   * farOppositeFail ) );
      justOneNegativeLeptonPass
      = ( nearNegativeLeptonSameHandednessTimesBr
          * ( nearSamePass * farSameFail + farOppositePass * nearOppositeFail )
          + ( cascadeBr - nearNegativeLeptonSameHandednessTimesBr )
            * ( farSamePass * nearSameFail
                + nearOppositePass * farOppositeFail ) );
      justOnePositiveLeptonPass
      = ( nearNegativeLeptonSameHandednessTimesBr
          * ( farSamePass * nearSameFail + nearOppositePass * farOppositeFail )
          + ( cascadeBr - nearNegativeLeptonSameHandednessTimesBr )
            * ( nearSamePass * farSameFail
                + farOppositePass * nearOppositeFail ) );

      if( muonsNotElectrons )
      {
        currentAcceptance->setMuonPlusAntimuon(
                                       currentAcceptance->getOssfMinusOsdf() );
        currentAcceptance->setZeroJetsOneNegativeMuon(
                                                   justOneNegativeLeptonPass );
        currentAcceptance->setZeroJetsOnePositiveMuon(
                                                   justOnePositiveLeptonPass );
      }
      else
      {
        currentAcceptance->setElectronPlusAntielectron(
                                       currentAcceptance->getOssfMinusOsdf() );
        currentAcceptance->setZeroJetsOneNegativeElectron(
                                                   justOneNegativeLeptonPass );
        currentAcceptance->setZeroJetsOnePositiveElectron(
                                                   justOnePositiveLeptonPass );
      }
    }  // end of if cascade branching ratio was large enough to care about.
    // debugging:
    /**std::cout
    << std::endl
    << "cascadeBr = " << cascadeBr;
    std::cout << std::endl;**/
  }



  chargeSummedNeutralinoToSemuCascade::chargeSummedNeutralinoToSemuCascade(
                                leptonAcceptanceParameterSet* const kinematics,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                                    bool const coloredDecayerIsNotAntiparticle,
                CppSLHA::particle_property_set const* const electroweakDecayer,
                                bool const electroweakDecayerIsNotAntiparticle,
               CppSLHA::particle_property_set const* const intermediateDecayer,
                                        inputHandler const* const shortcut ) :
    electroweakCascade( kinematics,
                        effectiveSquarkMass,
                        coloredDecayer,
                        coloredDecayerIsNotAntiparticle,
                        electroweakDecayer,
                        electroweakDecayerIsNotAntiparticle,
                        intermediateDecayer,
                        true,
                        shortcut ),
    nearDistribution( NULL ),
    farDistribution( NULL )
  {
    if( ( CppSLHA::PDG_code::smuon_L == intermediateDecayer->get_PDG_code() )
        ||
        ( CppSLHA::PDG_code::smuon_R == intermediateDecayer->get_PDG_code() ) )
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
    nearDistribution
    = new same_chirality_near_muon( shortcut->getReadier(),
                                    shortcut->getCppSlha(),
                                    coloredDecayer,
                                    effectiveSquarkMass,
                                    electroweakDecayer,
                                    intermediateDecayer );
    farDistribution
    = new same_chirality_far_muon( shortcut->getReadier(),
                                   shortcut->getCppSlha(),
                                   coloredDecayer,
                                   effectiveSquarkMass,
                                   electroweakDecayer,
                                   intermediateDecayer,
                                   shortcut->getNeutralinoOne() );
    activeDistributions.push_back( nearDistribution );
    activeDistributions.push_back( farDistribution );
  }

  chargeSummedNeutralinoToSemuCascade::~chargeSummedNeutralinoToSemuCascade()
  {
    // does nothing.
  }


  bool
  chargeSummedNeutralinoToSemuCascade::validSignal( int const numberOfJets,
                                           int const numberOfNegativeElectrons,
                                           int const numberOfPositiveElectrons,
                                               int const numberOfNegativeMuons,
                                              int const numberOfPositiveMuons )
  // this returns true if numberOfJets == 0, & either 1 of or both signs of
  // electron OR muon are asked for, or zero leptons are asked for.
  {
    if( ( 0 == numberOfJets )
        // if we are asked for 0 jets...
        &&
        ( 0 <= numberOfNegativeElectrons )
        &&
        ( 0 <= numberOfPositiveElectrons )
        &&
        ( 0 <= numberOfNegativeMuons )
        &&
        ( 0 <= numberOfPositiveMuons )
        &&
        ( 1 >= ( numberOfNegativeElectrons + numberOfNegativeMuons ) )
        &&
        ( 1 >= ( numberOfNegativeElectrons + numberOfPositiveMuons ) )
        &&
        ( 1 >= ( numberOfPositiveElectrons + numberOfNegativeMuons ) )
        &&
        ( 1 >= ( numberOfPositiveElectrons + numberOfPositiveMuons ) ) )
        // AND if we are asked for EITHER electrons OR muons but NOT both...
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  void
  chargeSummedNeutralinoToSemuCascade::calculateAcceptance(
                                                  acceptanceCutSet* const cuts,
                                    acceptanceValues* const currentAcceptance )
  // this returns the appropriate acceptances multiplied by branching ratios
  // from the electroweakino through the selectron or smuon to the LSP.
  {
    cascadeBr = ( first_BR->getBr() * second_BR->getBr() );
    // this is just for *1* of the possibilities (e.g. negatively-charged
    // slepton) - the charge-conjugate is accounted for with factors of 2.
    if( lhcFaserGlobal::negligibleBr < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      nearPass = integrateAcceptance( nearDistribution,
                                      cuts->getPrimaryLeptonCut() );
      nearFail = ( 1.0 - integrateAcceptance( nearDistribution,
                                             cuts->getSecondaryLeptonCut() ) );
      farPass = integrateAcceptance( farDistribution,
                                     cuts->getPrimaryLeptonCut() );
      farFail = ( 1.0 - integrateAcceptance( farDistribution,
                                             cuts->getSecondaryLeptonCut() ) );
      bothPass = ( cascadeBr * ( 2.0 * nearPass * farPass ) );
      currentAcceptance->setOssfMinusOsdf( bothPass );
      // the factor of 2.0 accounts for the charge-conjugate cascade.
      currentAcceptance->setZeroJetsZeroLeptons(
                                ( cascadeBr * ( 2.0 * nearFail * farFail ) ) );
      // the factor of 2.0 accounts for the charge-conjugate cascade.
      onePassOneFail
      = ( cascadeBr * ( nearPass * farFail + nearFail * farPass ) );

      if( muonsNotElectrons )
      {
        currentAcceptance->setMuonPlusAntimuon( bothPass );
        currentAcceptance->setZeroJetsOneNegativeMuon( onePassOneFail );
        currentAcceptance->setZeroJetsOnePositiveMuon( onePassOneFail );
      }
      else
      {
        currentAcceptance->setElectronPlusAntielectron( bothPass );
        currentAcceptance->setZeroJetsOneNegativeElectron( onePassOneFail );
        currentAcceptance->setZeroJetsOnePositiveElectron( onePassOneFail );
      }
    }  // end of if cascade branching ratio was large enough to care about.
    // debugging:
    /**std::cout
    << std::endl
    << "cascade_BR = " << cascade_BR;
    std::cout << std::endl;**/
  }



  neutralinoToStauCascade::neutralinoToStauCascade(
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
                        shortcut )
  {
    // we set up the BR of the (only-ish) channel:
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
    nearSameTauDistribution
    = new same_chirality_near_muon( shortcut->getReadier(),
                                    shortcut->getCppSlha(),
                                    coloredDecayer,
                                    effectiveSquarkMass,
                                    electroweakDecayer,
                                    intermediateDecayer );
    nearOppositeTauDistribution
    = new opposite_chirality_near_muon( shortcut->getReadier(),
                                        shortcut->getCppSlha(),
                                        coloredDecayer,
                                        effectiveSquarkMass,
                                        electroweakDecayer,
                                        intermediateDecayer );
    farSameTauDistribution
    = new same_chirality_far_muon( shortcut->getReadier(),
                                   shortcut->getCppSlha(),
                                   coloredDecayer,
                                   effectiveSquarkMass,
                                   electroweakDecayer,
                                   intermediateDecayer,
                                   shortcut->getNeutralinoOne() );
    farOppositeTauDistribution
    = new opposite_chirality_far_muon( shortcut->getReadier(),
                                       shortcut->getCppSlha(),
                                       coloredDecayer,
                                       effectiveSquarkMass,
                                       electroweakDecayer,
                                       intermediateDecayer,
                                       shortcut->getNeutralinoOne() );
    activeDistributions.push_back( nearSameTauDistribution );
    activeDistributions.push_back( nearOppositeTauDistribution );
    activeDistributions.push_back( farSameTauDistribution );
    activeDistributions.push_back( farOppositeTauDistribution );

    nearSameHardMuonDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     nearSameTauDistribution,
                                     shortcut->getHardMuonFromTau() );
    nearSameSoftMuonDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     nearSameTauDistribution,
                                     shortcut->getSoftMuonFromTau() );
    nearSameHardPionDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     nearSameTauDistribution,
                                     shortcut->getHardPionFromTau() );
    nearSameSoftPionDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     nearSameTauDistribution,
                                     shortcut->getSoftPionFromTau() );
    nearOppositeHardMuonDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     nearOppositeTauDistribution,
                                     shortcut->getHardMuonFromTau() );
    nearOppositeSoftMuonDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     nearOppositeTauDistribution,
                                     shortcut->getSoftMuonFromTau() );
    nearOppositeHardPionDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     nearOppositeTauDistribution,
                                     shortcut->getHardPionFromTau() );
    nearOppositeSoftPionDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     nearOppositeTauDistribution,
                                     shortcut->getSoftPionFromTau() );
    farSameHardMuonDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     farSameTauDistribution,
                                     shortcut->getHardMuonFromTau() );
    farSameSoftMuonDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     farSameTauDistribution,
                                     shortcut->getSoftMuonFromTau() );
    farSameHardPionDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     farSameTauDistribution,
                                     shortcut->getHardPionFromTau() );
    farSameSoftPionDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     farSameTauDistribution,
                                     shortcut->getSoftPionFromTau() );
    farOppositeHardMuonDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     farOppositeTauDistribution,
                                     shortcut->getHardMuonFromTau() );
    farOppositeSoftMuonDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     farOppositeTauDistribution,
                                     shortcut->getSoftMuonFromTau() );
    farOppositeHardPionDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     farOppositeTauDistribution,
                                     shortcut->getHardPionFromTau() );
    farOppositeSoftPionDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     farOppositeTauDistribution,
                                     shortcut->getSoftPionFromTau() );
    activeDistributions.push_back( nearSameHardMuonDistribution );
    activeDistributions.push_back( nearSameSoftMuonDistribution );
    activeDistributions.push_back( nearSameHardPionDistribution );
    activeDistributions.push_back( nearSameSoftPionDistribution );
    activeDistributions.push_back( nearOppositeHardMuonDistribution );
    activeDistributions.push_back( nearOppositeSoftMuonDistribution );
    activeDistributions.push_back( nearOppositeHardPionDistribution );
    activeDistributions.push_back( nearOppositeSoftPionDistribution );
    activeDistributions.push_back( farSameHardMuonDistribution );
    activeDistributions.push_back( farSameSoftMuonDistribution );
    activeDistributions.push_back( farSameHardPionDistribution );
    activeDistributions.push_back( farSameSoftPionDistribution );
    activeDistributions.push_back( farOppositeHardMuonDistribution );
    activeDistributions.push_back( farOppositeSoftMuonDistribution );
    activeDistributions.push_back( farOppositeHardPionDistribution );
    activeDistributions.push_back( farOppositeSoftPionDistribution );
  }

  neutralinoToStauCascade::~neutralinoToStauCascade()
  {
    // does nothing.
  }


  bool
  neutralinoToStauCascade::validSignal( int const numberOfJets,
                                        int const numberOfNegativeElectrons,
                                        int const numberOfPositiveElectrons,
                                        int const numberOfNegativeMuons,
                                        int const numberOfPositiveMuons )
  /* this returns true if a configuration where each of the signs of tau lepton
   * decayed either into a detected jet, detected lepton, or undetected
   * particle, & false otherwise.
   */
  {
    if( 2 == numberOfJets )
    {
      if( ( 0 == numberOfNegativeElectrons )
          &&
          ( 0 == numberOfPositiveElectrons )
          &&
          ( 0 == numberOfNegativeMuons )
          &&
          ( 0 == numberOfPositiveMuons ) )
      {
        return true;
      }
      else
      {
        return false;
      }
    }
    else if( ( 0 <= numberOfJets )
             &&
             ( 0 <= numberOfNegativeElectrons )
             &&
             ( 0 <= numberOfPositiveElectrons )
             &&
             ( 0 <= numberOfNegativeMuons )
             &&
             ( 0 <= numberOfPositiveMuons )
             &&
             ( 1 >= ( numberOfNegativeElectrons + numberOfNegativeMuons ) )
             &&
             ( 1 >= ( numberOfPositiveElectrons + numberOfPositiveMuons ) ) )
    {
      if( 2 >=
          ( numberOfJets
            + numberOfNegativeElectrons + numberOfPositiveElectrons
            + numberOfNegativeMuons + numberOfPositiveMuons ) )
      {
        return true;
      }
      else
      {
        return false;
      }
    }
    else
    {
      return false;
    }
  }

  void
  neutralinoToStauCascade::calculateAcceptance( acceptanceCutSet* const cuts,
                                    acceptanceValues* const currentAcceptance )
  // this returns the appropriate acceptances multiplied by branching ratios
  // from the electroweakino through the stau to the LSP.
  {
    currentAcceptance->setTwoJets( 0.0 );
    currentAcceptance->setOneJetOneNegativeElectron( 0.0 );
    currentAcceptance->setOneJetOnePositiveElectron( 0.0 );
    currentAcceptance->setOneJetOneNegativeMuon( 0.0 );
    currentAcceptance->setOneJetOnePositiveMuon( 0.0 );
    currentAcceptance->setElectronPlusAntielectron( 0.0 );
    currentAcceptance->setNegativeElectronPlusPositiveMuon( 0.0 );
    currentAcceptance->setNegativeMuonPlusPositiveElectron( 0.0 );
    currentAcceptance->setMuonPlusAntimuon( 0.0 );
    currentAcceptance->setOneJetZeroLeptons( 0.0 );
    currentAcceptance->setZeroJetsOneNegativeElectron( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveElectron( 0.0 );
    currentAcceptance->setZeroJetsOneNegativeMuon( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveMuon( 0.0 );
    currentAcceptance->setZeroJetsZeroLeptons( 0.0 );
    cascadeBr = ( first_BR->getBr() * second_BR->getBr() );
    if( lhcFaserGlobal::negligibleBr < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      currentCuts = cuts;
      this->currentAcceptance = currentAcceptance;
      jetLeftHandedness
      = shortcut->quarkOrLeptonLeftHandedness(
                                                coloredDecayer->get_PDG_code(),
                                          electroweakDecayer->get_PDG_code() );
      jetRightHandedness = ( 1.0 - jetLeftHandedness );
      nearNegativeTauLeftHandedness
      = shortcut->quarkOrLeptonLeftHandedness(
                                            mediating_particle->get_PDG_code(),
                                              decaying_EWino->get_PDG_code() );
      nearNegativeTauRightHandedness = ( 1.0 - nearNegativeTauLeftHandedness );
      farNegativeTauLeftHandedness
      = shortcut->quarkOrLeptonLeftHandedness(
                                            mediating_particle->get_PDG_code(),
                              shortcut->getNeutralinoOne()->get_PDG_code() );
      farNegativeTauRightHandedness = ( 1.0 - farNegativeTauLeftHandedness );


      /* in the following, XY(bar)Z(bar) is quark handedness, near tau
       * lepton handedness, far tau lepton handedness, YbarZ means that
       * the near tau lepton is positive while YZbar means that the far
       * tau lepton is positive.
       */
      // near lepton same handedness as quark:
      // hard near muons, hard far muons:
      currentNearMuonDistribution = nearSameHardMuonDistribution;
      currentNearPionDistribution = nearSameSoftPionDistribution;
      currentFarMuonDistribution = farSameHardMuonDistribution;
      currentFarPionDistribution = farSameSoftPionDistribution;
      // LLRbar
      negativeNearConfigurationBr = ( cascadeBr * jetLeftHandedness
                                                * nearNegativeTauLeftHandedness
                                              * farNegativeTauLeftHandedness );
      // RRbarL
      positiveNearConfigurationBr = ( cascadeBr * jetRightHandedness
                                                * nearNegativeTauLeftHandedness
                                              * farNegativeTauLeftHandedness );
      calculateForCurrentConfiguration();

      // hard near muons, soft far muons:
      currentNearMuonDistribution = nearSameHardMuonDistribution;
      currentNearPionDistribution = nearSameSoftPionDistribution;
      currentFarMuonDistribution = farSameSoftMuonDistribution;
      currentFarPionDistribution = farSameHardPionDistribution;
      // LLLbar
      negativeNearConfigurationBr = ( cascadeBr * jetLeftHandedness
                                                * nearNegativeTauLeftHandedness
                                             * farNegativeTauRightHandedness );
      // RRbarR
      positiveNearConfigurationBr = ( cascadeBr * jetRightHandedness
                                                * nearNegativeTauLeftHandedness
                                             * farNegativeTauRightHandedness );
      calculateForCurrentConfiguration();

      // soft near muons, hard far muons:
      currentNearMuonDistribution = nearSameSoftMuonDistribution;
      currentNearPionDistribution = nearSameHardPionDistribution;
      currentFarMuonDistribution = farSameHardMuonDistribution;
      currentFarPionDistribution = farSameSoftPionDistribution;
      // RRRbar
      negativeNearConfigurationBr = ( cascadeBr * jetRightHandedness
                                               * nearNegativeTauRightHandedness
                                              * farNegativeTauLeftHandedness );
      // LLbarL
      positiveNearConfigurationBr = ( cascadeBr * jetLeftHandedness
                                               * nearNegativeTauRightHandedness
                                              * farNegativeTauLeftHandedness );
      calculateForCurrentConfiguration();

      // soft near muons, soft far muons:
      currentNearMuonDistribution = nearSameSoftMuonDistribution;
      currentNearPionDistribution = nearSameHardPionDistribution;
      currentFarMuonDistribution = farSameSoftMuonDistribution;
      currentFarPionDistribution = farSameHardPionDistribution;
      // RRLbar
      negativeNearConfigurationBr = ( cascadeBr * jetRightHandedness
                                               * nearNegativeTauRightHandedness
                                             * farNegativeTauRightHandedness );
      // LLbarR
      positiveNearConfigurationBr = ( cascadeBr * jetLeftHandedness
                                               * nearNegativeTauRightHandedness
                                             * farNegativeTauRightHandedness );
      calculateForCurrentConfiguration();

      // near lepton opposite handedness to quark:
      // hard near muons, hard far muons:
      currentNearMuonDistribution = nearOppositeHardMuonDistribution;
      currentNearPionDistribution = nearOppositeSoftPionDistribution;
      currentFarMuonDistribution = farOppositeHardMuonDistribution;
      currentFarPionDistribution = farOppositeSoftPionDistribution;
      // RLRbar
      negativeNearConfigurationBr = ( cascadeBr * jetRightHandedness
                                                * nearNegativeTauLeftHandedness
                                              * farNegativeTauLeftHandedness );
      // LRbarL
      positiveNearConfigurationBr = ( cascadeBr * jetLeftHandedness
                                                * nearNegativeTauLeftHandedness
                                              * farNegativeTauLeftHandedness );
      calculateForCurrentConfiguration();

      // hard near muons, soft far muons:
      currentNearMuonDistribution = nearOppositeHardMuonDistribution;
      currentNearPionDistribution = nearOppositeSoftPionDistribution;
      currentFarMuonDistribution = farOppositeSoftMuonDistribution;
      currentFarPionDistribution = farOppositeHardPionDistribution;
      // RLLbar
      negativeNearConfigurationBr = ( cascadeBr * jetRightHandedness
                                                * nearNegativeTauLeftHandedness
                                             * farNegativeTauRightHandedness );
      // LRbarR
      positiveNearConfigurationBr = ( cascadeBr * jetLeftHandedness
                                                * nearNegativeTauLeftHandedness
                                             * farNegativeTauRightHandedness );
      calculateForCurrentConfiguration();

      // soft near muons, hard far muons:
      currentNearMuonDistribution = nearOppositeSoftMuonDistribution;
      currentNearPionDistribution = nearOppositeHardPionDistribution;
      currentFarMuonDistribution = farOppositeHardMuonDistribution;
      currentFarPionDistribution = farOppositeSoftPionDistribution;
      // LRRbar
      negativeNearConfigurationBr = ( cascadeBr * jetLeftHandedness
                                               * nearNegativeTauRightHandedness
                                              * farNegativeTauLeftHandedness );
      // RLbarL
      positiveNearConfigurationBr = ( cascadeBr * jetRightHandedness
                                               * nearNegativeTauRightHandedness
                                              * farNegativeTauLeftHandedness );
      calculateForCurrentConfiguration();

      // soft near muons, soft far muons:
      currentNearMuonDistribution = nearOppositeSoftMuonDistribution;
      currentNearPionDistribution = nearOppositeHardPionDistribution;
      currentFarMuonDistribution = farOppositeSoftMuonDistribution;
      currentFarPionDistribution = farOppositeHardPionDistribution;

      // LRLbar
      negativeNearConfigurationBr = ( cascadeBr * jetLeftHandedness
                                               * nearNegativeTauRightHandedness
                                             * farNegativeTauRightHandedness );
      // RLbarR
      positiveNearConfigurationBr = ( cascadeBr * jetRightHandedness
                                               * nearNegativeTauRightHandedness
                                             * farNegativeTauRightHandedness );
      calculateForCurrentConfiguration();
    }
    // debugging:
    /**std::cout
    << std::endl
    << "cascadeBr = " << cascadeBr;
    std::cout << std::endl;**/
  }

  void
  neutralinoToStauCascade::calculateForCurrentConfiguration()
  {
    nearMuonPass = integrateAcceptance( currentNearMuonDistribution,
                                        currentCuts->getPrimaryLeptonCut() );
    nearMuonFail = ( 1.0 - integrateAcceptance( currentNearMuonDistribution,
                                      currentCuts->getSecondaryLeptonCut() ) );
    nearPionPass = integrateAcceptance( currentNearPionDistribution,
                                        currentCuts->getJetCut() );
    nearPionFail = ( 1.0 - nearPionPass );
    farMuonPass = integrateAcceptance( currentFarMuonDistribution,
                                       currentCuts->getPrimaryLeptonCut() );
    farMuonFail = ( 1.0 - integrateAcceptance( currentFarMuonDistribution,
                                      currentCuts->getSecondaryLeptonCut() ) );
    farPionPass = integrateAcceptance( currentFarPionDistribution,
                                       currentCuts->getJetCut() );
    farPionFail = ( 1.0 - farPionPass );

    currentAcceptance->addToTwoJets( ( negativeNearConfigurationBr
                                       + positiveNearConfigurationBr )
                                     * tauPairToPionPairBr
                                     * nearPionPass * farPionPass );

    currentPass = ( negativeNearConfigurationBr * nearMuonPass * farPionPass
                  + positiveNearConfigurationBr * nearPionPass * farMuonPass );
    currentAcceptance->addToOneJetOneNegativeElectron(
                                 tauToPionTimesTauToElectronBr * currentPass );
    currentAcceptance->addToOneJetOneNegativeMuon( tauToPionTimesTauToMuonBr
                                                   * currentPass );

    currentPass = ( positiveNearConfigurationBr * nearMuonPass * farPionPass
                  + negativeNearConfigurationBr * nearPionPass * farMuonPass );
    currentAcceptance->addToOneJetOnePositiveElectron(
                                 tauToPionTimesTauToElectronBr * currentPass );
    currentAcceptance->addToOneJetOnePositiveMuon( tauToPionTimesTauToMuonBr
                                                   * currentPass );

    currentPass
    = ( ( negativeNearConfigurationBr + positiveNearConfigurationBr )
        * nearMuonPass * farMuonPass );
    currentAcceptance->addToElectronPlusAntielectron( tauPairToElectronPairBr
                                                      * currentPass );
    currentAcceptance->addToMuonPlusAntimuon( tauPairToMuonPairBr
                                              * currentPass );
    currentPass *= tauToElectronTimesTauToMuonBr;
    currentAcceptance->addToNegativeElectronPlusPositiveMuon( currentPass );
    currentAcceptance->addToNegativeMuonPlusPositiveElectron( currentPass );

    currentAcceptance->addToOneJetZeroLeptons( ( negativeNearConfigurationBr
                                                + positiveNearConfigurationBr )
                                               * ( tauPairToPionPairBr
                                                 * ( nearPionPass * farPionFail
                                                 + nearPionFail * farPionPass )
                                              + ( tauToPionTimesTauToElectronBr
                                                  + tauToPionTimesTauToMuonBr )
                                                 * ( nearPionPass * farMuonFail
                                            + nearMuonFail * farPionPass ) ) );

    currentNearFailSum
    = ( tauToPionTimesTauToElectronBr * nearPionFail
        + ( tauPairToElectronPairBr + tauToElectronTimesTauToMuonBr )
          * nearMuonFail );
    currentFarFailSum
    = ( tauToPionTimesTauToElectronBr * farPionFail
        + ( tauPairToElectronPairBr + tauToElectronTimesTauToMuonBr )
          * farMuonFail );
    currentAcceptance->addToZeroJetsOneNegativeElectron(
                                 ( negativeNearConfigurationBr * nearMuonPass )
                                 * currentFarFailSum
                                + ( positiveNearConfigurationBr * farMuonPass )
                                  * currentNearFailSum );
    currentAcceptance->addToZeroJetsOnePositiveElectron(
                                 ( positiveNearConfigurationBr * nearMuonPass )
                                 * currentFarFailSum
                                + ( negativeNearConfigurationBr * farMuonPass )
                                  * currentNearFailSum );

    currentNearFailSum
    = ( tauToPionTimesTauToMuonBr * nearPionFail
        + ( tauToElectronTimesTauToMuonBr + tauPairToMuonPairBr )
        * nearMuonFail );
    currentFarFailSum
    = ( tauToPionTimesTauToMuonBr * farPionFail
        + ( tauToElectronTimesTauToMuonBr + tauPairToMuonPairBr )
          * farMuonFail );
    currentAcceptance->addToZeroJetsOneNegativeMuon(
                                 ( negativeNearConfigurationBr * nearMuonPass )
                                 * currentFarFailSum
                                + ( positiveNearConfigurationBr * farMuonPass )
                                  * currentNearFailSum );
    currentAcceptance->addToZeroJetsOnePositiveMuon(
                                 ( positiveNearConfigurationBr * nearMuonPass )
                                 * currentFarFailSum
                                + ( negativeNearConfigurationBr * farMuonPass )
                                  * currentNearFailSum );

    currentAcceptance->addToZeroJetsZeroLeptons( ( negativeNearConfigurationBr
                                                + positiveNearConfigurationBr )
                                                 * ( tauPairToPionPairBr
                                                   * nearPionFail * farPionFail
                                              + ( tauToPionTimesTauToElectronBr
                                                  + tauToPionTimesTauToMuonBr )
                                                * ( nearPionFail * farMuonFail
                                                 + nearMuonFail * farPionFail )
                                                 + ( tauPairToElectronPairBr
                                          + 2.0 * tauToElectronTimesTauToMuonBr
                                                     + tauPairToMuonPairBr )
                                              * nearMuonFail * farMuonFail ) );
  }



  chargeSummedNeutralinoToStauCascade::chargeSummedNeutralinoToStauCascade(
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
                        shortcut )
  {
    // we set up the BR of the (only-ish) channel:
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
    nearTauDistribution
    = new flat_near_muon_plus_antimuon( shortcut->getReadier(),
                                        shortcut->getCppSlha(),
                                        coloredDecayer,
                                        effectiveSquarkMass,
                                        electroweakDecayer,
                                        intermediateDecayer );
    farTauDistribution
    = new flat_far_muon_plus_antimuon( shortcut->getReadier(),
                                       shortcut->getCppSlha(),
                                       coloredDecayer,
                                       effectiveSquarkMass,
                                       electroweakDecayer,
                                       intermediateDecayer,
                                       shortcut->getNeutralinoOne() );
    activeDistributions.push_back( nearTauDistribution );
    activeDistributions.push_back( farTauDistribution );

    nearHardMuonDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     nearTauDistribution,
                                     shortcut->getHardMuonFromTau() );
    nearSoftMuonDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     nearTauDistribution,
                                     shortcut->getSoftMuonFromTau() );
    nearHardPionDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     nearTauDistribution,
                                     shortcut->getHardPionFromTau() );
    nearSoftPionDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     nearTauDistribution,
                                     shortcut->getSoftPionFromTau() );
    farHardMuonDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     farTauDistribution,
                                     shortcut->getHardMuonFromTau() );
    farSoftMuonDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     farTauDistribution,
                                     shortcut->getSoftMuonFromTau() );
    farHardPionDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     farTauDistribution,
                                     shortcut->getHardPionFromTau() );
    farSoftPionDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     farTauDistribution,
                                     shortcut->getSoftPionFromTau() );
    activeDistributions.push_back( nearHardMuonDistribution );
    activeDistributions.push_back( nearSoftMuonDistribution );
    activeDistributions.push_back( nearHardPionDistribution );
    activeDistributions.push_back( nearSoftPionDistribution );
    activeDistributions.push_back( farHardMuonDistribution );
    activeDistributions.push_back( farSoftMuonDistribution );
    activeDistributions.push_back( farHardPionDistribution );
    activeDistributions.push_back( farSoftPionDistribution );
  }

  chargeSummedNeutralinoToStauCascade::~chargeSummedNeutralinoToStauCascade()
  {
    // does nothing.
  }


  void
  chargeSummedNeutralinoToStauCascade::calculateAcceptance(
                                                  acceptanceCutSet* const cuts,
                                    acceptanceValues* const currentAcceptance )
  // this returns the appropriate acceptances multiplied by branching ratios
  // from the electroweakino through the stau to the LSP.
  {
    currentAcceptance->setTwoJets( 0.0 );
    currentAcceptance->setOneJetOneNegativeElectron( 0.0 );
    currentAcceptance->setOneJetOnePositiveElectron( 0.0 );
    currentAcceptance->setOneJetOneNegativeMuon( 0.0 );
    currentAcceptance->setOneJetOnePositiveMuon( 0.0 );
    currentAcceptance->setElectronPlusAntielectron( 0.0 );
    currentAcceptance->setNegativeElectronPlusPositiveMuon( 0.0 );
    currentAcceptance->setNegativeMuonPlusPositiveElectron( 0.0 );
    currentAcceptance->setMuonPlusAntimuon( 0.0 );
    currentAcceptance->setOneJetZeroLeptons( 0.0 );
    currentAcceptance->setZeroJetsOneNegativeElectron( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveElectron( 0.0 );
    currentAcceptance->setZeroJetsOneNegativeMuon( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveMuon( 0.0 );
    currentAcceptance->setZeroJetsZeroLeptons( 0.0 );
    cascadeBr = ( first_BR->getBr() * second_BR->getBr() );
    if( lhcFaserGlobal::negligibleBr < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      currentCuts = cuts;
      this->currentAcceptance = currentAcceptance;
      jetLeftHandedness
      = shortcut->quarkOrLeptonLeftHandedness(
                                                coloredDecayer->get_PDG_code(),
                                          electroweakDecayer->get_PDG_code() );
      jetRightHandedness = ( 1.0 - jetLeftHandedness );
      nearNegativeTauLeftHandedness
      = shortcut->quarkOrLeptonLeftHandedness(
                                            mediating_particle->get_PDG_code(),
                                              decaying_EWino->get_PDG_code() );
      nearNegativeTauRightHandedness = ( 1.0 - nearNegativeTauLeftHandedness );
      farNegativeTauLeftHandedness
      = shortcut->quarkOrLeptonLeftHandedness(
                                            mediating_particle->get_PDG_code(),
                              shortcut->getNeutralinoOne()->get_PDG_code() );
      farNegativeTauRightHandedness = ( 1.0 - farNegativeTauLeftHandedness );

      /* 1st, left-handed quark, left-handed near negative tau lepton
       * with left-handed far positive tau lepton + right-handed near
       * positive tau lepton with right-handed far negative tau lepton,
       * plus the all-handedness-flipped versions:
       */
      currentNearMuonDistribution = nearHardMuonDistribution;
      currentNearPionDistribution = nearSoftPionDistribution;
      currentFarMuonDistribution = farSoftMuonDistribution;
      currentFarPionDistribution = farHardPionDistribution;
      configurationBr = ( cascadeBr * nearNegativeTauLeftHandedness
                                    * farNegativeTauRightHandedness );
      calculateForCurrentConfiguration();

      /* 2nd, left-handed quark, left-handed near negative tau lepton
       * with right-handed far positive tau lepton + right-handed near
       * positive tau lepton with left-handed far negative tau lepton,
       * plus the all-handedness-flipped versions:
       */
      currentNearMuonDistribution = nearHardMuonDistribution;
      currentNearPionDistribution = nearSoftPionDistribution;
      currentFarMuonDistribution = farHardMuonDistribution;
      currentFarPionDistribution = farSoftPionDistribution;
      configurationBr = ( cascadeBr * nearNegativeTauLeftHandedness
                                    * farNegativeTauLeftHandedness );
      calculateForCurrentConfiguration();

      /* 3rd, left-handed quark, right-handed near negative tau lepton
       * with left-handed far positive tau lepton + left-handed near
       * positive tau lepton with right-handed far negative tau lepton,
       * plus the all-handedness-flipped versions:
       */
      currentNearMuonDistribution = nearSoftMuonDistribution;
      currentNearPionDistribution = nearHardPionDistribution;
      currentFarMuonDistribution = farSoftMuonDistribution;
      currentFarPionDistribution = farHardPionDistribution;
      configurationBr = ( cascadeBr * nearNegativeTauRightHandedness
                                    * farNegativeTauRightHandedness );
      calculateForCurrentConfiguration();

      /* 4th, left-handed quark, right-handed near negative tau lepton
       * with right-handed far positive tau lepton + left-handed near
       * positive tau lepton with left-handed far negative tau lepton,
       * plus the all-handedness-flipped versions:
       */
      currentNearMuonDistribution = nearSoftMuonDistribution;
      currentNearPionDistribution = nearHardPionDistribution;
      currentFarMuonDistribution = farHardMuonDistribution;
      currentFarPionDistribution = farSoftPionDistribution;
      configurationBr = ( cascadeBr * nearNegativeTauRightHandedness
                                    * farNegativeTauLeftHandedness );
      calculateForCurrentConfiguration();
    }
    // debugging:
    /**std::cout
    << std::endl
    << "cascadeBr = " << cascadeBr;
    std::cout << std::endl;**/
  }

  void
  chargeSummedNeutralinoToStauCascade::calculateForCurrentConfiguration()
  {
    nearMuonPass = integrateAcceptance( currentNearMuonDistribution,
                                        currentCuts->getPrimaryLeptonCut() );
    nearMuonFail = ( 1.0 - integrateAcceptance( currentNearMuonDistribution,
                                      currentCuts->getSecondaryLeptonCut() ) );
    nearPionPass = integrateAcceptance( currentNearPionDistribution,
                                        currentCuts->getJetCut() );
    nearPionFail = ( 1.0 - nearPionPass );
    farMuonPass = integrateAcceptance( currentFarMuonDistribution,
                                       currentCuts->getPrimaryLeptonCut() );
    farMuonFail = ( 1.0 - integrateAcceptance( currentFarMuonDistribution,
                                      currentCuts->getSecondaryLeptonCut() ) );
    farPionPass = integrateAcceptance( currentFarPionDistribution,
                                       currentCuts->getJetCut() );
    farPionFail = ( 1.0 - farPionPass );

    currentAcceptance->addToTwoJets( 2.0 * configurationBr
                                         * tauPairToPionPairBr
                                         * nearPionPass * farPionPass );
    // the factor of 2.0 accounts for the charge-conjugate decay.

    currentPass = ( configurationBr * ( nearMuonPass * farPionPass
                                        + nearPionPass * farMuonPass ) );
    currentAcceptance->addToOneJetOneNegativeElectron(
                                 tauToPionTimesTauToElectronBr * currentPass );
    currentAcceptance->addToOneJetOnePositiveElectron(
                                 tauToPionTimesTauToElectronBr * currentPass );
    currentAcceptance->addToOneJetOnePositiveMuon(
                                     tauToPionTimesTauToMuonBr * currentPass );
    currentAcceptance->addToOneJetOneNegativeMuon(
                                     tauToPionTimesTauToMuonBr * currentPass );
    /* the charges are split evenly (chargeSummedNeutralinoToStauCascade should
     * only be used in situations where overall the charges are even, such as
     * in the case of a gluino decaying, which has equal branching ratio into
     * each charge version, so all the charge asymmetry cancels out, so we
     * might as well not bother calculating the asymmetry in the 1st place.
     */

    currentPass = ( 2.0 * configurationBr * nearMuonPass * farMuonPass );
    // the factor of 2.0 accounts for the charge-conjugate decay.
    currentAcceptance->addToElectronPlusAntielectron( tauPairToElectronPairBr
                                                      * currentPass );
    currentAcceptance->addToMuonPlusAntimuon( tauPairToMuonPairBr
                                              * currentPass );
    currentPass *= tauToElectronTimesTauToMuonBr;
    currentAcceptance->addToNegativeElectronPlusPositiveMuon( currentPass );
    currentAcceptance->addToNegativeMuonPlusPositiveElectron( currentPass );
    // as above, the charges are split evenly.

    currentAcceptance->addToOneJetZeroLeptons( 2.0 * configurationBr
                                                   * ( tauPairToPionPairBr
                                                 * ( nearPionPass * farPionFail
                                                 + nearPionFail * farPionPass )
                                              + ( tauToPionTimesTauToElectronBr
                                                  + tauToPionTimesTauToMuonBr )
                                                 * ( nearPionPass * farMuonFail
                                            + nearMuonFail * farPionPass ) ) );
    // the factor of 2.0 accounts for the charge-conjugate decay.

    currentPass
    = ( configurationBr
        * ( nearMuonPass
            * ( tauToPionTimesTauToElectronBr * nearPionFail
                + ( tauPairToElectronPairBr + tauToElectronTimesTauToMuonBr )
                  * nearMuonFail )
            + farMuonPass
              * ( tauToPionTimesTauToElectronBr * farPionFail
                  + ( tauPairToElectronPairBr + tauToElectronTimesTauToMuonBr )
                    * farMuonFail ) ) );
    currentAcceptance->addToZeroJetsOneNegativeElectron( currentPass );
    currentAcceptance->addToZeroJetsOnePositiveElectron( currentPass );
    currentPass
    = ( configurationBr
        * ( nearMuonPass
            * ( tauToPionTimesTauToMuonBr * nearPionFail
                + ( tauToElectronTimesTauToMuonBr + tauPairToMuonPairBr )
                  * nearMuonFail )
            + farMuonPass
              * ( tauToPionTimesTauToMuonBr * farPionFail
                  + ( tauToElectronTimesTauToMuonBr + tauPairToMuonPairBr )
                    * farMuonFail ) ) );
    currentAcceptance->addToZeroJetsOneNegativeMuon( currentPass );
    currentAcceptance->addToZeroJetsOnePositiveMuon( currentPass );

    currentAcceptance->addToZeroJetsZeroLeptons( configurationBr
                                                 * ( tauPairToPionPairBr
                                                   * nearPionFail * farPionFail
                                              + ( tauToPionTimesTauToElectronBr
                                                  + tauToPionTimesTauToMuonBr )
                                                 * ( nearPionFail * farMuonFail
                                                 + nearMuonFail * farPionFail )
                                                 + ( tauPairToElectronPairBr
                                          + 2.0 * tauToElectronTimesTauToMuonBr
                                                 + tauPairToMuonPairBr )
                                              * nearMuonFail * farMuonFail ) );
  }



  neutralinoToZCascade::neutralinoToZCascade(
                                leptonAcceptanceParameterSet* const kinematics,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                CppSLHA::particle_property_set const* const electroweakDecayer,
                                        inputHandler const* const shortcut ) :
    electroweakCascade( kinematics,
                        effectiveSquarkMass,
                        coloredDecayer,
                        electroweakDecayer,
                        shortcut->getZ(),
                        true,
                        shortcut )
  {
    double
    twiceSquareOfWeakCosine
    = ( 2.0 * shortcut->getWeakCosine() * shortcut->getWeakCosine() );
    negativeTauLeftHandedness
    = ( ( 1 - 2.0 * twiceSquareOfWeakCosine
          + twiceSquareOfWeakCosine * twiceSquareOfWeakCosine )
        / ( 5.0 - 2.0 * twiceSquareOfWeakCosine
            + 2.0 * twiceSquareOfWeakCosine * twiceSquareOfWeakCosine ) );
    // here firstBr & secondBr don't work as well as just getting the branching
    // ratios from the CppSLHA::particle_property_set pointers.
    directMuonDistribution = new Z_handed_muon( shortcut->getReadier(),
                                                shortcut->getCppSlha(),
                                                coloredDecayer,
                                                effectiveSquarkMass,
                                                electroweakDecayer,
                                                intermediateDecayer,
                                                shortcut->getNeutralinoOne(),
                                                true,
                                                true );
    sameHandedTauDistribution = new Z_handed_muon( shortcut->getReadier(),
                                                   shortcut->getCppSlha(),
                                                   coloredDecayer,
                                                   effectiveSquarkMass,
                                                   electroweakDecayer,
                                                   intermediateDecayer,
                                                shortcut->getNeutralinoOne(),
                                                   true,
                                                   false );
    oppositeHandedTauDistribution = new Z_handed_muon( shortcut->getReadier(),
                                                       shortcut->getCppSlha(),
                                                       coloredDecayer,
                                                       effectiveSquarkMass,
                                                       electroweakDecayer,
                                                       intermediateDecayer,
                                                shortcut->getNeutralinoOne(),
                                                       false,
                                                       false );
    activeDistributions.push_back( directMuonDistribution );
    activeDistributions.push_back( sameHandedTauDistribution );
    activeDistributions.push_back( oppositeHandedTauDistribution );

    sameHardMuonDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     sameHandedTauDistribution,
                                     shortcut->getHardMuonFromTau() );
    sameSoftMuonDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     sameHandedTauDistribution,
                                     shortcut->getSoftMuonFromTau() );
    sameHardPionDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     sameHandedTauDistribution,
                                     shortcut->getHardPionFromTau() );
    sameSoftPionDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     sameHandedTauDistribution,
                                     shortcut->getSoftPionFromTau() );
    oppositeHardMuonDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     oppositeHandedTauDistribution,
                                     shortcut->getHardMuonFromTau() );
    oppositeSoftMuonDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     oppositeHandedTauDistribution,
                                     shortcut->getSoftMuonFromTau() );
    oppositeHardPionDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     oppositeHandedTauDistribution,
                                     shortcut->getHardPionFromTau() );
    oppositeSoftPionDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     oppositeHandedTauDistribution,
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

  neutralinoToZCascade::~neutralinoToZCascade()
  {
    // does nothing.
  }


  bool
  neutralinoToZCascade::validSignal( int const numberOfJets,
                                     int const numberOfNegativeElectrons,
                                     int const numberOfPositiveElectrons,
                                     int const numberOfNegativeMuons,
                                     int const numberOfPositiveMuons )
  /* this returns true if a configuration where each of the signs of tau lepton
   * decayed either into a detected jet, detected lepton, or undetected
   * particle, & false otherwise.
   */
  {
    if( 2 == numberOfJets )
    {
      if( ( 0 == numberOfNegativeElectrons )
          &&
          ( 0 == numberOfPositiveElectrons )
          &&
          ( 0 == numberOfNegativeMuons )
          &&
          ( 0 == numberOfPositiveMuons ) )
      {
        return true;
      }
      else
      {
        return false;
      }
    }
    else if( ( 0 <= numberOfJets )
             &&
             ( 0 <= numberOfNegativeElectrons )
             &&
             ( 0 <= numberOfPositiveElectrons )
             &&
             ( 0 <= numberOfNegativeMuons )
             &&
             ( 0 <= numberOfPositiveMuons )
             &&
             ( 1 >= ( numberOfNegativeElectrons + numberOfNegativeMuons ) )
             &&
             ( 1 >= ( numberOfPositiveElectrons + numberOfPositiveMuons ) ) )
    {
      if( 2 >=
          ( numberOfJets
            + numberOfNegativeElectrons + numberOfPositiveElectrons
            + numberOfNegativeMuons + numberOfPositiveMuons ) )
      {
        return true;
      }
      else
      {
        return false;
      }
    }
    else
    {
      return false;
    }
  }

  void
  neutralinoToZCascade::calculateAcceptance( acceptanceCutSet* const cuts,
                                    acceptanceValues* const currentAcceptance )
  // this returns the appropriate acceptances multiplied by branching ratios..
  {
    currentAcceptance->setTwoJets( 0.0 );
    currentAcceptance->setOneJetOneNegativeElectron( 0.0 );
    currentAcceptance->setOneJetOnePositiveElectron( 0.0 );
    currentAcceptance->setOneJetOneNegativeMuon( 0.0 );
    currentAcceptance->setOneJetOnePositiveMuon( 0.0 );
    currentAcceptance->setElectronPlusAntielectron( 0.0 );
    currentAcceptance->setNegativeElectronPlusPositiveMuon( 0.0 );
    currentAcceptance->setNegativeMuonPlusPositiveElectron( 0.0 );
    currentAcceptance->setMuonPlusAntimuon( 0.0 );
    currentAcceptance->setOneJetZeroLeptons( 0.0 );
    currentAcceptance->setZeroJetsOneNegativeElectron( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveElectron( 0.0 );
    currentAcceptance->setZeroJetsOneNegativeMuon( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveMuon( 0.0 );
    currentAcceptance->setZeroJetsZeroLeptons( 0.0 );
    cascadeBr
    = ( electroweakDecayer->inspect_direct_decay_handler(
                  )->get_branching_ratio_for_exact_match( CppSLHA::PDG_code::Z,
                                         CppSLHA::PDG_code::neutralino_one ) );
    // the branching ratios of the Z into the various SM fermions are covered
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

      directMuonPass = integrateAcceptance( directMuonDistribution,
                                            cuts->getPrimaryLeptonCut() );
      directMuonFail = ( 1.0 - integrateAcceptance( directMuonDistribution,
                                             cuts->getSecondaryLeptonCut() ) );
      directJetPass = integrateAcceptance( directMuonDistribution,
                                           cuts->getJetCut() );
      directJetFail = ( 1.0 - directJetPass );

      configurationBr
      = ( cascadeBr * CppSLHA::PDG_data::Z_to_hadrons_BR * directJetPass );
      currentAcceptance->addToTwoJets( configurationBr * directJetPass );
      currentAcceptance->addToOneJetZeroLeptons( 2.0 * configurationBr
                                                     * directJetFail );
      currentAcceptance->addToOssfMinusOsdf( cascadeBr
                           * ( CppSLHA::PDG_data::Z_to_electron_antielectron_BR
                               + CppSLHA::PDG_data::Z_to_muon_antimuon_BR )
                             * directMuonPass * directMuonPass );
      configurationBr
      = ( cascadeBr * CppSLHA::PDG_data::Z_to_electron_antielectron_BR
          * directMuonPass );
      currentAcceptance->addToElectronPlusAntielectron( configurationBr
                                                        * directMuonPass );
      configurationBr *= directMuonFail;
      currentAcceptance->addToZeroJetsOneNegativeElectron( configurationBr );
      currentAcceptance->addToZeroJetsOnePositiveElectron( configurationBr );
      configurationBr
      = ( cascadeBr * CppSLHA::PDG_data::Z_to_muon_antimuon_BR
          * directMuonPass );
      currentAcceptance->addToMuonPlusAntimuon( configurationBr
                                                * directMuonPass );
      configurationBr *= directMuonFail;
      currentAcceptance->addToZeroJetsOneNegativeMuon( configurationBr );
      currentAcceptance->addToZeroJetsOnePositiveMuon( configurationBr );
      currentAcceptance->addToZeroJetsZeroLeptons( cascadeBr
                                         * ( CppSLHA::PDG_data::Z_to_hadrons_BR
                                             * directJetFail * directJetFail
                           + ( CppSLHA::PDG_data::Z_to_electron_antielectron_BR
                                   + CppSLHA::PDG_data::Z_to_muon_antimuon_BR )
                             * directMuonFail * directMuonFail ) );

      // now Z decays to tau-antitau pairs, followed by the decays of the taus:
      cascadeBr *= CppSLHA::PDG_data::Z_to_tau_lepton_antilepton_BR;
      jetLeftHandedness
      = shortcut->quarkOrLeptonLeftHandedness(
                                                coloredDecayer->get_PDG_code(),
                                          electroweakDecayer->get_PDG_code() );
      // left-handed quark, left-handed negative tau lepton:
      configurationBr
      = ( cascadeBr * jetLeftHandedness * negativeTauLeftHandedness );
      currentNegativeMuonDistribution = sameHardMuonDistribution;
      currentNegativePionDistribution = sameSoftPionDistribution;
      currentPositiveMuonDistribution = oppositeHardMuonDistribution;
      currentPositivePionDistribution = oppositeSoftPionDistribution;
      calculateForCurrentConfiguration();

      // left-handed quark, right-handed negative tau lepton:
      configurationBr
      = ( cascadeBr * jetLeftHandedness
                    * ( 1.0 - negativeTauLeftHandedness ) );
      currentNegativeMuonDistribution = oppositeSoftMuonDistribution;
      currentNegativePionDistribution = oppositeHardPionDistribution;
      currentPositiveMuonDistribution = sameSoftMuonDistribution;
      currentPositivePionDistribution = sameHardPionDistribution;
      calculateForCurrentConfiguration();

      // right-handed quark, left-handed negative tau lepton:
      configurationBr
      = ( cascadeBr * ( 1.0 - jetLeftHandedness )
          * negativeTauLeftHandedness );
      currentNegativeMuonDistribution = oppositeHardMuonDistribution;
      currentNegativePionDistribution = oppositeSoftPionDistribution;
      currentPositiveMuonDistribution = sameHardMuonDistribution;
      currentPositivePionDistribution = sameSoftPionDistribution;
      calculateForCurrentConfiguration();

      // right-handed quark, right-handed negative tau lepton:
      configurationBr
      = ( cascadeBr * ( 1.0 - jetLeftHandedness )
                    * ( 1.0 - negativeTauLeftHandedness ) );
      currentNegativeMuonDistribution = sameSoftMuonDistribution;
      currentNegativePionDistribution = sameHardPionDistribution;
      currentPositiveMuonDistribution = oppositeSoftMuonDistribution;
      currentPositivePionDistribution = oppositeHardPionDistribution;
      calculateForCurrentConfiguration();
    }
    // debugging:
    /**std::cout
    << std::endl
    << "cascadeBr = " << cascadeBr;
    std::cout << std::endl;**/
  }

  void
  neutralinoToZCascade::calculateForCurrentConfiguration()
  {
    negativeTauMuonPass = integrateAcceptance( currentNegativeMuonDistribution,
                                          currentCuts->getPrimaryLeptonCut() );
    negativeTauMuonFail
    = ( 1.0 - integrateAcceptance( currentNegativeMuonDistribution,
                                   currentCuts->getSecondaryLeptonCut() ) );
    negativeTauPionPass = integrateAcceptance( currentNegativePionDistribution,
                                               currentCuts->getJetCut() );
    negativeTauPionFail = ( 1.0 - negativeTauPionPass );
    positiveTauMuonPass = integrateAcceptance( currentPositiveMuonDistribution,
                                          currentCuts->getPrimaryLeptonCut() );
    positiveTauMuonFail
    = ( 1.0 - integrateAcceptance( currentPositiveMuonDistribution,
                                   currentCuts->getSecondaryLeptonCut() ) );
    positiveTauPionPass = integrateAcceptance( currentPositivePionDistribution,
                                               currentCuts->getJetCut() );
    positiveTauPionFail = ( 1.0 - positiveTauPionPass );

    currentAcceptance->addToTwoJets( configurationBr
                                     * tauPairToPionPairBr
                                     * negativeTauPionPass
                                     * positiveTauPionPass );
    currentPass = ( configurationBr * tauToPionTimesTauToElectronBr );
    currentAcceptance->addToOneJetOneNegativeElectron( currentPass
                                                       * positiveTauPionPass
                                                       * negativeTauMuonPass );
    currentAcceptance->addToOneJetOnePositiveElectron( currentPass
                                                       * negativeTauPionPass
                                                       * positiveTauMuonPass );
    currentPass = ( configurationBr * tauToPionTimesTauToMuonBr );
    currentAcceptance->addToOneJetOneNegativeMuon( currentPass
                                                   * positiveTauPionPass
                                                   * negativeTauMuonPass );
    currentAcceptance->addToOneJetOnePositiveMuon( currentPass
                                                   * negativeTauPionPass
                                                   * positiveTauMuonPass );
    currentPass
    = ( configurationBr * negativeTauMuonPass * positiveTauMuonPass );
    currentAcceptance->addToElectronPlusAntielectron( tauPairToElectronPairBr
                                                      * currentPass );
    currentAcceptance->addToMuonPlusAntimuon( tauPairToMuonPairBr
                                              * currentPass );
    currentPass *= tauToElectronTimesTauToMuonBr;
    currentAcceptance->addToNegativeElectronPlusPositiveMuon( currentPass );
    currentAcceptance->addToNegativeMuonPlusPositiveElectron( currentPass );

    currentAcceptance->addToOneJetZeroLeptons( configurationBr
                                               * ( tauPairToPionPairBr
                                                   * ( negativeTauPionPass
                                                       * positiveTauPionFail
                                                       + positiveTauPionPass
                                                        * negativeTauPionFail )
                                              + ( tauToPionTimesTauToElectronBr
                                                  + tauToPionTimesTauToMuonBr )
                                                * ( negativeTauPionPass
                                                    * positiveTauMuonFail
                                                    + positiveTauPionPass
                                                   * negativeTauMuonFail ) ) );

    currentAcceptance->addToZeroJetsOneNegativeElectron( configurationBr
                                                         * negativeTauMuonPass
                                              * ( tauToPionTimesTauToElectronBr
                                                  * positiveTauPionFail
                                                  + ( tauPairToElectronPairBr
                                              + tauToElectronTimesTauToMuonBr )
                                                     * positiveTauMuonFail ) );
    currentAcceptance->addToZeroJetsOnePositiveElectron( configurationBr
                                                         * positiveTauMuonPass
                                              * ( tauToPionTimesTauToElectronBr
                                                  * negativeTauPionFail
                                                  + ( tauPairToElectronPairBr
                                              + tauToElectronTimesTauToMuonBr )
                                                     * negativeTauMuonFail ) );
    currentAcceptance->addToZeroJetsOneNegativeMuon( configurationBr
                                                     * negativeTauMuonPass
                                                  * ( tauToPionTimesTauToMuonBr
                                                      * positiveTauPionFail
                                              + ( tauToElectronTimesTauToMuonBr
                                                  + tauPairToMuonPairBr )
                                                     * positiveTauMuonFail ) );
    currentAcceptance->addToZeroJetsOnePositiveMuon( configurationBr
                                                     * positiveTauMuonPass
                                                  * ( tauToPionTimesTauToMuonBr
                                                      * negativeTauPionFail
                                              + ( tauToElectronTimesTauToMuonBr
                                                  + tauPairToMuonPairBr )
                                                     * negativeTauMuonFail ) );

    currentAcceptance->addToZeroJetsZeroLeptons( configurationBr
                                                 * ( tauPairToPionPairBr
                                                     * negativeTauPionFail
                                                     * positiveTauPionFail
                                              + ( tauToPionTimesTauToElectronBr
                                                  + tauToPionTimesTauToMuonBr )
                                                  * ( negativeTauPionFail
                                                      * positiveTauMuonFail
                                                        + positiveTauPionFail
                                                        * negativeTauMuonFail )
                                                 + ( tauPairToElectronPairBr
                                          + 2.0 * tauToElectronTimesTauToMuonBr
                                                     + tauPairToMuonPairBr )
                                                   * negativeTauMuonFail
                                                   * negativeTauMuonFail ) );
  }



  chargeSummedNeutralinoToZCascade::chargeSummedNeutralinoToZCascade(
                                leptonAcceptanceParameterSet* const kinematics,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                CppSLHA::particle_property_set const* const electroweakDecayer,
                                        inputHandler const* const shortcut ) :
    electroweakCascade( kinematics,
                        effectiveSquarkMass,
                        coloredDecayer,
                        electroweakDecayer,
                        shortcut->getZ(),
                        true,
                        shortcut )
  {
    double
    twiceSquareOfWeakCosine
    = ( 2.0 * shortcut->getWeakCosine() * shortcut->getWeakCosine() );
    negativeTauLeftHandedness
    = ( ( 1 - 2.0 * twiceSquareOfWeakCosine
          + twiceSquareOfWeakCosine * twiceSquareOfWeakCosine )
        / ( 5.0 - 2.0 * twiceSquareOfWeakCosine
            + 2.0 * twiceSquareOfWeakCosine * twiceSquareOfWeakCosine ) );
    // here firstBr & secondBr don't work as well as just getting the branching
    // ratios from the CppSLHA::particle_property_set pointers.
    directMuonDistribution = new Z_handed_muon( shortcut->getReadier(),
                                                shortcut->getCppSlha(),
                                                coloredDecayer,
                                                effectiveSquarkMass,
                                                electroweakDecayer,
                                                intermediateDecayer,
                                                shortcut->getNeutralinoOne(),
                                                true,
                                                true );
    activeDistributions.push_back( directMuonDistribution );

    hardMuonDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     directMuonDistribution,
                                     shortcut->getHardMuonFromTau() );
    softMuonDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     directMuonDistribution,
                                     shortcut->getSoftMuonFromTau() );
    hardPionDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     directMuonDistribution,
                                     shortcut->getHardPionFromTau() );
    softPionDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     directMuonDistribution,
                                     shortcut->getSoftPionFromTau() );
    activeDistributions.push_back( hardMuonDistribution );
    activeDistributions.push_back( softMuonDistribution );
    activeDistributions.push_back( hardPionDistribution );
    activeDistributions.push_back( softPionDistribution );
  }

  chargeSummedNeutralinoToZCascade::~chargeSummedNeutralinoToZCascade()
  {
    // does nothing.
  }


  bool
  chargeSummedNeutralinoToZCascade::validSignal( int const numberOfJets,
                                           int const numberOfNegativeElectrons,
                                           int const numberOfPositiveElectrons,
                                               int const numberOfNegativeMuons,
                                              int const numberOfPositiveMuons )
  /* this returns true if a configuration where each of the signs of tau lepton
   * decayed either into a detected jet, detected lepton, or undetected
   * particle, & false otherwise.
   */
  {
    if( 2 == numberOfJets )
    {
      if( ( 0 == numberOfNegativeElectrons )
          &&
          ( 0 == numberOfPositiveElectrons )
          &&
          ( 0 == numberOfNegativeMuons )
          &&
          ( 0 == numberOfPositiveMuons ) )
      {
        return true;
      }
      else
      {
        return false;
      }
    }
    else if( ( 0 <= numberOfJets )
             &&
             ( 0 <= numberOfNegativeElectrons )
             &&
             ( 0 <= numberOfPositiveElectrons )
             &&
             ( 0 <= numberOfNegativeMuons )
             &&
             ( 0 <= numberOfPositiveMuons )
             &&
             ( 1 >= ( numberOfNegativeElectrons + numberOfNegativeMuons ) )
             &&
             ( 1 >= ( numberOfPositiveElectrons + numberOfPositiveMuons ) ) )
    {
      if( 2 >=
          ( numberOfJets
            + numberOfNegativeElectrons + numberOfPositiveElectrons
            + numberOfNegativeMuons + numberOfPositiveMuons ) )
      {
        return true;
      }
      else
      {
        return false;
      }
    }
    else
    {
      return false;
    }
  }

  void
  chargeSummedNeutralinoToZCascade::calculateAcceptance(
                                                  acceptanceCutSet* const cuts,
                                    acceptanceValues* const currentAcceptance )
  // this returns the appropriate acceptances multiplied by branching ratios.
  {
    currentAcceptance->setTwoJets( 0.0 );
    currentAcceptance->setOneJetOneNegativeElectron( 0.0 );
    currentAcceptance->setOneJetOnePositiveElectron( 0.0 );
    currentAcceptance->setOneJetOneNegativeMuon( 0.0 );
    currentAcceptance->setOneJetOnePositiveMuon( 0.0 );
    currentAcceptance->setElectronPlusAntielectron( 0.0 );
    currentAcceptance->setNegativeElectronPlusPositiveMuon( 0.0 );
    currentAcceptance->setNegativeMuonPlusPositiveElectron( 0.0 );
    currentAcceptance->setMuonPlusAntimuon( 0.0 );
    currentAcceptance->setOneJetZeroLeptons( 0.0 );
    currentAcceptance->setZeroJetsOneNegativeElectron( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveElectron( 0.0 );
    currentAcceptance->setZeroJetsOneNegativeMuon( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveMuon( 0.0 );
    currentAcceptance->setZeroJetsZeroLeptons( 0.0 );
    cascadeBr
    = ( electroweakDecayer->inspect_direct_decay_handler(
                  )->get_branching_ratio_for_exact_match( CppSLHA::PDG_code::Z,
                                         CppSLHA::PDG_code::neutralino_one ) );
    // the branching ratios of the Z into the various SM fermions are covered
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

      directMuonPass = integrateAcceptance( directMuonDistribution,
                                            cuts->getPrimaryLeptonCut() );
      directMuonFail = ( 1.0 - integrateAcceptance( directMuonDistribution,
                                             cuts->getSecondaryLeptonCut() ) );
      directJetPass = integrateAcceptance( directMuonDistribution,
                                           cuts->getJetCut() );
      directJetFail = ( 1.0 - directJetPass );

      configurationBr
      = ( cascadeBr * CppSLHA::PDG_data::Z_to_hadrons_BR * directJetPass );
      currentAcceptance->addToTwoJets( configurationBr * directJetPass );
      currentAcceptance->addToOneJetZeroLeptons( 2.0 * configurationBr
                                                     * directJetFail );
      currentAcceptance->addToOssfMinusOsdf( cascadeBr
                           * ( CppSLHA::PDG_data::Z_to_electron_antielectron_BR
                               + CppSLHA::PDG_data::Z_to_muon_antimuon_BR )
                             * directMuonPass * directMuonPass );
      configurationBr
      = ( cascadeBr * CppSLHA::PDG_data::Z_to_electron_antielectron_BR
          * directMuonPass );
      currentAcceptance->addToElectronPlusAntielectron( configurationBr
                                                        * directMuonPass );
      configurationBr *= directMuonFail;
      currentAcceptance->addToZeroJetsOneNegativeElectron( configurationBr );
      currentAcceptance->addToZeroJetsOnePositiveElectron( configurationBr );
      configurationBr
      = ( cascadeBr * CppSLHA::PDG_data::Z_to_muon_antimuon_BR
          * directMuonPass );
      currentAcceptance->addToMuonPlusAntimuon( configurationBr
                                                * directMuonPass );
      configurationBr *= directMuonFail;
      currentAcceptance->addToZeroJetsOneNegativeMuon( configurationBr );
      currentAcceptance->addToZeroJetsOnePositiveMuon( configurationBr );
      currentAcceptance->addToZeroJetsZeroLeptons( cascadeBr
                                         * ( CppSLHA::PDG_data::Z_to_hadrons_BR
                                             * directJetFail * directJetFail
                           + ( CppSLHA::PDG_data::Z_to_electron_antielectron_BR
                                   + CppSLHA::PDG_data::Z_to_muon_antimuon_BR )
                             * directMuonFail * directMuonFail ) );

      // now Z decays to tau-antitau pairs, followed by the decays of the taus:
      cascadeBr *= CppSLHA::PDG_data::Z_to_tau_lepton_antilepton_BR;
      // left-handed negative tau lepton:
      configurationBr
      = ( cascadeBr * negativeTauLeftHandedness );
      currentMuonDistribution = hardMuonDistribution;
      currentPionDistribution = softPionDistribution;
      calculateForCurrentConfiguration();

      // right-handed negative tau lepton:
      configurationBr
      = ( cascadeBr * ( 1.0 - negativeTauLeftHandedness ) );
      currentMuonDistribution = softMuonDistribution;
      currentPionDistribution = hardPionDistribution;
      calculateForCurrentConfiguration();
    }
    // debugging:
    /**std::cout
    << std::endl
    << "cascadeBr = " << cascadeBr;
    std::cout << std::endl;**/
  }

  void
  chargeSummedNeutralinoToZCascade::calculateForCurrentConfiguration()
  {
    tauMuonPass = integrateAcceptance( currentMuonDistribution,
                                       currentCuts->getPrimaryLeptonCut() );
    tauMuonFail = ( 1.0 - integrateAcceptance( currentMuonDistribution,
                                      currentCuts->getSecondaryLeptonCut() ) );
    tauPionPass = integrateAcceptance( currentPionDistribution,
                                       currentCuts->getJetCut() );
    tauPionFail = ( 1.0 - tauPionPass );

    currentAcceptance->addToTwoJets( configurationBr
                                     * tauPairToPionPairBr
                                     * tauPionPass * tauPionPass );
    currentPass = ( configurationBr * tauToPionTimesTauToElectronBr
                    * tauPionPass * tauMuonPass );
    currentAcceptance->addToOneJetOneNegativeElectron( currentPass );
    currentAcceptance->addToOneJetOnePositiveElectron( currentPass );
    currentPass = ( configurationBr * tauToPionTimesTauToMuonBr
                    * tauPionPass * tauMuonPass );
    currentAcceptance->addToOneJetOneNegativeMuon( currentPass );
    currentAcceptance->addToOneJetOnePositiveMuon( currentPass );
    currentPass = ( configurationBr * tauMuonPass * tauMuonPass );
    currentAcceptance->addToElectronPlusAntielectron( tauPairToElectronPairBr
                                                      * currentPass );
    currentAcceptance->addToMuonPlusAntimuon( tauPairToMuonPairBr
                                              * currentPass );
    currentPass *= tauToElectronTimesTauToMuonBr;
    currentAcceptance->addToNegativeElectronPlusPositiveMuon( currentPass );
    currentAcceptance->addToNegativeMuonPlusPositiveElectron( currentPass );

    currentAcceptance->addToOneJetZeroLeptons( 2.0 * configurationBr
                                                   * tauPionPass
                                                   * ( tauPairToPionPairBr
                                                       * tauPionFail
                                              + ( tauToPionTimesTauToElectronBr
                                                  + tauToPionTimesTauToMuonBr )
                                                * tauMuonFail ) );
    currentPass = ( configurationBr * tauMuonPass
                    * ( tauToPionTimesTauToElectronBr * tauPionFail
                  + ( tauPairToElectronPairBr + tauToElectronTimesTauToMuonBr )
                    * tauMuonFail ) );
    currentAcceptance->addToZeroJetsOneNegativeElectron( currentPass );
    currentAcceptance->addToZeroJetsOnePositiveElectron( currentPass );
    currentPass = ( configurationBr * tauMuonPass
                    * ( tauToPionTimesTauToMuonBr * tauPionFail
                      + ( tauToElectronTimesTauToMuonBr + tauPairToMuonPairBr )
                        * tauMuonFail ) );
    currentAcceptance->addToZeroJetsOneNegativeMuon( currentPass );
    currentAcceptance->addToZeroJetsOnePositiveMuon( currentPass );

    currentAcceptance->addToZeroJetsZeroLeptons( configurationBr
                                                 * ( tauPairToPionPairBr
                                                   * nearPionFail * farPionFail
                                        + 2.0 * ( tauToPionTimesTauToElectronBr
                                                  + tauToPionTimesTauToMuonBr )
                                                  * tauPionFail * tauMuonFail
                                                 + ( tauPairToElectronPairBr
                                          + 2.0 * tauToElectronTimesTauToMuonBr
                                                     + tauPairToMuonPairBr )
                                              * nearMuonFail * farMuonFail ) );
  }



  neutralinoToHiggsCascade::neutralinoToHiggsCascade(
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
                        true,
                        shortcut ),
    intermediateDecayerDecays(
                          intermediateDecayer->inspect_direct_decay_handler() )
  {
    negativeTauLeftHandedness = 0.5;
    // here firstBr & secondBr don't work as well as just getting the branching
    // ratios from the CppSLHA::particle_property_set pointers.
    directMuonDistribution
    = new Higgs_muon_plus_antimuon( shortcut->getReadier(),
                                    shortcut->getCppSlha(),
                                    coloredDecayer,
                                    effectiveSquarkMass,
                                    electroweakDecayer,
                                    intermediateDecayer,
                                    shortcut->getNeutralinoOne() );
    activeDistributions.push_back( directMuonDistribution );

    hardMuonDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     directMuonDistribution,
                                     shortcut->getHardMuonFromTau() );
    softMuonDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     directMuonDistribution,
                                     shortcut->getSoftMuonFromTau() );
    hardPionDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     directMuonDistribution,
                                     shortcut->getHardPionFromTau() );
    softPionDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     directMuonDistribution,
                                     shortcut->getSoftPionFromTau() );
    activeDistributions.push_back( hardMuonDistribution );
    activeDistributions.push_back( softMuonDistribution );
    activeDistributions.push_back( hardPionDistribution );
    activeDistributions.push_back( softPionDistribution );
  }

  neutralinoToHiggsCascade::~neutralinoToHiggsCascade()
  {
    // does nothing.
  }


  bool
  neutralinoToHiggsCascade::validSignal( int const numberOfJets,
                                         int const numberOfNegativeElectrons,
                                         int const numberOfPositiveElectrons,
                                         int const numberOfNegativeMuons,
                                         int const numberOfPositiveMuons )
  /* this returns true if a configuration where each of the signs of tau lepton
   * decayed either into a detected jet, detected lepton, or undetected
   * particle, & false otherwise.
   */
  {
    if( 2 == numberOfJets )
    {
      if( ( 0 == numberOfNegativeElectrons )
          &&
          ( 0 == numberOfPositiveElectrons )
          &&
          ( 0 == numberOfNegativeMuons )
          &&
          ( 0 == numberOfPositiveMuons ) )
      {
        return true;
      }
      else
      {
        return false;
      }
    }
    else if( ( 0 <= numberOfJets )
             &&
             ( 0 <= numberOfNegativeElectrons )
             &&
             ( 0 <= numberOfPositiveElectrons )
             &&
             ( 0 <= numberOfNegativeMuons )
             &&
             ( 0 <= numberOfPositiveMuons )
             &&
             ( 1 >= ( numberOfNegativeElectrons + numberOfNegativeMuons ) )
             &&
             ( 1 >= ( numberOfPositiveElectrons + numberOfPositiveMuons ) ) )
    {
      if( 2 >=
          ( numberOfJets
            + numberOfNegativeElectrons + numberOfPositiveElectrons
            + numberOfNegativeMuons + numberOfPositiveMuons ) )
      {
        return true;
      }
      else
      {
        return false;
      }
    }
    else
    {
      return false;
    }
  }

  void
  neutralinoToHiggsCascade::calculateAcceptance( acceptanceCutSet* const cuts,
                                    acceptanceValues* const currentAcceptance )
  // this returns the appropriate acceptances multiplied by branching ratios.
  {
    currentAcceptance->setTwoJets( 0.0 );
    currentAcceptance->setOneJetOneNegativeElectron( 0.0 );
    currentAcceptance->setOneJetOnePositiveElectron( 0.0 );
    currentAcceptance->setOneJetOneNegativeMuon( 0.0 );
    currentAcceptance->setOneJetOnePositiveMuon( 0.0 );
    currentAcceptance->setElectronPlusAntielectron( 0.0 );
    currentAcceptance->setNegativeElectronPlusPositiveMuon( 0.0 );
    currentAcceptance->setNegativeMuonPlusPositiveElectron( 0.0 );
    currentAcceptance->setMuonPlusAntimuon( 0.0 );
    currentAcceptance->setOneJetZeroLeptons( 0.0 );
    currentAcceptance->setZeroJetsOneNegativeElectron( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveElectron( 0.0 );
    currentAcceptance->setZeroJetsOneNegativeMuon( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveMuon( 0.0 );
    currentAcceptance->setZeroJetsZeroLeptons( 0.0 );
    cascadeBr
    = ( electroweakDecayer->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                           intermediateDecayer->get_PDG_code(),
                                         CppSLHA::PDG_code::neutralino_one ) );
    // the branching ratios of the scalar/pseudoscalar into the various SM
    // fermions are covered by the following code.
    if( lhcFaserGlobal::negligibleBr < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      currentBrToHadrons
      = ( intermediateDecayerDecays->get_branching_ratio_for_exact_match(
                                                       CppSLHA::PDG_code::down,
                                                  -(CppSLHA::PDG_code::down ) )
          + intermediateDecayerDecays->get_branching_ratio_for_exact_match(
                                                         CppSLHA::PDG_code::up,
                                                    -(CppSLHA::PDG_code::up ) )
          + intermediateDecayerDecays->get_branching_ratio_for_exact_match(
                                                    CppSLHA::PDG_code::strange,
                                               -(CppSLHA::PDG_code::strange ) )
          + intermediateDecayerDecays->get_branching_ratio_for_exact_match(
                                                      CppSLHA::PDG_code::charm,
                                                 -(CppSLHA::PDG_code::charm ) )
          + intermediateDecayerDecays->get_branching_ratio_for_exact_match(
                                                     CppSLHA::PDG_code::bottom,
                                             -(CppSLHA::PDG_code::bottom ) ) );
      currentBrToElectrons
      = intermediateDecayerDecays->get_branching_ratio_for_exact_match(
                                                   CppSLHA::PDG_code::electron,
                                             -(CppSLHA::PDG_code::electron ) );
      currentBrToMuons
      = intermediateDecayerDecays->get_branching_ratio_for_exact_match(
                                                       CppSLHA::PDG_code::muon,
                                                 -(CppSLHA::PDG_code::muon ) );
      currentBrToTaus
      = intermediateDecayerDecays->get_branching_ratio_for_exact_match(
                                                 CppSLHA::PDG_code::tau_lepton,
                                           -(CppSLHA::PDG_code::tau_lepton ) );

      directMuonPass = integrateAcceptance( directMuonDistribution,
                                            cuts->getPrimaryLeptonCut() );
      directMuonFail = ( 1.0 - integrateAcceptance( directMuonDistribution,
                                             cuts->getSecondaryLeptonCut() ) );
      directJetPass = integrateAcceptance( directMuonDistribution,
                                           cuts->getJetCut() );
      directJetFail = ( 1.0 - directJetPass );

      configurationBr
      = ( cascadeBr * currentBrToHadrons * directJetPass );
      currentAcceptance->addToTwoJets( configurationBr * directJetPass );
      currentAcceptance->addToOneJetZeroLeptons( 2.0 * configurationBr
                                                     * directJetFail );
      currentAcceptance->addToOssfMinusOsdf( cascadeBr
                                             * ( currentBrToElectrons
                                                 + currentBrToMuons )
                                             * directMuonPass
                                             * directMuonPass );
      configurationBr
      = ( cascadeBr * currentBrToElectrons * directMuonPass );
      currentAcceptance->addToElectronPlusAntielectron( configurationBr
                                                        * directMuonPass );
      configurationBr *= directMuonFail;
      currentAcceptance->addToZeroJetsOneNegativeElectron( configurationBr );
      currentAcceptance->addToZeroJetsOnePositiveElectron( configurationBr );
      configurationBr
      = ( cascadeBr * currentBrToMuons * directMuonPass );
      currentAcceptance->addToMuonPlusAntimuon( configurationBr
                                                * directMuonPass );
      configurationBr *= directMuonFail;
      currentAcceptance->addToZeroJetsOneNegativeMuon( configurationBr );
      currentAcceptance->addToZeroJetsOnePositiveMuon( configurationBr );
      currentAcceptance->addToZeroJetsZeroLeptons( cascadeBr
                                                   * ( currentBrToHadrons
                                                       * directJetFail
                                                       * directJetFail
                                                       + ( currentBrToElectrons
                                                           + currentBrToMuons )
                                                         * directMuonFail
                                                         * directMuonFail ) );

      // now the decays to tau-antitau pairs, followed by the decays of the
      // taus:
      cascadeBr *= currentBrToTaus;
      // left-handed negative tau lepton:
      configurationBr
      = ( cascadeBr * negativeTauLeftHandedness );
      currentNegativeMuonDistribution = hardMuonDistribution;
      currentNegativePionDistribution = softPionDistribution;
      currentPositiveMuonDistribution = softMuonDistribution;
      currentPositivePionDistribution = hardPionDistribution;
      calculateForCurrentConfiguration();

      // right-handed negative tau lepton:
      //configurationBr
      //= ( cascadeBr * ( 1.0 - negativeTauLeftHandedness ) );
      // equal numbers of left- & right-handed tau leptons are produced.
      currentNegativeMuonDistribution = softMuonDistribution;
      currentNegativePionDistribution = hardPionDistribution;
      currentPositiveMuonDistribution = hardMuonDistribution;
      currentPositivePionDistribution = softPionDistribution;
      calculateForCurrentConfiguration();
    }
    // debugging:
    /**std::cout
    << std::endl
    << "cascadeBr = " << cascadeBr;
    std::cout << std::endl;**/
  }

  void
  neutralinoToHiggsCascade::calculateForCurrentConfiguration()
  {
    negativeTauMuonPass = integrateAcceptance( currentNegativeMuonDistribution,
                                          currentCuts->getPrimaryLeptonCut() );
    negativeTauMuonFail
    = ( 1.0 - integrateAcceptance( currentNegativeMuonDistribution,
                                   currentCuts->getSecondaryLeptonCut() ) );
    negativeTauPionPass = integrateAcceptance( currentNegativePionDistribution,
                                               currentCuts->getJetCut() );
    negativeTauPionFail = ( 1.0 - negativeTauPionPass );
    positiveTauMuonPass = integrateAcceptance( currentPositiveMuonDistribution,
                                          currentCuts->getPrimaryLeptonCut() );
    positiveTauMuonFail
    = ( 1.0 - integrateAcceptance( currentPositiveMuonDistribution,
                                   currentCuts->getSecondaryLeptonCut() ) );
    positiveTauPionPass = integrateAcceptance( currentPositivePionDistribution,
                                               currentCuts->getJetCut() );
    positiveTauPionFail = ( 1.0 - positiveTauPionPass );

    currentAcceptance->addToTwoJets( configurationBr
                                     * tauPairToPionPairBr
                                     * negativeTauPionPass
                                     * positiveTauPionPass );
    currentPass = ( configurationBr * tauToPionTimesTauToElectronBr );
    currentAcceptance->addToOneJetOneNegativeElectron( currentPass
                                                       * positiveTauPionPass
                                                       * negativeTauMuonPass );
    currentAcceptance->addToOneJetOnePositiveElectron( currentPass
                                                       * negativeTauPionPass
                                                       * positiveTauMuonPass );
    currentPass = ( configurationBr * tauToPionTimesTauToMuonBr );
    currentAcceptance->addToOneJetOneNegativeMuon( currentPass
                                                   * positiveTauPionPass
                                                   * negativeTauMuonPass );
    currentAcceptance->addToOneJetOnePositiveMuon( currentPass
                                                   * negativeTauPionPass
                                                   * positiveTauMuonPass );
    currentPass
    = ( configurationBr * negativeTauMuonPass * positiveTauMuonPass );
    currentAcceptance->addToElectronPlusAntielectron( tauPairToElectronPairBr
                                                      * currentPass );
    currentAcceptance->addToMuonPlusAntimuon( tauPairToMuonPairBr
                                              * currentPass );
    currentPass *= tauToElectronTimesTauToMuonBr;
    currentAcceptance->addToNegativeElectronPlusPositiveMuon( currentPass );
    currentAcceptance->addToNegativeMuonPlusPositiveElectron( currentPass );

    currentAcceptance->addToOneJetZeroLeptons( configurationBr
                                               * ( tauPairToPionPairBr
                                                   * ( negativeTauPionPass
                                                       * positiveTauPionFail
                                                       + positiveTauPionPass
                                                        * negativeTauPionFail )
                                              + ( tauToPionTimesTauToElectronBr
                                                  + tauToPionTimesTauToMuonBr )
                                                * ( negativeTauPionPass
                                                    * positiveTauMuonFail
                                                    + positiveTauPionPass
                                                   * negativeTauMuonFail ) ) );

    currentAcceptance->addToZeroJetsOneNegativeElectron( configurationBr
                                                         * negativeTauMuonPass
                                              * ( tauToPionTimesTauToElectronBr
                                                  * positiveTauPionFail
                                                  + ( tauPairToElectronPairBr
                                              + tauToElectronTimesTauToMuonBr )
                                                     * positiveTauMuonFail ) );
    currentAcceptance->addToZeroJetsOnePositiveElectron( configurationBr
                                                         * positiveTauMuonPass
                                              * ( tauToPionTimesTauToElectronBr
                                                  * negativeTauPionFail
                                                  + ( tauPairToElectronPairBr
                                              + tauToElectronTimesTauToMuonBr )
                                                     * negativeTauMuonFail ) );
    currentAcceptance->addToZeroJetsOneNegativeMuon( configurationBr
                                                     * negativeTauMuonPass
                                                  * ( tauToPionTimesTauToMuonBr
                                                      * positiveTauPionFail
                                              + ( tauToElectronTimesTauToMuonBr
                                                  + tauPairToMuonPairBr )
                                                     * positiveTauMuonFail ) );
    currentAcceptance->addToZeroJetsOnePositiveMuon( configurationBr
                                                     * positiveTauMuonPass
                                                  * ( tauToPionTimesTauToMuonBr
                                                      * negativeTauPionFail
                                              + ( tauToElectronTimesTauToMuonBr
                                                  + tauPairToMuonPairBr )
                                                     * negativeTauMuonFail ) );

    currentAcceptance->addToZeroJetsZeroLeptons( configurationBr
                                                 * ( tauPairToPionPairBr
                                                     * negativeTauPionFail
                                                     * positiveTauPionFail
                                              + ( tauToPionTimesTauToElectronBr
                                                  + tauToPionTimesTauToMuonBr )
                                                  * ( negativeTauPionFail
                                                      * positiveTauMuonFail
                                                        + positiveTauPionFail
                                                        * negativeTauMuonFail )
                                                 + ( tauPairToElectronPairBr
                                          + 2.0 * tauToElectronTimesTauToMuonBr
                                                     + tauPairToMuonPairBr )
                                                   * negativeTauMuonFail
                                                   * negativeTauMuonFail ) );
  }



  chargeSummedNeutralinoVirtualCascade::chargeSummedNeutralinoVirtualCascade(
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
                        true,
                        shortcut ),
    electroweakDecayerDecays(
                           electroweakDecayer->inspect_direct_decay_handler() )
  {
    negativeTauLeftHandedness = 0.5;
    // this will probably have to be changed if we ever do the off-shell decays
    // properly.

    // here firstBr & secondBr don't work as well as just getting the branching
    // ratios from the CppSLHA::particle_property_set pointers.

    directDownDistribution
    = new neutralino_three_body_decay( shortcut->getReadier(),
                                       shortcut->getCppSlha(),
                                       coloredDecayer,
                                       effectiveSquarkMass,
                                       electroweakDecayer,
                                       shortcut->getNeutralinoOne(),
                                       shortcut->getSdownL(),
                                       shortcut->getSdownR() );
    directUpDistribution
    = new neutralino_three_body_decay( shortcut->getReadier(),
                                       shortcut->getCppSlha(),
                                       coloredDecayer,
                                       effectiveSquarkMass,
                                       electroweakDecayer,
                                       shortcut->getNeutralinoOne(),
                                       shortcut->getSupL(),
                                       shortcut->getSupR() );
    directStrangeDistribution
    = new neutralino_three_body_decay( shortcut->getReadier(),
                                       shortcut->getCppSlha(),
                                       coloredDecayer,
                                       effectiveSquarkMass,
                                       electroweakDecayer,
                                       shortcut->getNeutralinoOne(),
                                       shortcut->getSstrangeL(),
                                       shortcut->getSstrangeR() );
    directCharmDistribution
    = new neutralino_three_body_decay( shortcut->getReadier(),
                                       shortcut->getCppSlha(),
                                       coloredDecayer,
                                       effectiveSquarkMass,
                                       electroweakDecayer,
                                       shortcut->getNeutralinoOne(),
                                       shortcut->getScharmL(),
                                       shortcut->getScharmR() );
    directBottomDistribution
    = new neutralino_three_body_decay( shortcut->getReadier(),
                                       shortcut->getCppSlha(),
                                       coloredDecayer,
                                       effectiveSquarkMass,
                                       electroweakDecayer,
                                       shortcut->getNeutralinoOne(),
                                       shortcut->get_sbottom_L(),
                                       shortcut->get_sbottom_R() );
    directElectronDistribution
    = new neutralino_three_body_decay( shortcut->getReadier(),
                                       shortcut->getCppSlha(),
                                       coloredDecayer,
                                       effectiveSquarkMass,
                                       electroweakDecayer,
                                       shortcut->getNeutralinoOne(),
                                       shortcut->getSelectronL(),
                                       shortcut->getSelectronR() );
    directMuonDistribution
    = new neutralino_three_body_decay( shortcut->getReadier(),
                                       shortcut->getCppSlha(),
                                       coloredDecayer,
                                       effectiveSquarkMass,
                                       electroweakDecayer,
                                       shortcut->getNeutralinoOne(),
                                       shortcut->getSmuonL(),
                                       shortcut->getSmuonR() );
    directTauDistribution
    = new neutralino_three_body_decay( shortcut->getReadier(),
                                       shortcut->getCppSlha(),
                                       coloredDecayer,
                                       effectiveSquarkMass,
                                       electroweakDecayer,
                                       shortcut->getNeutralinoOne(),
                                       shortcut->getStauOne(),
                                       shortcut->getStauTwo() );
    activeDistributions.push_back( directDownDistribution );
    activeDistributions.push_back( directUpDistribution );
    activeDistributions.push_back( directStrangeDistribution );
    activeDistributions.push_back( directCharmDistribution );
    activeDistributions.push_back( directBottomDistribution );
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

  chargeSummedNeutralinoVirtualCascade::~chargeSummedNeutralinoVirtualCascade()
  {
    // does nothing.
  }


  bool
  chargeSummedNeutralinoVirtualCascade::validSignal( int const numberOfJets,
                                           int const numberOfNegativeElectrons,
                                           int const numberOfPositiveElectrons,
                                               int const numberOfNegativeMuons,
                                              int const numberOfPositiveMuons )
  /* this returns true if a configuration where each of the signs of tau lepton
   * decayed either into a detected jet, detected lepton, or undetected
   * particle, & false otherwise.
   */
  {
    if( 2 == numberOfJets )
    {
      if( ( 0 == numberOfNegativeElectrons )
          &&
          ( 0 == numberOfPositiveElectrons )
          &&
          ( 0 == numberOfNegativeMuons )
          &&
          ( 0 == numberOfPositiveMuons ) )
      {
        return true;
      }
      else
      {
        return false;
      }
    }
    else if( ( 0 <= numberOfJets )
             &&
             ( 0 <= numberOfNegativeElectrons )
             &&
             ( 0 <= numberOfPositiveElectrons )
             &&
             ( 0 <= numberOfNegativeMuons )
             &&
             ( 0 <= numberOfPositiveMuons )
             &&
             ( 1 >= ( numberOfNegativeElectrons + numberOfNegativeMuons ) )
             &&
             ( 1 >= ( numberOfPositiveElectrons + numberOfPositiveMuons ) ) )
    {
      if( 2 >=
          ( numberOfJets
            + numberOfNegativeElectrons + numberOfPositiveElectrons
            + numberOfNegativeMuons + numberOfPositiveMuons ) )
      {
        return true;
      }
      else
      {
        return false;
      }
    }
    else
    {
      return false;
    }
  }

  void
  chargeSummedNeutralinoVirtualCascade::calculateAcceptance(
                                                  acceptanceCutSet* const cuts,
                                    acceptanceValues* const currentAcceptance )
  // this returns the appropriate acceptances multiplied by branching ratios.
  {
    currentAcceptance->setTwoJets( 0.0 );
    currentAcceptance->setOneJetOneNegativeElectron( 0.0 );
    currentAcceptance->setOneJetOnePositiveElectron( 0.0 );
    currentAcceptance->setOneJetOneNegativeMuon( 0.0 );
    currentAcceptance->setOneJetOnePositiveMuon( 0.0 );
    currentAcceptance->setElectronPlusAntielectron( 0.0 );
    currentAcceptance->setNegativeElectronPlusPositiveMuon( 0.0 );
    currentAcceptance->setNegativeMuonPlusPositiveElectron( 0.0 );
    currentAcceptance->setMuonPlusAntimuon( 0.0 );
    currentAcceptance->setOneJetZeroLeptons( 0.0 );
    currentAcceptance->setZeroJetsOneNegativeElectron( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveElectron( 0.0 );
    currentAcceptance->setZeroJetsOneNegativeMuon( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveMuon( 0.0 );
    currentAcceptance->setZeroJetsZeroLeptons( 0.0 );
    currentCuts = cuts;
    this->currentAcceptance = currentAcceptance;
    cascadeBr
    =  electroweakDecayerDecays->get_branching_ratio_for_exact_match(
                                             CppSLHA::PDG_code::neutralino_one,
                                                       CppSLHA::PDG_code::down,
                                                  -(CppSLHA::PDG_code::down) );
    if( lhcFaserGlobal::negligibleBr < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      currentJetDistribution = directDownDistribution;
      calculateForCurrentJetConfiguration();
    }
    cascadeBr
    =  electroweakDecayerDecays->get_branching_ratio_for_exact_match(
                                             CppSLHA::PDG_code::neutralino_one,
                                                         CppSLHA::PDG_code::up,
                                                    -(CppSLHA::PDG_code::up) );
    if( lhcFaserGlobal::negligibleBr < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      currentJetDistribution = directUpDistribution;
      calculateForCurrentJetConfiguration();
    }
    cascadeBr
    =  electroweakDecayerDecays->get_branching_ratio_for_exact_match(
                                             CppSLHA::PDG_code::neutralino_one,
                                                    CppSLHA::PDG_code::strange,
                                               -(CppSLHA::PDG_code::strange) );
    if( lhcFaserGlobal::negligibleBr < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      currentJetDistribution = directStrangeDistribution;
      calculateForCurrentJetConfiguration();
    }
    cascadeBr
    =  electroweakDecayerDecays->get_branching_ratio_for_exact_match(
                                             CppSLHA::PDG_code::neutralino_one,
                                                      CppSLHA::PDG_code::charm,
                                                 -(CppSLHA::PDG_code::charm) );
    if( lhcFaserGlobal::negligibleBr < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      currentJetDistribution = directCharmDistribution;
      calculateForCurrentJetConfiguration();
    }
    cascadeBr
    =  electroweakDecayerDecays->get_branching_ratio_for_exact_match(
                                             CppSLHA::PDG_code::neutralino_one,
                                                     CppSLHA::PDG_code::bottom,
                                                -(CppSLHA::PDG_code::bottom) );
    if( lhcFaserGlobal::negligibleBr < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      currentJetDistribution = directBottomDistribution;
      calculateForCurrentJetConfiguration();
    }
    cascadeBr
    =  electroweakDecayerDecays->get_branching_ratio_for_exact_match(
                                             CppSLHA::PDG_code::neutralino_one,
                                                   CppSLHA::PDG_code::electron,
                                              -(CppSLHA::PDG_code::electron) );
    if( lhcFaserGlobal::negligibleBr < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      directPass = integrateAcceptance( directElectronDistribution,
                                        cuts->getPrimaryLeptonCut() );
      directFail = ( 1.0 - integrateAcceptance( directElectronDistribution,
                                             cuts->getSecondaryLeptonCut() ) );
      currentPass = ( cascadeBr * directPass * directPass );
      currentAcceptance->addToOssfMinusOsdf( currentPass );
      currentAcceptance->addToElectronPlusAntielectron( currentPass );
      currentPass = ( cascadeBr * directPass * directFail );
      currentAcceptance->addToZeroJetsOneNegativeElectron( currentPass );
      currentAcceptance->addToZeroJetsOnePositiveElectron( currentPass );
      currentAcceptance->addToZeroJetsZeroLeptons( cascadeBr
                                                   * directFail * directFail );
    }
    cascadeBr
    =  electroweakDecayerDecays->get_branching_ratio_for_exact_match(
                                             CppSLHA::PDG_code::neutralino_one,
                                                       CppSLHA::PDG_code::muon,
                                                  -(CppSLHA::PDG_code::muon) );
    if( lhcFaserGlobal::negligibleBr < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      directPass = integrateAcceptance( directMuonDistribution,
                                        cuts->getPrimaryLeptonCut() );
      directFail = ( 1.0 - integrateAcceptance( directMuonDistribution,
                                             cuts->getSecondaryLeptonCut() ) );
      currentPass = ( cascadeBr * directPass * directPass );
      currentAcceptance->addToOssfMinusOsdf( currentPass );
      currentAcceptance->addToMuonPlusAntimuon( currentPass );
      currentPass = ( cascadeBr * directPass * directFail );
      currentAcceptance->addToZeroJetsOneNegativeMuon( currentPass );
      currentAcceptance->addToZeroJetsOnePositiveMuon( currentPass );
      currentAcceptance->addToZeroJetsZeroLeptons( cascadeBr
                                                   * directFail * directFail );
    }
    cascadeBr
    =  electroweakDecayerDecays->get_branching_ratio_for_exact_match(
                                             CppSLHA::PDG_code::neutralino_one,
                                                 CppSLHA::PDG_code::tau_lepton,
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

      currentAcceptance->addToTwoJets( cascadeBr * tauPairToPionPairBr
                                       * averagePionPass * averagePionPass );
      // until I can think of a better way to do this, we take an average of
      // the pass rates.
      currentPass = ( cascadeBr * tauToPionTimesTauToElectronBr
                      * ( hardPionPass * hardMuonPass
                          + softPionPass * softMuonPass ) );
      currentAcceptance->addToOneJetOneNegativeElectron( currentPass );
      currentAcceptance->addToOneJetOnePositiveElectron( currentPass );
      currentPass = ( cascadeBr * tauToPionTimesTauToMuonBr
                      * ( 2.0 * averagePionPass * averageMuonPass ) );
      currentAcceptance->addToOneJetOneNegativeMuon( currentPass );
      currentAcceptance->addToOneJetOnePositiveMuon( currentPass );
      currentPass = ( cascadeBr * 2.0 * averageMuonPass * averageMuonPass );
      currentAcceptance->addToElectronPlusAntielectron( tauPairToElectronPairBr
                                                        * currentPass );
      currentAcceptance->addToMuonPlusAntimuon( tauPairToMuonPairBr
                                                * currentPass );
      currentPass *= tauToElectronTimesTauToMuonBr;
      currentAcceptance->addToNegativeElectronPlusPositiveMuon( currentPass );
      currentAcceptance->addToNegativeMuonPlusPositiveElectron( currentPass );

      currentAcceptance->addToOneJetZeroLeptons( cascadeBr
                                                 * ( tauPairToPionPairBr
                                                     * ( 2.0 * averagePionPass
                                                         * averagePionFail )
                                              + ( tauToPionTimesTauToElectronBr
                                                  + tauToPionTimesTauToMuonBr )
                                                * ( 2.0 * averagePionPass
                                                       * averageMuonFail ) ) );
      currentPass = ( cascadeBr * averageMuonPass
                      * ( tauToPionTimesTauToElectronBr * averagePionFail
                  + ( tauPairToElectronPairBr + tauToElectronTimesTauToMuonBr )
                    * averageMuonFail ) );
      currentAcceptance->addToZeroJetsOneNegativeElectron( currentPass );
      currentAcceptance->addToZeroJetsOnePositiveElectron( currentPass );
      currentPass = ( cascadeBr * averageMuonPass
                      * ( tauToPionTimesTauToMuonBr * averagePionFail
                  + ( tauToElectronTimesTauToMuonBr + tauPairToMuonPairBr )
                    * averageMuonFail ) );
      currentAcceptance->addToZeroJetsOneNegativeMuon( currentPass );
      currentAcceptance->addToZeroJetsOnePositiveMuon( currentPass );

      currentAcceptance->addToZeroJetsZeroLeptons( configurationBr
                                                   * ( tauPairToPionPairBr
                                                       * averagePionFail
                                                       * averagePionFail
                                              + ( tauToPionTimesTauToElectronBr
                                                  + tauToPionTimesTauToMuonBr )
                                                    * ( 2.0 * averagePionFail
                                                            * averageMuonFail )
                                                   + ( tauPairToElectronPairBr
                                          + 2.0 * tauToElectronTimesTauToMuonBr
                                                       + tauPairToMuonPairBr )
                                                     * averageMuonFail
                                                     * averageMuonFail ) );
    }
    // debugging:
    /**std::cout
    << std::endl
    << "cascadeBr = " << cascadeBr;
    std::cout << std::endl;**/
  }

  void
  chargeSummedNeutralinoVirtualCascade::calculateForCurrentJetConfiguration()
  {
    directPass = integrateAcceptance( currentJetDistribution,
                                      cuts->getJetCut() );
    directFail = ( 1.0 - directPass );
    currentAcceptance->addToTwoJets( cascadeBr * directPass * directPass );
    currentAcceptance->addToOneJetZeroLeptons( 2.0 * cascadeBr * directPass
                                                               * directFail );
    currentAcceptance->addToZeroJetsZeroLeptons( cascadeBr * directFail
                                                           * directFail );
  }



  neutralinoThreeBodyCascade::neutralinoThreeBodyCascade(
                                leptonAcceptanceParameterSet* const kinematics,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                CppSLHA::particle_property_set const* const electroweakDecayer,
                                        inputHandler const* const shortcut ) :
    chargeSummedNeutralinoVirtualCascade( kinematics,
                                            coloredDecayer,
                                            electroweakDecayer,
                                            shortcut )
  {
    // just substituting in a chargeSummedNeutralinoThreeBodyCascade instead.
  }

  neutralinoThreeBodyCascade::~neutralinoThreeBodyCascade()
  {
    // does nothing.
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



  scoloredToZPlusScoloredCascade::scoloredToZPlusScoloredCascade(
                                leptonAcceptanceParameterSet* const kinematics,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const coloredDecayer,
                   CppSLHA::particle_property_set const* const lighterScolored,
                                        inputHandler const* const shortcut ) :
    electroweakCascade( kinematics,
                        effectiveSquarkMass,
                        coloredDecayer,
                        lighterScolored,
                        shortcut->getZ(),
                        true,
                        shortcut )
  {
    double
    twiceSquareOfWeakCosine
    = ( 2.0 * shortcut->getWeakCosine() * shortcut->getWeakCosine() );
    negativeTauLeftHandedness
    = ( ( 1 - 2.0 * twiceSquareOfWeakCosine
          + twiceSquareOfWeakCosine * twiceSquareOfWeakCosine )
        / ( 5.0 - 2.0 * twiceSquareOfWeakCosine
            + 2.0 * twiceSquareOfWeakCosine * twiceSquareOfWeakCosine ) );
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

    hardMuonDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     directMuonDistribution,
                                     shortcut->getHardMuonFromTau() );
    softMuonDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     directMuonDistribution,
                                     shortcut->getSoftMuonFromTau() );
    hardPionDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     directMuonDistribution,
                                     shortcut->getHardPionFromTau() );
    softPionDistribution
    = new visible_tau_decay_product( shortcut->getReadier(),
                                     directMuonDistribution,
                                     shortcut->getSoftPionFromTau() );
    activeDistributions.push_back( hardMuonDistribution );
    activeDistributions.push_back( softMuonDistribution );
    activeDistributions.push_back( hardPionDistribution );
    activeDistributions.push_back( softPionDistribution );
  }

  scoloredToZPlusScoloredCascade::~scoloredToZPlusScoloredCascade()
  {
    // does nothing.
  }


  bool
  scoloredToZPlusScoloredCascade::validSignal( int const numberOfJets,
                                           int const numberOfNegativeElectrons,
                                           int const numberOfPositiveElectrons,
                                               int const numberOfNegativeMuons,
                                              int const numberOfPositiveMuons )
  /* this returns true if a configuration where each of the signs of tau lepton
   * decayed either into a detected jet, detected lepton, or undetected
   * particle, & false otherwise.
   */
  {
    if( 2 == numberOfJets )
    {
      if( ( 0 == numberOfNegativeElectrons )
          &&
          ( 0 == numberOfPositiveElectrons )
          &&
          ( 0 == numberOfNegativeMuons )
          &&
          ( 0 == numberOfPositiveMuons ) )
      {
        return true;
      }
      else
      {
        return false;
      }
    }
    else if( ( 0 <= numberOfJets )
             &&
             ( 0 <= numberOfNegativeElectrons )
             &&
             ( 0 <= numberOfPositiveElectrons )
             &&
             ( 0 <= numberOfNegativeMuons )
             &&
             ( 0 <= numberOfPositiveMuons )
             &&
             ( 1 >= ( numberOfNegativeElectrons + numberOfNegativeMuons ) )
             &&
             ( 1 >= ( numberOfPositiveElectrons + numberOfPositiveMuons ) ) )
    {
      if( 2 >=
          ( numberOfJets
            + numberOfNegativeElectrons + numberOfPositiveElectrons
            + numberOfNegativeMuons + numberOfPositiveMuons ) )
      {
        return true;
      }
      else
      {
        return false;
      }
    }
    else
    {
      return false;
    }
  }

  void
  scoloredToZPlusScoloredCascade::calculateAcceptance(
                                                  acceptanceCutSet* const cuts,
                                    acceptanceValues* const currentAcceptance )
  // this returns the appropriate acceptances multiplied by branching ratios.
  {
    currentAcceptance->setTwoJets( 0.0 );
    currentAcceptance->setOneJetOneNegativeElectron( 0.0 );
    currentAcceptance->setOneJetOnePositiveElectron( 0.0 );
    currentAcceptance->setOneJetOneNegativeMuon( 0.0 );
    currentAcceptance->setOneJetOnePositiveMuon( 0.0 );
    currentAcceptance->setElectronPlusAntielectron( 0.0 );
    currentAcceptance->setNegativeElectronPlusPositiveMuon( 0.0 );
    currentAcceptance->setNegativeMuonPlusPositiveElectron( 0.0 );
    currentAcceptance->setMuonPlusAntimuon( 0.0 );
    currentAcceptance->setOneJetZeroLeptons( 0.0 );
    currentAcceptance->setZeroJetsOneNegativeElectron( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveElectron( 0.0 );
    currentAcceptance->setZeroJetsOneNegativeMuon( 0.0 );
    currentAcceptance->setZeroJetsOnePositiveMuon( 0.0 );
    currentAcceptance->setZeroJetsZeroLeptons( 0.0 );
    cascadeBr
    = ( electroweakDecayer->inspect_direct_decay_handler(
                  )->get_branching_ratio_for_exact_match( CppSLHA::PDG_code::Z,
                                             electroweakDecayer->get_PDG_code()
        /* the lighter scolored got put into electroweakDecayer, while the Z
         * boson got put into intermediateDecayer. */ ) );
    // the branching ratios of the Z into the various SM fermions are covered
    // by the following code.
    if( lhcFaserGlobal::negligibleBr < cascadeBr )
      // if the branching ratio into this channel is not negligible...
    {
      /* at some point, it'd be nice to break this down into transverse
       * contributions & longitudinal contributions, since actually the
       * Goldstone boson contribution could both be large (for sbottom -> stop,
       * for example) & has a different (& easier to calculate) distribution.
       *
       * however, ( ( m_top * m_tau ) / m_Z^2 ) is still pretty small, so it
       * should be fine to neglect the longitudinal contributions.
       */

      directMuonPass = integrateAcceptance( directMuonDistribution,
                                            cuts->getPrimaryLeptonCut() );
      directMuonFail = ( 1.0 - integrateAcceptance( directMuonDistribution,
                                             cuts->getSecondaryLeptonCut() ) );
      directJetPass = integrateAcceptance( directMuonDistribution,
                                           cuts->getJetCut() );
      directJetFail = ( 1.0 - directJetPass );

      configurationBr
      = ( cascadeBr * CppSLHA::PDG_data::Z_to_hadrons_BR * directJetPass );
      currentAcceptance->addToTwoJets( configurationBr * directJetPass );
      currentAcceptance->addToOneJetZeroLeptons( 2.0 * configurationBr
                                                     * directJetFail );
      currentAcceptance->addToOssfMinusOsdf( cascadeBr
                           * ( CppSLHA::PDG_data::Z_to_electron_antielectron_BR
                               + CppSLHA::PDG_data::Z_to_muon_antimuon_BR )
                             * directMuonPass * directMuonPass );
      configurationBr
      = ( cascadeBr * CppSLHA::PDG_data::Z_to_electron_antielectron_BR
          * directMuonPass );
      currentAcceptance->addToElectronPlusAntielectron( configurationBr
                                                        * directMuonPass );
      configurationBr *= directMuonFail;
      currentAcceptance->addToZeroJetsOneNegativeElectron( configurationBr );
      currentAcceptance->addToZeroJetsOnePositiveElectron( configurationBr );
      configurationBr
      = ( cascadeBr * CppSLHA::PDG_data::Z_to_muon_antimuon_BR
          * directMuonPass );
      currentAcceptance->addToMuonPlusAntimuon( configurationBr
                                                * directMuonPass );
      configurationBr *= directMuonFail;
      currentAcceptance->addToZeroJetsOneNegativeMuon( configurationBr );
      currentAcceptance->addToZeroJetsOnePositiveMuon( configurationBr );
      currentAcceptance->addToZeroJetsZeroLeptons( cascadeBr
                                         * ( CppSLHA::PDG_data::Z_to_hadrons_BR
                                             * directJetFail * directJetFail
                           + ( CppSLHA::PDG_data::Z_to_electron_antielectron_BR
                                   + CppSLHA::PDG_data::Z_to_muon_antimuon_BR )
                             * directMuonFail * directMuonFail ) );

      // now Z decays to tau-antitau pairs, followed by the decays of the taus:
      cascadeBr *= CppSLHA::PDG_data::Z_to_tau_lepton_antilepton_BR;
      // left-handed negative tau lepton:
      configurationBr
      = ( cascadeBr * negativeTauLeftHandedness );
      currentMuonDistribution = hardMuonDistribution;
      currentPionDistribution = softPionDistribution;
      calculateForCurrentConfiguration();

      // right-handed negative tau lepton:
      configurationBr
      = ( cascadeBr * ( 1.0 - negativeTauLeftHandedness ) );
      currentMuonDistribution = softMuonDistribution;
      currentPionDistribution = hardPionDistribution;
      calculateForCurrentConfiguration();
    }
    // debugging:
    /**std::cout
    << std::endl
    << "cascadeBr = " << cascadeBr;
    std::cout << std::endl;**/
  }

  void
  scoloredToZPlusScoloredCascade::calculateForCurrentConfiguration()
  {
    tauMuonPass = integrateAcceptance( currentMuonDistribution,
                                       currentCuts->getPrimaryLeptonCut() );
    tauMuonFail = ( 1.0 - integrateAcceptance( currentMuonDistribution,
                                      currentCuts->getSecondaryLeptonCut() ) );
    tauPionPass = integrateAcceptance( currentPionDistribution,
                                       currentCuts->getJetCut() );
    tauPionFail = ( 1.0 - tauPionPass );

    currentAcceptance->addToTwoJets( configurationBr
                                     * tauPairToPionPairBr
                                     * tauPionPass * tauPionPass );
    currentPass = ( configurationBr * tauToPionTimesTauToElectronBr
                    * tauPionPass * tauMuonPass );
    currentAcceptance->addToOneJetOneNegativeElectron( currentPass );
    currentAcceptance->addToOneJetOnePositiveElectron( currentPass );
    currentPass = ( configurationBr * tauToPionTimesTauToMuonBr
                    * tauPionPass * tauMuonPass );
    currentAcceptance->addToOneJetOneNegativeMuon( currentPass );
    currentAcceptance->addToOneJetOnePositiveMuon( currentPass );
    currentPass = ( configurationBr * tauMuonPass * tauMuonPass );
    currentAcceptance->addToElectronPlusAntielectron( tauPairToElectronPairBr
                                                      * currentPass );
    currentAcceptance->addToMuonPlusAntimuon( tauPairToMuonPairBr
                                              * currentPass );
    currentPass *= tauToElectronTimesTauToMuonBr;
    currentAcceptance->addToNegativeElectronPlusPositiveMuon( currentPass );
    currentAcceptance->addToNegativeMuonPlusPositiveElectron( currentPass );

    currentAcceptance->addToOneJetZeroLeptons( 2.0 * configurationBr
                                                   * tauPionPass
                                                   * ( tauPairToPionPairBr
                                                       * tauPionFail
                                              + ( tauToPionTimesTauToElectronBr
                                                  + tauToPionTimesTauToMuonBr )
                                                * tauMuonFail ) );
    currentPass = ( configurationBr * tauMuonPass
                    * ( tauToPionTimesTauToElectronBr * tauPionFail
                  + ( tauPairToElectronPairBr + tauToElectronTimesTauToMuonBr )
                    * tauMuonFail ) );
    currentAcceptance->addToZeroJetsOneNegativeElectron( currentPass );
    currentAcceptance->addToZeroJetsOnePositiveElectron( currentPass );
    currentPass = ( configurationBr * tauMuonPass
                    * ( tauToPionTimesTauToMuonBr * tauPionFail
                      + ( tauToElectronTimesTauToMuonBr + tauPairToMuonPairBr )
                        * tauMuonFail ) );
    currentAcceptance->addToZeroJetsOneNegativeMuon( currentPass );
    currentAcceptance->addToZeroJetsOnePositiveMuon( currentPass );

    currentAcceptance->addToZeroJetsZeroLeptons( configurationBr
                                                 * ( tauPairToPionPairBr
                                                   * nearPionFail * farPionFail
                                        + 2.0 * ( tauToPionTimesTauToElectronBr
                                                  + tauToPionTimesTauToMuonBr )
                                                  * tauPionFail * tauMuonFail
                                                 + ( tauPairToElectronPairBr
                                          + 2.0 * tauToElectronTimesTauToMuonBr
                                                     + tauPairToMuonPairBr )
                                              * nearMuonFail * farMuonFail ) );
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
