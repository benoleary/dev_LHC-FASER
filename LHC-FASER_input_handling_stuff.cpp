/*
 * LHC-FASER_input_handling_stuff.cpp
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

#include "LHC-FASER_input_handling_stuff.hpp"

namespace LHC_FASER
{

  input_handler::input_handler(
                          CppSLHA::CppSLHA0 const* const given_CppSLHA_pointer,
                CppSLHA::EW_scale_spectrum const* const given_spectrum_pointer,
                                std::string const given_path_to_grids,
                                readier_for_new_point* const given_readier ) :
    readier( given_readier ),
    verbosity_flag( false ),
    using_NLO_flag( true ),
    CppSLHA_pointer( given_CppSLHA_pointer ),
    spectrum_pointer( given_spectrum_pointer ),
    path_to_grids( given_path_to_grids ),
    path_to_cross_section_grids( given_path_to_grids ),
    path_to_kinematics_grids( given_path_to_grids ),
    down_pointer( given_spectrum_pointer->inspect_particle_property_set(
                                                   CppSLHA::PDG_code::down ) ),
    up_pointer( given_spectrum_pointer->inspect_particle_property_set(
                                                     CppSLHA::PDG_code::up ) ),
    strange_pointer( given_spectrum_pointer->inspect_particle_property_set(
                                                CppSLHA::PDG_code::strange ) ),
    charm_pointer( given_spectrum_pointer->inspect_particle_property_set(
                                                  CppSLHA::PDG_code::charm ) ),
    bottom_pointer( given_spectrum_pointer->inspect_particle_property_set(
                                                 CppSLHA::PDG_code::bottom ) ),
    top_pointer( given_spectrum_pointer->inspect_particle_property_set(
                                                    CppSLHA::PDG_code::top ) ),
    electron_pointer( given_spectrum_pointer->inspect_particle_property_set(
                                               CppSLHA::PDG_code::electron ) ),
    electron_neutrino_pointer(
                         given_spectrum_pointer->inspect_particle_property_set(
                                      CppSLHA::PDG_code::electron_neutrino ) ),
    muon_pointer( given_spectrum_pointer->inspect_particle_property_set(
                                                   CppSLHA::PDG_code::muon ) ),
    muon_neutrino_pointer(
                         given_spectrum_pointer->inspect_particle_property_set(
                                          CppSLHA::PDG_code::muon_neutrino ) ),
    tau_lepton_pointer( given_spectrum_pointer->inspect_particle_property_set(
                                             CppSLHA::PDG_code::tau_lepton ) ),
    tau_neutrino_pointer(
                         given_spectrum_pointer->inspect_particle_property_set(
                                           CppSLHA::PDG_code::tau_neutrino ) ),
    gluon_pointer( given_spectrum_pointer->inspect_particle_property_set(
                                                  CppSLHA::PDG_code::gluon ) ),
    photon_pointer( given_spectrum_pointer->inspect_particle_property_set(
                                                 CppSLHA::PDG_code::photon ) ),
    Z_pointer( given_spectrum_pointer->inspect_particle_property_set(
                                                      CppSLHA::PDG_code::Z ) ),
    W_plus_pointer( given_spectrum_pointer->inspect_particle_property_set(
                                                 CppSLHA::PDG_code::W_plus ) ),
    light_neutral_EWSB_scalar_pointer(
                         given_spectrum_pointer->inspect_particle_property_set(
                              CppSLHA::PDG_code::light_neutral_EWSB_scalar ) ),
    heavy_neutral_EWSB_scalar_pointer(
                         given_spectrum_pointer->inspect_particle_property_set(
                              CppSLHA::PDG_code::heavy_neutral_EWSB_scalar ) ),
    neutral_EWSB_pseudoscalar_pointer(
                         given_spectrum_pointer->inspect_particle_property_set(
                              CppSLHA::PDG_code::neutral_EWSB_pseudoscalar ) ),
    charged_EWSB_scalar_pointer(
                         given_spectrum_pointer->inspect_particle_property_set(
                                    CppSLHA::PDG_code::charged_EWSB_scalar ) ),
    gluino_pointer( given_spectrum_pointer->inspect_particle_property_set(
                                                 CppSLHA::PDG_code::gluino ) ),
    sdown_L_pointer( given_spectrum_pointer->inspect_particle_property_set(
                                                CppSLHA::PDG_code::sdown_L ) ),
    sup_L_pointer( given_spectrum_pointer->inspect_particle_property_set(
                                                  CppSLHA::PDG_code::sup_L ) ),
    sstrange_L_pointer( given_spectrum_pointer->inspect_particle_property_set(
                                             CppSLHA::PDG_code::sstrange_L ) ),
    scharm_L_pointer( given_spectrum_pointer->inspect_particle_property_set(
                                               CppSLHA::PDG_code::scharm_L ) ),
    sbottom_one_pointer( given_spectrum_pointer->inspect_particle_property_set(
                                            CppSLHA::PDG_code::sbottom_one ) ),
    stop_one_pointer( given_spectrum_pointer->inspect_particle_property_set(
                                               CppSLHA::PDG_code::stop_one ) ),
    sdown_R_pointer( given_spectrum_pointer->inspect_particle_property_set(
                                                CppSLHA::PDG_code::sdown_R ) ),
    sup_R_pointer( given_spectrum_pointer->inspect_particle_property_set(
                                                  CppSLHA::PDG_code::sup_R ) ),
    sstrange_R_pointer( given_spectrum_pointer->inspect_particle_property_set(
                                             CppSLHA::PDG_code::sstrange_R ) ),
    scharm_R_pointer( given_spectrum_pointer->inspect_particle_property_set(
                                               CppSLHA::PDG_code::scharm_R ) ),
    sbottom_two_pointer( given_spectrum_pointer->inspect_particle_property_set(
                                            CppSLHA::PDG_code::sbottom_two ) ),
    stop_two_pointer( given_spectrum_pointer->inspect_particle_property_set(
                                               CppSLHA::PDG_code::stop_two ) ),
    selectron_L_pointer( given_spectrum_pointer->inspect_particle_property_set(
                                            CppSLHA::PDG_code::selectron_L ) ),
    electron_sneutrino_L_pointer(
                         given_spectrum_pointer->inspect_particle_property_set(
                                   CppSLHA::PDG_code::electron_sneutrino_L ) ),
    smuon_L_pointer( given_spectrum_pointer->inspect_particle_property_set(
                                                CppSLHA::PDG_code::smuon_L ) ),
    muon_sneutrino_L_pointer(
                         given_spectrum_pointer->inspect_particle_property_set(
                                       CppSLHA::PDG_code::muon_sneutrino_L ) ),
    stau_one_pointer( given_spectrum_pointer->inspect_particle_property_set(
                                               CppSLHA::PDG_code::stau_one ) ),
    tau_sneutrino_L_pointer(
                         given_spectrum_pointer->inspect_particle_property_set(
                                        CppSLHA::PDG_code::tau_sneutrino_L ) ),
    selectron_R_pointer( given_spectrum_pointer->inspect_particle_property_set(
                                            CppSLHA::PDG_code::selectron_R ) ),
    electron_sneutrino_R_pointer(
                         given_spectrum_pointer->inspect_particle_property_set(
                                   CppSLHA::PDG_code::electron_sneutrino_R ) ),
    smuon_R_pointer( given_spectrum_pointer->inspect_particle_property_set(
                                                CppSLHA::PDG_code::smuon_R ) ),
    muon_sneutrino_R_pointer(
                         given_spectrum_pointer->inspect_particle_property_set(
                                       CppSLHA::PDG_code::muon_sneutrino_R ) ),
    stau_two_pointer( given_spectrum_pointer->inspect_particle_property_set(
                                               CppSLHA::PDG_code::stau_two ) ),
    tau_sneutrino_R_pointer(
                         given_spectrum_pointer->inspect_particle_property_set(
                                        CppSLHA::PDG_code::tau_sneutrino_R ) ),
    neutralino_one_pointer(
                         given_spectrum_pointer->inspect_particle_property_set(
                                         CppSLHA::PDG_code::neutralino_one ) ),
    neutralino_two_pointer(
                         given_spectrum_pointer->inspect_particle_property_set(
                                         CppSLHA::PDG_code::neutralino_two ) ),
    neutralino_three_pointer(
                         given_spectrum_pointer->inspect_particle_property_set(
                                       CppSLHA::PDG_code::neutralino_three ) ),
    neutralino_four_pointer(
                         given_spectrum_pointer->inspect_particle_property_set(
                                        CppSLHA::PDG_code::neutralino_four ) ),
    chargino_one_pointer(
                         given_spectrum_pointer->inspect_particle_property_set(
                                           CppSLHA::PDG_code::chargino_one ) ),
    chargino_two_pointer(
                         given_spectrum_pointer->inspect_particle_property_set(
                                           CppSLHA::PDG_code::chargino_two ) ),
                  decay_checkers( new decay_checker_handler( given_readier ) ),
    NMIX_pointer( given_CppSLHA_pointer->inspect_BLOCK( "NMIX" ) ),
    UMIX_pointer( given_CppSLHA_pointer->inspect_BLOCK( "UMIX" ) ),
    VMIX_pointer( given_CppSLHA_pointer->inspect_BLOCK( "VMIX" ) ),
    SBOTMIX_pointer( given_CppSLHA_pointer->inspect_BLOCK( "SBOTMIX" ) ),
    STOPMIX_pointer( given_CppSLHA_pointer->inspect_BLOCK( "STOPMIX" ) ),
    STAUMIX_pointer( given_CppSLHA_pointer->inspect_BLOCK( "STAUMIX" ) ),
    YD_pointer( given_CppSLHA_pointer->inspect_BLOCK( "YD" ) ),
    YU_pointer( given_CppSLHA_pointer->inspect_BLOCK( "YU" ) ),
    YE_pointer( given_CppSLHA_pointer->inspect_BLOCK( "YE" ) ),
    hard_muon_from_tau_function(),
    soft_muon_from_tau_function(),
    hard_pion_from_tau_function(),
    soft_pion_from_tau_function()
  {

    update_dependent_input_values
    = new update_dependent_inputs( given_CppSLHA_pointer->inspect_BLOCK(
                                                                  "SMINPUTS" ),
                                   W_plus_pointer,
                                   Z_pointer,
                                   given_readier );

    // the default cross-sections are those generated with the MSTW2008 PDF
    // set, & the default kinematics grids are those from PYTHIA 8.0:
    path_to_cross_section_grids.append( "/cross-sections/MSTW2008" );
    path_to_kinematics_grids.append( "/kinematics/PYTHIA8" );


    // find all the shortcut pointers to particles:

    colored_sparticles.push_back( gluino_pointer );
    sparticle_production_combinations.push_back(
                             new signed_particle_shortcut_pair( gluino_pointer,
                                                                          true,
                                                                gluino_pointer,
                                                                      true ) );
    add_squark_as_colored_sparticle_and_in_production( sdown_L_pointer );
    squarks4.push_back( sdown_L_pointer );
    squarks5.push_back( sdown_L_pointer );
    sdown_types.push_back( sdown_L_pointer );
    add_squark_as_colored_sparticle_and_in_production( sdown_R_pointer );
    squarks4.push_back( sdown_R_pointer );
    squarks5.push_back( sdown_R_pointer );
    sdown_types.push_back( sdown_R_pointer );
    add_squark_as_colored_sparticle_and_in_production( sup_L_pointer );
    squarks4.push_back( sup_L_pointer );
    squarks5.push_back( sup_L_pointer );
    sup_types.push_back( sup_L_pointer );
    add_squark_as_colored_sparticle_and_in_production( sup_R_pointer );
    squarks4.push_back( sup_R_pointer );
    squarks5.push_back( sup_R_pointer );
    sup_types.push_back( sup_R_pointer );
    add_squark_as_colored_sparticle_and_in_production( sstrange_L_pointer );
    squarks4.push_back( sstrange_L_pointer );
    squarks5.push_back( sstrange_L_pointer );
    sdown_types.push_back( sstrange_L_pointer );
    add_squark_as_colored_sparticle_and_in_production( sstrange_R_pointer );
    squarks4.push_back( sstrange_R_pointer );
    squarks5.push_back( sstrange_R_pointer );
    sdown_types.push_back( sstrange_R_pointer );
    add_squark_as_colored_sparticle_and_in_production( scharm_L_pointer );
    squarks4.push_back( scharm_L_pointer );
    squarks5.push_back( scharm_L_pointer );
    sup_types.push_back( scharm_L_pointer );
    add_squark_as_colored_sparticle_and_in_production( scharm_R_pointer );
    squarks4.push_back( scharm_R_pointer );
    squarks5.push_back( scharm_R_pointer );
    sup_types.push_back( scharm_R_pointer );
    add_squark_as_colored_sparticle_and_in_production( sbottom_one_pointer );
    squarks5.push_back( sbottom_one_pointer );
    sdown_types.push_back( sbottom_one_pointer );
    add_squark_as_colored_sparticle_and_in_production( sbottom_two_pointer );
    squarks5.push_back( sbottom_two_pointer );
    sdown_types.push_back( sbottom_two_pointer );
    add_squark_as_colored_sparticle_and_in_production( stop_one_pointer );
    sup_types.push_back( stop_one_pointer );
    add_squark_as_colored_sparticle_and_in_production( stop_two_pointer );
    sup_types.push_back( stop_two_pointer );

    sparticle_production_combinations.push_back(
                            new signed_particle_shortcut_pair( sdown_L_pointer,
                                                               true,
                                                               sdown_R_pointer,
                                                               false ) );
    sparticle_production_combinations.push_back(
                            new signed_particle_shortcut_pair( sdown_R_pointer,
                                                               true,
                                                               sdown_L_pointer,
                                                               false ) );
    sparticle_production_combinations.push_back(
                            new signed_particle_shortcut_pair( sdown_L_pointer,
                                                               true,
                                                               sup_R_pointer,
                                                               false ) );
    sparticle_production_combinations.push_back(
                            new signed_particle_shortcut_pair( sdown_R_pointer,
                                                               true,
                                                               sup_L_pointer,
                                                               false ) );
    sparticle_production_combinations.push_back(
                            new signed_particle_shortcut_pair( sup_L_pointer,
                                                               true,
                                                               sdown_R_pointer,
                                                               false ) );
    sparticle_production_combinations.push_back(
                            new signed_particle_shortcut_pair( sup_R_pointer,
                                                               true,
                                                               sdown_L_pointer,
                                                               false ) );
    sparticle_production_combinations.push_back(
                            new signed_particle_shortcut_pair( sup_L_pointer,
                                                               true,
                                                               sup_R_pointer,
                                                               false ) );
    sparticle_production_combinations.push_back(
                            new signed_particle_shortcut_pair( sup_R_pointer,
                                                               true,
                                                               sup_L_pointer,
                                                               false ) );
    sparticle_production_combinations.push_back(
                            new signed_particle_shortcut_pair( sdown_L_pointer,
                                                               true,
                                                               sup_L_pointer,
                                                               true ) );
    sparticle_production_combinations.push_back(
                            new signed_particle_shortcut_pair( sdown_R_pointer,
                                                               true,
                                                               sup_R_pointer,
                                                               true ) );

    exclusive_BRs = new exclusive_BR_handler( &sdown_types,
                                              &sup_types,
                                              given_readier );

    colored_cascades = new colored_cascade_handler( gluino_pointer,
                                                    &colored_sparticles,
                                                    &electroweakinos,
                                                    given_readier );


    // make all the shortcut sparticle lists:

    electroweakinos.push_back( neutralino_one_pointer );
    neutralinos.push_back( neutralino_one_pointer );
    electroweakinos.push_back( neutralino_two_pointer );
    unstable_electroweakinos.push_back( neutralino_two_pointer );
    neutralinos.push_back( neutralino_two_pointer );
    unstable_neutralinos.push_back( neutralino_two_pointer );
    electroweakinos.push_back( neutralino_three_pointer );
    unstable_electroweakinos.push_back( neutralino_three_pointer );
    neutralinos.push_back( neutralino_three_pointer );
    unstable_neutralinos.push_back( neutralino_three_pointer );
    electroweakinos.push_back( neutralino_four_pointer );
    unstable_electroweakinos.push_back( neutralino_four_pointer );
    neutralinos.push_back( neutralino_four_pointer );
    unstable_neutralinos.push_back( neutralino_four_pointer );
    electroweakinos.push_back( chargino_one_pointer );
    unstable_electroweakinos.push_back( chargino_one_pointer );
    charginos.push_back( chargino_one_pointer );
    electroweakinos.push_back( chargino_two_pointer );
    unstable_electroweakinos.push_back( chargino_two_pointer );
    charginos.push_back( chargino_two_pointer );

    charged_sleptons.push_back( selectron_L_pointer );
    charged_sleptons.push_back( selectron_R_pointer );
    charged_sleptons.push_back( smuon_L_pointer );
    charged_sleptons.push_back( smuon_R_pointer );
    charged_sleptons.push_back( stau_one_pointer );
    charged_sleptons.push_back( stau_two_pointer );

    sneutrinos.push_back( electron_sneutrino_L_pointer );
    sneutrinos.push_back( electron_sneutrino_R_pointer );
    sneutrinos.push_back( muon_sneutrino_L_pointer );
    sneutrinos.push_back( muon_sneutrino_R_pointer );
    sneutrinos.push_back( tau_sneutrino_L_pointer );
    sneutrinos.push_back( tau_sneutrino_R_pointer );

    neutral_EWSB_scalars_and_pseudoscalars.push_back(
                                           light_neutral_EWSB_scalar_pointer );
    neutral_EWSB_scalars_and_pseudoscalars.push_back(
                                           heavy_neutral_EWSB_scalar_pointer );
    neutral_EWSB_scalars_and_pseudoscalars.push_back(
                                           neutral_EWSB_pseudoscalar_pointer );

    // make the jet particle code lists:

    for( std::vector< CppSLHA::particle_property_set* >::const_iterator
           particle_iterator
             = spectrum_pointer->get_particle_property_sets()->begin();
         spectrum_pointer->get_particle_property_sets()->end()
             > particle_iterator;
         ++particle_iterator )
      // go through each particle in the spectrum...
      {

        if( (*particle_iterator)->counts_as_jet() )
          // if we find a jet...
          {

            jets.push_back( (*particle_iterator)->get_PDG_code() );
            // note the jet particle.

            if( !((*particle_iterator)->counts_as_self_conjugate()) )
              {

                jets.push_back( -((*particle_iterator)->get_PDG_code()) );

              }

          }

      }

    not_in_jets5.push_back( CppSLHA::PDG_code::top );
    not_in_jets5.push_back( -CppSLHA::PDG_code::top );
    not_in_jets4.assign( not_in_jets5.begin(), not_in_jets5.end() );
    not_in_jets4.push_back( CppSLHA::PDG_code::bottom );
    not_in_jets4.push_back( -CppSLHA::PDG_code::bottom );

    bool not_in_not_in_jets;

    for( std::list< int >::iterator jet_iterator = jets.begin();
         jets.end() != jet_iterator;
         ++jet_iterator )
      {

        // assume that this jet is not in not_in_jets5:
        not_in_not_in_jets = true;

        for( std::list< int >::const_iterator
               not_in_jets_iterator = not_in_jets5.begin();
             ( ( not_in_jets5.end() != not_in_jets_iterator )
               && not_in_not_in_jets );
             ++not_in_jets_iterator )
          {

            if( *not_in_jets_iterator == *jet_iterator )
              {

                not_in_not_in_jets = false;

              }

          }

        if( not_in_not_in_jets )
          {

            jets5.push_back( *jet_iterator );

          }

        // assume that this jet is not in not_in_jets4:
        not_in_not_in_jets = true;

        for( std::list< int >::const_iterator
               not_in_jets_iterator = not_in_jets4.begin();
             ( ( not_in_jets4.end() != not_in_jets_iterator )
               && not_in_not_in_jets );
             ++not_in_jets_iterator )
          {

            if( *not_in_jets_iterator == *jet_iterator )
              {

                not_in_not_in_jets = false;

              }

          }

        if( not_in_not_in_jets )
          {

            jets4.push_back( *jet_iterator );

          }

      }

  }

  input_handler::~input_handler()
  {

    for( std::vector< signed_particle_shortcut_pair* >::iterator
           deletion_iterator = sparticle_production_combinations.begin();
         sparticle_production_combinations.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

    delete decay_checkers;
    delete exclusive_BRs;
    delete colored_cascades;
    delete update_dependent_input_values;

  }


  void
  input_handler::add_squark_as_colored_sparticle_and_in_production(
                     CppSLHA::particle_property_set const* const given_squark )
  /* this adds the given squark to the list of squarks & the list of colored
   * sparticles, & also adds its relevant production combinations (with a
   * gluino, with its antisquark, with itself) to the list.
   */
  {

    squarks.push_back( given_squark );
    colored_sparticles.push_back( given_squark );

    sparticle_production_combinations.push_back(
                             new signed_particle_shortcut_pair( given_squark,
                                                                true,
                                                                gluino_pointer,
                                                                true ) );
    sparticle_production_combinations.push_back(
                             new signed_particle_shortcut_pair( given_squark,
                                                                true,
                                                                given_squark,
                                                                false ) );
    sparticle_production_combinations.push_back(
                             new signed_particle_shortcut_pair( given_squark,
                                                                true,
                                                                given_squark,
                                                                true ) );

  }


  bool
  input_handler::is_in( int const given_PDG_code,
              std::vector< CppSLHA::particle_property_set const* > const* const
                        given_collection )
  const
  /* this checks to see if the particle with the absolute value of the given
   * PDG particle code is in the given collection. (hence it checks to see if
   * the given particle or its antiparticle are in the collection.)
   */
   {

    bool return_bool = false;

    for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
           set_iterator = given_collection->begin();
         given_collection->end() > set_iterator;
         ++set_iterator )
      {

        if( ( given_PDG_code == (*set_iterator)->get_PDG_code() )
            ||
            ( given_PDG_code == -(*set_iterator)->get_PDG_code() ) )
          {

            return_bool = true;

          }

      }

    return return_bool;

   }


  double
  input_handler::quark_or_lepton_left_handedness( int const sfermion_code,
                                                  int const EWino_code )
  const
  /* this returns the square of the coupling associated with the left-handed
   * projection operator of the Feynman rule for the vertex of the given
   * sfermion with the given electroweakino & the appropriate quark, divided
   * by the sum of the squares of the couplings for the left- & right-handed
   * projection operators, based on Rosiek's Feynman rules (hence always for
   * incoming quarks & leptons). this has to be changed if general squark
   * flavor mixing is incorporated. also, it ignores the possibility of
   * right-handed neutrinos.
   */
  {

    /* I might as well keep some notes here:
     * Rosiek's paper with all the MSSM Feynman rules is great, but it has
     * a rather annoying convention on the mixing matrices, in that he uses
     * a different convention to most others, most notably the convention
     * used in the SLHA format.
     * neutralinos:
     * Rosiek's Z_N(j,k) is NMIX(k,j) if neutralino_k has positive mass in
     * the SLHA file. if it's negative in the SLHA file, then Z_N(j,k) is
     * instead -i * NMIX(k,j).
     * Rosiek's Z_+(j,k) is VMIX(k,j), & Rosiek's Z_-(j,k) is UMIX(k,j).
     * in the absence of CP violation, the U & V mixing matrices can be
     * chosen to be real.
     */

    // by the way, if we get a gluino, we assume it has the same handedness as
    // if we got a sbottom_1.

    if( CppSLHA::PDG_code::gluino == sfermion_code )
      {

        return quark_or_lepton_left_handedness( CppSLHA::PDG_code::sbottom_one,
                                                EWino_code );

      }
    else if( ( CppSLHA::PDG_code::chargino_one == EWino_code )
             ||
             ( CppSLHA::PDG_code::chargino_one == -EWino_code )
             ||
             ( CppSLHA::PDG_code::chargino_two == EWino_code )
             ||
             ( CppSLHA::PDG_code::chargino_two == -EWino_code ) )
      {

        // if it's an up-type squark, we assume that it's a positive chargino;
        // if it's a down-type squark, we assume that it's a negative chargino:
        if( ( CppSLHA::PDG_code::sup_R == sfermion_code )
            ||
            ( CppSLHA::PDG_code::scharm_R == sfermion_code )
            ||
            ( CppSLHA::PDG_code::sdown_R == sfermion_code )
            ||
            ( CppSLHA::PDG_code::sstrange_R == sfermion_code )
            ||
            ( CppSLHA::PDG_code::selectron_L == sfermion_code )
            ||
            ( CppSLHA::PDG_code::smuon_L == sfermion_code )
            ||
            ( CppSLHA::PDG_code::stau_one == sfermion_code )
            ||
            ( CppSLHA::PDG_code::selectron_R == sfermion_code )
            ||
            ( CppSLHA::PDG_code::smuon_R == sfermion_code )
            ||
            ( CppSLHA::PDG_code::stau_two == sfermion_code )
            ||
            ( CppSLHA::PDG_code::electron_sneutrino_R == sfermion_code )
            ||
            ( CppSLHA::PDG_code::muon_sneutrino_R == sfermion_code )
            ||
            ( CppSLHA::PDG_code::tau_sneutrino_R == sfermion_code ) )
          {

            return 1.0;
            /* the right-chiral squarks only interact with the Higgsino
             * component of charginos, so produce only *left*-handed quarks,
             * while charged sleptons produce only neutrinos in their
             * interactions with charginos, which are assumed to always be
             * left-handed. (right-chiral sneutrinos also only would interact
             * with the Higgsino components.)
             */

          }
        else
          {

            double left_coupling;
            double right_coupling;
            double left_squared;
            int which_chargino = 1;
            if( CppSLHA::PDG_code::chargino_two == EWino_code )
              {

                which_chargino = 2;

              }
            if( ( NULL == UMIX_pointer )
                ||
                ( NULL == VMIX_pointer )
                ||
                ( NULL == YD_pointer )
                ||
                ( NULL == YU_pointer )
                ||
                ( NULL == SBOTMIX_pointer )
                ||
                ( NULL == STOPMIX_pointer )
                ||
                ( NULL == STAUMIX_pointer ) )
              {

                if( is_verbose() )
                  {

                    std::cout
                    << std::endl << "LHC-FASER::warning!"
                    << " input_handler::quark_or_lepton_left_handedness"
                    << " couldn't find at least one of UMIX, VMIX, YD, YU,"
                    << " SBOTMIX, STOPMIX or STAUMIX. returning "
                    << CppSLHA::CppSLHA_global::really_wrong_value_string;
                    std::cout << std::endl;

                  }

                return CppSLHA::CppSLHA_global::really_wrong_value;

              }
            else if( CppSLHA::PDG_code::sup_L == sfermion_code )
              {

                left_coupling
                = ( get_g_two() * VMIX_pointer->get_entry( which_chargino,
                                                           1 ) );
                // actually the Feynman rule has an additional minus sign, but
                // it doesn't matter because this is getting squared.

                right_coupling = ( YD_pointer->get_entry( 1,
                                                          1 )
                                   * UMIX_pointer->get_entry( which_chargino,
                                                              2 ) );
                // actually, the Feynman rule has an additional minus sign, but
                // it doesn't matter because this is getting squared.

                left_squared = ( left_coupling * left_coupling );

                return ( left_squared
                         / ( left_squared
                             + right_coupling * right_coupling ) );

              }
            else if( CppSLHA::PDG_code::scharm_L == sfermion_code )
              {

                left_coupling
                = ( get_g_two() * VMIX_pointer->get_entry( which_chargino,
                                                           1 ) );
                // actually, the Feynman rule has an additional minus sign, but
                // it doesn't matter because this is getting squared.

                right_coupling = ( YD_pointer->get_entry( 2,
                                                          2 )
                                   * UMIX_pointer->get_entry( which_chargino,
                                                              2 ) );
                // actually, the Feynman rule has an additional minus sign, but
                // it doesn't matter because this is getting squared.

                left_squared = ( left_coupling * left_coupling );

                return ( left_squared
                         / ( left_squared
                             + right_coupling * right_coupling ) );

              }
            else if( CppSLHA::PDG_code::stop_one == sfermion_code )
              {

                left_coupling
                = ( YU_pointer->get_entry( 3,
                                           3 )
                    * STOPMIX_pointer->get_entry( 1,
                                                  2 )
                    * VMIX_pointer->get_entry( which_chargino,
                                               2 )
                    - get_g_two()
                    * STOPMIX_pointer->get_entry( 1,
                                                  1 )
                    * VMIX_pointer->get_entry( which_chargino,
                                               1 ) );
                // the (relative) signs of the terms do matter in this case.

                right_coupling
                = ( YD_pointer->get_entry( 3,
                                           3 )
                    * STOPMIX_pointer->get_entry( 1,
                                                  1 )
                   * UMIX_pointer->get_entry( which_chargino,
                                              2 ) );
                // actually, the Feynman rule has an additional minus sign, but
                // it doesn't matter because this is getting squared.

                left_squared = ( left_coupling * left_coupling );

                return ( left_squared
                         / ( left_squared
                             + right_coupling * right_coupling ) );

              }
            else if( CppSLHA::PDG_code::stop_two == sfermion_code )
              {

                left_coupling
                = ( YU_pointer->get_entry( 3,
                                           3 )
                    * STOPMIX_pointer->get_entry( 2,
                                                  2 )
                    * VMIX_pointer->get_entry( which_chargino,
                                               2 )
                    - get_g_two()
                      * STOPMIX_pointer->get_entry( 2,
                                                    1 )
                      * VMIX_pointer->get_entry( which_chargino,
                                                 1 ) );
                // the (relative) signs of the terms do matter in this case.

                right_coupling
                = ( YD_pointer->get_entry( 3,
                                           3 )
                    * STOPMIX_pointer->get_entry( 2,
                                                  1 )
                    * UMIX_pointer->get_entry( which_chargino,
                                               2 ) );
                // actually, the Feynman rule has an additional minus sign, but
                // it doesn't matter because this is getting squared.

                left_squared = ( left_coupling * left_coupling );

                return ( left_squared
                         / ( left_squared
                             + right_coupling * right_coupling ) );

              }
            else if( CppSLHA::PDG_code::sdown_L == sfermion_code )
              {

                left_coupling
                = ( get_g_two() * UMIX_pointer->get_entry( which_chargino,
                                                           1 ) );
                // actually the Feynman rule has an additional minus sign, but
                // it doesn't matter because this is getting squared.

                right_coupling = ( YU_pointer->get_entry( 1,
                                                          1 )
                                   * VMIX_pointer->get_entry( which_chargino,
                                                              2 ) );
                // actually, the Feynman rule has an additional minus sign, but
                // it doesn't matter because this is getting squared.

                left_squared = ( left_coupling * left_coupling );

                return ( left_squared
                         / ( left_squared
                             + right_coupling * right_coupling ) );

              }
            else if( CppSLHA::PDG_code::sstrange_L == sfermion_code )
              {

                left_coupling
                = ( get_g_two() * UMIX_pointer->get_entry( which_chargino,
                                                           1 ) );
                // actually, the Feynman rule has an additional minus sign, but
                // it doesn't matter because this is getting squared.

                right_coupling = ( YU_pointer->get_entry( 2,
                                                          2 )
                                   * VMIX_pointer->get_entry( which_chargino,
                                                              2 ) );
                // actually, the Feynman rule has an additional minus sign, but
                // it doesn't matter because this is getting squared.

                left_squared = ( left_coupling * left_coupling );

                return ( left_squared
                         / ( left_squared
                             + right_coupling * right_coupling ) );

              }
            else if( CppSLHA::PDG_code::sbottom_one == sfermion_code )
              {

                left_coupling
                = ( YD_pointer->get_entry( 3,
                                           3 )
                    * SBOTMIX_pointer->get_entry( 1,
                                                  2 )
                    * UMIX_pointer->get_entry( which_chargino,
                                               2 )
                    - get_g_two()
                      * SBOTMIX_pointer->get_entry( 1,
                                                    1 )
                      * UMIX_pointer->get_entry( which_chargino,
                                                 1 ) );
                // the (relative) signs of the terms do matter in this case.

                right_coupling
                = ( YU_pointer->get_entry( 3,
                                           3 )
                    * SBOTMIX_pointer->get_entry( 1,
                                                  1 )
                    * VMIX_pointer->get_entry( which_chargino,
                                               2 ) );
                // actually, the Feynman rule has an additional minus sign, but
                // it doesn't matter because this is getting squared.

                left_squared = ( left_coupling * left_coupling );

                return ( left_squared
                         / ( left_squared
                             + right_coupling * right_coupling ) );

              }
            else if( CppSLHA::PDG_code::sbottom_two == sfermion_code )
              {

                left_coupling
                = ( YD_pointer->get_entry( 3,
                                           3 )
                    * SBOTMIX_pointer->get_entry( 2,
                                                  2 )
                    * UMIX_pointer->get_entry( which_chargino,
                                               2 )
                    - get_g_two()
                      * SBOTMIX_pointer->get_entry( 2,
                                                    1 )
                      * UMIX_pointer->get_entry( which_chargino,
                                                 1 ) );
                // the (relative) signs of the terms do matter in this case.

                right_coupling
                = ( YU_pointer->get_entry( 3,
                                           3 )
                    * SBOTMIX_pointer->get_entry( 2,
                                                  1 )
                    * VMIX_pointer->get_entry( which_chargino,
                                               2 ) );
                // actually, the Feynman rule has an additional minus sign, but
                // it doesn't matter because this is getting squared.

                left_squared = ( left_coupling * left_coupling );

                return ( left_squared
                         / ( left_squared
                             + right_coupling * right_coupling ) );

              }
            else if( CppSLHA::PDG_code::electron_sneutrino_L == sfermion_code )
              {

                left_coupling
                = ( get_g_two() * VMIX_pointer->get_entry( which_chargino,
                                                           1 ) );
                // actually the Feynman rule has an additional minus sign, but
                // it doesn't matter because this is getting squared.

                right_coupling = ( YE_pointer->get_entry( 1,
                                                          1 )
                                   * UMIX_pointer->get_entry( which_chargino,
                                                              2 ) );
                // actually, the Feynman rule has an additional minus sign, but
                // it doesn't matter because this is getting squared.

                left_squared = ( left_coupling * left_coupling );

                return ( left_squared
                         / ( left_squared
                             + right_coupling * right_coupling ) );

              }
            else if( CppSLHA::PDG_code::muon_sneutrino_L == sfermion_code )
              {

                left_coupling
                = ( get_g_two() * VMIX_pointer->get_entry( which_chargino,
                                                           1 ) );
                // actually, the Feynman rule has an additional minus sign, but
                // it doesn't matter because this is getting squared.

                right_coupling = ( YE_pointer->get_entry( 2,
                                                          2 )
                                   * UMIX_pointer->get_entry( which_chargino,
                                                              2 ) );
                // actually, the Feynman rule has an additional minus sign, but
                // it doesn't matter because this is getting squared.

                left_squared = ( left_coupling * left_coupling );

                return ( left_squared
                         / ( left_squared
                             + right_coupling * right_coupling ) );

              }
            else if( CppSLHA::PDG_code::tau_sneutrino_L == sfermion_code )
              {

                left_coupling
                = ( get_g_two() * VMIX_pointer->get_entry( which_chargino,
                                                           1 ) );
                // actually, the Feynman rule has an additional minus sign, but
                // it doesn't matter because this is getting squared.

                right_coupling = ( YE_pointer->get_entry( 3,
                                                          3 )
                                   * UMIX_pointer->get_entry( which_chargino,
                                                              2 ) );
                // actually, the Feynman rule has an additional minus sign, but
                // it doesn't matter because this is getting squared.

                left_squared = ( left_coupling * left_coupling );

                return ( left_squared
                         / ( left_squared
                             + right_coupling * right_coupling ) );

              }
            else
              {

                if( is_verbose() )
                  {

                    std::cout
                    << std::endl << "LHC-FASER::warning!"
                    << " input_handler::quark_or_lepton_left_handedness"
                    << " didn't recognize the sfermion PDG code \""
                    << sfermion_code << "\". returning "
                    << CppSLHA::CppSLHA_global::really_wrong_value_string;
                    std::cout << std::endl;

                  }

                return CppSLHA::CppSLHA_global::really_wrong_value;

              }

          }  // end of if-else over the types of sfermion.

      }  // end of if it's a chargino.
    else
      // otherwise we assume that it's a neutralino...
      {

        int which_neutralino = 1;
        // we start by assuming that it's neutralino_1.
        if( CppSLHA::PDG_code::neutralino_two == EWino_code )
          {

            which_neutralino = 2;

          }
        else if( CppSLHA::PDG_code::neutralino_three == EWino_code )
          {

            which_neutralino = 3;

          }
        else if( CppSLHA::PDG_code::neutralino_three == EWino_code )
          {

            which_neutralino = 4;

          }

        double left_coupling;
        double right_coupling;
        double left_squared;
        if( ( NULL == NMIX_pointer )
            ||
            ( NULL == YD_pointer )
            ||
            ( NULL == YU_pointer )
            ||
            ( NULL == SBOTMIX_pointer )
            ||
            ( NULL == STOPMIX_pointer )
            ||
            ( NULL == STAUMIX_pointer ) )
          {

            if( is_verbose() )
              {

                std::cout
                << std::endl << "LHC-FASER::warning!"
                << " input_handler::quark_or_lepton_left_handedness"
                << " couldn't find at least one of NMIX, YD, YU,"
                << " SBOTMIX, STOPMIX or STAUMIX. returning "
                << CppSLHA::CppSLHA_global::really_wrong_value_string;
                std::cout << std::endl;

              }

            return CppSLHA::CppSLHA_global::really_wrong_value;

          }
        else if( CppSLHA::PDG_code::sdown_L == sfermion_code )
          {

            left_coupling
            = ( M_SQRT1_2 * ( get_g_two()
                              * NMIX_pointer->get_entry( which_neutralino,
                                                         2 )
                              - ( get_g_one() / 3.0 )
                                * NMIX_pointer->get_entry( which_neutralino,
                                                           1 ) ) );
            // by the way, M_SQRT1_2 is 1 over the square root of 2, defined in
            // math.h.

            right_coupling = ( YD_pointer->get_entry( 1,
                                                      1 )
                               * NMIX_pointer->get_entry( which_neutralino,
                                                          3 ) );

            left_squared = ( left_coupling * left_coupling );

            return ( left_squared
                     / ( left_squared
                         + right_coupling * right_coupling ) );

          }
        else if( CppSLHA::PDG_code::sdown_R == sfermion_code )
          {

            left_coupling
            = ( YD_pointer->get_entry( 1,
                                       1 )
                * NMIX_pointer->get_entry( which_neutralino,
                                           3 ) );

            right_coupling = ( M_SQRT1_2 * ( get_g_one() / 1.5 )
                               * NMIX_pointer->get_entry( which_neutralino,
                                                          1 ) );
            // actually, the Feynman rule has an additional minus sign, but it
            // doesn't matter because this is getting squared.

            left_squared = ( left_coupling * left_coupling );

            return ( left_squared
                     / ( left_squared
                         + right_coupling * right_coupling ) );

          }
        else if( CppSLHA::PDG_code::sup_L == sfermion_code )
          {

            left_coupling
            = ( M_SQRT1_2 * ( get_g_two()
                              * NMIX_pointer->get_entry( which_neutralino,
                                                         2 )
                              - ( get_g_one() / 3.0 )
                                * NMIX_pointer->get_entry( which_neutralino,
                                                           1 ) ) );

            right_coupling = ( YU_pointer->get_entry( 1,
                                                      1 )
                               * NMIX_pointer->get_entry( which_neutralino,
                                                          4 ) );
            // actually, the Feynman rule has an additional minus sign, but it
            // doesn't matter because this is getting squared.

            left_squared = ( left_coupling * left_coupling );

            return ( left_squared
                     / ( left_squared
                         + right_coupling * right_coupling ) );

          }
        else if( CppSLHA::PDG_code::sup_R == sfermion_code )
          {

            left_coupling
            = ( YU_pointer->get_entry( 1,
                                       1 )
                * NMIX_pointer->get_entry( which_neutralino,
                                           4 ) );
            // actually, the Feynman rule has an additional minus sign, but it
            // doesn't matter because this is getting squared.

            right_coupling = ( M_SQRT1_2 * ( get_g_one() / 0.75 )
                               * NMIX_pointer->get_entry( which_neutralino,
                                                          1 ) );

            left_squared = ( left_coupling * left_coupling );

            return ( left_squared
                     / ( left_squared
                         + right_coupling * right_coupling ) );

          }
        else if( CppSLHA::PDG_code::sstrange_L == sfermion_code )
          {

            left_coupling
            = ( M_SQRT1_2 * ( get_g_two()
                              * NMIX_pointer->get_entry( which_neutralino,
                                                         2 )
                              - ( get_g_one() / 3.0 )
                                * NMIX_pointer->get_entry( which_neutralino,
                                                           1 ) ) );

            right_coupling = ( YD_pointer->get_entry( 2,
                                                      2 )
                               * NMIX_pointer->get_entry( which_neutralino,
                                                          3 ) );

            left_squared = ( left_coupling * left_coupling );

            return ( left_squared
                     / ( left_squared
                         + right_coupling * right_coupling ) );

          }
        else if( CppSLHA::PDG_code::sstrange_R == sfermion_code )
          {

            left_coupling
            = ( YD_pointer->get_entry( 2,
                                       2 )
                * NMIX_pointer->get_entry( which_neutralino,
                                           3 ) );

            right_coupling = ( M_SQRT1_2 * ( get_g_one() / 1.5 )
                               * NMIX_pointer->get_entry( which_neutralino,
                                                          1 ) );
            // actually, the Feynman rule has an additional minus sign, but it
            // doesn't matter because this is getting squared.

            left_squared = ( left_coupling * left_coupling );

            return ( left_squared
                     / ( left_squared
                         + right_coupling * right_coupling ) );

          }
        else if( CppSLHA::PDG_code::scharm_L == sfermion_code )
          {

            left_coupling
            = ( M_SQRT1_2 * ( get_g_two()
                              * NMIX_pointer->get_entry( which_neutralino,
                                                         2 )
                              - ( get_g_one() / 3.0 )
                                * NMIX_pointer->get_entry( which_neutralino,
                                                           1 ) ) );

            right_coupling = ( YU_pointer->get_entry( 2,
                                                      2 )
                               * NMIX_pointer->get_entry( which_neutralino,
                                                          4 ) );
            // actually, the Feynman rule has an additional minus sign, but it
            // doesn't matter because this is getting squared.

            left_squared = ( left_coupling * left_coupling );

            return ( left_squared
                     / ( left_squared
                         + right_coupling * right_coupling ) );

          }
        else if( CppSLHA::PDG_code::scharm_R == sfermion_code )
          {

            left_coupling
            = ( YU_pointer->get_entry( 2,
                                       2 )
                * NMIX_pointer->get_entry( which_neutralino,
                                           4 ) );
            // actually, the Feynman rule has an additional minus sign, but it
            // doesn't matter because this is getting squared.

            right_coupling = ( M_SQRT1_2 * ( get_g_one() / 0.75 )
                               * NMIX_pointer->get_entry( which_neutralino,
                                                          1 ) );

            left_squared = ( left_coupling * left_coupling );

            return ( left_squared
                     / ( left_squared
                         + right_coupling * right_coupling ) );

          }
        else if( CppSLHA::PDG_code::sbottom_one == sfermion_code )
          {

            left_coupling
            = ( ( M_SQRT1_2
                  * ( get_g_two() * NMIX_pointer->get_entry( which_neutralino,
                                                             2 )
                      - ( get_g_one() / 3.0 )
                        * NMIX_pointer->get_entry( which_neutralino,
                                                   1 ) ) )
                    * SBOTMIX_pointer->get_entry( 1,
                                                  1 )
                    + ( YD_pointer->get_entry( 3,
                                               3 )
                        * NMIX_pointer->get_entry( which_neutralino,
                                                   3 ) )
                      * SBOTMIX_pointer->get_entry( 1,
                                                    2 ) );
            // the (relative) signs of the terms do matter in this case.

            right_coupling = ( ( YD_pointer->get_entry( 3,
                                                        3 )
                                 * NMIX_pointer->get_entry( which_neutralino,
                                                            3 ) )
                               * SBOTMIX_pointer->get_entry( 1,
                                                             1 )
                               - ( M_SQRT1_2 * ( get_g_one() / 1.5 )
                                   * NMIX_pointer->get_entry( which_neutralino,
                                                              1 ) )
                                 * SBOTMIX_pointer->get_entry( 1,
                                                               2 ) );
            // the (relative) signs of the terms do matter in this case.

            left_squared = ( left_coupling * left_coupling );

            return ( left_squared
                     / ( left_squared
                         + right_coupling * right_coupling ) );

          }
        else if( CppSLHA::PDG_code::sbottom_two == sfermion_code )
          {

            left_coupling
            = ( ( M_SQRT1_2
                  * ( get_g_two() * NMIX_pointer->get_entry( which_neutralino,
                                                             2 )
                      - ( get_g_one()  / 3.0 )
                        * NMIX_pointer->get_entry( which_neutralino,
                                                   1 ) ) )
                * SBOTMIX_pointer->get_entry( 2,
                                              1 )
                + ( YD_pointer->get_entry( 3,
                                           3 )
                    * NMIX_pointer->get_entry( which_neutralino,
                                               3 ) )
                  * SBOTMIX_pointer->get_entry( 2,
                                                2 ) );
            // the (relative) signs of the terms do matter in this case.

            right_coupling = ( ( YD_pointer->get_entry( 3,
                                                        3 )
                                 * NMIX_pointer->get_entry( which_neutralino,
                                                            3 ) )
                               * SBOTMIX_pointer->get_entry( 2,
                                                             1 )
                               - ( M_SQRT1_2 * ( get_g_one() / 1.5 )
                                   * NMIX_pointer->get_entry( which_neutralino,
                                                              1 ) )
                                 * SBOTMIX_pointer->get_entry( 2,
                                                               2 ) );
            // the (relative) signs of the terms do matter in this case.

            left_squared = ( left_coupling * left_coupling );

            return ( left_squared
                     / ( left_squared
                         + right_coupling * right_coupling ) );

          }
        else if( CppSLHA::PDG_code::stop_one == sfermion_code )
          {

            left_coupling
            = ( ( M_SQRT1_2
                  * ( get_g_two() * NMIX_pointer->get_entry( which_neutralino,
                                                             2 )
                      - ( get_g_one() / 3.0 )
                        * NMIX_pointer->get_entry( which_neutralino,
                                                   1 ) ) )
                * STOPMIX_pointer->get_entry( 1,
                                              1 )
                - ( YU_pointer->get_entry( 3,
                                           3 )
                    * NMIX_pointer->get_entry( which_neutralino,
                                               4 ) )
                  * STOPMIX_pointer->get_entry( 1,
                                                2 ) );
            // the (relative) signs of the terms do matter in this case.

            right_coupling = ( ( M_SQRT1_2 * ( get_g_one() / 0.75 )
                                 * NMIX_pointer->get_entry( which_neutralino,
                                                            1 ) )
                               * STOPMIX_pointer->get_entry( 1,
                                                         2 )
                                   - ( YU_pointer->get_entry( 3,
                                                      3 )
                                       * NMIX_pointer->get_entry( which_neutralino,
                                                          4 ) )
                                     * STOPMIX_pointer->get_entry( 1,
                                                           1 ) );
            // the (relative) signs of the terms do matter in this case.

            left_squared = ( left_coupling * left_coupling );

            return ( left_squared
                     / ( left_squared
                         + right_coupling * right_coupling ) );

          }
        else if( CppSLHA::PDG_code::stop_two == sfermion_code )
          {

            left_coupling
            = ( ( M_SQRT1_2
                  * ( get_g_two() * NMIX_pointer->get_entry( which_neutralino,
                                                             2 )
                      - ( get_g_one() / 3.0 )
                        * NMIX_pointer->get_entry( which_neutralino,
                                                   1 ) ) )
                * STOPMIX_pointer->get_entry( 2,
                                              1 )
                - ( YU_pointer->get_entry( 3,
                                           3 )
                    * NMIX_pointer->get_entry( which_neutralino,
                                               4 ) )
                  * STOPMIX_pointer->get_entry( 2,
                                                2 ) );
            // the (relative) signs of the terms do matter in this case.

            right_coupling = ( ( M_SQRT1_2 * ( get_g_one() / 0.75 )
                                 * NMIX_pointer->get_entry( which_neutralino,
                                                            1 ) )
                               * STOPMIX_pointer->get_entry( 2,
                                                             2 )
                              - ( YU_pointer->get_entry( 3,
                                                         3 )
                                  * NMIX_pointer->get_entry( which_neutralino,
                                                             4 ) )
                                * STOPMIX_pointer->get_entry( 2,
                                                              1 ) );
            // the (relative) signs of the terms do matter in this case.

            left_squared = ( left_coupling * left_coupling );

            return ( left_squared
                     / ( left_squared
                         + right_coupling * right_coupling ) );

          }
        else if( CppSLHA::PDG_code::selectron_L == sfermion_code )
          {

            left_coupling
            = ( M_SQRT1_2 * ( get_g_two()
                              * NMIX_pointer->get_entry( which_neutralino,
                                                         2 )
                              + get_g_one()
                                * NMIX_pointer->get_entry( which_neutralino,
                                                           1 ) ) );

            right_coupling = ( YE_pointer->get_entry( 1,
                                                      1 )
                               * NMIX_pointer->get_entry( which_neutralino,
                                                          3 ) );

            left_squared = ( left_coupling * left_coupling );

            return ( left_squared
                     / ( left_squared
                         + right_coupling * right_coupling ) );

          }
        else if( CppSLHA::PDG_code::selectron_R == sfermion_code )
          {

            left_coupling
            = ( YE_pointer->get_entry( 1,
                                       1 )
                * NMIX_pointer->get_entry( which_neutralino,
                                           3 ) );

            right_coupling = ( M_SQRT1_2 * ( get_g_one() * 2.0 )
                               * NMIX_pointer->get_entry( which_neutralino,
                                                          1 ) );
            // actually, the Feynman rule has an additional minus sign, but it
            // doesn't matter because this is getting squared.

            left_squared = ( left_coupling * left_coupling );

            return ( left_squared
                     / ( left_squared
                         + right_coupling * right_coupling ) );

          }
        else if( CppSLHA::PDG_code::smuon_L == sfermion_code )
          {

            left_coupling
            = ( M_SQRT1_2 * ( get_g_two()
                              * NMIX_pointer->get_entry( which_neutralino,
                                                         2 )
                              + get_g_one()
                                * NMIX_pointer->get_entry( which_neutralino,
                                                           1 ) ) );

            right_coupling = ( YE_pointer->get_entry( 2,
                                                      2 )
                               * NMIX_pointer->get_entry( which_neutralino,
                                                          3 ) );

            left_squared = ( left_coupling * left_coupling );

            return ( left_squared
                     / ( left_squared
                         + right_coupling * right_coupling ) );

          }
        else if( CppSLHA::PDG_code::smuon_R == sfermion_code )
          {

            left_coupling
            = ( YE_pointer->get_entry( 2,
                                       2 )
                * NMIX_pointer->get_entry( which_neutralino,
                                           3 ) );

            right_coupling = ( M_SQRT1_2 * ( get_g_one() * 2.0 )
                               * NMIX_pointer->get_entry( which_neutralino,
                                                          1 ) );
            // actually, the Feynman rule has an additional minus sign, but it
            // doesn't matter because this is getting squared.

            left_squared = ( left_coupling * left_coupling );

            return ( left_squared
                     / ( left_squared
                         + right_coupling * right_coupling ) );

          }
        else if( CppSLHA::PDG_code::stau_one == sfermion_code )
          {

            left_coupling
            = ( ( M_SQRT1_2
                  * ( get_g_two() * NMIX_pointer->get_entry( which_neutralino,
                                                       2 )
                      + get_g_one()
                        * NMIX_pointer->get_entry( which_neutralino,
                                                   1 ) ) )
                    * STAUMIX_pointer->get_entry( 1,
                                                  1 )
                    + ( YE_pointer->get_entry( 3,
                                               3 )
                        * NMIX_pointer->get_entry( which_neutralino,
                                                   3 ) )
                      * STAUMIX_pointer->get_entry( 1,
                                                    2 ) );
            // the (relative) signs of the terms do matter in this case.

            right_coupling = ( ( YE_pointer->get_entry( 3,
                                                        3 )
                                 * NMIX_pointer->get_entry( which_neutralino,
                                                            3 ) )
                               * STAUMIX_pointer->get_entry( 1,
                                                             1 )
                               - ( M_SQRT1_2 * ( get_g_one() * 2.0 )
                                   * NMIX_pointer->get_entry( which_neutralino,
                                                              1 ) )
                                 * STAUMIX_pointer->get_entry( 1,
                                                               2 ) );
            // the (relative) signs of the terms do matter in this case.

            left_squared = ( left_coupling * left_coupling );

            return ( left_squared
                     / ( left_squared
                         + right_coupling * right_coupling ) );

          }
        else if( CppSLHA::PDG_code::stau_two == sfermion_code )
          {

            left_coupling
            = ( ( M_SQRT1_2
                  * ( get_g_two() * NMIX_pointer->get_entry( which_neutralino,
                                                             2 )
                      + get_g_one()
                        * NMIX_pointer->get_entry( which_neutralino,
                                                   1 ) ) )
                    * STAUMIX_pointer->get_entry( 2,
                                                  1 )
                    + ( YE_pointer->get_entry( 3,
                                               3 )
                        * NMIX_pointer->get_entry( which_neutralino,
                                                   3 ) )
                      * STAUMIX_pointer->get_entry( 2,
                                                    2 ) );
            // the (relative) signs of the terms do matter in this case.

            right_coupling = ( ( YE_pointer->get_entry( 3,
                                                        3 )
                                 * NMIX_pointer->get_entry( which_neutralino,
                                                            3 ) )
                               * STAUMIX_pointer->get_entry( 2,
                                                             1 )
                               - ( M_SQRT1_2 * ( get_g_one() * 2.0 )
                                   * NMIX_pointer->get_entry( which_neutralino,
                                                              1 ) )
                                 * STAUMIX_pointer->get_entry( 2,
                                                               2 ) );
            // the (relative) signs of the terms do matter in this case.

            left_squared = ( left_coupling * left_coupling );

            return ( left_squared
                     / ( left_squared
                         + right_coupling * right_coupling ) );

          }
        else
          {

            if( is_verbose() )
              {

                std::cout
                << std::endl << "LHC-FASER::warning!"
                << " input_handler::quark_or_lepton_left_handedness"
                << " didn't recognize the sfermion PDG code \""
                << sfermion_code << "\". returning "
                << CppSLHA::CppSLHA_global::really_wrong_value_string;
                std::cout << std::endl;

              }

            return CppSLHA::CppSLHA_global::really_wrong_value;

          }

      }  // end of if it's a neutralino.

  }

}  // end of LHC_FASER namespace.

