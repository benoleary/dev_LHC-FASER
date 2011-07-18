/*
 * LHC-FASER_input_handling_stuff.hpp
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

#ifndef LHC_FASER_INPUT_HANDLING_STUFF_HPP_
#define LHC_FASER_INPUT_HANDLING_STUFF_HPP_

#include "CppSLHA/CppSLHA.hpp"
#include "LHC-FASER_global_stuff.hpp"
#include "LHC-FASER_derived_lepton_distributions.hpp"
#include "LHC-FASER_sparticle_decay_stuff.hpp"

namespace LHC_FASER
{

  /* this class holds the various gauge couplings which are re-calculated for
   * every new point. it is kept separate from the input_handler class for ease
   * of (abusing) const correctness.
   */
  class update_dependent_inputs : public readied_for_new_point
  {

  public:

    update_dependent_inputs( CppSLHA::SLHA_BLOCK const* given_SM_inputs,
              CppSLHA::particle_property_set const* const given_W_plus_pointer,
                   CppSLHA::particle_property_set const* const given_Z_pointer,
                             readier_for_new_point* const given_readier ) :
      readied_for_new_point( given_readier ),
      SM_inputs( given_SM_inputs ),
      W_plus_pointer( given_W_plus_pointer ),
      Z_pointer( given_Z_pointer )
    {

      // just an initialization list.

    }

    ~update_dependent_inputs()
    {

      // does nothing.

    }


    inline double
    get_electromagnetic_coupling()
    {

      if( needs_to_prepare_for_this_point() )
        {

          calculate_gauge_couplings();
          finish_preparing_for_this_point();

        }

      return electromagnetic_coupling;

    }

    inline double
    get_g_one()
    {

      if( needs_to_prepare_for_this_point() )
        {

          calculate_gauge_couplings();
          finish_preparing_for_this_point();

        }

      return g_one;

    }

    inline double
    get_g_two()
    {

      if( needs_to_prepare_for_this_point() )
        {

          calculate_gauge_couplings();
          finish_preparing_for_this_point();

        }

      return g_two;

    }

    inline double
    get_g_three()
    {

      if( needs_to_prepare_for_this_point() )
        {

          calculate_gauge_couplings();
          finish_preparing_for_this_point();

        }

      return g_three;

    }

    inline double
    get_weak_sine()
    {

      if( needs_to_prepare_for_this_point() )
        {

          calculate_gauge_couplings();
          finish_preparing_for_this_point();

        }

      return weak_sine;

    }

    inline double
    get_weak_cosine()
    {

      if( needs_to_prepare_for_this_point() )
        {

          calculate_gauge_couplings();
          finish_preparing_for_this_point();

        }

      return weak_cosine;

    }



  protected:

    CppSLHA::SLHA_BLOCK const* SM_inputs;
    CppSLHA::particle_property_set const* const W_plus_pointer;
    CppSLHA::particle_property_set const* const Z_pointer;

    // the gauge couplings:
    double electromagnetic_coupling;
    double g_one;
    // the gauge coupling for U(1) at m_Z.
    double g_two;
    // the gauge coupling for SU(2) at m_Z.
    double g_three;
    // the gauge coupling for SU(3) at m_Z.
    double weak_sine;
    // the sine of the weak mixing angle.
    double weak_cosine;
    // the cosine of the weak mixing angle.

    inline void
    calculate_gauge_couplings()
    {

      // we work out the gauge couplings from the SMINPUTS BLOCK:
      if( NULL != SM_inputs )
        {

          electromagnetic_coupling
          = sqrt( ( ( 4.0 * M_PI ) / SM_inputs->get_entry( 1 ) ) );
          weak_cosine = ( W_plus_pointer->get_absolute_mass()
                          / Z_pointer->get_absolute_mass() );
          weak_sine = sqrt( ( 1.0 - weak_cosine * weak_cosine ) );
          g_one = ( electromagnetic_coupling / weak_cosine );
          g_two = ( electromagnetic_coupling / weak_sine );
          g_three = sqrt( ( 4.0 * M_PI * SM_inputs->get_entry( 3 ) ) );

        }
      else
        {

          std::cout
          << std::endl
          << "LHC-FASER::warning! input_handler::input_handler() did not"
          <<" provide a valid pointer to SMINPUTS. gauge couplings are set to "
          << CppSLHA::CppSLHA_global::really_wrong_value_string;
          std::cout << std::endl;

          electromagnetic_coupling
          = CppSLHA::CppSLHA_global::really_wrong_value;
          weak_sine = CppSLHA::CppSLHA_global::really_wrong_value;
          weak_cosine = CppSLHA::CppSLHA_global::really_wrong_value;
          g_one = CppSLHA::CppSLHA_global::really_wrong_value;
          g_two = CppSLHA::CppSLHA_global::really_wrong_value;
          g_three = CppSLHA::CppSLHA_global::really_wrong_value;

        }

    }

  };


  // this class holds a set of pointers to various objects for ease of
  // reference to particles.
  class input_handler
  {

  public:

    input_handler( CppSLHA::CppSLHA0 const* const given_CppSLHA_pointer,
                CppSLHA::EW_scale_spectrum const* const given_spectrum_pointer,
                   std::string const given_path_to_grids,
                   readier_for_new_point* const given_readier )
    /* code in .cpp file. */;

    ~input_handler()
    /* code in .cpp file. */;


    void
    add_squark_as_colored_sparticle_and_in_production(
                     CppSLHA::particle_property_set const* const given_squark )
    /* this adds the given squark to the list of squarks & the list of colored
     * sparticles, & also adds its relevant production combinations (with
     * a gluino, with its antisquark, with itself) to the list.
     */
    /* code in .cpp file. */;


    inline bool
    is_verbose()
    const
    // if this is true, lots of warnings will be printed. if false, the code
    // should stay quiet without warning of the impending seg faults...
    {

      return verbosity_flag;

    }

    inline void
    set_verbosity( bool given_verbosity )
    {

      verbosity_flag = given_verbosity;

    }


    inline bool
    using_NLO()
    const
    {

      return using_NLO_flag;

    }

    inline void
    set_using_NLO( bool given_using_NLO_flag )
    {

      using_NLO_flag = given_using_NLO_flag;

    }


    inline std::string const*
    get_path_to_grids()
    const
    {

      return &path_to_grids;

    }

    inline std::string*
    get_path_to_cross_section_grids()
    {

      return &path_to_cross_section_grids;

    }

    inline std::string*
    get_path_to_kinematics_grids()
    {

      return &path_to_kinematics_grids;

    }

    inline std::string const*
    inspect_path_to_cross_section_grids()
    const
    {

      return &path_to_cross_section_grids;

    }

    inline std::string const*
    inspect_path_to_kinematics_grids()
    const
    {

      return &path_to_kinematics_grids;

    }


    inline CppSLHA::CppSLHA0 const*
    get_CppSLHA()
    const
    {

      return CppSLHA_pointer;

    }

    inline CppSLHA::EW_scale_spectrum const*
    get_spectrum()
    const
    {

      return spectrum_pointer;

    }


    inline double
    get_electromagnetic_coupling()
    const
    {

      return update_dependent_input_values->get_electromagnetic_coupling();

    }

    inline double
    get_g_one()
    const
    {

      return update_dependent_input_values->get_g_one();

    }

    inline double
    get_g_two()
    const
    {

      return update_dependent_input_values->get_g_two();

    }

    inline double
    get_g_three()
    const
    {

      return update_dependent_input_values->get_g_three();

    }

    inline double
    get_weak_sine()
    const
    {

      return update_dependent_input_values->get_weak_sine();

    }

    inline double
    get_weak_cosine()
    const
    {

      return update_dependent_input_values->get_weak_cosine();

    }


    inline readier_for_new_point*
    get_readier()
    const
    {

      return readier;

    }


    inline CppSLHA::particle_property_set const*
    get_down()
    const
    {

      return down_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_up()
    const
    {

      return up_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_strange()
    const
    {

      return strange_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_charm()
    const
    {

      return charm_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_bottom()
    const
    {

      return bottom_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_top()
    const
    {

      return top_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_electron()
    const
    {

      return electron_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_electron_neutrino()
    const
    {

      return electron_neutrino_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_muon()
    const
    {

      return muon_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_muon_neutrino()
    const
    {

      return muon_neutrino_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_tau_lepton()
    const
    {

      return tau_lepton_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_tau_neutrino()
    const
    {

      return tau_neutrino_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_gluon()
    const
    {

      return gluon_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_photon()
    const
    {

      return photon_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_Z()
    const
    {

      return Z_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_W_plus()
    const
    {

      return W_plus_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_light_neutral_EWSB_scalar()
    const
    {

      return light_neutral_EWSB_scalar_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_heavy_neutral_EWSB_scalar()
    const
    {

      return heavy_neutral_EWSB_scalar_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_neutral_EWSB_pseudoscalar()
    const
    {

      return neutral_EWSB_pseudoscalar_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_charged_EWSB_scalar()
    const
    {

      return charged_EWSB_scalar_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_gluino()
    const
    {

      return gluino_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_sdown_L()
    const
    {

      return sdown_L_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_sup_L()
    const
    {

      return sup_L_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_sstrange_L()
    const
    {

      return sstrange_L_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_scharm_L()
    const
    {

      return scharm_L_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_sbottom_one()
    const
    {

      return sbottom_one_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_stop_one()
    const
    {

      return stop_one_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_sdown_R()
    const
    {

      return sdown_R_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_sup_R()
    const
    {

      return sup_R_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_sstrange_R()
    const
    {

      return sstrange_R_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_scharm_R()
    const
    {

      return scharm_R_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_sbottom_two()
    const
    {

      return sbottom_two_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_stop_two()
    const
    {

      return stop_two_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_selectron_L()
    const
    {

      return selectron_L_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_electron_sneutrino_L()
    const
    {

      return electron_sneutrino_L_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_smuon_L()
    const
    {

      return smuon_L_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_muon_sneutrino_L()
    const
    {

      return muon_sneutrino_L_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_stau_one()
    const
    {

      return stau_one_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_tau_sneutrino_L()
    const
    {

      return tau_sneutrino_L_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_selectron_R()
    const
    {

      return selectron_R_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_electron_sneutrino_R()
    const
    {

      return electron_sneutrino_R_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_smuon_R()
    const
    {

      return smuon_R_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_muon_sneutrino_R()
    const
    {

      return muon_sneutrino_R_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_stau_two()
    const
    {

      return stau_two_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_tau_sneutrino_R()
    const
    {

      return tau_sneutrino_R_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_neutralino_one()
    const
    {

      return neutralino_one_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_neutralino_two()
    const
    {

      return neutralino_two_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_neutralino_three()
    const
    {

      return neutralino_three_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_neutralino_four()
    const
    {

      return neutralino_four_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_chargino_one()
    const
    {

      return chargino_one_pointer;

    }

    inline CppSLHA::particle_property_set const*
    get_chargino_two()
    const
    {

      return chargino_two_pointer;

    }


    bool
    is_in( int const given_PDG_code,
           std::vector< CppSLHA::particle_property_set const* > const* const
           given_collection )
    const
    /* this checks to see if the particle with the absolute value of the given
     * PDG particle code is in the given collection. (hence it checks to see if
     * the given particle or its antiparticle are in the collection.)
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    inline direct_decay_checker*
    get_direct_decay_checker(
                     CppSLHA::particle_property_set const* const given_decayer,
              std::vector< CppSLHA::particle_property_set const* > const* const
                              given_products )
    const
    {

       return decay_checkers->get_decay_checker( given_decayer,
                                                 given_products );

    }

    inline exclusive_BR_calculator*
    get_exclusive_BR_calculator(
                     CppSLHA::particle_property_set const* const given_decayer,
                     CppSLHA::particle_property_set const* const given_product,
                                 bool const product_is_not_antiparticle_flag,
                                 std::list< int > const* const given_excluded )
    const
    {

      return exclusive_BRs->get_BR_calculator( given_decayer,
                                               given_product,
                                              product_is_not_antiparticle_flag,
                                               given_excluded );

    }

    inline colored_cascade_set*
    get_colored_cascade_set(
                    CppSLHA::particle_property_set const* const given_decayer )
    const
    {

      return colored_cascades->get_colored_cascade_set( given_decayer );

    }


    std::list< int > const*
    get_empty_list()
    const
    {

      return &empty_list;

    }


    std::list< int > const*
    get_jets()
    const
    {

      return &jets;

    }

    std::list< int > const*
    get_jets5()
    const
    {

      return &jets5;

    }

    std::list< int > const*
    get_jets4()
    const
    {

      return &jets4;

    }

    std::list< int > const*
    get_not_in_jets5()
    const
    {

      return &not_in_jets5;

    }

    std::list< int > const*
    get_not_in_jets4()
    const
    {

      return &not_in_jets4;

    }


    inline std::vector< CppSLHA::particle_property_set const* > const*
    get_colored_sparticles()
    const
    {

      return &colored_sparticles;

    }

    inline std::vector< CppSLHA::particle_property_set const* > const*
    get_squarks()
    const
    {

      return &squarks;

    }

    inline std::vector< CppSLHA::particle_property_set const* > const*
    get_squarks5()
    const
    {

      return &squarks5;

    }

    inline std::vector< CppSLHA::particle_property_set const* > const*
    get_squarks4()
    const
    {

      return &squarks4;

    }

    inline std::vector< CppSLHA::particle_property_set const* > const*
    get_sdown_types()
    const
    {

      return &sdown_types;

    }

    inline std::vector< CppSLHA::particle_property_set const* > const*
    get_sup_types()
    const
    {

      return &sup_types;

    }

    inline std::vector< signed_particle_shortcut_pair* > const*
    get_sparticle_production_combinations()
    const
    {

      return &sparticle_production_combinations;

    }


    inline std::vector< CppSLHA::particle_property_set const* > const*
    get_electroweakinos()
    const
    {

      return &electroweakinos;

    }

    inline std::vector< CppSLHA::particle_property_set const* > const*
    get_unstable_electroweakinos()
    const
    {

      return &unstable_electroweakinos;

    }

    inline std::vector< CppSLHA::particle_property_set const* > const*
    get_neutralinos()
    const
    {

      return &neutralinos;

    }

    inline std::vector< CppSLHA::particle_property_set const* > const*
    get_unstable_neutralinos()
    const
    {

      return &unstable_neutralinos;

    }

    inline std::vector< CppSLHA::particle_property_set const* > const*
    get_charginos()
    const
    {

      return &charginos;

    }

    inline std::vector< CppSLHA::particle_property_set const* > const*
    get_charged_sleptons()
    const
    {

      return &charged_sleptons;

    }

    inline std::vector< CppSLHA::particle_property_set const* > const*
    get_sneutrinos()
    const
    {

      return &sneutrinos;

    }

    inline std::vector< CppSLHA::particle_property_set const* > const*
    get_neutral_EWSB_scalars_and_pseudoscalars()
    const
    {

      return &neutral_EWSB_scalars_and_pseudoscalars;

    }

    inline std::vector< CppSLHA::particle_property_set const* > const*
    get_EW_veved_and_vector_bosons()
    const
    {

      return &EW_veved_and_vector_bosons;

    }


    hard_muon_from_tau const*
    get_hard_muon_from_tau()
    const
    {

      return &hard_muon_from_tau_function;

    }

    soft_muon_from_tau const*
    get_soft_muon_from_tau()
    const
    {

      return &soft_muon_from_tau_function;

    }

    hard_pion_from_tau const*
    get_hard_pion_from_tau()
    const
    {

      return &hard_pion_from_tau_function;

    }

    soft_pion_from_tau const*
    get_soft_pion_from_tau()
    const
    {

      return &soft_pion_from_tau_function;

    }


    inline double
    get_average_squarks4_mass()
    const
    {

      return ( 0.125 * ( get_sdown_L()->get_absolute_mass()
                         + get_sdown_R()->get_absolute_mass()
                         + get_sup_L()->get_absolute_mass()
                         + get_sup_R()->get_absolute_mass()
                         + get_sstrange_L()->get_absolute_mass()
                         + get_sstrange_R()->get_absolute_mass()
                         + get_scharm_L()->get_absolute_mass()
                         + get_scharm_R()->get_absolute_mass() ) );

    }


    double
    quark_or_lepton_left_handedness( int const sfermion_code,
                                     int const EWino_code )
    const
    /* this returns the square of the coupling associated with the left-handed
     * projection operator of the Feynman rule for the vertex of the given
     * sfermion with the given electroweakino & the appropriate quark, divided
     * by the sum of the squares of the couplings for the left- & right-handed
     * projection operators, based on Rosiek's Feynman rules. this has to be
     * changed if general squark flavor mixing is incorporated. also, it
     * ignores the possibility of right-handed neutrinos.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:

    bool verbosity_flag;
    bool using_NLO_flag;

    readier_for_new_point* const readier;

    CppSLHA::CppSLHA0 const* const CppSLHA_pointer;
    CppSLHA::EW_scale_spectrum const* const spectrum_pointer;

    std::string path_to_grids;
    std::string path_to_cross_section_grids;
    std::string path_to_kinematics_grids;

    CppSLHA::particle_property_set const* const down_pointer;
    CppSLHA::particle_property_set const* const up_pointer;
    CppSLHA::particle_property_set const* const strange_pointer;
    CppSLHA::particle_property_set const* const charm_pointer;
    CppSLHA::particle_property_set const* const bottom_pointer;
    CppSLHA::particle_property_set const* const top_pointer;
    CppSLHA::particle_property_set const* const electron_pointer;
    CppSLHA::particle_property_set const* const electron_neutrino_pointer;
    CppSLHA::particle_property_set const* const muon_pointer;
    CppSLHA::particle_property_set const* const muon_neutrino_pointer;
    CppSLHA::particle_property_set const* const tau_lepton_pointer;
    CppSLHA::particle_property_set const* const tau_neutrino_pointer;
    CppSLHA::particle_property_set const* const gluon_pointer;
    CppSLHA::particle_property_set const* const photon_pointer;
    CppSLHA::particle_property_set const* const Z_pointer;
    CppSLHA::particle_property_set const* const W_plus_pointer;

    CppSLHA::particle_property_set const* const
    light_neutral_EWSB_scalar_pointer;
    CppSLHA::particle_property_set const* const
    heavy_neutral_EWSB_scalar_pointer;
    CppSLHA::particle_property_set const* const
    neutral_EWSB_pseudoscalar_pointer;
    CppSLHA::particle_property_set const* const charged_EWSB_scalar_pointer;
    CppSLHA::particle_property_set const* const gluino_pointer;
    CppSLHA::particle_property_set const* const sdown_L_pointer;
    CppSLHA::particle_property_set const* const sup_L_pointer;
    CppSLHA::particle_property_set const* const sstrange_L_pointer;
    CppSLHA::particle_property_set const* const scharm_L_pointer;
    CppSLHA::particle_property_set const* const sbottom_one_pointer;
    CppSLHA::particle_property_set const* const stop_one_pointer;
    CppSLHA::particle_property_set const* const sdown_R_pointer;
    CppSLHA::particle_property_set const* const sup_R_pointer;
    CppSLHA::particle_property_set const* const sstrange_R_pointer;
    CppSLHA::particle_property_set const* const scharm_R_pointer;
    CppSLHA::particle_property_set const* const sbottom_two_pointer;
    CppSLHA::particle_property_set const* const stop_two_pointer;
    CppSLHA::particle_property_set const* const selectron_L_pointer;
    CppSLHA::particle_property_set const* const electron_sneutrino_L_pointer;
    CppSLHA::particle_property_set const* const smuon_L_pointer;
    CppSLHA::particle_property_set const* const muon_sneutrino_L_pointer;
    CppSLHA::particle_property_set const* const stau_one_pointer;
    CppSLHA::particle_property_set const* const tau_sneutrino_L_pointer;
    CppSLHA::particle_property_set const* const selectron_R_pointer;
    CppSLHA::particle_property_set const* const electron_sneutrino_R_pointer;
    CppSLHA::particle_property_set const* const smuon_R_pointer;
    CppSLHA::particle_property_set const* const muon_sneutrino_R_pointer;
    CppSLHA::particle_property_set const* const stau_two_pointer;
    CppSLHA::particle_property_set const* const tau_sneutrino_R_pointer;
    CppSLHA::particle_property_set const* const neutralino_one_pointer;
    CppSLHA::particle_property_set const* const neutralino_two_pointer;
    CppSLHA::particle_property_set const* const neutralino_three_pointer;
    CppSLHA::particle_property_set const* const neutralino_four_pointer;
    CppSLHA::particle_property_set const* const chargino_one_pointer;
    CppSLHA::particle_property_set const* const chargino_two_pointer;


    // a few useful SLHA BLOCK shortcuts:
    CppSLHA::SLHA_BLOCK const* NMIX_pointer;
    CppSLHA::SLHA_BLOCK const* UMIX_pointer;
    CppSLHA::SLHA_BLOCK const* VMIX_pointer;
    CppSLHA::SLHA_BLOCK const* SBOTMIX_pointer;
    CppSLHA::SLHA_BLOCK const* STOPMIX_pointer;
    CppSLHA::SLHA_BLOCK const* STAUMIX_pointer;
    CppSLHA::SLHA_BLOCK const* YD_pointer;
    CppSLHA::SLHA_BLOCK const* YU_pointer;
    CppSLHA::SLHA_BLOCK const* YE_pointer;

    update_dependent_inputs* update_dependent_input_values;


    std::list< int > empty_list;
    // empty_list is there for when we don't care about excluding particles.

    /* the following std::list< int > are lists of PDG codes for particles:
     * jets refers to all the particles in the spectrum which return true for
     * counts_as_jet(), as well as their antiparticles, unless they are
     * self-conjugate.
     * not_in_jets5 is just top & antitop, & not_in_jets4 is not_in_jets5 plus
     * bottom & antibottom.
     * jets5 is jets minus everything in not_int_jets5, likewise jets4 is jets
     * minus everything in not_int_jets4.
     * other exotic heavy quarks should either be labelled as not counting as
     * a jet in the spectrum, or put in not_in_jets5 in the same manner as top.
     * hence, without additional light jet particles, jets5 is
     * { gluon + the 5 light quarks (down, up, strange, charm, bottom)
     * + the 5 light antiquarks }, & jets4 is jets5 with bottom & antibottom
     * removed.
     */
    std::list< int > jets;
    std::list< int > jets5;
    std::list< int > jets4;
    std::list< int > not_in_jets5;
    std::list< int > not_in_jets4;

    std::vector< CppSLHA::particle_property_set const* >
    colored_sparticles;

    // squarks4 refers to { sdown, sup, sstrange, scharm }, while squarks5
    // is squarks4 plus sbottom, & squarks is squarks5 plus stop.
    std::vector< CppSLHA::particle_property_set const* > squarks;
    std::vector< CppSLHA::particle_property_set const* > squarks5;
    std::vector< CppSLHA::particle_property_set const* > squarks4;

    std::vector< CppSLHA::particle_property_set const* > sdown_types;
    std::vector< CppSLHA::particle_property_set const* > sup_types;

    std::vector< signed_particle_shortcut_pair* >
    sparticle_production_combinations;

    std::vector< CppSLHA::particle_property_set const* > electroweakinos;
    std::vector< CppSLHA::particle_property_set const* >
    unstable_electroweakinos;
    std::vector< CppSLHA::particle_property_set const* > neutralinos;
    std::vector< CppSLHA::particle_property_set const* > unstable_neutralinos;
    std::vector< CppSLHA::particle_property_set const* > charginos;
    std::vector< CppSLHA::particle_property_set const* > charged_sleptons;
    std::vector< CppSLHA::particle_property_set const* > sneutrinos;

    std::vector< CppSLHA::particle_property_set const* >
    neutral_EWSB_scalars_and_pseudoscalars;
    std::vector< CppSLHA::particle_property_set const* >
    EW_veved_and_vector_bosons;

    decay_checker_handler* decay_checkers;
    exclusive_BR_handler* exclusive_BRs;
    colored_cascade_handler* colored_cascades;

    hard_muon_from_tau hard_muon_from_tau_function;
    soft_muon_from_tau soft_muon_from_tau_function;
    hard_pion_from_tau hard_pion_from_tau_function;
    soft_pion_from_tau soft_pion_from_tau_function;

  };  // end of input_handler class.




  // inline functions:


  inline CppSLHA::particle_property_set const*
  signed_particle_shortcut_pair::get_first_pointer()
  const
  {

    return first_particle_pointer;

  }

  inline CppSLHA::particle_property_set const*
  signed_particle_shortcut_pair::get_second_pointer()
  const
  {

    return second_particle_pointer;

  }


  inline bool
  signed_particle_shortcut_pair::first_is_not_antiparticle()
  const
  {

    return first_is_not_antiparticle_flag;

  }

  inline bool
  signed_particle_shortcut_pair::second_is_not_antiparticle()
  const
  {

    return second_is_not_antiparticle_flag;

  }

}

#endif /* LHC_FASER_INPUT_HANDLING_STUFF_HPP_ */
