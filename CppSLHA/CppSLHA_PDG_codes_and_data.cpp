/*
 * CppSLHA_particle_codes.cpp
 *
 *  Created on: 19 Jan 2011
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2011 Ben O'Leary
 *
 *      This file is part of CppSLHA.
 *
 *      CppSLHA is free software: you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation, either version 3 of the License, or
 *      (at your option) any later version.
 *
 *      CppSLHA is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *
 *      You should have received a copy of the GNU General Public License
 *      along with CppSLHA.  If not, see <http://www.gnu.org/licenses/>.
 *
 *      The GNU General Public License should be in GNU_public_license.txt
 *      the files of CppSLHA are:
 *      CppSLHA.hpp
 *      CppSLHA.cpp
 *      CppSLHA_BLOCK.hpp
 *      CppSLHA_BLOCK.cpp
 *      CppSLHA_EW_scale_spectrum.hpp
 *      CppSLHA_EW_scale_spectrum.cpp
 *      CppSLHA_file_reading_stuff.hpp
 *      CppSLHA_file_reading_stuff.cpp
 *      CppSLHA_global_stuff.hpp
 *      CppSLHA_global_stuff.cpp
 *      CppSLHA_PDG_codes_and_data.hpp
 *      CppSLHA_PDG_codes_and_data.cpp
 *      CppSLHA_particle_data.hpp
 *      CppSLHA_particle_data.cpp
 *      CppSLHA_waiting_on_subprocess_executor.hpp
 *      CppSLHA_waiting_on_subprocess_executor.cpp
 *      and ./README.CppSLHA.txt which describes the package.
 *
 */

#include "CppSLHA_PDG_codes_and_data.hpp"

namespace CppSLHA
{

  /* here is where all the defining for the static consts happens, in this
   * source file because at least 1 compiler tries to include the definitions
   * multiple times if they are put outside the class but in the header file,
   * regardless of whether they are within the #ifndef bit.
   */

  // Standard Model particle codes:
  int const PDG_code::down = 1;
  int const PDG_code::up = 2;
  int const PDG_code::strange = 3;
  int const PDG_code::charm = 4;
  int const PDG_code::bottom = 5;
  int const PDG_code::top = 6;
  int const PDG_code::electron = 11;
  int const PDG_code::electron_neutrino = 12;
  int const PDG_code::muon = 13;
  int const PDG_code::muon_neutrino = 14;
  int const PDG_code::tau_lepton = 15;
  int const PDG_code::tau_neutrino = 16;
  int const PDG_code::gluon = 21;
  int const PDG_code::photon = 22;
  int const PDG_code::Z = 23;
  int const PDG_code::W_plus = 24;
  int const PDG_code::EWSB_scalar = 25;

  // 2 Higgs Doublet Model particle codes:
  int const PDG_code::light_neutral_EWSB_scalar = 25;
  int const PDG_code::heavy_neutral_EWSB_scalar = 35;
  int const PDG_code::neutral_EWSB_pseudoscalar = 36;
  int const PDG_code::charged_EWSB_scalar = 37;

  // almost-SM particle codes:
  int const PDG_code::graviton = 39;

  // Minimal Supersymmetric Standard Model particle codes:
  int const PDG_code::gluino = 1000021;
  int const PDG_code::sdown_L = 1000001;
  int const PDG_code::sup_L = 1000002;
  int const PDG_code::sstrange_L = 1000003;
  int const PDG_code::scharm_L = 1000004;
  int const PDG_code::sbottom_one = 1000005;
  int const PDG_code::stop_one = 1000006;
  int const PDG_code::sdown_R = 2000001;
  int const PDG_code::sup_R = 2000002;
  int const PDG_code::sstrange_R = 2000003;
  int const PDG_code::scharm_R = 2000004;
  int const PDG_code::sbottom_two = 2000005;
  int const PDG_code::stop_two = 2000006;
  int const PDG_code::selectron_L = 1000011;
  int const PDG_code::smuon_L = 1000013;
  int const PDG_code::stau_one = 1000015;
  int const PDG_code::selectron_R = 2000011;
  int const PDG_code::smuon_R = 2000013;
  int const PDG_code::stau_two = 2000015;
  int const PDG_code::electron_sneutrino_L = 1000012;
  int const PDG_code::muon_sneutrino_L = 1000014;
  int const PDG_code::tau_sneutrino_L = 1000016;
  int const PDG_code::neutralino_one = 1000022;
  int const PDG_code::neutralino_two = 1000023;
  int const PDG_code::neutralino_three = 1000025;
  int const PDG_code::neutralino_four = 1000035;
  int const PDG_code::chargino_one = 1000024;
  int const PDG_code::chargino_two = 1000037;

  // flavor-violation codes for mass-ordered sfermions:
  int const PDG_code::sdown_one = 1000001;
  int const PDG_code::sdown_two = 1000003;
  int const PDG_code::sdown_three = 1000005;
  int const PDG_code::sdown_four = 2000001;
  int const PDG_code::sdown_five = 2000003;
  int const PDG_code::sdown_six = 2000005;
  int const PDG_code::sup_one = 1000002;
  int const PDG_code::sup_two = 1000004;
  int const PDG_code::sup_three = 1000006;
  int const PDG_code::sup_four = 2000002;
  int const PDG_code::sup_five = 2000004;
  int const PDG_code::sup_six = 2000006;
  int const PDG_code::selectron_one = 1000011;
  int const PDG_code::selectron_two = 1000013;
  int const PDG_code::selectron_three = 1000015;
  int const PDG_code::selectron_four = 2000011;
  int const PDG_code::selectron_five = 2000013;
  int const PDG_code::selectron_six = 2000015;
  int const PDG_code::sneutrino_one = 1000012;
  int const PDG_code::sneutrino_two = 1000014;
  int const PDG_code::sneutrino_three = 1000016;

  // sneutrino CP-eigenstate codes:
  int const PDG_code::sneutrino_scalar_one = 1000012;
  int const PDG_code::sneutrino_scalar_two = 1000014;
  int const PDG_code::sneutrino_scalar_three = 1000016;
  int const PDG_code::sneutrino_pseudoscalar_one = 1000017;
  int const PDG_code::sneutrino_pseudoscalar_two = 1000018;
  int const PDG_code::sneutrino_pseudoscalar_three = 1000019;

  // Next-to-Minimal Supersymmetric Standard Model particle codes:
  int const PDG_code::singlet_scalar = 45;
  int const PDG_code::singlet_pseudoscalar = 46;
  int const PDG_code::neutralino_five = 1000045;

  /* extra sneutrino codes for 3 generations of right-handed neutrinos (these
   * are NOT official, they are just extrapolated from the pattern of the other
   * sfermions), including flavor-violation codes:
   */
  int const PDG_code::electron_sneutrino_R = 2000012;
  int const PDG_code::muon_sneutrino_R = 2000014;
  int const PDG_code::tau_sneutrino_R = 2000016;
  int const PDG_code::sneutrino_four = 2000012;
  int const PDG_code::sneutrino_five = 2000014;
  int const PDG_code::sneutrino_six = 2000016;
  int const PDG_code::sneutrino_scalar_four = 2000012;
  int const PDG_code::sneutrino_scalar_five = 2000014;
  int const PDG_code::sneutrino_scalar_six = 2000016;
  int const PDG_code::sneutrino_pseudoscalar_four = 2000017;
  int const PDG_code::sneutrino_pseudoscalar_five = 2000018;
  int const PDG_code::sneutrino_pseudoscalar_six = 2000019;

  // supergravity particlecodes:
  int const PDG_code::gravitino = 1000039;

  // Standard Model particle masses:
  double const PDG_data::down_mass = 0.00504;
  double const PDG_data::up_mass = 0.00255;
  double const PDG_data::strange_mass = 0.105;
  double const PDG_data::charm_mass = 1.27;
  double const PDG_data::bottom_mass = 4.20;
  double const PDG_data::top_mass = 171.3;
  double const PDG_data::electron_mass = 0.000510998910;
  double const PDG_data::electron_neutrino_mass = 0.0;
  double const PDG_data::muon_mass = 0.105658367;
  double const PDG_data::muon_neutrino_mass = 0.0;
  double const PDG_data::tau_lepton_mass = 1.77684;
  double const PDG_data::tau_neutrino_mass = 0.0;
  double const PDG_data::gluon_mass = 0.0;
  double const PDG_data::photon_mass = 0.0;
  double const PDG_data::Z_mass = 91.1876;
  double const PDG_data::W_plus_mass = 80.398;

  double const PDG_data::Z_decay_width = 2.4952;
  double const PDG_data::Z_to_electron_antielectron_BR = 0.03363;
  double const PDG_data::Z_to_muon_antimuon_BR = 0.03363;
  double const PDG_data::Z_to_tau_lepton_antilepton_BR = 0.03370;
  double const PDG_data::Z_to_invisible_BR = 0.20;

  // this code assumes that the invisible decays of Z bosons are equally
  // divided between the 3 flavors of neutrino:
  double const
  PDG_data::Z_to_electron_neutrino_antineutrino_BR
    = ( PDG_data::Z_to_invisible_BR / 3.0 );
  double const
  PDG_data::Z_to_muon_neutrino_antineutrino_BR
    = PDG_data::Z_to_electron_neutrino_antineutrino_BR;
  double const
  PDG_data::Z_to_tau_neutrino_antineutrino_BR
    = PDG_data::Z_to_electron_neutrino_antineutrino_BR;

  // this code assumes that the rest of the hadronic decay width of the Z boson
  // is equally divided between down, up and strange:
  double const PDG_data::Z_to_charm_anticharm_BR = 0.1203;
  double const PDG_data::Z_to_bottom_antibottom_BR = 0.1512;
  double const
  PDG_data::Z_to_down_antidown_BR
    = ( ( 0.6991 - PDG_data::Z_to_charm_anticharm_BR
                 - PDG_data::Z_to_bottom_antibottom_BR ) / 3.0 );
  double const
  PDG_data::Z_to_up_antiup_BR = PDG_data::Z_to_down_antidown_BR;
  double const
  PDG_data::Z_to_strange_antistrange_BR = PDG_data::Z_to_down_antidown_BR;

  double const PDG_data::W_plus_decay_width = 2.141;
  double const PDG_data::W_plus_to_neutrino_antielectron_BR = 0.1075;
  double const PDG_data::W_plus_to_neutrino_antimuon_BR = 0.1057;
  double const PDG_data::W_plus_to_neutrino_tau_antilepton_BR = 0.1125;
  double const PDG_data::W_plus_to_hadrons_BR = 0.6760;

  double const PDG_data::V_cd = 0.23;
  double const PDG_data::V_cs = 1.04;
  double const PDG_data::V_cd_squared = ( PDG_data::V_cd * PDG_data::V_cd );
  double const PDG_data::V_cs_squared = ( PDG_data::V_cs * PDG_data::V_cs );
  double const
  PDG_data::one_over_V_cd_squared_plus_V_cs_squared
    = ( 1.0 / ( PDG_data::V_cd_squared + PDG_data::V_cs_squared ) );
  double const
  PDG_data::V_cd_squared_fraction
    = ( PDG_data::V_cd_squared
        * PDG_data::one_over_V_cd_squared_plus_V_cs_squared );
  double const
  PDG_data::V_cs_squared_fraction
    = ( PDG_data::V_cs_squared
        * PDG_data::one_over_V_cd_squared_plus_V_cs_squared );
  double const PDG_data::V_ud = 0.97418;
  double const PDG_data::V_us = 0.2255;
  double const PDG_data::V_ud_squared = ( PDG_data::V_ud * PDG_data::V_ud );
  double const PDG_data::V_us_squared = ( PDG_data::V_us * PDG_data::V_us );
  double const
  PDG_data::one_over_V_ud_squared_plus_V_us_squared
    = ( 1.0 / ( PDG_data::V_ud_squared + PDG_data::V_us_squared ) );
  double const
  PDG_data::V_ud_squared_fraction
    = ( PDG_data::V_ud_squared
        * PDG_data::one_over_V_ud_squared_plus_V_us_squared );
  double const
  PDG_data::V_us_squared_fraction
    = ( PDG_data::V_us_squared
        * PDG_data::one_over_V_ud_squared_plus_V_us_squared );

  /* this code assumes that the BRs of the W^+ into charm + antidown &
   * charm + antistrange account for all the BR of the charm + X,
   * in the ratio ( |V_cd|^2 ) to ( |V_cs|^2 ):
   */
  double const PDG_data::W_plus_to_charm_X_BR = 0.334;
  double const
  PDG_data::W_plus_to_charm_antidown_BR
    = ( PDG_data::V_cd_squared_fraction * PDG_data::W_plus_to_charm_X_BR );
  double const PDG_data::W_plus_to_charm_antistrange_BR
    = ( PDG_data::V_cs_squared_fraction * PDG_data::W_plus_to_charm_X_BR );

  /* this code assumes that the rest of the hadronic decay width of the W is
   * divided between up + antidown and up + antistrange
   * in the ratio ( |V_ud|^2 ) to ( |V_us|^2 ):
   */
  double const
  PDG_data::W_plus_to_charmless_X_BR
    = ( PDG_data::W_plus_to_hadrons_BR - PDG_data::W_plus_to_charm_X_BR );
  double const
  PDG_data::W_plus_to_up_antidown_BR
    = ( PDG_data::V_ud_squared_fraction * PDG_data::W_plus_to_charmless_X_BR );
  double const
  PDG_data::W_plus_to_up_antistrange_BR
    = ( PDG_data::V_us_squared_fraction * PDG_data::W_plus_to_charmless_X_BR );

  /* currently tau leptons are treated as stable by the constructor, but if
   * they are to be implemented, these values should be used.
   * all these values were taken from the PDG on 2009-11-10.
   */
  double const PDG_data::tau_lepton_to_neutrinos_electron_BR = 0.1785;
  double const PDG_data::tau_lepton_to_neutrinos_muon_BR = 0.1736;
  /* this code assumes that the rest of the decay width of the tau lepton is
   * divided between down + antiup and strange + antiup in the ratio
   * ( |V_ud|^2 ) to ( |V_us|^2 ) as in the case of the decays of the W^+:
   */
  double const
  PDG_data::tau_lepton_to_neutrino_hadron_BR
    = ( 1.0 - PDG_data::tau_lepton_to_neutrinos_electron_BR
            - PDG_data::tau_lepton_to_neutrinos_muon_BR );
  double const
  PDG_data::tau_lepton_to_neutrino_down_antiup_BR
    = ( PDG_data::V_ud_squared_fraction
        * PDG_data::tau_lepton_to_neutrino_hadron_BR );
  double const PDG_data::tau_lepton_to_neutrino_strange_antiup_BR
    = ( PDG_data::V_us_squared_fraction
        * PDG_data::tau_lepton_to_neutrino_hadron_BR );

}  // end of CppSLHA namespace.
