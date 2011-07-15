/*
 * CppSLHA_particle_codes.hpp
 *
 *  Created on: 15 Dec 2010
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2010 Ben O'Leary
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

/*********************************************************************
 * PLEASE NOTE that I have just used the PDG particle codes as they  *
 * are given in the SLHA2 paper:                                     *
 * B. C. Allanach et al., Comp.Phys.Commun.180:8-25,2009,            *
 * http://arxiv.org/abs/0801.0045 online.                            *
 * there the authors have given the 2 extra scalar boson states      *
 * numbers 45 & 46 (neutral scalar & pseudoscalar), & the 5th        *
 * neutralino number 1000045.                                        *
 * please also note that there are separate classes for the          *
 * flavor violation sparticle numbering scheme, called               *
 * EW_scale_FVMSSM_spectrum & EW_scale_FVNMSSM_spectrum which have   *
 * the MSSM & NMSSM particle content respectively, but relabelled    *
 * sfermions.                                                        *
 *********************************************************************/

#ifndef CPPSLHA_PDG_CODES_AND_DATA_HPP_
#define CPPSLHA_PDG_CODES_AND_DATA_HPP_

namespace CppSLHA
{

  class PDG_code
  {

  public:

    // Standard Model particle codes:
    static int const down;
    static int const up;
    static int const strange;
    static int const charm;
    static int const bottom;
    static int const top;
    static int const electron;
    static int const electron_neutrino;
    static int const muon;
    static int const muon_neutrino;
    static int const tau_lepton;
    static int const tau_neutrino;
    static int const gluon;
    static int const photon;
    static int const Z;
    static int const W_plus;
    static int const EWSB_scalar;

    // 2 Higgs Doublet Model particle codes:
    static int const light_neutral_EWSB_scalar;
    static int const heavy_neutral_EWSB_scalar;
    static int const neutral_EWSB_pseudoscalar;
    static int const charged_EWSB_scalar;

    // almost-SM particle codes:
    static int const graviton;

    // Minimal Supersymmetric Standard Model particle codes:
    static int const gluino;
    static int const sdown_L;
    static int const sup_L;
    static int const sstrange_L;
    static int const scharm_L;
    static int const sbottom_one;
    static int const stop_one;
    static int const sdown_R;
    static int const sup_R;
    static int const sstrange_R;
    static int const scharm_R;
    static int const sbottom_two;
    static int const stop_two;
    static int const selectron_L;
    static int const smuon_L;
    static int const stau_one;
    static int const selectron_R;
    static int const smuon_R;
    static int const stau_two;
    static int const electron_sneutrino_L;
    static int const muon_sneutrino_L;
    static int const tau_sneutrino_L;
    static int const neutralino_one;
    static int const neutralino_two;
    static int const neutralino_three;
    static int const neutralino_four;
    static int const chargino_one;
    static int const chargino_two;

    // flavor-violation codes for mass-ordered sfermions:
    static int const sdown_one;
    static int const sdown_two;
    static int const sdown_three;
    static int const sdown_four;
    static int const sdown_five;
    static int const sdown_six;
    static int const sup_one;
    static int const sup_two;
    static int const sup_three;
    static int const sup_four;
    static int const sup_five;
    static int const sup_six;
    static int const selectron_one;
    static int const selectron_two;
    static int const selectron_three;
    static int const selectron_four;
    static int const selectron_five;
    static int const selectron_six;
    static int const sneutrino_one;
    static int const sneutrino_two;
    static int const sneutrino_three;

    // sneutrino CP-eigenstate codes:
    static int const sneutrino_scalar_one;
    static int const sneutrino_scalar_two;
    static int const sneutrino_scalar_three;
    static int const sneutrino_pseudoscalar_one;
    static int const sneutrino_pseudoscalar_two;
    static int const sneutrino_pseudoscalar_three;

    // Next-to-Minimal Supersymmetric Standard Model particle codes:
    static int const singlet_scalar;
    static int const singlet_pseudoscalar;
    static int const neutralino_five;

    /* extra sneutrino codes for 3 generations of right-handed neutrinos (these
     * are NOT official, they are just extrapolated from the pattern of the
     * other sfermions), including flavor-violation codes:
     */
    static int const electron_sneutrino_R;
    static int const muon_sneutrino_R;
    static int const tau_sneutrino_R;
    static int const sneutrino_four;
    static int const sneutrino_five;
    static int const sneutrino_six;
    static int const sneutrino_scalar_four;
    static int const sneutrino_scalar_five;
    static int const sneutrino_scalar_six;
    static int const sneutrino_pseudoscalar_four;
    static int const sneutrino_pseudoscalar_five;
    static int const sneutrino_pseudoscalar_six;

    // supergravity particlecodes:
    static int const gravitino;

  };  // end of PDG_code class.


  class PDG_data
   // all data here taken from the PDG on 2009-11-08.
  {

  public:

    // Standard Model particle masses:
    static double const down_mass;
    static double const up_mass;
    static double const strange_mass;
    static double const charm_mass;
    static double const bottom_mass;
    static double const top_mass;
    static double const electron_mass;
    static double const electron_neutrino_mass;
    static double const muon_mass;
    static double const muon_neutrino_mass;
    static double const tau_lepton_mass;
    static double const tau_neutrino_mass;
    static double const gluon_mass;
    static double const photon_mass;
    static double const Z_mass;
    static double const W_plus_mass;

    // some SM particle decay widths & branching ratios:
    static double const Z_decay_width;
    static double const Z_to_electron_antielectron_BR;
    static double const Z_to_muon_antimuon_BR;
    static double const Z_to_tau_lepton_antilepton_BR;
    static double const Z_to_invisible_BR;
    static double const Z_to_electron_neutrino_antineutrino_BR;
    static double const Z_to_muon_neutrino_antineutrino_BR;
    static double const Z_to_tau_neutrino_antineutrino_BR;
    static double const Z_to_charm_anticharm_BR;
    static double const Z_to_bottom_antibottom_BR;
    static double const Z_to_down_antidown_BR;
    static double const Z_to_up_antiup_BR;
    static double const Z_to_strange_antistrange_BR;

    static double const W_plus_decay_width;
    static double const W_plus_to_neutrino_antielectron_BR;
    static double const W_plus_to_neutrino_antimuon_BR;
    static double const W_plus_to_neutrino_tau_antilepton_BR;
    static double const W_plus_to_hadrons_BR;
    static double const W_plus_to_charm_X_BR;
    static double const W_plus_to_charm_antidown_BR;
    static double const W_plus_to_charm_antistrange_BR;
    static double const W_plus_to_charmless_X_BR;
    static double const W_plus_to_up_antidown_BR;
    static double const W_plus_to_up_antistrange_BR;

    // CKM matrix elements, required for dividing up the hadronic decays:
    static double const V_ud;
    static double const V_us;
    static double const V_ud_squared;
    static double const V_us_squared;
    static double const one_over_V_ud_squared_plus_V_us_squared;
    static double const V_ud_squared_fraction;
    static double const V_us_squared_fraction;
    static double const V_cd;
    static double const V_cs;
    static double const V_cd_squared;
    static double const V_cs_squared;
    static double const one_over_V_cd_squared_plus_V_cs_squared;
    static double const V_cd_squared_fraction;
    static double const V_cs_squared_fraction;

    /* currently tau leptons are treated as stable in EW_scale_spectrum, but if
     * they are to be implemented, these values should be used.
     * all these values were taken from the PDG on 2009-11-10.
     */
    static double const tau_lepton_to_neutrinos_electron_BR;
    static double const tau_lepton_to_neutrinos_muon_BR;
    /* this code assumes that the rest of the decay width of the tau lepton is
     * divided between down + antiup and strange + antiup in the ratio
     * ( |V_ud|^2 ) to ( |V_us|^2 ) as in the case of the decays of the W^+:
     */
    static double const tau_lepton_to_neutrino_hadron_BR;
    static double const tau_lepton_to_neutrino_down_antiup_BR;
    static double const tau_lepton_to_neutrino_strange_antiup_BR;

  };

}  // end of CppSLHA namespace.

#endif /* CPPSLHA_PDG_CODES_AND_DATA_HPP_ */
