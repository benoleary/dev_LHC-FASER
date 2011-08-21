/*
 * LHC-FASER_derived_lepton_distributions.hpp
 *
 *  Created on: 05 Apr 5 2011
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
 *      LHC-FASER_base_lepton_distribution_stuff.hpp
 *      LHC-FASER_base_lepton_distribution_stuff.cpp
 *      LHC-FASER_derived_lepton_distributions.hpp
 *      LHC-FASER_derived_lepton_distributions.cpp
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

#ifndef LHC_FASER_DERIVED_LEPTON_DISTRIBUTIONS_HPP_
#define LHC_FASER_DERIVED_LEPTON_DISTRIBUTIONS_HPP_


#include "CppSLHA/CppSLHA.hpp"
#include "LHC-FASER_global_stuff.hpp"
#include "LHC-FASER_base_lepton_distribution_stuff.hpp"

namespace LHC_FASER
{

  /* this derived class sets up the energy distribution for a light lepton
   * which is the "near" lepton in a cascade decay, summing over charges, which
   * averages over the chirality effects & thus produces a flat distribution.
   */
  class flat_near_muon_plus_antimuon : public leptonEnergyDistribution
  {

  public:

    flat_near_muon_plus_antimuon( readier_for_new_point* const given_readier,
                                 CppSLHA::CppSLHA0 const* const given_spectrum,
              CppSLHA::particle_property_set const* const given_first_particle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const given_second_particle,
             CppSLHA::particle_property_set const* const given_third_particle )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~flat_near_muon_plus_antimuon()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    segment_term_set MINtoMAX_segment;
    // the terms between minimum_energy & maximum_energy.

    lepton_distribution_expansion_term MINtoMAX_const;
    // the term constant with respect to input_energy in the above segment.


    void
    calculateCoefficients()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };

  /* this derived class sets up the energy distribution for a light lepton
   * which is the "near" lepton in a cascade decay from a squark to a slepton,
   * where the near lepton or antilepton has the same helicity to the jet (thus
   * the distibution peaks at higher energies, hence "hard" ).  Actually, in
   * practice, this is from a cascade decay of a squark to a chargino & then
   * the chargino to a light lepton plus sneutrino.
   */
  class same_chirality_near_muon : public leptonEnergyDistribution
  {

  public:

    same_chirality_near_muon( readier_for_new_point* const given_readier,
                              CppSLHA::CppSLHA0 const* const given_spectrum,
              CppSLHA::particle_property_set const* const given_first_particle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const given_second_particle,
             CppSLHA::particle_property_set const* const given_third_particle )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~same_chirality_near_muon()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    segment_term_set MINtoMAX_segment;
    // the terms between minimum_energy & maximum_energy.

    lepton_distribution_expansion_term MINtoMAX_const;
    // the term constant with respect to input_energy in the above segment.
    lepton_distribution_expansion_term MINtoMAX_lin;
    // the term linear in input_energy in the above segment.


    void
    calculateCoefficients()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };

  /* this derived class sets up the energy distribution for a light lepton
   * which is the "near" lepton in a cascade decay from a squark to a slepton,
   * where the near lepton or antilepton has the opposite helicity to the jet
   * (thus the distibution peaks at lower energies, hence "soft" ).  Actually,
   * in practice, this is from a cascade decay of a squark to a chargino & then
   * the chargino to a light lepton plus sneutrino.
   */
  class opposite_chirality_near_muon : public leptonEnergyDistribution
  {


  public:

    opposite_chirality_near_muon( readier_for_new_point* const given_readier,
                                 CppSLHA::CppSLHA0 const* const given_spectrum,
              CppSLHA::particle_property_set const* const given_first_particle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const given_second_particle,
             CppSLHA::particle_property_set const* const given_third_particle )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~opposite_chirality_near_muon()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    segment_term_set MINtoMAX_segment;
    // the terms between minimum_energy & maximum_energy.

    lepton_distribution_expansion_term MINtoMAX_const;
    // the term constant with respect to input_energy in the above segment.
    lepton_distribution_expansion_term MINtoMAX_lin;
    // the term linear in input_energy in the above segment.


    void
    calculateCoefficients()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };

  /* this derived class sets up the energy distribution for a light lepton
   * which is the "far" lepton in a cascade decay, summing over charges, which
   * averages over the chirality effects & thus produces a relatively simple
   * distribution.
   */
  class flat_far_muon_plus_antimuon : public leptonEnergyDistribution
  {

  public:

    flat_far_muon_plus_antimuon( readier_for_new_point* const given_readier,
                                 CppSLHA::CppSLHA0 const* const given_spectrum,
              CppSLHA::particle_property_set const* const given_first_particle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const given_second_particle,
              CppSLHA::particle_property_set const* const given_third_particle,
            CppSLHA::particle_property_set const* const given_fourth_particle )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~flat_far_muon_plus_antimuon()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    double Elk;
    double Ehk;


    segment_term_set MINtoLK_segment;
    // the terms between minimum_energy & Elk.

    lepton_distribution_expansion_term MINtoLK_const;
    // the term constant with respect to input_energy in the above segment.
    lepton_distribution_expansion_term MINtoLK_log;
    // the term linear in the logarithm of input_energy in the above segment.


    segment_term_set LKtoHK_segment;
    // the terms between Elk & Ehk.

    lepton_distribution_expansion_term LKtoHK_const;
    // the term constant with respect to input_energy in the above segment.


    segment_term_set HKtoMAX_segment;
    // the terms between Ehk & maximum_energy.

    lepton_distribution_expansion_term HKtoMAX_const;
    // the term constant with respect to input_energy in the above segment.
    lepton_distribution_expansion_term HKtoMAX_log;
    // the term linear in the logarithm of input_energy in the above segment.


    void
    calculateCoefficients()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };

  /* this derived class sets up the energy distribution for a light lepton
   * which is the "far" lepton in a cascade decay from a squark to a slepton,
   * where the near lepton or antilepton has the same helicity to the jet (thus
   * the distibution peaks at lower energies, hence "soft" ).  Actually, in
   * practice, this is from a cascade decay of a squark to a chargino & then
   * the chargino to a neutrino plus slepton which then decays to the light
   * lepton.
   */
  class same_chirality_far_muon : public leptonEnergyDistribution
  {

  public:

    same_chirality_far_muon( readier_for_new_point* const given_readier,
                             CppSLHA::CppSLHA0 const* const given_spectrum,
              CppSLHA::particle_property_set const* const given_first_particle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const given_second_particle,
              CppSLHA::particle_property_set const* const given_third_particle,
            CppSLHA::particle_property_set const* const given_fourth_particle )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~same_chirality_far_muon()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    double Elk;
    double Ehk;


    segment_term_set MINtoLK_segment;
    // the terms between minimum_energy & Elk.

    lepton_distribution_expansion_term MINtoLK_inv;
    // the term linear in the inverse power of input_energy in the above
    // segment.
    lepton_distribution_expansion_term MINtoLK_const;
    // the term constant with respect to input_energy in the above segment.
    lepton_distribution_expansion_term MINtoLK_log;
    // the term linear in the logarithm of input_energy in the above segment.
    lepton_distribution_expansion_term MINtoLK_lin;
    // the term linear in input_energy in the above segment.


    segment_term_set LKtoHK_segment;
    // the terms between Elk & Ehk.

    lepton_distribution_expansion_term LKtoHK_inv;
    // the term linear in the inverse power of input_energy in the above
    // segment.
    lepton_distribution_expansion_term LKtoHK_const;
    // the term constant with respect to input_energy in the above segment.


    segment_term_set HKtoMAX_segment;
    // the terms between Ehk & maximum_energy.

    lepton_distribution_expansion_term HKtoMAX_inv;
    // the term linear in the inverse power of input_energy in the above
    // segment.
    lepton_distribution_expansion_term HKtoMAX_const;
    // the term constant with respect to input_energy in the above segment.
    lepton_distribution_expansion_term HKtoMAX_log;
    // the term linear in the logarithm of input_energy in the above segment.
    lepton_distribution_expansion_term HKtoMAX_lin;
    // the term linear in input_energy in the above segment.


    void
    calculateCoefficients()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };

  /* this derived class sets up the energy distribution for a light lepton
   * which is the "far" lepton in a cascade decay from a squark to a slepton,
   * where the near lepton or antilepton has the opposite helicity to the jet
   * (thus the distibution peaks at higher energies, hence "hard" ).  Actually,
   * in practice, this is from a cascade decay of a squark to a chargino & then
   * the chargino to a neutrino plus slepton which then decays to the light
   * lepton.
   */
  class opposite_chirality_far_muon : public leptonEnergyDistribution
  {

  public:

    opposite_chirality_far_muon( readier_for_new_point* const given_readier,
                                 CppSLHA::CppSLHA0 const* const given_spectrum,
              CppSLHA::particle_property_set const* const given_first_particle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const given_second_particle,
              CppSLHA::particle_property_set const* const given_third_particle,
            CppSLHA::particle_property_set const* const given_fourth_particle )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~opposite_chirality_far_muon()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    double Elk;
    double Ehk;


    segment_term_set MINtoLK_segment;
    // the terms between minimum_energy & Elk.

    lepton_distribution_expansion_term MINtoLK_inv;
    // the term linear in the inverse power of input_energy in the above
    // segment.
    lepton_distribution_expansion_term MINtoLK_const;
    // the term constant with respect to input_energy in the above segment.
    lepton_distribution_expansion_term MINtoLK_log;
    // the term linear in the logarithm of input_energy in the above segment.
    lepton_distribution_expansion_term MINtoLK_lin;
    // the term linear in input_energy in the above segment.


    segment_term_set LKtoHK_segment;
    // the terms between Elk & Ehk.

    lepton_distribution_expansion_term LKtoHK_inv;
    // the term linear in the inverse power of input_energy in the above
    // segment.
    lepton_distribution_expansion_term LKtoHK_const;
    // the term constant with respect to input_energy in the above segment.


    segment_term_set HKtoMAX_segment;
    // the terms between Ehk & maximum_energy.

    lepton_distribution_expansion_term HKtoMAX_inv;
    // the term linear in the inverse power of input_energy in the above
    // segment.
    lepton_distribution_expansion_term HKtoMAX_const;
    // the term constant with respect to input_energy in the above segment.
    lepton_distribution_expansion_term HKtoMAX_log;
    // the term linear in the logarithm of input_energy in the above segment.
    lepton_distribution_expansion_term HKtoMAX_lin;
    // the term linear in input_energy in the above segment.


    void
    calculateCoefficients()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };

  /* the negatively- & positively-charged leptons from neutral EWSB scalar
   * decays have the same distribution thanks to the spin-0 nature of the
   * bosons.
   */
  /* this derived class sets up the energy distribution for a light lepton
   * which is from the decay of an on-shell Higgs boson from the cascade decay
   * of a squark to a jet of either chirality plus a light lepton-antilepton
   * pair plus a lightest neutralino, summing over charges, which averages over
   * the chirality effects & thus produces a relatively simple distribution.
   */
  class Higgs_muon_plus_antimuon : public leptonEnergyDistribution
  {

  public:

      Higgs_muon_plus_antimuon( readier_for_new_point* const given_readier,
                                CppSLHA::CppSLHA0 const* const given_spectrum,
              CppSLHA::particle_property_set const* const given_first_particle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const given_second_particle,
              CppSLHA::particle_property_set const* const given_third_particle,
            CppSLHA::particle_property_set const* const given_fourth_particle )
      /* code after the classes in this .hpp file, or in the .cpp file. */;

      virtual
      ~Higgs_muon_plus_antimuon()
      /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

      double Elk;
      double Ehk;


      segment_term_set MINtoLK_segment;
      // the terms between minimum_energy & Elk.

      lepton_distribution_expansion_term MINtoLK_const;
      // the term constant with respect to input_energy in the above segment.
      lepton_distribution_expansion_term MINtoLK_log;
      // the term linear in the logarithm of input_energy in the above segment.


      segment_term_set LKtoHK_segment;
      // the terms between Elk & Ehk.

      lepton_distribution_expansion_term LKtoHK_const;
      // the term constant with respect to input_energy in the above segment.


      segment_term_set HKtoMAX_segment;
      // the terms between Ehk & maximum_energy.

      lepton_distribution_expansion_term HKtoMAX_const;
      // the term constant with respect to input_energy in the above segment.
      lepton_distribution_expansion_term HKtoMAX_log;
      // the term linear in the logarithm of input_energy in the above segment.


      void
      calculateCoefficients()
      /* code after the classes in this .hpp file, or in the .cpp file. */;

  };

  /* specific lepton or antilepton versions of the Z-channel distribution could
   * be done, but at the moment there is no need for them.
   */
  /* this derived class sets up the energy distribution for a light lepton
   * which is from the decay of an on-shell Z boson from the cascade decay of a
   * squark to a jet of either chirality plus a light lepton-antilepton pair
   * plus a lightest neutralino, keeping chirality effects.
   */
  class Z_handed_muon : public leptonEnergyDistribution
  {

  public:

    Z_handed_muon( readier_for_new_point* const given_readier,
                   CppSLHA::CppSLHA0 const* const given_spectrum,
              CppSLHA::particle_property_set const* const given_first_particle,
                   effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const given_second_particle,
              CppSLHA::particle_property_set const* const given_third_particle,
            CppSLHA::particle_property_set const* const given_fourth_particle,
                   bool const negative_muon_is_same_handedness,
                   bool const should_sum_over_handedness )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~Z_handed_muon()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    double Elk;
    double Ehk;
    double same_handedness_factor;
    double axial_coupling_factor;
    bool energy_limited_cosines;

    /* these are for ease of calculating the coefficients (referring to the
     * squark as Q, the neutralino_2 as C, the neutralino_1 as X, the jet as j,
     * the Z boson as V, the negative muon as l, & the positive muon as v; it
     * makes more sense to know that the notation came from a calculation with
     * a chargino decaying to a W boson...):
     */
    double ElMin;
    // the minimum muon energy in the neutralino_2 rest frame.
    double lnElMin;
    // ln( ElMin ).
    double ElMinsq;
    // ElMin^2.
    double mQ;
    // the mass of the squark.
    double mQsq;
    // mQ^2.
    double mQcu;
    // mQ^3.
    double mC;
    // the mass of the neutralino_2.
    double mCsq;
    // mC^2.
    double mQCsqSum;
    // mQ^2 + mC^2.
    double mQCsqDiff;
    // mQ^2 - mC^2.
    double lnmQC;
    // ln( ( mQ / mC ) ).
    double mCmX;
    // mC * the mass of the neutralino_1.
    double mXsq;
    // mX^2.
    double mVsq;
    // the square of the mass of the Z boson.
    double EV;
    // the energy of the Z boson in the rest frame of neutralino_2.
    double EightmCEVElMin;
    // 8 * mC * EV * ElMin.
    double lnVsqOverFourElMinsq;
    // ln( ( mVsq / ( 4 * ElMin^2 ) ).
    double EightmCXElMinmQsq;
    // 8 * mC * mX * ElMin * mQ^2.
    double lnEmin;
    // ln( minimum_energy ).
    double lnElk;
    // ln( Elk ).
    double lnEhk;
    // ln( Ehk ).
    double lnEmax;
    // ln( maximum_energy ).


    double current_term_coefficient;
    // this is for holding the coefficient being calculated.

    segment_term_set MINtoLK_segment;
    // the terms between minimum_energy & Elk.

    lepton_distribution_expansion_term MINtoLK_inv;
    // the term linear in the inverse power of input_energy in the above
    // segment.
    lepton_distribution_expansion_term MINtoLK_const;
    // the term constant with respect to input_energy in the above segment.
    lepton_distribution_expansion_term MINtoLK_log;
    // the term linear in the logarithm of input_energy in the above segment.
    lepton_distribution_expansion_term MINtoLK_lin;
    // the term linear in input_energy in the above segment.
    lepton_distribution_expansion_term MINtoLK_linlog;
    // the term linear in input_energy & in the logarithm of input_energy in
    // the above segment.
    lepton_distribution_expansion_term MINtoLK_sq;
    // the term quadratic in input_energy in the above segment.


    segment_term_set LKtoHK_segment;
    // the terms between Elk & Ehk.

    lepton_distribution_expansion_term LKtoHK_const;
    // the term constant with respect to input_energy in the above segment.
    lepton_distribution_expansion_term LKtoHK_lin;
    // the term linear in input_energy in the above
    // segment.
    lepton_distribution_expansion_term LKtoHK_sq;
    // the term quadratic in input_energy in the above segment.


    segment_term_set HKtoMAX_segment;
    // the terms between Ehk & maximum_energy.

    lepton_distribution_expansion_term HKtoMAX_inv;
    // the term linear in the inverse power of input_energy in the above
    // segment.
    lepton_distribution_expansion_term HKtoMAX_const;
    // the term constant with respect to input_energy in the above segment.
    lepton_distribution_expansion_term HKtoMAX_log;
    // the term linear in the logarithm of input_energy in the above segment.
    lepton_distribution_expansion_term HKtoMAX_lin;
    // the term linear in input_energy in the above segment.
    lepton_distribution_expansion_term HKtoMAX_linlog;
    // the term linear in input_energy & in the logarithm of input_energy in
    // the above segment.
    lepton_distribution_expansion_term HKtoMAX_sq;
    // the term quadratic in input_energy in the above segment.


    void
    calculateCoefficients()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };


  /* this derived class sets up the energy distribution for a light lepton
   * which is from the decay of an on-shell negatively-charged W boson from the
   * cascade decay of a down-type squark or an up-type antisquark.
   * unfortunately, the couplings are very model-specific, so this assumes just
   * the MSSM.
   */
  class W_minus_handed_muon : public leptonEnergyDistribution
  {

  public:

    W_minus_handed_muon( readier_for_new_point* const given_readier,
                         CppSLHA::CppSLHA0 const* const given_spectrum,
              CppSLHA::particle_property_set const* const given_first_particle,
                         effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const given_second_particle,
              CppSLHA::particle_property_set const* const given_third_particle,
            CppSLHA::particle_property_set const* const given_fourth_particle,
                         bool const given_jet_is_left_handed )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~W_minus_handed_muon()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    double Elk;
    double Ehk;
    bool const jet_is_left_handed;
    bool energy_limited_cosines;

    int which_chargino;
    int which_neutralino;
    CppSLHA::SLHA_BLOCK const* NMIX;
    CppSLHA::SLHA_BLOCK const* UMIX;
    CppSLHA::SLHA_BLOCK const* VMIX;

    /* these are for ease of calculating the coefficients (referring to the
     * squark as Q, the chargino as C, the neutralino as X, the jet as j,
     * the W boson as V, the negative muon as l, & the neutrino as v):
     */
    double Vv;
    // the vector coupling of the W to the electroweakinos.
    double Vvsq;
    // Vv^sq.
    double Aa;
    // the axial coupling of the W to the electroweakinos.
    double Aasq;
    // Aa^2.
    double TwiceVvAa;
    // 2 * Vv * Aa.
    // Jj is +1 if the jet is left-handed, -1 if right-handed:
    double VvPlusJjAaAllSq;
    // ( Vv + Jj * Aa )^2.
    double VvSqMinusAaSq;
    // ( Vv^2 - Aa^2 ).
    double MinusTwiceVvJjAa;
    // -2.0 * Vv * Jj * Aa.
    double JjVvPlusJjAaAllSq;
    // Jj * ( Vv + Jj * Aa )^2.
    double JjVvSqMinusAaSq;
    // Jj * ( Vv^2 - Aa^2 ).
    double JjVvSqPlusAaSq;
    // Jj * ( Vv^2 + Aa^2 ).
    double ElMin;
    // the minimum muon energy in the chargino rest frame.
    double lnElMin;
    // ln( ElMin ).
    double ElMinsq;
    // ElMin^2.
    double mQ;
    // the mass of the squark.
    double mQsq;
    // mQ^2.
    double mQcu;
    // mQ^3.
    double mC;
    // the mass of the chargino.
    double mCsq;
    // mC^2.
    double mQCsqSum;
    // mQ^2 + mC^2.
    double mQCsqDiff;
    // mQ^2 - mC^2.
    double mQCsqDiffsq;
    // mQCsqDiff^2.
    double lnmQC;
    // ln( ( mQ / mC ) ).
    double mCmX;
    // mC * the mass of the neutralino.
    double mXsq;
    // mX^2.
    double mVsq;
    // the square of the mass of the W boson.
    double mVqu;
    // mVsq^2.
    double EV;
    // the energy of the W boson in the rest frame of the chargino.
    double EightmCEVElMin;
    // 8 * mC * EV * ElMin.
    double lnVsqOverFourElMinsq;
    // ln( ( mVsq / ( 4 * ElMin^2 ) ).
    double EightmCXElMinmQsq;
    // 8 * mC * mX * ElMin * mQ^2.
    double lnEmin;
    // ln( minimum_energy ).
    double lnElk;
    // ln( Elk ).
    double lnEhk;
    // ln( Ehk ).
    double lnEmax;
    // ln( maximum_energy ).


    double current_term_coefficient;
    // this is for holding the coefficient being calculated.

    segment_term_set MINtoLK_segment;
    // the terms between minimum_energy & Elk.

    lepton_distribution_expansion_term MINtoLK_inv;
    // the term linear in the inverse power of input_energy in the above
    // segment.
    lepton_distribution_expansion_term MINtoLK_const;
    // the term constant with respect to input_energy in the above segment.
    lepton_distribution_expansion_term MINtoLK_log;
    // the term linear in the logarithm of input_energy in the above segment.
    lepton_distribution_expansion_term MINtoLK_lin;
    // the term linear in input_energy in the above segment.
    lepton_distribution_expansion_term MINtoLK_linlog;
    // the term linear in input_energy & in the logarithm of input_energy in
    // the above segment.
    lepton_distribution_expansion_term MINtoLK_sq;
    // the term quadratic in input_energy in the above segment.


    segment_term_set LKtoHK_segment;
    // the terms between Elk & Ehk.

    lepton_distribution_expansion_term LKtoHK_const;
    // the term constant with respect to input_energy in the above segment.
    lepton_distribution_expansion_term LKtoHK_lin;
    // the term linear in input_energy in the above
    // segment.
    lepton_distribution_expansion_term LKtoHK_sq;
    // the term quadratic in input_energy in the above segment.


    segment_term_set HKtoMAX_segment;
    // the terms between Ehk & maximum_energy.

    lepton_distribution_expansion_term HKtoMAX_inv;
    // the term linear in the inverse power of input_energy in the above
    // segment.
    lepton_distribution_expansion_term HKtoMAX_const;
    // the term constant with respect to input_energy in the above segment.
    lepton_distribution_expansion_term HKtoMAX_log;
    // the term linear in the logarithm of input_energy in the above segment.
    lepton_distribution_expansion_term HKtoMAX_lin;
    // the term linear in input_energy in the above segment.
    lepton_distribution_expansion_term HKtoMAX_linlog;
    // the term linear in input_energy & in the logarithm of input_energy in
    // the above segment.
    lepton_distribution_expansion_term HKtoMAX_sq;
    // the term quadratic in input_energy in the above segment.


    void
    calculateCoefficients()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };

  /* this derived class sets up the energy distribution for a light lepton
   * which is from the decay of an on-shell gauge vector boson from the
   * cascade decay of a squark to a lighter squark.
   */
  class vector_from_squark_to_muon : public leptonEnergyDistribution
  {

  public:

    vector_from_squark_to_muon( readier_for_new_point* const given_readier,
                                CppSLHA::CppSLHA0 const* const given_spectrum,
              CppSLHA::particle_property_set const* const given_first_particle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const given_second_particle,
             CppSLHA::particle_property_set const* const given_third_particle )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~vector_from_squark_to_muon()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    // these are for ease of calculating the coefficients (referring to the
    // decaying squark as Qh & the product squark as Ql):
    double mQhSq;
    // the square of the mass of the heavier squark.
    double mQlSq;
    // the square of the mass of the lighter squark.
    double mVSq;
    // the square of the mass of the vector boson.
    double gammaV;
    // the gamma factor of the boost from the Qh rest frame to the vector boson
    // rest frame.
    double betaV;
    // the beta factor of the boost from the Qh rest frame to the vector boson
    // rest frame.
    double energyDifference;
    // maximum_energy - minimum_energy

    segment_term_set MINtoMAX_segment;
    // the terms between minimum_energy & maximum_energy.

    lepton_distribution_expansion_term MINtoMAX_const;
    // the term constant with respect to input_energy in the above segment.
    lepton_distribution_expansion_term MINtoMAX_lin;
    // the term linear in input_energy in the above segment.
    lepton_distribution_expansion_term MINtoMAX_sq;
    // the term quadratic in input_energy in the above segment.


    void
    calculateCoefficients()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };

  // averaging over both jets from the Z boson produces the same distribution
  // as directly-produced light leptons.
  class Z_direct_jet : public Z_handed_muon
  {

  public:

    Z_direct_jet( readier_for_new_point* const given_readier,
                  CppSLHA::CppSLHA0 const* const given_spectrum,
            CppSLHA::particle_property_set const* const given_first_particle,
                  effectiveSquarkMassHolder* const effectiveSquarkMass,
            CppSLHA::particle_property_set const* const given_second_particle,
            CppSLHA::particle_property_set const* const given_third_particle,
            CppSLHA::particle_property_set const* const given_fourth_particle )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~Z_direct_jet()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  //protected:

    // nothing.

  };

  // averaging over both jets from the W boson produces the same distribution
  // which Z bosons produce.
  class W_minus_direct_jet : public Z_handed_muon

  {

  public:

    W_minus_direct_jet( readier_for_new_point* const given_readier,
                        CppSLHA::CppSLHA0 const* const given_spectrum,
              CppSLHA::particle_property_set const* const given_first_particle,
                        effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const given_second_particle,
              CppSLHA::particle_property_set const* const given_third_particle,
            CppSLHA::particle_property_set const* const given_fourth_particle )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~W_minus_direct_jet()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    //protected:

      // nothing.

  };


  /* FROM THIS POINT ON ARE DERIVED CLASSES WHICH STILL NEED WORK!
   * ALL ARE CURRENTLY ROUGH GUESSES!
   * direct light leptons: charged EWSB scalar decay
   *   - equations in progress...
   * direct light leptons: 3-body neutralino decay
   *   - equations in progress...
   * direct light leptons: 3-body chargino decay
   *   - equations in progress...
   * direct jets: charged EWSB scalar decay
   *   - should be a quick modification of the light lepton case.
   * direct jets: 3-body neutralino decay
   *   - should be a quick modification of the light lepton case.
   * direct jets: 3-body chargino decay
   *   - should be a quick modification of the light lepton case.
   */

  /* placeholder of neutral Higgs distribution. it's different, because the
   * chiral couplings of the charged Higgs with the chargino & neutralino make
   * the neutralino prefer to be aligned or anti-aligned (depending on the
   * couplings) with the jet, & thus the charged Higgs is preferentially
   * boosted more or less compared to the neutral Higgs case, & thus the
   * muon distribution is different.
   */
  class negatively_charged_Higgs_muon : public Higgs_muon_plus_antimuon
  {

  public:

    negatively_charged_Higgs_muon( readier_for_new_point* const given_readier,
                                 CppSLHA::CppSLHA0 const* const given_spectrum,
              CppSLHA::particle_property_set const* const given_first_particle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const given_second_particle,
              CppSLHA::particle_property_set const* const given_third_particle,
          CppSLHA::particle_property_set const* const given_fourth_particle )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~negatively_charged_Higgs_muon()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    //protected:

      // nothing.

  };

  /* placeholder of neutral Higgs distribution. it's different, because the
   * chiral couplings of the charged Higgs with the chargino & neutralino make
   * the neutralino prefer to be aligned or anti-aligned (depending on the
   * couplings) with the jet, & thus the charged Higgs is preferentially
   * boosted more or less compared to the neutral Higgs case, & thus the
   * jet distribution is different.
   */
  class negatively_charged_Higgs_jet : public Higgs_muon_plus_antimuon
  {

  public:

    negatively_charged_Higgs_jet( readier_for_new_point* const given_readier,
                                 CppSLHA::CppSLHA0 const* const given_spectrum,
              CppSLHA::particle_property_set const* const given_first_particle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const given_second_particle,
              CppSLHA::particle_property_set const* const given_third_particle,
          CppSLHA::particle_property_set const* const given_fourth_particle )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~negatively_charged_Higgs_jet()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    //protected:

      // nothing.

  };



  /* this needs to know the 2 sfermions which appear off-shell in the process.
   * they are given as "left" & "right", but the code checks for 3rd-generation
   * sfermions & treats them accordingly. it combines the contribution from the
   * off-shell sfermions with the contributions from off-shell Z & spin-0
   * bosons, though without doing interference between the channels properly.
   *
   * currently this is a placeholder guess of a flat distribution from 0.0 GeV
   * to the neutralino mass difference.
   */
  class neutralino_three_body_decay : public leptonEnergyDistribution
  {

  public:

    neutralino_three_body_decay( readier_for_new_point* const given_readier,
                                 CppSLHA::CppSLHA0 const* const given_spectrum,
              CppSLHA::particle_property_set const* const given_first_particle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const given_second_particle,
              CppSLHA::particle_property_set const* const given_third_particle,
               CppSLHA::particle_property_set const* const given_left_sfermion,
             CppSLHA::particle_property_set const* const given_right_sfermion )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~neutralino_three_body_decay()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    CppSLHA::particle_property_set const* const left_sfermion;
    CppSLHA::particle_property_set const* const right_sfermion;

    segment_term_set MINtoMAX_segment;
    // the terms between minimum_energy & maximum_energy.

    lepton_distribution_expansion_term MINtoMAX_const;
    // the term constant with respect to input_energy in the above segment.


    void
    calculateCoefficients()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };

  /* this needs to know the 2 sfermions which appear off-shell in the process.
   * they are given as "left" & "right", but the code checks for 3rd-generation
   * sfermions & treats them accordingly. it combines the contribution from the
   * off-shell sfermions with the contributions from off-shell W & spin-0
   * bosons, though without doing interference between the channels properly.
   *
   * currently this is a placeholder guess of a flat distribution from 0.0 GeV
   * to the mass difference between the chargino & the neutralino.
   */
  class chargino_three_body_decay : public leptonEnergyDistribution
  {

  public:

    chargino_three_body_decay( readier_for_new_point* const given_readier,
                               CppSLHA::CppSLHA0 const* const given_spectrum,
              CppSLHA::particle_property_set const* const given_first_particle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const given_second_particle,
              CppSLHA::particle_property_set const* const given_third_particle,
             CppSLHA::particle_property_set const* const leftUpIsospinSfermion,
            CppSLHA::particle_property_set const* const rightUpIsospinSfermion,
           CppSLHA::particle_property_set const* const leftDownIsospinSfermion,
         CppSLHA::particle_property_set const* const rightDownIsospinSfermion )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~chargino_three_body_decay()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    CppSLHA::particle_property_set const* const leftUpIsospinSfermion;
    CppSLHA::particle_property_set const* const rightUpIsospinSfermion;
    CppSLHA::particle_property_set const* const leftDownIsospinSfermion;
    CppSLHA::particle_property_set const* const rightDownIsospinSfermion;

    segment_term_set MINtoMAX_segment;
    // the terms between minimum_energy & maximum_energy.

    lepton_distribution_expansion_term MINtoMAX_const;
    // the term constant with respect to input_energy in the above segment.


    void
    calculateCoefficients()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };





  // inline functions:

  // none.

}  // end of LHC_FASER namespace.

#endif /* LHC_FASER_DERIVED_LEPTON_DISTRIBUTIONS_HPP_ */
