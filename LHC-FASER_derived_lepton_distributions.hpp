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

#include "LHC-FASER_base_lepton_distribution_stuff.hpp"

namespace LHC_FASER
{

  /* this derived class sets up the energy distribution for a light lepton
   * which is the "near" lepton in a cascade decay, summing over charges, which
   * averages over the chirality effects & thus produces a flat distribution.
   */
  class flatNearMuonPlusAntimuon : public leptonEnergyDistribution
  {
  public:
    flatNearMuonPlusAntimuon( readierForNewPoint* const readierPointer,
                              CppSLHA::CppSLHA0 const* const spectrumData,
                     CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                    CppSLHA::particle_property_set const* const thirdParticle )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    virtual
    ~flatNearMuonPlusAntimuon()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    segmentTermSet minToMaxSegment;
    // the terms between minimumEnergy & maximumEnergy.
    leptonDistributionExpansionTerm* const minToMaxConst;
    // the term constant with respect to inputEnergy in the above segment.

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
  class sameChiralityNearMuon : public leptonEnergyDistribution
  {
  public:
    sameChiralityNearMuon( readierForNewPoint* const readierPointer,
                           CppSLHA::CppSLHA0 const* const spectrumData,
                     CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                    CppSLHA::particle_property_set const* const thirdParticle )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    virtual
    ~sameChiralityNearMuon()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    segmentTermSet minToMaxSegment;
    // the terms between minimumEnergy & maximumEnergy.
    leptonDistributionExpansionTerm* const minToMaxConst;
    // the term constant with respect to inputEnergy in the above segment.
    leptonDistributionExpansionTerm* const minToMaxLin;
    // the term linear in inputEnergy in the above segment.

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
  class oppositeChiralityNearMuon : public leptonEnergyDistribution
  {
  public:
    oppositeChiralityNearMuon( readierForNewPoint* const readierPointer,
                               CppSLHA::CppSLHA0 const* const spectrumData,
                     CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                    CppSLHA::particle_property_set const* const thirdParticle )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    virtual
    ~oppositeChiralityNearMuon()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    segmentTermSet minToMaxSegment;
    // the terms between minimumEnergy & maximumEnergy.
    leptonDistributionExpansionTerm* const minToMaxConst;
    // the term constant with respect to inputEnergy in the above segment.
    leptonDistributionExpansionTerm* const minToMaxLin;
    // the term linear in inputEnergy in the above segment.

    void
    calculateCoefficients()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };

  /* this derived class sets up the energy distribution for a light lepton
   * which is the "far" lepton in a cascade decay, summing over charges, which
   * averages over the chirality effects & thus produces a relatively simple
   * distribution.
   */
  class flatFarMuonPlusAntimuon : public leptonEnergyDistribution
  {
  public:
    flatFarMuonPlusAntimuon( readierForNewPoint* const readierPointer,
                             CppSLHA::CppSLHA0 const* const spectrumData,
              CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                     CppSLHA::particle_property_set const* const thirdParticle,
                   CppSLHA::particle_property_set const* const fourthParticle )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    virtual
    ~flatFarMuonPlusAntimuon()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    double Elk;
    double Ehk;

    segmentTermSet minToLkSegment;
    // the terms between minimumEnergy & Elk.
    leptonDistributionExpansionTerm* const minToLkConst;
    // the term constant with respect to inputEnergy in the above segment.
    leptonDistributionExpansionTerm* const minToLkLog;
    // the term linear in the logarithm of inputEnergy in the above segment.

    segmentTermSet lkToHkSegment;
    // the terms between Elk & Ehk.
    leptonDistributionExpansionTerm* const lkToHkConst;
    // the term constant with respect to inputEnergy in the above segment.

    segmentTermSet hkToMaxSegment;
    // the terms between Ehk & maximumEnergy.
    leptonDistributionExpansionTerm* const hkToMaxConst;
    // the term constant with respect to inputEnergy in the above segment.
    leptonDistributionExpansionTerm* const hkToMaxLog;
    // the term linear in the logarithm of inputEnergy in the above segment.

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
  class sameChiralityFarMuon : public leptonEnergyDistribution
  {
  public:
    sameChiralityFarMuon( readierForNewPoint* const readierPointer,
                          CppSLHA::CppSLHA0 const* const spectrumData,
              CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                     CppSLHA::particle_property_set const* const thirdParticle,
                   CppSLHA::particle_property_set const* const fourthParticle )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    virtual
    ~sameChiralityFarMuon()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    double Elk;
    double Ehk;

    segmentTermSet minToLkSegment;
    // the terms between minimumEnergy & Elk.
    leptonDistributionExpansionTerm* const minToLkInv;
    // the term linear in the inverse power of inputEnergy in the above
    // segment.
    leptonDistributionExpansionTerm* const minToLkConst;
    // the term constant with respect to inputEnergy in the above segment.
    leptonDistributionExpansionTerm* const minToLkLog;
    // the term linear in the logarithm of inputEnergy in the above segment.
    leptonDistributionExpansionTerm* const minToLkLin;
    // the term linear in inputEnergy in the above segment.

    segmentTermSet lkToHkSegment;
    // the terms between Elk & Ehk.
    leptonDistributionExpansionTerm* const lkToHkInv;
    // the term linear in the inverse power of inputEnergy in the above
    // segment.
    leptonDistributionExpansionTerm* const lkToHkConst;
    // the term constant with respect to inputEnergy in the above segment.

    segmentTermSet hkToMaxSegment;
    // the terms between Ehk & maximumEnergy.
    leptonDistributionExpansionTerm* const hkToMaxInv;
    // the term linear in the inverse power of inputEnergy in the above
    // segment.
    leptonDistributionExpansionTerm* const hkToMaxConst;
    // the term constant with respect to inputEnergy in the above segment.
    leptonDistributionExpansionTerm* const hkToMaxLog;
    // the term linear in the logarithm of inputEnergy in the above segment.
    leptonDistributionExpansionTerm* const hkToMaxLin;
    // the term linear in inputEnergy in the above segment.

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
  class oppositeChiralityFarMuon : public leptonEnergyDistribution
  {
  public:
    oppositeChiralityFarMuon( readierForNewPoint* const readierPointer,
                              CppSLHA::CppSLHA0 const* const spectrumData,
              CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                     CppSLHA::particle_property_set const* const thirdParticle,
                   CppSLHA::particle_property_set const* const fourthParticle )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    virtual
    ~oppositeChiralityFarMuon()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:
    double Elk;
    double Ehk;

    segmentTermSet minToLkSegment;
    // the terms between minimumEnergy & Elk.
    leptonDistributionExpansionTerm* const minToLkInv;
    // the term linear in the inverse power of inputEnergy in the above
    // segment.
    leptonDistributionExpansionTerm* const minToLkConst;
    // the term constant with respect to inputEnergy in the above segment.
    leptonDistributionExpansionTerm* const minToLkLog;
    // the term linear in the logarithm of inputEnergy in the above segment.
    leptonDistributionExpansionTerm* const minToLkLin;
    // the term linear in inputEnergy in the above segment.

    segmentTermSet lkToHkSegment;
    // the terms between Elk & Ehk.
    leptonDistributionExpansionTerm* const lkToHkInv;
    // the term linear in the inverse power of inputEnergy in the above
    // segment.
    leptonDistributionExpansionTerm* const lkToHkConst;
    // the term constant with respect to inputEnergy in the above segment.

    segmentTermSet hkToMaxSegment;
    // the terms between Ehk & maximumEnergy.
    leptonDistributionExpansionTerm* const hkToMaxInv;
    // the term linear in the inverse power of inputEnergy in the above
    // segment.
    leptonDistributionExpansionTerm* const hkToMaxConst;
    // the term constant with respect to inputEnergy in the above segment.
    leptonDistributionExpansionTerm* const hkToMaxLog;
    // the term linear in the logarithm of inputEnergy in the above segment.
    leptonDistributionExpansionTerm* const hkToMaxLin;
    // the term linear in inputEnergy in the above segment.

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
  class HiggsMuonPlusAntimuon : public leptonEnergyDistribution
  {
  public:
      HiggsMuonPlusAntimuon( readierForNewPoint* const readierPointer,
                             CppSLHA::CppSLHA0 const* const spectrumData,
              CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                     CppSLHA::particle_property_set const* const thirdParticle,
                   CppSLHA::particle_property_set const* const fourthParticle )
      /* code after the classes in this .hpp file, or in the .cpp file. */;
      virtual
      ~HiggsMuonPlusAntimuon()
      /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
      double Elk;
      double Ehk;

      segmentTermSet minToLkSegment;
      // the terms between minimumEnergy & Elk.
      leptonDistributionExpansionTerm* const minToLkConst;
      // the term constant with respect to inputEnergy in the above segment.
      leptonDistributionExpansionTerm* const minToLkLog;
      // the term linear in the logarithm of inputEnergy in the above segment.

      segmentTermSet lkToHkSegment;
      // the terms between Elk & Ehk.
      leptonDistributionExpansionTerm* const lkToHkConst;
      // the term constant with respect to inputEnergy in the above segment.

      segmentTermSet hkToMaxSegment;
      // the terms between Ehk & maximumEnergy.
      leptonDistributionExpansionTerm* const hkToMaxConst;
      // the term constant with respect to inputEnergy in the above segment.
      leptonDistributionExpansionTerm* const hkToMaxLog;
      // the term linear in the logarithm of inputEnergy in the above segment.

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
  class zHandedMuon : public leptonEnergyDistribution
  {
  public:
    zHandedMuon( readierForNewPoint* const readierPointer,
                 CppSLHA::CppSLHA0 const* const spectrumData,
                 CppSLHA::particle_property_set const* const firstParticle,
                 effectiveSquarkMassHolder* const effectiveSquarkMass,
                 CppSLHA::particle_property_set const* const secondParticle,
                 CppSLHA::particle_property_set const* const thirdParticle,
                 CppSLHA::particle_property_set const* const fourthParticle,
                 bool const negativeMuonIsSameHandednessAsJet,
                 bool const shouldSumOverHandedness )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    virtual
    ~zHandedMuon()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    bool const negativeMuonIsSameHandednessAsJet;
    bool const shouldNotSumOverHandedness;
    double Elk;
    double Ehk;
    double sameHandednessFactor;
    double axialCouplingFactor;
    bool couplesAsVectorNotAxial;
    bool cosinesLimitedByEnergy;

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
    // ln( minimumEnergy ).
    double lnElk;
    // ln( Elk ).
    double lnEhk;
    // ln( Ehk ).
    double lnEmax;
    // ln( maximumEnergy ).

    // these are for holding the coefficients being calculated:
    double currentMinToLkInvCoefficient;
    double currentMinToLkConstCoefficient;
    double currentMinToLkLogCoefficient;
    double currentMinToLkLinCoefficient;
    double currentMinToLkLinLogCoefficient;
    double currentMinToLkSqCoefficient;
    double currentLkToHkConstCoefficient;
    double currentLkToHkLinCoefficient;
    double currentLkToHkSqCoefficient;
    double currentHkToMaxInvCoefficient;
    double currentHkToMaxConstCoefficient;
    double currentHkToMaxLogCoefficient;
    double currentHkToMaxLinCoefficient;
    double currentHkToMaxLinLogCoefficient;
    double currentHkToMaxSqCoefficient;

    segmentTermSet minToLkSegment;
    // the terms between minimumEnergy & Elk.
    leptonDistributionExpansionTerm* const minToLkInv;
    // the term linear in the inverse power of inputEnergy in the above
    // segment.
    leptonDistributionExpansionTerm* const minToLkConst;
    // the term constant with respect to inputEnergy in the above segment.
    leptonDistributionExpansionTerm* const minToLkLog;
    // the term linear in the logarithm of inputEnergy in the above segment.
    leptonDistributionExpansionTerm* const minToLkLin;
    // the term linear in inputEnergy in the above segment.
    leptonDistributionExpansionTerm* const minToLkLinlog;
    // the term linear in inputEnergy & in the logarithm of inputEnergy in
    // the above segment.
    leptonDistributionExpansionTerm* const minToLkSq;
    // the term quadratic in inputEnergy in the above segment.

    segmentTermSet lkToHkSegment;
    // the terms between Elk & Ehk.
    leptonDistributionExpansionTerm* const lkToHkConst;
    // the term constant with respect to inputEnergy in the above segment.
    leptonDistributionExpansionTerm* const lkToHkLin;
    // the term linear in inputEnergy in the above segment.
    leptonDistributionExpansionTerm* const lkToHkSq;
    // the term quadratic in inputEnergy in the above segment.

    segmentTermSet hkToMaxSegment;
    // the terms between Ehk & maximumEnergy.
    leptonDistributionExpansionTerm* const hkToMaxInv;
    // the term linear in the inverse power of inputEnergy in the above
    // segment.
    leptonDistributionExpansionTerm* const hkToMaxConst;
    // the term constant with respect to inputEnergy in the above segment.
    leptonDistributionExpansionTerm* const hkToMaxLog;
    // the term linear in the logarithm of inputEnergy in the above segment.
    leptonDistributionExpansionTerm* const hkToMaxLin;
    // the term linear in inputEnergy in the above segment.
    leptonDistributionExpansionTerm* const hkToMaxLinlog;
    // the term linear in inputEnergy & in the logarithm of inputEnergy in
    // the above segment.
    leptonDistributionExpansionTerm* const hkToMaxSq;
    // the term quadratic in inputEnergy in the above segment.

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
    W_minus_handed_muon( readierForNewPoint* const readierPointer,
                         CppSLHA::CppSLHA0 const* const spectrumData,
              CppSLHA::particle_property_set const* const firstParticle,
                         effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                     CppSLHA::particle_property_set const* const thirdParticle,
            CppSLHA::particle_property_set const* const fourthParticle,
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
    // ln( minimumEnergy ).
    double lnElk;
    // ln( Elk ).
    double lnEhk;
    // ln( Ehk ).
    double lnEmax;
    // ln( maximumEnergy ).

    double current_term_coefficient;
    // this is for holding the coefficientValue being calculated.

    segmentTermSet minToLkSegment;
    // the terms between minimumEnergy & Elk.
    leptonDistributionExpansionTerm* const minToLkInv;
    // the term linear in the inverse power of inputEnergy in the above
    // segment.
    leptonDistributionExpansionTerm* const minToLkConst;
    // the term constant with respect to inputEnergy in the above segment.
    leptonDistributionExpansionTerm* const minToLkLog;
    // the term linear in the logarithm of inputEnergy in the above segment.
    leptonDistributionExpansionTerm* const minToLkLin;
    // the term linear in inputEnergy in the above segment.
    leptonDistributionExpansionTerm* const minToLkLinlog;
    // the term linear in inputEnergy & in the logarithm of inputEnergy in
    // the above segment.
    leptonDistributionExpansionTerm* const minToLkSq;
    // the term quadratic in inputEnergy in the above segment.

    segmentTermSet lkToHkSegment;
    // the terms between Elk & Ehk.
    leptonDistributionExpansionTerm* const lkToHkConst;
    // the term constant with respect to inputEnergy in the above segment.
    leptonDistributionExpansionTerm* const lkToHkLin;
    // the term linear in inputEnergy in the above
    // segment.
    leptonDistributionExpansionTerm* const lkToHkSq;
    // the term quadratic in inputEnergy in the above segment.

    segmentTermSet hkToMaxSegment;
    // the terms between Ehk & maximumEnergy.
    leptonDistributionExpansionTerm* const hkToMaxInv;
    // the term linear in the inverse power of inputEnergy in the above
    // segment.
    leptonDistributionExpansionTerm* const hkToMaxConst;
    // the term constant with respect to inputEnergy in the above segment.
    leptonDistributionExpansionTerm* const hkToMaxLog;
    // the term linear in the logarithm of inputEnergy in the above segment.
    leptonDistributionExpansionTerm* const hkToMaxLin;
    // the term linear in inputEnergy in the above segment.
    leptonDistributionExpansionTerm* const hkToMaxLinlog;
    // the term linear in inputEnergy & in the logarithm of inputEnergy in
    // the above segment.
    leptonDistributionExpansionTerm* const hkToMaxSq;
    // the term quadratic in inputEnergy in the above segment.

    void
    calculateCoefficients()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };

  /* this derived class sets up the energy distribution for a light lepton
   * which is from the decay of an on-shell gauge vector boson from the
   * cascade decay of a squark to a lighter squark.
   */
  class vectorFromSquarkToMuon : public leptonEnergyDistribution
  {
  public:
    vectorFromSquarkToMuon( readierForNewPoint* const readierPointer,
                            CppSLHA::CppSLHA0 const* const spectrumData,
              CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
             CppSLHA::particle_property_set const* const thirdParticle )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    virtual
    ~vectorFromSquarkToMuon()
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
    // maximumEnergy - minimumEnergy

    segmentTermSet minToMaxSegment;
    // the terms between minimumEnergy & maximumEnergy.
    leptonDistributionExpansionTerm* const minToMaxConst;
    // the term constant with respect to inputEnergy in the above segment.
    leptonDistributionExpansionTerm* const minToMaxLin;
    // the term linear in inputEnergy in the above segment.
    leptonDistributionExpansionTerm* const minToMaxSq;
    // the term quadratic in inputEnergy in the above segment.

    void
    calculateCoefficients()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };

  // averaging over both jets from the Z boson produces the same distribution
  // as directly-produced light leptons.
  class zDirectJet : public zHandedMuon
  {
  public:
    zDirectJet( readierForNewPoint* const readierPointer,
                CppSLHA::CppSLHA0 const* const spectrumData,
                CppSLHA::particle_property_set const* const firstParticle,
                effectiveSquarkMassHolder* const effectiveSquarkMass,
                CppSLHA::particle_property_set const* const secondParticle,
                CppSLHA::particle_property_set const* const thirdParticle,
                CppSLHA::particle_property_set const* const fourthParticle )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    virtual
    ~zDirectJet()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  //protected:
    // nothing.
  };

  // averaging over both jets from the W boson produces the same distribution
  // which Z bosons produce.
  class wMinusDirectJet : public zHandedMuon
  {
  public:
    wMinusDirectJet( readierForNewPoint* const readierPointer,
                     CppSLHA::CppSLHA0 const* const spectrumData,
                     CppSLHA::particle_property_set const* const firstParticle,
                     effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                     CppSLHA::particle_property_set const* const thirdParticle,
                   CppSLHA::particle_property_set const* const fourthParticle )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    virtual
    ~wMinusDirectJet()
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
  class negativelyChargedHiggsMuon : public HiggsMuonPlusAntimuon
  {
  public:
    negativelyChargedHiggsMuon( readierForNewPoint* const readierPointer,
                                CppSLHA::CppSLHA0 const* const spectrumData,
                     CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                     CppSLHA::particle_property_set const* const thirdParticle,
                   CppSLHA::particle_property_set const* const fourthParticle )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    virtual
    ~negativelyChargedHiggsMuon()
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
  class negativelyChargedHiggsJet : public HiggsMuonPlusAntimuon
  {
  public:
    negativelyChargedHiggsJet( readierForNewPoint* const readierPointer,
                               CppSLHA::CppSLHA0 const* const spectrumData,
                     CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                     CppSLHA::particle_property_set const* const thirdParticle,
                   CppSLHA::particle_property_set const* const fourthParticle )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    virtual
    ~negativelyChargedHiggsJet()
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
  class neutralinoThreeBodyDecay : public leptonEnergyDistribution
  {
  public:
    neutralinoThreeBodyDecay( readierForNewPoint* const readierPointer,
                              CppSLHA::CppSLHA0 const* const spectrumData,
                     CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                     CppSLHA::particle_property_set const* const thirdParticle,
                      CppSLHA::particle_property_set const* const leftSfermion,
                    CppSLHA::particle_property_set const* const rightSfermion )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    virtual
    ~neutralinoThreeBodyDecay()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    CppSLHA::particle_property_set const* const leftSfermion;
    CppSLHA::particle_property_set const* const rightSfermion;

    segmentTermSet minToMaxSegment;
    // the terms between minimumEnergy & maximumEnergy.
    leptonDistributionExpansionTerm* const minToMaxConst;
    // the term constant with respect to inputEnergy in the above segment.

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
  class charginoThreeBodyDecay : public leptonEnergyDistribution
  {
  public:
    charginoThreeBodyDecay( readierForNewPoint* const readierPointer,
                            CppSLHA::CppSLHA0 const* const spectrumData,
                     CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                     CppSLHA::particle_property_set const* const thirdParticle,
             CppSLHA::particle_property_set const* const leftUpIsospinSfermion,
            CppSLHA::particle_property_set const* const rightUpIsospinSfermion,
           CppSLHA::particle_property_set const* const leftDownIsospinSfermion,
         CppSLHA::particle_property_set const* const rightDownIsospinSfermion )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    virtual
    ~charginoThreeBodyDecay()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    CppSLHA::particle_property_set const* const leftUpIsospinSfermion;
    CppSLHA::particle_property_set const* const rightUpIsospinSfermion;
    CppSLHA::particle_property_set const* const leftDownIsospinSfermion;
    CppSLHA::particle_property_set const* const rightDownIsospinSfermion;

    segmentTermSet minToMaxSegment;
    // the terms between minimumEnergy & maximumEnergy.
    leptonDistributionExpansionTerm* const minToMaxConst;
    // the term constant with respect to inputEnergy in the above segment.

    void
    calculateCoefficients()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };





  // inline functions:


  inline void
  flatNearMuonPlusAntimuon::calculateCoefficients()
  {
    productionFrameEnergy
    = ( ( ( secondMass * secondMass ) - ( thirdMass * thirdMass ) )
        / ( 2.0 * secondMass ) );
    minimumEnergy = ( ( secondMass / firstMass ) * productionFrameEnergy );
    maximumEnergy = ( ( firstMass / secondMass ) * productionFrameEnergy );

    // now we set up the coefficients of various terms:
    // actually, it just stays constant at 1.0.

    // MIN to MAX segment:
    minToMaxSegment.setSegmentRange( minimumEnergy,
                                     maximumEnergy );
  }



  inline void
  sameChiralityNearMuon::calculateCoefficients()
  {
    productionFrameEnergy
    = ( ( ( secondMass * secondMass ) - ( thirdMass * thirdMass ) )
        / ( 2.0 * secondMass ) );
    minimumEnergy
    = ( ( secondMass / firstMass ) * productionFrameEnergy );
    maximumEnergy
    = ( ( firstMass / secondMass ) * productionFrameEnergy );

    // now we set up the coefficients of various terms:

    // MIN to MAX segment:
    minToMaxSegment.setSegmentRange( minimumEnergy,
                                     maximumEnergy );
    minToMaxConst->setCoefficient( -minimumEnergy );
    minToMaxLin->setCoefficient( 1.0 );
  }



  inline void
  oppositeChiralityNearMuon::calculateCoefficients()
  {
    productionFrameEnergy
    = ( ( ( secondMass * secondMass ) - ( thirdMass * thirdMass ) )
        / ( 2.0 * secondMass ) );
    minimumEnergy
    = ( ( secondMass / firstMass ) * productionFrameEnergy );
    maximumEnergy
    = ( ( firstMass / secondMass ) * productionFrameEnergy );

    // now we set up the coefficients of various terms:

    // MIN to MAX segment:
    minToMaxSegment.setSegmentRange( minimumEnergy,
                                     maximumEnergy );
    minToMaxConst->setCoefficient( maximumEnergy );
    minToMaxLin->setCoefficient( -1.0 );
  }



  inline void
  vectorFromSquarkToMuon::calculateCoefficients()
  {
    mQhSq = ( firstMass * firstMass );
    mQlSq = ( secondMass * secondMass );
    mVSq = ( thirdMass * thirdMass );
    gammaV = ( ( mQhSq - mQlSq + mVSq ) / ( 2.0 * firstMass * thirdMass ) );
    betaV = sqrt( ( 1.0 - ( 1.0 / ( gammaV * gammaV ) ) ) );
    productionFrameEnergy = ( 0.5 * thirdMass );
    minimumEnergy = ( gammaV * ( 1.0 - betaV ) * productionFrameEnergy );
    maximumEnergy = ( gammaV * ( 1.0 + betaV ) * productionFrameEnergy );
    energyDifference = ( maximumEnergy - minimumEnergy );

    // now we set up the coefficients of various terms:
    // MIN to MAX segment:
    minToMaxSegment.setSegmentRange( minimumEnergy,
                                      maximumEnergy );
    minToMaxConst->setCoefficient( -minimumEnergy * maximumEnergy );
    minToMaxLin->setCoefficient( minimumEnergy + maximumEnergy );
    minToMaxSq->setCoefficient( -1.0 );
    //normalization
    //= ( ( energyDifference * energyDifference * energyDifference ) / 6.0 );
  }

}  // end of LHC_FASER namespace.

#endif /* LHC_FASER_DERIVED_LEPTON_DISTRIBUTIONS_HPP_ */
