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
 *      LHC-FASER_base_electroweak_cascade_stuff.hpp
 *      LHC-FASER_base_electroweak_cascade_stuff.cpp
 *      LHC-FASER_base_kinematics_stuff.hpp
 *      LHC-FASER_base_kinematics_stuff.cpp
 *      LHC-FASER_base_lepton_distribution_stuff.hpp
 *      LHC-FASER_base_lepton_distribution_stuff.cpp
 *      LHC-FASER_charged_electroweak_cascade_stuff.hpp
 *      LHC-FASER_charged_electroweak_cascade_stuff.cpp
 *      LHC-FASER_cross-section_stuff.hpp
 *      LHC-FASER_cross-section_stuff.cpp
 *      LHC-FASER_derived_lepton_distributions.hpp
 *      LHC-FASER_derived_lepton_distributions.cpp
 *      LHC-FASER_electroweak_cascade_collection_stuff.hpp
 *      LHC-FASER_electroweak_cascade_collection_stuff.cpp
 *      LHC-FASER_full_cascade_stuff.hpp
 *      LHC-FASER_full_cascade_stuff.cpp
 *      LHC-FASER_global_stuff.hpp
 *      LHC-FASER_global_stuff.cpp
 *      LHC-FASER_input_handling_stuff.hpp
 *      LHC-FASER_input_handling_stuff.cpp
 *      LHC-FASER_jet_kinematics_stuff.hpp
 *      LHC-FASER_jet_kinematics_stuff.cpp
 *      LHC-FASER_lepton_kinematics_stuff.hpp
 *      LHC-FASER_lepton_kinematics_stuff.cpp
 *      LHC-FASER_neutral_electroweak_cascade_stuff.hpp
 *      LHC-FASER_neutral_electroweak_cascade_stuff.cpp
 *      LHC-FASER_signal_calculator_stuff.hpp
 *      LHC-FASER_signal_calculator_stuff.cpp
 *      LHC-FASER_signal_data_collection_stuff.hpp
 *      LHC-FASER_signal_data_collection_stuff.cpp
 *      LHC-FASER_sparticle_decay_stuff.hpp
 *      LHC-FASER_sparticle_decay_stuff.cpp
 *      LHC-FASER_template_classes.hpp
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
                   CppSLHA::particle_property_set const* const thirdParticle );
    virtual
    ~flatNearMuonPlusAntimuon();

  protected:
    segmentTermSet minToMaxSegment;
    // the terms between minimumEnergy & maximumEnergy.
    leptonDistributionExpansionTerm* const minToMaxConst;
    // the term constant with respect to inputEnergy in the above segment.

    void
    calculateCoefficients();
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
                   CppSLHA::particle_property_set const* const thirdParticle );
    virtual
    ~sameChiralityNearMuon();

  protected:
    segmentTermSet minToMaxSegment;
    // the terms between minimumEnergy & maximumEnergy.
    leptonDistributionExpansionTerm* const minToMaxConst;
    // the term constant with respect to inputEnergy in the above segment.
    leptonDistributionExpansionTerm* const minToMaxLin;
    // the term linear in inputEnergy in the above segment.

    void
    calculateCoefficients();
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
                   CppSLHA::particle_property_set const* const thirdParticle );
    virtual
    ~oppositeChiralityNearMuon();

  protected:
    segmentTermSet minToMaxSegment;
    // the terms between minimumEnergy & maximumEnergy.
    leptonDistributionExpansionTerm* const minToMaxConst;
    // the term constant with respect to inputEnergy in the above segment.
    leptonDistributionExpansionTerm* const minToMaxLin;
    // the term linear in inputEnergy in the above segment.

    void
    calculateCoefficients();
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
                  CppSLHA::particle_property_set const* const fourthParticle );
    virtual
    ~flatFarMuonPlusAntimuon();

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
    calculateCoefficients();
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
                  CppSLHA::particle_property_set const* const fourthParticle );
    virtual
    ~sameChiralityFarMuon();

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
    calculateCoefficients();
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
                  CppSLHA::particle_property_set const* const fourthParticle );
    virtual
    ~oppositeChiralityFarMuon();


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
    calculateCoefficients();
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
                  CppSLHA::particle_property_set const* const fourthParticle );
      virtual
      ~HiggsMuonPlusAntimuon();

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
      calculateCoefficients();
  };

  /* this derived class to act as a base class for both Z & W^- boson-based
   * distributions. it is mainly so that I can modify the common elements (if
   * they ever need to be fixed again) without duplication.
   */
  class weakVectorBosonHandedMuon : public leptonEnergyDistribution
  {
  public:
    weakVectorBosonHandedMuon( readierForNewPoint* const readierPointer,
                               CppSLHA::CppSLHA0 const* const spectrumData,
                     CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                     CppSLHA::particle_property_set const* const thirdParticle,
                    CppSLHA::particle_property_set const* const fourthParticle,
                               bool const negativeMuonIsSameHandednessAsJet );
    virtual
    ~weakVectorBosonHandedMuon();

  protected:
    bool const negativeMuonIsSameHandednessAsJet;
    double Elk;
    double Ehk;
    bool cosinesLimitedByEnergy;

    /* these are for ease of calculating the coefficients (referring to the
     * squark as Q, the neutralino_2 as C, the neutralino_1 as X, the jet as j,
     * the Z boson as V, the negative muon as l, & the positive muon as v; it
     * makes more sense to know that the notation came from a calculation with
     * a chargino decaying to a W boson...):
     */
    double gammaCQ;
    // the gamma boost factor from going from the decaying electroweakino rest
    // frame to the squark rest frame.
    double gammaCQSq;
    // gammaCQ^2.
    double gammaCQCu;
    // gammaCQ^3.
    double gammaCQQu;
    // gammaCQ^4.
    double betaCQ;
    // the beta boost factor from going from the decaying electroweakino rest
    // frame to the squark rest frame.
    double betaCQSq;
    // betaCQ^2.
    double oneMinusBetaCQ;
    // ( 1 - betaCQ ).
    double oneMinusBetaCQSq;
    // oneMinusBetaCQ^2.
    double oneMinusBetaCQCu;
    // oneMinusBetaCQ^3.
    double oneMinusBetaCQQu;
    // oneMinusBetaCQ^4.
    double onePlusBetaCQ;
    double lnOnePlusBetaCQOverOneMinusBetaCQ;
    // ln( ( 1 + betaCQ ) / ( 1 - betaCQ ) ).
    double ElMin;
    // the minimum muon energy in the decaying electroweakino rest frame.
    double lnElMin;
    // ln( ElMin ).
    double ElMinSq;
    // ElMin^2.
    double ElMinCu;
    // ElMin^3.
    double ElMinQu;
    // ElMin^4.
    double mQ;
    // the mass of the squark.
    //double mQSq;
    // mQ^2.
    //double mQCu;
    // mQ^3.
    double mC;
    // the mass of the neutralino_2.
    double mCSq;
    // mC^2.
    double mCCu;
    // mC^3.
    //double mQCsqSum;
    // mQ^2 + mC^2.
    //double mQCsqDiff;
    // mQ^2 - mC^2.
    //double lnmQC;
    // ln( ( mQ / mC ) ).
    double mCX;
    // mC * the mass of the neutralino_1.
    double mXSq;
    // mX^2.
    double mVBSq;
    // the square of the mass of the weak vector boson.
    double mVBQu;
    // mVBSq^2.
    double EVB;
    // the energy of the Z boson in the rest frame of neutralino_2.
    //double EightmCEVElMin;
    // 8 * mC * EV * ElMin.
    double lnmVBSqOverFourElMinSq;
    // ln( ( mVBSq / ( 4 * ElMin^2 ) ).
    //double EightmCXElMinmQsq;
    // 8 * mC * mX * ElMin * mQ^2.
    //double lnEmin;
    // ln( minimumEnergy ).
    //double lnElk;
    // ln( Elk ).
    //double lnEhk;
    // ln( Ehk ).
    //double lnEmax;
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
    calculateEnergiesAndFactors();
    void
    calculateVvPlusJjAaAllSqSymCoefficients();
    void
    calculateTwiceVvAaSymCoefficients();
    void
    calculateVvSqMinusAaSqSymCoefficients();
    void
    calculateVvPlusJjAaAllSqAntiCoefficients();
    void
    calculateVvSqPlusAaSqAntiCoefficients();
    void
    calculateVvSqMinusAaSqAntiCoefficients();
    void
    setCurrentCoefficientsAsTotals();
    void
    addCurrentCoefficientsToTotals();
  };


  /* specific lepton or antilepton versions of the Z-channel distribution could
   * be done, but at the moment there is no need for them.
   */
  /* this derived class sets up the energy distribution for a light lepton
   * which is from the decay of an on-shell Z boson from the cascade decay of a
   * squark to a jet of either chirality plus a light lepton-antilepton pair
   * plus a lightest neutralino, keeping chirality effects.
   */
  class zHandedMuon : public weakVectorBosonHandedMuon
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
                 bool const shouldSumOverHandedness );
    virtual
    ~zHandedMuon();

  protected:
    bool const shouldNotSumOverHandedness;
    double sameHandednessFactor;
    double axialCouplingFactor;
    bool couplesAsAxialNotVector;

    void
    flipSignsOfCurrentCoefficients();

    void
    calculateCoefficients();
  };


  /* this derived class sets up the energy distribution for a light lepton
   * which is from the decay of an on-shell negatively-charged W boson from the
   * cascade decay of a down-type squark or an up-type antisquark.
   * unfortunately, the couplings are very model-specific, so this assumes just
   * the MSSM.
   */
  class wMinusHandedMuon : public weakVectorBosonHandedMuon
  {
  public:
    wMinusHandedMuon( readierForNewPoint* const readierPointer,
                      CppSLHA::CppSLHA0 const* const spectrumData,
                     CppSLHA::particle_property_set const* const firstParticle,
                      effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                     CppSLHA::particle_property_set const* const thirdParticle,
                    CppSLHA::particle_property_set const* const fourthParticle,
                      bool const jetIsLeftHanded );
    virtual
    ~wMinusHandedMuon();

  protected:
    int whichChargino;
    int whichNeutralino;
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

    void
    scaleCurrentCoefficients( double const scalingFactor );

    void
    calculateCoefficients();
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
             CppSLHA::particle_property_set const* const thirdParticle );
    virtual
    ~vectorFromSquarkToMuon();

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
    //double energyDifference;
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
    calculateCoefficients();
  };


  /* this derived class sets up the energy distribution for a light lepton
   * which is from the decay of an on-shell spin-0 boson from the
   * cascade decay of a squark to a lighter squark.
   */
  class scalarFromSquarkToMuon : public leptonEnergyDistribution
  {
  public:
    scalarFromSquarkToMuon( readierForNewPoint* const readierPointer,
                            CppSLHA::CppSLHA0 const* const spectrumData,
              CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
             CppSLHA::particle_property_set const* const thirdParticle );
    virtual
    ~scalarFromSquarkToMuon();

  protected:
    // these are for ease of calculating the coefficients (referring to the
    // decaying squark as Qh & the product squark as Ql):
    double mQhSq;
    // the square of the mass of the heavier squark.
    double mQlSq;
    // the square of the mass of the lighter squark.
    double mSSq;
    // the square of the mass of the spin-0 boson.
    double gammaS;
    // the gamma factor of the boost from the Qh rest frame to the spin-0 boson
    // rest frame.
    double betaS;
    // the beta factor of the boost from the Qh rest frame to the spin-0 boson
    // rest frame.

    segmentTermSet minToMaxSegment;
    // the terms between minimumEnergy & maximumEnergy.
    leptonDistributionExpansionTerm* const minToMaxConst;
    // the term constant with respect to inputEnergy in the above segment.

    void
    calculateCoefficients();
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
                CppSLHA::particle_property_set const* const fourthParticle );
    virtual
    ~zDirectJet();

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
                  CppSLHA::particle_property_set const* const fourthParticle );
    virtual
    ~wMinusDirectJet();

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
                  CppSLHA::particle_property_set const* const fourthParticle );
    virtual
    ~negativelyChargedHiggsMuon();

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
                  CppSLHA::particle_property_set const* const fourthParticle );
    virtual
    ~negativelyChargedHiggsJet();

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
                   CppSLHA::particle_property_set const* const rightSfermion );
    virtual
    ~neutralinoThreeBodyDecay();

  protected:
    CppSLHA::particle_property_set const* const leftSfermion;
    CppSLHA::particle_property_set const* const rightSfermion;

    segmentTermSet minToMaxSegment;
    // the terms between minimumEnergy & maximumEnergy.
    leptonDistributionExpansionTerm* const minToMaxConst;
    // the term constant with respect to inputEnergy in the above segment.

    void
    calculateCoefficients();
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
        CppSLHA::particle_property_set const* const rightDownIsospinSfermion );
    virtual
    ~charginoThreeBodyDecay();

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
    calculateCoefficients();
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
    //energyDifference = ( maximumEnergy - minimumEnergy );

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



  inline void
  scalarFromSquarkToMuon::calculateCoefficients()
  {
    mQhSq = ( firstMass * firstMass );
    mQlSq = ( secondMass * secondMass );
    mSSq = ( thirdMass * thirdMass );
    gammaS = ( ( mQhSq - mQlSq + mSSq ) / ( 2.0 * firstMass * thirdMass ) );
    betaS = sqrt( ( 1.0 - ( 1.0 / ( gammaS * gammaS ) ) ) );
    productionFrameEnergy = ( 0.5 * thirdMass );
    minimumEnergy = ( gammaS * ( 1.0 - betaS ) * productionFrameEnergy );
    maximumEnergy = ( gammaS * ( 1.0 + betaS ) * productionFrameEnergy );

    // now we set up the coefficients of various terms:
    // MIN to MAX segment:
    minToMaxSegment.setSegmentRange( minimumEnergy,
                                      maximumEnergy );
    minToMaxConst->setCoefficient( 1.0 );
    //normalization
    //= ( maximumEnergy - minimumEnergy );
  }



  inline void
  weakVectorBosonHandedMuon::setCurrentCoefficientsAsTotals()
  {
    minToLkInv->setCoefficient( currentMinToLkInvCoefficient );
    minToLkConst->setCoefficient( currentMinToLkConstCoefficient );
    minToLkLog->setCoefficient( currentMinToLkLogCoefficient );
    minToLkLin->setCoefficient( currentMinToLkLinCoefficient );
    minToLkLinlog->setCoefficient( currentMinToLkLinLogCoefficient );
    minToLkSq->setCoefficient( currentMinToLkSqCoefficient );

    lkToHkConst->setCoefficient( currentLkToHkConstCoefficient );
    lkToHkLin->setCoefficient( currentLkToHkLinCoefficient );
    lkToHkSq->setCoefficient( currentLkToHkSqCoefficient );

    hkToMaxInv->setCoefficient( currentHkToMaxInvCoefficient );
    hkToMaxConst->setCoefficient( currentHkToMaxConstCoefficient );
    hkToMaxLog->setCoefficient( currentHkToMaxLogCoefficient );
    hkToMaxLin->setCoefficient( currentHkToMaxLinCoefficient );
    hkToMaxLinlog->setCoefficient( currentHkToMaxLinLogCoefficient );
    hkToMaxSq->setCoefficient( currentHkToMaxSqCoefficient );
  }

  inline void
  weakVectorBosonHandedMuon::addCurrentCoefficientsToTotals()
  {
    minToLkInv->addToCoefficient( currentMinToLkInvCoefficient );
    minToLkConst->addToCoefficient( currentMinToLkConstCoefficient );
    minToLkLog->addToCoefficient( currentMinToLkLogCoefficient );
    minToLkLin->addToCoefficient( currentMinToLkLinCoefficient );
    minToLkLinlog->addToCoefficient( currentMinToLkLinLogCoefficient );
    minToLkSq->addToCoefficient( currentMinToLkSqCoefficient );

    lkToHkConst->addToCoefficient( currentLkToHkConstCoefficient );
    lkToHkLin->addToCoefficient( currentLkToHkLinCoefficient );
    lkToHkSq->addToCoefficient( currentLkToHkSqCoefficient );

    hkToMaxInv->addToCoefficient( currentHkToMaxInvCoefficient );
    hkToMaxConst->addToCoefficient( currentHkToMaxConstCoefficient );
    hkToMaxLog->addToCoefficient( currentHkToMaxLogCoefficient );
    hkToMaxLin->addToCoefficient( currentHkToMaxLinCoefficient );
    hkToMaxLinlog->addToCoefficient( currentHkToMaxLinLogCoefficient );
    hkToMaxSq->addToCoefficient( currentHkToMaxSqCoefficient );
  }



  inline void
  zHandedMuon::flipSignsOfCurrentCoefficients()
  {
    currentMinToLkInvCoefficient = -currentMinToLkInvCoefficient;
    currentMinToLkConstCoefficient = -currentMinToLkConstCoefficient;
    currentMinToLkLogCoefficient = -currentMinToLkLogCoefficient;
    currentMinToLkLinCoefficient = -currentMinToLkLinCoefficient;
    currentMinToLkLinLogCoefficient = -currentMinToLkLinLogCoefficient;
    currentMinToLkSqCoefficient = -currentMinToLkSqCoefficient;
    currentLkToHkConstCoefficient = -currentLkToHkConstCoefficient;
    currentLkToHkLinCoefficient = -currentLkToHkLinCoefficient;
    currentLkToHkSqCoefficient = -currentLkToHkSqCoefficient;
    currentHkToMaxInvCoefficient = -currentHkToMaxInvCoefficient;
    currentHkToMaxConstCoefficient = -currentHkToMaxConstCoefficient;
    currentHkToMaxLogCoefficient = -currentHkToMaxLogCoefficient;
    currentHkToMaxLinCoefficient = -currentHkToMaxLinCoefficient;
    currentHkToMaxLinLogCoefficient = -currentHkToMaxLinLogCoefficient;
    currentHkToMaxSqCoefficient = -currentHkToMaxSqCoefficient;
  }



  inline void
  wMinusHandedMuon::scaleCurrentCoefficients( double const scalingFactor )
  {
    currentMinToLkInvCoefficient *= scalingFactor;
    currentMinToLkConstCoefficient *= scalingFactor;
    currentMinToLkLogCoefficient *= scalingFactor;
    currentMinToLkLinCoefficient *= scalingFactor;
    currentMinToLkLinLogCoefficient *= scalingFactor;
    currentMinToLkSqCoefficient *= scalingFactor;
    currentLkToHkConstCoefficient *= scalingFactor;
    currentLkToHkLinCoefficient *= scalingFactor;
    currentLkToHkSqCoefficient *= scalingFactor;
    currentHkToMaxInvCoefficient *= scalingFactor;
    currentHkToMaxConstCoefficient *= scalingFactor;
    currentHkToMaxLogCoefficient *= scalingFactor;
    currentHkToMaxLinCoefficient *= scalingFactor;
    currentHkToMaxLinLogCoefficient *= scalingFactor;
    currentHkToMaxSqCoefficient *= scalingFactor;
  }

}  // end of LHC_FASER namespace.

#endif /* LHC_FASER_DERIVED_LEPTON_DISTRIBUTIONS_HPP_ */
