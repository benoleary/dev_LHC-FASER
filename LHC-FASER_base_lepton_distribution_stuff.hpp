/*
 * LHC-FASER_base_lepton_distribution_stuff.hpp
 *
 *  Created on: 17 Mar 2010
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

/* this is a set of classes used for the calculation of the fraction of leptons
 * that pass various cuts at the moment. it also includes jets by pretending
 * that they're basically leptons.
 */


#ifndef LHC_FASER_BASE_LEPTON_DISTRIBUTION_STUFF_HPP
#define LHC_FASER_BASE_LEPTON_DISTRIBUTION_STUFF_HPP

#include "CppSLHA/CppSLHA.hpp"
#include "LHC-FASER_global_stuff.hpp"

namespace LHC_FASER
{
  // this class is to make it easier to update the effective squark mass.
  class effectiveSquarkMassHolder
  {
  public:
    effectiveSquarkMassHolder();
    virtual
    ~effectiveSquarkMassHolder();

    virtual double
    getEffectiveSquarkMass()
    = 0;

  //protected:
    // nothing
  };


  // this derived class is for the case of on-shell squarks:
  class onShellSquarkMassHolder : public effectiveSquarkMassHolder
  {
  public:
    onShellSquarkMassHolder(
                   CppSLHA::particle_property_set const* const onShellSquark );
    virtual
    ~onShellSquarkMassHolder();

    virtual double
    getEffectiveSquarkMass();
    CppSLHA::particle_property_set const*
    getOnShellSquark()
    const;

  protected:
    CppSLHA::particle_property_set const* onShellSquark;
  };


  // this derived class is for the case of splitting an on-shell squark's mass
  // between an on-shell boson and a virtual squark:
  class squarkMinusBosonMassHolder : public effectiveSquarkMassHolder
  {
  public:
    squarkMinusBosonMassHolder(
                     CppSLHA::particle_property_set const* const onShellSquark,
                      CppSLHA::particle_property_set const* const onShellBoson,
                    CppSLHA::particle_property_set const* const onShellEwino );
    virtual
    ~squarkMinusBosonMassHolder();

    virtual double
    getEffectiveSquarkMass();
    bool
    isEquivalent( CppSLHA::particle_property_set const* const onShellSquark,
                  CppSLHA::particle_property_set const* const onShellBoson,
                  CppSLHA::particle_property_set const* const onShellEwino )
    const;

  protected:
    CppSLHA::particle_property_set const* onShellSquark;
    CppSLHA::particle_property_set const* onShellBoson;
    CppSLHA::particle_property_set const* onShellEwino;
  };


  /* this derived class is for the case of fudging splitting an on-shell
   * squark's mass between an on-shell boson and a virtual squark, using
   * instead the mass it should have to have approximately the right energy for
   * the boson for the case of the virtual squark being on-shell (meant to be
   * used along with giving the on-shell squark in place of the virtual squark
   * for the decay to an electroweakino, so that the masses work out):
   */
  class squarkPlusBosonMassHolder : public effectiveSquarkMassHolder
  {
  public:
    squarkPlusBosonMassHolder(
                     CppSLHA::particle_property_set const* const virtualSquark,
                      CppSLHA::particle_property_set const* const onShellBoson,
                    CppSLHA::particle_property_set const* const onShellEwino );
    virtual
    ~squarkPlusBosonMassHolder();

    virtual double
    getEffectiveSquarkMass();
    bool
    isEquivalent( CppSLHA::particle_property_set const* const virtualSquark,
                  CppSLHA::particle_property_set const* const onShellBoson,
                  CppSLHA::particle_property_set const* const onShellEwino )
    const;

  protected:
    CppSLHA::particle_property_set const* virtualSquark;
    CppSLHA::particle_property_set const* onShellBoson;
    CppSLHA::particle_property_set const* onShellEwino;
  };


  // this class is for emulating doing minor symbolic algebra.
  class leptonDistributionExpansionTerm
  {
  public:
    static double const minimumInputEnergy;
    // this is used to prevent trying to take logarithms of 0.0 or less.

    leptonDistributionExpansionTerm( int const powerOfEnergy,
                                     int const powerOfLogarithm,
                                     double const coefficientValue,
                 leptonDistributionExpansionTerm* const referenceTerm = NULL );
    virtual
    ~leptonDistributionExpansionTerm();

    leptonDistributionExpansionTerm const*
    getReferenceTerm()
    const;
    int
    getPowerOfEnergy()
    const;
    int
    getPowerOfLogarithm()
    const;
    double
    getCoefficient()
    const;
    void
    setCoefficient( double const coefficientValue );
    void
    addToCoefficient( double const coefficientValue );
    void
    normalizeCoefficient( double const normalizingFactor );
    // this just multiplies coefficientValue by normalizingFactor.
    virtual double
    evaluate( double const inputEnergy )
    const;
    // this evaluates this term for the given input energy.
    virtual double
    getArea( double startEnergy,
             double endEnergy )
    const;
    // this gives the definite integral of the term from startEnergy to
    // endEnergy.

  protected:
    int powerOfEnergy;
    int powerOfLogarithm;
    double coefficientValue;
    leptonDistributionExpansionTerm* const referenceTerm;
    // this is for ease of updating coefficients in tauon-decay distributions.
  };

  class leptonDistributionInverseTerm : public leptonDistributionExpansionTerm
  {
  public:
    leptonDistributionInverseTerm( double const coefficientValue,
                 leptonDistributionExpansionTerm* const referenceTerm = NULL );
    virtual
    ~leptonDistributionInverseTerm();

    virtual double
    evaluate( double const inputEnergy )
    const;
    // this evaluates this term for the given input energy.
    virtual double
    getArea( double startEnergy,
             double endEnergy )
    const;
    // this gives the definite integral of the term from startEnergy to
    // endEnergy.


  //protected:
    // nothing.
  };

  class leptonDistributionConstantTerm : public leptonDistributionExpansionTerm
  {
  public:
    leptonDistributionConstantTerm( double const coefficientValue,
                 leptonDistributionExpansionTerm* const referenceTerm = NULL );
    virtual
    ~leptonDistributionConstantTerm();

    virtual double
    evaluate( double const inputEnergy )
    const;
    // this evaluates this term for the given input energy.
    virtual double
    getArea( double startEnergy,
             double endEnergy )
    const;
    // this gives the definite integral of the term from startEnergy to
    // endEnergy.


  //protected:
    // nothing.
  };

  class leptonDistributionPowerTerm : public leptonDistributionExpansionTerm
  {
  public:
    leptonDistributionPowerTerm( int const powerOfEnergy,
                                 double const coefficientValue,
                 leptonDistributionExpansionTerm* const referenceTerm = NULL );
    virtual
    ~leptonDistributionPowerTerm();

    virtual double
    evaluate( double const inputEnergy )
    const;
    // this evaluates this term for the given input energy.
    virtual double
    getArea( double startEnergy,
             double endEnergy )
    const;
    // this gives the definite integral of the term from startEnergy to
    // endEnergy.


  //protected:
    // nothing.
  };


  // this is a class to hold a set of leptonDistributionExpansionTerm
  // objects relevant for an energy range for the distribution.
  class segmentTermSet
  {
  public:
    segmentTermSet();
    ~segmentTermSet();

    double
    getSegmentStart()
    const;
    void
    setSegmentStart( double const inputValue );
    double
    getSegmentEnd()
    const;
    void
    setSegmentEnd( double const inputValue );
    void
    setSegmentRange( double const rangeStart,
                     double const rangeEnd );
    std::vector< leptonDistributionExpansionTerm* > const*
    getTerms()
    const;
    int
    getNumberOfTerms()
    const;
    leptonDistributionExpansionTerm*
    addTerm( int const powerOfEnergy,
             int const powerOfLogarithm,
             leptonDistributionExpansionTerm* const referenceTerm = NULL,
             double const
             coefficientValue = CppSLHA::CppSLHA_global::really_wrong_value );
    void
    normalizeCoefficients( double const normalizingFactor );
    // this goes through each term & gets it to multiply its coefficient
    // by normalizingFactor.
    double
    getCoefficient( int const powerOfEnergy,
                    int const powerOfLogarithm )
    const;
    double
    evaluate( double const inputEnergy )
    const;
    // this evaluates the sum of the set of terms for the given input energy.
    double
    getArea()
    const;
    // this evaluates the sum of the areas of the set of terms for the
    // segment's range.


  protected:
    std::vector< leptonDistributionExpansionTerm* > segmentTerms;
    // the segmentTermSet is responsible for calling the destructors of the
    // leptonDistributionExpansionTerms held in segmentTerms.
    leptonDistributionExpansionTerm* termAdder;
    double segmentStart;
    double segmentEnd;
  };



  /* this class is the base abstract class for the functors for obtaining
   * energy distributions for the visible decay productSet of tau leptons with
   * given energy distributions.
   * (yes, left-handed tau leptons / right-handed tau antileptons give harder
   * muons / electrons but softer pions than the other helicity.)
   */
  class tauDecayCoefficient
  {
  public:
    tauDecayCoefficient();
    virtual
    ~tauDecayCoefficient();

    virtual double
    operator()( int const visibleProductEnergyPower
                 /* the power of the visible product's energy for the term. */,
                int const visibleProductLogPower
/* the power of the logarithm of the visible product's energy for the term. */,
                int const tauLeptonEnergyPower
                /* the power of the tau lepton's energy for the term. */,
                int const tauLeptonLogPower
     /* the power of the logarithm of the tau lepton's energy for the term. */,
                double tauMinEnergy
         /* the minimum of the range of the tau lepton energies considered. */,
                double const tauMaxEnergy
         /* the maximum of the range of the tau lepton energies considered. */,
                bool const isInsideRange )
    const
    = 0;
    /* this is overridden by functions which return the value of the
     * coefficient for the specified term for the tauon decay product from the
     * specified term for the tauon.
     */

  //protected:
    // nothing.
  };


  /* this functor returns the coefficient for the muon distribution which would
   * come from a left-handed tau lepton distribution of the given power of the
   * tau lepton's energy & power of logarithm thereof, for requested integer
   * powers of the muon's energy or logarithm thereof, which also depends on
   * whether the muon's energy is inside the range of the tau lepton's energy
   * for this segment.
   */
  class hardMuonFromTau : public tauDecayCoefficient
  {
  public:
    hardMuonFromTau();
    virtual
    ~hardMuonFromTau();

    virtual double
    operator()( int const visibleProductEnergyPower
                /* the power of the visible product's energy for the term. */,
                int const visibleProductLogPower
/* the power of the logarithm of the visible product's energy for the term. */,
                int const tauLeptonEnergyPower
                /* the power of the tau lepton's energy for the term. */,
                int const tauLeptonLogPower
    /* the power of the logarithm of the tau lepton's energy for the term. */,
                double tauMinEnergy
         /* the minimum of the range of the tau lepton energies considered. */,
                double const tauMaxEnergy
         /* the maximum of the range of the tau lepton energies considered. */,
                bool const isInsideRange )
    const;


  //protected:
    // nothing.
  };


  /* this functor returns the coefficient for the muon distribution which would
   * come from a right-handed tau lepton distribution of the given power of the
   * tau lepton's energy & power of logarithm thereof, for requested integer
   * powers of the muon's energy or logarithm thereof, which also depends on
   * whether the muon's energy is inside the range of the tau lepton's energy
   * for this segment.
   */
  class softMuonFromTau : public tauDecayCoefficient
  {
  public:
    softMuonFromTau();
    virtual
    ~softMuonFromTau();

    virtual double
    operator()( int const visibleProductEnergyPower
                /* the power of the visible product's energy for the term. */,
                int const visibleProductLogPower
/* the power of the logarithm of the visible product's energy for the term. */,
                int const tauLeptonEnergyPower
                /* the power of the tau lepton's energy for the term. */,
                int const tauLeptonLogPower
     /* the power of the logarithm of the tau lepton's energy for the term. */,
                double tauMinEnergy
         /* the minimum of the range of the tau lepton energies considered. */,
                double const tauMaxEnergy
         /* the maximum of the range of the tau lepton energies considered. */,
                bool const isInsideRange )
    const;


  //protected:
    // nothing.
  };


  /* this functor returns the coefficient for the pion distribution which would
   * come from a right-handed tau lepton distribution of the given power of the
   * tau lepton's energy & power of logarithm thereof, for requested integer
   * powers of the pion's energy or logarithm thereof, which also depends on
   * whether the pion's energy is inside the range of the tau lepton's energy
   * for this segment.
   */
  class hardPionFromTau : public tauDecayCoefficient
  {
  public:
    hardPionFromTau();
    virtual
    ~hardPionFromTau();

    virtual double
    operator()( int const visibleProductEnergyPower
                /* the power of the visible product's energy for the term. */,
                int const visibleProductLogPower
/* the power of the logarithm of the visible product's energy for the term. */,
                int const tauLeptonEnergyPower
                /* the power of the tau lepton's energy for the term. */,
                int const tauLeptonLogPower
     /* the power of the logarithm of the tau lepton's energy for the term. */,
                double tauMinEnergy
         /* the minimum of the range of the tau lepton energies considered. */,
                double const tauMaxEnergy
         /* the maximum of the range of the tau lepton energies considered. */,
                bool const isInsideRange )
    const;


  //protected:
    // nothing.
  };


  /* this functor returns the coefficient for the pion distribution which would
   * come from a left-handed tau lepton distribution of the given power of the
   * tau lepton's energy & power of logarithm thereof, for requested integer
   * powers of the pion's energy or logarithm thereof, which also depends on
   * whether the pion's energy is inside the range of the tau lepton's energy
   * for this segment.
   */
  class softPionFromTau : public tauDecayCoefficient
  {
  public:
    softPionFromTau();
    virtual
    ~softPionFromTau();

    virtual double
    operator()( int const visibleProductEnergyPower
                /* the power of the visible product's energy for the term. */,
                int const visibleProductLogPower
/* the power of the logarithm of the visible product's energy for the term. */,
                int const tauLeptonEnergyPower
                /* the power of the tau lepton's energy for the term. */,
                int const tauLeptonLogPower
     /* the power of the logarithm of the tau lepton's energy for the term. */,
                double tauMinEnergy
         /* the minimum of the range of the tau lepton energies considered. */,
                double const tauMaxEnergy
         /* the maximum of the range of the tau lepton energies considered. */,
                bool const isInsideRange )
    const;


  //protected:
    // nothing.
  };


  /* this holds together a pair of segmentTermSet pointers for the 2 parts of
   * of the distribution for the visible decay product of a tau lepton with
   * distribution given by the 3rd segmentTermSet pointer the
   * tauSegmentTriple has. it performs the role of updating the visible decay
   * product terms based on the reference terms.
   */
  class tauSegmentTriple
  {
  public:
    tauSegmentTriple( std::vector< segmentTermSet* >* segmentSetToPopulate,
                      segmentTermSet* referenceSegment,
                      tauDecayCoefficient const* const tauDecay );
    ~tauSegmentTriple();

    void
    updateSegments();

  protected:
    segmentTermSet* referenceSegment;
    segmentTermSet* lowSegment;
    segmentTermSet* highSegment;
    tauDecayCoefficient const* const tauDecay;
  };


  /* this class is just an abstract base class for some polymorphic derived
   * classes. However, it does do most of the work, the derived classes just
   * set up the leptonDistributionExpansionTerm objects properly.
   */
  class leptonEnergyDistribution : public getsReadiedForNewPoint
  {
  public:
    leptonEnergyDistribution( readierForNewPoint* const readierPointer,
                              CppSLHA::CppSLHA0 const* const spectrumPointer,
                     CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                     CppSLHA::particle_property_set const* const thirdParticle,
                  CppSLHA::particle_property_set const* const fourthParticle );
    virtual
    ~leptonEnergyDistribution();

    double
    getMinimumEnergy();
    double
    getMaximumEnergy();
    std::vector< segmentTermSet* >*
    getSegments();
    std::vector< segmentTermSet* > const*
    inspectSegments()
    const;
    double
    valueAt( double const inputEnergy );
    // this returns the value of the distribution at the given value, with
    // arbitrary normalization (different for each distribution).


  protected:
    CppSLHA::CppSLHA0 const* spectrumPointer;
    // this holds all the information about the MSSM point. care must be taken
    // to ensure that it has the appropriate BLOCKs for the distribution.
    double minimumEnergy;
    double maximumEnergy;

    std::vector< segmentTermSet* > segmentSet;

    CppSLHA::particle_property_set const* const firstParticle;
    effectiveSquarkMassHolder* const effectiveSquarkMass;
    CppSLHA::particle_property_set const* const secondParticle;
    CppSLHA::particle_property_set const* const thirdParticle;
    CppSLHA::particle_property_set const* const fourthParticle;
    /* these masses could be pulled out of spectrumData every time they are
     * needed, but it's probably worthwhile to have them as quickly-accessible
     * doubles:
     */
    double firstMass;
    // the squark mass in a standard decay, but let's pretend that I'm being
    // general.
    double secondMass;
    // the intermediate neutralino mass in the standard decay, but could be
    // e.g. the chargino mass in a certain decay.
    double thirdMass;
    // the slepton mass in the standard decay, but could be e.g.
    // the Z boson mass or sneutrino mass.
    double fourthMass;
    // the lightest neutralino mass in the standard decay, but could be
    // something else in some crazy cascade.
    double normalizingDivisor;
    // this is *an overall denominator*, rather than a multiplicative
    // normalization.
    double normalizingFactor;
    // this *is* an overall multiplicative normalization
    // (= 1.0 / normalizingDivisor ).
    // this is not strictly necessary, but is handy for distributions from
    // on-shell cascades:
    double productionFrameEnergy;

    void
    reset();
    // this updates the masses then calls calculate_coefficients().
    virtual void
    calculateCoefficients()
    = 0;
    // this is overridden by functions which calculate & set the coefficients
    // of all the terms of the distribution.
    void
    normalizeCoefficients();
    // this goes through each segment & gets it to divide its terms'
    // coefficients by normalizingDivisor.
  };


  /* this is a derived class which sets up the energy distribution for a muon
   * or pion which is from the decay of a left- or right-handed tau lepton or
   * antilepton with the provided energy distribution.
   */
  class visibleTauDecayProduct : public leptonEnergyDistribution
  {
  public:
    visibleTauDecayProduct( readierForNewPoint* const readierPointer,
                            leptonEnergyDistribution* const tauDistribution,
                            tauDecayCoefficient const* const tauDecay );
    virtual
    ~visibleTauDecayProduct();


  protected:
    leptonEnergyDistribution* const tauDistribution;
    std::vector< tauSegmentTriple* > tauTriples;

    virtual void
    calculateCoefficients();
    // this is overridden by functions which calculate & set the coefficients
    // of all the terms of the distribution.
  };




  // inline functions


  inline double
  onShellSquarkMassHolder::getEffectiveSquarkMass()
  {
    return onShellSquark->get_absolute_mass();
  }

  inline CppSLHA::particle_property_set const*
  onShellSquarkMassHolder::getOnShellSquark()
  const
  {
    return onShellSquark;
  }



  inline double
  squarkMinusBosonMassHolder::getEffectiveSquarkMass()
  {
    // return the on-shell squark's mass, minus the boson's mass, plus a third
    // of the energy from the rest of the squark's mass.
    return ( ( 4.0 * onShellSquark->get_absolute_mass()
               - 4.0 * onShellBoson->get_absolute_mass()
               - onShellEwino->get_absolute_mass() ) / 3.0 );
  }

  inline bool
  squarkMinusBosonMassHolder::isEquivalent(
                     CppSLHA::particle_property_set const* const onShellSquark,
                      CppSLHA::particle_property_set const* const onShellBoson,
                     CppSLHA::particle_property_set const* const onShellEwino )
  const
  {
    if( ( onShellSquark == this->onShellSquark )
        &&
        ( onShellBoson == this->onShellBoson )
        &&
        ( onShellEwino == this->onShellEwino ) )
    {
      return true;
    }
    else
    {
      return false;
    }
  }



  inline double
  squarkPlusBosonMassHolder::getEffectiveSquarkMass()
  {
    // return the on-shell squark's mass, plus the boson's mass, plus a third
    // of the energy from the rest of the squark's mass.
    return ( ( 4.0 * virtualSquark->get_absolute_mass()
               + 2.0 * onShellBoson->get_absolute_mass()
               - onShellEwino->get_absolute_mass() ) / 3.0 );
  }

  inline bool
  squarkPlusBosonMassHolder::isEquivalent(
                     CppSLHA::particle_property_set const* const virtualSquark,
                      CppSLHA::particle_property_set const* const onShellBoson,
                     CppSLHA::particle_property_set const* const onShellEwino )
  const
  {
    if( ( virtualSquark == this->virtualSquark )
        &&
        ( onShellBoson == this->onShellBoson )
        &&
        ( onShellEwino == this->onShellEwino ) )
    {
      return true;
    }
    else
    {
      return false;
    }
  }



  inline leptonDistributionExpansionTerm const*
  leptonDistributionExpansionTerm::getReferenceTerm()
  const
  {
    return referenceTerm;
  }

  inline int
  leptonDistributionExpansionTerm::getPowerOfEnergy()
  const
  {
    return powerOfEnergy;
  }

  inline int
  leptonDistributionExpansionTerm::getPowerOfLogarithm()
  const
  {
    return powerOfLogarithm;
  }

  inline double
  leptonDistributionExpansionTerm::getCoefficient()
  const
  {
    return coefficientValue;
  }

  inline void
  leptonDistributionExpansionTerm::setCoefficient(
                                                double const coefficientValue )
  {
     this->coefficientValue = coefficientValue;
  }

  inline void
  leptonDistributionExpansionTerm::addToCoefficient(
                                                double const coefficientValue )
  {
     this->coefficientValue += coefficientValue;
  }

  inline void
  leptonDistributionExpansionTerm::normalizeCoefficient(
                                               double const normalizingFactor )
  // this just multiplies coefficientValue by normalizingFactor.
  {
    coefficientValue *= normalizingFactor;
  }



  inline double
  leptonDistributionInverseTerm::evaluate( double const inputEnergy )
  const
  // this evaluates this term for the given input energy.
  {
    return ( coefficientValue / inputEnergy );
  }

  inline double
  leptonDistributionInverseTerm::getArea( double startEnergy,
                                          double endEnergy )
  const
  // this gives the definite integral of the term from startEnergy to
  // endEnergy.
  {
    double returnValue( 0.0 );
    if( 0.0 != coefficientValue )
    {
      if( minimumInputEnergy > startEnergy )
      {
        startEnergy = minimumInputEnergy;
      }

      if( endEnergy > startEnergy )
      {
        returnValue
        = ( coefficientValue * log( ( endEnergy / startEnergy ) ) );
      }
    }
    return returnValue;
  }



  inline double
  leptonDistributionConstantTerm::evaluate( double const inputEnergy )
  const
  // this evaluates this term for the given input energy.
  {
    return coefficientValue;
  }

  inline double
  leptonDistributionConstantTerm::getArea( double startEnergy,
                                           double endEnergy )
  const
  // this gives the definite integral of the term from startEnergy to
  // endEnergy.
  {
    double returnValue( 0.0 );
    if( 0.0 != coefficientValue )
    {
      if( minimumInputEnergy > startEnergy )
      {
        startEnergy = minimumInputEnergy;
      }

      if( endEnergy > startEnergy )
      {
        returnValue
        = ( coefficientValue * ( endEnergy - startEnergy ) );
      }
    }
    return returnValue;
  }



  inline double
  leptonDistributionPowerTerm::evaluate( double const inputEnergy )
  const
  // this evaluates this term for the given input energy.
  {
    if( 0.0 != coefficientValue )
    {
      double returnValue( coefficientValue );

      for( int powerCount( powerOfEnergy );
           0 < powerCount;
           --powerCount )
      {
        returnValue *= inputEnergy;
      }
      return returnValue;
    }
    else
    {
      return 0.0;
    }
  }



  inline double
  segmentTermSet::getSegmentStart()
  const
  {
    return segmentStart;
  }

  inline void
  segmentTermSet::setSegmentStart( double const inputValue )
  {
    segmentStart = inputValue;
  }

  inline double
  segmentTermSet::getSegmentEnd()
  const
  {
    return segmentEnd;
  }

  inline void
  segmentTermSet::setSegmentEnd( double const inputValue )
  {
    segmentEnd = inputValue;
  }

  inline void
  segmentTermSet::setSegmentRange( double const rangeStart,
                                   double const rangeEnd )
  {
    segmentStart = rangeStart;
    segmentEnd = rangeEnd;
  }

  inline std::vector< leptonDistributionExpansionTerm* > const*
  segmentTermSet::getTerms()
  const
  {
    return &segmentTerms;
  }

  inline int
  segmentTermSet::getNumberOfTerms()
  const
  {
    return (int)(segmentTerms.size());
  }

  inline void
  segmentTermSet::normalizeCoefficients( double const normalizingFactor )
  // this goes through each term & gets it to multiply its coefficientValue by
  // normalizingFactor.
  {
    for( std::vector< leptonDistributionExpansionTerm* >::iterator
         termIterator( segmentTerms.begin() );
         segmentTerms.end() > termIterator;
         ++termIterator )
    {
      (*termIterator)->normalizeCoefficient( normalizingFactor );
    }
  }

  inline double
  segmentTermSet::getCoefficient( int const powerOfEnergy,
                                  int const powerOfLogarithm )
  const
  {
    double returnValue( 0.0 );
    // start by assuming that the coefficientValue is zero.
    for( std::vector< leptonDistributionExpansionTerm* >::const_iterator
         termIterator( segmentTerms.begin() );
         segmentTerms.end() > termIterator;
         ++termIterator )
      // go through all the terms...
      {
        if( ( powerOfEnergy == (*termIterator)->getPowerOfEnergy() )
            &&
            ( powerOfLogarithm == (*termIterator)->getPowerOfLogarithm() ) )
          // if we find the requested term with a non-zero coefficientValue...
        {
          returnValue = (*termIterator)->getCoefficient();
          // note the coefficientValue to return.
          termIterator = segmentTerms.end();
          // stop looking.
        }
      }
    return returnValue;
  }

  inline double
  segmentTermSet::getArea()
  const
  // this evaluates the sum of the areas of the set of terms for the segment's
  // range.
  {
    double returnValue( 0.0 );
    for( std::vector< leptonDistributionExpansionTerm* >::const_iterator
        summationIterator( segmentTerms.begin() );
        summationIterator < segmentTerms.end();
        ++summationIterator )
    {
      returnValue += (*summationIterator)->getArea( segmentStart,
                                                    segmentEnd );
    }
    return returnValue;
  }



  inline double
  leptonEnergyDistribution::valueAt( double const inputEnergy )
  // this returns the value of the distribution at the given value, with
  // arbitrary normalization (different for each distribution).
  {
    if( needsToPrepareForThisPoint() )
    {
      reset();
    }
    double returnValue( 0.0 );
    for( std::vector< segmentTermSet* >::const_iterator
         summationIterator( segmentSet.begin() );
         segmentSet.end() > summationIterator;
         ++summationIterator )
    {
      returnValue += (*summationIterator)->evaluate( inputEnergy );
    }
    return returnValue;
  }

  inline double
  leptonEnergyDistribution::getMinimumEnergy()
  {
    if( needsToPrepareForThisPoint() )
    {
      reset();
    }
    return minimumEnergy;
  }

  inline double
  leptonEnergyDistribution::getMaximumEnergy()
  {
    if( needsToPrepareForThisPoint() )
    {
      reset();
    }
    return maximumEnergy;
  }

  inline std::vector< segmentTermSet* >*
  leptonEnergyDistribution::getSegments()
  {
    return &segmentSet;
  }

  inline std::vector< segmentTermSet* > const*
  leptonEnergyDistribution::inspectSegments()
  const
  {
    return &segmentSet;
  }

  inline void
  leptonEnergyDistribution::normalizeCoefficients()
  // this goes through each segment & gets it to divide its terms' coefficients
  // by normalizingDivisor.
  {
    normalizingDivisor = 0.0;
    for( std::vector< segmentTermSet* >::iterator
         segmentIterator( segmentSet.begin() );
         segmentSet.end() > segmentIterator;
         ++segmentIterator )
    {
      normalizingDivisor += (*segmentIterator)->getArea();
    }
    normalizingFactor = ( 1.0 / normalizingDivisor );
    for( std::vector< segmentTermSet* >::iterator
         segmentIterator( segmentSet.begin() );
         segmentSet.end() > segmentIterator;
         ++segmentIterator )
    {
      (*segmentIterator)->normalizeCoefficients( normalizingFactor );
    }
  }

}  // end of LHC_FASER namespace.

#endif  // LHC_FASER_BASE_LEPTON_DISTRIBUTION_STUFF_HPP
