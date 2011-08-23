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
    effectiveSquarkMassHolder()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~effectiveSquarkMassHolder()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual double
    getEffectiveSquarkMass()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  //protected:
    // nothing
  };

  // this class is for emulating doing minor symbolic algebra.
  class lepton_distribution_expansion_term
  {

  public:

    lepton_distribution_expansion_term( int const given_power_of_energy,
                                        int const given_power_of_logarithm,
                                        double const given_coefficient,
               lepton_distribution_expansion_term* const given_reference_term )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    // if not given a reference term, this sets its pointer to NULL.
    lepton_distribution_expansion_term( int const given_power_of_energy,
                                        int const given_power_of_logarithm,
                                        double const given_coefficient )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~lepton_distribution_expansion_term()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    lepton_distribution_expansion_term const*
    get_reference_term()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    int
    get_power_of_energy()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    int
    get_power_of_logarithm()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_coefficient()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    set_coefficient( double const given_value )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    normalize_coefficient( double const normalization )
    // this just divides coefficient by normalization.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    inline double
    evaluate( double const given_energy )
    const
    // this evaluates this term for the given input energy.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    int power_of_energy;
    int power_of_logarithm;
    double coefficient;
    lepton_distribution_expansion_term* const reference_term;
    // this is for ease of updating coefficients in tauon-decay distributions.

  };


  // this is a class to hold a set of lepton_distribution_expansion_term
  // objects relevant for an energy range for the distribution.
  class segment_term_set
  {

  public:

    segment_term_set()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~segment_term_set()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    double
    get_segment_start()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    set_segment_start( double const given_value )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_segment_end()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    set_segment_end( double const given_value )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    set_segment_range( double const given_start,
                       double const given_end )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    std::vector< lepton_distribution_expansion_term* > const*
    get_terms()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    int
    get_number_of_terms()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    lepton_distribution_expansion_term*
    add_term( lepton_distribution_expansion_term* const given_term )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    normalize_coefficients( double const normalization )
    // this goes through each term & gets it to divide its coefficient by
    // normalization.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    double
    evaluate( double const given_energy )
    const
    // this evaluates the sum of the set of terms for the given input energy.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    std::vector< lepton_distribution_expansion_term* > segment_terms;
    double segment_start;
    double segment_end;

  };



  /* this class is the base abstract class for the functors for obtaining
   * energy distributions for the visible decay products of tau leptons with
   * given energy distributions.
   * (yes, left-handed tau leptons / right-handed tau antileptons give harder
   * muons / electrons but softer pions than the other helicity.)
   */
  class tau_decay_coefficient
  {

  public:

    tau_decay_coefficient()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~tau_decay_coefficient()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual double
    operator()( int const visible_product_energy_power
                /* the power of the visible product's energy for the term. */,
                int const visible_product_log_power
/* the power of the logarithm of the visible product's energy for the term. */,
                int const tau_lepton_energy_power
                /* the power of the tau lepton's energy for the term. */,
                int const tau_lepton_log_power
     /* the power of the logarithm of the tau lepton's energy for the term. */,
                double const Emin
         /* the minimum of the range of the tau lepton energies considered. */,
                double const Emax
         /* the maximum of the range of the tau lepton energies considered. */,
                bool const is_inside_range )
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
  class hard_muon_from_tau : public tau_decay_coefficient
  {

  public:

    hard_muon_from_tau()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~hard_muon_from_tau()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    virtual double
    operator()( int const muon_energy_power
                /* the power of the muon's energy for the term. */,
                int const muon_log_power
           /* the power of the logarithm of the muon's energy for the term. */,
                int const tau_energy_power
                /* the power of the tau lepton's energy for the term. */,
                int const tau_log_power
     /* the power of the logarithm of the tau lepton's energy for the term. */,
                double const Emin
         /* the minimum of the range of the tau lepton energies considered. */,
                double const Emax
         /* the maximum of the range of the tau lepton energies considered. */,
                bool const is_inside_range )
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;


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
  class soft_muon_from_tau : public tau_decay_coefficient
  {

  public:

    soft_muon_from_tau()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~soft_muon_from_tau()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    virtual double
    operator()( int const muon_energy_power
                /* the power of the muon's energy for the term. */,
                int const muon_log_power
           /* the power of the logarithm of the muon's energy for the term. */,
                int const tau_energy_power
                /* the power of the tau lepton's energy for the term. */,
                int const tau_log_power
     /* the power of the logarithm of the tau lepton's energy for the term. */,
                double const Emin
         /* the minimum of the range of the tau lepton energies considered. */,
                double const Emax
         /* the maximum of the range of the tau lepton energies considered. */,
                bool const is_inside_range )
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;


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
  class hard_pion_from_tau : public tau_decay_coefficient
  {

  public:

    hard_pion_from_tau()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~hard_pion_from_tau()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    virtual double
    operator()( int const pion_energy_power
                /* the power of the pion's energy for the term. */,
                int const pion_log_power
           /* the power of the logarithm of the pion's energy for the term. */,
                int const tau_energy_power
                /* the power of the tau lepton's energy for the term. */,
                int const tau_log_power
     /* the power of the logarithm of the tau lepton's energy for the term. */,
                double const Emin
         /* the minimum of the range of the tau lepton energies considered. */,
                double const Emax
         /* the maximum of the range of the tau lepton energies considered. */,
                bool const is_inside_range )
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;


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
  class soft_pion_from_tau : public tau_decay_coefficient
  {

  public:

    soft_pion_from_tau()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~soft_pion_from_tau()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    virtual double
    operator()( int const pion_energy_power
                /* the power of the pion's energy for the term. */,
                int const pion_log_power
           /* the power of the logarithm of the pion's energy for the term. */,
                int const tau_energy_power
                /* the power of the tau lepton's energy for the term. */,
                int const tau_log_power
     /* the power of the logarithm of the tau lepton's energy for the term. */,
                double const Emin
         /* the minimum of the range of the tau lepton energies considered. */,
                double const Emax
         /* the maximum of the range of the tau lepton energies considered. */,
                bool const is_inside_range )
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    //protected:

      // nothing.

  };


  /* this holds together a pair of segment_term_set pointers for the 2 parts of
   * of the distribution for the visible decay product of a tau lepton with
   * distribution given by the 3rd segment_term_set pointer the
   * tau_segment_triple has. it performs the role of updating the visible decay
   * product terms based on the reference terms.
   */
  class tau_segment_triple
  {

  public:

    tau_segment_triple( std::vector< segment_term_set* >* given_segment_set,
                        segment_term_set* given_segment,
                        tau_decay_coefficient const* const given_tau_decay )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~tau_segment_triple()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    update_segments()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    segment_term_set* reference_segment;
    segment_term_set* low_segment;
    segment_term_set* high_segment;
    tau_decay_coefficient const* const tau_decay;

  };

  /* this class is just an abstract base class for some polymorphic derived
   * classes. However, it does do most of the work, the derived classes just
   * set up the lepton_distribution_expansion_term objects properly.
   */
  class leptonEnergyDistribution : public getsReadiedForNewPoint
  {
  public:
    leptonEnergyDistribution( readierForNewPoint* const readier,
                              CppSLHA::CppSLHA0 const* const spectrum,
                     CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                     CppSLHA::particle_property_set const* const thirdParticle,
                   CppSLHA::particle_property_set const* const fourthParticle )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~leptonEnergyDistribution()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    getMinimumEnergy()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getMaximumEnergy()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    std::vector< segment_term_set* >*
    getSegments()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    std::vector< segment_term_set* > const*
    inspectSegments()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    valueAt( double const givenEnergy )
    // this returns the value of the distribution at the given value, with
    // arbitrary normalization (different for each distribution).
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    CppSLHA::CppSLHA0 const* spectrum;
    // this holds all the information about the MSSM point. care must be taken
    // to ensure that it has the appropriate BLOCKs for the distribution.
    double minimumEnergy;
    double maximumEnergy;

    std::vector< segment_term_set* > segments;

    CppSLHA::particle_property_set const* const firstParticle;
    effectiveSquarkMassHolder* const effectiveSquarkMass;
    CppSLHA::particle_property_set const* const secondParticle;
    CppSLHA::particle_property_set const* const thirdParticle;
    CppSLHA::particle_property_set const* const fourthParticle;
    /* these masses could be pulled out of spectrum_data every time they are
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
    double normalization;
    // this is *an overall denominator*, rather than multiplicative
    // normalization.
    double normalizationInverse;
    // this is *an overall multiplicative normalization*
    // (= 1.0 / normalization).
    // this is not strictly necessary, but is handy for distributions from
    // on-shell cascades:
    double productionFrameEnergy;

    void
    reset()
    // this updates the masses then calls calculate_coefficients().
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual void
    calculateCoefficients()
    = 0;
    // this is overridden by functions which calculate & set the coefficients
    // of all the terms of the distribution.

    void
    normalizeCoefficients()
    // this goes through each segment & gets it to divide its terms'
    // coefficients by normalization.
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };


  /* this is a derived class which sets up the energy distribution for a muon
   * or pion which is from the decay of a left- or right-handed tau lepton or
   * antilepton with the provided energy distribution.
   */
  class visible_tau_decay_product : public leptonEnergyDistribution
  {

  public:

    visible_tau_decay_product( readierForNewPoint* const given_readier,
                      leptonEnergyDistribution* const given_tau_distribution,
                               tau_decay_coefficient const* const tau_decay )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~visible_tau_decay_product()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    leptonEnergyDistribution* const tau_distribution;
    std::vector< tau_segment_triple* > tau_triples;

    virtual void
    calculateCoefficients()
    // this is overridden by functions which calculate & set the coefficients
    // of all the terms of the distribution.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };




  // inline functions


  inline lepton_distribution_expansion_term const*
  lepton_distribution_expansion_term::get_reference_term()
  const
  {

    return reference_term;

  }

  inline int
  lepton_distribution_expansion_term::get_power_of_energy()
  const
  {

    return power_of_energy;

  }

  inline int
  lepton_distribution_expansion_term::get_power_of_logarithm()
  const
  {

    return power_of_logarithm;

  }

  inline double
  lepton_distribution_expansion_term::get_coefficient()
  const
  {

    return coefficient;

  }

  inline void
  lepton_distribution_expansion_term::set_coefficient(
                                                     double const given_value )
  {

     coefficient = given_value;

  }

  inline void
  lepton_distribution_expansion_term::normalize_coefficient(
                                           double const inverse_normalization )
  // this just divides coefficient by normalization.
  {

    coefficient *= inverse_normalization;

  }



  inline double
  segment_term_set::get_segment_start()
  const
  {

    return segment_start;

  }

  inline void
  segment_term_set::set_segment_start( double const given_value )
  {

    segment_start = given_value;

  }

  inline double
  segment_term_set::get_segment_end()
  const
  {

    return segment_end;

  }

  inline void
  segment_term_set::set_segment_end( double const given_value )
  {

    segment_end = given_value;

  }

  inline void
  segment_term_set::set_segment_range( double const given_start,
                                       double const given_end )
  {

    segment_start = given_start;
    segment_end = given_end;

  }

  inline std::vector< lepton_distribution_expansion_term* > const*
  segment_term_set::get_terms()
  const
  {

    return &segment_terms;

  }

  inline int
  segment_term_set::get_number_of_terms()
  const
  {

    return segment_terms.size();

  }

  inline lepton_distribution_expansion_term*
  segment_term_set::add_term(
                         lepton_distribution_expansion_term* const given_term )
  {

    segment_terms.push_back( given_term );

    return given_term;

  }

  inline void
  segment_term_set::normalize_coefficients(
                                           double const inverse_normalization )
  // this goes through each term & gets it to multiply its coefficient by
  // inverse_normalization.
  {

    for( std::vector< lepton_distribution_expansion_term* >::iterator
         term_iterator = segment_terms.begin();
         segment_terms.end() > term_iterator;
         ++term_iterator )
      {

        (*term_iterator)->normalize_coefficient( inverse_normalization );

      }

  }


  inline double
  leptonEnergyDistribution::valueAt( double const givenEnergy )
  // this returns the value of the distribution at the given value, with
  // arbitrary normalization (different for each distribution).
  {

    if( needsToPrepareForThisPoint() )
      {

        reset();
        finishPreparingForThisPoint();

      }

    double return_value = 0.0;

    for( std::vector< segment_term_set* >::const_iterator
         summation_iterator = segments.begin();
         summation_iterator != segments.end();
         ++summation_iterator )
      {

        return_value += (*summation_iterator)->evaluate( givenEnergy );

      }

    return return_value;

  }

  inline double
  leptonEnergyDistribution::getMinimumEnergy()
  {

    if( needsToPrepareForThisPoint() )
      {

        reset();
        finishPreparingForThisPoint();

      }

    return minimumEnergy;

  }

  inline double
  leptonEnergyDistribution::getMaximumEnergy()
  {

    if( needsToPrepareForThisPoint() )
      {

        reset();
        finishPreparingForThisPoint();

      }

    return maximumEnergy;

  }

  inline std::vector< segment_term_set* >*
  leptonEnergyDistribution::getSegments()
  {

    return &segments;

  }

  inline std::vector< segment_term_set* > const*
  leptonEnergyDistribution::inspectSegments()
  const
  {

    return &segments;

  }

  inline void
  leptonEnergyDistribution::normalizeCoefficients()
  // this goes through each segment & gets it to divide its terms' coefficients
  // by normalization.
  {

    normalizationInverse = ( 1.0 / normalization );

    for( std::vector< segment_term_set* >::iterator
         segment_iterator = segments.begin();
         segments.end() > segment_iterator;
         ++segment_iterator )
      {

        (*segment_iterator)->normalize_coefficients( normalizationInverse );

      }

  }

}  // end of LHC_FASER namespace.

#endif  // LHC_FASER_BASE_LEPTON_DISTRIBUTION_STUFF_HPP
