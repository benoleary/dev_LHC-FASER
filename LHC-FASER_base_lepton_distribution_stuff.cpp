/*
 * LHC-FASER_base_lepton_distribution_stuff.cpp
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


#include "LHC-FASER_base_lepton_distribution_stuff.hpp"


namespace LHC_FASER
{
  effectiveSquarkMassHolder::effectiveSquarkMassHolder()
  {
    // does nothing.
  }

  effectiveSquarkMassHolder::~effectiveSquarkMassHolder()
  {
    // does nothing.
  }



  onshellSquarkMassHolder::onshellSquarkMassHolder(
                  CppSLHA::particle_property_set const* const onshellSquark ) :
    onshellSquark( onshellSquark )
  {
    // just an initialization list.
  }

  onshellSquarkMassHolder::~onshellSquarkMassHolder()
  {
    // does nothing.
  }



  lepton_distribution_expansion_term::lepton_distribution_expansion_term(
                                               int const given_power_of_energy,
                                            int const given_power_of_logarithm,
                                             double const given_coefficient,
             lepton_distribution_expansion_term* const given_reference_term ) :
    power_of_energy( given_power_of_energy ),
    power_of_logarithm( given_power_of_logarithm ),
    coefficient( given_coefficient ),
    reference_term( given_reference_term )
  {

    // just an initialization list.

  }

  lepton_distribution_expansion_term::lepton_distribution_expansion_term(
                                               int const given_power_of_energy,
                                            int const given_power_of_logarithm,
                                             double const given_coefficient ) :
    power_of_energy( given_power_of_energy ),
    power_of_logarithm( given_power_of_logarithm ),
    coefficient( given_coefficient ),
    reference_term( NULL )
  {

    // just an initialization list.

  }

  lepton_distribution_expansion_term::~lepton_distribution_expansion_term()
  {

    // does nothing.

  }


  double
  lepton_distribution_expansion_term::evaluate( double const given_energy )
  const
  // this evaluates this term for the given input energy.
  {

    // debugging:
    /**std::cout << std::endl << "debugging: "
    << "lepton_distribution_expansion_term::evaluate( " << given_energy
    << " ) called";
    std::cout << std::endl;**/


    double return_value = 0.0;

    if( ( 0.0 < given_energy )
        // only give a value if the energy is above 0.0 since otherwise the
        // logarithms are not going to work.
        &&
        ( 0.0 != coefficient ) )
      // & don't bother if the coefficient is 0.0 anyway.
      {

        if( -1 == power_of_energy )
          {

            return_value = ( coefficient / given_energy );

          }
        else if( 0 == power_of_energy )
          {

            return_value = coefficient;

          }
        else if( 1 == power_of_energy )
          {

            return_value = ( coefficient * given_energy );

          }
        else if( 2 == power_of_energy )
          {

            return_value = ( coefficient * given_energy * given_energy );

          }
        else if( 3 == power_of_energy )
          {

            return_value = ( coefficient * given_energy
                             * given_energy * given_energy );

          }
        // I could set up something recursive, but there shouldn't be any
        // higher powers than 3.

        if( 0 < power_of_logarithm )
          {

            double log_energy = log( given_energy );

            /*for( int logarithm_counter = 0;
            logarithm_counter < power_of_logarithm;
            logarithm_counter++ )
            {

            return_value *= log_energy;

            }*/
            // go go super-irrelevant, premature optimization!
            if( 1 == power_of_logarithm )
              {

                return_value *= log_energy;

              }
            else if( 2 == power_of_logarithm )
              {

                return_value *= ( log_energy * log_energy );

              }
            else if( 3 == power_of_logarithm )
              {

                return_value *= ( log_energy * log_energy * log_energy );

              }

          }

      }

    return return_value;

  }


  segment_term_set::segment_term_set() :
    segment_start( CppSLHA::CppSLHA_global::really_wrong_value ),
    segment_end( CppSLHA::CppSLHA_global::really_wrong_value )
  {

    // just an initialization list.

  }

  segment_term_set::~segment_term_set()
  {

    // does nothing.

  }


  double
  segment_term_set::evaluate( double const given_energy )
  const
  // this evaluates the sum of the set of terms for the given input energy.
  {

    // debugging:
    /**std::cout << std::endl << "debugging: "
    << "segment_term_set::evaluate( " << given_energy
    << " ) called";
    std::cout << std::endl;**/

    double return_value = 0.0;

    // only evaluate if the input is within the segment's range:
    if( ( given_energy >= segment_start )
        &&
        ( given_energy <= segment_end ) )
      {

        // debugging:
        /**int segment_term_counter = 0;
        std::cout << std::endl << "debugging: "
        << "given_energy within range, segment_terms.size() = "
        << segment_terms.size();
        std::cout << std::endl;**/

        for( std::vector< lepton_distribution_expansion_term* >::const_iterator
             summation_iterator = segment_terms.begin();
             summation_iterator < segment_terms.end();
             ++summation_iterator )
          {

            // debugging:
            /**std::cout << std::endl << "debugging: "
            << "segment_term_counter = " << ++segment_term_counter;
            std::cout << std::endl;**/

            return_value += (*summation_iterator)->evaluate( given_energy );

          }

      }

    return return_value;

  }


  tau_decay_coefficient::tau_decay_coefficient()
  {

    // just an initialization.

  }

  tau_decay_coefficient::~tau_decay_coefficient()
  {

    // does nothing.

  }


  hard_muon_from_tau::hard_muon_from_tau() :
    tau_decay_coefficient()
  {

    // just an initialization.

  }

  hard_muon_from_tau::~hard_muon_from_tau()
  {

    // does nothing.

  }


  double
  hard_muon_from_tau::operator()( int const muon_energy_power
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
  /* this returns the coefficient for the muon distribution which would come
   * from a left-handed tau lepton distribution of the given power of the tau
   * lepton's energy & power of logarithm thereof, for requested integer powers
   * of the muon's energy or logarithm thereof, which also depends on whether
   * the muon's energy is inside the range of the tau lepton's energy for this
   * segment.
   */
  {

    // debugging:
    /**std::cout << std::endl << "debugging: "
    << "hard_muon_from_tau::operator()( " << muon_energy_power << ", "
    << muon_log_power << ", "
    << tau_energy_power << ", "
    << tau_log_power << ", "
    << Emin << ", "
    << Emax << ", "
    << is_inside_range
    << " ) called";
    std::cout << std::endl;**/

    double return_value = 0.0;

    if( ( -1 == tau_energy_power )
        &&
        ( 0 == tau_log_power ) )
      {

        if( ( -1 == muon_energy_power )
            &&
            ( 0 == muon_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = 1.0;

              }

          }
        else if( ( 0 == muon_energy_power )
                 &&
                 ( 0 == muon_log_power ) )
          {

            return_value = ( -4.0 / ( 3.0 * Emax ) );

            if( !is_inside_range )
              {

                return_value += ( 4.0 / ( 3.0 * Emin ) );

              }

          }
        else if( ( 3 == muon_energy_power )
                 &&
                 ( 0 == muon_log_power ) )
          {

            return_value = ( 1.0 / ( 3.0 * Emax * Emax * Emax * Emax ) );

            if( !is_inside_range )
              {

                return_value -= ( 1.0 / ( 3.0 * Emin * Emin * Emin * Emin ) );

              }

          }

      }  // end of if( -1 == tau_energy_power ) && ( 0 == tau_log_power )
    else if( ( 0 == tau_energy_power )
             &&
             ( 0 == tau_log_power ) )
      {

        if( ( 0 == muon_energy_power )
            &&
            ( 0 == muon_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( ( 4.0 * log_max ) / 3.0 );

            if( is_inside_range )
              {

                return_value -= ( 4.0 / 9.0 );

              }
            else
              {

                double log_min = log( Emin );

                return_value -= ( ( 4.0 * log_min ) / 3.0 );

              }

          }
        else if( ( 0 == muon_energy_power )
                 &&
                 ( 1 == muon_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = ( -4.0 / 3.0 );

              }

          }
        else if( ( 3 == muon_energy_power )
                 &&
                 ( 0 == muon_log_power ) )
          {

            return_value = ( 4.0 / ( 9.0 * Emax * Emax * Emax ) );

            if( !is_inside_range )
              {

                return_value -= ( 4.0 / ( 9.0 * Emin * Emin * Emin ) );

              }

          }

      }  // end of if( 0 == tau_energy_power ) && ( 0 == tau_log_power )
    else if( ( 0 == tau_energy_power )
             &&
             ( 1 == tau_log_power ) )
      {

        if( ( 0 == muon_energy_power )
            &&
            ( 0 == muon_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( ( 2.0 * log_max * log_max ) / 3.0 );

            if( is_inside_range )
              {

                return_value -= ( 4.0 / 27.0 );

              }
            else
              {

                double log_min = log( Emin );

                return_value -= ( ( 2.0 * log_min * log_min ) / 3.0 );

              }

          }

        else if( ( 0 == muon_energy_power )
                 &&
                 ( 1 == muon_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = ( -4.0 / 9.0 );

              }

          }
        else if( ( 0 == muon_energy_power )
                 &&
                 ( 2 == muon_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = ( -2.0 / 3.0 );

              }

          }
        else if( ( 3 == muon_energy_power )
                 &&
                 ( 0 == muon_log_power ) )
          {

            double log_max = log( Emax );

            return_value
            = ( ( 4.0 * ( 1.0 + 3.0 * log_max ) )
                / ( 27.0 * Emax * Emax * Emax ) );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value -= ( ( 4.0 * ( 1.0 + 3.0 * log_min ) )
                                  / ( 27.0 * Emin * Emin * Emin ) );

              }

          }

      }  // end of if( 0 == tau_energy_power ) && ( 1 == tau_log_power )
    else if( ( 1 == tau_energy_power )
             &&
             ( 0 == tau_log_power ) )
      {

        if( ( 0 == muon_energy_power )
            &&
            ( 0 == muon_log_power ) )
          {

            return_value = ( ( 4.0 * Emax ) / 3.0 );

            if( !is_inside_range )
              {

                return_value -= ( ( 4.0 * Emin ) / 3.0 );

              }

          }
        else if( ( 1 == muon_energy_power )
                 &&
                 ( 0 == muon_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = -2.0;

              }

          }
        else if( ( 3 == muon_energy_power )
                 &&
                 ( 0 == muon_log_power ) )
          {

            return_value = ( 2.0 / ( 3.0 * Emax * Emax ) );

            if( !is_inside_range )
              {

                return_value -= ( 2.0 / ( 3.0 * Emin * Emin ) );

              }

          }

      }  // end of if( 1 == tau_energy_power ) && ( 0 == tau_log_power )
    else if( ( 1 == tau_energy_power )
             &&
             ( 1 == tau_log_power ) )
      {

        if( ( 0 == muon_energy_power )
            &&
            ( 0 == muon_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( ( 4.0 * Emax * ( log_max - 1.0 ) ) / 3.0 );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value += ( ( 4.0 * Emin * ( 1.0 - log_min ) ) / 3.0 );

              }

          }
        else if( ( 1 == muon_energy_power )
                 &&
                 ( 0 == muon_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = 1.0;

              }

          }
        else if( ( 1 == muon_energy_power )
                 &&
                 ( 1 == muon_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = -2.0;

              }

          }
        else if( ( 3 == muon_energy_power )
                 &&
                 ( 0 == muon_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( ( 1.0 + 2.0 * log_max ) / ( 3.0 * Emax * Emax ) );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value
                -= ( ( 1.0 + 2.0 * log_min ) / ( 3.0 * Emin * Emin ) );

              }

          }

      }  // end of if( 1 == tau_energy_power ) && ( 1 == tau_log_power )
    else if( ( 2 == tau_energy_power )
             &&
             ( 0 == tau_log_power ) )
      {

        if( ( 0 == muon_energy_power )
            &&
            ( 0 == muon_log_power ) )
          {

            return_value = ( ( 2.0 * Emax * Emax ) / 3.0 );

            if( !is_inside_range )
              {

                return_value -= ( ( 2.0 * Emin * Emin ) / 3.0 );

              }

          }
        else if( ( 2 == muon_energy_power )
                 &&
                 ( 0 == muon_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = -2.0;

              }

          }
        else if( ( 3 == muon_energy_power )
                 &&
                 ( 0 == muon_log_power ) )
          {

            return_value = ( 4.0 / ( 3.0 * Emax ) );

            if( !is_inside_range )
              {

                return_value -= ( 4.0 / ( 3.0 * Emin ) );

              }

          }

      }  // end of if( 2 == tau_energy_power ) && ( 0 == tau_log_power )

    return return_value;

  }


  soft_muon_from_tau::soft_muon_from_tau() :
    tau_decay_coefficient()
  {

    // just an initialization.

  }

  soft_muon_from_tau::~soft_muon_from_tau()
  {

    // does nothing.

  }


  double
  soft_muon_from_tau::operator()( int const muon_energy_power
                            /* the power of the muon's energy for the term. */,
                                  int const muon_log_power
           /* the power of the logarithm of the muon's energy for the term. */,
                                  int const tau_energy_power
                       /* the power of the tau lepton's energyfor the term. */,
                                  int const tau_log_power
     /* the power of the logarithm of the tau lepton's energy for the term. */,
                                  double const Emin
         /* the minimum of the range of the tau lepton energies considered. */,
                                  double const Emax
         /* the maximum of the range of the tau lepton energies considered. */,
                                  bool const is_inside_range )
  const
  /* this returns the coefficient for the muon distribution which would come
   * from a right-handed tau lepton distribution of the given power of the tau
   * lepton's energy & power of logarithm thereof, for requested integer powers
   * of the muon's energy or logarithm thereof, which also depends on whether
   * the muon's energy is inside the range of the tau lepton's energy for this
   * segment.
   */
  {

    double return_value = 0.0;

    if( ( -1 == tau_energy_power )
        &&
        ( 0 == tau_log_power ) )
      {

        if( ( -1 == muon_energy_power )
            &&
            ( 0 == muon_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = 1.0;

              }

          }
        else if( ( 0 == muon_energy_power )
                 &&
                 ( 0 == muon_log_power ) )
          {

            return_value = ( -2.0 / Emax );

            if( !is_inside_range )
              {

                return_value += ( 2.0 / Emin );

              }

          }
        else if( ( 2 == muon_energy_power )
                 &&
                 ( 0 == muon_log_power ) )
          {

            return_value = ( 2.0 / ( Emax * Emax * Emax ) );

            if( !is_inside_range )
              {

                return_value -= ( 2.0 / ( Emin * Emin * Emin ) );

              }

          }
        else if( ( 3 == muon_energy_power )
                 &&
                 ( 0 == muon_log_power ) )
          {

            return_value = ( -1.0 / ( Emax * Emax * Emax * Emax ) );

            if( !is_inside_range )
              {

                return_value += ( 1.0 / ( Emin * Emin * Emin * Emin ) );

              }

          }

      }  // end of if( -1 == tau_energy_power ) && ( 0 == tau_log_power )
    else if( ( 0 == tau_energy_power )
             &&
             ( 0 == tau_log_power ) )
      {

        if( ( 0 == muon_energy_power )
            &&
            ( 0 == muon_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( 2.0 * log_max );

            if( is_inside_range )
              {

                return_value -= ( 5.0 / 3.0 );

              }
            else
              {

                double log_min = log( Emin );

                return_value -= ( 2.0 * log_min );

              }

          }
        else if( ( 0 == muon_energy_power )
                 &&
                 ( 1 == muon_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = -2.0;

              }

          }
        else if( ( 2 == muon_energy_power )
                 &&
                 ( 0 == muon_log_power ) )
          {

            return_value = ( 3.0 / ( Emax * Emax ) );

            if( !is_inside_range )
              {

                return_value -= ( 3.0 / ( Emin * Emin ) );

              }

          }
        else if( ( 3 == muon_energy_power )
                 &&
                 ( 0 == muon_log_power ) )
          {

            return_value = ( -4.0 / ( 3.0 * Emax * Emax * Emax ) );

            if( !is_inside_range )
              {

                return_value += ( 4.0 / ( 3.0 * Emin * Emin * Emin ) );

              }

          }

      }  // end of if( 0 == tau_energy_power ) && ( 0 == tau_log_power )
    else if( ( 0 == tau_energy_power )
             &&
             ( 1 == tau_log_power ) )
      {

        if( ( 0 == muon_energy_power )
            &&
            ( 0 == muon_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( log_max * log_max );

            if( is_inside_range )
              {

                return_value -= ( 19.0 / 18.0 );

              }
            else
              {

                double log_min = log( Emin );

                return_value -= ( log_min * log_min );

              }

          }

        else if( ( 0 == muon_energy_power )
                 &&
                 ( 1 == muon_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = ( -5.0 / 3.0 );

              }

          }
        else if( ( 0 == muon_energy_power )
                 &&
                 ( 2 == muon_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = -1.0;

              }

          }
        else if( ( 2 == muon_energy_power )
                 &&
                 ( 0 == muon_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( ( 1.5 + 3.0 * log_max ) / ( Emax * Emax ) );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value -= ( ( 1.5 + 3.0 * log_min ) / ( Emin * Emin ) );

              }

          }
        else if( ( 3 == muon_energy_power )
                 &&
                 ( 0 == muon_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( ( -4.0 * ( 1.0 + 3.0 * log_max ) )
                             / ( 9.0 * Emax * Emax * Emax ) );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value += ( ( 4.0 * ( 1.0 + 3.0 * log_min ) )
                                  / ( 9.0 * Emin * Emin * Emin ) );

              }

          }

      }  // end of if( 0 == tau_energy_power ) && ( 1 == tau_log_power )
    else if( ( 1 == tau_energy_power )
             &&
             ( 0 == tau_log_power ) )
      {

        if( ( 0 == muon_energy_power )
            &&
            ( 0 == muon_log_power ) )
          {

            return_value = ( 2.0 * Emax );

            if( !is_inside_range )
              {

                return_value -= ( 2.0 * Emin );

              }

          }
        else if( ( 1 == muon_energy_power )
                 &&
                 ( 0 == muon_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = -6.0;

              }

          }
        else if( ( 2 == muon_energy_power )
                 &&
                 ( 0 == muon_log_power ) )
          {

            return_value = ( 6.0 / Emax );

            if( !is_inside_range )
              {

                return_value -= ( 6.0 / Emin );

              }

          }
        else if( ( 3 == muon_energy_power )
                 &&
                 ( 0 == muon_log_power ) )
          {

            return_value = ( -2.0 / ( Emax * Emax ) );

            if( !is_inside_range )
              {

                return_value += ( 2.0 / ( Emin * Emin ) );

              }

          }

      }  // end of if( 1 == tau_energy_power ) && ( 0 == tau_log_power )
    else if( ( 1 == tau_energy_power )
             &&
             ( 1 == tau_log_power ) )
      {

        if( ( 0 == muon_energy_power )
            &&
            ( 0 == muon_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( 2.0 * Emax * ( log_max - 1.0 ) );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value += ( 2.0 * Emin * ( 1.0 - log_min ) );

              }

          }
        else if( ( 1 == muon_energy_power )
                 &&
                 ( 0 == muon_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = -3.0;

              }

          }
        else if( ( 1 == muon_energy_power )
                 &&
                 ( 1 == muon_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = -6.0;

              }

          }
        else if( ( 2 == muon_energy_power )
                 &&
                 ( 0 == muon_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( ( 6.0 * ( 1.0 + log_max ) ) / Emax );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value -= ( ( 6.0 * ( 1.0 + log_min ) ) / Emin );

              }

          }
        else if( ( 3 == muon_energy_power )
                 &&
                 ( 0 == muon_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( ( -1.0 * ( 1.0 + 2.0 * log_max ) )
                             / ( Emax * Emax ) );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value += ( ( 1.0 + 2.0 * log_min ) / ( Emin * Emin ) );

              }

          }

      }  // end of if( 1 == tau_energy_power ) && ( 1 == tau_log_power )
    else if( ( 2 == tau_energy_power )
             &&
             ( 0 == tau_log_power ) )
      {

        if( ( 0 == muon_energy_power )
            &&
            ( 0 == muon_log_power ) )
          {

            return_value = ( Emax * Emax );

            if( !is_inside_range )
              {

                return_value -= ( Emin * Emin );

              }

          }
        else if( ( 2 == muon_energy_power )
                 &&
                 ( 0 == muon_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( -6.0 * log_max );

            if( is_inside_range )
              {

                return_value += 3.0;

              }
            else
              {

                double log_min = log( Emin );

                return_value += ( 6.0 * log_min );


              }

          }
        else if( ( 2 == muon_energy_power )
                 &&
                 ( 1 == muon_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = 6.0;

              }

          }
        else if( ( 3 == muon_energy_power )
                 &&
                 ( 0 == muon_log_power ) )
          {

            return_value = ( -4.0 / ( Emax ) );

            if( !is_inside_range )
              {

                return_value += ( 4.0 / ( Emin ) );

              }

          }

      }  // end of if( 2 == tau_energy_power ) && ( 0 == tau_log_power )

    return return_value;

  }


  hard_pion_from_tau::hard_pion_from_tau() :
    tau_decay_coefficient()
  {

    // just an initialization.

  }

  hard_pion_from_tau::~hard_pion_from_tau()
  {

    // does nothing.

  }


  double
  hard_pion_from_tau::operator()( int const pion_energy_power
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
  /* this returns the coefficient for the pion distribution which would come
   * from a right-handed tau lepton distribution of the given power of the
   * tau lepton's energy & power of logarithm thereof, for requested integer
   * powers of the pion's energy or logarithm thereof, which also depends on
   * whether the pion's energy is inside the range of the tau lepton's energy
   * for this segment.
   */
  {

    double return_value = 0.0;

    if( ( -1 == tau_energy_power )
        &&
        ( 0 == tau_log_power ) )
      {

        if( ( -1 == pion_energy_power )
            &&
            ( 0 == pion_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = 1.0;

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 0 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value
            = ( ( -1.0 * ( 1.0 + 2.0 * log_max ) ) / ( Emax * Emax ) );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value
                += ( ( ( 1.0 + 2.0 * log_min ) ) / ( Emin * Emin ) );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 1 == pion_log_power ) )
          {

            return_value = ( 2.0 / ( Emax * Emax ) );

            if( !is_inside_range )
              {

                return_value -= ( 2.0 / ( Emin * Emin ) );

              }

          }

      }  // end of if( -1 == tau_energy_power ) && ( 0 == tau_log_power )
    else if( ( 0 == tau_energy_power )
             &&
             ( 0 == tau_log_power ) )
      {

        if( ( 0 == pion_energy_power )
            &&
            ( 0 == pion_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = 4.0;

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 0 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( ( -4.0 * ( 1.0 + log_max ) ) / Emax );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value += ( ( 4.0 * ( 1.0 + log_min ) ) / Emin );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 1 == pion_log_power ) )
          {

            return_value = ( 4.0 / Emax );

            if( !is_inside_range )
              {

                return_value -= ( 4.0 / Emin );

              }

          }

      }  // end of if( 0 == tau_energy_power ) && ( 0 == tau_log_power )
    else if( ( 0 == tau_energy_power )
             &&
             ( 1 == tau_log_power ) )
      {

        if( ( 0 == pion_energy_power )
            &&
            ( 0 == pion_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = 8.0;

              }

          }
        else if( ( 0 == pion_energy_power )
                 &&
                 ( 1 == pion_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = 4.0;

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 0 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value
            = ( ( -4.0 * ( 2.0 + 2.0 * log_max + log_max * log_max ) )
                / Emax );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value
                += ( ( 4.0 * ( 2.0 + 2.0 * log_min + log_min * log_min ) )
                     / Emin );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 1 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( ( 4.0 * ( 1.0 + log_max ) ) / Emax );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value -= ( ( 4.0 * ( 1.0 + log_min ) ) / Emin );

              }

          }

      }  // end of if( 0 == tau_energy_power ) && ( 1 == tau_log_power )
    else if( ( 1 == tau_energy_power )
             &&
             ( 0 == tau_log_power ) )
      {

        if( ( 1 == pion_energy_power )
            &&
            ( 0 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( 2.0 * log_max * log_max );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value -= ( 2.0 * log_min * log_min );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 1 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( -4.0 * log_max );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value += ( 4.0 * log_min );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 2 == pion_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = 2.0;

              }

          }

      }  // end of if( 1 == tau_energy_power ) && ( 0 == tau_log_power )
    else if( ( 1 == tau_energy_power )
             &&
             ( 1 == tau_log_power ) )
      {

        if( ( 1 == pion_energy_power )
            &&
            ( 0 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( ( 4.0 * log_max * log_max * log_max ) / 3.0 );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value
                -= ( ( 4.0 * log_min * log_min * log_min ) / 3.0 );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 1 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( -2.0 * log_max * log_max );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value += ( 2.0 * log_min * log_min );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 3 == pion_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = ( 2.0 / 3.0 );

              }

          }

      }  // end of if( 1 == tau_energy_power ) && ( 1 == tau_log_power )
    else if( ( 2 == tau_energy_power )
             &&
             ( 0 == tau_log_power ) )
      {

        if( ( 1 == pion_energy_power )
            &&
            ( 0 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( 4.0 * Emax * ( log_max - 1.0 ) );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value += ( 4.0 * Emin * ( 1.0 - log_min ) );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 1 == pion_log_power ) )
          {

            return_value = ( -4.0 * Emax );

            if( !is_inside_range )
              {

                return_value += ( 4.0 * Emin );

              }

          }
        else if( ( 2 == pion_energy_power )
                 &&
                 ( 0 == pion_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = 4.0;

              }

          }

      }  // end of if( 2 == tau_energy_power ) && ( 0 == tau_log_power )

    return return_value;

  }


  soft_pion_from_tau::soft_pion_from_tau() :
    tau_decay_coefficient()
  {

    // just an initialization.

  }

  soft_pion_from_tau::~soft_pion_from_tau()
  {

    // does nothing.

  }


  double
  soft_pion_from_tau::operator()( int const pion_energy_power
                            /* the power of the pion's energy for the term. */,
                                   int const pion_log_power
           /* the power of the logarithm of the pion's energy for the term. */,
                                   int const tau_energy_power
                       /* the power of the tau lepton's energyfor the term. */,
                                   int const tau_log_power
     /* the power of the logarithm of the tau lepton's energy for the term. */,
                                   double const Emin
         /* the minimum of the range of the tau lepton energies considered. */,
                                   double const Emax
         /* the maximum of the range of the tau lepton energies considered. */,
                                   bool const is_inside_range )
  const
  /* this returns the coefficient for the pion distribution which would come
   * from a left-handed tau lepton distribution of the given power of the tau
   * lepton's energy & power of logarithm thereof, for requested integer
   * powers of the pion's energy or logarithm thereof, which also depends on
   * whether the pion's energy is inside the range of the tau lepton's energy
   * for this segment.
   */
  {

    double return_value = 0.0;

    if( ( -1 == tau_energy_power )
        &&
        ( 0 == tau_log_power ) )
      {

        if( ( -1 == pion_energy_power )
            &&
            ( 0 == pion_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = 1.0;

              }

          }
        else if( ( 0 == pion_energy_power )
                 &&
                 ( 0 == pion_log_power ) )
          {

            return_value = ( -4.0 / Emax );

            if( !is_inside_range )
              {

                return_value += ( 4.0 / Emin );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 0 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( ( 2.0 * log_max + 3.0 ) / ( Emax * Emax ) );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value -= ( ( 2.0 * log_min + 3.0 ) / ( Emin * Emin ) );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 1 == pion_log_power ) )
          {

            return_value = ( -2.0 / ( Emax * Emax ) );

            if( !is_inside_range )
              {

                return_value += ( 2.0 / ( Emin * Emin ) );

              }

          }

      }  // end of if( -1 == tau_energy_power ) && ( 0 == tau_log_power )
    else if( ( 0 == tau_energy_power )
             &&
             ( 0 == tau_log_power ) )
      {

        if( ( 0 == pion_energy_power )
            &&
            ( 0 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( 4.0 * log_max );

            if( is_inside_range )
              {

                return_value -= 8.0;

              }
            else
              {
                double log_min = log( Emin );

                return_value -= ( 4.0 * log_min );

              }

          }
        else if( ( 0 == pion_energy_power )
                 &&
                 ( 1 == pion_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = -4.0;

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 0 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( ( 4.0 * ( 2.0 + log_max ) ) / Emax );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value -= ( ( 4.0 * ( 2.0 + log_min ) ) / Emin );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 1 == pion_log_power ) )
          {

            return_value = ( -4.0 / Emax );

            if( !is_inside_range )
              {

                return_value += ( 4.0 / Emin );

              }

          }

      }  // end of if( 0 == tau_energy_power ) && ( 0 == tau_log_power )
    else if( ( 0 == tau_energy_power )
             &&
             ( 1 == tau_log_power ) )
      {

        if( ( 0 == pion_energy_power )
            &&
            ( 0 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( 2.0 * log_max * log_max );

            if( is_inside_range )
              {

                return_value -= 12.0;

              }
            else
              {

                double log_min = log( Emin );

                return_value -= ( 2.0 * log_min * log_min );

              }

          }
        else if( ( 0 == pion_energy_power )
                 &&
                 ( 1 == pion_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = -8.0;

              }

          }
        else if( ( 0 == pion_energy_power )
                 &&
                 ( 2 == pion_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = -2.0;

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 0 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value
            = ( ( 4.0 * ( 3.0 + 3.0 * log_max + log_max * log_max ) ) / Emax );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value
                -= ( ( 4.0 * ( 3.0 + 3.0 * log_min + log_min * log_min ) )
                     / Emin );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 1 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( ( -4.0 * ( 1.0 + log_max ) ) / Emax );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value += ( ( 4.0 * ( 1.0 + log_min ) ) / Emin );

              }

          }

      }  // end of if( 0 == tau_energy_power ) && ( 1 == tau_log_power )
    else if( ( 1 == tau_energy_power )
             &&
             ( 0 == tau_log_power ) )
      {

        if( ( 0 == pion_energy_power )
            &&
            ( 0 == pion_log_power ) )
          {

            return_value = ( 4.0 * Emax );

            if( !is_inside_range )
              {

                return_value -= ( 4.0 * Emin );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 0 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( -2.0 * log_max * ( log_max + 2.0 ) );

            if( is_inside_range )
              {

                return_value -= 4.0;

              }
            else
              {

                double log_min = log( Emin );

                return_value += ( 2.0 * log_min * ( log_min + 2.0 ) );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 1 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( 4.0 * log_max );

            if( is_inside_range )
              {

                return_value += 4.0;

              }
            else
              {

                double log_min = log( Emin );

                return_value -= ( 4.0 * log_min );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 2 == pion_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = -2.0;

              }

          }

      }  // end of if( 1 == tau_energy_power ) && ( 0 == tau_log_power )
    else if( ( 1 == tau_energy_power )
             &&
             ( 1 == tau_log_power ) )
      {

        if( ( 0 == pion_energy_power )
            &&
            ( 0 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( 4.0 * Emax * ( log_max - 1.0 ) );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value += ( 4.0 * Emin * ( 1.0 - log_min ) );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 0 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value
            = ( ( -2.0 * log_max * log_max * ( 3.0 + 2.0 * log_max ) ) / 3.0 );

            if( is_inside_range )
              {

                return_value += 4.0;

              }
            else
              {

                double log_min = log( Emin );

                return_value
                += ( ( 2.0 * log_min * log_min * ( 3.0 + 2.0 * log_min ) )
                     / 3.0 );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 1 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( 2.0 * log_max * log_max );

            if( is_inside_range )
              {

                return_value -= 4.0;

              }
            else
              {

                double log_min = log( Emin );

                return_value -= ( 2.0 * log_min * log_min );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 2 == pion_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = 2.0;

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 3 == pion_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = ( -2.0 / 3.0 );

              }

          }

      }  // end of if( 1 == tau_energy_power ) && ( 1 == tau_log_power )
    else if( ( 2 == tau_energy_power )
             &&
             ( 0 == tau_log_power ) )
      {

        if( ( 0 == pion_energy_power )
            &&
            ( 0 == pion_log_power ) )
          {

            return_value = ( 2.0 * Emax * Emax );

            if( !is_inside_range )
              {

                return_value -= ( 2.0 * Emin * Emin );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 0 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( -4.0 * Emax * log_max );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value += ( 4.0 * Emin * log_min );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 1 == pion_log_power ) )
          {

            return_value = ( 4.0 * Emax );

            if( !is_inside_range )
              {

                return_value -= ( 4.0 * Emin );

              }

          }
        else if( ( 2 == pion_energy_power )
                 &&
                 ( 0 == pion_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = -2.0;

              }

          }

      }  // end of if( 2 == tau_energy_power ) && ( 0 == tau_log_power )

    return return_value;

  }



  tau_segment_triple::tau_segment_triple(
                           std::vector< segment_term_set* >* given_segment_set,
                                          segment_term_set* given_segment,
                         tau_decay_coefficient const* const given_tau_decay ) :
    reference_segment( given_segment ),
    tau_decay( given_tau_decay )
  {

    double coefficient_checker;

    low_segment = new segment_term_set();
    low_segment->set_segment_start( 0.0 );

    high_segment = new segment_term_set();

    for( std::vector< lepton_distribution_expansion_term* >::const_iterator
         term_iterator = given_segment->get_terms()->begin();
         given_segment->get_terms()->end() > term_iterator;
         ++term_iterator )
      {

        for( int energy_power_counter = -1;
             3 >= energy_power_counter;
             ++energy_power_counter )
          {

            for( int log_power_counter = 0;
                 3 >= log_power_counter;
                 ++log_power_counter )
              {

                /* each term is checked to see if it's non-zero. "3.1" & "4.3"
                 * were chosen since there should be no accidental
                 * cancellations with these values, since there are just
                 * logarithms, powers of the end energies, & some rational
                 * numbers, but none are multiples of 3.1 or 4.3 so there can
                 * be no cancellation.
                 */

                coefficient_checker
                = (*tau_decay)( energy_power_counter,
                                log_power_counter,
                                (*term_iterator)->get_power_of_energy(),
                                (*term_iterator)->get_power_of_logarithm(),
                                3.1,
                                4.3,
                                false );

                if( 0.0 != coefficient_checker )
                  {

                    low_segment->add_term(
                                        new lepton_distribution_expansion_term(
                                                          energy_power_counter,
                                                             log_power_counter,
                                   CppSLHA::CppSLHA_global::really_wrong_value,
                                                            *term_iterator ) );

                  }

                coefficient_checker
                = (*tau_decay)( energy_power_counter,
                                log_power_counter,
                                (*term_iterator)->get_power_of_energy(),
                                (*term_iterator)->get_power_of_logarithm(),
                                3.1,
                                4.3,
                                true );

                if( 0.0 != coefficient_checker )
                  {

                    high_segment->add_term(
                                        new lepton_distribution_expansion_term(
                                                          energy_power_counter,
                                                             log_power_counter,
                                   CppSLHA::CppSLHA_global::really_wrong_value,
                                                            *term_iterator ) );

                  }

              }  // end of loop over powers of log( energy ).

          }  // end of loop over powers of energy.

      }  // end of loop over lepton_distribution_expansion_term pointers.

    // now that all the low_segment & high_segment terms have been added to the
    // segments, we add them to the vector of segment_term_sets.
    given_segment_set->push_back( low_segment );
    given_segment_set->push_back( high_segment );

  }

  tau_segment_triple::~tau_segment_triple()
  {

    // does nothing.

  }


  void
  tau_segment_triple::update_segments()
  {

    low_segment->set_segment_end( reference_segment->get_segment_start() );
    high_segment->set_segment_range( reference_segment->get_segment_start(),
                                     reference_segment->get_segment_end() );

    for( std::vector< lepton_distribution_expansion_term* >::const_iterator
         term_iterator = low_segment->get_terms()->begin();
         low_segment->get_terms()->end() > term_iterator;
         ++term_iterator )
      {

        (*term_iterator)->set_coefficient(
                         (*tau_decay)( (*term_iterator)->get_power_of_energy(),
                                    (*term_iterator)->get_power_of_logarithm(),
                                       (*term_iterator)->get_reference_term(
                                                      )->get_power_of_energy(),
                                       (*term_iterator)->get_reference_term(
                                                   )->get_power_of_logarithm(),
                                       0.0,
                                       low_segment->get_segment_end(),
                                       false ) );

      }
    for( std::vector< lepton_distribution_expansion_term* >::const_iterator
         term_iterator = high_segment->get_terms()->begin();
         high_segment->get_terms()->end() > term_iterator;
         ++term_iterator )
      {

        (*term_iterator)->set_coefficient(
                         (*tau_decay)( (*term_iterator)->get_power_of_energy(),
                                    (*term_iterator)->get_power_of_logarithm(),
                                       (*term_iterator)->get_reference_term(
                                                      )->get_power_of_energy(),
                                       (*term_iterator)->get_reference_term(
                                                   )->get_power_of_logarithm(),
                                       high_segment->get_segment_start(),
                                       high_segment->get_segment_end(),
                                       true ) );

      }

  }



  leptonEnergyDistribution::leptonEnergyDistribution(
                                          readierForNewPoint* const readier,
                                       CppSLHA::CppSLHA0 const* const spectrum,
                     CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                     CppSLHA::particle_property_set const* const thirdParticle,
                 CppSLHA::particle_property_set const* const fourthParticle ) :
    getsReadiedForNewPoint( readier ),
    spectrum( spectrum ),
    firstParticle( firstParticle ),
    effectiveSquarkMass( effectiveSquarkMass ),
    secondParticle( secondParticle ),
    thirdParticle( thirdParticle ),
    fourthParticle( fourthParticle )
  {
    // just an initialization list.
  }

  leptonEnergyDistribution::~leptonEnergyDistribution()
  {
    for( std::vector< segment_term_set* >::iterator
         deletionIterator = segments.begin();
         segments.end() > deletionIterator;
         ++deletionIterator )
      {
        delete *deletionIterator;
      }
  }


  void
  leptonEnergyDistribution::reset()
  {
    if( NULL != effectiveSquarkMass )
      {

        firstMass = effectiveSquarkMass->getEffectiveSquarkMass();

      }
    else
      {

        firstMass = CppSLHA::CppSLHA_global::really_wrong_value;

      }
    if( NULL != secondParticle )
      {

        secondMass = secondParticle->get_absolute_mass();

      }
    else
      {

        secondMass = CppSLHA::CppSLHA_global::really_wrong_value;

      }
    if( NULL != thirdParticle )
      {

        thirdMass = thirdParticle->get_absolute_mass();

      }
    else
      {

        thirdMass = CppSLHA::CppSLHA_global::really_wrong_value;

      }
    if( NULL != fourthParticle )
      {

        fourthMass = fourthParticle->get_absolute_mass();

      }
    else
      {

        fourthMass = CppSLHA::CppSLHA_global::really_wrong_value;

      }
    calculateCoefficients();
  }



  visible_tau_decay_product::visible_tau_decay_product(
                                    readierForNewPoint* const given_readier,
                      leptonEnergyDistribution* const given_tau_distribution,
                               tau_decay_coefficient const* const tau_decay ) :
    leptonEnergyDistribution( given_readier,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                NULL ),
    tau_distribution( given_tau_distribution )
  {

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: "
    << "visible_tau_decay_product( "
    << tau_distribution << ", " << tau_decay << " )";
    std::cout
    << std::endl;**/

    minimumEnergy = 0.0;

    for( std::vector< segment_term_set* >::iterator
         segment_iterator = given_tau_distribution->getSegments()->begin();
         given_tau_distribution->getSegments()->end() > segment_iterator;
         ++segment_iterator )
      {

        tau_triples.push_back( new tau_segment_triple( &segments,
                                                       *segment_iterator,
                                                       tau_decay ) );

      }  // end of loop through segments.

  }

  visible_tau_decay_product::~visible_tau_decay_product()
  {

    for( std::vector< tau_segment_triple* >::iterator
         deletion_iterator = tau_triples.begin();
         tau_triples.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

  }


  void
  visible_tau_decay_product::calculateCoefficients()
  {

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: "
    << "visible_tau_decay_product::calculate_coefficients()";
    std::cout
    << std::endl;**/

    maximumEnergy = tau_distribution->getMaximumEnergy();
    // now all the underlying stuff should have been updated for the point.

    for( std::vector< tau_segment_triple* >::iterator
         triple_iterator = tau_triples.begin();
         tau_triples.end() > triple_iterator;
         ++triple_iterator )
      {

        (*triple_iterator)->update_segments();

      }  // end of loop through segments.

  }

}  // end of LHC_FASER namespace

