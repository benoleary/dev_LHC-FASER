/*
 * CppSLHA_particle_data.hpp
 *
 *  Created on: 12 Sep 2009
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2009 Ben O'Leary
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

#ifndef CPPSLHA_PARTICLE_DATA_HPP_
#define CPPSLHA_PARTICLE_DATA_HPP_

// the codes assigned by the PDG to the various particles can be found in
// the CppSLHA_particle_codes.hpp file.

#include "CppSLHA_global_stuff.hpp"

namespace CppSLHA
{

  /* this is a class to hold a variable-length set of PDG codes that represent
   * the possible decay products of a particle. instances of this class do not
   * know which particle decays into these sets, rather they belong to a
   * particle_property_set instance.
   */
  class particle_decay
  {

  public:

    // constructor for a 2-body decay:
    particle_decay( int const product_PDG_code_one,
                    int const product_PDG_code_two,
                    double const input_branching_ratio )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    // constructor for a 3-body decay:
    particle_decay( int const product_PDG_code_one,
                    int const product_PDG_code_two,
                    int const product_PDG_code_three,
                    double const input_branching_ratio )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    // further constructors (4-body, 5-body, & so forth) are a bit of a waste
    // of time.

    particle_decay( particle_decay const* particle_decay_to_copy )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~particle_decay()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    std::list< int > const*
    inspect_decay_product_PDG_codes()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    std::list< int >*
    get_decay_product_PDG_codes()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    double
    get_branching_ratio()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_decay*
    set_branching_ratio( double const input_branching_ratio )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    bool
    has_as_subset( std::list< int > const* sought_subset_PDG_code_list )
    const
    // compares the decay products against a variable amount of given
    // PDG codes.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    bool
    has_this_PDG_code( int const input_PDG_code )
    const
    // checks to see if the sought particle is among the decay products.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    bool
    has_as_subset( std::list< int > const* sought_subset_PDG_code_list,
                   std::list< int > const* excluded_PDG_code_list )
    const
    // compares the decay products against a variable amount of given PDG
    // codes, along with a list of excluded PDG codes.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    bool
    is_same_as( int const product_PDG_code_one,
                int const product_PDG_code_two )
    const
    // compares the decay products against a pair of given PDG codes.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    bool
    is_same_as( int const product_PDG_code_one,
                int const product_PDG_code_two,
                int const product_PDG_code_three )
    const
    // compares the decay products against a set of 3 given PDG codes.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    bool
    is_same_as( std::list< int > const* input_PDG_code_list )
    const
    // compares the decay products against a variable amount of given PDG
    // codes.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    std::list< int > decay_product_PDG_codes;
    // it is most convenient to have a std::list to hold the decays, as shown
    // in EW_scale_spectrum::find_all_cascade_decays().

    double branching_ratio;

  };

  /* an instance of this class holds a list of particle_decay
   * instances. it handles the addition of decays to the
   * particle_property_set instance, & might gain more functions if I
   * think of any which would be useful.
   */
  class particle_decay_set_handler
  {

  public:

    particle_decay_set_handler()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~particle_decay_set_handler()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    particle_decay*
    add_decay( int const product_PDG_code_one,
               int const product_PDG_code_two,
               double const input_branching_ratio )
    /* this creates a new particle_decay with the input, then adds a
     * pointer to it to the vector of pointers to the particle_decay
     * instances of this particle_property_set.
     * N.B.: this does not check that the decay does not already exist for
     * this particle!
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_decay*
    add_decay( int const product_PDG_code_one,
               int const product_PDG_code_two,
               int const product_PDG_code_three,
               double const input_branching_ratio )
    /* this creates a new particle_decay with the input, then adds a
     * pointer to it to the vector of pointers to the particle_decay
     * instances of this particle_property_set.
     * N.B.: this does not check that the decay does not already exist for
     * this particle!
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_decay*
    add_decay( particle_decay const* const decay_to_be_copied )
    /* this creates a new particle_decay with the input, then adds a
     * pointer to it to the vector of pointers to the particle_decay
     * instances of this particle_property_set.
     * N.B.: this does not check that the decay does not already exist for
     * this particle!
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    // add_decay for 4- or 5-body decays or even more complicated decays is a
    // bit of a waste of time.


    std::list< particle_decay* > const*
    inspect_decay_set()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    std::list< particle_decay* >*
    get_decay_set()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    particle_decay_set_handler*
    become_copy_of(
                particle_decay_set_handler const* const decay_handler_to_copy )
    // this copies a particle_decay_set_handler instance.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    std::list< particle_decay* >*
    copy_decay_set(
           std::list< particle_decay* > const* const decay_sets_to_copy )
    // this copies a particle_decay_set_handler instance.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    double
    get_branching_ratio_for_exact_match( int const product_PDG_code_one,
                                         int const product_PDG_code_two )
    const
    /* this looks for the given set of PDG codes amongst this instance's
     * decays, & returns the branching ratio associated with the decay if
     * found, or else 0.0 if not.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_branching_ratio_for_exact_match( int const product_PDG_code_one,
                                         int const product_PDG_code_two,
                                         int const product_PDG_code_three )
    const
    /* this looks for the given set of PDG codes amongst this instance's
     * decays, & returns the branching ratio associated with the decay if
     * found, or else 0.0 if not.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_inclusive_branching_ratio_for_particle( int const sought_PDG_code )
    const
    /* this looks for the given PDG code amongst this instance's decays,
     * without caring what the remainder (if any) of the particles in any
     * decay is, & returns the aggregate branching ratio of all these decays.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_branching_ratio_for_exact_match(
                            std::list< int > const* const input_PDG_code_list )
    const
    /* this looks for an exact match to the given set of PDG codes amongst
     * this instance's decays, & returns the branching ratio associated with
     * the decay if found, or else 0.0 if not.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_branching_ratio_for_subset(
                    std::list< int > const* const sought_subset_PDG_code_list )
    const
    /* this looks for the given set of PDG codes for matches to subsets
     * amongst this instance's decays, without caring what the remainder (if
     * any) of the particles in any decay is, & returns the aggregate
     * branching ratio of all these decays.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_branching_ratio_for_subset(
                     std::list< int > const* const sought_subset_PDG_code_list,
                         std::list< int > const* const excluded_PDG_code_list )
    const
    /* this looks for the given set of PDG codes for matches to subsets
     * amongst this instance's decays, checking if any of the particles in the
     * remainder (if any) in each decay appear on the exclusion list, &
     * returns the aggregate branching ratio of all the decays which have the
     * subset & do not have any excluded particles in the remainder.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    int
    merge_identical_decays()
    /* this checks all the decays in the decay handler to see if there are
     * identical sets of decay products, & then if there are, it adds the BR
     * of every identical set after the 1st encountered to the 1st set's BR,
     * & erases the other.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    std::list< particle_decay* > decay_set;
    // I think a std::list is best here, for merge_identical_decays().

  };

  /* instances of this class hold all the relevant data about a type of
   * particle, typically gleaned from a SLHA-format output file from a
   * spectrum generator, such as SPheno. this code assumes CP invariance, so
   * antiparticles have the same masses as their particles, and the same decay
   * widths and branching ratios into the conjugate decays.
   */
  class particle_property_set
  {

  public:

    particle_property_set( int const given_particle_PDG_code,
                           std::string const* const given_particle_name,
                           std::string const* const given_antiparticle_name,
                           bool const given_counts_as_stable_invisible,
                           bool const given_counts_as_jet,
                           bool const given_counts_as_light_lepton,
                           bool const given_counts_as_self_conjugate )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set( int const given_particle_PDG_code,
                           std::string const* const given_particle_name,
                           std::string const* const given_antiparticle_name,
                           std::string const* const given_LaTeX_name,
                           std::string const* const given_LaTeX_antiname,
                           bool const given_counts_as_stable_invisible,
                           bool const given_counts_as_jet,
                           bool const given_counts_as_light_lepton,
                           bool const given_counts_as_self_conjugate )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~particle_property_set()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    become_copy_of( particle_property_set const* const copy_source )
    // this copies all the data from another particle_property_set.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    int
    get_PDG_code()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    bool
    has_mass_been_recorded()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_mass()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_absolute_mass()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set*
    set_mass( double const input_mass )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    bool
    has_decay_width_been_recorded()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_decay_width()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set*
    set_decay_width( double const input_decay_width )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    std::string const*
    get_name()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set*
    set_name( std::string const input_particle_name )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    std::string const*
    get_antiname()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set*
    set_antiname( std::string const input_antiparticle_name )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    std::string const*
    get_LaTeX_name()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set*
    set_LaTeX_name( std::string const input_LaTeX_name )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    std::string const*
    get_LaTeX_antiname()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set*
    set_LaTeX_antiname( std::string const input_LaTeX_antiname )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    std::string const*
    get_name_or_antiname( bool const name_not_antiname )
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    std::list< particle_decay* > const*
    get_direct_decays()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_decay_set_handler const*
    inspect_direct_decay_handler()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_decay_set_handler*
    get_direct_decay_handler()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    std::list< particle_decay* > const*
    get_cascade_decays()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_decay_set_handler const*
    inspect_cascade_decay_handler()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    std::list< particle_decay* >*
    get_cascade_decays()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_decay_set_handler*
    get_cascade_handler()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    bool
    make_count_as_stable_invisible( bool const input_flag )
    // this might have to be set because a sneutrino might be the LSP instead
    // of the lightest neutralino.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    bool
    counts_as_stable_invisible()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    bool
    counts_as_jet()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    bool
    counts_as_light_lepton()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    bool
    counts_as_self_conjugate()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    particle_decay*
    add_direct_decay( int const product_PDG_code_one,
                      int const product_PDG_code_two,
                      double const input_branching_ratio )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_decay*
    add_cascade_decay( int const product_PDG_code_one,
                       int const product_PDG_code_two,
                       double const input_branching_ratio )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_decay*
    add_direct_decay( int const product_PDG_code_one,
                      int const product_PDG_code_two,
                      int const product_PDG_code_three,
                      double const input_branching_ratio )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_decay*
    add_cascade_decay( int const product_PDG_code_one,
                       int const product_PDG_code_two,
                       int const product_PDG_code_three,
                       double const input_branching_ratio )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_decay*
    add_direct_decay( particle_decay const* const decay_to_be_copied )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_decay*
    add_cascade_decay( particle_decay const* const decay_to_be_copied )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    // add_decay for 4- or 5-body decays or even more complicated decays is a
    // bit of a waste of time.


  protected:

    int PDG_code;
    // see codes at the start of this file.  note that antiparticles get
    // negative versions of the particle codes.

    double mass;
    bool mass_has_been_recorded;

    double decay_width;
    bool decay_width_has_been_recorded;

    std::string particle_name;
    // this should be the name of the particle in human-readable form.
    std::string antiparticle_name;
    // this should be the name of the antiparticle in human-readable form.
    std::string LaTeX_name;
    // this is to allow the particle to have a name in LaTeX format.
    std::string LaTeX_antiname;
    // this is to allow the antiparticle to have a name in LaTeX format.

    bool counts_as_stable_invisible_flag;
    // if true, escapes from the LHC detectors, such as neutrinos or
    // lightest neutralinos.

    bool counts_as_jet_flag;
    // if true, basically "is a quark or a gluon"

    bool counts_as_light_lepton_flag;
    // if true, basically "is an electron or a muon"

    bool counts_as_self_conjugate_flag;
    // if true, the particle is also its own antiparticle as far as this code
    // is concerned.

    particle_decay_set_handler direct_decays;
    // this handles all the direct decays of this particle.

    particle_decay_set_handler cascade_decays;
    // this handles all the cascade decays of this particle.

  };



  // inline functions:

  inline std::list< int > const*
  particle_decay::inspect_decay_product_PDG_codes()
  const
  {

    return &decay_product_PDG_codes;

  }

  inline std::list< int >*
  particle_decay::get_decay_product_PDG_codes()
  {

    return &decay_product_PDG_codes;

  }

  inline double
  particle_decay::get_branching_ratio()
  const
  {

    return branching_ratio;

  }

  inline particle_decay*
  particle_decay::set_branching_ratio( double const input_branching_ratio )
  {

    branching_ratio = input_branching_ratio;

    return this;

  }

  inline bool
  particle_decay::has_as_subset(
                          std::list< int > const* sought_subset_PDG_code_list )
  const
  // compares the decay products against a variable amount of given
  // PDG codes.
  {

    std::list< int > empty_exclusion_list;

    return has_as_subset( sought_subset_PDG_code_list,
                          &empty_exclusion_list );

  }



  inline particle_decay*
  particle_decay_set_handler::add_decay( int const product_PDG_code_one,
             int const product_PDG_code_two,
             double const input_branching_ratio )
  /* this creates a new particle_decay with the input, then adds a
   * pointer to it to the vector of pointers to the particle_decay
   * instances of this particle_property_set.
   * N.B.: this does not check that the decay does not already exist for
   * this particle!
   */
  {

    particle_decay*
    pointer_to_decay_to_be_added = new particle_decay( product_PDG_code_one,
                                                       product_PDG_code_two,
                                                       input_branching_ratio );

    decay_set.push_back( pointer_to_decay_to_be_added );

    return pointer_to_decay_to_be_added;

  }

  inline particle_decay*
  particle_decay_set_handler::add_decay( int const product_PDG_code_one,
                                         int const product_PDG_code_two,
                                         int const product_PDG_code_three,
                                         double const input_branching_ratio )
  /* this creates a new particle_decay with the input, then adds a
   * pointer to it to the vector of pointers to the particle_decay
   * instances of this particle_property_set.
   * N.B.: this does not check that the decay does not already exist for
   * this particle!
   */
  {

    particle_decay*
    pointer_to_decay_to_be_added = new particle_decay( product_PDG_code_one,
                                                       product_PDG_code_two,
                                                       product_PDG_code_three,
                                                       input_branching_ratio );

    decay_set.push_back( pointer_to_decay_to_be_added );

    return pointer_to_decay_to_be_added;

  }

  inline particle_decay*
  particle_decay_set_handler::add_decay(
                               particle_decay const* const decay_to_be_copied )
  /* this creates a new particle_decay with the input, then adds a
   * pointer to it to the vector of pointers to the particle_decay
   * instances of this particle_property_set.
   * N.B.: this does not check that the decay does not already exist for
   * this particle!
   */
  {

    particle_decay*
    pointer_to_decay_to_be_added = new particle_decay( decay_to_be_copied );

    decay_set.push_back( pointer_to_decay_to_be_added );

    return pointer_to_decay_to_be_added;

  }

  // add_decay for 4- or 5-body decays or even more complicated decays is a
  // bit of a waste of time.


  inline std::list< particle_decay* > const*
  particle_decay_set_handler::inspect_decay_set()
  const
  {

    return &decay_set;

  }

  inline std::list< particle_decay* >*
  particle_decay_set_handler::get_decay_set()
  {

    return &decay_set;

  }



  inline int
  particle_property_set::get_PDG_code()
  const
  {

    return PDG_code;

  }


  inline bool
  particle_property_set::has_mass_been_recorded()
  const
  {

    return mass_has_been_recorded;

  }

  inline double
  particle_property_set::get_mass()
  const
  {

    if( mass_has_been_recorded )
      {

        return mass;

      }
    else
      {

        std::cout
        << std::endl
        << "CppSLHA::error! a particle mass was asked for which has not been"
        << " recorded.";
        std::cout
        << std::endl <<  CppSLHA_global::really_wrong_value_string
        << " has been returned.";
        std::cout << std::endl;

        return CppSLHA_global::really_wrong_value;

      }

  }

  inline double
  particle_property_set::get_absolute_mass()
  const
  {

    if( mass_has_been_recorded )
      {

        if( 0.0 > mass )
          {

            return -mass;

          }
        else
          {

            return mass;

          }

      }
    else
      {

        std::cout
        << std::endl
        << "CppSLHA::error! a particle mass was asked for which has not been"
        << " recorded.";
        std::cout
        << std::endl <<  CppSLHA_global::really_wrong_value_string
        << " has been returned.";
        std::cout << std::endl;

        return CppSLHA_global::really_wrong_value;

      }

  }

  inline particle_property_set*
  particle_property_set::set_mass( double const input_mass )
  {

    mass = input_mass;
    mass_has_been_recorded = true;

    return this;

  }


  inline bool
  particle_property_set::has_decay_width_been_recorded()
  const
  {

    return decay_width_has_been_recorded;

  }

  inline double
  particle_property_set::get_decay_width()
  const
  {

    if( decay_width_has_been_recorded )
      {

        return decay_width;

      }
    else
      {

        std::cout
        << std::endl
        << "CppSLHA::error! a particle decay width was asked for which has"
        << " not been recorded.";
        std::cout
        << std::endl <<  CppSLHA_global::really_wrong_value_string
        << " has been returned.";
        std::cout << std::endl;

        return CppSLHA_global::really_wrong_value;

      }

  }

  inline particle_property_set*
  particle_property_set::set_decay_width( double const input_decay_width )
  {

    decay_width = input_decay_width;
    decay_width_has_been_recorded = true;

    return this;

  }


  inline std::string const*
  particle_property_set::get_name()
  const
  {

    return &particle_name;

  }

  inline particle_property_set*
  particle_property_set::set_name( std::string const input_particle_name )
  {

    particle_name.assign( input_particle_name );

    return this;

  }

  inline std::string const*
  particle_property_set::get_antiname()
  const
  {

    return &antiparticle_name;

  }

  inline particle_property_set*
  particle_property_set::set_antiname(
                                    std::string const input_antiparticle_name )
  {

    antiparticle_name.assign( input_antiparticle_name );

    return this;

  }

  inline std::string const*
  particle_property_set::get_LaTeX_name()
  const
  {

    return &LaTeX_name;

  }

  inline particle_property_set*
  particle_property_set::set_LaTeX_name( std::string const input_LaTeX_name )
  {

    LaTeX_name.assign( input_LaTeX_name );

    return this;

  }

  inline std::string const*
  particle_property_set::get_LaTeX_antiname()
  const
  {

    return &LaTeX_antiname;

  }

  inline particle_property_set*
  particle_property_set::set_LaTeX_antiname(
                                       std::string const input_LaTeX_antiname )
  {

    LaTeX_antiname.assign( input_LaTeX_antiname );

    return this;

  }

  inline std::string const*
  particle_property_set::get_name_or_antiname( bool const name_not_antiname )
  const
  {

    if( name_not_antiname )
      {

        return &particle_name;

      }
    else
      {

        return &antiparticle_name;

      }

  }

  inline std::list< particle_decay* > const*
  particle_property_set::get_direct_decays()
  const
  {

    return direct_decays.inspect_decay_set();

  }

  inline particle_decay_set_handler const*
  particle_property_set::inspect_direct_decay_handler()
  const
  {

    return &direct_decays;

  }

  inline particle_decay_set_handler*
  particle_property_set::get_direct_decay_handler()
  {

    return &direct_decays;

  }

  inline std::list< particle_decay* > const*
  particle_property_set::get_cascade_decays()
  const
  {

    return cascade_decays.inspect_decay_set();

  }

  inline particle_decay_set_handler const*
  particle_property_set::inspect_cascade_decay_handler()
  const
  {

    return &cascade_decays;

  }

  inline std::list< particle_decay* >*
  particle_property_set::get_cascade_decays()
  {

    return cascade_decays.get_decay_set();

  }

  inline particle_decay_set_handler*
  particle_property_set::get_cascade_handler()
  {

    return &cascade_decays;

  }


  inline bool
  particle_property_set::make_count_as_stable_invisible(
                                                        bool const input_flag )
  // this might have to be set because a sneutrino might be the LSP instead
  // of the lightest neutralino.
  {

    counts_as_stable_invisible_flag = input_flag;

    return counts_as_stable_invisible_flag;

  }

  inline bool
  particle_property_set::counts_as_stable_invisible()
  const
  {

    return counts_as_stable_invisible_flag;

  }

  inline bool
  particle_property_set::counts_as_jet()
  const
  {

    return counts_as_jet_flag;

  }

  inline bool
  particle_property_set::counts_as_light_lepton()
  const
  {

    return counts_as_light_lepton_flag;

  }

  inline bool
  particle_property_set::counts_as_self_conjugate()
  const
  {

    return counts_as_self_conjugate_flag;

  }


  inline particle_decay*
  particle_property_set::add_direct_decay( int const product_PDG_code_one,
                                           int const product_PDG_code_two,
                                           double const input_branching_ratio )
  {

    return direct_decays.add_decay( product_PDG_code_one,
                                    product_PDG_code_two,
                                    input_branching_ratio );

  }

  inline particle_decay*
  particle_property_set::add_cascade_decay( int const product_PDG_code_one,
                                            int const product_PDG_code_two,
                                           double const input_branching_ratio )
  {

    return cascade_decays.add_decay( product_PDG_code_one,
                                     product_PDG_code_two,
                                     input_branching_ratio );

  }

  inline particle_decay*
  particle_property_set::add_direct_decay( int const product_PDG_code_one,
                                           int const product_PDG_code_two,
                                           int const product_PDG_code_three,
                                           double const input_branching_ratio )
  {

    return direct_decays.add_decay( product_PDG_code_one,
                                    product_PDG_code_two,
                                    product_PDG_code_three,
                                    input_branching_ratio );

  }

  inline particle_decay*
  particle_property_set::add_cascade_decay( int const product_PDG_code_one,
                                            int const product_PDG_code_two,
                                            int const product_PDG_code_three,
                                           double const input_branching_ratio )
  {

    return cascade_decays.add_decay( product_PDG_code_one,
                                     product_PDG_code_two,
                                     product_PDG_code_three,
                                     input_branching_ratio );

  }

  inline particle_decay*
  particle_property_set::add_direct_decay(
                               particle_decay const* const decay_to_be_copied )
  {

    return direct_decays.add_decay( decay_to_be_copied );

  }

  inline particle_decay*
  particle_property_set::add_cascade_decay(
                               particle_decay const* const decay_to_be_copied )
  {

    return cascade_decays.add_decay( decay_to_be_copied );

  }

} // end of CppSLHA namespace.

#endif /* CPPSLHA_PARTICLE_DATA_HPP_ */
