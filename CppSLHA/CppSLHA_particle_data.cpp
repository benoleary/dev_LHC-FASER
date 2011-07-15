/*
 * CppSLHA_particle_data.cpp
 *
 *  Created on: 12 Jan 2010
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

#include "CppSLHA_particle_data.hpp"

namespace CppSLHA
{

  // constructor for a 2-body decay:
  particle_decay::particle_decay( int const product_PDG_code_one,
                                  int const product_PDG_code_two,
                                  double const input_branching_ratio ) :
    branching_ratio( input_branching_ratio )
  {

    decay_product_PDG_codes.push_back( product_PDG_code_one );
    decay_product_PDG_codes.push_back( product_PDG_code_two );

  }

  // constructor for a 3-body decay:
  particle_decay::particle_decay( int const product_PDG_code_one,
                                  int const product_PDG_code_two,
                                  int const product_PDG_code_three,
                                  double const input_branching_ratio ) :
    branching_ratio( input_branching_ratio )
  {

    decay_product_PDG_codes.push_back( product_PDG_code_one );
    decay_product_PDG_codes.push_back( product_PDG_code_two );
    decay_product_PDG_codes.push_back( product_PDG_code_three );

  }

  // further constructors (4-body, 5-body, & so forth) are a bit of a waste of
  // time.

  particle_decay::particle_decay(
                         particle_decay const* const particle_decay_to_copy ) :
    branching_ratio( particle_decay_to_copy->get_branching_ratio() )
  {

    for( std::list< int >::const_iterator
         decay_PDG_codes_iterator
         = particle_decay_to_copy->inspect_decay_product_PDG_codes()->begin();
         particle_decay_to_copy->inspect_decay_product_PDG_codes()->end()
         != decay_PDG_codes_iterator;
         ++decay_PDG_codes_iterator )
      {

        decay_product_PDG_codes.push_back( *decay_PDG_codes_iterator );

      }

  }

  particle_decay::~particle_decay()
  {

    // does nothing.

  }


  bool
  particle_decay::has_this_PDG_code( int const input_PDG_code )
  const
  // checks to see if the sought particle is among the decay products.
  {

    bool return_bool = false;
    // assume that the particle is not in the set of decay products.

    // check each of the decay products, & if it has the sought PDG code,
    // note that.
    for( std::list< int >::const_iterator
           list_searcher = decay_product_PDG_codes.begin();
         decay_product_PDG_codes.end() != list_searcher;
         ++list_searcher )
      {

        if( input_PDG_code == *list_searcher )
          {

            return_bool = true;

          }

      }

    return return_bool;

  }


  bool
  particle_decay::is_same_as( int const product_PDG_code_one,
                              int const product_PDG_code_two )
  const
  // compares the decay products against a pair of given PDG codes.
  {

    bool comparison_to_return = false;
    // assume that the input list does not match the decay product list.

    if( 2 == decay_product_PDG_codes.size() )
    // if at least the decay product list is the right length...
      {

        // make a small list of the input codes & pass it to the general case
        // is_same_as():
        std::list< int > input_PDG_codes_as_list;
        input_PDG_codes_as_list.push_back( product_PDG_code_one );
        input_PDG_codes_as_list.push_back( product_PDG_code_two );

        comparison_to_return = is_same_as( &input_PDG_codes_as_list );

      }

    return comparison_to_return;

  }

  bool
  particle_decay::is_same_as( int const product_PDG_code_one,
                              int const product_PDG_code_two,
                              int const product_PDG_code_three )
  const
  // compares the decay products against a set of 3 given PDG codes.
  {

    bool comparison_to_return = false;
    // assume that the input list does not match the decay product list.

    if( 3 == decay_product_PDG_codes.size() )
      // if the decay product list is at least the right length...
      {

        // make a small list of the input codes & pass it to the general case
        // is_same_as().
        std::list< int > input_PDG_codes_as_list;
        input_PDG_codes_as_list.push_back( product_PDG_code_one );
        input_PDG_codes_as_list.push_back( product_PDG_code_two );
        input_PDG_codes_as_list.push_back( product_PDG_code_three );

        comparison_to_return = is_same_as( &input_PDG_codes_as_list );

      }

    return comparison_to_return;

  }

  bool
  particle_decay::is_same_as(
                            std::list< int > const* const input_PDG_code_list )
  const
  // compares the decay products against a variable amount of given PDG codes.
  {

    bool comparison_to_return = false;
    // assume that the input list does not match the decay product list.

    if( ( input_PDG_code_list->size() == decay_product_PDG_codes.size() )
        &&
        has_as_subset( input_PDG_code_list ) )
    // if the numbers match up & the given set is a subset, they are the same.
      {

        comparison_to_return = true;

      }

    return comparison_to_return;

  }


  bool
  particle_decay::has_as_subset(
                     std::list< int > const* const sought_subset_PDG_code_list,
                         std::list< int > const* const excluded_PDG_code_list )
  const
  // compares the decay products against a variable amount of given PDG codes,
  // along with a list of excluded PDG codes.
  {

    bool comparison_to_return = false;
    // assume that the input list does not match the decay product list.

    if( sought_subset_PDG_code_list->size()
          <= decay_product_PDG_codes.size() )
    // check that the numbers of PDG codes to check match up.
      {

        // if the input list is at least not longer than this decay, start by
        // assuming that the subset will be found:
        comparison_to_return = true;

        // create a list of PDG codes which can have codes removed, to prevent
        // confusion over comparison of duplicates in either list.
        std::list< int > list_of_PDG_codes_to_check_off;
        for( std::list< int >::const_iterator
               list_filler = decay_product_PDG_codes.begin();
             decay_product_PDG_codes.end() != list_filler;
             ++list_filler )
          {

            list_of_PDG_codes_to_check_off.push_back( *list_filler );

          }

        /* loop (backwards) through all the provided PDG codes, checking that
         * each is accounted for once, removing checked PDG codes from the
         * list.
         */
        std::list< int >::const_iterator
        sought_subset_PDG_code_list_iterator
          = sought_subset_PDG_code_list->begin();
        while( comparison_to_return
               &&
               ( sought_subset_PDG_code_list->end()
                   != sought_subset_PDG_code_list_iterator ) )
          {

            // search through the list until either the sought PDG code is
            // found, or the end is reached.
            bool PDG_code_unfound_flag = true;
            std::list< int >::iterator
            list_searcher = list_of_PDG_codes_to_check_off.begin();
            while( PDG_code_unfound_flag
                   &&
                   ( list_of_PDG_codes_to_check_off.end()
                       != list_searcher ) )
              {

                // if the PDG code is found in the list, remove it & note that
                // it was found:
                if( *list_searcher == *sought_subset_PDG_code_list_iterator )
                  {

                    PDG_code_unfound_flag = false;

                    list_searcher
                    = list_of_PDG_codes_to_check_off.erase( list_searcher );

                  }
                else // otherwise, go on to the next element.
                  {

                    ++list_searcher;

                  }

              }

            /* if the previous loop did not check off the sought PDG code,
             * note that it was not found so the sets of PDG codes do not
             * match, so this method should return false.
             */
            if( PDG_code_unfound_flag )
              {

                comparison_to_return = false;

              }

            ++sought_subset_PDG_code_list_iterator;
            // go on to the next sought PDG code.

          }

        /* now go through the list of remaining PDG codes (the excess to the
         * list sought), checking to see if any are on the list of
         * excluded PDG codes.
         */
        if( comparison_to_return
            &&
            !( excluded_PDG_code_list->empty() ) )
        // if the subset was found, & there are PDG codes to exclude...
          {

            std::list< int >::const_iterator
            excluded_PDG_code_list_iterator = excluded_PDG_code_list->begin();
            while( comparison_to_return
                   &&
                   ( excluded_PDG_code_list->end()
                       != excluded_PDG_code_list_iterator ) )
              {

                // search through the list until either the excluded PDG code
                // is found, or the end is reached:
                std::list< int >::iterator
                list_searcher = list_of_PDG_codes_to_check_off.begin();
                while( comparison_to_return
                       &&
                       ( list_of_PDG_codes_to_check_off.end()
                           != list_searcher ) )
                  {

                    if( *excluded_PDG_code_list_iterator == *list_searcher )
                    // if the remainder list has an excluded particle...
                      {

                        comparison_to_return = false;
                        // note this.

                      }
                    else
                    // otherwise, go on to the next element.
                      {

                        ++list_searcher;

                      }

                  }

                ++excluded_PDG_code_list_iterator;
                // go on to the next sought PDG code.

              }

          }

      }

    return comparison_to_return;

  }



  particle_decay_set_handler::particle_decay_set_handler()
  {

    // does nothing.

  }

  particle_decay_set_handler::~particle_decay_set_handler()
  {

    for( std::list< particle_decay* >::iterator
           list_destructor = decay_set.begin();
         decay_set.end() != list_destructor;
         ++list_destructor )
      // go through each pointer in the std::list.
      {

        delete *list_destructor;

      }

  }


  particle_decay_set_handler*
  particle_decay_set_handler::become_copy_of(
                particle_decay_set_handler const* const decay_handler_to_copy )
  // this copies a particle_decay_set_handler instance.
  {

    // 1st, empty the existing decays:
    for( std::list< particle_decay* >::iterator
           list_destructor = decay_set.begin();
         decay_set.end() != list_destructor;
         ++list_destructor )
      // go through each pointer in the std::list.
      {

        // delete the allocated memory pointed at by each pointer:
        delete *list_destructor;

      }
    decay_set.clear();

    particle_decay* pointer_to_decay_to_be_added;

    for( std::list< particle_decay* >::const_iterator
           decay_copying_iterator
             = decay_handler_to_copy->inspect_decay_set()->begin();
         decay_handler_to_copy->inspect_decay_set()->end()
           != decay_copying_iterator;
         ++decay_copying_iterator )
      {

        pointer_to_decay_to_be_added
        = new particle_decay( *decay_copying_iterator );

        decay_set.push_back( pointer_to_decay_to_be_added );

      }

    return this;

  }

  std::list< particle_decay* >*
  particle_decay_set_handler::copy_decay_set(
                 std::list< particle_decay* > const* const decay_sets_to_copy )
  // this copies a particle_decay_set_handler instance.
  {

    particle_decay* pointer_to_decay_to_be_added;

    for( std::list< particle_decay* >::const_iterator
           decay_copying_iterator = decay_sets_to_copy->begin();
         decay_sets_to_copy->end()
           != decay_copying_iterator;
         ++decay_copying_iterator )
      {

        pointer_to_decay_to_be_added
        = new particle_decay( *decay_copying_iterator );

        decay_set.push_back( pointer_to_decay_to_be_added );

      }

    return &decay_set;

  }


  double
  particle_decay_set_handler::get_branching_ratio_for_exact_match(
                                               int const product_PDG_code_one,
                                               int const product_PDG_code_two )
  const
  /* this looks for the given set of PDG codes amongst this instance's decays,
   * & returns the branching ratio associated with the decay if found, or else
   * 0.0 if not.
   */
  {

    double branching_ratio_to_return = 0.0;
    // start by assuming that the decay won't be found.

    for( std::list< particle_decay* >::const_iterator
           kept_decay_iterator = decay_set.begin();
         decay_set.end() != kept_decay_iterator;
         ++kept_decay_iterator )
      // go through each stored decay.
      {

        if( (*kept_decay_iterator)->is_same_as( product_PDG_code_one,
                                                product_PDG_code_two ) )
          // if this decay matches the sought decay...
          {

            branching_ratio_to_return
              += (*kept_decay_iterator)->get_branching_ratio();
            // note the branching ratio to return (allowing for multiple
            // occurrences of the decay in the set of decays).

          }

      }

    return branching_ratio_to_return;

  }

  double
  particle_decay_set_handler::get_branching_ratio_for_exact_match(
                                             int const product_PDG_code_one,
                                             int const product_PDG_code_two,
                                             int const product_PDG_code_three )
  const
  /* this looks for the given set of PDG codes amongst this instance's decays,
   * & returns the branching ratio associated with the decay if found, or else
   * 0.0 if not.
   */
  {

    double branching_ratio_to_return = 0.0;
    // start by assuming that the decay won't be found.

    for( std::list< particle_decay* >::const_iterator
           kept_decay_iterator = decay_set.begin();
         decay_set.end() != kept_decay_iterator;
         ++kept_decay_iterator )
      // go through each stored decay.
      {

        if( (*kept_decay_iterator)->is_same_as( product_PDG_code_one,
                                                product_PDG_code_two,
                                                product_PDG_code_three ) )
          // if this decay matches the sought decay...
          {

            branching_ratio_to_return
                += (*kept_decay_iterator)->get_branching_ratio();
            // note the branching ratio to return (allowing for multiple
            // occurrences of the decay in the set of decays).

          }

      }

    return branching_ratio_to_return;

  }

  double
  particle_decay_set_handler::get_branching_ratio_for_exact_match(
                            std::list< int > const* const input_PDG_code_list )
  const
  /* this looks for an exact match to the given set of PDG codes amongst this
   * instance's decays, & returns the branching ratio associated with the
   * decay if found, or else 0.0 if not.
   */
  {

    double branching_ratio_to_return = 0.0;
    // start by assuming that the decay won't be found.

    for( std::list< particle_decay* >::const_iterator
           kept_decay_iterator = decay_set.begin();
         decay_set.end() != kept_decay_iterator;
         ++kept_decay_iterator )
      // go through each stored decay.
      {

        if( (*kept_decay_iterator)->is_same_as( input_PDG_code_list ) )
          // if this decay matches the sought decay...
          {

            branching_ratio_to_return
                += (*kept_decay_iterator)->get_branching_ratio();
            // note the branching ratio to return (allowing for multiple
            // occurrences of the decay in the set of decays).

          }

      }

    return branching_ratio_to_return;

  }

  double
  particle_decay_set_handler::get_inclusive_branching_ratio_for_particle(
                                                    int const sought_PDG_code )
  const
  /* this looks for the given PDG code amongst this instance's decays, without
   * caring what the remainder (if any) of the particles in any decay is, &
   * returns the aggregate branching ratio of all these decays.
   */
  {

    double branching_ratio_to_return = 0.0;
    // start by assuming that the decay won't be found.

    for( std::list< particle_decay* >::const_iterator
           kept_decay_iterator = decay_set.begin();
         decay_set.end() != kept_decay_iterator;
         ++kept_decay_iterator )
    // go through each stored decay.
      {

        if( (*kept_decay_iterator)->has_this_PDG_code( sought_PDG_code ) )
        // if this decay has the sought PDG code...
          {

            branching_ratio_to_return
                += (*kept_decay_iterator)->get_branching_ratio();
            // update the branching ratio to return.

          }

      }

    return branching_ratio_to_return;

  }

  double
  particle_decay_set_handler::get_branching_ratio_for_subset(
                    std::list< int > const* const sought_subset_PDG_code_list )
  const
  /* this looks for the given set of PDG codes for matches to subsets amongst
   * this instance's decays, without caring what the remainder (if any) of the
   * particles in any decay is, & returns the aggregate branching ratio of all
   * these decays.
   */
  {

    double branching_ratio_to_return = 0.0;
    // start by assuming that the decay won't be found.

    for( std::list< particle_decay* >::const_iterator
           kept_decay_iterator = decay_set.begin();
         decay_set.end() != kept_decay_iterator;
         ++kept_decay_iterator )
      // go through each stored decay.
      {

        if( (*kept_decay_iterator)->has_as_subset(
                                      sought_subset_PDG_code_list ) )
          // if this decay has the sought PDG codes as a subset...
          {

            branching_ratio_to_return
                += (*kept_decay_iterator)->get_branching_ratio();
            // update the branching ratio to return.

          }

      }

    return branching_ratio_to_return;

  }

  double
  particle_decay_set_handler::get_branching_ratio_for_subset(
                     std::list< int > const* const sought_subset_PDG_code_list,
                         std::list< int > const* const excluded_PDG_code_list )
  const
  /* this looks for the given set of PDG codes for matches to subsets amongst
   * this instance's decays, checking if any of the particles in the remainder
   * (if any) in each decay appear on the exclusion list, & returns the
   * aggregate branching ratio of all the decays which have the subset & do
   * not have any excluded particles in the remainder.
   */
  {

    double branching_ratio_to_return = 0.0;
    // start by assuming that the decay won't be found.

    for( std::list< particle_decay* >::const_iterator
           kept_decay_iterator = decay_set.begin();
         decay_set.end() != kept_decay_iterator;
         ++kept_decay_iterator )
      // go through each stored decay.
      {

        if( (*kept_decay_iterator)->has_as_subset( sought_subset_PDG_code_list,
                                                   excluded_PDG_code_list ) )
          // if this decay has the sought PDG codes as a subset & none of the
          // excluded PDG codes...
          {

            branching_ratio_to_return
                += (*kept_decay_iterator)->get_branching_ratio();
            // update the branching ratio to return.

          }

      }

    return branching_ratio_to_return;

  }


  int
  particle_decay_set_handler::merge_identical_decays()
  /* this checks all the decays in the decay handler to see if there are
   * identical sets of decay products, & then if there are, it adds the BR of
   * every identical set after the 1st encountered to the 1st set's BR, &
   * erases the other.
   */
  {

    int number_of_mergers = 0;

    for( std::list< particle_decay* >::iterator
           kept_decay_iterator = decay_set.begin();
         decay_set.end() != kept_decay_iterator;
         ++kept_decay_iterator )
      // go through each stored decay.
      {

        /* go through each decay after (backwards, so actually before) the one
         * that will be kept. this would be a for loop, except we have to
         * initialize an iterator & advance it once before starting the loop.
         */
        std::list< particle_decay* >::iterator
        checked_decay_iterator = kept_decay_iterator;
        if( decay_set.end() != checked_decay_iterator )
          // if there are still decays to check...
          {

            ++checked_decay_iterator;
            // don't check the decay against itself.

          }
        while( decay_set.end() != checked_decay_iterator )
          // before we run out of decays to check against...
          {

            if( (*checked_decay_iterator)->is_same_as(
                  (*kept_decay_iterator)->get_decay_product_PDG_codes() ) )
              // if the decays have the same decay products...
              {

                (*kept_decay_iterator)->set_branching_ratio( (
                    (*kept_decay_iterator)->get_branching_ratio()
                    + (*checked_decay_iterator)->get_branching_ratio() ) );
                // update the BR of the decay that will be kept.

                checked_decay_iterator
                    = decay_set.erase( checked_decay_iterator );
                // erase the other decay.

                ++number_of_mergers;

              }
            else
              {

                ++checked_decay_iterator;
                // go on to the next decay in the list to be checked.

              }

          }

      }

    return number_of_mergers;

  }



  particle_property_set::particle_property_set( int const given_PDG_code,
                                  std::string const* const given_particle_name,
                              std::string const* const given_antiparticle_name,
                              bool const given_counts_as_stable_invisible_flag,
                                           bool const given_counts_as_jet_flag,
                                  bool const given_counts_as_light_lepton_flag,
                             bool const given_counts_as_self_conjugate_flag ) :
    PDG_code( given_PDG_code ),
    mass( CppSLHA_global::really_wrong_value ),
    mass_has_been_recorded( false ),
    /* this is to note that a mass has not been given for this particle
     * property set, as opposed to guessing whether it should be really
     * massless, or if its mass was never given.
     */
    decay_width( CppSLHA_global::really_wrong_value ),
    decay_width_has_been_recorded( false ),
    /* this is to note that a decay width has not been given for this particle
     * property set, as opposed to guessing whether it should be really
     * stable, or if its decay width was never given.
     */
    particle_name( *given_particle_name ),
    antiparticle_name( *given_antiparticle_name ),
    counts_as_stable_invisible_flag( given_counts_as_stable_invisible_flag ),
    counts_as_jet_flag( given_counts_as_jet_flag ),
    counts_as_light_lepton_flag( given_counts_as_light_lepton_flag ),
    counts_as_self_conjugate_flag( given_counts_as_self_conjugate_flag )
  {

    if( 0 > given_PDG_code )
      // if the particle was defined with a negative PDG code, record the
      // antiparticle with positive PDG code.
      {

        PDG_code = -given_PDG_code;
        particle_name.assign( *given_antiparticle_name );
        antiparticle_name.assign( *given_particle_name );

      }
    std::stringstream code_to_string;
    code_to_string << PDG_code;
    if( 0 > given_PDG_code )
      // if the particle was defined with a negative PDG code...
      {

        // record the antiparticle name as the positive PDG code.
        LaTeX_antiname.assign( code_to_string.str() );
        LaTeX_name.assign( "-" );
        LaTeX_name.append( code_to_string.str() );

      }
    else
      {

        LaTeX_name.assign( code_to_string.str() );
        LaTeX_antiname.assign( "-" );
        LaTeX_antiname.append( code_to_string.str() );

      }

  }

  particle_property_set::particle_property_set( int const given_PDG_code,
                                  std::string const* const given_particle_name,
                              std::string const* const given_antiparticle_name,
                                     std::string const* const given_LaTeX_name,
                                 std::string const* const given_LaTeX_antiname,
                              bool const given_counts_as_stable_invisible_flag,
                                           bool const given_counts_as_jet_flag,
                                  bool const given_counts_as_light_lepton_flag,
                             bool const given_counts_as_self_conjugate_flag ) :
    PDG_code( given_PDG_code ),
    mass( CppSLHA_global::really_wrong_value ),
    mass_has_been_recorded( false ),
    /* this is to note that a mass has not been given for this particle
     * property set, as opposed to guessing whether it should be really
     * massless, or if its mass was never given.
     */
    decay_width( CppSLHA_global::really_wrong_value ),
    decay_width_has_been_recorded( false ),
    /* this is to note that a decay width has not been given for this particle
     * property set, as opposed to guessing whether it should be really
     * stable, or if its decay width was never given.
     */
    particle_name( *given_particle_name ),
    antiparticle_name( *given_antiparticle_name ),
    LaTeX_name( *given_LaTeX_name ),
    LaTeX_antiname( *given_LaTeX_antiname ),
    counts_as_stable_invisible_flag( given_counts_as_stable_invisible_flag ),
    counts_as_jet_flag( given_counts_as_jet_flag ),
    counts_as_light_lepton_flag( given_counts_as_light_lepton_flag ),
    counts_as_self_conjugate_flag( given_counts_as_self_conjugate_flag )
  {

    if( 0 > given_PDG_code )
      // if the particle was defined with a negative PDG code, record the
      // antiparticle with positive PDG code.
      {

        PDG_code = -given_PDG_code;
        particle_name.assign( *given_antiparticle_name );
        antiparticle_name.assign( *given_particle_name );
        LaTeX_name.assign( *given_LaTeX_antiname );
        LaTeX_antiname.assign( *given_LaTeX_name );

      }

  }

  particle_property_set::~particle_property_set()
  {

    // does nothing.

  }


  void
  particle_property_set::become_copy_of(
                               particle_property_set const* const copy_source )
  {

    PDG_code = copy_source->get_PDG_code();

    mass_has_been_recorded = copy_source->has_mass_been_recorded();
    if( copy_source->has_mass_been_recorded() )
      {

        mass = copy_source->get_mass();

      }
    else
      {

        mass = CppSLHA_global::really_wrong_value;

      }

    particle_name.assign( *(copy_source->get_name()) );
    antiparticle_name.assign( *(copy_source->get_antiname()) );
    LaTeX_name.assign( *(copy_source->get_LaTeX_name()) );
    LaTeX_antiname.assign( *(copy_source->get_LaTeX_antiname()) );

    decay_width_has_been_recorded
    = copy_source->has_decay_width_been_recorded();
    if( copy_source->has_decay_width_been_recorded() )
      {

        decay_width = copy_source->get_decay_width();

      }
    else
      {

        decay_width = CppSLHA_global::really_wrong_value;

      }

    counts_as_stable_invisible_flag
    = copy_source->counts_as_stable_invisible();
    counts_as_jet_flag = copy_source->counts_as_jet();
    counts_as_light_lepton_flag = copy_source->counts_as_light_lepton();
    counts_as_self_conjugate_flag = copy_source->counts_as_self_conjugate();

    direct_decays.become_copy_of(
                                 copy_source->inspect_direct_decay_handler() );
    cascade_decays.become_copy_of(
                                copy_source->inspect_cascade_decay_handler() );

  }

} // end of CppSLHA namespace.

