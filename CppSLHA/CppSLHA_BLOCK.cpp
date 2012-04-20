/*
 * CppSLHA_BLOCK.cpp
 *
 *  Created on: 07 Oct 2010
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

#include "CppSLHA_BLOCK.hpp"


namespace CppSLHA
{

  indexed_double::indexed_double( int const given_index_int,
                                  double const given_value_double ) :
    index_int( given_index_int ),
    value_double( given_value_double )
  {

    // just an initialization list.

  }

  indexed_double::~indexed_double()
  {

    // does nothing.

  }



  indexed_double_set::indexed_double_set()
  {

    // does nothing.

  }

  indexed_double_set::~indexed_double_set()
  {

    for( std::vector< indexed_double* >::iterator
           deletion_iterator = entries.begin();
         entries.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

  }


  bool
  indexed_double_set::has_entry( int const requested_index )
  const
  // this checks to see if the entry exists in this set.
  {

    bool return_bool = false;
    // assume that the requested entry will not be found.

    for( std::vector< indexed_double* >::const_iterator
           entry_iterator = entries.begin();
         entries.end() > entry_iterator;
         ++entry_iterator )
      {

        if( requested_index == (*entry_iterator)->get_index() )
        // if the requested number is found...
          {

            return_bool = true;
            // note that the associated number has been found.

            entry_iterator = entries.end();
            // stop looking for the number now that it has been found.

          }

      }

    return return_bool;
    // return true if requested_first_number was found, false otherwise.

  }

  bool
  indexed_double_set::has_entry( int const requested_index,
                                 double* const double_for_entry )
  const
  // this checks to see if the entry exists in this set & puts its value into
  // the provided memory address for a double.
  {

    bool return_bool = false;
    // assume that the requested entry will not be found.

    for( std::vector< indexed_double* >::const_iterator
           entry_iterator = entries.begin();
         entries.end() > entry_iterator;
         ++entry_iterator )
      {

        if( requested_index == (*entry_iterator)->get_index() )
        // if the requested number is found...
          {

            return_bool = true;
            // note that the associated number has been found.

            *double_for_entry = (*entry_iterator)->get_double();
            // record the entry in the given memory address.

            entry_iterator = entries.end();
            // stop looking for the number now that it has been found.

          }

      }

    return return_bool;
    // return true if requested_first_number was found, false otherwise.

  }

  double
  indexed_double_set::get_entry( int const requested_index )
  const
  {

    bool entry_found = false;

    double return_number = CppSLHA_global::really_wrong_value;
    // assume that the requested entry will not be found.

    for( std::vector< indexed_double* >::const_iterator
           entry_iterator = entries.begin();
         entries.end() > entry_iterator;
         ++entry_iterator )
      {

        if( requested_index == (*entry_iterator)->get_index() )
        // if the requested number is found...
          {

            return_number = (*entry_iterator)->get_double();
            // note the found associated number.

            entry_found = true;
            // note that the entry was found.

            entry_iterator = entries.end();
            // stop looking for the number now that it has been found.

          }

      }

    if( entry_found )
    // if the entry was found...
      {

        return return_number;

      }
    else
    // otherwise the entry was not found...
      {

      /** I've decided to turn off these messages.
        std::cout
        << std::endl
        << "CppSLHA::error! an SLHA_BLOCK was asked for an entry which it had"
        << " not recorded.";
        std::cout
        << std::endl
        << CppSLHA_global::really_wrong_value_string
        << " has been returned.";
        std::cout << std::endl;**/

        return CppSLHA_global::really_wrong_value;

      }

  }

  void
  indexed_double_set::set_entry( int const requested_index,
                                 double const given_double )
  {

    bool entry_needs_to_be_added = true;

    for( std::vector< indexed_double* >::iterator
           entry_iterator = entries.begin();
         entries.end() > entry_iterator;
         ++entry_iterator )
      {

        if( requested_index == (*entry_iterator)->get_index() )
        // if the requested number is found...
          {

            (*entry_iterator)->set_double( given_double );
            // note the found associated number.

            entry_needs_to_be_added = false;
            // note that the entry does not need to be added.

            entry_iterator = entries.end();
            // stop looking for the number now that it has been found.

          }

      }

    if( entry_needs_to_be_added )
    // if requested_first_number was not found...
      {

        indexed_double*
        indexed_double_adder = new indexed_double( requested_index,
                                                   given_double );

        entries.push_back( indexed_double_adder );

      }

    if( 1 >= entries.size() )
      {

        lowest_index = requested_index;
        highest_index = requested_index;

      }
    else
      {

        if( requested_index < lowest_index )
          {

            lowest_index = requested_index;

          }

        if( requested_index > highest_index )
          {

            highest_index = requested_index;

          }

      }

  }



  single_indexed_comment::single_indexed_comment( int const given_index_int,
                              std::string const* const given_comment_string ) :
    index_int( given_index_int ),
    comment_string( *given_comment_string )
  {

    // just an initialization list.

  }

  single_indexed_comment::~single_indexed_comment()
  {

    // does nothing.

  }



  pair_indexed_comment::pair_indexed_comment( int const given_first_index_int,
                                              int const given_second_index_int,
                              std::string const* const given_comment_string ) :
    first_index_int( given_first_index_int ),
    second_index_int( given_second_index_int ),
    comment_string( *given_comment_string )
  {

    // just an initialization list.

  }

  pair_indexed_comment::~pair_indexed_comment()
  {

    // does nothing.

  }



  triplet_indexed_comment::triplet_indexed_comment(
                                               int const given_first_index_int,
                                              int const given_second_index_int,
                                               int const given_third_index_int,
                                std::string const* const given_comment_string ) :
    first_index_int( given_first_index_int ),
    second_index_int( given_second_index_int ),
    third_index_int( given_third_index_int ),
    comment_string( *given_comment_string )
  {

    // just an initialization list.

  }

  triplet_indexed_comment::~triplet_indexed_comment()
  {

    // does nothing.

  }



  comment_organizer::comment_organizer()
  {

    // does nothing.

  }

  comment_organizer::~comment_organizer()
  {

    // does nothing.

  }



  one_index_comment_organizer::one_index_comment_organizer() :
    comment_organizer()
  {

    // does nothing.

  }

  one_index_comment_organizer::~one_index_comment_organizer()
  {

    for( std::vector< single_indexed_comment* >::iterator
           deletion_iterator = comment_strings.begin();
         comment_strings.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

  }


  bool
  one_index_comment_organizer::has_comment( int const requested_index,
                                        std::string* const string_for_comment )
  const
  // this looks for the comment for the requested index, & if it is found,
  // writes it to the provided string, while returning true.
  {

    bool return_bool = false;
    // assume that the requested entry will not be found.

    for( std::vector< single_indexed_comment* >::const_iterator
           comment_iterator = comment_strings.begin();
         comment_strings.end() > comment_iterator;
         ++comment_iterator )
      {

        if( requested_index == (*comment_iterator)->get_index() )
        // if the requested index is found...
          {

            return_bool = true;
            // note that the associated comment has been found.

            string_for_comment->assign(
                                       *((*comment_iterator)->get_comment()) );
            // record the comment in the given memory address.

            comment_iterator = comment_strings.end();
            // stop looking for the index now that it has been found.

          }

      }

    return return_bool;
    // return true if requested_first_number was found, false otherwise.

  }

  void
  one_index_comment_organizer::set_comment( int const requested_index,
                                       std::string const* const given_comment )
  {

    // debugging:
    /*std::cout
     << std::endl
     << "debugging: one_index_comment_organizer::set_comment( "
     << requested_index << ", " <<  given_comment << " ) called";
     std::cout << std::endl;
     */

    bool comment_needs_to_be_added = true;

    for( std::vector< single_indexed_comment* >::iterator
           comment_iterator = comment_strings.begin();
         comment_strings.end() > comment_iterator;
         ++comment_iterator )
      {

        if( requested_index == (*comment_iterator)->get_index() )
        // if the requested index is found...
          {

            (*comment_iterator)->set_comment( given_comment );
            // set the associated comment.

            comment_needs_to_be_added = false;
            // note that the comment does not need to be added.

            comment_iterator = comment_strings.end();
            // stop looking for the index now that it has been found.

          }

      }

    if( comment_needs_to_be_added )
    // if requested_index was not found...
      {

        // debugging:
        /*std::cout << std::endl << "debugging: have to add a comment.";
         std::cout << std::endl;
         */

        // create an indexed comment & add it to the set:
        single_indexed_comment*
        indexed_comment_adder
        = new single_indexed_comment( requested_index,
                                      given_comment );

        comment_strings.push_back( indexed_comment_adder );

      }

  }



  two_index_comment_organizer::two_index_comment_organizer() :
    comment_organizer()
  {

    // does nothing.

  }

  two_index_comment_organizer::~two_index_comment_organizer()
  {

    for( std::vector< pair_indexed_comment* >::iterator
           deletion_iterator = comment_strings.begin();
         comment_strings.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;
        *deletion_iterator = NULL;  // just to be safe.

      }

  }


  bool
  two_index_comment_organizer::has_comment( int const requested_first_index,
                                            int const requested_second_index,
                                        std::string* const string_for_comment )
  const
  // this looks for the comment for the requested indices, & if it is found,
  // writes it to the provided string, while returning true.
  {

    bool return_bool = false;
    // assume that the requested entry will not be found.

    for( std::vector< pair_indexed_comment* >::const_iterator
           comment_iterator = comment_strings.begin();
         comment_strings.end() > comment_iterator;
         ++comment_iterator )
      {

        if( ( requested_first_index == (*comment_iterator)->get_first_index() )
            &&
            ( requested_second_index
                == (*comment_iterator)->get_second_index() ) )
        // if the requested indices are found...
          {

            return_bool = true;
            // note that the associated comment has been found.

            string_for_comment->assign(
                                       *((*comment_iterator)->get_comment()) );
            // record the comment in the given memory address.

            comment_iterator = comment_strings.end();
            // stop looking for the indices now that they has been found.

          }

      }

    return return_bool;
    // return true if requested_first_number was found, false otherwise.

  }

  void
  two_index_comment_organizer::set_comment( int const requested_first_index,
                                            int const requested_second_index,
                                       std::string const* const given_comment )
  {

    bool comment_needs_to_be_added = true;

    for( std::vector< pair_indexed_comment* >::iterator
           comment_iterator = comment_strings.begin();
         comment_strings.end() > comment_iterator;
         ++comment_iterator )
      {

        if( ( requested_first_index == (*comment_iterator)->get_first_index() )
            &&
            ( requested_second_index
                == (*comment_iterator)->get_second_index() ) )
        // if the requested indices are found...
          {

            (*comment_iterator)->set_comment( given_comment );
            // set the requested comment.

            comment_needs_to_be_added = false;
            // note that the comment does not need to be added.

            comment_iterator = comment_strings.end();
            // stop looking for the index now that it has been found.

          }

      }

    if( comment_needs_to_be_added )
    // if the requested indices were not found...
      {

        // create an indexed comment & add it to the set:
        pair_indexed_comment*
        indexed_comment_adder =
                new pair_indexed_comment( requested_first_index,
                                          requested_second_index,
                                          given_comment );

        comment_strings.push_back( indexed_comment_adder );

      }

  }



  three_index_comment_organizer::three_index_comment_organizer() :
    comment_organizer()
  {

    // does nothing.

  }

  three_index_comment_organizer::~three_index_comment_organizer()
  {

    for( std::vector< triplet_indexed_comment* >::iterator
           deletion_iterator = comment_strings.begin();
         comment_strings.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

  }


  bool
  three_index_comment_organizer::has_comment( int const requested_first_index,
                                              int const requested_second_index,
                                              int const requested_third_index,
                                        std::string* const string_for_comment )
  const
  // this looks for the comment for the requested indices, & if it is found,
  // writes it to the provided string, while returning true.
  {

    bool return_bool = false;
    // assume that the requested entry will not be found.

    for( std::vector< triplet_indexed_comment* >::const_iterator
           comment_iterator = comment_strings.begin();
         comment_strings.end() > comment_iterator;
         ++comment_iterator )
      {

        if( ( requested_first_index == (*comment_iterator)->get_first_index() )
            &&
            ( requested_second_index
                == (*comment_iterator)->get_second_index() )
            &&
            ( requested_third_index
                == (*comment_iterator)->get_third_index() ) )
        // if the requested indices are found...
          {

            return_bool = true;
            // note that the associated comment has been found.

            string_for_comment->assign(
                                       *((*comment_iterator)->get_comment()) );
            // record the comment in the given memory address.

            comment_iterator = comment_strings.end();
            // stop looking for the indices now that they has been found.

          }

      }

    return return_bool;
    // return true if requested_first_number was found, false otherwise.

  }

  void
  three_index_comment_organizer::set_comment( int const requested_first_index,
                                              int const requested_second_index,
                                              int const requested_third_index,
                                       std::string const* const given_comment )
  {

    bool comment_needs_to_be_added = true;

    for( std::vector< triplet_indexed_comment* >::iterator
           comment_iterator = comment_strings.begin();
         comment_strings.end() > comment_iterator;
         ++comment_iterator )
      {

        if( ( requested_first_index == (*comment_iterator)->get_first_index() )
            &&
            ( requested_second_index
                == (*comment_iterator)->get_second_index() )
            &&
            ( requested_third_index
                == (*comment_iterator)->get_third_index() ) )
        // if the requested indices are found...
          {

            (*comment_iterator)->set_comment( given_comment );
            // set the requested comment.

            comment_needs_to_be_added = false;
            // note that the comment does not need to be added.

            comment_iterator = comment_strings.end();
            // stop looking for the index now that it has been found.

          }

      }

    if( comment_needs_to_be_added )
    // if the requested indices were not found...
      {

        // create an indexed comment & add it to the set:
        triplet_indexed_comment*
        indexed_comment_adder
        = new triplet_indexed_comment( requested_first_index,
                                       requested_second_index,
                                       requested_third_index,
                                       given_comment );

        comment_strings.push_back( indexed_comment_adder );

      }

  }



  SLHA_BLOCK::SLHA_BLOCK( std::string const given_BLOCK_name,
                          std::string const given_name_comment ) :
    BLOCK_name( given_BLOCK_name ),
    scale_value( CppSLHA_global::really_wrong_value ),
    scale_recorded( false ),
    name_comment_string( given_name_comment ),
    has_recorded_flag( false ),
    has_extra_BLOCK( false ),
    write_as_int_flag( false )
  {

    CppSLHA_global::transform_to_uppercase( &BLOCK_name );
    // this is to ensure all BLOCKs have names in all uppercase.

    // debugging:
    /*std::cout << std::endl << "debugging: SLHA_BLOCK::SLHA_BLOCK( "
    << given_BLOCK_name << ", " << given_name_comment << " ) called.";
    std::cout << std::endl << "this = " << this << ", comment_strings = "
    << comment_strings;
    std::cout << std::endl;
    */

  }

  SLHA_BLOCK::~SLHA_BLOCK()
  {

    if( has_extra_BLOCK )
      {

        delete extra_BLOCK;

      }

  }


  void
  SLHA_BLOCK::make_BLOCK_as_string()
  // this writes the BLOCK in SLHA-format ASCII into BLOCK_as_string.
  {

    std::stringstream return_string_builder;

    return_string_builder << "BLOCK " << BLOCK_name;
    if( scale_recorded )
      {

        return_string_builder << " Q= "
                              << CppSLHA_global::SLHA_double( scale_value );

      }
    return_string_builder << " #";
    if( !( name_comment_string.empty() ) )
      {

        return_string_builder << " " << name_comment_string;

      }
    return_string_builder << std::endl;

    return_string_builder
    << "# nothing recorded from original file for some reason."
    << std::endl;

    BLOCK_as_string.assign( return_string_builder.str() );

  }



  zero_index_BLOCK::zero_index_BLOCK( std::string const given_BLOCK_name,
                                      std::string const given_name_comment ) :
    SLHA_BLOCK( given_BLOCK_name,
                given_name_comment ),
    entry_value( CppSLHA_global::really_wrong_value ),
    entry_recorded( false )
  {

    comment_strings = new one_index_comment_organizer();
    // an entire organizer seems like overkill, but it saves on extra functions
    // for the SLHA_BLOCK derived classes.

  }

  zero_index_BLOCK::~zero_index_BLOCK()
  {

    delete comment_strings;

  }


  void
  zero_index_BLOCK::make_BLOCK_as_string()
  // this writes the BLOCK in SLHA-format ASCII into BLOCK_as_string.
  {

    std::stringstream return_string_builder;

    return_string_builder << "BLOCK " << BLOCK_name;
    if( scale_recorded )
      {

        return_string_builder
        << " Q= " << CppSLHA_global::SLHA_double( scale_value );

      }
    return_string_builder << " #";
    if( !( name_comment_string.empty() ) )
      {

        return_string_builder << " " << name_comment_string;

      }
    return_string_builder << std::endl;

    if( entry_recorded )
      {

        return_string_builder
        << "             "
        << CppSLHA_global::SLHA_int_or_double( entry_value,
                                               write_as_int_flag ) << "   #";
        std::string comment_writer;
        if( comment_strings->has_comment( 1,
                                          &comment_writer ) )
          {

            return_string_builder << " " << comment_writer;

          }
        return_string_builder << std::endl;

      }
    else
      {

        return_string_builder
        << "# was not recorded from original file for some reason."
        << std::endl;

      }

    if( has_extra_BLOCK )
      {

        return_string_builder << *(extra_BLOCK->get_as_string());

      }

    BLOCK_as_string.assign( return_string_builder.str() );

  }

  void
  zero_index_BLOCK::read_input( std::stringstream* const input_line )
  {

    double entry_recorder;

    *input_line >> entry_recorder;

    set_entry( entry_recorder );

  }



  one_index_BLOCK::one_index_BLOCK( std::string const given_BLOCK_name,
                                    std::string const given_name_comment ) :
    SLHA_BLOCK( given_BLOCK_name,
                given_name_comment )
  {

    comment_strings = new one_index_comment_organizer();

    // debugging:
    /*std::cout
     << std::endl
     << "debugging:"
     << " one_index_BLOCK::one_index_BLOCK( "
     << given_BLOCK_name << ", " << given_name_comment << " ) called.";
     std::cout << std::endl << "this = " << this << ", comment_strings = "
     << comment_strings;
     std::cout << std::endl;
     */

  }

  one_index_BLOCK::~one_index_BLOCK()
  {

    delete comment_strings;

  }


  bool
  one_index_BLOCK::has_entry( int const requested_element,
                              double* const double_for_entry )
  const
  // this checks to see if the entry exists in this set & puts its value into
  // the provided memory address for a double.
  {

    if( ( 0 < requested_element )
        &&
        ( (unsigned int)requested_element <= recorded.size() ) )
      {

        *double_for_entry = entries.at( ( requested_element - 1 ) );

        return recorded.at( ( requested_element - 1 ) );

      }
    else
      {

        return false;

      }

  }

  double
  one_index_BLOCK::get_entry( int const requested_element )
  const
  // this returns the entry in non-insane notation (i.e. the first entry is
  // sought with 1, not 0) or returns CppSLHA_global::really_wrong_value.
  {

    if( ( 0 < requested_element )
        &&
        ( (unsigned int)requested_element <= entries.size() ) )
      {

        if( recorded.at( ( requested_element - 1 ) ) )
          {

            return entries.at( ( requested_element - 1 ) );

          }
        else
        // otherwise the entry was not found...
          {

          /** I've decided to turn off these messages.
            std::cout
            << std::endl
            << "CppSLHA::error! an SLHA BLOCK was asked for an entry which it"
            << " had not recorded.";
            std::cout
            << std::endl
            << CppSLHA_global::really_wrong_value_string
            << " has been returned.";
            std::cout << std::endl;**/

            return CppSLHA_global::really_wrong_value;

          }

      }
    else
    // otherwise an entry out of range was requested.
      {

      /** I've decided to turn off these messages.
        std::cout
        << std::endl
        << "CppSLHA::error! an SLHA BLOCK was asked for an entry which it"
        << " had not recorded.";
        std::cout
        << std::endl
        << CppSLHA_global::really_wrong_value_string
        << " has been returned.";
        std::cout << std::endl;**/

        return CppSLHA_global::really_wrong_value;

      }

  }

  void
  one_index_BLOCK::set_entry( int const requested_element,
                              double const input_entry )
  /* this either overwrites the value at entry requested_element, or adds it
   * to the std::vector if it didn't already have it, possibly filling in
   * CppSLHA_global::really_wrong_value entries along the way.
   */
  {

    BLOCK_as_string.clear();  // the string has to be redone now.

    if( ( 0 < requested_element )
        &&
        ( (unsigned int)requested_element <= entries.size() ) )
      // if the entry exists...
      {

        entries.at( ( requested_element - 1 ) ) = input_entry;
        // overwrite the entry.
        recorded.at( ( requested_element - 1 ) ) = true;
        // note that a value for this entry has been recorded (so if it's
        // CppSLHA_global::really_wrong_value, then it's a coincidence).

      }
    else if( 0 < requested_element )
      // otherwise...
      {

        while( (unsigned int)( requested_element - 1 ) > entries.size() )
          // fill out intermediate entries with
          // CppSLHA_global::really_wrong_value.
          {

            entries.push_back( CppSLHA_global::really_wrong_value );
            recorded.push_back( false );
            // note that a value for this entry has been NOT recorded.

          }

        entries.push_back( input_entry );
        // add the entry.
        recorded.push_back( true );
        // note that a value for this entry has been recorded (so if it's
        // CppSLHA_global::really_wrong_value, then it's a coincidence).

      }
    else
      {

      /** I've decided to turn off these messages.
        std::cout
        << std::endl
        << "CppSLHA::error! this SLHA file has entries with negative indices"
        << " for " << BLOCK_name;
        std::cout
        << std::endl
        << " which it SHOULD NOT have! calling exit( EXIT_FAILURE ).";
        std::cout << std::endl;**/

        exit( EXIT_FAILURE );

      }

  }

  void
  one_index_BLOCK::make_BLOCK_as_string()
  {

    // debugging:
    /*std::cout
    << std::endl
    << "debugging: one_index_BLOCK::make_BLOCK_as_string called"
    << " with pointer " << this;
    std::cout << std::endl;
    */

    std::stringstream return_string_builder;

    int entries_written = 0;

    return_string_builder << "BLOCK " << BLOCK_name;
    if( scale_recorded )
      {

        return_string_builder
        << " Q= " << CppSLHA_global::SLHA_double( scale_value );

      }
    return_string_builder << " #";
    if( !( name_comment_string.empty() ) )
      {

        return_string_builder << " " << name_comment_string;

      }
    return_string_builder << std::endl;

    for( unsigned int entry_counter = 1;
         entries.size() >= entry_counter;
         ++entry_counter )
      {

        double entry_writer;

        if( has_entry( entry_counter,
                       &entry_writer ) )
          {

            ++entries_written;

            return_string_builder
            << " " << CppSLHA_global::SLHA_int( entry_counter,
                                                5 )
            << "   "
            << CppSLHA_global::SLHA_int_or_double( entry_writer,
                                                   write_as_int_flag )
            << "   #";
            std::string comment_writer;
            if( comment_strings->has_comment( entry_counter,
                                              &comment_writer ) )
              {

                return_string_builder << " " << comment_writer;

              }
            return_string_builder << std::endl;

          } // end of if there is an entry.

      } // end of for loop over entries.

    if( 1 > entries_written )
      {

        return_string_builder
        << "# nothing recorded from original file for some reason."
        << std::endl;

      }

    if( has_extra_BLOCK )
      {

        return_string_builder << *(extra_BLOCK->get_as_string());

      }

    BLOCK_as_string.assign( return_string_builder.str() );

  }

  void
  one_index_BLOCK::read_input( std::stringstream* const input_line )
  {

    int index_one;
    double entry_recorder;

    *input_line >> index_one
                >> entry_recorder;

    set_entry( index_one,
               entry_recorder );

  }



  sparse_one_index_BLOCK::sparse_one_index_BLOCK(
                                            std::string const given_BLOCK_name,
                                       std::string const given_name_comment ) :
    SLHA_BLOCK( given_BLOCK_name,
                given_name_comment )
  {

    comment_strings = new one_index_comment_organizer();

    // debugging:
    /*std::cout
     << std::endl
     << "debugging:"
     << " sparse_one_index_BLOCK::"
     << "sparse_one_index_BLOCK( "
     << given_BLOCK_name << ", " <<  given_name_comment << " ) called";
     std::cout << std::endl << "comment_strings pointer = " << comment_strings;
     std::cout << std::endl;
     */

  }

  sparse_one_index_BLOCK::~sparse_one_index_BLOCK()
  {

    delete comment_strings;

  }


  void
  sparse_one_index_BLOCK::make_BLOCK_as_string()
  {

    std::stringstream return_string_builder;

    int entries_written = 0;

    return_string_builder << "BLOCK " << BLOCK_name;
    if( scale_recorded )
      {

        return_string_builder
        << " Q= " << CppSLHA_global::SLHA_double( scale_value );

      }
    return_string_builder << " #";
    if( !( name_comment_string.empty() ) )
      {

        return_string_builder << " " << name_comment_string;

      }
    return_string_builder << std::endl;

    for( int entry_counter = entries.get_lowest_index();
         entries.get_highest_index() >= entry_counter;
         ++entry_counter )
      {

        double entry_writer;

        if( has_entry( entry_counter,
                       &entry_writer ) )
          {

            ++entries_written;

            return_string_builder
            << " " << CppSLHA_global::SLHA_int( entry_counter,
                                                5 )
            << "   "
            << CppSLHA_global::SLHA_int_or_double( entry_writer,
                                                   write_as_int_flag )
            << "   #";
            std::string comment_writer;
            if( comment_strings->has_comment( entry_counter,
                                              &comment_writer ) )
              {

                return_string_builder << " " << comment_writer;

              }
            return_string_builder << std::endl;

          } // end of if there is an entry.

      } // end of for loop over entries.

    if( 1 > entries_written )
      {

        return_string_builder
        << "# nothing recorded from original file for some reason."
        << std::endl;

      }

    if( has_extra_BLOCK )
      {

        return_string_builder << *(extra_BLOCK->get_as_string());

      }

    BLOCK_as_string.assign( return_string_builder.str() );

  }

  void
  sparse_one_index_BLOCK::read_input( std::stringstream* const input_line )
  {

    int index_one;
    double entry_recorder;

    *input_line >> index_one
                >> entry_recorder;

    set_entry( index_one,
               entry_recorder );

  }



  sparse_one_index_string_BLOCK::sparse_one_index_string_BLOCK(
                                            std::string const given_BLOCK_name,
                                       std::string const given_name_comment ) :
    SLHA_BLOCK( given_BLOCK_name,
                given_name_comment )
  {

    entries = new one_index_comment_organizer();
    comment_strings = new one_index_comment_organizer();

  }

  sparse_one_index_string_BLOCK::~sparse_one_index_string_BLOCK()
  {

    delete entries;
    entries = NULL;  // just to be safe.
    delete comment_strings;
    comment_strings = NULL;  // just to be safe.

  }


  void
  sparse_one_index_string_BLOCK::make_BLOCK_as_string()
  {

    std::stringstream return_string_builder;

    int entries_written = 0;

    return_string_builder << "BLOCK " << BLOCK_name;
    if( scale_recorded )
      {

        return_string_builder
        << " Q= " << CppSLHA_global::SLHA_double( scale_value );

      }
    return_string_builder << " #";
    if( !( name_comment_string.empty() ) )
      {

        return_string_builder << " " << name_comment_string;

      }
    return_string_builder << std::endl;

    for( int entry_counter = indices.get_lowest_index();
         indices.get_highest_index() >= entry_counter;
         ++entry_counter )
      {

        std::string entry_writer;

        if( has_string( entry_counter,
                        &entry_writer ) )
          {

            ++entries_written;

            return_string_builder
            << " " << CppSLHA_global::SLHA_int( entry_counter,
                                                5 )
            << "   " << entry_writer << "   #";
            std::string comment_writer;
            if( comment_strings->has_comment( entry_counter,
                                              &comment_writer ) )
              {

                return_string_builder << " " << comment_writer;

              }
            return_string_builder << std::endl;

          } // end of if there is an entry.

      } // end of for loop over entries.

    if( 1 > entries_written )
      {

        return_string_builder
        << "# nothing recorded from original file for some reason."
        << std::endl;

      }

    if( has_extra_BLOCK )
      {

        return_string_builder << *(extra_BLOCK->get_as_string());

      }

    BLOCK_as_string.assign( return_string_builder.str() );

  }

  void
  sparse_one_index_string_BLOCK::read_input(
                                          std::stringstream* const input_line )
  {

    int index_one;

    // find the index:
    *input_line >> index_one;


    // now strip the index, any whitespace & any comment away from the string:

    std::string input_line_as_string( input_line->str() );
    unsigned int entry_start;
    unsigned int entry_end;
    unsigned int comment_start;
    unsigned int comment_end;

    // start at the beginning of the string:
    unsigned int position_of_character_in_line = 0;
    char character_being_read
    = input_line_as_string[ position_of_character_in_line ];

    // read along until the 1st non-space character is found or the end of the
    // line is reached:
    while( ( ' ' == character_being_read )
           &&
           ( input_line_as_string.size() > position_of_character_in_line ) )
      {

        ++position_of_character_in_line;

        character_being_read
        = input_line_as_string[ position_of_character_in_line ];

      }

    if( input_line_as_string.size() > position_of_character_in_line )
      // if we didn't run out of line yet...
      {

        // read along until the 1st non-digit character is found or the end of
        // the line is reached:
        while( ( ( '0' == character_being_read )
                 ||
                 ( '1' == character_being_read )
                 ||
                 ( '2' == character_being_read )
                 ||
                 ( '3' == character_being_read )
                 ||
                 ( '4' == character_being_read )
                 ||
                 ( '5' == character_being_read )
                 ||
                 ( '6' == character_being_read )
                 ||
                 ( '7' == character_being_read )
                 ||
                 ( '8' == character_being_read )
                 ||
                 ( '9' == character_being_read ) )
               &&
               ( input_line_as_string.size()
                   > position_of_character_in_line ) )
          {

            ++position_of_character_in_line;

            character_being_read
            = input_line_as_string[ position_of_character_in_line ];

          }

      }

    if( input_line_as_string.size() > position_of_character_in_line )
      // if we didn't run out of line yet...
      {

        // read along until the 1st non-space character is found or the end of
        // the line is reached:
        while( ( ' ' == character_being_read )
               &&
               ( input_line_as_string.size()
                   > position_of_character_in_line ) )
          {

            ++position_of_character_in_line;

            character_being_read
            = input_line_as_string[ position_of_character_in_line ];

          }

      }
    // at this point we should have read past initial whitespace, the index,
    // then whitespace up to the entry string.
    entry_start = position_of_character_in_line;


    // read along until the comment character is found or the end of the line
    // is reached:
    while( ( '#' != character_being_read )
           &&
           ( input_line_as_string.size() > position_of_character_in_line ) )
      {

        ++position_of_character_in_line;

        character_being_read
        = input_line_as_string[ position_of_character_in_line ];

      }
    if( '#' == character_being_read )
    // if we did find a comment character...
      {

        entry_end = ( position_of_character_in_line - 1 );
        // note that the entry ends before the comment.

        // strip away the whitespace before & after the comment:
        ++position_of_character_in_line;
        while( ( ' ' == character_being_read )
               &&
               ( input_line_as_string.size()
                   > position_of_character_in_line ) )
          {

            ++position_of_character_in_line;

            character_being_read
            = input_line_as_string[ position_of_character_in_line ];

          }
        // at this point we should have read past initial whitespace.
        comment_start = position_of_character_in_line;

        // read back through any whitespace at the end of the line:
        position_of_character_in_line = ( input_line_as_string.size() - 1 );
        while( ( ' ' == character_being_read )
               &&
               ( comment_start < position_of_character_in_line ) )
          {

            --position_of_character_in_line;

            character_being_read
            = input_line_as_string[ position_of_character_in_line ];

          }
        // after this while loop, we have stripped away any whitespace at the
        // end of the line, so now this is where the comment ends.
        comment_end = position_of_character_in_line;

        if( comment_start < comment_end )
          // if the comment exists...
          {

            std::string comment_string( input_line_as_string.substr(
                                                                 comment_start,
                                           ( comment_end - comment_start ) ) );
            comment_strings->set_comment( index_one,
                                          &comment_string );
            // record the comment.

          }

      }  // end of if # character was found.
    else
      {

        // read back through any whitespace at the end of the line:
        position_of_character_in_line = ( input_line_as_string.size() - 1 );
        while( ( ' ' == character_being_read )
               &&
               ( entry_start < position_of_character_in_line ) )
          {

            --position_of_character_in_line;

            character_being_read
            = input_line_as_string[ position_of_character_in_line ];

          }
        entry_end = position_of_character_in_line;

      }
    set_entry( index_one,
               input_line_as_string.substr( entry_start,
                                            ( entry_end - entry_start ) ) );

  }



  multiple_same_one_index_string_BLOCK::multiple_same_one_index_string_BLOCK(
                                            std::string const given_BLOCK_name,
                                       std::string const given_name_comment ) :
    SLHA_BLOCK( given_BLOCK_name,
                given_name_comment ),
    lowest_index( 0 ),
    highest_index( 0 ),
    largest_hidden_index( 0 )
  {

    entries = new two_index_comment_organizer();
    comment_strings = new two_index_comment_organizer();

  }

  multiple_same_one_index_string_BLOCK::~multiple_same_one_index_string_BLOCK()
  {

    delete entries;
    delete comment_strings;

  }


  void
  multiple_same_one_index_string_BLOCK::make_BLOCK_as_string()
  {

    std::stringstream return_string_builder;

    int entries_written = 0;

    return_string_builder << "BLOCK " << BLOCK_name;
    if( scale_recorded )
      {

        return_string_builder
        << " Q= " << CppSLHA_global::SLHA_double( scale_value );

      }
    return_string_builder << " #";
    if( !( name_comment_string.empty() ) )
      {

        return_string_builder << " " << name_comment_string;

      }
    return_string_builder << std::endl;

    for( int entry_counter = lowest_index;
         highest_index >= entry_counter;
         ++entry_counter )
      {

        for( int hidden_counter = 1;
             largest_hidden_index >= hidden_counter;
             ++hidden_counter )
          {

            std::string entry_writer;

            if( has_string( entry_counter,
                            hidden_counter,
                            &entry_writer ) )
              {

                ++entries_written;

                return_string_builder
                << " " << CppSLHA_global::SLHA_int( entry_counter,
                                                    5 )
                << "   " << entry_writer << "   #";
                std::string comment_writer;
                if( comment_strings->has_comment( entry_counter,
                                                  hidden_counter,
                                                  &comment_writer ) )
                  {

                    return_string_builder << " " << comment_writer;

                  }
                return_string_builder << std::endl;

              } // end of if there is an entry.

          } // end of loop over hidden entries.

      } // end of for loop over shown entries.

    if( 1 > entries_written )
      {

        return_string_builder
        << "# nothing recorded from original file for some reason."
        << std::endl;

      }

    if( has_extra_BLOCK )
      {

        return_string_builder << *(extra_BLOCK->get_as_string());

      }

    BLOCK_as_string.assign( return_string_builder.str() );

  }

  void
  multiple_same_one_index_string_BLOCK::read_input(
                                          std::stringstream* const input_line )
  {

    // debugging:
    /*std::cout
    << std::endl
    << "debugging: multiple_same_one_index_string_BLOCK::read_input()"
    << " called with:";
    std::cout << std::endl << input_line->str();
    std::cout << std::endl;
    */


    int index_one;

    // find the index:
    *input_line >> index_one;

    int hidden_index = 1;
    std::string dummy_string;
    while( has_string( index_one, hidden_index, &dummy_string ) )
      {

        ++hidden_index;

      }
    // after the while loop, hidden index should be 1 larger than the largest
    // recorded hidden index for that shown index.

    // update the index endpoints:
    if( largest_hidden_index < hidden_index )
      {

        largest_hidden_index = hidden_index;

      }
    if( lowest_index > index_one )
      {

        lowest_index = index_one;

      }
    if( highest_index < index_one )
      {

        highest_index = index_one;

      }


    // now strip the index, any whitespace & any comment away from the string:

    std::string input_line_as_string( input_line->str() );
    unsigned int entry_start;
    unsigned int entry_end;
    unsigned int comment_start;
    unsigned int comment_end;

    // start at the beginning of the string:
    unsigned int position_of_character_in_line = 0;
    char character_being_read
    = input_line_as_string[ position_of_character_in_line ];

    // read along until the 1st non-space character is found or the end of the
    // line is reached:
    while( ( ' ' == character_being_read )
           &&
           ( input_line_as_string.size() > position_of_character_in_line ) )
      {

        ++position_of_character_in_line;

        character_being_read
        = input_line_as_string[ position_of_character_in_line ];

      }

    if( input_line_as_string.size() > position_of_character_in_line )
      // if we didn't run out of line yet...
      {

        // read along until the 1st non-digit character is found or the end of
        // the line is reached:
        while( ( ( '0' == character_being_read )
                 ||
                 ( '1' == character_being_read )
                 ||
                 ( '2' == character_being_read )
                 ||
                 ( '3' == character_being_read )
                 ||
                 ( '4' == character_being_read )
                 ||
                 ( '5' == character_being_read )
                 ||
                 ( '6' == character_being_read )
                 ||
                 ( '7' == character_being_read )
                 ||
                 ( '8' == character_being_read )
                 ||
                 ( '9' == character_being_read ) )
               &&
               ( input_line_as_string.size()
                   > position_of_character_in_line ) )
          {

            ++position_of_character_in_line;

            character_being_read
            = input_line_as_string[ position_of_character_in_line ];

          }

      }

    if( input_line_as_string.size() > position_of_character_in_line )
      // if we didn't run out of line yet...
      {

        // read along until the 1st non-space character is found or the end of
        // the line is reached:
        while( ( ' ' == character_being_read )
               &&
               ( input_line_as_string.size()
                   > position_of_character_in_line ) )
          {

            ++position_of_character_in_line;

            character_being_read
            = input_line_as_string[ position_of_character_in_line ];

          }

      }
    // at this point we should have read past initial whitespace, the index,
    // then whitespace up to the entry string.
    entry_start = position_of_character_in_line;


    // read along until the comment character is found or the end of the line
    // is reached:
    while( ( '#' != character_being_read )
           &&
           ( input_line_as_string.size() > position_of_character_in_line ) )
      {

        ++position_of_character_in_line;

        character_being_read
        = input_line_as_string[ position_of_character_in_line ];

      }
    if( '#' == character_being_read )
    // if we did find a comment character...
      {

        entry_end = ( position_of_character_in_line - 1 );
        // note that the entry ends before the comment.

        // strip away the whitespace before & after the comment:
        ++position_of_character_in_line;
        while( ( ' ' == character_being_read )
               &&
               ( input_line_as_string.size()
                   > position_of_character_in_line ) )
          {

            ++position_of_character_in_line;

            character_being_read
            = input_line_as_string[ position_of_character_in_line ];

          }
        // at this point we should have read past initial whitespace.
        comment_start = position_of_character_in_line;

        // read back through any whitespace at the end of the line:
        position_of_character_in_line = ( input_line_as_string.size() - 1 );
        while( ( ' ' == character_being_read )
               &&
               ( comment_start < position_of_character_in_line ) )
          {

            --position_of_character_in_line;

            character_being_read
            = input_line_as_string[ position_of_character_in_line ];

          }
        // at this point we should have read back past final whitespace.
        comment_end = ( position_of_character_in_line + 1 );

        if( comment_start < comment_end )
          {

            std::string comment_string( input_line_as_string.substr(
                                                                 comment_start,
                                           ( comment_end - comment_start ) ) );
            comment_strings->set_comment( index_one,
                                          hidden_index,
                                          &comment_string );

          }

      }  // end of if # character was found.
    else
      {

        // read back through any whitespace at the end of the line:
        position_of_character_in_line = ( input_line_as_string.size() - 1 );
        while( ( ' ' == character_being_read )
               &&
               ( entry_start < position_of_character_in_line ) )
          {

            --position_of_character_in_line;

            character_being_read
            = input_line_as_string[ position_of_character_in_line ];

          }
        entry_end = position_of_character_in_line;

      }
    set_entry( index_one,
               hidden_index,
               input_line_as_string.substr( entry_start,
                                            ( entry_end - entry_start ) ) );

  }



  two_index_BLOCK::two_index_BLOCK( std::string const given_BLOCK_name,
                                    std::string const given_name_comment ) :
    SLHA_BLOCK( given_BLOCK_name,
                given_name_comment )
  {

    comment_strings = new two_index_comment_organizer();

  }

  two_index_BLOCK::~two_index_BLOCK()
  {

    for( std::vector< std::vector< double >* >::iterator
           deletion_iterator = entries.begin();
         entries.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

    for( std::vector< std::vector< bool >* >::iterator
           deletion_iterator = recorded.begin();
         recorded.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

    delete comment_strings;

  }


  void
  two_index_BLOCK::specific_reset()
  // this drops all held data apart from the name.
  {

    for( std::vector< std::vector< double >* >::iterator
           deletion_iterator = entries.begin();
         entries.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;
        *deletion_iterator = NULL;  // just to be safe.

      }

    for( std::vector< std::vector< bool >* >::iterator
           deletion_iterator = recorded.begin();
         recorded.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;
        *deletion_iterator = NULL;  // just to be safe.

      }

    entries.clear();
    recorded.clear();

    has_recorded_flag = false;

    if( has_extra_BLOCK )
      {

        delete extra_BLOCK;
        extra_BLOCK = NULL;  // just to be safe.
        has_extra_BLOCK = false;

      }

  }


  bool
  two_index_BLOCK::has_entry( int const requested_first_element,
                              int const requested_second_element,
                              double* const double_for_entry )
  const
  {

    if( ( 0 < requested_first_element )
        &&
        ( (unsigned int)requested_first_element <= recorded.size() ) )
      {

        if( ( 0 < requested_second_element )
            &&
            ( (unsigned int)requested_second_element
              <= recorded.at( ( requested_first_element - 1 ) )->size() ) )
          {

            *double_for_entry
            = entries.at( ( requested_first_element - 1 )
                                     )->at( ( requested_second_element - 1 ) );

            return recorded.at( ( requested_first_element - 1 )
                                     )->at( ( requested_second_element - 1 ) );

          }
        else
          {

            return false;

          }

      }
    else
      {

        return false;

      }

  }

  double
  two_index_BLOCK::get_entry( int const requested_first_element,
                              int const requested_second_element )
  const
  // this returns the entry in non-insane notation (i.e. the first entry is
  // sought with 1, not 0) or returns CppSLHA_global::really_wrong_value.
  {

    if( ( 0 < requested_first_element )
        &&
        ( (unsigned int)requested_first_element <= entries.size() ) )
      {

        if( ( 0 < requested_second_element )
            &&
            ( (unsigned int)requested_second_element
              <= entries.at( ( requested_first_element - 1 ) )->size() ) )
          {

            if( recorded.at( ( requested_first_element - 1 )
                                    )->at( ( requested_second_element - 1 ) ) )
              {

                return entries.at( ( requested_first_element - 1 )
                                     )->at( ( requested_second_element - 1 ) );

              }
            else
            // otherwise the entry was not recorded...
              {

              /** I've decided to turn off these messages.
                std::cout
                << std::endl
                << "CppSLHA::error! an SLHA BLOCK was asked for an entry"
                << " which it had not recorded.";
                std::cout
                << std::endl
                << CppSLHA_global::really_wrong_value_string
                << " has been returned.";
                std::cout << std::endl;**/

                return CppSLHA_global::really_wrong_value;

              }

          }
        else
        // otherwise an entry out of the column range was requested.
          {

          /** I've decided to turn off these messages.
            std::cout
            << std::endl
            << "CppSLHA::error! an SLHA BLOCK was asked for an entry which it"
            << " had not recorded.";
            std::cout
            << std::endl
            << CppSLHA_global::really_wrong_value_string
            << " has been returned.";
            std::cout << std::endl;**/

            return CppSLHA_global::really_wrong_value;

          }

      }
    else
    // otherwise the entry out of the row range was requested.
      {

      /** I've decided to turn off these messages.
        std::cout
        << std::endl
        << "CppSLHA::error! an SLHA BLOCK was asked for an entry which it"
        << " had not recorded.";
        std::cout
        << std::endl
        << CppSLHA_global::really_wrong_value_string
        << " has been returned.";
        std::cout << std::endl;**/

        return CppSLHA_global::really_wrong_value;

      }

  }

  void
  two_index_BLOCK::set_entry( int const requested_first_element,
                              int const requested_second_element,
                              double const input_entry )
  /* this either overwrites the value at entry requested_element, or adds it
   * to the std::vector if it didn't already have it, possibly filling in
   * CppSLHA_global::really_wrong_value entries along the way (in the 2nd
   * argument, just empty std::vectors in the 1st argument).
   */
  {

    BLOCK_as_string.clear();  // the string has to be redone now.

    if( ( 0 < requested_first_element )
        &&
        ( (unsigned int)requested_first_element <= entries.size() ) )
      {

        if( ( 0 < requested_second_element )
            &&
            ( (unsigned int)requested_second_element
              <= entries.at( ( requested_first_element - 1 ) )->size() ) )
          // if the entry exists...
          {

            entries.at( ( requested_first_element - 1 )
                                      )->at( ( requested_second_element - 1 ) )
            = input_entry;
            // overwrite the entry with the given value.
            recorded.at( ( requested_first_element - 1 )
                                      )->at( ( requested_second_element - 1 ) )
            = true;
            // note that a value for this entry has been recorded (so if it's
            // CppSLHA_global::really_wrong_value, then it's a coincidence).

          }
        else if( 0 < requested_second_element )
          // otherwise...
          {

            while( (unsigned int)( requested_second_element - 1 )
                     > entries.at( ( requested_first_element - 1 ) )->size() )
              // fill out intermediate entries with
              // CppSLHA_global::really_wrong_value.
              {

                entries.at( ( requested_first_element - 1 )
                            )->push_back( CppSLHA_global::really_wrong_value );
                recorded.at( ( requested_first_element - 1 )
                                                         )->push_back( false );
                // note that a value for this entry has been NOT recorded.

              }

            entries.at( ( requested_first_element - 1 ) )->push_back(
                                                                 input_entry );
            // add the entry.
            recorded.at( ( requested_first_element - 1 ) )->push_back( true );
            // note that a value for this entry has been recorded (so if it's
            // CppSLHA_global::really_wrong_value, then it's a coincidence).

          }
        else
          {

            std::cout
            << std::endl
            << "CppSLHA::error! this SLHA file has entries with negative"
            << " indices for " << BLOCK_name;
            std::cout
            << std::endl
            << " which it SHOULD NOT have! calling exit( EXIT_FAILURE ).";
            std::cout << std::endl;

            exit( EXIT_FAILURE );

          }

      }
    else if( 0 < requested_first_element )
    // otherwise the entry needs a new row or rows.
      {

        while( (unsigned int)( requested_first_element - 1 ) > entries.size() )
        // fill out intermediate entries with NULL vectors.
          {

            std::vector< double >* vector_adder = new std::vector< double >;
            // make a new empty std::vector for the entry.

            entries.push_back( vector_adder );

            std::vector< bool >* empty_record_adder = new std::vector< bool >;
            // make a new empty std::vector for the entry record.

            recorded.push_back( empty_record_adder );

          }

        std::vector< double >* vector_adder = new std::vector< double >;
        // make a new std::vector for the entry.

        std::vector< bool >* record_adder = new std::vector< bool >;
        // make a new empty std::vector for the entry record.

        while( (unsigned int)( requested_second_element - 1 )
               > vector_adder->size() )
          // fill out intermediate entries with
          // CppSLHA_global::really_wrong_value.
          {

            vector_adder->push_back( CppSLHA_global::really_wrong_value );
            record_adder->push_back( false );
            // note that no value for this entry has actually been recorded.

          }

        vector_adder->push_back( input_entry );
        // add the entry to this std::vector.
        record_adder->push_back( true );
        // note that a value for this entry has been recorded (so if it's
        // CppSLHA_global::really_wrong_value, then it's a coincidence).

        // add the pointers to the std::vectors to the std::vectors of
        // pointers:
        entries.push_back( vector_adder );
        recorded.push_back( record_adder );

      }
    else
      {

        std::cout
        << std::endl
        << "CppSLHA::error! this SLHA file has entries with negative indices"
        << " for " << BLOCK_name;
        std::cout
        << std::endl
        << " which it SHOULD NOT have! calling exit( EXIT_FAILURE ).";
        std::cout << std::endl;

        exit( EXIT_FAILURE );

      }

  }


  void
  two_index_BLOCK::make_BLOCK_as_string()
  {

    std::stringstream return_string_builder;

    int entries_written = 0;

    return_string_builder << "BLOCK " << BLOCK_name;
    if( scale_recorded )
      {

        return_string_builder
        << " Q= " << CppSLHA_global::SLHA_double( scale_value );

      }
    return_string_builder << " #";
    if( !( name_comment_string.empty() ) )
      {

        return_string_builder << " " << name_comment_string;

      }
    return_string_builder << std::endl;

    for( unsigned int first_index_counter = 1;
         entries.size() >= first_index_counter;
         ++first_index_counter )
      {

        for( unsigned int second_index_counter = 1;
             entries.at( ( first_index_counter - 1 ) )->size()
             >= second_index_counter;
             ++second_index_counter )
          {

            double entry_writer;

            if( has_entry( first_index_counter,
                           second_index_counter,
                           &entry_writer ) )
              {

                ++entries_written;

                return_string_builder
                << " " << CppSLHA_global::SLHA_int( first_index_counter,
                                                    2 )
                << " " << CppSLHA_global::SLHA_int( second_index_counter,
                                                    2 )
                << "   "
                << CppSLHA_global::SLHA_int_or_double( entry_writer,
                                                       write_as_int_flag )
                << "   #";
                std::string comment_writer;

                if( comment_strings->has_comment( first_index_counter,
                                                  second_index_counter,
                                                  &comment_writer ) )
                  {

                    return_string_builder << " " << comment_writer;

                  }
                return_string_builder << std::endl;

              } // end of if there is an entry.

          } // end of for loop over columns.

      } // end of for loop over rows.

    if( 1 > entries_written )
      {

        return_string_builder
        << "# nothing recorded from original file for some reason."
        << std::endl;

      }

    if( has_extra_BLOCK )
      {

        return_string_builder << *(extra_BLOCK->get_as_string());

      }

    BLOCK_as_string.assign( return_string_builder.str() );

  }

  void
  two_index_BLOCK::read_input( std::stringstream* const input_line )
  {

    int index_one;
    int index_two;
    double entry_recorder;

    *input_line >> index_one
                >> index_two
                >> entry_recorder;

    set_entry( index_one,
               index_two,
               entry_recorder );

  }



  three_index_BLOCK::three_index_BLOCK( std::string const given_BLOCK_name,
                                       std::string const given_name_comment ) :
    SLHA_BLOCK( given_BLOCK_name,
                given_name_comment )
  {

    comment_strings = new three_index_comment_organizer();

  }

  three_index_BLOCK::~three_index_BLOCK()
  {

    for( std::vector< std::vector< std::vector< double >* >* >::iterator
           vector_iterator = entries.begin();
         entries.end() > vector_iterator;
         ++vector_iterator )
      {

        for( std::vector< std::vector< double >* >::iterator
               deletion_iterator = (*vector_iterator)->begin();
             (*vector_iterator)->end() > deletion_iterator;
             ++deletion_iterator )
          {

            delete *deletion_iterator;

          }

        delete *vector_iterator;

      }

    for( std::vector< std::vector< std::vector< bool >* >* >::iterator
           vector_iterator = recorded.begin();
         recorded.end() > vector_iterator;
         ++vector_iterator )
      {

        for( std::vector< std::vector< bool >* >::iterator
               deletion_iterator = (*vector_iterator)->begin();
             (*vector_iterator)->end() > deletion_iterator;
             ++deletion_iterator )
          {

            delete *deletion_iterator;

          }

        delete *vector_iterator;

      }

    delete comment_strings;

  }


  void
  three_index_BLOCK::specific_reset()
  // this drops all held data apart from the name.
  {

    for( std::vector< std::vector< std::vector< double >* >* >::iterator
           vector_iterator = entries.begin();
         entries.end() > vector_iterator;
         ++vector_iterator )
      {

        for( std::vector< std::vector< double >* >::iterator
               deletion_iterator = (*vector_iterator)->begin();
             (*vector_iterator)->end() > deletion_iterator;
             ++deletion_iterator )
          {

            delete *deletion_iterator;
            *deletion_iterator = NULL;  // just to be safe.

          }

      }

    for( std::vector< std::vector< std::vector< bool >* >* >::iterator
           vector_iterator = recorded.begin();
         recorded.end() > vector_iterator;
         ++vector_iterator )
      {

        for( std::vector< std::vector< bool >* >::iterator
               deletion_iterator = (*vector_iterator)->begin();
             (*vector_iterator)->end() > deletion_iterator;
             ++deletion_iterator )
          {

            delete *deletion_iterator;
            *deletion_iterator = NULL;  // just to be safe.

          }

      }

    entries.clear();
    recorded.clear();

    has_recorded_flag = false;

    if( has_extra_BLOCK )
      {

        delete extra_BLOCK;
        extra_BLOCK = NULL;  // just to be safe.
        has_extra_BLOCK = false;

      }

  }


  bool
  three_index_BLOCK::has_entry( int const requested_first_element,
                                int const requested_second_element,
                                int const requested_third_element,
                                double* const double_for_entry )
  const
  {

    if( ( 0 < requested_first_element )
        &&
        ( (unsigned int)requested_first_element <= recorded.size() ) )
      {

        if( ( 0 < requested_second_element )
            &&
            ( (unsigned int)requested_second_element
              <= recorded.at( ( requested_first_element - 1 ) )->size() ) )
          {

            if( ( 0 < requested_third_element )
                &&
                ( (unsigned int)requested_third_element
                  <= recorded.at( ( requested_first_element - 1 )
                          )->at( ( requested_second_element - 1 ) )->size() ) )
              {

                *double_for_entry
                = entries.at( ( requested_first_element - 1 )
                                        )->at( ( requested_second_element - 1 )
                                      )->at( ( requested_third_element - 1 ) );

                return recorded.at( ( requested_first_element - 1 )
                                        )->at( ( requested_second_element - 1 )
                                      )->at( ( requested_third_element - 1 ) );

              }
            else
              {

                return false;

              }

          }
        else
          {

            return false;

          }

      }
    else
      {

        return false;

      }

  }

  double
  three_index_BLOCK::get_entry( int const requested_first_element,
                                int const requested_second_element,
                                int const requested_third_element )
  const
  // this returns the entry in non-insane notation (i.e. the first entry is
  // sought with 1, not 0) or returns CppSLHA_global::really_wrong_value.
  {

    if( ( 0 < requested_first_element )
        &&
        ( (unsigned int)requested_first_element <= entries.size() ) )
      {

        if( ( 0 < requested_second_element )
            &&
            ( (unsigned int)requested_second_element
              <= entries.at( ( requested_first_element - 1 ) )->size() ) )
          {

            if( ( 0 < requested_third_element )
                &&
                ( (unsigned int)requested_third_element
                  <= entries.at( ( requested_first_element - 1 )
                          )->at( ( requested_second_element - 1 ) )->size() ) )
              {

                if( recorded.at( ( requested_first_element - 1 )
                                        )->at( ( requested_second_element - 1 )
                                     )->at( ( requested_third_element - 1 ) ) )
                  {

                    return entries.at( ( requested_first_element - 1 )
                                        )->at( ( requested_second_element - 1 )
                                      )->at( ( requested_third_element - 1 ) );

                  }
                else
                // otherwise the entry was not recorded...
                  {

                  /** I've decided to turn off these messages.
                    std::cout
                    << std::endl
                    << "CppSLHA::error! an SLHA BLOCK was asked for an entry"
                    << " which it had not recorded.";
                    std::cout
                    << std::endl
                    << CppSLHA_global::really_wrong_value_string
                    << " has been returned.";
                    std::cout << std::endl;**/

                    return CppSLHA_global::really_wrong_value;

                  }

              }
            else
            // otherwise the entry out of the
            // thing-after-rows-and-columns-have-been-used range was requested.
              {

              /** I've decided to turn off these messages.
                std::cout
                << std::endl
                << "CppSLHA::error! an SLHA BLOCK was asked for an entry"
                << " which it had not recorded.";
                std::cout
                << std::endl
                << CppSLHA_global::really_wrong_value_string
                << " has been returned.";
                std::cout << std::endl;**/

                return CppSLHA_global::really_wrong_value;

              }

          }
        else
        // otherwise an entry out of the column range was requested.
          {

          /** I've decided to turn off these messages.
            std::cout
            << std::endl
            << "CppSLHA::error! an SLHA BLOCK was asked for an entry which it"
            << " had not recorded.";
            std::cout
            << std::endl
            << CppSLHA_global::really_wrong_value_string
            << " has been returned.";
            std::cout << std::endl;**/

            return CppSLHA_global::really_wrong_value;

          }

      }
    else
    // otherwise the entry out of the row range was requested.
      {

      /** I've decided to turn off these messages.
        std::cout
        << std::endl
        << "CppSLHA::error! an SLHA BLOCK was asked for an entry which it"
        << " had not recorded.";
        std::cout
        << std::endl
        << CppSLHA_global::really_wrong_value_string
        << " has been returned.";
        std::cout << std::endl;**/

        return CppSLHA_global::really_wrong_value;

      }

  }

  void
  three_index_BLOCK::set_entry( int const requested_first_element,
                                int const requested_second_element,
                                int const requested_third_element,
                                double const input_entry )
  /* this either overwrites the value at entry requested_element, or adds it
   * to the std::vector if it didn't already have it, possibly filling in
   * CppSLHA_global::really_wrong_value entries along the way (in the 3rd
   * argument, just empty std::vectors in the 1st & 2nd arguments).
   */
  {

    BLOCK_as_string.clear();  // the string has to be redone now.

    if( ( 0 < requested_first_element )
        &&
        ( (unsigned int)requested_first_element <= entries.size() ) )
      {

        if( ( 0 < requested_second_element )
            &&
            ( (unsigned int)requested_second_element
              <= entries.at( ( requested_first_element - 1 ) )->size() ) )
        // if the entry exists...
          {

            if( ( 0 < requested_third_element )
                &&
                ( (unsigned int)requested_third_element
                  <= entries.at( ( requested_first_element - 1 )
                          )->at( ( requested_second_element - 1 ) )->size() ) )
              // if the entry exists...
              {

                entries.at( ( requested_first_element - 1 )
                                        )->at( ( requested_second_element - 1 )
                                       )->at( ( requested_third_element - 1 ) )
                = input_entry;
                // overwrite the entry with the given value.
                recorded.at( ( requested_first_element - 1 )
                                        )->at( ( requested_second_element - 1 )
                                       )->at( ( requested_third_element - 1 ) )
                = true;
                /* note that a value for this entry has been recorded (so if
                 * it's CppSLHA_global::really_wrong_value, then it's a
                 * coincidence).
                 */

              }
            else if( 0 < requested_third_element )
              // otherwise...
              {

                while( (unsigned int)( requested_third_element - 1 )
                       > entries.at( ( requested_first_element - 1 )
                            )->at( ( requested_second_element - 1 ) )->size() )
                  // fill out intermediate entries with
                  // CppSLHA_global::really_wrong_value.
                  {

                    entries.at( ( requested_first_element - 1 )
                                        )->at( ( requested_second_element - 1 )
                            )->push_back( CppSLHA_global::really_wrong_value );
                    recorded.at( ( requested_first_element - 1 )
                                        )->at( ( requested_second_element - 1 )
                                                         )->push_back( false );
                    // note that a value for this entry has been NOT recorded.

                  }

                entries.at( ( requested_first_element - 1 )
                                        )->at( ( requested_second_element - 1 )
                                                   )->push_back( input_entry );
                // add the entry.
                recorded.at( ( requested_first_element - 1 )
                                        )->at( ( requested_second_element - 1 )
                                                          )->push_back( true );
                /* note that a value for this entry has been recorded (so if
                 * it's CppSLHA_global::really_wrong_value, then it's a
                 * coincidence).
                 */

              }
            else
              {

                std::cout
                << std::endl
                << "CppSLHA::error! this SLHA file has entries with negative"
                << " indices for " << BLOCK_name;
                std::cout
                << std::endl
                << " which it SHOULD NOT have! calling exit( EXIT_FAILURE ).";
                std::cout << std::endl;

                exit( EXIT_FAILURE );

              }

          }
        else if( 0 < requested_second_element )
        // otherwise the entry needs a new column or new columns.
          {

            while( (unsigned int)( requested_second_element - 1 )
                   > entries.at( ( requested_first_element - 1 ) )->size() )
            // fill out intermediate entries with NULL vectors.
              {

                std::vector< double >*
                vector_adder = new std::vector< double >;
                // make a new empty std::vector for the entry.

                entries.at( ( requested_first_element - 1 )
                                                  )->push_back( vector_adder );

                std::vector< bool >*
                empty_record_adder = new std::vector< bool >;
                // make a new empty std::vector for the entry record.

                recorded.at( ( requested_first_element - 1 )
                                            )->push_back( empty_record_adder );

              }

            std::vector< double >* vector_adder = new std::vector< double >;
            // make a new std::vector for the entry.

            std::vector< bool >* record_adder = new std::vector< bool >;
            // make a new empty std::vector for the entry record.

            while( (unsigned int)( requested_third_element - 1 )
                   > vector_adder->size() )
              // fill out intermediate entries with
              // CppSLHA_global::really_wrong_value.
              {

                vector_adder->push_back( CppSLHA_global::really_wrong_value );
                record_adder->push_back( false );
                // note that no value for this entry has actually been
                // recorded.

              }

            vector_adder->push_back( input_entry );
            // add the entry to this std::vector.
            record_adder->push_back( true );
            // note that a value for this entry has been recorded (so if it's
            // CppSLHA_global::really_wrong_value, then it's a coincidence).

            // add the pointers to the std::vectors to the std::vectors of
            // pointers:
            entries.at( ( requested_first_element - 1 )
                                                  )->push_back( vector_adder );
            recorded.at( ( requested_first_element - 1 )
                                                  )->push_back( record_adder );
          }
        else
          {

            std::cout
            << std::endl
            << "CppSLHA::error! this SLHA file has entries with negative"
            << " indices for " << BLOCK_name;
            std::cout
            << std::endl
            << " which it SHOULD NOT have! calling exit( EXIT_FAILURE ).";
            std::cout << std::endl;

            exit( EXIT_FAILURE );

          }

      }
    else if( 0 < requested_first_element )
    // otherwise the entry needs a new row or new rows.
      {

        std::vector< std::vector< double >* >* requested_double_column;
        std::vector< std::vector< bool >* >* requested_bool_column;

        while( (unsigned int)requested_first_element > entries.size() )
        /* fill out intermediate entries with NULL vectors as well as the
         * vector of the requested index (hence not
         * "( requested_first_element - 1 )").
         */
          {

            requested_double_column
            = new std::vector< std::vector< double >* >;
            // make a new empty std::vector for the entry.

            entries.push_back( requested_double_column );

            requested_bool_column = new std::vector< std::vector< bool >* >;
            // make a new empty std::vector for the entry record.

            recorded.push_back( requested_bool_column );

          }
        while( (unsigned int)( requested_second_element - 1 )
               > requested_double_column->size() )
        // fill out intermediate entries with NULL vectors.
          {

            std::vector< double >* vector_adder = new std::vector< double >;
            // make a new empty std::vector for the entry.

            requested_double_column->push_back( vector_adder );

            std::vector< bool >* empty_record_adder = new std::vector< bool >;
            // make a new empty std::vector for the entry record.

            requested_bool_column->push_back( empty_record_adder );

          }

        std::vector< double >* vector_adder = new std::vector< double >;
        // make a new std::vector for the entry.

        std::vector< bool >* record_adder = new std::vector< bool >;
        // make a new empty std::vector for the entry record.

        while( (unsigned int)( requested_third_element - 1 )
               > vector_adder->size() )
          // fill out intermediate entries with
          // CppSLHA_global::really_wrong_value.
          {

            vector_adder->push_back( CppSLHA_global::really_wrong_value );
            record_adder->push_back( false );
            // note that no value for this entry has actually been recorded.

          }

        vector_adder->push_back( input_entry );
        // add the entry to this std::vector.
        record_adder->push_back( true );
        // note that a value for this entry has been recorded (so if it's
        // CppSLHA_global::really_wrong_value, then it's a coincidence).

        // add the pointers to the std::vectors to the std::vectors of
        // pointers:
        requested_double_column->push_back( vector_adder );
        requested_bool_column->push_back( record_adder );

      }
    else
      {

        std::cout
        << std::endl
        << "CppSLHA::error! this SLHA file has entries with negative indices"
        << " for " << BLOCK_name;
        std::cout
        << std::endl
        << " which it SHOULD NOT have! calling exit( EXIT_FAILURE ).";
        std::cout << std::endl;

        exit( EXIT_FAILURE );

      }

  }


  void
  three_index_BLOCK::make_BLOCK_as_string()
  {

    std::stringstream return_string_builder;

    int entries_written = 0;

    return_string_builder << "BLOCK " << BLOCK_name;
    if( scale_recorded )
      {

        return_string_builder
        << " Q= " << CppSLHA_global::SLHA_double( scale_value );

      }
    return_string_builder << " #";
    if( !( name_comment_string.empty() ) )
      {

        return_string_builder << " " << name_comment_string;

      }
    return_string_builder << std::endl;

    for( unsigned int first_index_counter = 1;
         entries.size() >= first_index_counter;
         ++first_index_counter )
      {

        for( unsigned int second_index_counter = 1;
             entries.at( ( first_index_counter - 1 ) )->size()
                 >= second_index_counter;
             ++second_index_counter )
          {

            for( unsigned int third_index_counter = 1;
                 entries.at( ( first_index_counter - 1 )
                             )->at( ( second_index_counter - 1 ) )->size()
                   >= third_index_counter;
                 ++third_index_counter )
              {

                double entry_writer;

                if( has_entry( first_index_counter,
                               second_index_counter,
                               third_index_counter,
                               &entry_writer ) )
                  {

                    ++entries_written;

                    return_string_builder
                    << " " << CppSLHA_global::SLHA_int( first_index_counter,
                                                        2 )
                    << " " << CppSLHA_global::SLHA_int( second_index_counter,
                                                        2 )
                    << " " << CppSLHA_global::SLHA_int( third_index_counter,
                                                        2 )
                    << "   "
                    << CppSLHA_global::SLHA_int_or_double( entry_writer,
                                                           write_as_int_flag )
                    << "   #";
                    std::string comment_writer;
                    if( comment_strings->has_comment( first_index_counter,
                                                      second_index_counter,
                                                      third_index_counter,
                                                      &comment_writer ) )
                      {

                        return_string_builder << " " << comment_writer;

                      }
                    return_string_builder << std::endl;

                  } // end of if there is an entry.

              }
            // end of loop over thing-after-rows-and-columns-have-been-used.

          } // end of for loop over columns.

      } // end of for loop over rows.

    if( 1 > entries_written )
      {

        return_string_builder
        << "# nothing recorded from original file for some reason."
        << std::endl;

      }

    if( has_extra_BLOCK )
      {

        return_string_builder << *(extra_BLOCK->get_as_string());

      }

    BLOCK_as_string.assign( return_string_builder.str() );

  }

  void
  three_index_BLOCK::read_input( std::stringstream* const input_line )
  {

    int index_one;
    int index_two;
    int index_three;
    double entry_recorder;

    *input_line >> index_one
                >> index_two
                >> index_three
                >> entry_recorder;

    set_entry( index_one,
               index_two,
               index_three,
               entry_recorder );

  }



  MASS_BLOCK::MASS_BLOCK( EW_scale_spectrum* given_spectrum ) :
    SLHA_BLOCK( "MASS",
                "particle masses" ),
    reference_spectrum( given_spectrum )
  {

    // just an initialization list.

  }

  MASS_BLOCK::~MASS_BLOCK()
  {

    // does nothing.

  }

  void
  MASS_BLOCK::make_BLOCK_as_string()
  {

    std::stringstream return_string_builder;

    int entries_written = 0;

    return_string_builder << "BLOCK MASS   #" << std::endl;

    for( std::vector< particle_property_set* >::const_iterator
         particle_set_iterator
         = reference_spectrum->get_particle_property_sets()->begin();
         reference_spectrum->get_particle_property_sets()->end()
         > particle_set_iterator;
         ++particle_set_iterator )
    // for each stored particle property set...
      {

        if( (*particle_set_iterator)->has_mass_been_recorded() )
          {

            ++entries_written;

            return_string_builder
            << " "
            << CppSLHA_global::SLHA_int(
                                      (*particle_set_iterator)->get_PDG_code(),
                                         9 )
            << "   "
            << CppSLHA_global::SLHA_double(
                                         (*particle_set_iterator)->get_mass() )
            << "   # " << *((*particle_set_iterator)->get_name()) << std::endl;

          }

      } // end of for loop over entries.

    if( 1 > entries_written )
      {

        return_string_builder
            << "# nothing recorded from original file for some reason."
            << std::endl;

      }

    BLOCK_as_string.assign( return_string_builder.str() );

  }


  void
  MASS_BLOCK::read_input( std::stringstream* const input_line )
  {

    int index_one;
    double entry_recorder;

    *input_line >> index_one
                >> entry_recorder;

    reference_spectrum->record_PDG_coded_mass( index_one,
                                               entry_recorder );

  }



  SPheno_input_BLOCK::SPheno_input_BLOCK( std::string const given_BLOCK_name,
                                       std::string const given_name_comment ) :
    SLHA_BLOCK( given_BLOCK_name,
                given_name_comment ),
    lowest_index( 0 ),
    highest_index( 0 ),
    largest_hidden_index( 0 )
  {

    comment_strings = new two_index_comment_organizer();

  }

  SPheno_input_BLOCK::~SPheno_input_BLOCK()
  {

    delete comment_strings;

  }


  bool
  SPheno_input_BLOCK::has_entry( int const requested_element,
                                 double* const double_for_entry )
  const
  {

    std::pair< int, int > requested_pair( requested_element,
                                          0 );
    std::map< std::pair< int, int >, double >::const_iterator
    entry_iterator = entries.find( requested_pair );
    // map::find() returns the end of the container if the element is not
    // found.
    if( entries.end() != entry_iterator )
      {

        *double_for_entry = entry_iterator->second;

        return true;

      }
    else
      {

        return false;

      }

  }

  bool
  SPheno_input_BLOCK::has_entry( int const requested_first_element,
                                 int const requested_second_element,
                                 double* const double_for_entry )
  const
  {

    std::pair< int, int > requested_pair( requested_first_element,
                                          requested_second_element );
    std::map< std::pair< int, int >, double >::const_iterator
    entry_iterator = entries.find( requested_pair );
    // map::find() returns the end of the container if the element is not
    // found.
    if( entries.end() != entry_iterator )
      {

        *double_for_entry = entry_iterator->second;

        return true;

      }
    else
      {

        return false;

      }

  }

  void
  SPheno_input_BLOCK::set_entry( int const requested_element,
                                 double const input_entry )
  // this either overwrites the value at entry ( requested_element, 1 ),
  // or adds it to the set if it didn't have it.
  {

    BLOCK_as_string.clear();  // the string has to be redone now.

    std::pair< int, int > requested_pair( requested_element,
                                          1 );
    entries[ requested_pair ] = input_entry;
    /* map::operator[] is great. if the entry exists, it gives a reference
     * to the entry, otherwise it pushes back in a new entry with the
     * sought key.
     */

    // update the index endpoints:
    if( largest_hidden_index < 1 )
      {

        largest_hidden_index = 1;

      }
    if( lowest_index > requested_element )
      {

        lowest_index = requested_element;

      }
    if( highest_index < requested_element )
      {

        highest_index = requested_element;

      }

  }

  void
  SPheno_input_BLOCK::set_entry( int const requested_first_element,
                                 int const requested_second_element,
                                 double const input_entry )
  /* this either overwrites the value at entry
   * ( requested_first_element, requested_second_element ), or adds it to
   * the set if it didn't have it.
   */
  {

    BLOCK_as_string.clear();  // the string has to be redone now.

    std::pair< int, int > requested_pair( requested_first_element,
                                          requested_second_element );
    entries[ requested_pair ] = input_entry;
    /* map::operator[] is great. if the entry exists, it gives a reference
     * to the entry, otherwise it pushes back in a new entry with the
     * sought key.
     */

    // update the index endpoints:
    if( largest_hidden_index < requested_second_element )
      {

        largest_hidden_index = requested_second_element;

      }
    if( lowest_index > requested_first_element )
      {

        lowest_index = requested_first_element;

      }
    if( highest_index < requested_first_element )
      {

        highest_index = requested_first_element;

      }

  }


  void
  SPheno_input_BLOCK::make_BLOCK_as_string()
  {

    // debugging:
    /*std::cout
    << std::endl
    << "debugging! SPheno_input_BLOCK::make_BLOCK_as_string() called.";
    std::cout << std::endl;*/

    std::stringstream return_string_builder;

    int entries_written = 0;

    return_string_builder << "BLOCK " << BLOCK_name;
    if( scale_recorded )
      {

        return_string_builder
        << " Q= " << CppSLHA_global::SLHA_double( scale_value );

      }
    return_string_builder << " #";
    if( !( name_comment_string.empty() ) )
      {

        return_string_builder << " " << name_comment_string;

      }
    return_string_builder << std::endl;

    for( int entry_counter = lowest_index;
         highest_index >= entry_counter;
         ++entry_counter )
      {

        // debugging:
        /*std::cout << std::endl << "checking visible index " << entry_counter;
        std::cout << std::endl;*/

        double entry_writer;

        if( 22 == entry_counter )
          // if we find the set of repeatable acceptanceValues...
          {

            entry_counter = 26;
            /* note that we're writing out entries 22 to 26 as soon as we find
             * the 1st entry 22. (at the end of this statement, the for-loop
             * ++ incrementation should kick in.)
             */

            for( int hidden_counter = 1;
                 largest_hidden_index >= hidden_counter;
                 ++hidden_counter )
              // go through all repetitions of entries 22 to 26.
              {

                // debugging:
                /*std::cout << std::endl << "checking hidden index "
                                       << hidden_counter;
                std::cout << std::endl;*/

                double entry_writer;
                for( int repetition_counter = 22;
                     26 >= repetition_counter;
                     ++repetition_counter )
                  {

                    if( has_entry( repetition_counter,
                                   hidden_counter,
                                   &entry_writer ) )
                      {

                        ++entries_written;

                        return_string_builder
                        << " "
                        << CppSLHA_global::SLHA_int( repetition_counter,
                                                     5 )
                        << "   "
                        << CppSLHA_global::SLHA_int_or_double( entry_writer,
                                                            write_as_int_flag )
                        << "   #";
                        std::string comment_writer;
                        if( comment_strings->has_comment( repetition_counter,
                                                          hidden_counter,
                                                          &comment_writer ) )
                          {

                            return_string_builder << " " << comment_writer;

                          }
                          return_string_builder << std::endl;

                      }  // end of if statement constructing line.

                  }
                // end of for loop over repeatable entries of given hidden
                // index.

              } // end of looping over hidden indices for repeatable entries.

          }  // end of if( 22 == entry_counter )
        else
          {

            if( has_entry( entry_counter,
                           1,
                           &entry_writer ) )
              {

                // debugging:
                /*std::cout << std::endl << "has entry " << entry_counter;
                std::cout << std::endl;*/

                ++entries_written;

                return_string_builder
                << " "
                << CppSLHA_global::SLHA_int( entry_counter,
                                             5 )
                << "   "
                << CppSLHA_global::SLHA_int_or_double( entry_writer,
                                                       write_as_int_flag )
                << "   #";
                std::string comment_writer;
                if( comment_strings->has_comment( entry_counter,
                                                  1,
                                                  &comment_writer ) )
                  {

                    return_string_builder << " " << comment_writer;

                  }
                  return_string_builder << std::endl;

              } // end of if there is an entry.

          } // end of loop over hidden entries.

      } // end of for loop over shown entries.

    if( 1 > entries_written )
      {

        return_string_builder
        << "# nothing recorded from original file for some reason."
        << std::endl;

      }

    if( has_extra_BLOCK )
      {

        return_string_builder << *(extra_BLOCK->get_as_string());

      }


    // debugging:
    /*std::cout << std::endl << "about to return";
    std::cout << std::endl << return_string_builder.str();
    std::cout << std::endl;*/

    BLOCK_as_string.assign( return_string_builder.str() );

  }


  void
  SPheno_input_BLOCK::read_input( std::stringstream* const input_line )
  {

    int index_one;
    double entry_value;

    // find the index:
    *input_line >> index_one >> entry_value;

    int hidden_index = 1;
    double dummy_double;
    while( has_entry( index_one, hidden_index, &dummy_double ) )
      {

        ++hidden_index;

      }
    // after the while loop, hidden index should be 1 larger than the largest
    // recorded hidden index for that shown index.

    // update the index endpoints:
    if( largest_hidden_index < hidden_index )
      {

        largest_hidden_index = hidden_index;

      }
    if( lowest_index > index_one )
      {

        lowest_index = index_one;

      }
    if( highest_index < index_one )
      {

        highest_index = index_one;

      }

    set_entry( index_one,
               hidden_index,
               entry_value );

  }



  SPheno_cross_section_line::SPheno_cross_section_line(
                                         std::stringstream* const given_input )
  {

    *given_input >> cross_section
                 >> number_of_final_state_particles
                 >> first_particle_code
                 >> second_particle_code;

    if( 3 == number_of_final_state_particles )
      {

        *given_input >> third_particle_code;

      }
    else
      {

        third_particle_code = (int)CppSLHA_global::really_wrong_value;

      }

  }

  SPheno_cross_section_line::~SPheno_cross_section_line()
  {

    // does nothing.

  }


  SPheno_cross_section_set::SPheno_cross_section_set(
                                 EW_scale_spectrum const* const given_spectrum,
                                      std::stringstream* const given_header ) :
    reference_spectrum( given_spectrum )
  {

    std::string dummy_string;

    *given_header >> dummy_string
                  // ignore the initial XS.
                  >> first_incoming_particle_code
                  >> second_incoming_particle_code
                  >> center_of_momentum_energy
                  >> first_incoming_particle_polarization
                  >> second_incoming_particle_polarization
                  >> include_ISR_flag;

  }

  SPheno_cross_section_set::~SPheno_cross_section_set()
  {

    for( std::vector< SPheno_cross_section_line* >::iterator
           deletion_iterator = cross_sections.begin();
         cross_sections.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

  }


  double
  SPheno_cross_section_set::get_cross_section_for_PDG_codes(
                                                      int const first_PDG_code,
                                                    int const second_PDG_code )
  const
  {

    double return_value = CppSLHA_global::really_wrong_value;

    for( std::vector< SPheno_cross_section_line* >::const_iterator
           line_iterator = cross_sections.begin();
         cross_sections.end() > line_iterator;
         ++line_iterator )
      {

        if( 2 == (*line_iterator)->get_number_of_final_state_particles() )
          {

            if( ( ( first_PDG_code
                    == (*line_iterator)->get_first_particle_code() )
                  &&
                  ( second_PDG_code
                      == (*line_iterator)->get_second_particle_code() ) )
                  ||
                  ( ( first_PDG_code
                        == (*line_iterator)->get_second_particle_code() )
                    &&
                    ( second_PDG_code
                        == (*line_iterator)->get_first_particle_code() ) ) )
              // if the requested codes were found, in either order...
              {

                return_value = (*line_iterator)->get_cross_section();
                // note the cross-section.
                line_iterator = cross_sections.end();
                // stop looking.

              }

          }

      }

    return return_value;

  }

  double
  SPheno_cross_section_set::get_cross_section_for_PDG_codes(
                                                      int const first_PDG_code,
                                                     int const second_PDG_code,
                                                     int const third_PDG_code )
  const
  {

    double return_value = CppSLHA_global::really_wrong_value;

    for( std::vector< SPheno_cross_section_line* >::const_iterator
           line_iterator = cross_sections.begin();
         cross_sections.end() > line_iterator;
         ++line_iterator )
      {

        if( 3 == (*line_iterator)->get_number_of_final_state_particles() )
          {

            if( ( ( first_PDG_code
                    == (*line_iterator)->get_first_particle_code() )
                  &&
                  ( second_PDG_code
                      == (*line_iterator)->get_second_particle_code() )
                  &&
                  ( third_PDG_code
                      == (*line_iterator)->get_third_particle_code() ) )
                // if the codes were found in the order 1, 2, 3...
                ||
                ( ( first_PDG_code
                      == (*line_iterator)->get_first_particle_code() )
                  &&
                  ( second_PDG_code
                      == (*line_iterator)->get_third_particle_code() )
                  &&
                  ( third_PDG_code
                      == (*line_iterator)->get_second_particle_code() ) )
                // or if the codes were found in the order 1, 3, 2...
                ||
                ( ( first_PDG_code
                      == (*line_iterator)->get_second_particle_code() )
                  &&
                  ( second_PDG_code
                      == (*line_iterator)->get_first_particle_code() )
                  &&
                  ( third_PDG_code
                      == (*line_iterator)->get_third_particle_code() ) )
                // or if the codes were found in the order 2, 1, 3...
                ||
                ( ( first_PDG_code
                      == (*line_iterator)->get_second_particle_code() )
                  &&
                  ( second_PDG_code
                      == (*line_iterator)->get_third_particle_code() )
                  &&
                  ( third_PDG_code
                      == (*line_iterator)->get_first_particle_code() ) )
                // or if the codes were found in the order 2, 3, 1...
                ||
                ( ( first_PDG_code
                       == (*line_iterator)->get_third_particle_code() )
                  &&
                  ( second_PDG_code
                      == (*line_iterator)->get_first_particle_code() )
                  &&
                  ( third_PDG_code
                      == (*line_iterator)->get_second_particle_code() ) )
                // or if the codes were found in the order 3, 1, 2...
                ||
                ( ( first_PDG_code
                      == (*line_iterator)->get_third_particle_code() )
                  &&
                  ( second_PDG_code
                      == (*line_iterator)->get_second_particle_code() )
                  &&
                  ( third_PDG_code
                      == (*line_iterator)->get_first_particle_code() ) ) )
                      // or if the codes were found in the order 3, 2, 1...
              {

                return_value = (*line_iterator)->get_cross_section();
                // note the cross-section.
                line_iterator = cross_sections.end();
                // stop looking.

              }

          }

      }

    return return_value;

  }


  std::string
  SPheno_cross_section_set::get_as_string()
  const
  {

    // debugging:
    /*std::cout
    << std::endl
    << "debugging! SPheno_cross_section_set::get_as_string() called.";
    std::cout << std::endl;*/

    std::stringstream return_string_builder;

    int entries_written = 0;

    return_string_builder
    << "XS "
    << CppSLHA_global::SLHA_int( first_incoming_particle_code,
                                 9 ) << " "
    << CppSLHA_global::SLHA_int( second_incoming_particle_code,
                                 9 ) << " "
    << CppSLHA_global::SLHA_double( center_of_momentum_energy ) << " "
    << CppSLHA_global::SLHA_double( first_incoming_particle_polarization )
    << " "
    << CppSLHA_global::SLHA_double( second_incoming_particle_polarization )
    << ""
    << CppSLHA_global::SLHA_int( include_ISR_flag,
                                 2 ) << "  # "
    << *(reference_spectrum->get_PDG_coded_name(
                                                first_incoming_particle_code ))
    << ", "
    << *(reference_spectrum->get_PDG_coded_name(
                                               second_incoming_particle_code ))
    << ", c.o.m. energy, polarizations, ISR on/off";
    return_string_builder << std::endl;
    return_string_builder << "# cross-section/fb   NDA   product PDG codes";
    return_string_builder << std::endl;

    for( std::vector< SPheno_cross_section_line* >::const_iterator
           line_iterator = cross_sections.begin();
         cross_sections.end() > line_iterator;
         ++line_iterator )
      {

        ++entries_written;

        int number_of_products
        = (*line_iterator)->get_number_of_final_state_particles();
        int first_product_code = (*line_iterator)->get_first_particle_code();
        int second_product_code
        = (*line_iterator)->get_second_particle_code();
        int third_product_code = (*line_iterator)->get_third_particle_code();
        return_string_builder
        << "   "
        << CppSLHA_global::SLHA_double( (*line_iterator)->get_cross_section() )
        << "   "
        << CppSLHA_global::SLHA_int( number_of_products,
                                     2 ) << "   "
        << CppSLHA_global::SLHA_int( first_product_code,
                                     9 ) << " "
        << CppSLHA_global::SLHA_int( second_product_code,
                                     9 ) << " ";
        if( 3 == number_of_products )
          {

            return_string_builder
            << CppSLHA_global::SLHA_int( third_product_code,
                                         9 ) << " ";

          }
        return_string_builder
        << "  # "
        << *(reference_spectrum->get_PDG_coded_name( first_product_code ))
        << " "
        << *(reference_spectrum->get_PDG_coded_name( second_product_code ));
        if( 3 == number_of_products )
          {

            return_string_builder
            << " "
            << *(reference_spectrum->get_PDG_coded_name( third_product_code ));

          }
        return_string_builder << std::endl;

      } // end of for loop over shown entries.

    if( 1 > entries_written )
      {

        return_string_builder
        << "# nothing recorded from original file for some reason."
        << std::endl;

      }

    // debugging:
    /*std::cout << std::endl << "about to return";
    std::cout << std::endl << return_string_builder.str();
    std::cout << std::endl;*/

    return return_string_builder.str();

  }


  void
  SPheno_cross_section_set::add_cross_section(
                                          std::stringstream* const input_line )
  {

    cross_sections.push_back(
                     ( new SPheno_cross_section_line( input_line ) ) );

  }



  SPheno_cross_sections_BLOCK::SPheno_cross_sections_BLOCK(
                                            std::string const given_BLOCK_name,
                                          std::string const given_name_comment,
                              EW_scale_spectrum const* const given_spectrum ) :
    SLHA_BLOCK( given_BLOCK_name,
                given_name_comment ),
    current_cross_section( NULL ),
    reference_spectrum( given_spectrum )
  {

    comment_strings = NULL;

  }

  SPheno_cross_sections_BLOCK::~SPheno_cross_sections_BLOCK()
  {

    for( std::vector< SPheno_cross_section_set* >::iterator
           deletion_iterator = entries.begin();
         entries.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

  }


  void
  SPheno_cross_sections_BLOCK::make_BLOCK_as_string()
  {

    std::stringstream return_string_builder;

    int entries_written = 0;

    return_string_builder << "BLOCK " << BLOCK_name;
    return_string_builder << " #";
    if( !( name_comment_string.empty() ) )
      {

        return_string_builder << " " << name_comment_string;

      }
    return_string_builder << std::endl;

    for( std::vector< SPheno_cross_section_set* >::const_iterator
           set_iterator = entries.begin();
         entries.end() > set_iterator;
         ++set_iterator )
      {

        ++entries_written;

        return_string_builder << (*set_iterator)->get_as_string();

      } // end of for loop over shown entries.

    if( 1 > entries_written )
      {

        return_string_builder
        << "# nothing recorded from original file for some reason."
        << std::endl;

      }

    BLOCK_as_string.assign( return_string_builder.str() );

  }


  void
  SPheno_cross_sections_BLOCK::read_input(
                                          std::stringstream* const input_line )
  {

    std::string input_line_as_string( input_line->str() );
    std::size_t substring_found_position = input_line_as_string.find( "XS" );
    // look to see if this line is a header for a SPheno cross-section.

    if( substring_found_position != std::string::npos )
    // if it is a header...
      {

        current_cross_section
        = new SPheno_cross_section_set( reference_spectrum,
                                        input_line );
        entries.push_back( current_cross_section );

      }
    else
      // otherwise, we are recording a line for current_cross_section:
      {

        if( NULL == current_cross_section )
          // if somehow there was not a header yet, & thus
          // current_cross_section does not point at allocated memory...
          {

            std::cout
            << std::endl
            << "CppSLHA::error! tried to read in a SPheno cross-section, but"
            << " did not find a line beginning with XS directly after the"
            << " declaration of the SPHENOCROSSSECTIONS BLOCK! ignoring lines"
            << " until the XS header line is found (or until another BLOCK is"
            << " found).";
            std::cout << std::endl;

          }
        else
          {

            current_cross_section->add_cross_section( input_line );

          }

      }

  }

} // end of CppSLHA namespace.

