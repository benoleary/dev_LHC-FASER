/*
 * CppSLHA_BLOCK.hpp
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

#ifndef CPPSLHA_BLOCK_HPP_
#define CPPSLHA_BLOCK_HPP_
// this should prevent this class being included multiple times.

// this is the header file for a C++ class to store SLHA BLOCKs for CppSLHA.

#include "CppSLHA_EW_scale_spectrum.hpp"
#include "CppSLHA_global_stuff.hpp"

/* this class stores SLHA BLOCKs for CppSLHA.
 */

/* I tried to do this with fancy templating, but
 * 1) templates cannot be split across header & source files.
 * 2) std::vectors don't seem to like taking template types.
 * 3) I was only going to use it for < int, int > & < int, double > anyway.
 */


namespace CppSLHA
{

  // this class is for ease of pairing doubles with ints.
  class indexed_double
  {

  public:

    indexed_double( int const given_index_int,
                    double const given_value_double )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~indexed_double()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    int
    get_index()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    set_index( int const input_index_int )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    double
    get_double()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    set_double( double const input_value_double )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    int index_int;
    double value_double;

  };  // end of indexed_double class.

  // this class is for ease of collecting indexed_double objects.
  class indexed_double_set
  {

  public:

    indexed_double_set()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~indexed_double_set()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    reset()
    // this drops all held data.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    int
    get_lowest_index()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    int
    get_highest_index()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    bool
    has_entry( int const requested_index )
    const
    // this checks to see if the entry exists in this set.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    bool
    has_entry( int const requested_index,
               double* const double_for_entry )
    const
    // this checks to see if the entry exists in this set & puts its value
    // into the provided memory address for a double.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_entry( int const requested_index )
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    set_entry( int const requested_index,
               double const given_double )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    std::vector< indexed_double* > entries;

    int lowest_index;
    int highest_index;

  }; // end of indexed_double_set class.


  // this class is for ease of pairing strings with single ints.
  class single_indexed_comment
  {

  public:

    single_indexed_comment( int const given_index_int,
                            std::string const* const given_comment_string )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~single_indexed_comment()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    int
    get_index()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    std::string const*
    get_comment()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    set_comment( std::string const* const given_comment_string )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    int index_int;
    std::string comment_string;

  };  // end of single_indexed_comment class.


  // this class is for ease of pairing strings with pairs of ints.
  class pair_indexed_comment
  {

  public:

    pair_indexed_comment( int const given_first_index_int,
                          int const given_second_index_int,
                          std::string const* const given_comment_string )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~pair_indexed_comment()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    int
    get_first_index()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    int
    get_second_index()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    std::string const*
    get_comment()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    set_comment( std::string const* const given_comment_string )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    int first_index_int;
    int second_index_int;
    std::string comment_string;

  };  // end of pair_indexed_comment class.


  // this class is for ease of pairing strings with triplets of ints.
  class triplet_indexed_comment
  {

  public:

    triplet_indexed_comment( int const given_first_index_int,
                             int const given_second_index_int,
                             int const given_third_index_int,
                             std::string const* const given_comment_string )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~triplet_indexed_comment()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    int
    get_first_index()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    int
    get_second_index()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    int
    get_third_index()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    std::string const*
    get_comment()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    set_comment( std::string const* const given_comment_string )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    int first_index_int;
    int second_index_int;
    int third_index_int;
    std::string comment_string;

  };  // end of triplet_indexed_comment class.


  // this abstract base class is for ease of matching comments to entries of
  // a SLHA_BLOCK for printing out.
  class comment_organizer
  {

  public:

    comment_organizer()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~comment_organizer()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    virtual bool
    has_comment( int const requested_index,
                 std::string* const string_for_comment )
    const
    // this looks for the comment for the requested index, & if it is found,
    // writes it to the provided string, while returning true.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual void
    set_comment( int const requested_index,
                 std::string const* const given_comment )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual bool
    has_comment( int const requested_first_index,
                 int const requested_second_index,
                 std::string* const string_for_comment )
    const
    // this looks for the comment for the requested index, & if it is found,
    // writes it to the provided string, while returning true.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual void
    set_comment( int const requested_first_index,
                 int const requested_second_index,
                 std::string const* const given_comment )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual bool
    has_comment( int const requested_first_index,
                 int const requested_second_index,
                 int const requested_third_index,
                 std::string* const string_for_comment )
    const
    // this looks for the comment for the requested index, & if it is found,
    // writes it to the provided string, while returning true.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual void
    set_comment( int const requested_first_index,
                 int const requested_second_index,
                 int const requested_third_index,
                 std::string const* const given_comment )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  //protected:

    // nothing.

  }; // end of comment_organizer class.

  // this class is for ease of matching comments to entries of a
  // one_index_BLOCK for printing out.
  class one_index_comment_organizer : public comment_organizer
  {

  public:

    one_index_comment_organizer()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~one_index_comment_organizer()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    bool
    has_comment( int const requested_index,
                 std::string* const string_for_comment )
    const
    // this looks for the comment for the requested index, & if it is found,
    // writes it to the provided string, while returning true.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    set_comment( int const requested_index,
                 std::string const* const given_comment )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    std::vector< single_indexed_comment* > comment_strings;

  }; // end of one_index_comment_organizer class.

  // this class is for ease of matching comments to entries of a
  // two_index_BLOCK for printing out.
  class two_index_comment_organizer : public comment_organizer
  {

  public:

    two_index_comment_organizer()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~two_index_comment_organizer()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    bool
    has_comment( int const requested_first_index,
                 int const requested_second_index,
                 std::string* const string_for_comment )
    const
    // this looks for the comment for the requested indices, & if it is found,
    // writes it to the provided string, while returning true.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    set_comment( int const requested_first_index,
                 int const requested_second_index,
                 std::string const* const given_comment )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    std::vector< pair_indexed_comment* > comment_strings;

  }; // end of two_index_comment_organizer class

  // this class is for ease of matching comments to entries of a
  // two_index_BLOCK for printing out.
  class three_index_comment_organizer : public comment_organizer
  {

  public:

    three_index_comment_organizer()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~three_index_comment_organizer()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    bool
    has_comment( int const requested_first_index,
                 int const requested_second_index,
                 int const requested_third_index,
                 std::string* const string_for_comment )
    const
    // this looks for the comment for the requested indices, & if it is found,
    // writes it to the provided string, while returning true.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    set_comment( int const requested_first_index,
                 int const requested_second_index,
                 int const requested_third_index,
                 std::string const* const given_comment )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    std::vector< triplet_indexed_comment* > comment_strings;

  }; // end of three_index_comment_organizer class


  // this is an abstract base class for storing all BLOCKs polymorphically.
  class SLHA_BLOCK
  {

  public:

    SLHA_BLOCK( std::string const given_BLOCK_name,
                std::string const given_name_comment )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~SLHA_BLOCK()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    bool
    has_recorded()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    note_recorded()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    SLHA_BLOCK*
    get_extra_BLOCK()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    SLHA_BLOCK const*
    inspect_extra_BLOCK()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    std::string const*
    get_name()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    double
    get_Q()
    const
    // this returns the scale of the BLOCK or returns
    // CppSLHA_global::really_wrong_value.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    set_Q( double const input_entry )
    // this sets the scale & records that it was set, if it had not been set.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    set_comment( int const requested_index,
                 std::string const given_comment )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    set_comment( int const requested_first_index,
                 int const requested_second_index,
                 std::string const given_comment )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    set_comment( int const requested_first_index,
                 int const requested_second_index,
                 int const requested_third_index,
                 std::string const given_comment )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    reset()
    // this wipes any saved string representation of the BLOCK, then calls the
    // specific_reset function which is specific to each derived class.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual SLHA_BLOCK*
    make_extra_BLOCK()
    = 0;

    virtual void
    read_input( std::stringstream* const input_line )
    = 0;


    virtual bool
    has_entry( double* const double_for_entry )
    const
    // this checks to see if the entry exists in this set & puts its value
    // into the provided memory address for a double.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual bool
    has_entry( int const requested_index,
               double* const double_for_entry )
    const
    // this checks to see if the entry exists in this set & puts its value
    // into the provided memory address for a double.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual bool
    has_entry( int const requested_first_index,
               int const requested_second_index,
               double* const double_for_entry )
    const
    // this checks to see if the entry exists in this set & puts its value
    // into the provided memory address for a double.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual bool
    has_entry( int const requested_first_index,
               int const requested_second_index,
               int const requested_third_index,
               double* const double_for_entry )
    const
    // this checks to see if the entry exists in this set & puts its value
    // into the provided memory address for a double.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    virtual double
    get_entry()
    const
    // this should return the entry or return
    // CppSLHA_global::really_wrong_value.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual double
    get_entry( int const requested_index )
    const
    // this should return the entry or return
    // CppSLHA_global::really_wrong_value.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual double
    get_entry( int const requested_first_index,
               int const requested_second_index )
    const
    // this should return the entry or return
    // CppSLHA_global::really_wrong_value.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual double
    get_entry( int const requested_first_index,
               int const requested_second_index,
               int const requested_third_index )
    const
    // this should return the entry or return
    // CppSLHA_global::really_wrong_value.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    virtual bool
    has_string( std::string* const string_for_entry )
    const
    // this looks for the entry for the requested indices, & if it is found,
    // writes it to the provided string, while returning true.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual bool
    has_string( int const requested_index,
                std::string* const string_for_entry )
    const
    // this looks for the entry for the requested indices, & if it is found,
    // writes it to the provided string, while returning true.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual bool
    has_string( int const requested_first_index,
                int const requested_second_index,
                std::string* const string_for_entry )
    const
    // this looks for the entry for the requested indices, & if it is found,
    // writes it to the provided string, while returning true.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual bool
    has_string( int const requested_first_index,
                int const requested_second_index,
                int const requested_third_index,
                std::string* const string_for_entry )
    const
    // this looks for the entry for the requested indices, & if it is found,
    // writes it to the provided string, while returning true.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    virtual void
    set_entry( double const input_entry )
    // this sets the entry & records that it was set, if it had not been set.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual void
    set_entry( int const requested_index,
               double const input_entry )
    // this sets the entry & records that it was set, if it had not been set.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual void
    set_entry( int const requested_first_index,
               int const requested_second_index,
               double const input_entry )
    // this sets the entry & records that it was set, if it had not been set.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual void
    set_entry( int const requested_first_index,
               int const requested_second_index,
               int const requested_third_index,
               double const input_entry )
    // this sets the entry & records that it was set, if it had not been set.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual void
    set_string( std::string const* const input_string )
    // this sets the entry & records that it was set, if it had not been set.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual void
    set_string( int const requested_first_index,
                std::string const* const input_string )
    // this sets the entry & records that it was set, if it had not been set.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual void
    set_string( int const requested_first_index,
                int const requested_second_index,
                std::string const* const input_string )
    // this sets the entry & records that it was set, if it had not been set.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual void
    set_string( int const requested_first_index,
                int const requested_second_index,
                int const requested_third_index,
                std::string const* const input_string )
    // this sets the entry & records that it was set, if it had not been set.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    virtual std::string const*
    get_as_string()
    // this returns the BLOCK in SLHA-format ASCII.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    virtual void
    make_BLOCK_as_string()
    // this writes the BLOCK in SLHA-format ASCII into BLOCK_as_string.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    write_entry_as_integer( bool const given_write_as_int_flag )
    /* if given true, subsequent writing of the output as a string will use
     * ints for the entries rather than doubles. hence this also resets
     * BLOCK_as_string.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    std::string BLOCK_name;

    double scale_value;
    bool scale_recorded;

    std::string name_comment_string;

    comment_organizer* comment_strings;

    // the following are to try to cope with multiple BLOCKs of the same name
    // with different Qs:
    bool has_recorded_flag;
    bool has_extra_BLOCK;

    SLHA_BLOCK* extra_BLOCK;
    std::string BLOCK_as_string;

    bool write_as_int_flag;

    virtual void
    specific_reset()
    // this is to be over-written in each derived class to do whatever needs to
    // be done to properly specific_reset the BLOCK.
    = 0;

  };  // end of SLHA_BLOCK class.

  // this stores the acceptanceValues for a SLHA BLOCK that can be represented as
  // a single number.
  class zero_index_BLOCK : public SLHA_BLOCK
  {

  public:

    zero_index_BLOCK( std::string const given_BLOCK_name,
                      std::string const given_name_comment )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~zero_index_BLOCK()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    SLHA_BLOCK*
    make_extra_BLOCK()
    // this makes a blank copy of this BLOCK.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    bool
    has_entry()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    bool
    has_entry( double* const double_for_entry )
    const
    // this checks to see if the entry exists in this set & puts its value
    // into the provided memory address for a double.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_entry()
    const
    // this returns the entry or returns CppSLHA_global::really_wrong_value.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    set_entry( double const input_entry )
    // this sets the entry & records that it was set, if it had not been set.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    make_BLOCK_as_string()
    // this writes the BLOCK in SLHA-format ASCII into BLOCK_as_string.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    read_input( std::stringstream* const input_line )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    double entry_value;
    bool entry_recorded;


    void
    specific_reset()
    // this drops all held data apart from the name & comments.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  }; // end of zero_index_BLOCK class.

  // this stores the acceptanceValues for a SLHA BLOCK that can be represented as
  // a 1-dimensional array.
  class one_index_BLOCK : public SLHA_BLOCK
  {

  public:

    one_index_BLOCK( std::string const given_BLOCK_name,
                     std::string const given_name_comment )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~one_index_BLOCK()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    SLHA_BLOCK*
    make_extra_BLOCK()
    // this makes a blank copy of this BLOCK.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    bool
    has_entry( int const requested_element,
               double* const double_for_entry )
    const
    // this checks to see if the entry exists in this set & puts its value
    // into the provided memory address for a double.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    double
    get_entry( int const requested_element )
    const
    // this returns the entry in non-insane notation (i.e. the first entry is
    // sought with 1, not 0) or returns CppSLHA_global::really_wrong_value.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    set_entry( int const requested_element,
               double const input_entry )
    /* this either overwrites the value at entry requested_element, or adds it
     * to the std::vector if it didn't already have it, possibly filling in
     * CppSLHA_global::really_wrong_value entries along the way.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    make_BLOCK_as_string()
    // this writes the BLOCK in SLHA-format ASCII into BLOCK_as_string.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    read_input( std::stringstream* const input_line )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    std::vector< double > entries;
    std::vector< bool > recorded;


    void
    specific_reset()
    // this drops all held data apart from the name & comments.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  }; // end of one_index_BLOCK class.

  // this is for storing doubles when there would be large std::vectors
  // full of zero entries but for a few.
  class sparse_one_index_BLOCK : public SLHA_BLOCK
  {

  public:

    sparse_one_index_BLOCK( std::string const given_BLOCK_name,
                            std::string const given_name_comment )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~sparse_one_index_BLOCK()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    SLHA_BLOCK*
    make_extra_BLOCK()
    // this makes a blank copy of this BLOCK.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    bool
    has_entry( int const requested_element,
               double* const double_for_entry )
    const
    // this checks to see if the entry exists in this set & puts its value
    // into the provided memory address for a double.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_entry( int const requested_element )
    const
    // this either returns the value associated with requested_element
    // or CppSLHA_global::really_wrong_value if there is none.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    set_entry( int const requested_element,
               double const input_entry )
    // this either overwrites the value at entry requested_element,
    // or adds it to the number_pair_set if it didn't have it.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    make_BLOCK_as_string()
    // this writes the BLOCK in SLHA-format ASCII into BLOCK_as_string.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    read_input( std::stringstream* const input_line )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    indexed_double_set entries;


    void
    specific_reset()
    // this drops all held data apart from the name & comments.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };  // end of sparse_one_index_BLOCK class.

  // this is for storing string entries when there would be large
  // std::vectors full of zero entries but for a few.
  class sparse_one_index_string_BLOCK : public SLHA_BLOCK
  {

  public:

    sparse_one_index_string_BLOCK( std::string const given_BLOCK_name,
                                   std::string const given_name_comment )
      /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~sparse_one_index_string_BLOCK()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    SLHA_BLOCK*
    make_extra_BLOCK()
    // this makes a blank copy of this BLOCK.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    bool
    has_string( int const requested_element,
                std::string* const string_for_entry )
    const
    // this looks for the entry for the requested indices, & if it is found,
    // writes it to the provided string, while returning true.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    set_entry( int const requested_element,
               std::string const input_entry )
    // this either overwrites the value at entry requested_element,
    // or adds it to the set if it didn't have it.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    make_BLOCK_as_string()
    // this writes the BLOCK in SLHA-format ASCII into BLOCK_as_string.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    read_input( std::stringstream* const input_line )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    one_index_comment_organizer* entries;
    // we use the comment organizer thing as a means of storing the strings
    // with their indices.

    indexed_double_set indices;
    // this is because I am lazy and just want the functionality of recording
    // the lowest & highest indices for the entries.


    void
    specific_reset()
    // this drops all held data apart from the name & comments.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };  // end of sparse_one_index_string_BLOCK class.

  // this is for storing string entries when there can be multiple lines with
  // the same index.
  class multiple_same_one_index_string_BLOCK : public SLHA_BLOCK
  {

  public:

    multiple_same_one_index_string_BLOCK( std::string const given_BLOCK_name,
                                         std::string const given_name_comment )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~multiple_same_one_index_string_BLOCK()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    SLHA_BLOCK*
    make_extra_BLOCK()
    // this makes a blank copy of this BLOCK.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    bool
    has_string( int const requested_element,
                std::string* const string_for_entry )
    const
    // this looks for the entry for the requested indices, & if it is found,
    // writes it to the provided string, while returning true.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    bool
    has_string( int const requested_first_element,
                int const requested_second_element,
                std::string* const string_for_entry )
    const
    // this looks for the entry for the requested indices, & if it is found,
    // writes it to the provided string, while returning true.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    set_entry( int const requested_first_element,
               int const requested_second_element,
               std::string const input_entry )
    // this either overwrites the value at entry requested_element,
    // or adds it to the set if it didn't have it.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    make_BLOCK_as_string()
    // this writes the BLOCK in SLHA-format ASCII into BLOCK_as_string.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    read_input( std::stringstream* const input_line )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    two_index_comment_organizer* entries;
    // we use the comment organizer thing as a means of storing the strings
    // with their indices.

    int lowest_index;
    int highest_index;
    int largest_hidden_index;
    // there is only 1 dimension to the indices, though each value for the
    // index can have multiple entries.


    void
    specific_reset()
    // this drops all held data apart from the name & comments.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  }; // end of multiple_same_one_index_string_BLOCK class.

  // this stores the acceptanceValues for a SLHA BLOCK that can be represented as
  // a 2-dimensional array.
  class two_index_BLOCK : public SLHA_BLOCK
  {

  public:

    two_index_BLOCK( std::string const given_BLOCK_name,
                     std::string const given_name_comment )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~two_index_BLOCK()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    SLHA_BLOCK*
    make_extra_BLOCK()
    // this makes a blank copy of this BLOCK.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    bool
    has_entry( int const requested_first_element,
               int const requested_second_element,
               double* const double_for_entry )
    const
    // this checks to see if the entry exists in this set & puts its value
    // into the provided memory address for a double.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_entry( int const requested_first_element,
               int const requested_second_element )
    const
    // this returns the entry in non-insane notation (i.e. the first entry is
    // sought with 1, not 0) or returns CppSLHA_global::really_wrong_value.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    set_entry( int const requested_first_element,
               int const requested_second_element,
               double const input_entry )
    /* this either overwrites the value at entry requested_element, or adds
     * it to the std::vector if it didn't already have it, possibly filling
     * in CppSLHA_global::really_wrong_value entries along the way (in the 2nd
     * argument, just empty std::vectors in the 1st argument).
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    make_BLOCK_as_string()
    // this writes the BLOCK in SLHA-format ASCII into BLOCK_as_string.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    read_input( std::stringstream* const input_line )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    std::vector< std::vector< double >* > entries;
    std::vector< std::vector< bool >* > recorded;


    void
    specific_reset()
    // this drops all held data apart from the name & comments.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  }; // end of two_index_BLOCK class.


  // this stores the acceptanceValues for a SLHA BLOCK that can be represented as a
  // 3-dimensional array.
  class three_index_BLOCK : public SLHA_BLOCK
  {

  public:

    three_index_BLOCK( std::string const given_BLOCK_name,
                       std::string const given_name_comment )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~three_index_BLOCK()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    bool
    has_entry( int const requested_first_element,
               int const requested_second_element,
               int const requested_third_element,
               double* const double_for_entry )
    const
    // this checks to see if the entry exists in this set & puts its value
    // into the provided memory address for a double.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    SLHA_BLOCK*
    make_extra_BLOCK()
    // this makes a blank copy of this BLOCK.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    double
    get_entry( int const requested_first_element,
               int const requested_second_element,
               int const requested_third_element )
    const
    // this returns the entry in non-insane notation (i.e. the first entry is
    // sought with 1, not 0) or returns CppSLHA_global::really_wrong_value.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    set_entry( int const requested_first_element,
               int const requested_second_element,
               int const requested_third_element,
               double const input_entry )
    /* this either overwrites the value at entry requested_element, or adds
     * it to the std::vector if it didn't already have it, possibly filling
     * in CppSLHA_global::really_wrong_value entries along the way (in the 3rd
     * argument, just empty std::vectors in the 1st & 2nd arguments).
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    make_BLOCK_as_string()
    // this writes the BLOCK in SLHA-format ASCII into BLOCK_as_string.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    read_input( std::stringstream* const input_line )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    std::vector< std::vector< std::vector< double >* >* > entries;
    std::vector< std::vector< std::vector< bool >* >* > recorded;


    void
    specific_reset()
    // this drops all held data apart from the name & comments.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  }; // end of three_index_BLOCK class.


  // this is a BLOCK type for accessing the masses in the manner of a normal
  // BLOCK
  class MASS_BLOCK : public SLHA_BLOCK
  {

  public:

    MASS_BLOCK( EW_scale_spectrum* given_spectrum )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~MASS_BLOCK()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    bool
    has_entry( int const requested_element,
               double* const double_for_entry )
    const
    // this checks to see if the entry exists in this set & puts its value
    // into the provided memory address for a double.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    SLHA_BLOCK*
    make_extra_BLOCK()
    // this makes a blank copy of this BLOCK.
    /* actually, while I *could* re-write it all such that each MASS BLOCK has
     * its own EW_scale_spectrum, or re-write the EW_scale_spectrum class so
     * that each particle can have a scale-dependent mass, I'm not going to
     * bother right now. I think probably changing the EW_scale_spectrum class
     * would be the best way to do it.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    double
    get_entry( int const requested_element )
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    set_entry( int const requested_element,
               double const input_entry )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    make_BLOCK_as_string()
    // this writes the BLOCK in SLHA-format ASCII into BLOCK_as_string.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    read_input( std::stringstream* const input_line )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    EW_scale_spectrum* reference_spectrum;


    void
    specific_reset()
    // this drops all held data apart from the name & comments.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };  // end of MASS_BLOCK class.

  /* this is for storing entries for the SPHENOINPUT BLOCK, which has the
   * annoying property of a set of lines which can be repeated, so this is
   * implemented with the "hidden index" method of
   * multiple_same_one_index_string_BLOCK.
   */
  class SPheno_input_BLOCK : public SLHA_BLOCK
  {

  public:

    SPheno_input_BLOCK( std::string const given_BLOCK_name,
                        std::string const given_name_comment )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~SPheno_input_BLOCK()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    SLHA_BLOCK*
    make_extra_BLOCK()
    // this makes a blank copy of this BLOCK.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    bool
    has_entry( int const requested_element,
               double* const double_for_entry )
    const
    // this checks to see if the entry exists in this set & puts its value
    // into the provided memory address for a double.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    bool
    has_entry( int const requested_first_element,
               int const requested_second_element,
               double* const double_for_entry )
    const
    // this checks to see if the entry exists in this set & puts its value
    // into the provided memory address for a double.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    set_entry( int const requested_element,
               double const input_entry )
    // this either overwrites the value at entry ( requested_element, 1 ),
    // or adds it to the set if it didn't have it.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    set_entry( int const requested_first_element,
               int const requested_second_element,
               double const input_entry )
    /* this either overwrites the value at entry
     * ( requested_first_element, requested_second_element ), or adds it to
     * the set if it didn't have it.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    make_BLOCK_as_string()
    // this writes the BLOCK in SLHA-format ASCII into BLOCK_as_string.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    read_input( std::stringstream* const input_line )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    std::map< std::pair< int, int >, double > entries;
    // the 1st int of the pair is the "visible" index, while the 2nd is the
    // "hidden" index.

    int lowest_index;
    int highest_index;
    int largest_hidden_index;
    // there is only 1 dimension to the indices, though each value for the
    // index can have multiple entries.


    void
    specific_reset()
    // this drops all held data apart from the name & comments.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  }; // end of SPheno_input_BLOCK class.

  /* this is for storing entries for the SPHENOCROSSSECTIONS BLOCK, which has
   * an annoying DECAY-like structure. it is a header set of ints & doubles,
   * & a std::vector of lines of sets of ints & doubles.
   */
  class SPheno_cross_section_line
  {

  public:

    SPheno_cross_section_line( std::stringstream* const given_input )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~SPheno_cross_section_line()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    int
    get_number_of_final_state_particles()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    int
    get_first_particle_code()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    int
    get_second_particle_code()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    int
    get_third_particle_code()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_cross_section()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    int number_of_final_state_particles;
    int first_particle_code;
    int second_particle_code;
    int third_particle_code;
    double cross_section;

  }; // end of SPheno_cross_section_line class.

  /* this is for storing entries for the SPHENOCROSSSECTIONS BLOCK, which has
   * an annoying DECAY-like structure. it is a header set of ints & doubles,
   * & a std::vector of SPheno_cross_section_line pointers.
   */
  class SPheno_cross_section_set
  {

  public:

    SPheno_cross_section_set( EW_scale_spectrum const* const given_spectrum,
                              std::stringstream* const given_header )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~SPheno_cross_section_set()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    int
    get_first_incoming_particle_code()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    int
    get_second_incoming_particle_code()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_center_of_momentum_energy()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_first_incoming_particle_polarization()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_second_incoming_particle_polarization()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    int
    get_include_ISR_flag()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    std::vector< SPheno_cross_section_line* >*
    get_cross_section_lines()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    double
    get_cross_section_for_PDG_codes( int const first_PDG_code,
                                     int const second_PDG_code )
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_cross_section_for_PDG_codes( int const first_PDG_code,
                                     int const second_PDG_code,
                                     int const third_PDG_code )
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    std::string
    get_as_string()
    const
    // this writes the SPheno_cross_section_set in SLHA-format ASCII.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    add_cross_section( std::stringstream* const input_line )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    int first_incoming_particle_code;
    int second_incoming_particle_code;
    double center_of_momentum_energy;
    double first_incoming_particle_polarization;
    double second_incoming_particle_polarization;
    int include_ISR_flag;

    std::vector< SPheno_cross_section_line* > cross_sections;

    EW_scale_spectrum const* const reference_spectrum;

  }; // end of SPheno_cross_section_set class.


  /* this is for storing entries for the SPHENOCROSSSECTIONS BLOCK, which has
   * an annoying DECAY-like structure. it is a derived class from
   * CppSLHA_BLOCK for purposes of reading the SLHA file, but most of its
   * functionality is not available to the base class, and as such can only be
   * used by CppSLHA_SPheno, which knows about this derived class.
   */
  class SPheno_cross_sections_BLOCK : public SLHA_BLOCK
  {

  public:

    SPheno_cross_sections_BLOCK( std::string const given_BLOCK_name,
                                 std::string const given_name_comment,
                                EW_scale_spectrum const* const given_spectrum )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~SPheno_cross_sections_BLOCK()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    SLHA_BLOCK*
    make_extra_BLOCK()
    // this makes a blank copy of this BLOCK.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    std::vector< SPheno_cross_section_set* >*
    get_SPheno_cross_section_vector()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    make_BLOCK_as_string()
    // this writes the BLOCK in SLHA-format ASCII into BLOCK_as_string.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    read_input( std::stringstream* const input_line )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    std::vector< SPheno_cross_section_set* > entries;
    SPheno_cross_section_set* current_cross_section;
    EW_scale_spectrum const* reference_spectrum;


    void
    specific_reset()
    // this drops all held data apart from the name & comments.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  }; // end of SPheno_cross_sections_BLOCK class.




  // inline functions:

  inline int
  indexed_double::get_index()
  const
  {

    return index_int;

  }

  inline void
  indexed_double::set_index( int const input_index_int )
  {

    index_int = input_index_int;

  }

  inline double
  indexed_double::get_double()
  const
  {

    return value_double;

  }

  inline void
  indexed_double::set_double( double const input_value_double )
  {

    value_double = input_value_double;

  }



  inline void
  indexed_double_set::reset()
  // this drops all held data.
  {

    for( std::vector< indexed_double* >::iterator
           deletion_iterator = entries.begin();
         entries.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }
    entries.clear();

    lowest_index = (int)CppSLHA_global::really_wrong_value;
    highest_index = (int)CppSLHA_global::really_wrong_value;

  }


  inline int
  indexed_double_set::get_lowest_index()
  const
  {

    return lowest_index;

  }

  inline int
  indexed_double_set::get_highest_index()
  const
  {

    return highest_index;

  }



  inline int
  single_indexed_comment::get_index()
  const
  {

    return index_int;

  }

  inline std::string const*
  single_indexed_comment::get_comment()
  const
  {

    return &comment_string;

  }

  inline void
  single_indexed_comment::set_comment(
                                std::string const* const given_comment_string )
  {

    comment_string.assign( *given_comment_string );

  }



  inline int
  pair_indexed_comment::get_first_index()
  const
  {

    return first_index_int;

  }

  inline int
  pair_indexed_comment::get_second_index()
  const
  {

    return second_index_int;

  }

  inline std::string const*
  pair_indexed_comment::get_comment()
  const
  {

    return &comment_string;

  }

  inline void
  pair_indexed_comment::set_comment(
                                std::string const* const given_comment_string )
  {

    comment_string.assign( *given_comment_string );

  }



  inline int
  triplet_indexed_comment::get_first_index()
  const
  {

    return first_index_int;

  }

  inline int
  triplet_indexed_comment::get_second_index()
  const
  {

    return second_index_int;

  }

  inline int
  triplet_indexed_comment::get_third_index()
  const
  {

    return third_index_int;

  }

  inline std::string const*
  triplet_indexed_comment::get_comment()
  const
  {

      return &comment_string;

  }

  inline void
  triplet_indexed_comment::set_comment(
                                std::string const* const given_comment_string )
  {

    comment_string.assign( *given_comment_string );

  }



  inline bool
  comment_organizer::has_comment( int const requested_index,
                                  std::string* const string_for_comment )
  const
  // this looks for the comment for the requested index, & if it is found,
  // writes it to the provided string, while returning true.
  {

    return false;

  }

  inline void
  comment_organizer::set_comment( int const requested_index,
                                  std::string const* const given_comment )
  {

    std::cout
    << std::endl
    << "CppSLHA::error! a BLOCK was asked to record a comment with the wrong"
    << " number of indices";
    std::cout << std::endl << "calling exit( EXIT_FAILURE ).";
    std::cout << std::endl;

    exit( EXIT_FAILURE );

  }

  inline bool
  comment_organizer::has_comment( int const requested_first_index,
                                  int const requested_second_index,
                                  std::string* const string_for_comment )
  const
  // this looks for the comment for the requested index, & if it is found,
  // writes it to the provided string, while returning true.
  {

    return false;

  }

  inline void
  comment_organizer::set_comment( int const requested_first_index,
                                  int const requested_second_index,
                                  std::string const* const given_comment )
  {

    std::cout
    << std::endl
    << "CppSLHA::error! a BLOCK was asked to record a comment with the wrong"
    << " number of indices";
    std::cout << std::endl << "calling exit( EXIT_FAILURE ).";
    std::cout << std::endl;

    exit( EXIT_FAILURE );

  }

  inline bool
  comment_organizer::has_comment( int const requested_first_index,
                                  int const requested_second_index,
                                  int const requested_third_index,
                                  std::string* const string_for_comment )
  const
  // this looks for the comment for the requested index, & if it is found,
  // writes it to the provided string, while returning true.
  {

    return false;

  }

  inline void
  comment_organizer::set_comment( int const requested_first_index,
                                  int const requested_second_index,
                                  int const requested_third_index,
                                  std::string const* const given_comment )
  {

    std::cout
    << std::endl
    << "CppSLHA::error! a BLOCK was asked to record a comment with the wrong"
    << " number of indices";
    std::cout << std::endl << "calling exit( EXIT_FAILURE ).";
    std::cout << std::endl;

    exit( EXIT_FAILURE );

  }



  inline bool
  SLHA_BLOCK::has_recorded()
  const
  {

    return has_recorded_flag;

  }

  inline void
  SLHA_BLOCK::note_recorded()
  {

    has_recorded_flag = true;

  }

  inline SLHA_BLOCK*
  SLHA_BLOCK::get_extra_BLOCK()
  {

    if( has_extra_BLOCK )
      {

        return extra_BLOCK;

      }
    else
      {

        return NULL;

      }

  }

  inline SLHA_BLOCK const*
  SLHA_BLOCK::inspect_extra_BLOCK()
  const
  {

    if( has_extra_BLOCK )
      {

        return extra_BLOCK;

      }
    else
      {

        return NULL;

      }

  }


  inline std::string const*
  SLHA_BLOCK::get_name()
  const
  {

    return &BLOCK_name;

  }


  inline double
  SLHA_BLOCK::get_Q()
  const
  // this returns the scale of the BLOCK or returns
  // CppSLHA_global::really_wrong_value.
  {

    return scale_value;

  }

  inline void
  SLHA_BLOCK::set_Q( double const input_entry )
  // this sets the scale & records that it was set, if it had not been set.
  {

    scale_value = input_entry;
    scale_recorded = true;

  }


  inline void
  SLHA_BLOCK::set_comment( int const requested_index,
                           std::string const given_comment )
  {

    comment_strings->set_comment( requested_index,
                                  &given_comment );

  }

  inline void
  SLHA_BLOCK::set_comment( int const requested_first_index,
                           int const requested_second_index,
                           std::string const given_comment )
  {

    comment_strings->set_comment( requested_first_index,
                                  requested_second_index,
                                  &given_comment );

  }

  inline void
  SLHA_BLOCK::set_comment( int const requested_first_index,
                           int const requested_second_index,
                           int const requested_third_index,
                           std::string const given_comment )
  {

    comment_strings->set_comment( requested_first_index,
                                  requested_second_index,
                                  requested_third_index,
                                  &given_comment );

  }


  inline void
  SLHA_BLOCK::reset()
  // this wipes any saved string representation of the BLOCK, then calls the
  // specific_reset function which is specific to each derived class.
  {

    BLOCK_as_string.clear();

    specific_reset();

  }


  inline bool
  SLHA_BLOCK::has_entry( double* const double_for_entry )
  const
  // this checks to see if the entry exists in this set & puts its value into
  // the provided memory address for a double.
  {

    return false;

  }

  inline bool
  SLHA_BLOCK::has_entry( int const requested_index,
                         double* const double_for_entry )
  const
  // this checks to see if the entry exists in this set & puts its value into
  // the provided memory address for a double.
  {

    return false;

  }

  inline bool
  SLHA_BLOCK::has_entry( int const requested_first_index,
                         int const requested_second_index,
                         double* const double_for_entry )
  const
  // this checks to see if the entry exists in this set & puts its value into
  // the provided memory address for a double.
  {

    return false;

  }

  inline bool
  SLHA_BLOCK::has_entry( int const requested_first_index,
                         int const requested_second_index,
                         int const requested_third_index,
                         double* const double_for_entry )
  const
  // this checks to see if the entry exists in this set & puts its value into
  // the provided memory address for a double.
  {

    return false;

  }


  inline double
  SLHA_BLOCK::get_entry()
  const
  // this should return the entry or return
  // CppSLHA_global::really_wrong_value.
  {

    std::cout
    << std::endl
    << "CppSLHA::error! an SLHA BLOCK was asked for an entry which it had not"
    << " recorded.";
    std::cout
    << std::endl << CppSLHA_global::really_wrong_value_string
    << " has been returned.";
    std::cout << std::endl;

    return CppSLHA_global::really_wrong_value;

  }

  inline double
  SLHA_BLOCK::get_entry( int const requested_index )
  const
  // this should return the entry or return
  // CppSLHA_global::really_wrong_value.
  {

    std::cout
    << std::endl
    << "CppSLHA::error! an SLHA BLOCK was asked for an entry which it had not"
    << " recorded.";
    std::cout
    << std::endl << CppSLHA_global::really_wrong_value_string
    << " has been returned.";
    std::cout << std::endl;

    return CppSLHA_global::really_wrong_value;

  }

  inline double
  SLHA_BLOCK::get_entry( int const requested_first_index,
                         int const requested_second_index )
  const
  // this should return the entry or return
  // CppSLHA_global::really_wrong_value.
  {

    std::cout
    << std::endl
    << "CppSLHA::error! an SLHA BLOCK was asked for an entry which it had not"
    << " recorded.";
    std::cout
    << std::endl << CppSLHA_global::really_wrong_value_string
    << " has been returned.";
    std::cout << std::endl;

    return CppSLHA_global::really_wrong_value;

  }

  inline double
  SLHA_BLOCK::get_entry( int const requested_first_index,
                         int const requested_second_index,
                         int const requested_third_index )
  const
  // this should return the entry or return
  // CppSLHA_global::really_wrong_value.
  {

    std::cout
    << std::endl
    << "CppSLHA::error! an SLHA BLOCK was asked for an entry which it had not"
    << " recorded.";
    std::cout
    << std::endl << CppSLHA_global::really_wrong_value_string
    << " has been returned.";
    std::cout << std::endl;

    return CppSLHA_global::really_wrong_value;

  }


  inline bool
  SLHA_BLOCK::has_string( std::string* const string_for_entry )
  const
  // this looks for the entry for the requested indices, & if it is found,
  // writes it to the provided string, while returning true.
  {

    return false;

  }

  inline bool
  SLHA_BLOCK::has_string( int const requested_index,
                          std::string* const string_for_entry )
  const
  // this looks for the entry for the requested indices, & if it is found,
  // writes it to the provided string, while returning true.
  {

    return false;

  }

  inline bool
  SLHA_BLOCK::has_string( int const requested_first_index,
                          int const requested_second_index,
                          std::string* const string_for_entry )
  const
  // this looks for the entry for the requested indices, & if it is found,
  // writes it to the provided string, while returning true.
  {

    return false;

  }

  inline bool
  SLHA_BLOCK::has_string( int const requested_first_index,
                          int const requested_second_index,
                          int const requested_third_index,
                          std::string* const string_for_entry )
  const
  // this looks for the entry for the requested indices, & if it is found,
  // writes it to the provided string, while returning true.
  {

    return false;

  }


  inline void
  SLHA_BLOCK::set_entry( double const input_entry )
  // this sets the entry & records that it was set, if it had not been set.
  {

    std::cout
    << std::endl
    << "CppSLHA::error! a BLOCK (\"" << BLOCK_name << "\")"
    << " was asked to record an entry with the wrong number of indices";
    std::cout << std::endl << "calling exit( EXIT_FAILURE ).";
    std::cout << std::endl;

    exit( EXIT_FAILURE );

  }

  inline void
  SLHA_BLOCK::set_entry( int const requested_index,
                         double const input_entry )
  // this sets the entry & records that it was set, if it had not been set.
  {

    std::cout
    << std::endl
    << "CppSLHA::error! a BLOCK (\"" << BLOCK_name << "\")"
    << " was asked to record an entry with the wrong number of indices";
    std::cout << std::endl << "calling exit( EXIT_FAILURE ).";
    std::cout << std::endl;

    exit( EXIT_FAILURE );

  }

  inline void
  SLHA_BLOCK::set_entry( int const requested_first_index,
                         int const requested_second_index,
                         double const input_entry )
  // this sets the entry & records that it was set, if it had not been set.
  {

    std::cout
    << std::endl
    << "CppSLHA::error! a BLOCK (\"" << BLOCK_name << "\")"
    << " was asked to record an entry with the wrong number of indices";
    std::cout << std::endl << "calling exit( EXIT_FAILURE ).";
    std::cout << std::endl;

    exit( EXIT_FAILURE );

  }

  inline void
  SLHA_BLOCK::set_entry( int const requested_first_index,
                         int const requested_second_index,
                         int const requested_third_index,
                         double const input_entry )
  // this sets the entry & records that it was set, if it had not been set.
  {

    std::cout
    << std::endl
    << "CppSLHA::error! a BLOCK (\"" << BLOCK_name << "\")"
    << " was asked to record an entry with the wrong number of indices";
    std::cout << std::endl << "calling exit( EXIT_FAILURE ).";
    std::cout << std::endl;

    exit( EXIT_FAILURE );

  }

  inline void
  SLHA_BLOCK::set_string( std::string const* const input_string )
  // this sets the entry & records that it was set, if it had not been set.
  {

    std::cout
    << std::endl
    << "CppSLHA::error! a BLOCK (\"" << BLOCK_name << "\")"
    << " was asked to record an entry with the wrong number of indices";
    std::cout << std::endl << "calling exit( EXIT_FAILURE ).";
    std::cout << std::endl;

    exit( EXIT_FAILURE );

  }

  inline void
  SLHA_BLOCK::set_string( int const requested_first_index,
                          std::string const* const input_string )
  // this sets the entry & records that it was set, if it had not been set.
  {

    std::cout
    << std::endl
    << "CppSLHA::error! a BLOCK (\"" << BLOCK_name << "\")"
    << " was asked to record an entry with the wrong number of indices";
    std::cout << std::endl << "calling exit( EXIT_FAILURE ).";
    std::cout << std::endl;

    exit( EXIT_FAILURE );

  }

  inline void
  SLHA_BLOCK::set_string( int const requested_first_index,
                          int const requested_second_index,
                          std::string const* const input_string )
  // this sets the entry & records that it was set, if it had not been set.
  {

    std::cout
    << std::endl
    << "CppSLHA::error! a BLOCK (\"" << BLOCK_name << "\")"
    << " was asked to record an entry with the wrong number of indices";
    std::cout << std::endl << "calling exit( EXIT_FAILURE ).";
    std::cout << std::endl;

    exit( EXIT_FAILURE );

  }

  inline void
  SLHA_BLOCK::set_string( int const requested_first_index,
                          int const requested_second_index,
                          int const requested_third_index,
                          std::string const* const input_string )
  // this sets the entry & records that it was set, if it had not been set.
  {

    std::cout
    << std::endl
    << "CppSLHA::error! a BLOCK (\"" << BLOCK_name << "\")"
    << " was asked to record an entry with the wrong number of indices";
    std::cout << std::endl << "calling exit( EXIT_FAILURE ).";
    std::cout << std::endl;

    exit( EXIT_FAILURE );

  }


  inline std::string const*
  SLHA_BLOCK::get_as_string()
  // this returns the BLOCK in SLHA-format ASCII.
  {

    if( BLOCK_as_string.empty() )
      {

        make_BLOCK_as_string();

      }

    return &BLOCK_as_string;

  }


  inline void
  SLHA_BLOCK::write_entry_as_integer( bool const given_write_as_int_flag )
  /* if given true, subsequent writing of the output as a string will use
   * ints for the entries rather than doubles. hence this also resets
   * BLOCK_as_string.
   */
  {

    BLOCK_as_string.clear();
    write_as_int_flag = given_write_as_int_flag;

  }


  inline void
  zero_index_BLOCK::specific_reset()
  // this drops all held data apart from the name & comments.
  {

    entry_value = CppSLHA_global::really_wrong_value;
    entry_recorded = false;

    scale_value = CppSLHA_global::really_wrong_value;
    scale_recorded = false;

    has_recorded_flag = false;

    if( has_extra_BLOCK )
      {

        delete extra_BLOCK;
        extra_BLOCK = NULL;  // just to be safe.
        has_extra_BLOCK = false;

      }

  }


  inline SLHA_BLOCK*
  zero_index_BLOCK::make_extra_BLOCK()
  // this makes a blank copy of this BLOCK.
  {

    extra_BLOCK = new zero_index_BLOCK( BLOCK_name,
                                        name_comment_string );

    has_extra_BLOCK = true;

    return extra_BLOCK;

  }


  inline bool
  zero_index_BLOCK::has_entry()
  const
  // this checks to see if the entry exists in this set & puts its value into
  // the provided memory address for a double.
  {

    return entry_recorded;

  }

  inline bool
  zero_index_BLOCK::has_entry( double* const double_for_entry )
  const
  // this checks to see if the entry exists in this set & puts its value into
  // the provided memory address for a double.
  {

    if( entry_recorded )
      {

        *double_for_entry = entry_value;

      }

    return entry_recorded;

  }

  inline double
  zero_index_BLOCK::get_entry()
  const
  // this returns the entry or returns CppSLHA_global::really_wrong_value.
  {

    return entry_value;

  }

  inline void
  zero_index_BLOCK::set_entry( double const input_entry )
  // this sets the entry & records that it was set, if it had not been set.
  {

    entry_value = input_entry;
    entry_recorded = true;
    BLOCK_as_string.clear();  // the string has to be redone now.

  }



  inline void
  one_index_BLOCK::specific_reset()
  // this drops all held data apart from the name & comments.
  {

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


  inline SLHA_BLOCK*
  one_index_BLOCK::make_extra_BLOCK()
  // this makes a blank copy of this BLOCK.
  {

    extra_BLOCK = new one_index_BLOCK( BLOCK_name,
                                       name_comment_string );

    has_extra_BLOCK = true;

    return extra_BLOCK;

  }



  inline void
  sparse_one_index_BLOCK::specific_reset()
  // this drops all held data apart from the name & comments.
  {

    entries.reset();

    has_recorded_flag = false;

    if( has_extra_BLOCK )
      {

        delete extra_BLOCK;
        extra_BLOCK = NULL;  // just to be safe.
        has_extra_BLOCK = false;

      }

  }


  inline SLHA_BLOCK*
  sparse_one_index_BLOCK::make_extra_BLOCK()
  // this makes a blank copy of this BLOCK.
  {

    extra_BLOCK = new sparse_one_index_BLOCK( BLOCK_name,
                                              name_comment_string );

    has_extra_BLOCK = true;

    return extra_BLOCK;

  }


  inline bool
  sparse_one_index_BLOCK::has_entry( int const requested_element,
                                     double* const double_for_entry )
  const
  // this checks to see if the entry exists in this set & puts its value into
  // the provided memory address for a double.
  {

    return entries.has_entry( requested_element,
                              double_for_entry );

  }

  inline double
  sparse_one_index_BLOCK::get_entry( int const requested_element )
  const
  // this either returns the value associated with requested_element or
  // CppSLHA_global::really_wrong_value if there is none.
  {

    return entries.get_entry( requested_element );

  }

  inline void
  sparse_one_index_BLOCK::set_entry( int const requested_element,
                                     double const input_entry )
  // this either overwrites the value at entry requested_element, or adds it to
  // the number_pair_set if it didn't have it.
  {

    entries.set_entry( requested_element,
                       input_entry );
    BLOCK_as_string.clear();  // the string has to be redone now.

  }



  inline void
  sparse_one_index_string_BLOCK::specific_reset()
  // this drops all held data apart from the name & comments.
  {

    delete entries;
    entries = new one_index_comment_organizer();

    has_recorded_flag = false;

    if( has_extra_BLOCK )
      {

        delete extra_BLOCK;
        extra_BLOCK = NULL;  // just to be safe.
        has_extra_BLOCK = false;

      }

  }


  inline SLHA_BLOCK*
  sparse_one_index_string_BLOCK::make_extra_BLOCK()
  // this makes a blank copy of this BLOCK.
  {

    extra_BLOCK = new sparse_one_index_string_BLOCK( BLOCK_name,
                                                     name_comment_string );

    has_extra_BLOCK = true;

    return extra_BLOCK;

  }


  inline bool
  sparse_one_index_string_BLOCK::has_string( int const requested_element,
                                          std::string* const string_for_entry )
  const
  // this looks for the entry for the requested indices, & if it is found,
  // writes it to the provided string, while returning true.
  {

    return entries->has_comment( requested_element,
                                 string_for_entry );

  }

  inline void
  sparse_one_index_string_BLOCK::set_entry( int const requested_element,
                                            std::string const input_entry )
  // this either overwrites the value at entry requested_element, or adds it to
  // the set if it didn't have it.
  {

    entries->set_comment( requested_element,
                          &input_entry );
    indices.set_entry( requested_element,
                       1.0 );
    BLOCK_as_string.clear();  // the string has to be redone now.

  }



  inline void
  multiple_same_one_index_string_BLOCK::specific_reset()
  // this drops all held data apart from the name & comments.
  {

    delete entries;
    entries = new two_index_comment_organizer();

    has_recorded_flag = false;

    if( has_extra_BLOCK )
      {

        delete extra_BLOCK;
        extra_BLOCK = NULL;  // just to be safe.
        has_extra_BLOCK = false;

      }

  }


  inline SLHA_BLOCK*
  multiple_same_one_index_string_BLOCK::make_extra_BLOCK()
  // this makes a blank copy of this BLOCK.
  {

    extra_BLOCK = new multiple_same_one_index_string_BLOCK( BLOCK_name,
                                                         name_comment_string );

    has_extra_BLOCK = true;

    return extra_BLOCK;

  }


  inline bool
  multiple_same_one_index_string_BLOCK::has_string(
                                                   int const requested_element,
                                          std::string* const string_for_entry )
  const
  // this looks for the entry for the requested indices, & if it is found,
  // writes it to the provided string, while returning true.
  {

    return entries->has_comment( requested_element,
                                 1,
                                 string_for_entry );

  }

  inline bool
  multiple_same_one_index_string_BLOCK::has_string(
                                             int const requested_first_element,
                                            int const requested_second_element,
                                          std::string* const string_for_entry )
  const
  // this looks for the entry for the requested indices, & if it is found,
  // writes it to the provided string, while returning true.
  {

    return entries->has_comment( requested_first_element,
                                 requested_second_element,
                                 string_for_entry );

  }


  inline void
  multiple_same_one_index_string_BLOCK::set_entry(
                                             int const requested_first_element,
                                            int const requested_second_element,
                                                std::string const input_entry )
  // this either overwrites the value at entry requested_element, or adds it to
  // the set if it didn't have it.
  {

    entries->set_comment( requested_first_element,
                          requested_second_element,
                          &input_entry );
    BLOCK_as_string.clear();  // the string has to be redone now.

  }



  inline SLHA_BLOCK*
  two_index_BLOCK::make_extra_BLOCK()
  // this makes a blank copy of this BLOCK.
  {

    extra_BLOCK = new two_index_BLOCK( BLOCK_name,
                                       name_comment_string );

    has_extra_BLOCK = true;

    return extra_BLOCK;

  }



  inline SLHA_BLOCK*
  three_index_BLOCK::make_extra_BLOCK()
  // this makes a blank copy of this BLOCK.
  {

    extra_BLOCK
    = new three_index_BLOCK( BLOCK_name,
                             name_comment_string );

    has_extra_BLOCK = true;

    return extra_BLOCK;

  }



  inline bool
  MASS_BLOCK::has_entry( int const requested_element,
                         double* const double_for_entry )
  const
  // this checks to see if the entry exists in this set & puts its value into
  // the provided memory address for a double.
  {

    particle_property_set const*
    sought_particle
    = reference_spectrum->inspect_particle_property_set( requested_element );

    if( NULL != sought_particle )
      {

        *double_for_entry = sought_particle->get_mass();

        return true;

      }
    else
      {

        return false;

      }

  }


  inline SLHA_BLOCK*
  MASS_BLOCK::make_extra_BLOCK()
  // this makes a blank copy of this BLOCK.
  {

    std::cout
    << std::endl
    << "CppSLHA::warning! unfortunately, this version of CppSLHA can NOT"
    << " properly manage multiple MASS BLOCKs at different scales. the latest"
    << " MASS BLOCK to be read overwrites any previous masses, in this"
    << " version.";
    std::cout << std::endl;

    extra_BLOCK
    = new MASS_BLOCK( reference_spectrum );

    has_extra_BLOCK = true;

    return extra_BLOCK;

  }


  inline void
  MASS_BLOCK::specific_reset()
  // this drops all held data apart from the name & comments.
  {

    // does nothing in this case.

  }


  inline double
  MASS_BLOCK::get_entry( int const requested_element )
  const
  {

    return reference_spectrum->get_PDG_coded_mass( requested_element );

  }

  inline void
  MASS_BLOCK::set_entry( int const requested_element,
                         double const input_entry )
  {

    reference_spectrum->record_PDG_coded_mass( requested_element,
                                               input_entry );
    BLOCK_as_string.clear();  // the string has to be redone now.

  }



  inline void
  SPheno_input_BLOCK::specific_reset()
  // this drops all held data apart from the name & comments.
  {

    entries.clear();

    has_recorded_flag = false;

    if( has_extra_BLOCK )
      {

        delete extra_BLOCK;
        extra_BLOCK = NULL;  // just to be safe.
        has_extra_BLOCK = false;

      }

  }


  inline SLHA_BLOCK*
  SPheno_input_BLOCK::make_extra_BLOCK()
  // this makes a blank copy of this BLOCK.
  {

    extra_BLOCK = new SPheno_input_BLOCK( BLOCK_name,
                                          name_comment_string );

    has_extra_BLOCK = true;

    return extra_BLOCK;

  }



  inline int
  SPheno_cross_section_line::get_number_of_final_state_particles()
  const
  {

    return number_of_final_state_particles;

  }

  inline int
  SPheno_cross_section_line::get_first_particle_code()
  const
  {

    return first_particle_code;

  }

  inline int
  SPheno_cross_section_line::get_second_particle_code()
  const
  {

    return second_particle_code;

  }

  inline int
  SPheno_cross_section_line::get_third_particle_code()
  const
  {

    return third_particle_code;

  }

  inline double
  SPheno_cross_section_line::get_cross_section()
  const
  {

    return cross_section;

  }



  inline int
  SPheno_cross_section_set::get_first_incoming_particle_code()
  const
  {

    return first_incoming_particle_code;

  }

  inline int
  SPheno_cross_section_set::get_second_incoming_particle_code()
  const
  {

    return second_incoming_particle_code;

  }

  inline double
  SPheno_cross_section_set::get_center_of_momentum_energy()
  const
  {

    return center_of_momentum_energy;

  }

  inline double
  SPheno_cross_section_set::get_first_incoming_particle_polarization()
  const
  {

    return first_incoming_particle_polarization;

  }

  inline double
  SPheno_cross_section_set::get_second_incoming_particle_polarization()
  const
  {

    return second_incoming_particle_polarization;

  }

  inline int
  SPheno_cross_section_set::get_include_ISR_flag()
  const
  {

    return include_ISR_flag;

  }

  inline std::vector< SPheno_cross_section_line* >*
  SPheno_cross_section_set::get_cross_section_lines()
  {

    return &cross_sections;

  }



  inline void
  SPheno_cross_sections_BLOCK::specific_reset()
  // this drops all held data apart from the name & comments.
  {

    has_recorded_flag = false;

    for( std::vector< SPheno_cross_section_set* >::iterator
           entry_iterator = entries.begin();
         entries.end() > entry_iterator;
         ++entry_iterator )
      {

        delete *entry_iterator;
        *entry_iterator = NULL;  // just to be safe.

      }

    entries.clear();

  }


  inline SLHA_BLOCK*
  SPheno_cross_sections_BLOCK::make_extra_BLOCK()
  // this makes a blank copy of this BLOCK.
  {

    extra_BLOCK = new SPheno_cross_sections_BLOCK( BLOCK_name,
                                                   name_comment_string,
                                                   reference_spectrum );

    has_extra_BLOCK = true;

    return extra_BLOCK;

  }


  inline std::vector< SPheno_cross_section_set* >*
  SPheno_cross_sections_BLOCK::get_SPheno_cross_section_vector()
  {

    return &entries;

  }

} // end of CppSLHA namespace.

#endif /* CPPSLHA_BLOCK_HPP_ */
