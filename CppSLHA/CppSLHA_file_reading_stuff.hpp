/*
 * CppSLHA_file_reading_stuff.hpp
 *
 *  Created on: 18 Jan 2011
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2011 Ben O'Leary
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

#ifndef CPPSLHA_FILE_READING_STUFF_HPP_
#define CPPSLHA_FILE_READING_STUFF_HPP_

#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

namespace CppSLHA
{

  // this class reads in a file & returns lines which are not empty but for
  // whitespace &/or a comment, using '#' as the commenting character.
  class hash_commented_file_reader
  {

  public:

    hash_commented_file_reader( std::string const* const input_file_name,
                                int const given_max_lines_to_read,
                                bool const lines_should_be_stored_or_not )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~hash_commented_file_reader()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    bool
    file_is_open()
    const
    // this returns true if the file to be read was successfully opened.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    bool
    read_next_line()
    /* this reads in the next line into line_being_read, then returns false if
     * the line is empty or just a comment, true otherwise. it also updates
     * whether the file is to continue to be read.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    bool
    read_to_next_valid_line()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    bool
    read_line( std::string* const line_to_fill )
    // this reads in the next valid line & stores it in the given string.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    copy_last_read_line( std::string* const line_to_fill )
    const
    // this copies the last valid read line into the given string.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    std::vector< std::string* > const*
    get_stored_lines()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    bool
    was_successful()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    std::ifstream input_file;
    // this is the file being read in.

    std::string line_being_read;
    // the string used for holding the line read.
    std::string word_being_read;
    // a string to be the word being read.
    char character_being_read;
    // a character for checking for the commenting character #.
    int position_of_character_in_line;
    // a number to track how far along a line we have gone looking for a #.
    int line_count;
    // a counter for the number of lines read.
    int const max_line_count;
    // the given maximum number of lines to be read.
    bool storing_lines_flag;
    // this sets whether lines read in are held in memory or not.
    std::vector< std::string* > stored_lines;
    // this holds the lines read in if they are to be stored.
    bool finished_reading_flag;
    // whether there are any more lines to read or not (accounting for
    // max_line_count).
    bool successful_reading_flag;
    // this is to note whether the whole file, or as much as was allowed, was
    // read or not.

  };



  // inline functions:

  inline bool
  hash_commented_file_reader::file_is_open()
  const
  // this returns true if the file to be read was successfully opened.
  {

    return input_file.is_open();

  }

  inline bool
  hash_commented_file_reader::read_to_next_valid_line()
  {

    if( file_is_open() )
      {

        while( !(read_next_line())
               &&
               !finished_reading_flag )
          {

            /* just keep reading in lines until either a valid line is found,
             * or we run out of file to read. when the loop is done,
             * finished_reading_flag is true if we didn't find a new valid
             * line, so we should return false, or it is false if we did, so we
             * should return true.
             */

          }

      }

    return !finished_reading_flag;

  }

  inline bool
  hash_commented_file_reader::read_line( std::string* const line_to_fill )
  // this reads in the next valid line & stores it in the given string.
  {

    if( read_to_next_valid_line() )
      // if we find a new valid line...
      {

        copy_last_read_line( line_to_fill );
        return true;

      }
    else
      {

        return false;

      }

  }

  inline void
  hash_commented_file_reader::copy_last_read_line(
                                              std::string* const line_to_fill )
  const
  // this copies the last valid read line into the given string.
  {

    if( NULL != line_to_fill )
      // if we don't have a definitely invalid pointer...
      {

        // take the given string & assign it the read line, stripped of
        // leading whitespace:
        line_to_fill->assign( line_being_read,
                              position_of_character_in_line,
                              ( line_being_read.size()
                                - position_of_character_in_line ) );

      }

  }

  inline std::vector< std::string* > const*
  hash_commented_file_reader::get_stored_lines()
  const
  {

    return &stored_lines;

  }


  inline bool
  hash_commented_file_reader::was_successful()
  const
  {

    return successful_reading_flag;

  }

} // end of CppSLHA namespace.

#endif /* CPPSLHA_FILE_READING_STUFF_HPP_ */
