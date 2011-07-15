/*
 * CppSLHA_file_reading_stuff.cpp
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

#include "CppSLHA_file_reading_stuff.hpp"

namespace CppSLHA
{

  hash_commented_file_reader::hash_commented_file_reader(
                                      std::string const* const input_file_name,
                                             int const given_max_lines_to_read,
                                   bool const lines_should_be_stored_or_not ) :
      line_being_read( "before input" ),
      word_being_read( "unread" ),
      character_being_read( 'x' ),
      position_of_character_in_line( 0 ),
      line_count( 0 ),
      max_line_count( given_max_lines_to_read ),
      storing_lines_flag( lines_should_be_stored_or_not ),
      finished_reading_flag( false ),
      successful_reading_flag( false )
  {

    input_file.open( input_file_name->c_str() );

    if( !(input_file.is_open()) )
      {

        finished_reading_flag = true;

        std::cout
        << std::endl
        << "CppSLHA::error! hash_commented_file_reader tried to open "
        << *input_file_name << " but could not!";
        std::cout << std::endl;

      }

    if( max_line_count <= line_count )
      {

        finished_reading_flag = true;

      }

  }

  hash_commented_file_reader::~hash_commented_file_reader()
  {

    if( input_file.is_open() )
      {

        input_file.close();

      }

    for( std::vector< std::string* >::iterator
           deletion_iterator = stored_lines.begin();
         stored_lines.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;
        *deletion_iterator = NULL;  // just to be safe.

      }

  }

  bool
  hash_commented_file_reader::read_next_line()
  /* this reads in the next line into line_being_read, then returns false if
   * the line is empty or just a comment, true otherwise. it also updates
   * whether the file is to continue to be read.
   */
  {

    bool return_bool = false;

    if( input_file.eof()
        ||
        ( max_line_count <= line_count ) )
      /* if we have not gone "too far" (hoping to stop corrupted files
       * from being read in forever, though any file missing EOF is
       * probably going to create much more serious problems anyway)...
       */
      {

        finished_reading_flag = true;
        // note that no more reading of this file should be attempted.
        successful_reading_flag = true;
        // also note that we got to EOF or the maximum allowed length of the
        // file, so the reading was "successful".

      }

    if( !finished_reading_flag )
      // if we can still read more of this file...
      {

        std::getline( input_file,
                      line_being_read );
        // take in a line (should only read up to a newline character, & not
        // include it) & assign it to line_being_read.

        // debugging:
        /*std::cout << std::endl << "debugging: read line:";
        std::cout << std::endl << line_being_read;
        std::cout << std::endl;
        */

        if( 0 != line_being_read.compare( "" ) )
          // ignoring blank lines... (string.compare() returns 0 if the
          // comparison is good.)
          {

            // check to see if the 1st non-whitespace character is a letter, a
            // number or a #:

            position_of_character_in_line = 0;

            character_being_read
                            = line_being_read[ position_of_character_in_line ];

            while( ' ' == character_being_read )
              {

                ++position_of_character_in_line;

                character_being_read
                            = line_being_read[ position_of_character_in_line ];

              }

            if( '#' != character_being_read )
              // if we have not found a comment mark before anything else but
              // whitespace...
              {

                // ... we have a valid line, so:

                return_bool = true;
                ++line_count;

                if( storing_lines_flag )
                  {

                    stored_lines.push_back(
                                          new std::string( line_being_read ) );

                  }

              }

          }  // end of if statement ignoring blank lines.

      }  // end of if it's valid to look for a new line.

    return return_bool;

  }

} // end of CppSLHA namespace.
