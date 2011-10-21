/*
 * CppSLHA_global_stuff.cpp
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

#include "CppSLHA_global_stuff.hpp"

namespace CppSLHA
{

  /* here is where all the defining for the static consts happens, in this
   * source file because at least 1 compiler tries to include the definitions
   * multiple times if they are put outside the class but in the header file,
   * regardless of whether they are within the #ifndef bit.
   */

  // unrecorded acceptanceValues are returned as really_wrong_value, & I think that
  // NaN is the best value to return to show up errors:
  double const CppSLHA_global::really_wrong_value = NAN;
  std::string const
  CppSLHA_global::really_wrong_value_string( "NaN (\"Not a Number\")" );

  int const CppSLHA_global::silly_amount_of_parameters = 100;


  void
  CppSLHA_global::transform_to_uppercase( std::string* string_to_transform )
  {

    for( unsigned int char_counter = 0;
         string_to_transform->size() > char_counter;
         ++char_counter )
      // go through each character in the string:
      {

        // if it's a lowercase character, replace it with its uppercase:
        if( 'a' == (*string_to_transform)[ char_counter ] )
          {

            (*string_to_transform)[ char_counter ] = 'A';

          }
        else if( 'b' == (*string_to_transform)[ char_counter ] )
          {

            (*string_to_transform)[ char_counter ] = 'B';

          }
        else if( 'c' == (*string_to_transform)[ char_counter ] )
          {

            (*string_to_transform)[ char_counter ] = 'C';

          }
        else if( 'd' == (*string_to_transform)[ char_counter ] )
          {

            (*string_to_transform)[ char_counter ] = 'D';

          }
        else if( 'e' == (*string_to_transform)[ char_counter ] )
          {

            (*string_to_transform)[ char_counter ] = 'E';

          }
        else if( 'f' == (*string_to_transform)[ char_counter ] )
          {

            (*string_to_transform)[ char_counter ] = 'F';

          }
        else if( 'g' == (*string_to_transform)[ char_counter ] )
          {

            (*string_to_transform)[ char_counter ] = 'G';

          }
        else if( 'h' == (*string_to_transform)[ char_counter ] )
          {

            (*string_to_transform)[ char_counter ] = 'H';

          }
        else if( 'i' == (*string_to_transform)[ char_counter ] )
          {

            (*string_to_transform)[ char_counter ] = 'I';

          }
        else if( 'j' == (*string_to_transform)[ char_counter ] )
          {

            (*string_to_transform)[ char_counter ] = 'J';

          }
        else if( 'k' == (*string_to_transform)[ char_counter ] )
          {

            (*string_to_transform)[ char_counter ] = 'K';

          }
        else if( 'l' == (*string_to_transform)[ char_counter ] )
          {

            (*string_to_transform)[ char_counter ] = 'L';

          }
        else if( 'm' == (*string_to_transform)[ char_counter ] )
          {

            (*string_to_transform)[ char_counter ] = 'M';

          }
        else if( 'n' == (*string_to_transform)[ char_counter ] )
          {

            (*string_to_transform)[ char_counter ] = 'N';

          }
        else if( 'o' == (*string_to_transform)[ char_counter ] )
          {

            (*string_to_transform)[ char_counter ] = 'O';

          }
        else if( 'p' == (*string_to_transform)[ char_counter ] )
          {

            (*string_to_transform)[ char_counter ] = 'P';

          }
        else if( 'q' == (*string_to_transform)[ char_counter ] )
          {

            (*string_to_transform)[ char_counter ] = 'Q';

          }
        else if( 'r' == (*string_to_transform)[ char_counter ] )
          {

            (*string_to_transform)[ char_counter ] = 'R';

          }
        else if( 's' == (*string_to_transform)[ char_counter ] )
          {

            (*string_to_transform)[ char_counter ] = 'S';

          }
        else if( 't' == (*string_to_transform)[ char_counter ] )
          {

            (*string_to_transform)[ char_counter ] = 'T';

          }
        else if( 'u' == (*string_to_transform)[ char_counter ] )
          {

            (*string_to_transform)[ char_counter ] = 'U';

          }
        else if( 'v' == (*string_to_transform)[ char_counter ] )
          {

            (*string_to_transform)[ char_counter ] = 'V';

          }
        else if( 'w' == (*string_to_transform)[ char_counter ] )
          {

            (*string_to_transform)[ char_counter ] = 'W';

          }
        else if( 'x' == (*string_to_transform)[ char_counter ] )
          {

            (*string_to_transform)[ char_counter ] = 'X';

          }
        else if( 'y' == (*string_to_transform)[ char_counter ] )
          {

            (*string_to_transform)[ char_counter ] = 'Y';

          }

        else if( 'z' == (*string_to_transform)[ char_counter ] )
          {

            (*string_to_transform)[ char_counter ] = 'Z';

          }

      }

  }


  std::string
  CppSLHA_global::SLHA_double( double const given_double )
  /* this prepares a string that is the ASCII version of a double, in the
   * form necessary for an SLHA file (Fortran E16.8).
   * 1st character: either "-" for negative numbers, or a blank space for
   *                positive numbers,
   * 2nd character: the 1st digit,
   * 3rd character: the decimal point,
   * 4th-11th characters: the 1st 8 digits following the decimal point,
   * 12th character: E
   * 13th character: + or -, depending on the sign of the exponent,
   * 14th, 15th & 16th characters: the absolute value of the exponent, with
   * a preceding 0s to fill to 3 characters.
   * SPheno almost uses this format, but only prints 2 digits of the
   * exponent.
   * NaNs are returned as 1 space followed by NaN followed by 12 spaces.
   */
  {

    std::stringstream return_string_builder( "" );

    double formatted_mantissa = given_double;

    int formatted_exponent = 0;

    if( 0.0 > given_double )
      {

        formatted_mantissa = -given_double;

        return_string_builder << "-";

      }
    else
      {

        return_string_builder << " ";

      }

    if( 0.0 <= formatted_mantissa )
    // at this point any negative numbers will now be positive, so any that
    // fail this comparison should be NaN.
      {

        while( formatted_mantissa >= 10.0 )
          {

            formatted_mantissa *= 0.1;

            ++formatted_exponent;

          }

        while( ( formatted_mantissa < 1.0 )
               &&
               ( formatted_mantissa > 0.0 ) )
          {

            formatted_mantissa *= 10.0;

            --formatted_exponent;

          }

        // this should set the stream to format numbers in fixed-digit
        // notation, then set it to having 8 digits after the decimal point:
        return_string_builder.setf( std::ios::fixed,
                                    std::ios::floatfield );
        return_string_builder.precision( 8 );

        return_string_builder << formatted_mantissa;
        // then write the mantissa.

        return_string_builder << "E";

        if( formatted_exponent < 0 )
          {

            return_string_builder << "-";

            formatted_exponent = -formatted_exponent;

          }
        else
          {

            return_string_builder << "+";

          }
        if( formatted_exponent < 100 )
          {

            return_string_builder << "0";

          }
        if( formatted_exponent < 10 )
          {

            return_string_builder << "0";

          }

        return_string_builder << formatted_exponent;
        // since the exponent is an integer, it shouldn't get formatted in
        // scientific format.

      }
    else
    // if it failed the comparison, it should be a NaN.
      {

        return_string_builder << "NaN            ";

      }

    return return_string_builder.str();

  }

  std::string
  CppSLHA_global::SLHA_int( int const given_int,
                            int const given_size )
  /* this prepares a string that is the ASCII version of an int, in the
   * form necessary for an SLHA file
   * (Fortran Iw, where w is an integer, which is given as given_size).
   * the string is given_size characters long.
   */
  {

    std::stringstream return_string_builder( "" );

    int absolute_value = given_int;

    if( 0 > given_int )
      {

        absolute_value = -absolute_value;

      }

    int actual_size = 0;

    if( 0 == absolute_value )
      {

        actual_size = 1;

      }
    else
      {

        int ten_to_digits = 1;
        while( absolute_value >= ten_to_digits )
          {

            ten_to_digits *= 10;
            ++actual_size;

          }

      }

    if( given_size < ( actual_size + 1 ) )
    /* if the size needed ( the number of digits + 1 for the
     * space-instead-of-plus-sign or the minus sign ) is bigger than the
     * requested string size...
     */
      {

        return_string_builder
        << "CppSLHA::error! SLHA_int( " << given_int << ", " << given_size
        << " ) could not fit the integer into the given size!";
        // report an error message.

      }
    else
      {

        for( int space_inserter = 0;
             ( given_size - actual_size - 1 ) > space_inserter;
             ++space_inserter )
          {

            return_string_builder << " ";

          }

        if( 0 > given_int )
          {

            return_string_builder << "-";

          }
        else
          {

            return_string_builder << " ";

          }

      }

    return_string_builder << absolute_value;

    return return_string_builder.str();

  }

}  // end of CppSLHA namespace.
