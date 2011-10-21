/*
 * CppSLHA_global_defines.hpp
 *
 *  Created on: 13 Oct 2010
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

#ifndef CPPSLHA_GLOBAL_STUFF_HPP_
#define CPPSLHA_GLOBAL_STUFF_HPP_

#include <stdlib.h>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <list>
#include <map>


// namespace:
namespace CppSLHA
{

  class CppSLHA_global
  {

  public:

    static void
    transform_to_uppercase( std::string* string_to_transform )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    static std::string
    SLHA_double( double const given_double )
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
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    static std::string
    SLHA_int( int const given_int,
              int const given_size )
    /* this prepares a string that is the ASCII version of an int, in the
     * form necessary for an SLHA file
     * (Fortran Iw, where w is an integer, which is given as given_size).
     * the string is given_size characters long.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    static std::string
    SLHA_int_or_double( double const given_double,
                        bool const write_as_int )
    // this prepares a string either using SLHA_double if write_as_int is false
    // or just a normally formatted int if given true.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    // unrecorded acceptanceValues are returned as really_wrong_value, & I think that
    // NaN is the best value to return to show up errors:
    static double const really_wrong_value;
    static std::string const really_wrong_value_string;

    // I am allowing lots of parameters which I have not yet thought of, but
    // there are limits for a reasonable number of parameters:
    static int const silly_amount_of_parameters;


  //protected:
    // nothing

  };


  inline std::string
  CppSLHA_global::SLHA_int_or_double( double const given_double,
                                      bool const write_as_int )
  // this prepares a string either using SLHA_double if write_as_int is false
  // or just a normally formatted int if given true.
  {

    if( write_as_int )
      {

        std::stringstream return_string_builder( "" );
        return_string_builder << (int)given_double;
        return return_string_builder.str();

      }
    else
      {

        return SLHA_double( given_double );

      }

  }

} // end of CppSLHA namespace.

#endif /* CPPSLHA_GLOBAL_STUFF_HPP_ */
