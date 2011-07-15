/*
 * CppSLHA_waiting_on_subprocess_executor.hpp
 *
 *  Created on: 20 Mar 2010
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

#ifndef CPPSLHA_WAITING_ON_SUBPROCESS_EXECUTOR_HPP_
#define CPPSLHA_WAITING_ON_SUBPROCESS_EXECUTOR_HPP_

// includes:

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <sys/types.h>
#include <sys/wait.h>
#include <signal.h>

namespace CppSLHA
{

  class waiting_on_subprocess_executor
  {

  public:

    waiting_on_subprocess_executor( std::string const given_executable_name,
                std::vector< std::string const* > const* const given_arguments,
                                    int const given_patience_ticks )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~waiting_on_subprocess_executor()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    set_executable_name( std::string const input_executable_name )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    set_patience_ticks( int const input_patience_ticks )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    set_arguments(
               std::vector< std::string const* > const* const given_arguments )
    /* this assigns memory to make argument_list point to an array of pointers
     * to chars which are stored internally, & thus can be passed to execv to
     * properly execute the subprocess with its arguments.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    bool
    fork_and_execv_and_wait()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    int process_id;
    // this is for checking whether the process is the parent or the forked
    // child subprocess.

    std::string subprocess_executable_name;
    // this is the name for the executable to be run by execv.

    char const** argument_list_as_array;
    /* this is to hold the arguments for the executable. it is a constant
     * pointer to a pointer to a char because execv is ancient & requires
     * arguments in an irritatingly silly format.
     */

    std::vector< std::string* > argument_list_as_strings;

    int patience_ticks;
    /* this is the number of milliseconds for which the parent process will
     * wait before losing patience & killing the subprocess, if the subprocess
     * hasn't finished.
     */

  };  // end of waiting_on_subprocess_executor class.





  // inline functions:

  inline void
  waiting_on_subprocess_executor::set_executable_name(
                                      std::string const input_executable_name )
  {

    subprocess_executable_name.assign( input_executable_name );

  }

  inline void
  waiting_on_subprocess_executor::set_patience_ticks(
                                               int const input_patience_ticks )
  {

    patience_ticks = input_patience_ticks;

  }

}  // end of CppSLHA namespace

#endif /* CPPSLHA_WAITING_ON_SUBPROCESS_EXECUTOR_HPP_ */
