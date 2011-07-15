/*
 * LHC-FASER.hpp
 *
 *  Created on: 06 Oct 2010
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
 *      LHC-FASER_lepton_distributions.hpp
 *      LHC-FASER_lepton_distributions.cpp
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

#ifndef LHC_FASER_HPP_
#define LHC_FASER_HPP_

/* this is a C++ class for calculating event rates for signals of the Minimal
 * Supersymmetric Standard Model at the LHC. it calculates only the signal
 * event rate, assuming that the background is known.
 */

/* this version of the code was written mainly by Ben O'Leary ("BOL",
 * benjamin.oleary@gmail.com), based on a setup decided upon by Michael Krämer
 * ("MK", mkraemer@physik.rwth-aachen.de), Jonas Lindert ("JL",
 * jonas.lindert@googlemail.com) and BOL. Quite a lot is a refactoring of
 * JL's code. the cross-section grids were created by JL, while the current
 * acceptance grids were written by Carsten Robens
 * ("CR", Carsten.Robens@rwth-aachen.de).
 */

/* the method of calculation is as laid out in the publication
 * "SUSY parameter determination at the LHC using cross sections and kinematic
 * edges" by Herbi Dreiner, Michael Krämer, Jonas Lindert and Ben O'Leary,
 * published in JHEP 1004 (2010) 109 available on the arXiv as number
 * 1003.2648 as well.
 * other publications *should* follow, I hope.
 */

// this particular file is the header file for the main class. there is an
// accompanying .cpp file with the functions.


// include the header file with useful global-ish stuff, as well as including
// other general headers such as iostream:
#include "LHC-FASER_global_stuff.hpp"
// include the classes for handling signals (hence also those for handling the
// data collections for the signals):
#include "LHC-FASER_signal_calculator_stuff.hpp"


namespace LHC_FASER
{

  /* this is a class for calculating MSSM signal event rates at the LHC.
   * it tries to be a Factory (creating appropriate subtypes of signals based
   * on input).
   */
  class LHC_FASER_UI
  {

  public:

    // I thought about the factory pattern to replace all these constructors,
    // but then I decided that it wouldn't have any advantage.

    // default constructor assuming default path to grids & default SLHA
    // spectrum file:
    LHC_FASER_UI( CppSLHA::CppSLHA0* const given_spectrum_data )
    // default constructor with specified CppSLHA.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    // constructor with specified SLHA file:
    LHC_FASER_UI( std::string const given_SLHA_file_name )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    // constructor with specified CppSLHA & user-defined path to grids:
    LHC_FASER_UI( CppSLHA::CppSLHA0* const given_spectrum_data,
                  std::string const given_path_to_grids )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    // constructor with specified SLHA file & user-defined path to grids:
    LHC_FASER_UI( std::string const given_SLHA_file_name,
                  std::string const given_path_to_grids )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    // constructor with specified CppSLHA & user-defined path to grids,
    // & whether returned event rates should be in nb, pb or fb:
    LHC_FASER_UI( CppSLHA::CppSLHA0* const given_spectrum_data,
                  std::string const given_path_to_grids,
                  std::string const given_cross_section_unit )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    // constructor with specified SLHA file & user-defined path to grids,
    // & whether returned event rates should be in nb, pb or fb:
    LHC_FASER_UI( std::string const given_SLHA_file_name,
                  std::string const given_path_to_grids,
                  std::string const given_cross_section_unit )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    /* constructor with specified CppSLHA & user-defined path to grids,
     * whether returned event rates should be in nb, pb or fb& whether to use
     * LO or NLO:
     */
    LHC_FASER_UI( CppSLHA::CppSLHA0* const given_spectrum_datas,
                  std::string const given_path_to_grids,
                  std::string const given_cross_section_unit,
                  bool const given_NLO_to_be_used )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    /* constructor with specified SLHA file & user-defined path to grids,
     * whether returned event rates should be in nb, pb or fb& whether to use
     * LO or NLO:
     */
    LHC_FASER_UI( std::string const given_SLHA_file_name,
                  std::string const given_path_to_grids,
                  std::string const given_cross_section_unit,
                  bool const given_NLO_to_be_used )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~LHC_FASER_UI()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    signalHandler*
    add_signal( std::string const signal_name )
    // this adds a new signal to the set of signals based on its name.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    signalHandler*
    get_signal( std::string const signal_name )
    // this returns the handler object for the requested signal name.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    update_for_updated_point()
    /* this assumes that the CppSLHA was updated & so sets each signal to be
     * recalculated next time its value is requested. I expect that I could do
     * this more elegantly with throwing exceptions, but I'll leave that for
     * another day.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    update_for_new_point()
    // this reads in the CppSLHA's target file & recalculates all required
    // signals.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    update_for_new_point( std::string const SLHA_file_name )
    // this reads in the new file & recalculates all required signals.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    CppSLHA::CppSLHA0* spectrum_data;
    // this holds the MSSM spectrum data.
    bool own_CppSLHA;
    // this notes whether the CppSLHA was supplied externally (in which case
    // the destructor does NOT delete the CppSLHA).
    std::string path_to_grids;
    // this is where the lookup tables live.
    double cross_section_unit_factor;
    // this is to allow for the user to specify event rates in fb, pb or nb.
    bool NLO_to_be_used;
    // this is to allow the user to use LO or NLO cross-sections.

    input_handler* input_handler_object;
    // this keeps const pointers to useful objects together for ease of passing
    // around & for neater code.
    cross_section_handler* cross_section_handler_object;
    // this holds the lookup tables for LHC colored sparticle pair production
    // cross-sections.
    kinematics_handler* kinematics_handler_object;
    // this holds the kinematic data lookup tables for colored sparticle
    // production at the LHC.
    cascade_handler* cascade_handler_object;
    // this handles the cascade decays of colored sparticles.
    std::vector< signalHandler* > signal_set;
    // this tracks the various signals of LHC supersymmetric events & their
    // rates.

    signal_shortcuts* shortcut;
    // this keeps const pointers to useful objects together for ease of passing
    // around & for neater code.

    readier_for_new_point* readier_instance;


    void
    initialize( CppSLHA::CppSLHA0* const given_spectrum_data,
                std::string const given_path_to_grids,
                std::string const given_cross_section_unit,
                bool const given_NLO_to_be_used )
    // this is used by all the constructors to do most of the construction.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };  // end of LHC_FASER_UI class.




  // inline functions:

  inline signalHandler*
  LHC_FASER_UI::add_signal( std::string const signal_name )
  // this adds a new signal to the set of signals based on its name.
  {

    signalHandler* return_pointer = new signalHandler( signal_name,
                                                         shortcut,
                                                   cross_section_unit_factor );

    signal_set.push_back( return_pointer );

    return return_pointer;

  }

  inline signalHandler*
  LHC_FASER_UI::get_signal( std::string const signal_name )
  // this returns the handler object for the requested signal name.
  {

    signalHandler* return_pointer = NULL;

    for( std::vector< signalHandler* >::iterator
         signal_iterator = signal_set.begin();
         signal_set.end() > signal_iterator;
         signal_iterator++ )
      // look through all the signals...
      {

        if( 0 == (*signal_iterator)->getName()->compare( signal_name ) )
          // if we find the requested signal...
          {

            return_pointer = *signal_iterator;

            signal_iterator = signal_set.end();
            // stop looking.

          }

      }

    return return_pointer;

  }


  inline void
  LHC_FASER_UI::update_for_updated_point()
  /* this assumes that the CppSLHA was updated & so sets each signal to be
   * recalculated next time its value is requested. I expect that I could do
   * this more elegantly with throwing exceptions, but I'll leave that for
   * another day.
   */
  {

    readier_instance->ready_for_new_point();

  }

  inline void
  LHC_FASER_UI::update_for_new_point()
  // this reads in the CppSLHA's target file & recalculates all required
  // signals.
  {

    spectrum_data->read_file();
    update_for_updated_point();

  }

  inline void
  LHC_FASER_UI::update_for_new_point( std::string const SLHA_file_name )
  // this reads in the new file & recalculates all required signals.
  {

    spectrum_data->read_file( SLHA_file_name );
    update_for_updated_point();

  }

}  // end of LHC_FASER namespace.

#endif /* LHC_FASER_HPP_ */
