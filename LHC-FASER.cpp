/*
 * LHC-FASER.cpp
 *
 *  Created on: 10 Nov 2010
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

#include "LHC-FASER.hpp"

namespace LHC_FASER
{

  // default constructor assuming default path to grids & default SLHA
  // spectrum file:
  LHC_FASER_UI::LHC_FASER_UI( CppSLHA::CppSLHA0* const given_spectrum_data ) :
    own_CppSLHA( false )
  {

    initialize( given_spectrum_data,
                "./grids/",
                "fb",
                true );

  }

  // constructor with specified SLHA file:
  LHC_FASER_UI::LHC_FASER_UI( std::string const given_SLHA_file_name ) :
    own_CppSLHA( true )
  {

    spectrum_data = new CppSLHA::CppSLHA2( given_SLHA_file_name );

    initialize( spectrum_data,
                "./grids/",
                "fb",
                true );

  }

  // constructor with specified CppSLHA and user-defined path to grids:
  LHC_FASER_UI::LHC_FASER_UI( CppSLHA::CppSLHA0* const given_spectrum_data,
                              std::string const given_path_to_grids ) :
    own_CppSLHA( false )
  {

    initialize( given_spectrum_data,
                given_path_to_grids,
                "fb",
                true );

  }

  // constructor with specified SLHA file and user-defined path to grids:
  LHC_FASER_UI::LHC_FASER_UI( std::string const given_SLHA_file_name,
                              std::string const given_path_to_grids ) :
    own_CppSLHA( true )
  {

    spectrum_data = new CppSLHA::CppSLHA2( given_SLHA_file_name );

    initialize( spectrum_data,
                given_path_to_grids,
                "fb",
                true );

  }

  // constructor with specified CppSLHA & user-defined path to grids,
  // & whether returned event rates should be in nb, pb or fb:
  LHC_FASER_UI::LHC_FASER_UI( CppSLHA::CppSLHA0* const given_spectrum_data,
                              std::string const given_path_to_grids,
                              std::string const given_cross_section_unit ) :
    own_CppSLHA( false )
  {

    initialize( given_spectrum_data,
                given_path_to_grids,
                given_cross_section_unit,
                true );

  }

  // constructor with specified SLHA file & user-defined path to grids,
  // & whether returned event rates should be in nb, pb or fb:
  LHC_FASER_UI::LHC_FASER_UI( std::string const given_SLHA_file_name,
                              std::string const given_path_to_grids,
                              std::string const given_cross_section_unit ) :
    own_CppSLHA( true )
  {

    spectrum_data = new CppSLHA::CppSLHA2( given_SLHA_file_name );

    initialize( spectrum_data,
                given_path_to_grids,
                given_cross_section_unit,
                true );

  }

  /* constructor with specified CppSLHA & user-defined path to grids,
   * whether returned event rates should be in nb, pb or fb& whether to use
   * LO or NLO:
   */
  LHC_FASER_UI::LHC_FASER_UI( CppSLHA::CppSLHA0* const given_spectrum_data,
                              std::string const given_path_to_grids,
                              std::string const given_cross_section_unit,
                              bool const given_NLO_to_be_used ) :
    own_CppSLHA( false )
  {

    initialize( given_spectrum_data,
                given_path_to_grids,
                given_cross_section_unit,
                given_NLO_to_be_used );

  }

  /* constructor with specified SLHA file & user-defined path to grids,
   * whether returned event rates should be in nb, pb or fb& whether to use
   * LO or NLO:
   */
  LHC_FASER_UI::LHC_FASER_UI( std::string const given_SLHA_file_name,
                              std::string const given_path_to_grids,
                              std::string const given_cross_section_unit,
                              bool const given_NLO_to_be_used ) :
    own_CppSLHA( true )
  {

    spectrum_data = new CppSLHA::CppSLHA2( given_SLHA_file_name );

    initialize( spectrum_data,
                given_path_to_grids,
                given_cross_section_unit,
                given_NLO_to_be_used );

  }

  LHC_FASER_UI::~LHC_FASER_UI()
  {

    if( own_CppSLHA )
      {

        delete spectrum_data;

      }

    for( std::vector< signalHandler* >::iterator
           deletion_iterator = signal_set.begin();
         signal_set.end() > deletion_iterator;
         deletion_iterator++ )
      {

        delete *deletion_iterator;

      }

    delete shortcut;
    delete input_handler_object;
    delete cross_section_handler_object;
    delete kinematics_handler_object;
    delete cascade_handler_object;
    delete readier_instance;

  }


  void
  LHC_FASER_UI::initialize( CppSLHA::CppSLHA0* const given_spectrum_data,
                            std::string const given_path_to_grids,
                            std::string const given_cross_section_unit,
                            bool const given_NLO_to_be_used )
  // this is used by all the constructors to do most of the construction.
  {

    spectrum_data = given_spectrum_data;

    path_to_grids.assign( given_path_to_grids );

    NLO_to_be_used = given_NLO_to_be_used;

    if( 0 == given_cross_section_unit.compare( "pb" ) )
      {

        cross_section_unit_factor = 1.0;

      }
    else if( 0 == given_cross_section_unit.compare( "fb" ) )
      {

        cross_section_unit_factor = 1000.0;

      }
    else if( 0 == given_cross_section_unit.compare( "nb" ) )
      {

        cross_section_unit_factor = 0.001;

      }
    else
      {

        std::cout
        << "LHC-FASER::error! given_unit has to be pb/fb/nb (all lowercase"
        << " characters). Carrying on, assuming fb.";
        std::cout << std::endl;

        cross_section_unit_factor = 1000.0;

      }

    readier_instance = new readierForNewPoint();

    input_handler_object = new input_handler( spectrum_data,
                                        spectrum_data->get_particle_spectrum(),
                                              path_to_grids,
                                              readier_instance );

    cross_section_handler_object
    = new crossSectionHandler( input_handler_object );
    kinematics_handler_object = new kinematics_handler( input_handler_object );
    cascade_handler_object = new cascade_handler( input_handler_object );

    shortcut = new signal_shortcuts( input_handler_object,
                                     cross_section_handler_object,
                                     kinematics_handler_object,
                                     cascade_handler_object );

  }

}  // end of LHC_FASER namespace.


