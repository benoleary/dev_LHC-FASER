/*
 * LHC-FASER_Light.cpp
 *
 *  Created on: 31 Aug 2012
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
 *      LHC-FASER_Light.cpp
 *      LHC-FASER_base_electroweak_cascade_stuff.hpp
 *      LHC-FASER_base_electroweak_cascade_stuff.cpp
 *      LHC-FASER_base_kinematics_stuff.hpp
 *      LHC-FASER_base_kinematics_stuff.cpp
 *      LHC-FASER_base_lepton_distribution_stuff.hpp
 *      LHC-FASER_base_lepton_distribution_stuff.cpp
 *      LHC-FASER_charged_electroweak_cascade_stuff.hpp
 *      LHC-FASER_charged_electroweak_cascade_stuff.cpp
 *      LHC-FASER_cross-section_stuff.hpp
 *      LHC-FASER_cross-section_stuff.cpp
 *      LHC-FASER_derived_lepton_distributions.hpp
 *      LHC-FASER_derived_lepton_distributions.cpp
 *      LHC-FASER_electroweak_cascade_collection_stuff.hpp
 *      LHC-FASER_electroweak_cascade_collection_stuff.cpp
 *      LHC-FASER_full_cascade_stuff.hpp
 *      LHC-FASER_full_cascade_stuff.cpp
 *      LHC-FASER_global_stuff.hpp
 *      LHC-FASER_global_stuff.cpp
 *      LHC-FASER_input_handling_stuff.hpp
 *      LHC-FASER_input_handling_stuff.cpp
 *      LHC-FASER_jet_kinematics_stuff.hpp
 *      LHC-FASER_jet_kinematics_stuff.cpp
 *      LHC-FASER_lepton_kinematics_stuff.hpp
 *      LHC-FASER_lepton_kinematics_stuff.cpp
 *      LHC-FASER_neutral_electroweak_cascade_stuff.hpp
 *      LHC-FASER_neutral_electroweak_cascade_stuff.cpp
 *      LHC-FASER_signal_calculator_stuff.hpp
 *      LHC-FASER_signal_calculator_stuff.cpp
 *      LHC-FASER_signal_data_collection_stuff.hpp
 *      LHC-FASER_signal_data_collection_stuff.cpp
 *      LHC-FASER_sparticle_decay_stuff.hpp
 *      LHC-FASER_sparticle_decay_stuff.cpp
 *      LHC-FASER_template_classes.hpp
 *      and README.LHC-FASER.txt which describes the package.
 *
 *      LHC-FASER also requires CppSLHA. It should be found in a subdirectory
 *      included with this package.
 *
 *      LHC-FASER also requires grids of lookup values. These should also be
 *      found in a subdirectory included with this package.
 */

/* this program is to print out the cross-sections for a given SLHA spectrum on
 * a channel-by-channel basis, in fb. it takes a single argument, which is the
 * name of the SLHA file.
 */


// includes:

#include <fstream>
#include "LHC-FASER.hpp"


// the main program:

int main( int argumentCount,
          char* argumentStrings[] )
{
  LHC_FASER::lhcFaserLight onceOffLhcFaserLight( "./grids/",
                                                 "fb",
                                                 true );
  std::string slhaFilename( "" );
  std::string appendArgument( "" );
  if( 3 == argumentCount )
    // if the arguments should have been the SLHA filename & "--append"...
  {
    appendArgument.assign( argumentStrings[ 1 ] );
    if( 0 != appendArgument.compare( "--append" ) )
    {
      slhaFilename.assign( appendArgument );
      appendArgument.assign( argumentStrings[ 2 ] );
    }
    if( 0 != appendArgument.compare( "--append" ) )
    {
      std::cout
      << std::endl
      << "error! if 2 arguments are given, 1 of them should be the name of"
      << " the SLHA file and the other should be --append (which tells"
      << " LHC-FASER_Light to append its results to the SLHA file instead of"
      << " printing to the console)!";
      std::cout << std::endl;  // let the user know the format.
      return EXIT_FAILURE;
    }
  }
  else if( 2 == argumentCount )
    // if the argument should have been just the SLHA filename...
  {
    slhaFilename.assign( argumentStrings[ 1 ] );
  }
  else
  {
    std::cout
    << std::endl
    << "error! at least 1 argument must be given (the name of the SLHA file"
    << " with the spectrum)! if 2 arguments are given, 1 of them should be the"
    << " name of the SLHA file and the other should be --append (which tells"
    << " LHC-FASER_Light to append its results to the SLHA file instead of"
    << " printing to the console)!";
    std::cout << std::endl;  // let the user know the format.
    return EXIT_FAILURE;
  }

  if( appendArgument.empty() )
  {
    std::cout
    << std::endl
    << onceOffLhcFaserLight.fullResultsForNewSlha( slhaFilename );
    std::cout
    << "# in the format [particle 1 PDG code]  [particle 2 PDG code]"
    << "   [7 TeV cross-section in fb]   [14 TeV cross-section in fb]";
    std::cout << std::endl;
    std::cout << std::endl;
  }
  else
  {
    std::string
    resultsString( onceOffLhcFaserLight.fullResultsForNewSlha( slhaFilename,
                                                               " " ) );
    std::ofstream slhaFile( slhaFilename.c_str(),
                            std::ios::app );
    slhaFile
    << std::endl
    << "BLOCK LHCFASERLIGHT"
    << std::endl
    << "# particle 1  particle 2       7 TeV sigma/fb      14 TeV sigma/fb"
    << std::endl
    << resultsString
    << std::endl;
    slhaFile.close();
  }

    // this was a triumph! I'm making a note here:
  return EXIT_SUCCESS;

}
