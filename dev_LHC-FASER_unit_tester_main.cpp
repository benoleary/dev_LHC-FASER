/*
 * dev_LHC-FASER_unit_tester_main.cpp
 *
 *  Created on: 21 Jan 2011
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

/* this program is to test various bits of LHC-FASER in isolation (or as much
 * isolation as practical).
 *
 * I should really test each bit as it is written, BUT the whole structure of
 * LHC-FASER is subject to changes, even major changes, leaving good chances
 * of testing units that do not find their way into even the 1st prototype
 * compilable version.
 *
 *
 * the current ordering for testing is:
 *
 * done so far:
 *      LHC-FASER_template_classes.hpp
 *      LHC-FASER_global_stuff.hpp
 *      LHC-FASER_global_stuff.cpp
 *      LHC-FASER_base_lepton_distribution_stuff.hpp
 *      LHC-FASER_base_lepton_distribution_stuff.cpp
 *      LHC-FASER_derived_lepton_distributions.hpp
 *      LHC-FASER_derived_lepton_distributions.cpp
 *      LHC-FASER_sparticle_decay_stuff.hpp
 *      LHC-FASER_sparticle_decay_stuff.cpp
 *      LHC-FASER_input_handling_stuff.hpp
 *      LHC-FASER_input_handling_stuff.cpp
 *      LHC-FASER_cross-section_stuff.hpp
 *      LHC-FASER_cross-section_stuff.cpp
 *      LHC-FASER_base_kinematics_stuff.hpp
 *      LHC-FASER_base_kinematics_stuff.cpp
 *      LHC-FASER_lepton_kinematics_stuff.hpp
 *      LHC-FASER_lepton_kinematics_stuff.cpp
 *      LHC-FASER_base_electroweak_cascade_stuff.hpp
 *      LHC-FASER_base_electroweak_cascade_stuff.cpp
 *      LHC-FASER_neutral_electroweak_cascade_stuff.hpp
 *      LHC-FASER_neutral_electroweak_cascade_stuff.cpp
 *      LHC-FASER_charged_electroweak_cascade_stuff.hpp
 *      LHC-FASER_charged_electroweak_cascade_stuff.cpp
 *      LHC-FASER_electroweak_cascade_collection_stuff.hpp
 *      LHC-FASER_electroweak_cascade_collection_stuff.cpp
 *      LHC-FASER_jet_kinematics_stuff.hpp
 *      LHC-FASER_jet_kinematics_stuff.cpp
 *      LHC-FASER_signal_data_collection_stuff.hpp
 *      LHC-FASER_signal_data_collection_stuff.cpp
 *      LHC-FASER_signal_calculator_stuff.hpp
 *      LHC-FASER_signal_calculator_stuff.cpp
 *      LHC-FASER_full_cascade_stuff.hpp
 *      LHC-FASER_full_cascade_stuff.cpp
 *      LHC-FASER.hpp
 *      LHC-FASER.cpp
 *
 *
 * still to test:
 *      nothing! we're done. well, the accuracy of the approximations has to be
 *      tested...
 */


// includes:

#include <sys/time.h>
#include "testing_header_files/basic_testing_stuff.hpp"
#include "testing_header_files/lepton_energy_distributions_testing.hpp"
#include "testing_header_files/input_testing.hpp"
#include "testing_header_files/cross-section_testing.hpp"
#include "testing_header_files/base_and_lepton_kinematics_testing.hpp"
#include "testing_header_files/base_electroweak_cascade_testing.hpp"
#include "testing_header_files/full_electroweak_cascade_testing.hpp"
#include "testing_header_files/full_cascade_testing.hpp"
#include "testing_header_files/jet_kinematics_testing.hpp"
#include "testing_header_files/signal_testing.hpp"
#include "testing_header_files/signal_timing_testing.hpp"
#include "testing_header_files/lhcFaser_testing.hpp"

// future includes:
// none, since we're done!



// the main program:

int main( int argumentCount,
          char* argumentStrings[] )
{
  int recalculationAmount( 10 );

  if( 2 != argumentCount )
    // if the input was definitely wrong...
    {
      std::cout
      << std::endl
      << "error! this program requires an integer as its argument, as the"
      << " number of times to re-calculate all the points. defaulting to 10.";
      std::cout << std::endl;  // let the user know the format.
    }
  else  // otherwise, do what is supposed to be done.
    {
      std::string argumentInterpretter( argumentStrings[ 1 ] );
      std::stringstream argumentStream( argumentInterpretter );
      argumentStream >> recalculationAmount;
    }

  LHC_FASER::basicStuffTesting basicStuffTester( "SPS1a_spectrum.out" );

  /* testing the lepton energy distributions in the squark rest frame:
   * (works) *//**
  LHC_FASER::leptonEnergyDistributionsTestingClass
  leptonEnergyDistributionsTester( &basicStuffTester );
  leptonEnergyDistributionsTester.performTest();
  **/


  /* testing the means of checking if something has already been prepared for a
   * new point: (works)*//**/
  LHC_FASER::inputTesting inputTester( &basicStuffTester );
  //inputTester.performTest();
  /**/


  /* testing the cross-section grid stuff:
   * (works)*//**
  LHC_FASER::crossSectionTesting crossSectionTester( &basicStuffTester,
                                                     &inputTester );
  crossSectionTester.performTest();
  **/


  /* testing the base and lepton kinematics grid stuff:
   * (works, but the testing needs BOL_useful_stuff, which I'm not including in
   * the repository)*//**
  LHC_FASER::baseAndLeptonKinematicsTesting
  baseAndLeptonKinematicsTester( &basicStuffTester,
                                 &inputTester );
  baseAndLeptonKinematicsTester.performTest();
  **/


  /* testing the basic electroweak cascade calculation stuff:
   * (works)*//**
  LHC_FASER::baseElectroweakCascadeTesting
  baseElectroweakCascadeTester( &basicStuffTester,
                                &inputTester );
  baseElectroweakCascadeTester.performTest();
  **/

  /* testing the electroweak cascade combination stuff:
   * (works)*//**
  LHC_FASER::fullElectroweakCascadeTesting
  fullElectroweakCascadeTester( &basicStuffTester,
                                &inputTester );
  fullElectroweakCascadeTester.performTest();
  **/


  /* testing the full cascade stuff:
   * (needs checking for new version)*//**
  LHC_FASER::fullCascadeTesting fullCascadeTester( &basicStuffTester,
                                                   &inputTester );
  fullCascadeTester.performTest();
  **/


  /* testing the jet kinematics grid stuff:
   * (works)*//**
  LHC_FASER::jetKinematicsTesting jetKinematicsTester( &basicStuffTester,
                                                       &inputTester );
  jetKinematicsTester.performTest();
  **/

	

  /* testing the signal stuff:
   * (works)*//**
  LHC_FASER::signalTesting signalTester( &basicStuffTester,
                                         &inputTester );
  signalTester.performTest();
  **/

  /* finally, testing the lhcFaser itself:
   * (works)*//**
  LHC_FASER::lhcFaserTesting lhcFaserTester( &basicStuffTester,
                                             &inputTester );
  lhcFaserTester.performTest();
  **/


  /* testing adding new signals & how long it takes to load & calculate:
   * (works) *//**/
  LHC_FASER::signalTimingTesting signalTimingTester( &basicStuffTester,
                                                     &inputTester );
  signalTimingTester.performTest( recalculationAmount );
  /**/


  std::cout << std::endl;
  std::cout << std::endl << "ended successfully, I hope";
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;

  return EXIT_SUCCESS;

}
