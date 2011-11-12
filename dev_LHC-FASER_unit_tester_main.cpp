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
 *      LHC-FASER_full_cascade_stuff.hpp
 *      LHC-FASER_full_cascade_stuff.cpp
 *      LHC-FASER_jet_kinematics_stuff.hpp
 *      LHC-FASER_jet_kinematics_stuff.cpp
 *
 *
 * still to test:
 *      LHC-FASER_signal_data_collection_stuff.hpp
 *      LHC-FASER_signal_data_collection_stuff.cpp
 *      LHC-FASER_signal_calculator_stuff.hpp
 *      LHC-FASER_signal_calculator_stuff.cpp
 *      LHC-FASER.hpp
 *      LHC-FASER.cpp
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

// future includes:
//#include "testing_header_files/lhcFaser_testing.hpp"


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
   * (works)*//**
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
   * (needs work)*//**/
  LHC_FASER::signalTesting signalTester( &basicStuffTester,
                                         &inputTester );
  signalTester.performTest();
  /**/



  /**


  LHC_FASER::crossSectionHandler test_cross_section_handler(
                                                         &testInputHandler );

  LHC_FASER::kinematics_handler test_kinematics_handler( &testInputHandler );

  LHC_FASER::cascade_handler test_cascade_handler( &testInputHandler );
  LHC_FASER::signalShortcuts test_signal_shortcuts( &testInputHandler,
                                                   &test_cross_section_handler,
                                                     &test_kinematics_handler,
                                                     &test_cascade_handler );

  LHC_FASER::signal_handler test_bad_signal( "testBadSignal",
                                             &test_signal_shortcuts );

  std::cout
  << std::endl
  << *(test_bad_signal.get_name()) << " has signal value "
  << test_bad_signal.getValue();
  std::cout << std::endl;

  LHC_FASER::signal_handler
  cross_section_breakdown_test( "sigmaBreakdownTest",
                                &test_signal_shortcuts );

  std::cout
  << std::endl
  << *(cross_section_breakdown_test.get_name()) << " has signal value "
  << cross_section_breakdown_test.getValue();
  std::cout << std::endl;

  LHC_FASER::signal_handler test_Atlas_signal( "Atlas4jMET0l7TeV",
                                               &test_signal_shortcuts );

  std::cout
  << std::endl
  << *(test_Atlas_signal.get_name()) << " constructed without crashing.";
  std::cout
  << std::endl;
  std::cout
  << std::endl
  << *(test_Atlas_signal.get_name()) << " has signal value "
  << test_Atlas_signal.getValue();
  std::cout << std::endl;

  slhaData.read_file( "SPS2_spectrum.out" );
  std::cout << std::endl;
  std::cout
  << std::endl << "without updating for new point, "
  << *(cross_section_breakdown_test.get_name()) << " has signal value "
  << cross_section_breakdown_test.getValue()
  << std::endl
  << *(test_Atlas_signal.get_name()) << " has signal value "
  << test_Atlas_signal.getValue();
  std::cout << std::endl;

  testReadier.readyObserversForNewPoint();
  std::cout
  << std::endl << "after updating for new point, "
  << *(cross_section_breakdown_test.get_name()) << " has signal value "
  << cross_section_breakdown_test.getValue()
  << std::endl
  << *(test_Atlas_signal.get_name()) << " has signal value "
  << test_Atlas_signal.getValue();
  std::cout << std::endl;


  std::cout
  << std::endl
  << "now trying alternating SPS1a & SPS2 4 times.";
  std::cout << std::endl;

  timeval starting_time;
  timeval ending_time;
  double elapsed_time;
  for( int reload_count = 1;
       4 >= reload_count;
       ++reload_count )
    {

      gettimeofday( &starting_time,
                    NULL );

      std::cout
      << std::endl << "count = " << reload_count;

      slhaData.read_file( "SPS1a_spectrum.out" );
      testReadier.readyObserversForNewPoint();
      std::cout
      << std::endl << "SPS1a: "
      << *(test_Atlas_signal.get_name()) << " has signal value "
      << test_Atlas_signal.getValue();
      std::cout << std::endl;

      slhaData.read_file( "SPS2_spectrum.out" );
      gettimeofday( &ending_time,
                    NULL );
      elapsed_time
      = ( (double)( ending_time.tv_sec - starting_time.tv_sec )
          + 0.000001 * (double)( ending_time.tv_usec
                                 - starting_time.tv_usec ) );
      std::cout
      << std::endl
      << "took " << seconds_taken << " seconds.";
      std::cout << std::endl;

      gettimeofday( &starting_time,
                    NULL );

      testReadier.readyObserversForNewPoint();
      std::cout
      << std::endl << "SPS2: "
      << *(test_Atlas_signal.get_name()) << " has signal value "
      << test_Atlas_signal.getValue();
      std::cout << std::endl;

      gettimeofday( &ending_time,
                    NULL );
      elapsed_time
      = ( (double)( ending_time.tv_sec - starting_time.tv_sec )
          + 0.000001 * (double)( ending_time.tv_usec
                                 - starting_time.tv_usec ) );
      std::cout
      << std::endl
      << "took " << elapsed_time << " seconds.";
      std::cout << std::endl;

    }


   **/


  /* finally, testing the lhcFaser itself: (not working)*//**

  LHC_FASER::lhcFaser SPS1a_UI( "SPS1a_spectrum.out" );
  LHC_FASER::lhcFaser SPS2_UI( "SPS2_spectrum.out" );

  SPS1a_UI.addSignal( "sigmaBreakdownTest" );
  SPS1a_UI.addSignal( "Atlas4jMET0l7TeV" );
  SPS2_UI.addSignal( "badSignalName" );
  SPS2_UI.addSignal( "Atlas4jMET0l7TeV" );

  LHC_FASER::signal_handler* sigma_SPS1a
  = SPS1a_UI.getSignal( "sigmaBreakdownTest" );
  LHC_FASER::signal_handler* Atlas_SPS1a
  = SPS1a_UI.getSignal( "Atlas4jMET0l7TeV" );
  LHC_FASER::signal_handler* Atlas_SPS2
  = SPS2_UI.getSignal( "Atlas4jMET0l7TeV" );
  LHC_FASER::signal_handler* bad_SPS2
  = SPS2_UI.getSignal( "badSignalName" );


  std::cout << std::endl;
  std::cout << std::endl << "all 4 test signals from test UIs have pointers.";
  std::cout << std::endl;

  SPS1a_UI.updateForNewSlha();
  std::cout << std::endl;
  std::cout << std::endl << "SPS1a_UI.updateForNewSlha() successful.";
  std::cout << std::endl;
  SPS2_UI.updateForNewSlha();
  std::cout << std::endl;
  std::cout << std::endl << "SPS2_UI.updateForNewSlha() successful.";
  std::cout << std::endl;
  std::cout
  << std::endl << "after updating for new point, "
  << *(sigma_SPS1a->get_name()) << " (sigma_SPS1a) has signal value "
  << sigma_SPS1a->getValue();
  std::cout << std::endl;
  std::cout
  << *(bad_SPS2->get_name()) << " (bad_SPS2) has signal value "
  << bad_SPS2->getValue();
  std::cout << std::endl;
  std::cout
  << *(Atlas_SPS1a->get_name()) << " (Atlas_SPS1a) has signal value "
  << Atlas_SPS1a->getValue();
  std::cout << std::endl;
  std::cout
  << *(Atlas_SPS2->get_name()) << " (Atlas_SPS2) has signal value "
  << Atlas_SPS2->getValue();
  std::cout << std::endl;

  **/


  /* testing adding new signals & how long it takes to load & calculate: *//**

  timeval start_time;
  timeval end_time;
  double seconds_taken;

  std::vector< LHC_FASER::signalHandler* > signals;
  std::string signal_name;

  //std::string grid_path( "/home/bol/BOL_work/projects/Eclipse_workspace/");
  //grid_path.append( "LHC-FASER_unit_tester/testing/grids/" );
  std::string grid_path( "./grids/");

  system( "cp LM1_spectrum.out test_spectrum.out" );

  gettimeofday( &start_time,
                NULL );
  std::cout
  << std::endl << "constructing test_UI";
  LHC_FASER::lhcFaser test_UI( "test_spectrum.out",
                                   grid_path,
                                   "fb" );
  gettimeofday( &end_time,
                NULL );
  seconds_taken
  = ( (double)( end_time.tv_sec - start_time.tv_sec )
      + 0.000001 * (double)( end_time.tv_usec - start_time.tv_usec ) );
  std::cout
  << " took " << seconds_taken << " seconds.";
  std::cout << std::endl;

  gettimeofday( &start_time,
                NULL );
  signal_name.assign( "Atlas4jMET0l7TeV" );
  std::cout
  << std::endl << "adding " << signal_name;
  test_UI.addSignal( signal_name );
  signals.push_back( test_UI.getSignal( signal_name ) );
  gettimeofday( &end_time,
                NULL );
  seconds_taken
  = ( (double)( end_time.tv_sec - start_time.tv_sec )
      + 0.000001 * (double)( end_time.tv_usec - start_time.tv_usec ) );
  std::cout
  << " took " << seconds_taken << " seconds.";
  std::cout << std::endl;

  gettimeofday( &start_time,
                NULL );
  signal_name.assign( "Atlas3jMET1l7TeV" );
  std::cout
  << std::endl << "adding " << signal_name;
  test_UI.addSignal( signal_name );
  signals.push_back( test_UI.getSignal( signal_name ) );
  gettimeofday( &end_time,
                NULL );
  seconds_taken
  = ( (double)( end_time.tv_sec - start_time.tv_sec )
      + 0.000001 * (double)( end_time.tv_usec - start_time.tv_usec ) );
  std::cout
  << " took " << seconds_taken << " seconds.";
  std::cout << std::endl;


  std::vector< std::string* > spectrum_files;
  spectrum_files.push_back( new std::string( "SPS1a_spectrum.out" ) );
  spectrum_files.push_back( new std::string( "SPS1b_spectrum.out" ) );
  spectrum_files.push_back( new std::string( "SPS2_spectrum.out" ) );
  spectrum_files.push_back( new std::string( "SPS3_spectrum.out" ) );
  //spectrum_files.push_back( new std::string( "SPS4_spectrum.out" ) );
  // SPS4 is borked.
  spectrum_files.push_back( new std::string( "SPS5_spectrum.out" ) );
  spectrum_files.push_back( new std::string( "SPS6_spectrum.out" ) );
  spectrum_files.push_back( new std::string( "SPS7_spectrum.out" ) );
  spectrum_files.push_back( new std::string( "SPS8_spectrum.out" ) );
  spectrum_files.push_back( new std::string( "SPS9_spectrum.out" ) );
  spectrum_files.push_back( new std::string( "LM1_spectrum.out" ) );
  spectrum_files.push_back( new std::string( "SU4_spectrum.out" ) );

  std::string system_call_string;
  for( int reload_count = 1;
       recalculation_amount >= reload_count;
       ++reload_count )
    {

      std::cout
      << std::endl << "count = " << reload_count;

      for( std::vector< std::string* >::const_iterator
           spectrum_iterator = spectrum_files.begin();
           spectrum_files.end() > spectrum_iterator;
           ++spectrum_iterator )
        {

          system_call_string.assign( "cp " );
          system_call_string.append( *(*spectrum_iterator) );
          system_call_string.append( " test_spectrum.out" );
          gettimeofday( &start_time,
                        NULL );
          std::cout
          << std::endl << "copying " << *(*spectrum_iterator);
          system( system_call_string.c_str() );
          gettimeofday( &end_time,
                        NULL );
          seconds_taken
          = ( (double)( end_time.tv_sec - start_time.tv_sec )
              + 0.000001 * (double)( end_time.tv_usec - start_time.tv_usec ) );
          std::cout
          << " took " << seconds_taken << " seconds.";
          std::cout << std::endl;

          gettimeofday( &start_time,
                        NULL );
          std::cout
          << std::endl << "reading it in & updating";
          test_UI.updateForNewSlha();
          gettimeofday( &end_time,
                        NULL );
          seconds_taken
          = ( (double)( end_time.tv_sec - start_time.tv_sec )
              + 0.000001 * (double)( end_time.tv_usec - start_time.tv_usec ) );
          std::cout
          << " took " << seconds_taken << " seconds.";
          std::cout << std::endl;

          for( std::vector< LHC_FASER::signalHandler* >::iterator
               signal_iterator = signals.begin();
               signals.end() > signal_iterator;
               ++signal_iterator )
            {

              gettimeofday( &start_time,
                            NULL );
              std::cout
              << std::endl << "calculating the point's value for "
              << *((*signal_iterator)->get_name()) << ": ";
              std::cout << (*signal_iterator)->getValue();
              gettimeofday( &end_time,
                            NULL );
              seconds_taken
              = ( (double)( end_time.tv_sec - start_time.tv_sec )
                  + 0.000001 * (double)( end_time.tv_usec
                                         - start_time.tv_usec ) );
              std::cout
              << "; took " << seconds_taken << " seconds.";
              std::cout << std::endl;

            }

        }

    }

  for( std::vector< std::string* >::const_iterator
       deletion_iterator = spectrum_files.begin();
       spectrum_files.end() > deletion_iterator;
       ++deletion_iterator )
    {

      delete *deletion_iterator;

    }

  **/


  std::cout << std::endl;
  std::cout << std::endl << "ended successfully, I hope";
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;

  return EXIT_SUCCESS;

}
