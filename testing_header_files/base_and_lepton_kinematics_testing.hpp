/*
 * base_and_lepton_kinematics_testing.hpp
 *
 *  Created on: Oct 30, 2011
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *              Jonas Lindert (jonas.lindert@googlemail.com)
 *              Carsten Robens (carsten.robens@gmx.de)
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

#ifndef BASE_AND_LEPTON_KINEMATICS_TESTING_HPP_
#define BASE_AND_LEPTON_KINEMATICS_TESTING_HPP_

#include "basic_testing_stuff.hpp"
#include "input_testing.hpp"
//#include "BOL_useful_stuff/BOL_useful_stuff.hpp"
#include "../LHC-FASER_lepton_kinematics_stuff.hpp"

namespace LHC_FASER
{
  class baseAndLeptonKinematicsTesting
  {
  public:
    baseAndLeptonKinematicsTesting( basicStuffTesting* basicStuff,
                                    inputTesting* inputTester ) :
        basicStuff( basicStuff ),
        inputTester( inputTester ),
        testInputHandler( inputTester->getInput() )
    {
      // just an initialization list.
    }
    ~baseAndLeptonKinematicsTesting()
    {
      // does nothing.
    }

    void
    performTest()
    {
      std::string
      acceptanceGridPathToFile(
                         *(testInputHandler->inspectPathToKinematicsGrids()) );
      std::cout
      << std::endl
      << "path to kinematics grids = " << acceptanceGridPathToFile;
      std::cout << std::endl;
      acceptanceGridPathToFile.append(
                        "/7TeV/jets/Atlas4jMET/squark+gluino_acceptance.dat" );
      acceptanceGrid jetsTestAcceptance( &acceptanceGridPathToFile,
                                         testInputHandler );
      std::cout
      << std::endl << "jetsTestAcceptance.testFunction( "
      << jetsTestAcceptance.getLowestSquarkMass() << ", "
      << jetsTestAcceptance.getLowestGluinoMass() << ", 0, 5 ) = "
      << jetsTestAcceptance.testFunction(
                                      jetsTestAcceptance.getLowestSquarkMass(),
                                      jetsTestAcceptance.getLowestGluinoMass(),
                                          0,
                                          5 )
      << std::endl;
      std::cout
      << "jetsTestAcceptance.testFunction( 123.4, 567.8, 0, 5 ) = "
      << jetsTestAcceptance.testFunction( 123.4,
                                          567.8,
                                          0,
                                          5 )
      << std::endl;
      std::cout
      << "jetsTestAcceptance.testFunction( 234.5, 678.9, 0, 3 ) = "
      << jetsTestAcceptance.testFunction( 234.5,
                                          678.9,
                                          0,
                                          5 )
      << std::endl;
      std::cout
      << "jetsTestAcceptance.testFunction( 876.5, 432.1, 1, 5 ) = "
      << jetsTestAcceptance.testFunction( 876.5,
                                          432.1,
                                          1,
                                          5 )
      << std::endl;
      std::cout
      << "jetsTestAcceptance.testFunction( 501.0, 501.0, 1, 5 ) = "
      << jetsTestAcceptance.testFunction( 501.0,
                                          501.0,
                                          1,
                                          5 );
      std::cout << std::endl;
      acceptanceGridPathToFile.assign(
                         *(testInputHandler->inspectPathToKinematicsGrids()) );
      std::cout
      << std::endl
      << "path to kinematics grids = " << acceptanceGridPathToFile;
      std::cout << std::endl;
      acceptanceGridPathToFile.append(
                                "/7TeV/leptons/squark+gluino_acceptance.dat" );
      acceptanceGrid*
      leptonsTestAcceptance = new acceptanceGrid( &acceptanceGridPathToFile,
                                                  testInputHandler );
      std::cout
      << std::endl << "leptonsTestAcceptance->testFunction( "
      << leptonsTestAcceptance->getLowestSquarkMass() << ", "
      << leptonsTestAcceptance->getLowestGluinoMass() << ", 0, 0 ) = "
      << leptonsTestAcceptance->testFunction(
                                  leptonsTestAcceptance->getLowestSquarkMass(),
                                  leptonsTestAcceptance->getLowestGluinoMass(),
                                              0,
                                              0 )
      << std::endl
      << "leptonsTestAcceptance->testFunction( 123.4, 567.8, 0, 0 ) = "
      << leptonsTestAcceptance->testFunction( 123.4,
                                              567.8,
                                              0,
                                              0 )
      << std::endl
      << "leptonsTestAcceptance->testFunction( 234.5, 678.9, 0, 0 ) = "
      << leptonsTestAcceptance->testFunction( 234.5,
                                              678.9,
                                              0,
                                              0 )
      << std::endl
      << "leptonsTestAcceptance->testFunction( 876.5, 432.1, 3, 10 ) = "
      << leptonsTestAcceptance->testFunction( 876.5,
                                              432.1,
                                              3,
                                              10 )
      << std::endl
      << "leptonsTestAcceptance->testFunction( 501.0, 501.0, 3, 10 ) = "
      << leptonsTestAcceptance->testFunction( 501.0,
                                              501.0,
                                              3,
                                              10 );
      std::cout << std::endl;
      delete leptonsTestAcceptance;


      std::string
      acceptanceGridPathToLeptons(
                         *(testInputHandler->inspectPathToKinematicsGrids()) );
      std::cout
      << std::endl
      << "path to kinematics grids = " << acceptanceGridPathToLeptons;
      std::cout << std::endl;
      acceptanceGridPathToLeptons.append(
                            "/7TeV/leptons/squark+antisquark_acceptance.dat" );
      leptonAcceptanceFromSquarkGrid
      testLeptonAcceptanceGrid( &acceptanceGridPathToLeptons,
                                testInputHandler );

      double testEffectiveSquarkMass;
      double testPseudorapidityAcceptance;
      std::vector< double > testEnergyAcceptances;
      testLeptonAcceptanceGrid.interpolateAcceptances(
                            testInputHandler->getSdownR()->get_absolute_mass(),
                                             testInputHandler->getGluinoMass(),
                     testInputHandler->getNeutralinoTwo()->get_absolute_mass(),
                                                      &testEffectiveSquarkMass,
                                                 &testPseudorapidityAcceptance,
                                                      &testEnergyAcceptances );
      std::cout
      << std::endl
      << "after testLeptonAcceptanceGrid.interpolateAcceptances( "
      << testInputHandler->getSdownR()->get_absolute_mass() << ", "
      << testInputHandler->getGluinoMass() << ", "
      << testInputHandler->getNeutralinoTwo()->get_absolute_mass() << ", ...),"
      << " testEffectiveSquarkMass = " << testEffectiveSquarkMass
      << " testPseudorapidityAcceptance = " << testPseudorapidityAcceptance;
      std::cout << std::endl;

      /**BOL::twoDimensionalDataPlotter*
      testPlotter
      = new BOL::twoDimensionalDataPlotter( "/opt/local/bin/gnuplot",
                              "lepton_acceptance_from_SPS1a_sdownRS_neu2.eps");
      double xValue;
      double yValue;
      testPlotter->addPoint( 0.0,
                             0.0 );
      for( unsigned int binCounter( 0 );
           testEnergyAcceptances.size() > binCounter;
           ++binCounter )
      {
        xValue = ( 2.0 * ( (double)( binCounter + 1 ) ) );
        yValue = testEnergyAcceptances.at( binCounter );
        testPlotter->addPoint( xValue,
                               yValue );
      }
      testPlotter->plotData();
      std::cout
      << std::endl << "testPlotter->plotData() seems to have been successful.";
      std::cout << std::endl;
      delete testPlotter;**/

      leptonAcceptanceParameterSet
      testLeptonAcceptanceParameterSet( testInputHandler,
                                        &testLeptonAcceptanceGrid,
                                        testInputHandler->getSdownR(),
                                        testInputHandler->getNeutralinoTwo() );
      std::cout
      << std::endl
      << "testLeptonAcceptanceParameterSet.getEffectiveSquarkMass() = "
      << testLeptonAcceptanceParameterSet.getEffectiveSquarkMass();
      std::cout
      << std::endl
      << "testLeptonAcceptanceParameterSet.acceptanceAt( 6.4, 10.0 ) = "
      << testLeptonAcceptanceParameterSet.acceptanceAt( 6.4,
                                                        10.0 );
      std::cout
      << std::endl
      << "testLeptonAcceptanceParameterSet.acceptanceAt( 7.4, 12.3 ) = "
      << testLeptonAcceptanceParameterSet.acceptanceAt( 7.4,
                                                        12.3 );
      std::cout
      << std::endl
      << "testLeptonAcceptanceParameterSet.acceptanceAt( 1007.4, 12.3 ) = "
      << testLeptonAcceptanceParameterSet.acceptanceAt( 1007.4,
                                                        12.3 );
      std::cout << std::endl;

      leptonAcceptanceHandler
      testLeptonAcceptanceHandler( testInputHandler,
                            testInputHandler->inspectPathToKinematicsGrids() );
      leptonAcceptanceParameterSet*
      testLeptonAcceptanceParameterSetPointer(
                testLeptonAcceptanceHandler.getLeptonAcceptanceParameterSet( 7,
                                                testInputHandler->getScharmR(),
                                            testInputHandler->getCharginoOne(),
                                                                           2.0,
                                                                       9.0 ) );
      std::cout
      << std::endl
      << "testLeptonAcceptanceParameterSetPointer uses scharm_R rather than"
      << " sdown_L, and importantly that the grid's lepton cut was 9.0 GeV"
      << " rather than 10.0 GeV.";
      std::cout
      << std::endl
      << "testLeptonAcceptanceParameterSetPointer->getEffectiveSquarkMass() = "
      << testLeptonAcceptanceParameterSetPointer->getEffectiveSquarkMass();
      std::cout
      << std::endl
      << "testLeptonAcceptanceParameterSetPointer->acceptanceAt( 6.4, 10.0 )"
      << " = " << testLeptonAcceptanceParameterSetPointer->acceptanceAt( 6.4,
                                                                         10.0 )
      << std::endl
      << "testLeptonAcceptanceParameterSetPointer->acceptanceAt( 7.4, 12.3 )"
      << " = " << testLeptonAcceptanceParameterSetPointer->acceptanceAt( 7.4,
                                                                        12.3 );
      std::cout << std::endl;
    }


  protected:
    basicStuffTesting* basicStuff;
    inputTesting* inputTester;
    inputHandler const* const testInputHandler;
  };

}  // end of LHC_FASER namespace

#endif /* BASE_AND_LEPTON_KINEMATICS_TESTING_HPP_ */
