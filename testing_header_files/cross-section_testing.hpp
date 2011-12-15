/*
 * cross-section_testing.hpp
 *
 *  Created on: Oct 30, 2011
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

#ifndef CROSS_SECTION_TESTING_HPP_
#define CROSS_SECTION_TESTING_HPP_

#include "input_testing.hpp"
#include "../LHC-FASER_cross-section_stuff.hpp"

namespace LHC_FASER
{
  class crossSectionTesting
  {
  public:
    crossSectionTesting( basicStuffTesting* basicStuff,
                         inputTesting* inputTester ) :
        basicStuff( basicStuff ),
        inputTester( inputTester ),
        testInputHandler( inputTester->getInput() )
    {
      // just an initialization list.
    }
    ~crossSectionTesting()
    {
      // does nothing.
    }

    void
    performTest()
    {
      std::cout
      << std::endl
      << "pathToGrids = "
      << *(testInputHandler->getPathToGrids())
      << std::endl;

      std::string
      pathToTestGrid( *(testInputHandler->inspectPathToCrossSectionGrids()) );
      pathToTestGrid.append( "/7TeV/gluino+gluino_LO.dat" );
      squareGrid testGrid( &pathToTestGrid,
                           basicStuff->getBluh(),
                           NULL,
                           testInputHandler );
      crossSectionTable
      gluinoGluinoTable( &testGrid,
                testInputHandler->getScoloredProductionCombinations()->at( 0 ),
                         1.0,
                         testInputHandler );
      std::cout
      << std::endl << "x: " << testGrid.getLowestX() << " - "
      << testGrid.getHighestX();
      std::cout
      << std::endl << "y: " << testGrid.getLowestY() << " - "
      << testGrid.getHighestY();
      std::cout
      << std::endl
      << "value at ( 200.0, 200.0 ) = " << testGrid.valueAt( 200.0,
                                                             200.0 );
      std::cout << std::endl;
      std::cout
      << std::endl
      << "value at ( 200.0, 205.0 ) = " << testGrid.valueAt( 200.0,
                                                             205.0 );
      std::cout << std::endl;
      std::cout
      << std::endl
      << "value at ( 205.0, 200.0 ) = " << testGrid.valueAt( 205.0,
                                                             200.0 );
      std::cout << std::endl;
      std::cout
      << std::endl
      << "value at ( 201.0, 201.0 ) = " << testGrid.valueAt( 201.0,
                                                             201.0 );
      std::cout << std::endl;
      std::cout
      << std::endl
      << "value at ( 210.0, 210.0 ) = " << testGrid.valueAt( 210.0,
                                                             210.0 );
      std::cout << std::endl;
      std::cout
      << std::endl
      << "value at ( 199.5, 201.0 ) = " << testGrid.valueAt( 199.5,
                                                             201.0 );
      std::cout << std::endl;
      std::cout
      << std::endl
      << "value at ( 199.5, 199.5 ) = " << testGrid.valueAt( 199.5,
                                                             199.5 );
      std::cout << std::endl;
      std::cout
      << std::endl
      << "value at ( 2000.0, 200.5 ) = " << testGrid.valueAt( 2000.0,
                                                              200.5 );
      std::cout << std::endl;
      std::cout << std::endl;
      std::cout
      << std::endl
      << "value at ( 2000.5, 200.5 ) = " << testGrid.valueAt( 2000.5,
                                                              200.5 );
      std::cout << std::endl;
      std::cout
      << std::endl
      << "value at ( 2020.0, 200.5 ) = " << testGrid.valueAt( 2020.0,
                                                              200.5 );
      std::cout << std::endl;
      std::cout
      << std::endl
      << "value at ( 2222.2, 200.5 ) = " << testGrid.valueAt( 2222.2,
                                                              200.5 );
      std::cout << std::endl;
      std::cout
      << std::endl
      << "value at ( 200.5, 2000.0 ) = " << testGrid.valueAt( 200.5,
                                                              2000.0 );
      std::cout << std::endl;
      std::cout
      << std::endl
      << "value at ( 200.5, 2020.0 ) = " << testGrid.valueAt( 200.5,
                                                              2020.0 );
      std::cout << std::endl;
      std::cout
      << std::endl
      << "value at ( 200.5, 2222.2 ) = " << testGrid.valueAt( 200.5,
                                                              2222.2 );
      std::cout << std::endl;
      std::cout
      << std::endl
      << "value at ( 1999.5, 1999.5 ) = " << testGrid.valueAt( 1999.5,
                                                               1999.5 );
      std::cout << std::endl;
      std::cout
      << std::endl
      << "value at ( 2000.0, 2000.0 ) = " << testGrid.valueAt( 2000.0,
                                                               2000.0 );
      std::cout << std::endl;
      std::cout
      << std::endl
      << "value at ( 2000.5, 2000.5 ) = " << testGrid.valueAt( 2000.5,
                                                               2000.5 );
      std::cout << std::endl;
      std::cout
      << std::endl
      << "value at ( 2020.0, 2020.0 ) = " << testGrid.valueAt( 2020.0,
                                                               2020.0 );
      std::cout << std::endl;
      std::cout
      << std::endl
      << "value at ( 2200.0, 2200.0 ) = " << testGrid.valueAt( 2200.0,
                                                               2200.0 );
      std::cout << std::endl;
      std::cout
      << std::endl
      << "value at ( 2200.5, 2200.2 ) = " << testGrid.valueAt( 2200.5,
                                                               2200.2 );
      std::cout << std::endl;
      std::cout
      << std::endl
      << "gluinoGluinoTable value = " << gluinoGluinoTable.getValue();

      pathToTestGrid.assign(
                       *(testInputHandler->inspectPathToCrossSectionGrids()) );
      pathToTestGrid.append( "/7TeV/gluino+gluino_K.dat" );
      squareGrid testKFactorGrid( &pathToTestGrid,
                                  basicStuff->getBluh(),
                                  NULL,
                                  testInputHandler );
      std::cout
      << std::endl << "K-factor value at ( 200.5, 200.5 ) = "
      << testKFactorGrid.valueAt( 200.5,
                                  200.5 );
      std::cout << std::endl;

      crossSectionHandler
      crossSectionHandlerTester( testInputHandler );
      crossSectionTableSet*
      testCrossSectionSet( crossSectionHandlerTester.getTableSet( 7 ) );
      std::cout
      << std::endl << "testCrossSectionSet = " << testCrossSectionSet;
      // N.B.: for the next bit to work, crossSectionHandler::getTable() has to
      // be moved back to public.
      crossSectionTable*
      sdGoSigma( testCrossSectionSet->getTable(
            testInputHandler->getScoloredProductionCombinations()->at( 1 ) ) );
      crossSectionTable*
      suGoSigma( testCrossSectionSet->getTable(
           testInputHandler->getScoloredProductionCombinations()->at( 10 ) ) );
      crossSectionTableSet*
      testCrossSectionSet2( crossSectionHandlerTester.getTableSet( 7 ) );
      std::cout
      << std::endl << "testCrossSectionSet2 = " << testCrossSectionSet2;
      crossSectionTable*
      sdGoSigma2( testCrossSectionSet2->getTable(
            testInputHandler->getScoloredProductionCombinations()->at( 1 ) ) );
      std::cout
      << std::endl << "no extra squareGrid loading?";

      pathToTestGrid.assign(
                       *(testInputHandler->inspectPathToCrossSectionGrids()) );
      pathToTestGrid.append( "/7TeV/sdown_L+gluino_LO.dat" );
      squareGrid testSdGoLo( &pathToTestGrid,
                                        basicStuff->getBluh(),
                                        NULL,
                                        testInputHandler );
      pathToTestGrid.assign(
                       *(testInputHandler->inspectPathToCrossSectionGrids()) );
      pathToTestGrid.append( "/7TeV/sdown_L+gluino_K.dat" );
      squareGrid testSdGoK( &pathToTestGrid,
                                       basicStuff->getBluh(),
                                       NULL,
                                       testInputHandler );
      double
      sdMass( testInputHandler->getScoloredProductionCombinations()->at( 1
                                  )->getFirstParticle()->get_absolute_mass() );
      double
      goMass( testInputHandler->getScoloredProductionCombinations()->at( 1
                                 )->getSecondParticle()->get_absolute_mass() );
      std::cout << std::endl;
      std::cout
      << std::endl << "sd mass = " << sdMass
      << std::endl << "go mass = " << goMass
      << std::endl << "testSdGoLo value = " << testSdGoLo.valueAt( sdMass,
                                                                   goMass )
      << std::endl << "testSdGoK value = " << testSdGoK.valueAt( sdMass,
                                                                 goMass )
      << std::endl << "hence NLO value should be "
      << ( testSdGoLo.valueAt( sdMass,
                               goMass )
           * testSdGoK.valueAt( sdMass,
                                goMass ) )
      << std::endl << "sdGoSigma = " << sdGoSigma
      << std::endl << "sdGoSigma value = " << sdGoSigma->getValue()
      << std::endl << "sdGoSigma2 = " << sdGoSigma2
      << std::endl << "sdGoSigma2 value = " << sdGoSigma->getValue()
      << std::endl << "suGoSigma = " << suGoSigma
      << std::endl << "suGoSigma value = " << suGoSigma->getValue();
    }


  protected:
    basicStuffTesting* basicStuff;
    inputTesting* inputTester;
    inputHandler const* const testInputHandler;
  };

}  // end of LHC_FASER namespace

#endif /* CROSS_SECTION_TESTING_HPP_ */
