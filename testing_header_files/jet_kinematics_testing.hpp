/*
 * jet_kinematics_testing.hpp
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

#ifndef JET_KINEMATICS_TESTING_HPP_
#define JET_KINEMATICS_TESTING_HPP_

#include "basic_testing_stuff.hpp"
#include "input_testing.hpp"
#include "../LHC-FASER_jet_kinematics_stuff.hpp"

namespace LHC_FASER
{
  class jetKinematicsTesting
  {
  public:
    jetKinematicsTesting( basicStuffTesting* basicStuff,
                          inputTesting* inputTester ) :
      slhaData( basicStuff->getSlha() ),
      testReadier( basicStuff->getReadier() ),
      testInputHandler( inputTester->getInput() )
    {
      // just an initialization list.
    }
    ~jetKinematicsTesting()
    {
      // does nothing.
    }

    void
    performTest()
    {
      particlePointer gluinoPointer( testInputHandler->getGluino() );
      particlePointer sdownLPointer( testInputHandler->getSdownL() );
      particlePointer sbottomOnePointer( testInputHandler->getSbottomOne() );
      particlePointer sbottomTwoPointer( testInputHandler->getSbottomTwo() );
      particlePointer stopOnePointer( testInputHandler->getStopOne() );
      particlePointer
      neutralinoOnePointer( testInputHandler->getNeutralinoOne() );
      particlePointer
      neutralinoTwoPointer( testInputHandler->getNeutralinoTwo() );
      particlePointer
      charginoOnePointer( testInputHandler->getCharginoOne() );
      particlePointer wPlusPointer( testInputHandler->getWPlus() );
      particlePointer firstPointer( gluinoPointer );
      particlePointer secondPointer( sdownLPointer );
      squarkMassForGridDecider*
      testSquarkMassForGridDecider( new heavierThanGluinoSquarkMassForGrid(
                                                          testInputHandler ) );
      std::cout
      << std::endl
      << "heavierThanGluinoSquarkMassForGrid( sdown_L, gluino ) = "
      << (*testSquarkMassForGridDecider)( sdownLPointer,
                                          gluinoPointer );
      std::cout
      << std::endl
      << "heavierThanGluinoSquarkMassForGrid( gluino, sbottom_1 ) = "
      << (*testSquarkMassForGridDecider)( gluinoPointer,
                                          sbottomOnePointer );
      std::cout
      << std::endl
      << "heavierThanGluinoSquarkMassForGrid( sdown_L, sbottom_1 ) = "
      << (*testSquarkMassForGridDecider)( sdownLPointer,
                                          sbottomOnePointer );
      std::cout << std::endl;
      delete testSquarkMassForGridDecider;
      testSquarkMassForGridDecider = new firstMassForGrid( testInputHandler );
      std::cout
      << std::endl
      << "firstMassForGrid( " << *(firstPointer->get_name()) << " ("
      << firstPointer->get_absolute_mass() << "), "
      << *(secondPointer->get_name()) << " ("
      << secondPointer->get_absolute_mass() << ") ) = "
      << (*testSquarkMassForGridDecider)( firstPointer,
                                          secondPointer );
      firstPointer = sdownLPointer;
      secondPointer = gluinoPointer;
      std::cout
      << std::endl
      << "firstMassForGrid( " << *(firstPointer->get_name()) << " ("
      << firstPointer->get_absolute_mass() << "), "
      << *(secondPointer->get_name()) << " ("
      << secondPointer->get_absolute_mass() << ") ) = "
      << (*testSquarkMassForGridDecider)( firstPointer,
                                          secondPointer );
      firstPointer = sbottomOnePointer;
      secondPointer = sbottomOnePointer;
      std::cout
      << std::endl
      << "firstMassForGrid( " << *(firstPointer->get_name()) << " ("
      << firstPointer->get_absolute_mass() << "), "
      << *(secondPointer->get_name()) << " ("
      << secondPointer->get_absolute_mass() << ") ) = "
      << (*testSquarkMassForGridDecider)( firstPointer,
                                          secondPointer );
      firstPointer = sdownLPointer;
      std::cout
      << std::endl
      << "firstMassForGrid( " << *(firstPointer->get_name()) << " ("
      << firstPointer->get_absolute_mass() << "), "
      << *(secondPointer->get_name()) << " ("
      << secondPointer->get_absolute_mass() << ") ) = "
      << (*testSquarkMassForGridDecider)( firstPointer,
                                          secondPointer );
      std::cout << std::endl;
      delete testSquarkMassForGridDecider;
      testSquarkMassForGridDecider = new secondMassForGrid( testInputHandler );
      std::cout
      << std::endl
      << "secondMassForGrid( " << *(firstPointer->get_name()) << " ("
      << firstPointer->get_absolute_mass() << "), "
      << *(secondPointer->get_name()) << " ("
      << secondPointer->get_absolute_mass() << ") ) = "
      << (*testSquarkMassForGridDecider)( firstPointer,
                                          secondPointer );
      firstPointer = sdownLPointer;
      secondPointer = gluinoPointer;
      std::cout
      << std::endl
      << "secondMassForGrid( " << *(firstPointer->get_name()) << " ("
      << firstPointer->get_absolute_mass() << "), "
      << *(secondPointer->get_name()) << " ("
      << secondPointer->get_absolute_mass() << ") ) = "
      << (*testSquarkMassForGridDecider)( firstPointer,
                                          secondPointer );
      firstPointer = sbottomOnePointer;
      secondPointer = sbottomOnePointer;
      std::cout
      << std::endl
      << "secondMassForGrid( " << *(firstPointer->get_name()) << " ("
      << firstPointer->get_absolute_mass() << "), "
      << *(secondPointer->get_name()) << " ("
      << secondPointer->get_absolute_mass() << ") ) = "
      << (*testSquarkMassForGridDecider)( firstPointer,
                                          secondPointer );
      firstPointer = sdownLPointer;
      std::cout
      << std::endl
      << "secondMassForGrid( " << *(firstPointer->get_name()) << " ("
      << firstPointer->get_absolute_mass() << "), "
      << *(secondPointer->get_name()) << " ("
      << secondPointer->get_absolute_mass() << ") ) = "
      << (*testSquarkMassForGridDecider)( firstPointer,
                                          secondPointer );
      std::cout << std::endl;
      delete testSquarkMassForGridDecider;
      testSquarkMassForGridDecider
      = new averageSquarkMassForGrid( testInputHandler );
      std::cout
      << std::endl
      << "averageSquarkMassForGrid( " << *(firstPointer->get_name()) << " ("
      << firstPointer->get_absolute_mass() << "), "
      << *(secondPointer->get_name()) << " ("
      << secondPointer->get_absolute_mass() << ") ) = "
      << (*testSquarkMassForGridDecider)( firstPointer,
                                          secondPointer );
      firstPointer = sdownLPointer;
      secondPointer = gluinoPointer;
      std::cout
      << std::endl
      << "averageSquarkMassForGrid( " << *(firstPointer->get_name()) << " ("
      << firstPointer->get_absolute_mass() << "), "
      << *(secondPointer->get_name()) << " ("
      << secondPointer->get_absolute_mass() << ") ) = "
      << (*testSquarkMassForGridDecider)( firstPointer,
                                          secondPointer );
      firstPointer = sbottomOnePointer;
      secondPointer = sbottomOnePointer;
      std::cout
      << std::endl
      << "averageSquarkMassForGrid( " << *(firstPointer->get_name()) << " ("
      << firstPointer->get_absolute_mass() << "), "
      << *(secondPointer->get_name()) << " ("
      << secondPointer->get_absolute_mass() << ") ) = "
      << (*testSquarkMassForGridDecider)( firstPointer,
                                          secondPointer );
      firstPointer = sdownLPointer;
      std::cout
      << std::endl
      << "averageSquarkMassForGrid( " << *(firstPointer->get_name()) << " ("
      << firstPointer->get_absolute_mass() << "), "
      << *(secondPointer->get_name()) << " ("
      << secondPointer->get_absolute_mass() << ") ) = "
      << (*testSquarkMassForGridDecider)( firstPointer,
                                          secondPointer );
      std::cout << std::endl;
      delete testSquarkMassForGridDecider;

      std::string
      pathToGrids( *(testInputHandler->inspectPathToKinematicsGrids()) );
      pathToGrids.append(
                        "/7TeV/jets/Atlas4jMET/squark+gluino_acceptance.dat" );
      jetAcceptanceGrid*
      testJetAcceptanceGrid( new jetAcceptanceGrid( &pathToGrids,
                                                    testInputHandler ) );
      double squarkMass( 123.4 );
      double gluinoMass( 567.8 );
      double firstEwinoMass( 32.1 );
      double secondEwinoMass( 43.2 );
      int acceptanceColumn( 10 );
      std::cout
      << std::endl
      << "testJetAcceptanceGrid->getAcceptance( " << squarkMass << ", "
      << gluinoMass << ", " << firstEwinoMass << ", " << secondEwinoMass
      << ", " << acceptanceColumn << " ) = "
      << testJetAcceptanceGrid->getAcceptance( squarkMass,
                                               gluinoMass,
                                               firstEwinoMass,
                                               secondEwinoMass,
                                               acceptanceColumn );
      squarkMass = 987.6;
      gluinoMass = 876.5;
      firstEwinoMass = 500.0;
      secondEwinoMass = 333.1;
      std::cout
      << std::endl
      << "testJetAcceptanceGrid->getAcceptance( " << squarkMass << ", "
      << gluinoMass << ", " << firstEwinoMass << ", " << secondEwinoMass
      << ", " << acceptanceColumn << " ) = "
      << testJetAcceptanceGrid->getAcceptance( squarkMass,
                                               gluinoMass,
                                               firstEwinoMass,
                                               secondEwinoMass,
                                               acceptanceColumn );
      std::cout << std::endl;
      secondEwinoMass = 500.0;
      firstEwinoMass = 333.1;
      std::cout
      << std::endl
      << "testJetAcceptanceGrid->getAcceptance( " << squarkMass << ", "
      << gluinoMass << ", " << firstEwinoMass << ", " << secondEwinoMass
      << ", " << acceptanceColumn << " ) = "
      << testJetAcceptanceGrid->getAcceptance( squarkMass,
                                               gluinoMass,
                                               firstEwinoMass,
                                               secondEwinoMass,
                                               acceptanceColumn );
      std::cout << std::endl;
      gluinoMass = 2876.5;
      std::cout
      << std::endl
      << "testJetAcceptanceGrid->getAcceptance( " << squarkMass << ", "
      << gluinoMass << ", " << firstEwinoMass << ", " << secondEwinoMass
      << ", " << acceptanceColumn << " ) = "
      << testJetAcceptanceGrid->getAcceptance( squarkMass,
                                               gluinoMass,
                                               firstEwinoMass,
                                               secondEwinoMass,
                                               acceptanceColumn );
      std::cout << std::endl;
      gluinoMass = 876.5;
      firstEwinoMass = 1500.0;
      std::cout
      << std::endl
      << "testJetAcceptanceGrid->getAcceptance( " << squarkMass << ", "
      << gluinoMass << ", " << firstEwinoMass << ", " << secondEwinoMass
      << ", " << acceptanceColumn << " ) = "
      << testJetAcceptanceGrid->getAcceptance( squarkMass,
                                               gluinoMass,
                                               firstEwinoMass,
                                               secondEwinoMass,
                                               acceptanceColumn );
      std::cout << std::endl;
      secondEwinoMass = 1333.1;
      std::cout
      << std::endl
      << "testJetAcceptanceGrid->getAcceptance( " << squarkMass << ", "
      << gluinoMass << ", " << firstEwinoMass << ", " << secondEwinoMass
      << ", " << acceptanceColumn << " ) = "
      << testJetAcceptanceGrid->getAcceptance( squarkMass,
                                               gluinoMass,
                                               firstEwinoMass,
                                               secondEwinoMass,
                                               acceptanceColumn );
      std::cout << std::endl;
      squarkMass = 500.1;
      gluinoMass = 600.2;
      firstEwinoMass = 400.3;
      secondEwinoMass = 400.4;
      std::cout
      << std::endl
      << "testJetAcceptanceGrid->getAcceptance( " << squarkMass << ", "
      << gluinoMass << ", " << firstEwinoMass << ", " << secondEwinoMass
      << ", " << acceptanceColumn << " ) = "
      << testJetAcceptanceGrid->getAcceptance( squarkMass,
                                               gluinoMass,
                                               firstEwinoMass,
                                               secondEwinoMass,
                                               acceptanceColumn );
      std::cout << std::endl;
      firstEwinoMass = 200.3;
      secondEwinoMass = 400.4;
      std::cout
      << std::endl
      << "testJetAcceptanceGrid->getAcceptance( " << squarkMass << ", "
      << gluinoMass << ", " << firstEwinoMass << ", " << secondEwinoMass
      << ", " << acceptanceColumn << " ) = "
      << testJetAcceptanceGrid->getAcceptance( squarkMass,
                                               gluinoMass,
                                               firstEwinoMass,
                                               secondEwinoMass,
                                               acceptanceColumn );
      std::cout << std::endl;
      firstEwinoMass = 0.3;
      secondEwinoMass = 400.4;
      std::cout
      << std::endl
      << "testJetAcceptanceGrid->getAcceptance( " << squarkMass << ", "
      << gluinoMass << ", " << firstEwinoMass << ", " << secondEwinoMass
      << ", " << acceptanceColumn << " ) = "
      << testJetAcceptanceGrid->getAcceptance( squarkMass,
                                               gluinoMass,
                                               firstEwinoMass,
                                               secondEwinoMass,
                                               acceptanceColumn );
      std::cout << std::endl;
      firstEwinoMass = 0.0;
      std::cout
      << std::endl
      << "testJetAcceptanceGrid->getAcceptance( " << squarkMass << ", "
      << gluinoMass << ", " << firstEwinoMass << ", " << secondEwinoMass
      << ", " << acceptanceColumn << " ) = "
      << testJetAcceptanceGrid->getAcceptance( squarkMass,
                                               gluinoMass,
                                               firstEwinoMass,
                                               secondEwinoMass,
                                               acceptanceColumn );
      std::cout << std::endl;
      firstEwinoMass = 200.3;
      secondEwinoMass = 200.4;
      std::cout
      << std::endl
      << "testJetAcceptanceGrid->getAcceptance( " << squarkMass << ", "
      << gluinoMass << ", " << firstEwinoMass << ", " << secondEwinoMass
      << ", " << acceptanceColumn << " ) = "
      << testJetAcceptanceGrid->getAcceptance( squarkMass,
                                               gluinoMass,
                                               firstEwinoMass,
                                               secondEwinoMass,
                                               acceptanceColumn );
      std::cout << std::endl;
      firstEwinoMass = 0.3;
      secondEwinoMass = 200.4;
      std::cout
      << std::endl
      << "testJetAcceptanceGrid->getAcceptance( " << squarkMass << ", "
      << gluinoMass << ", " << firstEwinoMass << ", " << secondEwinoMass
      << ", " << acceptanceColumn << " ) = "
      << testJetAcceptanceGrid->getAcceptance( squarkMass,
                                               gluinoMass,
                                               firstEwinoMass,
                                               secondEwinoMass,
                                               acceptanceColumn );
      std::cout << std::endl;
      firstEwinoMass = 0.0;
      std::cout
      << std::endl
      << "testJetAcceptanceGrid->getAcceptance( " << squarkMass << ", "
      << gluinoMass << ", " << firstEwinoMass << ", " << secondEwinoMass
      << ", " << acceptanceColumn << " ) = "
      << testJetAcceptanceGrid->getAcceptance( squarkMass,
                                               gluinoMass,
                                               firstEwinoMass,
                                               secondEwinoMass,
                                               acceptanceColumn );
      std::cout << std::endl;
      firstEwinoMass = 0.3;
      secondEwinoMass = 0.4;
      std::cout
      << std::endl
      << "testJetAcceptanceGrid->getAcceptance( " << squarkMass << ", "
      << gluinoMass << ", " << firstEwinoMass << ", " << secondEwinoMass
      << ", " << acceptanceColumn << " ) = "
      << testJetAcceptanceGrid->getAcceptance( squarkMass,
                                               gluinoMass,
                                               firstEwinoMass,
                                               secondEwinoMass,
                                               acceptanceColumn );
      std::cout << std::endl;
      firstEwinoMass = 0.0;
      std::cout
      << std::endl
      << "testJetAcceptanceGrid->getAcceptance( " << squarkMass << ", "
      << gluinoMass << ", " << firstEwinoMass << ", " << secondEwinoMass
      << ", " << acceptanceColumn << " ) = "
      << testJetAcceptanceGrid->getAcceptance( squarkMass,
                                               gluinoMass,
                                               firstEwinoMass,
                                               secondEwinoMass,
                                               acceptanceColumn );
      std::cout << std::endl;
      secondEwinoMass = 0.0;
      std::cout
      << std::endl
      << "testJetAcceptanceGrid->getAcceptance( " << squarkMass << ", "
      << gluinoMass << ", " << firstEwinoMass << ", " << secondEwinoMass
      << ", " << acceptanceColumn << " ) = "
      << testJetAcceptanceGrid->getAcceptance( squarkMass,
                                               gluinoMass,
                                               firstEwinoMass,
                                               secondEwinoMass,
                                               acceptanceColumn );
      std::cout << std::endl;
      delete testJetAcceptanceGrid;
      std::cout
      << std::endl
      << "jetAcceptanceGrid class seems OK.";
      std::cout << std::endl;
      jetPlusMetAcceptanceHandler
      testJetPlusMetAcceptanceHandler( testInputHandler );
      std::string jetCutName( "Atlas4jMET" );
      jetAcceptanceTable*
      testJetAcceptanceTable(
               testJetPlusMetAcceptanceHandler.getJetPlusMetAcceptanceTable( 7,
                                                                   &jetCutName,
                                                                        10 ) );
      jetAcceptanceTable*
      testJetAcceptanceTable2(
               testJetPlusMetAcceptanceHandler.getJetPlusMetAcceptanceTable( 7,
                                                                   &jetCutName,
                                                                        10 ) );
      std::cout
      << std::endl
      << "testJetAcceptanceTable = " << testJetAcceptanceTable << ", "
      << "testJetAcceptanceTable2 = " << testJetAcceptanceTable2;
      std::cout << std::endl;


      electroweakCascadeHandler
      testElectroweakCascadeHandler( testInputHandler );
      electroweakCascadesForOneBeamEnergy*
      cascadesForSevenTev(
          testElectroweakCascadeHandler.getElectroweakCascadesForOneBeamEnergy(
                                                                         7 ) );
      electroweakCascadeSet*
      testEwinoCascadeSet( cascadesForSevenTev->getCascadeSet( sdownLPointer,
                                                      neutralinoTwoPointer ) );
      fullCascadeType::squarkDirectlyToElectroweakType::sdownType
      testSxFullCascade( testInputHandler,
                         sdownLPointer,
                         7,
                         testEwinoCascadeSet );
      testEwinoCascadeSet
      = cascadesForSevenTev->getCascadeSet( gluinoPointer,
                                            neutralinoTwoPointer );
      fullCascadeType::gluinoDirectlyToElectroweak
      testGxFullCascade( testInputHandler,
                         7,
                         testEwinoCascadeSet );
      testEwinoCascadeSet = cascadesForSevenTev->getCascadeSet( stopOnePointer,
                                                          charginoOnePointer );

      effectiveSquarkMassHolder*
      stopMinusW( testInputHandler->getSquarkMinusBosonEffectiveMass(
                                                                stopOnePointer,
                                                                  wPlusPointer,
                                                      neutralinoOnePointer ) );
      effectiveSquarkMassHolder*
      sbottomPlusW( testInputHandler->getSquarkPlusBosonEffectiveMass(
                                                             sbottomOnePointer,
                                                                  wPlusPointer,
                                                      neutralinoOnePointer ) );
      fullCascadeType::squarkDirectlyToElectroweakType::supType
      testStXOneFullCascade( testInputHandler,
                             stopOnePointer,
                             7,
                             testEwinoCascadeSet,
                            cascadesForSevenTev->getCascadeSet( stopOnePointer,
                                                            charginoOnePointer,
                                                                stopMinusW ),
                            cascadesForSevenTev->getCascadeSet( stopOnePointer,
                                                                wPlusPointer,
                                                             sbottomOnePointer,
                                                              sbottomPlusW ) );
      electroweakCascadeSet*
      testBosonCascadeSet( cascadesForSevenTev->getCascadeSet(
                                                             sbottomTwoPointer,
                                                               wPlusPointer,
                                                            stopOnePointer ) );
      fullCascadeType::squarkByBosonToCompound testSbsxFullCascade;
      testSbsxFullCascade.setProperties( sbottomTwoPointer,
                                         testBosonCascadeSet,
                                         &testStXOneFullCascade );
      fullCascadeType::gluinoOrNeutralinoToCompound testGjsxFullCascade;
      testGjsxFullCascade.setProperties( gluinoPointer,
                                         &testSxFullCascade,
                                         cascadesForSevenTev );
      fullCascadeType::gluinoOrNeutralinoToCompound testGjsbsxFullCascade;
      testGjsbsxFullCascade.setProperties( gluinoPointer,
                                           &testSbsxFullCascade,
                                           cascadesForSevenTev );
      signedParticleShortcutPair gluinoWithSdownL( gluinoPointer,
                                                   true,
                                                   sdownLPointer,
                                                   true );
      signedParticleShortcutPair sbottomOneWithGluino( gluinoPointer,
                                                       true,
                                                       sbottomOnePointer,
                                                       true );
      signedParticleShortcutPair gluinoPair( gluinoPointer,
                                             true,
                                             gluinoPointer,
                                             true );
      signedParticleShortcutPair
      sbottomTwoWithAntisbottomTwo( sbottomTwoPointer,
                                    true,
                                    sbottomTwoPointer,
                                    false );
      signedParticleShortcutPair sdownLWithSdownL( sdownLPointer,
                                                   true,
                                                   sdownLPointer,
                                                   true );
      std::cout
      << std::endl
      << "testSxFullCascade is sdown_L -> neutralino_2"
      << std::endl
      << "testGxFullCascade is gluino -> neutralino_2"
      << std::endl
      << "testGjsxFullCascade is gluino -> sdown_L -> neutralino_2"
      << std::endl
      << "testSbsxFullCascade is sbottom_2 -> stop_1 -> chargino_1"
      << std::endl
      << "testGjsbsxFullCascade is gluino -> sbottom_2 -> stop_1 ->"
      << " chargino_1";
      std::cout << std::endl;
      std::cout
      << std::endl
      << "testJetAcceptanceTable->getAcceptance( &gluinoWithSdownL,"
      << " &testGxFullCascade, &testSxFullCascade ) = "
      << testJetAcceptanceTable->getAcceptance( &gluinoWithSdownL,
                                                &testGxFullCascade,
                                                &testSxFullCascade );
      std::cout
      << std::endl
      << "seems fine.";
      std::cout << std::endl;
      std::cout
      << std::endl
      << "testJetAcceptanceTable->getAcceptance( &gluinoWithSdownL,"
      << " &testGjsxFullCascade, &testSxFullCascade ) = "
      << testJetAcceptanceTable->getAcceptance( &gluinoWithSdownL,
                                                &testGjsxFullCascade,
                                                &testSxFullCascade );
      std::cout
      << std::endl
      << "seems fine.";
      std::cout << std::endl;
      std::cout
      << std::endl
      << "testJetAcceptanceTable->getAcceptance( &gluinoPair,"
      << " &testGxFullCascade, &testGxFullCascade ) = "
      << testJetAcceptanceTable->getAcceptance( &gluinoPair,
                                                &testGxFullCascade,
                                                &testGxFullCascade );
      std::cout
      << std::endl
      << "seems fine.";
      std::cout << std::endl;
      std::cout
      << std::endl
      << "testJetAcceptanceTable->getAcceptance( &gluinoPair,"
      << " &testGjsxFullCascade, &testGjsxFullCascade ) = "
      << testJetAcceptanceTable->getAcceptance( &gluinoPair,
                                                &testGjsxFullCascade,
                                                &testGjsxFullCascade );
      std::cout
      << std::endl
      << "seems fine.";
      std::cout << std::endl;
      std::cout
      << std::endl
      << "testJetAcceptanceTable->getAcceptance( &gluinoPair,"
      << " &testGxFullCascade, &testGjsxFullCascade ) = "
      << testJetAcceptanceTable->getAcceptance( &gluinoPair,
                                                &testGxFullCascade,
                                                &testGjsxFullCascade );
      std::cout
      << std::endl
      << "seems fine.";
      std::cout << std::endl;
      std::cout
      << std::endl
      << "testJetAcceptanceTable->getAcceptance( &sdownLWithSdownL,"
      << " &testSxFullCascade, &testSxFullCascade ) = "
      << testJetAcceptanceTable->getAcceptance( &sdownLWithSdownL,
                                                &testSxFullCascade,
                                                &testSxFullCascade );
      std::cout
      << std::endl
      << "seems fine.";
      std::cout << std::endl;
      std::cout
      << std::endl
      << "testJetAcceptanceTable->getAcceptance("
      << " &sbottomTwoWithAntisbottomTwo, &testSbsxFullCascade,"
      << " &testSbsxFullCascade ) = "
      << testJetAcceptanceTable->getAcceptance( &sbottomTwoWithAntisbottomTwo,
                                                &testSbsxFullCascade,
                                                &testSbsxFullCascade );
      std::cout
      << std::endl
      << "seems fine.";
      std::cout << std::endl;
      std::cout << std::endl;
    }


  protected:
    CppSLHA::CppSLHA2* slhaData;
    readierForNewPoint* testReadier;
    inputHandler const* const testInputHandler;
  };

}  // end of LHC_FASER namespace

#endif /* JET_KINEMATICS_TESTING_HPP_ */
