/*
 * base_electroweak_cascade_testing.hpp
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

#ifndef BASE_ELECTROWEAK_CASCADE_TESTING_HPP_
#define BASE_ELECTROWEAK_CASCADE_TESTING_HPP_

#include "basic_testing_stuff.hpp"
#include "input_testing.hpp"
#include "../LHC-FASER_neutral_electroweak_cascade_stuff.hpp"
#include "../LHC-FASER_charged_electroweak_cascade_stuff.hpp"

namespace LHC_FASER
{
  class baseElectroweakCascadeTesting
  {
  public:
    baseElectroweakCascadeTesting( basicStuffTesting* basicStuff,
                                   inputTesting* inputTester ) :
        basicStuff( basicStuff ),
        inputTester( inputTester ),
        testInputHandler( inputTester->getInput() ),
        testAcceptanceCutSet(),
        noLeptonGapAcceptanceCutSet()
    {
      testAcceptanceCutSet.setJetCut( 40.0 );
      testAcceptanceCutSet.setPrimaryLeptonCut( 20.0 );
      testAcceptanceCutSet.setSecondaryLeptonCut( 10.0 );
      testAcceptanceCutSet.setExcludedStandardModelProducts(
                                           testInputHandler->getNotInJets4() );
      noLeptonGapAcceptanceCutSet.setJetCut( 40.0 );
      noLeptonGapAcceptanceCutSet.setPrimaryLeptonCut( 10.0 );
      noLeptonGapAcceptanceCutSet.setSecondaryLeptonCut( 10.0 );
      noLeptonGapAcceptanceCutSet.setExcludedStandardModelProducts(
                                           testInputHandler->getNotInJets4() );
    }
    ~baseElectroweakCascadeTesting()
    {
      // does nothing.
    }

    inputHandler const* const
    getInput()
    {
      return testInputHandler;
    }

    void
    performTest()
    {
      lightestNeutralinoCascade
      testLightestNeutralinoCascade( testInputHandler );

      std::cout
      << std::endl
      << "testLightestNeutralinoCascade.getAcceptance( "
      << &testAcceptanceCutSet
      << ", 0, 0, 0, 0, 0 ) = "
      << testLightestNeutralinoCascade.getAcceptance( &testAcceptanceCutSet,
                                                      0,
                                                      0,
                                                      0,
                                                      0,
                                                      0 )
      << std::endl
      << "testLightestNeutralinoCascade.getAcceptance( "
      << &testAcceptanceCutSet
      << ", 1, 1, 0, 0, 0 ) = "
      << testLightestNeutralinoCascade.getAcceptance( &testAcceptanceCutSet,
                                                      1,
                                                      1,
                                                      0,
                                                      0,
                                                      0 )
      << std::endl
      << "testLightestNeutralinoCascade.getOssfMinusOsdf( "
      << &testAcceptanceCutSet << " ) = "
      << testLightestNeutralinoCascade.getOssfMinusOsdf(
                                                       &testAcceptanceCutSet );
      std::cout << std::endl;

      leptonAcceptanceHandler
      testLeptonAcceptanceHandler( testInputHandler,
                            testInputHandler->inspectPathToKinematicsGrids() );
      leptonAcceptanceParameterSet*
      testLeptonAcceptanceParameterSet(
                testLeptonAcceptanceHandler.getLeptonAcceptanceParameterSet( 7,
                                                testInputHandler->getScharmL(),
                                          testInputHandler->getNeutralinoTwo(),
                                                                           2.0,
                                                                       9.0 ) );
      CppSLHA::particle_decay_set_handler const*
      neutralinoTwoDecay(
        testInputHandler->getNeutralinoTwo()->inspect_direct_decay_handler() );
      double
      cascadeBr( 2.0
                 * neutralinoTwoDecay->get_branching_ratio_for_exact_match(
                                                 -(CppSLHA::PDG_code::smuon_R),
                                                   CppSLHA::PDG_code::muon ) );
      double remainderBr;
      electroweakCascade*
      testElectroweakCascade( new neutralinoToSemuCascade(
                                              testLeptonAcceptanceParameterSet,
                                              testLeptonAcceptanceParameterSet,
                                                testInputHandler->getScharmL(),
                                          testInputHandler->getNeutralinoTwo(),
                                                 testInputHandler->getSmuonR(),
                                                          testInputHandler ) );
      std::cout
      << std::endl
      << "neutralinoToSemuCascade (seems OK):"
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 1, 1, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                1,
                                                1,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 1, 1 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                1,
                                                1 )
      << std::endl
      << "->getOssfMinusOsdf( ... ) = "
      << testElectroweakCascade->getOssfMinusOsdf( &testAcceptanceCutSet );
      remainderBr = getRemainderBr( cascadeBr,
                                    testElectroweakCascade,
                                    &testAcceptanceCutSet );
      std::cout
      << std::endl
      << "20 GeV / 10 GeV remainderBr = " << remainderBr
      << ", ( remainderBr / cascadeBr ) = " << ( remainderBr / cascadeBr );
      remainderBr = getRemainderBr( cascadeBr,
                                    testElectroweakCascade,
                                    &noLeptonGapAcceptanceCutSet );
      std::cout
      << std::endl
      << "10 GeV / 10 GeV remainderBr = " << remainderBr
      << ", ( remainderBr / cascadeBr ) = " << ( remainderBr / cascadeBr );
      std::cout << std::endl;
      std::cout << std::endl;
      delete testElectroweakCascade;
      testElectroweakCascade = new chargeSummedNeutralinoToSemuCascade(
                                              testLeptonAcceptanceParameterSet,
                                              testLeptonAcceptanceParameterSet,
                                                testInputHandler->getScharmL(),
                                          testInputHandler->getNeutralinoTwo(),
                                                 testInputHandler->getSmuonR(),
                                                            testInputHandler );
      std::cout
      << std::endl
      << "chargeSummedNeutralinoToSemuCascade (seems OK):"
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 1, 1, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                1,
                                                1,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 1, 1 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                1,
                                                1 )
      << std::endl
      << "->getOssfMinusOsdf( ... ) = "
      << testElectroweakCascade->getOssfMinusOsdf( &testAcceptanceCutSet );
      remainderBr = getRemainderBr( cascadeBr,
                                    testElectroweakCascade,
                                    &testAcceptanceCutSet );
      std::cout
      << std::endl
      << "20 GeV / 10 GeV remainderBr = " << remainderBr
      << ", ( remainderBr / cascadeBr ) = " << ( remainderBr / cascadeBr );
      remainderBr = getRemainderBr( cascadeBr,
                                    testElectroweakCascade,
                                    &noLeptonGapAcceptanceCutSet );
      std::cout
      << std::endl
      << "10 GeV / 10 GeV remainderBr = " << remainderBr
      << ", ( remainderBr / cascadeBr ) = " << ( remainderBr / cascadeBr );
      std::cout << std::endl;
      std::cout << std::endl;
      delete testElectroweakCascade;
      cascadeBr = ( 2.0
                    * neutralinoTwoDecay->get_branching_ratio_for_exact_match(
                                                -(CppSLHA::PDG_code::stau_one),
                                             CppSLHA::PDG_code::tau_lepton ) );
      testElectroweakCascade = new neutralinoToStauCascade(
                                              testLeptonAcceptanceParameterSet,
                                              testLeptonAcceptanceParameterSet,
                                                testInputHandler->getScharmL(),
                                          testInputHandler->getNeutralinoTwo(),
                                                testInputHandler->getStauOne(),
                                                            testInputHandler );
      std::cout
      << std::endl
      << "neutralinoToStauCascade (seems OK):"
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 1, 1, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                1,
                                                1,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 1, 1 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                1,
                                                1 )
      << std::endl
      << "->getOssfMinusOsdf( ... ) = "
      << testElectroweakCascade->getOssfMinusOsdf( &testAcceptanceCutSet );
      remainderBr = getRemainderBr( cascadeBr,
                                    testElectroweakCascade,
                                    &testAcceptanceCutSet );
      std::cout
      << std::endl
      << "20 GeV / 10 GeV remainderBr = " << remainderBr
      << ", ( remainderBr / cascadeBr ) = " << ( remainderBr / cascadeBr );
      remainderBr = getRemainderBr( cascadeBr,
                                    testElectroweakCascade,
                                    &noLeptonGapAcceptanceCutSet );
      std::cout
      << std::endl
      << "10 GeV / 10 GeV remainderBr = " << remainderBr
      << ", ( remainderBr / cascadeBr ) = " << ( remainderBr / cascadeBr );
      std::cout << std::endl;
      std::cout << std::endl;
      delete testElectroweakCascade;
      testElectroweakCascade = new chargeSummedNeutralinoToStauCascade(
                                              testLeptonAcceptanceParameterSet,
                                              testLeptonAcceptanceParameterSet,
                                                testInputHandler->getScharmL(),
                                          testInputHandler->getNeutralinoTwo(),
                                                testInputHandler->getStauOne(),
                                                            testInputHandler );
      std::cout
      << std::endl
      << "chargeSummedNeutralinoToStauCascade (seems OK):"
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 1, 1, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                1,
                                                1,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 1, 1 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                1,
                                                1 )
      << std::endl
      << "->getOssfMinusOsdf( ... ) = "
      << testElectroweakCascade->getOssfMinusOsdf( &testAcceptanceCutSet );
      remainderBr = getRemainderBr( cascadeBr,
                                    testElectroweakCascade,
                                    &testAcceptanceCutSet );
      std::cout
      << std::endl
      << "20 GeV / 10 GeV remainderBr = " << remainderBr
      << ", ( remainderBr / cascadeBr ) = " << ( remainderBr / cascadeBr );
      remainderBr = getRemainderBr( cascadeBr,
                                    testElectroweakCascade,
                                    &noLeptonGapAcceptanceCutSet );
      std::cout
      << std::endl
      << "10 GeV / 10 GeV remainderBr = " << remainderBr
      << ", ( remainderBr / cascadeBr ) = " << ( remainderBr / cascadeBr );
      std::cout << std::endl;
      std::cout << std::endl;
      delete testElectroweakCascade;
      cascadeBr = neutralinoTwoDecay->get_branching_ratio_for_exact_match(
                                                          CppSLHA::PDG_code::Z,
                                           CppSLHA::PDG_code::neutralino_one );
      testElectroweakCascade = new neutralinoToZCascade(
                                              testLeptonAcceptanceParameterSet,
                                              testLeptonAcceptanceParameterSet,
                                                testInputHandler->getScharmL(),
                                          testInputHandler->getNeutralinoTwo(),
                                                            testInputHandler );
      std::cout
      << std::endl
      << "neutralinoToZCascade (seems OK):"
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 1, 1, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                1,
                                                1,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 1, 1 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                1,
                                                1 )
      << std::endl
      << "->getOssfMinusOsdf( ... ) = "
      << testElectroweakCascade->getOssfMinusOsdf( &testAcceptanceCutSet );
      remainderBr = getRemainderBr( cascadeBr,
                                    testElectroweakCascade,
                                    &testAcceptanceCutSet );
      std::cout
      << std::endl
      << "20 GeV / 10 GeV remainderBr = " << remainderBr
      << ", ( remainderBr / cascadeBr ) = " << ( remainderBr / cascadeBr );
      remainderBr = getRemainderBr( cascadeBr,
                                    testElectroweakCascade,
                                    &noLeptonGapAcceptanceCutSet );
      std::cout
      << std::endl
      << "10 GeV / 10 GeV remainderBr = " << remainderBr
      << ", ( remainderBr / cascadeBr ) = " << ( remainderBr / cascadeBr );
      std::cout << std::endl;
      std::cout << std::endl;
      delete testElectroweakCascade;
      testElectroweakCascade = new chargeSummedNeutralinoToZCascade(
                                              testLeptonAcceptanceParameterSet,
                                              testLeptonAcceptanceParameterSet,
                                                testInputHandler->getScharmL(),
                                          testInputHandler->getNeutralinoTwo(),
                                                            testInputHandler );
      std::cout
      << std::endl
      << "chargeSummedNeutralinoToZCascade (seems OK):"
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 1, 1, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                1,
                                                1,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 1, 1 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                1,
                                                1 )
      << std::endl
      << "->getOssfMinusOsdf( ... ) = "
      << testElectroweakCascade->getOssfMinusOsdf( &testAcceptanceCutSet );
      remainderBr = getRemainderBr( cascadeBr,
                                    testElectroweakCascade,
                                    &testAcceptanceCutSet );
      std::cout
      << std::endl
      << "20 GeV / 10 GeV remainderBr = " << remainderBr
      << ", ( remainderBr / cascadeBr ) = " << ( remainderBr / cascadeBr );
      remainderBr = getRemainderBr( cascadeBr,
                                    testElectroweakCascade,
                                    &noLeptonGapAcceptanceCutSet );
      std::cout
      << std::endl
      << "10 GeV / 10 GeV remainderBr = " << remainderBr
      << ", ( remainderBr / cascadeBr ) = " << ( remainderBr / cascadeBr );
      std::cout << std::endl;
      std::cout << std::endl;
      delete testElectroweakCascade;
      cascadeBr
      = testInputHandler->getNeutralinoFour()->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                                          CppSLHA::PDG_code::Z,
                                           CppSLHA::PDG_code::neutralino_one );
      testElectroweakCascade = new neutralinoToZCascade(
                                              testLeptonAcceptanceParameterSet,
                                              testLeptonAcceptanceParameterSet,
                                                testInputHandler->getScharmL(),
                                         testInputHandler->getNeutralinoFour(),
                                                            testInputHandler );
      std::cout
      << std::endl
      << "neutralinoToZCascade by chi^0_4 (seems OK):";
      std::cout
      << std::endl;
      std::cout
      << "->getAcceptance( ..., 0, 0, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 1, 1, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                1,
                                                1,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 1, 1 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                1,
                                                1 )
      << std::endl
      << "->getOssfMinusOsdf( ... ) = "
      << testElectroweakCascade->getOssfMinusOsdf( &testAcceptanceCutSet );
      remainderBr = getRemainderBr( cascadeBr,
                                    testElectroweakCascade,
                                    &testAcceptanceCutSet );
      std::cout
      << std::endl
      << "20 GeV / 10 GeV remainderBr = " << remainderBr
      << ", ( remainderBr / cascadeBr ) = " << ( remainderBr / cascadeBr );
      remainderBr = getRemainderBr( cascadeBr,
                                    testElectroweakCascade,
                                    &noLeptonGapAcceptanceCutSet );
      std::cout
      << std::endl
      << "10 GeV / 10 GeV remainderBr = " << remainderBr
      << ", ( remainderBr / cascadeBr ) = " << ( remainderBr / cascadeBr );
      std::cout << std::endl;
      std::cout << std::endl;
      delete testElectroweakCascade;
      testElectroweakCascade = new chargeSummedNeutralinoToZCascade(
                                              testLeptonAcceptanceParameterSet,
                                              testLeptonAcceptanceParameterSet,
                                                testInputHandler->getScharmL(),
                                         testInputHandler->getNeutralinoFour(),
                                                            testInputHandler );
      std::cout
      << std::endl
      << "chargeSummedNeutralinoToZCascade by chi^0_4 (seems OK):"
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 1, 1, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                1,
                                                1,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 1, 1 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                1,
                                                1 )
      << std::endl
      << "->getOssfMinusOsdf( ... ) = "
      << testElectroweakCascade->getOssfMinusOsdf( &testAcceptanceCutSet );
      remainderBr = getRemainderBr( cascadeBr,
                                    testElectroweakCascade,
                                    &testAcceptanceCutSet );
      std::cout
      << std::endl
      << "20 GeV / 10 GeV remainderBr = " << remainderBr
      << ", ( remainderBr / cascadeBr ) = " << ( remainderBr / cascadeBr );
      remainderBr = getRemainderBr( cascadeBr,
                                    testElectroweakCascade,
                                    &noLeptonGapAcceptanceCutSet );
      std::cout
      << std::endl
      << "10 GeV / 10 GeV remainderBr = " << remainderBr
      << ", ( remainderBr / cascadeBr ) = " << ( remainderBr / cascadeBr );
      std::cout << std::endl;
      std::cout << std::endl;
      delete testElectroweakCascade;
      cascadeBr
      = testInputHandler->getNeutralinoFour()->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                  CppSLHA::PDG_code::light_neutral_EWSB_scalar,
                                           CppSLHA::PDG_code::neutralino_one );
      testElectroweakCascade = new neutralinoToHiggsCascade(
                                              testLeptonAcceptanceParameterSet,
                                              testLeptonAcceptanceParameterSet,
                                                testInputHandler->getScharmL(),
                                         testInputHandler->getNeutralinoFour(),
                                 testInputHandler->getLightNeutralEwsbScalar(),
                                                            testInputHandler );
      std::cout
      << std::endl
      << "neutralinoToHiggsCascade"
      << " (seems OK, since 5% of the h BR is photon pair or 3-body):"
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 1, 1, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                1,
                                                1,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 1, 1 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                1,
                                                1 )
      << std::endl
      << "->getOssfMinusOsdf( ... ) = "
      << testElectroweakCascade->getOssfMinusOsdf( &testAcceptanceCutSet );
      remainderBr = getRemainderBr( cascadeBr,
                                    testElectroweakCascade,
                                    &testAcceptanceCutSet );
      std::cout
      << std::endl
      << "20 GeV / 10 GeV remainderBr = " << remainderBr
      << ", ( remainderBr / cascadeBr ) = " << ( remainderBr / cascadeBr );
      remainderBr = getRemainderBr( cascadeBr,
                                    testElectroweakCascade,
                                    &noLeptonGapAcceptanceCutSet );
      std::cout
      << std::endl
      << "10 GeV / 10 GeV remainderBr = " << remainderBr
      << ", ( remainderBr / cascadeBr ) = " << ( remainderBr / cascadeBr );
      std::cout << std::endl;
      std::cout << std::endl;
      delete testElectroweakCascade;
      testElectroweakCascade = new chargeSummedNeutralinoVirtualCascade(
                                              testLeptonAcceptanceParameterSet,
                                              testLeptonAcceptanceParameterSet,
                                                testInputHandler->getScharmL(),
                                          testInputHandler->getNeutralinoTwo(),
                                                            testInputHandler );
      std::cout
      << std::endl
      << "chargeSummedNeutralinoVirtualCascade"
      << " (probably fine, not really checked though):"
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 1, 1, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                1,
                                                1,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 1, 1 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                1,
                                                1 )
      << std::endl
      << "->getOssfMinusOsdf( ... ) = "
      << testElectroweakCascade->getOssfMinusOsdf( &testAcceptanceCutSet );
      delete testElectroweakCascade;
      testElectroweakCascade = new neutralinoVirtualCascade(
                                              testLeptonAcceptanceParameterSet,
                                              testLeptonAcceptanceParameterSet,
                                                testInputHandler->getScharmL(),
                                          testInputHandler->getNeutralinoTwo(),
                                                            testInputHandler );
      std::cout
      << std::endl
      << "neutralinoVirtualCascade"
      << " (probably fine, not really checked though):"
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 1, 1, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                1,
                                                1,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 1, 1 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                1,
                                                1 )
      << std::endl
      << "->getOssfMinusOsdf( ... ) = "
      << testElectroweakCascade->getOssfMinusOsdf( &testAcceptanceCutSet );
      std::cout << std::endl;
      std::cout << std::endl;
      delete testElectroweakCascade;
      cascadeBr = testInputHandler->getStopTwo()->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                                   CppSLHA::PDG_code::stop_one,
                                                        CppSLHA::PDG_code::Z );
      testElectroweakCascade = new scoloredToZPlusScoloredCascade(
                                              testLeptonAcceptanceParameterSet,
                                              testLeptonAcceptanceParameterSet,
                                                testInputHandler->getStopTwo(),
                                                testInputHandler->getStopOne(),
                                                            testInputHandler );
      std::cout
      << std::endl
      << "scoloredToZPlusScoloredCascade (seems OK):"
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 1, 1, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                1,
                                                1,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 1, 1 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                1,
                                                1 )
      << std::endl
      << "->getOssfMinusOsdf( ... ) = "
      << testElectroweakCascade->getOssfMinusOsdf( &testAcceptanceCutSet );
      remainderBr = getRemainderBr( cascadeBr,
                                    testElectroweakCascade,
                                    &testAcceptanceCutSet );
      std::cout
      << std::endl
      << "20 GeV / 10 GeV remainderBr = " << remainderBr
      << ", ( remainderBr / cascadeBr ) = " << ( remainderBr / cascadeBr );
      remainderBr = getRemainderBr( cascadeBr,
                                    testElectroweakCascade,
                                    &noLeptonGapAcceptanceCutSet );
      std::cout
      << std::endl
      << "10 GeV / 10 GeV remainderBr = " << remainderBr
      << ", ( remainderBr / cascadeBr ) = " << ( remainderBr / cascadeBr );
      std::cout << std::endl;
      std::cout << std::endl;
      delete testElectroweakCascade;
      cascadeBr = testInputHandler->getStopTwo()->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                                   CppSLHA::PDG_code::stop_one,
                                CppSLHA::PDG_code::light_neutral_EWSB_scalar );
      testElectroweakCascade = new scoloredToHiggsPlusScoloredCascade(
                                              testLeptonAcceptanceParameterSet,
                                              testLeptonAcceptanceParameterSet,
                                                testInputHandler->getStopTwo(),
                                                testInputHandler->getStopOne(),
                                 testInputHandler->getLightNeutralEwsbScalar(),
                                                            testInputHandler );
      std::cout
      << std::endl
      << "scoloredToHiggsPlusScoloredCascade"
      << " (seems OK, since 5% of the h BR is photon pair or 3-body):"
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 1, 1, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                1,
                                                1,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 1, 1 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                1,
                                                1 )
      << std::endl
      << "->getOssfMinusOsdf( ... ) = "
      << testElectroweakCascade->getOssfMinusOsdf( &testAcceptanceCutSet );
      remainderBr = getRemainderBr( cascadeBr,
                                    testElectroweakCascade,
                                    &testAcceptanceCutSet );
      std::cout
      << std::endl
      << "20 GeV / 10 GeV remainderBr = " << remainderBr
      << ", ( remainderBr / cascadeBr ) = " << ( remainderBr / cascadeBr );
      remainderBr = getRemainderBr( cascadeBr,
                                    testElectroweakCascade,
                                    &noLeptonGapAcceptanceCutSet );
      std::cout
      << std::endl
      << "10 GeV / 10 GeV remainderBr = " << remainderBr
      << ", ( remainderBr / cascadeBr ) = " << ( remainderBr / cascadeBr );
      std::cout << std::endl;
      std::cout << std::endl;
      delete testElectroweakCascade;
      CppSLHA::particle_decay_set_handler const*
      charginoOneDecay
      = testInputHandler->getCharginoOne()->inspect_direct_decay_handler();
      cascadeBr = charginoOneDecay->get_branching_ratio_for_exact_match(
                                                 -(CppSLHA::PDG_code::smuon_R),
                                            CppSLHA::PDG_code::muon_neutrino );
      testElectroweakCascade = new charginoToSemuOrEmuSnuCascade(
                                              testLeptonAcceptanceParameterSet,
                                              testLeptonAcceptanceParameterSet,
                                                testInputHandler->getScharmL(),
                                            testInputHandler->getCharginoOne(),
                                                 testInputHandler->getSmuonR(),
                                                              testInputHandler,
                                                                       false );
      std::cout
      << std::endl
      << "charginoToSemuOrEmuSnuCascade (seems OK):"
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 1, 1, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                1,
                                                1,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 1, 1 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                1,
                                                1 )
      << std::endl
      << "->getOssfMinusOsdf( ... ) = "
      << testElectroweakCascade->getOssfMinusOsdf( &testAcceptanceCutSet );
      remainderBr = getRemainderBr( cascadeBr,
                                    testElectroweakCascade,
                                    &testAcceptanceCutSet );
      std::cout
      << std::endl
      << "20 GeV / 10 GeV remainderBr = " << remainderBr
      << ", ( remainderBr / cascadeBr ) = " << ( remainderBr / cascadeBr );
      remainderBr = getRemainderBr( cascadeBr,
                                    testElectroweakCascade,
                                    &noLeptonGapAcceptanceCutSet );
      std::cout
      << std::endl
      << "10 GeV / 10 GeV remainderBr = " << remainderBr
      << ", ( remainderBr / cascadeBr ) = " << ( remainderBr / cascadeBr );
      std::cout << std::endl;
      std::cout << std::endl;
      delete testElectroweakCascade;
      cascadeBr = charginoOneDecay->get_branching_ratio_for_exact_match(
                                                 -(CppSLHA::PDG_code::smuon_R),
                                            CppSLHA::PDG_code::muon_neutrino );
      testElectroweakCascade = new charginoToSemuOrEmuSnuCascade(
                                              testLeptonAcceptanceParameterSet,
                                              testLeptonAcceptanceParameterSet,
                                                testInputHandler->getScharmL(),
                                            testInputHandler->getCharginoOne(),
                                                 testInputHandler->getSmuonR(),
                                                              testInputHandler,
                                                                  true );
      std::cout
      << std::endl
      << "charginoToSemuOrEmuSnuCascade pretending to be sneutrino cascade"
      << " (seems OK):"
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 1, 1, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                1,
                                                1,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 1, 1 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                1,
                                                1 )
      << std::endl
      << "->getOssfMinusOsdf( ... ) = "
      << testElectroweakCascade->getOssfMinusOsdf( &testAcceptanceCutSet );
      remainderBr = getRemainderBr( cascadeBr,
                                    testElectroweakCascade,
                                    &testAcceptanceCutSet );
      std::cout
      << std::endl
      << "20 GeV / 10 GeV remainderBr = " << remainderBr
      << ", ( remainderBr / cascadeBr ) = " << ( remainderBr / cascadeBr );
      remainderBr = getRemainderBr( cascadeBr,
                                    testElectroweakCascade,
                                    &noLeptonGapAcceptanceCutSet );
      std::cout
      << std::endl
      << "10 GeV / 10 GeV remainderBr = " << remainderBr
      << ", ( remainderBr / cascadeBr ) = " << ( remainderBr / cascadeBr );
      std::cout << std::endl;
      std::cout << std::endl;
      delete testElectroweakCascade;
      cascadeBr = charginoOneDecay->get_branching_ratio_for_exact_match(
                                                -(CppSLHA::PDG_code::stau_one),
                                             CppSLHA::PDG_code::tau_neutrino );
      testElectroweakCascade = new charginoToStauOrTauSnuCascade(
                                              testLeptonAcceptanceParameterSet,
                                              testLeptonAcceptanceParameterSet,
                                                testInputHandler->getScharmL(),
                                            testInputHandler->getCharginoOne(),
                                                testInputHandler->getStauOne(),
                                                             testInputHandler,
                                                                  false );
      std::cout
      << std::endl
      << "charginoToStauOrTauSnuCascade (seems OK):"
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 1, 1, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                1,
                                                1,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 1, 1 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                1,
                                                1 )
      << std::endl
      << "->getOssfMinusOsdf( ... ) = "
      << testElectroweakCascade->getOssfMinusOsdf( &testAcceptanceCutSet );
      remainderBr = getRemainderBr( cascadeBr,
                                    testElectroweakCascade,
                                    &testAcceptanceCutSet );
      std::cout
      << std::endl
      << "20 GeV / 10 GeV remainderBr = " << remainderBr
      << ", ( remainderBr / cascadeBr ) = " << ( remainderBr / cascadeBr );
      remainderBr = getRemainderBr( cascadeBr,
                                    testElectroweakCascade,
                                    &noLeptonGapAcceptanceCutSet );
      std::cout
      << std::endl
      << "10 GeV / 10 GeV remainderBr = " << remainderBr
      << ", ( remainderBr / cascadeBr ) = " << ( remainderBr / cascadeBr );
      std::cout << std::endl;
      std::cout << std::endl;
      delete testElectroweakCascade;
      cascadeBr = charginoOneDecay->get_branching_ratio_for_exact_match(
                                                     CppSLHA::PDG_code::W_plus,
                                           CppSLHA::PDG_code::neutralino_one );
      testElectroweakCascade = new charginoToWCascade(
                                              testLeptonAcceptanceParameterSet,
                                              testLeptonAcceptanceParameterSet,
                                                testInputHandler->getScharmL(),
                                            testInputHandler->getCharginoOne(),
                                                            testInputHandler );
      std::cout
      << std::endl
      << "charginoToWCascade (seems OK):"
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 1, 1, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                1,
                                                1,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 1, 1 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                1,
                                                1 )
      << std::endl
      << "->getOssfMinusOsdf( ... ) = "
      << testElectroweakCascade->getOssfMinusOsdf( &testAcceptanceCutSet );
      remainderBr = getRemainderBr( cascadeBr,
                                    testElectroweakCascade,
                                    &testAcceptanceCutSet );
      std::cout
      << std::endl
      << "20 GeV / 10 GeV remainderBr = " << remainderBr
      << ", ( remainderBr / cascadeBr ) = " << ( remainderBr / cascadeBr );
      remainderBr = getRemainderBr( cascadeBr,
                                    testElectroweakCascade,
                                    &noLeptonGapAcceptanceCutSet );
      std::cout
      << std::endl
      << "10 GeV / 10 GeV remainderBr = " << remainderBr
      << ", ( remainderBr / cascadeBr ) = " << ( remainderBr / cascadeBr );
      std::cout << std::endl;
      std::cout << std::endl;
      delete testElectroweakCascade;
      cascadeBr = charginoOneDecay->get_branching_ratio_for_exact_match(
                                                     CppSLHA::PDG_code::W_plus,
                                           CppSLHA::PDG_code::neutralino_one );
      testElectroweakCascade = new charginoToHiggsCascade(
                                              testLeptonAcceptanceParameterSet,
                                              testLeptonAcceptanceParameterSet,
                                                testInputHandler->getScharmL(),
                                            testInputHandler->getCharginoOne(),
                                      testInputHandler->getChargedEwsbScalar(),
                                                            testInputHandler );
      std::cout
      << std::endl
      << "charginoToHiggsCascade"
      << " (needs checking, but it's going to need a hacked spectrum):"
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 1, 1, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                1,
                                                1,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 1, 1 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                1,
                                                1 )
      << std::endl
      << "->getOssfMinusOsdf( ... ) = "
      << testElectroweakCascade->getOssfMinusOsdf( &testAcceptanceCutSet );
      remainderBr = getRemainderBr( cascadeBr,
                                    testElectroweakCascade,
                                    &testAcceptanceCutSet );
      std::cout
      << std::endl
      << "20 GeV / 10 GeV remainderBr = " << remainderBr
      << ", ( remainderBr / cascadeBr ) = " << ( remainderBr / cascadeBr );
      remainderBr = getRemainderBr( cascadeBr,
                                    testElectroweakCascade,
                                    &noLeptonGapAcceptanceCutSet );
      std::cout
      << std::endl
      << "10 GeV / 10 GeV remainderBr = " << remainderBr
      << ", ( remainderBr / cascadeBr ) = " << ( remainderBr / cascadeBr );
      std::cout << std::endl;
      std::cout << std::endl;
      delete testElectroweakCascade;
      testElectroweakCascade = new charginoVirtualCascade(
                                              testLeptonAcceptanceParameterSet,
                                              testLeptonAcceptanceParameterSet,
                                                testInputHandler->getScharmL(),
                                          testInputHandler->getNeutralinoTwo(),
                                                            testInputHandler );
      std::cout
      << std::endl
      << "charginoVirtualCascade"
      << " (probably fine, not really checked though):"
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 1, 1, 0, 0, 0 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                1,
                                                1,
                                                0,
                                                0,
                                                0 )
      << std::endl
      << "->getAcceptance( ..., 0, 0, 0, 1, 1 ) = "
      << testElectroweakCascade->getAcceptance( &testAcceptanceCutSet,
                                                0,
                                                0,
                                                0,
                                                1,
                                                1 )
      << std::endl
      << "->getOssfMinusOsdf( ... ) = "
      << testElectroweakCascade->getOssfMinusOsdf( &testAcceptanceCutSet );
      std::cout << std::endl;
      std::cout << std::endl;
      delete testElectroweakCascade;
    }


  protected:
    basicStuffTesting* basicStuff;
    inputTesting* inputTester;
    inputHandler const* const testInputHandler;
    acceptanceCutSet testAcceptanceCutSet;
    acceptanceCutSet noLeptonGapAcceptanceCutSet;

    double
    getRemainderBr( double const cascadeBr,
                    electroweakCascade* testElectroweakCascade,
                    acceptanceCutSet* testAcceptanceCutSet )
    {
      double previousRemainder( cascadeBr );
      double returnValue( cascadeBr );
      /**std::cout
      << std::endl
      << "remainder starts = " << cascadeBr;**/
      for( int jetCounter( -1 );
           3 > jetCounter;
           ++jetCounter )
      {
        for( int eMinusCounter( -1 );
             3 > eMinusCounter;
             ++eMinusCounter )
        {
          for( int ePlusCounter( -1 );
               3 > ePlusCounter;
               ++ePlusCounter )
          {
            for( int muMinusCounter( -1 );
                 3 > muMinusCounter;
                 ++muMinusCounter )
            {
              for( int muPlusCounter( -1 );
                   3 > muPlusCounter;
                   ++muPlusCounter )
              {
                returnValue
                -= testElectroweakCascade->getAcceptance( testAcceptanceCutSet,
                                                          jetCounter,
                                                          eMinusCounter,
                                                          ePlusCounter,
                                                          muMinusCounter,
                                                          muPlusCounter );
                if( previousRemainder != returnValue )
                {
                  /**std::cout
                  << std::endl
                  << jetCounter << eMinusCounter << ePlusCounter
                  << muMinusCounter << muPlusCounter << " remainder = "
                  << returnValue;**/
                  previousRemainder = returnValue;
                }
              }
            }
          }
        }
      }
      return returnValue;
    }
  };

}  // end of LHC_FASER namespace

#endif /* BASE_ELECTROWEAK_CASCADE_TESTING_HPP_ */
