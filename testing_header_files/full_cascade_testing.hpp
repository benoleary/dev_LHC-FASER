/*
 * full_cascade_testing.hpp
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

#ifndef FULL_CASCADE_TESTING_HPP_
#define FULL_CASCADE_TESTING_HPP_

#include "basic_testing_stuff.hpp"
#include "input_testing.hpp"
#include "../LHC-FASER_full_cascade_stuff.hpp"

namespace LHC_FASER
{
  class fullCascadeTesting
  {
  public:
    fullCascadeTesting( basicStuffTesting* basicStuff,
                        inputTesting* inputTester ) :
        basicStuff( basicStuff ),
        inputTester( inputTester ),
        testInputHandler( inputTester->getInput() ),
        testAcceptanceCutSet()
    {
      testAcceptanceCutSet.setJetCut( 40.0 );
      testAcceptanceCutSet.setPrimaryLeptonCut( 20.0 );
      testAcceptanceCutSet.setSecondaryLeptonCut( 10.0 );
      testAcceptanceCutSet.setExcludedStandardModelProducts(
                                           testInputHandler->getNotInJets4() );
    }
    ~fullCascadeTesting()
    {
      // does nothing.
    }


    void
    performTest()
    {
      particlePointer sdownL( testInputHandler->getSdownL() );
      particlePointer stopOne( testInputHandler->getStopOne() );
      particlePointer sbottomOne( testInputHandler->getSbottomOne() );
      particlePointer sbottomTwo( testInputHandler->getSbottomTwo() );
      particlePointer gluinoPointer( testInputHandler->getGluino() );
      particlePointer neutralinoOne( testInputHandler->getNeutralinoOne() );
      particlePointer neutralinoTwo( testInputHandler->getNeutralinoTwo() );
      particlePointer charginoOne( testInputHandler->getCharginoOne() );
      particlePointer wPlus( testInputHandler->getWPlus() );

      // 1st, SPS1a:

      electroweakCascadeHandler
      testElectroweakCascadeHandler( testInputHandler );
      electroweakCascadesForOneBeamEnergy*
      cascadesForSevenTev(
          testElectroweakCascadeHandler.getElectroweakCascadesForOneBeamEnergy(
                                                                         7 ) );
      electroweakCascadeSet*
      testEwinoCascadeSet( cascadesForSevenTev->getCascadeSet( sdownL,
                                                             neutralinoTwo ) );
      fullCascadeType::squarkDirectlyToElectroweakType::sdownType
      testSdNTwoFullCascade( testInputHandler,
                             sdownL,
                             7,
                             testEwinoCascadeSet );
      std::cout
      << std::endl
      << "testSdNTwoFullCascade set as sdown_L -> chi^0_2 ( "
      << testSdNTwoFullCascade.getAsString() << " )";
      std::cout << std::endl;
      if( testSdNTwoFullCascade.isOpen() )
      {
        std::cout
        << std::endl
        << "testSdNTwoFullCascade is open.";
        std::cout << std::endl;
      }
      else
      {
        std::cout
        << std::endl
        << "testSdNTwoFullCascade is not open?!";
        std::cout << std::endl;
      }
      std::list< int > testSmParticleExclusion;
      std::cout
      << std::endl
      << "BR with no excluded SM fermions = "
      << testSdNTwoFullCascade.getBrToEwino( &testSmParticleExclusion );
      testSmParticleExclusion.push_back( CppSLHA::PDG_code::down );
      std::cout
      << std::endl
      << "BR excluding downs (but not antidowns) in final state = "
      << testSdNTwoFullCascade.getBrToEwino( &testSmParticleExclusion );
      std::cout << std::endl;
      std::cout
      << std::endl
      << ".getAcceptance( [particle], &testAcceptanceCutSet, 0, 1, 1, 0, 0 )"
      << " = " << testSdNTwoFullCascade.getAcceptance( true,
                                                       &testAcceptanceCutSet,
                                                       0,
                                                       1,
                                                       1,
                                                       0,
                                                       0 );
      std::cout
      << std::endl
      << ".getAcceptance("
      << " [antiparticle], &testAcceptanceCutSet, 0, 1, 1, 0, 0 ) = "
      << testSdNTwoFullCascade.getAcceptance( false,
                                              &testAcceptanceCutSet,
                                              0,
                                              1,
                                              1,
                                              0,
                                              0 );
      std::cout
      << std::endl
      << "testSdNTwoFullCascade seems OK.";
      std::cout << std::endl;
      std::cout << std::endl;
      testEwinoCascadeSet = cascadesForSevenTev->getCascadeSet( sdownL,
                                                                charginoOne );
      fullCascadeType::squarkDirectlyToElectroweakType::sdownType
      testSdXOneFullCascade( testInputHandler,
                             sdownL,
                             7,
                             testEwinoCascadeSet );
      std::cout
      << std::endl
      << "testSdXOneFullCascade set as sdown_L -> chi^-_1 ( "
      << testSdXOneFullCascade.getAsString() << " )";
      std::cout << std::endl;
      if( testSdXOneFullCascade.isOpen() )
      {
        std::cout
        << std::endl
        << "testSdXOneFullCascade is open.";
        std::cout << std::endl;
      }
      else
      {
        std::cout
        << std::endl
        << "testSdXOneFullCascade is not open?!";
        std::cout << std::endl;
      }
      testSmParticleExclusion.clear();
      std::cout
      << std::endl
      << "BR with no excluded SM fermions = "
      << testSdXOneFullCascade.getBrToEwino( &testSmParticleExclusion );
      testSmParticleExclusion.push_back( CppSLHA::PDG_code::up );
      testSmParticleExclusion.push_back( -(CppSLHA::PDG_code::up) );
      std::cout
      << std::endl
      << "BR excluding ups & antiups in final state = "
      << testSdXOneFullCascade.getBrToEwino( &testSmParticleExclusion );
      std::cout << std::endl;
      std::cout
      << std::endl
      << ".getAcceptance( [particle], &testAcceptanceCutSet, 0, 1, 0, 0, 0 )"
      << " = " << testSdXOneFullCascade.getAcceptance( true,
                                                       &testAcceptanceCutSet,
                                                       0,
                                                       1,
                                                       0,
                                                       0,
                                                       0 );
      std::cout
      << std::endl
      << ".getAcceptance("
      << " [antiparticle], &testAcceptanceCutSet, 0, 1, 0, 0, 0 ) = "
      << testSdXOneFullCascade.getAcceptance( false,
                                              &testAcceptanceCutSet,
                                              0,
                                              1,
                                              0,
                                              0,
                                              0 );
      std::cout
      << std::endl
      << ".getAcceptance("
      << " [antiparticle], &testAcceptanceCutSet, 0, 0, 1, 0, 0 ) = "
      << testSdXOneFullCascade.getAcceptance( false,
                                              &testAcceptanceCutSet,
                                              0,
                                              0,
                                              1,
                                              0,
                                              0 );
      std::cout
      << std::endl
      << "testSdXOneFullCascade seems OK.";
      std::cout << std::endl;
      std::cout << std::endl;
      testEwinoCascadeSet = cascadesForSevenTev->getCascadeSet( sdownL,
                                          testInputHandler->getCharginoTwo() );
      fullCascadeType::squarkDirectlyToElectroweakType::sdownType
      testSdXTwoFullCascade( testInputHandler,
                             sdownL,
                             7,
                             testEwinoCascadeSet );
      std::cout
      << std::endl
      << "testSdXTwoFullCascade set as sdown_L -> chi^-_2 ( "
      << testSdXTwoFullCascade.getAsString() << " )";
      std::cout << std::endl;
      if( testSdXTwoFullCascade.isOpen() )
      {
        std::cout
        << std::endl
        << "testSdXTwoFullCascade is open.";
        std::cout << std::endl;
      }
      else
      {
        std::cout
        << std::endl
        << "testSdXTwoFullCascade is not open?!";
        std::cout << std::endl;
      }
      testSmParticleExclusion.clear();
      std::cout
      << std::endl
      << "BR with no excluded SM fermions = "
      << testSdXTwoFullCascade.getBrToEwino( &testSmParticleExclusion );
      testSmParticleExclusion.push_back( CppSLHA::PDG_code::up );
      testSmParticleExclusion.push_back( -(CppSLHA::PDG_code::up) );
      std::cout
      << std::endl
      << "BR excluding ups & antiups in final state = "
      << testSdXTwoFullCascade.getBrToEwino( &testSmParticleExclusion );
      std::cout << std::endl;
      std::cout
      << std::endl
      << ".getAcceptance( [particle], &testAcceptanceCutSet, 0, 1, 0, 0, 0 )"
      << " = " << testSdXTwoFullCascade.getAcceptance( true,
                                                       &testAcceptanceCutSet,
                                                       0,
                                                       1,
                                                       0,
                                                       0,
                                                       0 );
      std::cout
      << std::endl
      << ".getAcceptance("
      << " [antiparticle], &testAcceptanceCutSet, 0, 1, 0, 0, 0 ) = "
      << testSdXTwoFullCascade.getAcceptance( false,
                                              &testAcceptanceCutSet,
                                              0,
                                              1,
                                              0,
                                              0,
                                              0 );
      std::cout
      << std::endl
      << ".getAcceptance("
      << " [antiparticle], &testAcceptanceCutSet, 0, 0, 1, 0, 0 ) = "
      << testSdXTwoFullCascade.getAcceptance( false,
                                              &testAcceptanceCutSet,
                                              0,
                                              0,
                                              1,
                                              0,
                                              0 );
      std::cout
      << std::endl
      << "testSdXTwoFullCascade seems OK (large amount of BR lost by ignoring"
      << " channel to W^- + neutralino_2).";
      std::cout << std::endl;
      std::cout << std::endl;
      testEwinoCascadeSet = cascadesForSevenTev->getCascadeSet( gluinoPointer,
                                                               neutralinoTwo );
      fullCascadeType::gluinoDirectlyToElectroweak
      testGoNTwoFullCascade( testInputHandler,
                             7,
                             testEwinoCascadeSet );
      std::cout
      << std::endl
      << "testGoNTwoFullCascade set as gluino -> chi^0_2 ( "
      << testGoNTwoFullCascade.getAsString() << " )";
      std::cout << std::endl;
      if( testGoNTwoFullCascade.isOpen() )
      {
        std::cout
        << std::endl
        << "testGoNTwoFullCascade is open?1";
        std::cout << std::endl;
      }
      else
      {
        std::cout
        << std::endl
        << "testGoNTwoFullCascade is not open (BR should be too low, despite"
        << " masses checking out).";
        std::cout << std::endl;
      }
      testSmParticleExclusion.clear();
      std::cout
      << std::endl
      << "BR with no excluded SM fermions = "
      << testGoNTwoFullCascade.getBrToEwino( &testSmParticleExclusion );
      testSmParticleExclusion.push_back( CppSLHA::PDG_code::down );
      std::cout
      << std::endl
      << "BR excluding downs (but not antidowns) in final state = "
      << testGoNTwoFullCascade.getBrToEwino( &testSmParticleExclusion );
      std::cout
      << std::endl
      << "testGoNTwoFullCascade is kinematically allowed, & gives non-zero"
      << " acceptances, but the BR is so low that the code snaps it to zero.";
      std::cout << std::endl;
      std::cout
      << std::endl
      << ".getAcceptance( [particle], &testAcceptanceCutSet, 0, 1, 1, 0, 0 )"
      << " = " << testGoNTwoFullCascade.getAcceptance( true,
                                                       &testAcceptanceCutSet,
                                                       0,
                                                       1,
                                                       1,
                                                       0,
                                                       0 );
      std::cout
      << std::endl
      << ".getAcceptance("
      << " [antiparticle], &testAcceptanceCutSet, 0, 1, 1, 0, 0 ) = "
      << testGoNTwoFullCascade.getAcceptance( false,
                                              &testAcceptanceCutSet,
                                              0,
                                              1,
                                              1,
                                              0,
                                              0 );
      std::cout << std::endl;
      std::cout
      << std::endl
      << "testGoNTwoFullCascade seems OK.";
      std::cout << std::endl;
      std::cout << std::endl;
      testEwinoCascadeSet = cascadesForSevenTev->getCascadeSet( stopOne,
                                                                charginoOne );

      effectiveSquarkMassHolder*
      stopMinusW( testInputHandler->getSquarkMinusBosonEffectiveMass( stopOne,
                                                                     wPlus,
                                                             neutralinoOne ) );
      effectiveSquarkMassHolder*
      sbottomPlusW( testInputHandler->getSquarkPlusBosonEffectiveMass(
                                                                    sbottomOne,
                                                                       wPlus,
                                                             neutralinoOne ) );
      std::cout
      << std::endl
      << "stopMinusW = " << stopMinusW << ", sbottomPlusW = " << sbottomPlusW;
      std::cout << std::endl;
      std::cout
      << std::endl
      << "stopMinusW->getEffectiveSquarkMass() = "
      << stopMinusW->getEffectiveSquarkMass()
      << std::endl
      << "sbottomPlusW->getEffectiveSquarkMass() = "
      << sbottomPlusW->getEffectiveSquarkMass();
      std::cout << std::endl;

      fullCascadeType::squarkDirectlyToElectroweakType::supType
      testStXOneFullCascade( testInputHandler,
                             stopOne,
                             7,
                             testEwinoCascadeSet,
                             cascadesForSevenTev->getCascadeSet( stopOne,
                                                                 charginoOne,
                                                                 stopMinusW ),
                             cascadesForSevenTev->getCascadeSet( stopOne,
                                                                 wPlus,
                                                                 sbottomOne,
                                                              sbottomPlusW ) );
      std::cout
      << std::endl
      << "testStXOneFullCascade set as stop_1 -> chi^+_1 ( "
      << testStXOneFullCascade.getAsString() << " )";
      std::cout << std::endl;
      if( testStXOneFullCascade.isOpen() )
      {
        std::cout
        << std::endl
        << "testStXOneFullCascade is open.";
        std::cout << std::endl;
      }
      else
      {
        std::cout
        << std::endl
        << "testStXOneFullCascade is not open?!";
        std::cout << std::endl;
      }
      testSmParticleExclusion.clear();
      std::cout
      << std::endl
      << "BR with no excluded SM fermions = "
      << testStXOneFullCascade.getBrToEwino( &testSmParticleExclusion );
      testSmParticleExclusion.push_back( CppSLHA::PDG_code::up );
      testSmParticleExclusion.push_back( -(CppSLHA::PDG_code::up) );
      std::cout
      << std::endl
      << "BR excluding ups & antiups in final state = "
      << testStXOneFullCascade.getBrToEwino( &testSmParticleExclusion );
      std::cout << std::endl;
      std::cout
      << std::endl
      << ".getAcceptance( [particle], &testAcceptanceCutSet, 0, 1, 0, 0, 0 )"
      << " = " << testStXOneFullCascade.getAcceptance( true,
                                                       &testAcceptanceCutSet,
                                                       0,
                                                       1,
                                                       0,
                                                       0,
                                                       0 );
      std::cout
      << std::endl
      << ".getAcceptance("
      << " [antiparticle], &testAcceptanceCutSet, 0, 1, 0, 0, 0 ) = "
      << testStXOneFullCascade.getAcceptance( false,
                                              &testAcceptanceCutSet,
                                              0,
                                              1,
                                              0,
                                              0,
                                              0 );
      std::cout
      << std::endl
      << ".getAcceptance("
      << " [antiparticle], &testAcceptanceCutSet, 0, 0, 1, 0, 0 ) = "
      << testStXOneFullCascade.getAcceptance( false,
                                              &testAcceptanceCutSet,
                                              0,
                                              0,
                                              1,
                                              0,
                                              0 );
      std::cout
      << std::endl
      << "testStXOneFullCascade seems OK (large part of the example BR is from"
      << " taus, so a lighter squark mass has a big effect).";
      std::cout << std::endl;
      std::cout << std::endl;
      electroweakCascadeSet*
      testBosonCascadeSet( cascadesForSevenTev->getCascadeSet( sbottomTwo,
                                                               wPlus,
                                                               stopOne ) );
      fullCascadeType::squarkByBosonToCompound testSbWMFullCascade;
      testSbWMFullCascade.setProperties( sbottomTwo,
                                         testBosonCascadeSet,
                                         &testStXOneFullCascade );
      std::cout
      << std::endl
      << "testSbWMFullCascade set with sbottom_2 -> W^- + stop_1,"
      << " stop_1 -> chargino_1 ( " << testSbWMFullCascade.getAsString()
      << " )";
      std::cout << std::endl;
      testSmParticleExclusion.clear();
      std::cout
      << std::endl
      << "BR with no excluded SM fermions = "
      << testSbWMFullCascade.getBrToEwino( &testSmParticleExclusion );
      testSmParticleExclusion.push_back( CppSLHA::PDG_code::down );
      std::cout
      << std::endl
      << "BR excluding downs (but not antidowns) in final state = "
      << testSbWMFullCascade.getBrToEwino( &testSmParticleExclusion );
      std::cout << std::endl;
      std::cout
      << std::endl
      << ".getAcceptance( [particle], &testAcceptanceCutSet, 0, 1, 1, 0, 0 )"
      << " = " << testSbWMFullCascade.getAcceptance( true,
                                                     &testAcceptanceCutSet,
                                                     0,
                                                     1,
                                                     1,
                                                     0,
                                                     0 );
      std::cout
      << std::endl
      << ".getAcceptance("
      << " [antiparticle], &testAcceptanceCutSet, 0, 1, 1, 0, 0 ) = "
      << testSbWMFullCascade.getAcceptance( false,
                                            &testAcceptanceCutSet,
                                            0,
                                            1,
                                            1,
                                            0,
                                            0 );
      std::cout
      << std::endl
      << "testSbWMFullCascade seems OK (~10% for W^- to pass e^-, ~4%"
      << " for chargino as shown by testStXOneFullCascade).";
      std::cout << std::endl;
      std::cout << std::endl;
      fullCascadeType::gluinoOrNeutralinoToCompound testGoMFullCascade;
      testGoMFullCascade.setProperties( gluinoPointer,
                                        &testSdXOneFullCascade,
                                        cascadesForSevenTev );
      std::cout
      << std::endl
      << "testGoMFullCascade set with testSdXOneFullCascade (pointer is "
      << &testSdXOneFullCascade
      << "), so is gluino => sdown_L => chargino_1 ( "
      << testGoMFullCascade.getAsString() << " )";
      testSmParticleExclusion.clear();
      std::cout
      << std::endl
      << "BR with no excluded SM fermions = "
      << testGoMFullCascade.getBrToEwino( &testSmParticleExclusion );
      testSmParticleExclusion.push_back( CppSLHA::PDG_code::down );
      std::cout
      << std::endl
      << "BR excluding downs (but not antidowns) in final state = "
      << testGoMFullCascade.getBrToEwino( &testSmParticleExclusion );
      std::cout << std::endl;
      std::cout
      << std::endl
      << ".getAcceptance( [particle], &testAcceptanceCutSet, 0, 1, 0, 0, 0 )"
      << " = " << testGoMFullCascade.getAcceptance( true,
                                                    &testAcceptanceCutSet,
                                                    0,
                                                    1,
                                                    0,
                                                    0,
                                                    0 );
      std::cout
      << std::endl
      << ".getAcceptance("
      << " [antiparticle], &testAcceptanceCutSet, 0, 1, 0, 0, 0 ) = "
      << testGoMFullCascade.getAcceptance( false,
                                           &testAcceptanceCutSet,
                                           0,
                                           1,
                                           0,
                                           0,
                                           0 );
      std::cout
      << std::endl
      << "testGoMFullCascade seems OK.";
      std::cout << std::endl;
      std::cout << std::endl;
      testGoMFullCascade.setProperties( gluinoPointer,
                                        &testSbWMFullCascade,
                                        cascadesForSevenTev );
      basicStuff->getReadier()->readyObserversForNewPoint();
      std::cout
      << std::endl
      << "testGoMFullCascade set with testSbWMFullCascade"
      << " (& readiers reset) ( "
      << testGoMFullCascade.getAsString() << " )";
      std::cout << std::endl;
      testSmParticleExclusion.clear();
      std::cout
      << std::endl
      << "BR with no excluded SM fermions = "
      << testGoMFullCascade.getBrToEwino( &testSmParticleExclusion );
      testSmParticleExclusion.push_back( CppSLHA::PDG_code::down );
      std::cout
      << std::endl
      << "BR excluding downs (but not antidowns) in final state = "
      << testGoMFullCascade.getBrToEwino( &testSmParticleExclusion );
      std::cout << std::endl;
      std::cout
      << std::endl
      << ".getAcceptance( [particle], &testAcceptanceCutSet, 0, 1, 1, 0, 0 )"
      << " = " << testGoMFullCascade.getAcceptance( true,
                                                    &testAcceptanceCutSet,
                                                    0,
                                                    1,
                                                    1,
                                                    0,
                                                    0 );
      std::cout
      << std::endl
      << ".getAcceptance("
      << " [antiparticle], &testAcceptanceCutSet, 0, 1, 1, 0, 0 ) = "
      << testGoMFullCascade.getAcceptance( false,
                                           &testAcceptanceCutSet,
                                           0,
                                           1,
                                           1,
                                           0,
                                           0 );
      std::cout
      << std::endl
      << "testGoMFullCascade seems OK.";
      std::cout << std::endl;
      std::cout << std::endl;


      // now the convolutedly long cascade:

      basicStuff->getSlha()->read_file( "convoluted_spectrum.out" );
      basicStuff->getReadier()->readyObserversForNewPoint();

      std::cout
      << std::endl
      << "stopMinusW->getEffectiveSquarkMass() = "
      << stopMinusW->getEffectiveSquarkMass()
      << std::endl
      << "sbottomPlusW->getEffectiveSquarkMass() = "
      << sbottomPlusW->getEffectiveSquarkMass()
      << std::endl
      << "stopOne->get_absolute_mass() = "
      << stopOne->get_absolute_mass()
      << std::endl
      << "sbottomOne->get_absolute_mass() = "
      << sbottomOne->get_absolute_mass();
      std::cout << std::endl;
      testEwinoCascadeSet = cascadesForSevenTev->getCascadeSet( stopOne,
                                                               neutralinoOne );
      fullCascadeType::squarkDirectlyToElectroweakType::supType
      testStWNOneFullCascade( testInputHandler,
                              stopOne,
                              7,
                              testEwinoCascadeSet,
                              cascadesForSevenTev->getCascadeSet( stopOne,
                                                                 neutralinoOne,
                                                                  stopMinusW ),
                              cascadesForSevenTev->getCascadeSet( stopOne,
                                                                  wPlus,
                                                                  sbottomOne,
                                                              sbottomPlusW ) );
      std::cout
      << std::endl
      << "testStWNOneFullCascade set as stop_1 -> chi^0_1 + W^+ + b ( "
      << testStWNOneFullCascade.getAsString() << " )";
      std::cout << std::endl;
      if( testStWNOneFullCascade.isOpen() )
      {
        std::cout
        << std::endl
        << "testStWNOneFullCascade is open.";
        std::cout << std::endl;
      }
      else
      {
        std::cout
        << std::endl
        << "testStWNOneFullCascade is not open?!";
        std::cout << std::endl;
      }
      testSmParticleExclusion.clear();
      std::cout
      << std::endl
      << "BR with no excluded SM fermions = "
      << testStWNOneFullCascade.getBrToEwino( &testSmParticleExclusion );
      testSmParticleExclusion.push_back( CppSLHA::PDG_code::down );
      std::cout
      << std::endl
      << "BR excluding downs (but not antidowns) in final state = "
      << testStWNOneFullCascade.getBrToEwino( &testSmParticleExclusion );
      std::cout << std::endl;
      std::cout
      << std::endl
      << ".getAcceptance( [particle], &testAcceptanceCutSet, 0, 0, 1, 0, 0 )"
      << " = " << testStWNOneFullCascade.getAcceptance( true,
                                                        &testAcceptanceCutSet,
                                                        0,
                                                        0,
                                                        1,
                                                        0,
                                                        0 );
      std::cout
      << std::endl
      << ".getAcceptance("
      << " [antiparticle], &testAcceptanceCutSet, 0, 0, 1, 0, 0 ) = "
      << testStWNOneFullCascade.getAcceptance( false,
                                               &testAcceptanceCutSet,
                                               0,
                                               0,
                                               1,
                                               0,
                                               0 );
      std::cout
      << std::endl
      << "testStWNOneFullCascade seems OK.";
      std::cout
      << std::endl
      << std::endl
      << "testSuperlongFullCascade should only work if negligibleBr has been"
      << " set low enough (~10^-9).";
      std::cout << std::endl;
      std::cout << std::endl;
      testSbWMFullCascade.setProperties( sbottomOne,
                                         cascadesForSevenTev->getCascadeSet(
                                                                    sbottomOne,
                                                                         wPlus,
                                                                     stopOne ),
                                         &testStWNOneFullCascade );
      std::cout
      << std::endl
      << "testSbWMFullCascade set with sbottom_2 -> W^- + stop_1,"
      << " stop_1 -> W^+ + b + neutralino_1 ( "
      << testSbWMFullCascade.getAsString() << " )";
      std::cout << std::endl;
      testSmParticleExclusion.clear();
      std::cout
      << std::endl
      << "BR with no excluded SM fermions = "
      << testSbWMFullCascade.getBrToEwino( &testSmParticleExclusion );
      testSmParticleExclusion.push_back( CppSLHA::PDG_code::down );
      std::cout
      << std::endl
      << "BR excluding downs (but not antidowns) in final state = "
      << testSbWMFullCascade.getBrToEwino( &testSmParticleExclusion );
      std::cout << std::endl;
      std::cout
      << std::endl
      << ".getAcceptance( [particle], &testAcceptanceCutSet, 0, 1, 1, 0, 0 )"
      << " = " << testSbWMFullCascade.getAcceptance( true,
                                                     &testAcceptanceCutSet,
                                                     0,
                                                     1,
                                                     1,
                                                     0,
                                                     0 );
      std::cout
      << std::endl
      << ".getAcceptance("
      << " [antiparticle], &testAcceptanceCutSet, 0, 1, 1, 0, 0 ) = "
      << testSbWMFullCascade.getAcceptance( false,
                                            &testAcceptanceCutSet,
                                            0,
                                            1,
                                            1,
                                            0,
                                            0 );
      std::cout
      << std::endl
      << "testSbWMFullCascade seems OK (~10% for W^- to pass e^-, ~10% for W^+"
      << " to pass e^+ as shown by testStWNOneFullCascade).";
      std::cout << std::endl;
      std::cout << std::endl;
      fullCascadeType::gluinoOrNeutralinoToCompound testNTwoMFullCascade;
      testNTwoMFullCascade.setProperties( neutralinoTwo,
                                          &testSbWMFullCascade,
                                          cascadesForSevenTev );
      std::cout
      << std::endl
      << "testNTwoMFullCascade set to be neutralino_2 -> testSbWMFullCascade"
      << " ( " << testNTwoMFullCascade.getAsString() << " )";
      std::cout << std::endl;
      testSmParticleExclusion.clear();
      std::cout
      << std::endl
      << "BR with no excluded SM fermions = "
      << testNTwoMFullCascade.getBrToEwino( &testSmParticleExclusion );
      testSmParticleExclusion.push_back( CppSLHA::PDG_code::down );
      std::cout
      << std::endl
      << "BR excluding downs (but not antidowns) in final state = "
      << testNTwoMFullCascade.getBrToEwino( &testSmParticleExclusion );
      std::cout << std::endl;
      std::cout
      << std::endl
      << ".getAcceptance( [particle], &testAcceptanceCutSet, 0, 1, 1, 0, 0 )"
      << " = " << testNTwoMFullCascade.getAcceptance( true,
                                                      &testAcceptanceCutSet,
                                                      0,
                                                      1,
                                                      1,
                                                      0,
                                                      0 );
      std::cout
      << std::endl
      << ".getAcceptance("
      << " [antiparticle], &testAcceptanceCutSet, 0, 1, 1, 0, 0 ) = "
      << testNTwoMFullCascade.getAcceptance( false,
                                             &testAcceptanceCutSet,
                                             0,
                                             1,
                                             1,
                                             0,
                                             0 );
      std::cout
      << std::endl
      << "testNTwoMFullCascade seems OK.";
      std::cout << std::endl;
      std::cout << std::endl;
      fullCascadeType::squarkByJetToCompoundType::sdownType
      testSbTwoJMFullCascade;
      testSbTwoJMFullCascade.setProperties( sbottomTwo,
                                            &testNTwoMFullCascade,
                                            cascadesForSevenTev );
      std::cout
      << std::endl
      << "testSbTwoJMFullCascade set to be sbottom_2 -> testNTwoMFullCascade"
      << " ( " << testSbTwoJMFullCascade.getAsString() << " )";
      std::cout << std::endl;
      testSmParticleExclusion.clear();
      std::cout
      << std::endl
      << "BR with no excluded SM fermions = "
      << testSbTwoJMFullCascade.getBrToEwino( &testSmParticleExclusion );
      testSmParticleExclusion.push_back( CppSLHA::PDG_code::down );
      std::cout
      << std::endl
      << "BR excluding downs (but not antidowns) in final state = "
      << testSbTwoJMFullCascade.getBrToEwino( &testSmParticleExclusion );
      std::cout << std::endl;
      std::cout
      << std::endl
      << ".getAcceptance( [particle], &testAcceptanceCutSet, 0, 1, 1, 0, 0 )"
      << " = " << testSbTwoJMFullCascade.getAcceptance( true,
                                                        &testAcceptanceCutSet,
                                                        0,
                                                        1,
                                                        1,
                                                        0,
                                                        0 );
      std::cout
      << std::endl
      << ".getAcceptance("
      << " [antiparticle], &testAcceptanceCutSet, 0, 1, 1, 0, 0 ) = "
      << testSbTwoJMFullCascade.getAcceptance( false,
                                               &testAcceptanceCutSet,
                                               0,
                                               1,
                                               1,
                                               0,
                                               0 );
      std::cout
      << std::endl
      << "testSbTwoJMFullCascade seems OK.";
      std::cout << std::endl;
      std::cout << std::endl;
      testGoMFullCascade.setProperties( gluinoPointer,
                                        &testSbTwoJMFullCascade,
                                        cascadesForSevenTev );
      std::cout
      << std::endl
      << "testGoMFullCascade set to be gluino -> testSbTwoJMFullCascade ( "
      << testGoMFullCascade.getAsString() << " )";
      std::cout << std::endl;
      testSmParticleExclusion.clear();
      std::cout
      << std::endl
      << "BR with no excluded SM fermions = "
      << testGoMFullCascade.getBrToEwino( &testSmParticleExclusion );
      testSmParticleExclusion.push_back( CppSLHA::PDG_code::down );
      std::cout
      << std::endl
      << "BR excluding downs (but not antidowns) in final state = "
      << testGoMFullCascade.getBrToEwino( &testSmParticleExclusion );
      std::cout << std::endl;
      std::cout
      << std::endl
      << ".getAcceptance( [particle], &testAcceptanceCutSet, 0, 1, 1, 0, 0 )"
      << " = " << testGoMFullCascade.getAcceptance( true,
                                                    &testAcceptanceCutSet,
                                                    0,
                                                    1,
                                                    1,
                                                    0,
                                                    0 );
      std::cout
      << std::endl
      << ".getAcceptance("
      << " [antiparticle], &testAcceptanceCutSet, 0, 1, 1, 0, 0 ) = "
      << testGoMFullCascade.getAcceptance( false,
                                           &testAcceptanceCutSet,
                                           0,
                                           1,
                                           1,
                                           0,
                                           0 );
      std::cout
      << std::endl
      << "testGoMFullCascade seems OK.";
      std::cout << std::endl;
      std::cout << std::endl;
      fullCascadeType::squarkByJetToCompoundType::supType testSupJMFullCascade;
      testSupJMFullCascade.setProperties( testInputHandler->getSupL(),
                                          &testGoMFullCascade,
                                          cascadesForSevenTev );
      std::cout
      << std::endl
      << "testSupJMFullCascade set to be sup_L -> testGoMFullCascade ( "
      << testSupJMFullCascade.getAsString() << " )";
      std::cout << std::endl;
      testSmParticleExclusion.clear();
      std::cout
      << std::endl
      << "BR with no excluded SM fermions = "
      << testSupJMFullCascade.getBrToEwino( &testSmParticleExclusion );
      testSmParticleExclusion.push_back( CppSLHA::PDG_code::down );
      std::cout
      << std::endl
      << "BR excluding downs (but not antidowns) in final state = "
      << testSupJMFullCascade.getBrToEwino( &testSmParticleExclusion );
      std::cout << std::endl;
      std::cout
      << std::endl
      << ".getAcceptance( [particle], &testAcceptanceCutSet, 0, 1, 1, 0, 0 )"
      << " = " << testSupJMFullCascade.getAcceptance( true,
                                                      &testAcceptanceCutSet,
                                                      0,
                                                      1,
                                                      1,
                                                      0,
                                                      0 );
      std::cout
      << std::endl
      << ".getAcceptance("
      << " [antiparticle], &testAcceptanceCutSet, 0, 1, 1, 0, 0 ) = "
      << testSupJMFullCascade.getAcceptance( false,
                                             &testAcceptanceCutSet,
                                             0,
                                             1,
                                             1,
                                             0,
                                             0 );
      std::cout
      << std::endl
      << "testSupJMFullCascade seems OK.";
      std::cout << std::endl;
      std::cout << std::endl;


      // go back to SPS1a:
      basicStuff->getSlha()->read_file( "SPS1a_spectrum.out" );
      basicStuff->getReadier()->readyObserversForNewPoint();


      fullCascadeSetFactory testFullCascadeSetFactory( testInputHandler,
                                              &testElectroweakCascadeHandler );
      std::cout
      << std::endl
      << "testFullCascadeSetFactory constructed without problem.";
      std::cout << std::endl;
      fullCascadeSet*
      testFullCascadeSetFromFactory(
                                   testFullCascadeSetFactory.getFullCascadeSet(
                                                 testInputHandler->getSdownL(),
                                                                         7 ) );

      std::cout
      << std::endl
      << "testFullCascadeSetFromFactory = " << testFullCascadeSetFromFactory;
      std::cout << std::endl;
      std::vector< fullCascade* >*
      openCascades( testFullCascadeSetFromFactory->getOpenCascades() );
      testSmParticleExclusion.clear();
      testSmParticleExclusion.push_back( CppSLHA::PDG_code::top );
      testSmParticleExclusion.push_back( -(CppSLHA::PDG_code::top) );
      std::cout
      << std::endl
      << "has " << openCascades->size()
      << " open cascades. 1 is squarkToEwino, 2 is gluinoToEwino, 3 is"
      << " squarkToSquarkThenMore, 4 is squarkToGauginoThenMore, 5 is"
      << " gluinoOrElectroweakinoToSquarkThenMore. BR is for"
      << " ->getBrToEwino( &testSmParticleExclusion ), which excludes top &"
      << " antitop, \"acc\" is for ->getAcceptance("
      << " [particle], &testAcceptanceCutSet, 0, 0, 0, 1, 0 ).";
      for( std::vector< fullCascade* >::iterator
           cascadeIterator( openCascades->begin() );
           openCascades->end() > cascadeIterator;
           ++cascadeIterator )
      {
        std::cout
        << std::endl
        << *((*cascadeIterator)->getInitialSparticle()->get_name())
        << " => " << (*cascadeIterator)->getColorfulCascadeType()
        << " => "
        << *((*cascadeIterator)->getElectroweakinoAtEndOfScoloredness(
                                                                 )->get_name())
        << ", BR = "
        << (*cascadeIterator)->getBrToEwino( &testSmParticleExclusion )
        << ", acc = "
        << (*cascadeIterator)->getAcceptance( true,
                                              &testAcceptanceCutSet,
                                              0,
                                              0,
                                              0,
                                              1,
                                              0 );
      }
      testFullCascadeSetFromFactory
      = testFullCascadeSetFactory.getFullCascadeSet(
                                                 testInputHandler->getGluino(),
                                                     7 );
    std::cout
    << std::endl
    << "testFullCascadeSetFromFactory = " << testFullCascadeSetFromFactory;
    std::cout << std::endl;
    openCascades = testFullCascadeSetFromFactory->getOpenCascades();
    std::cout
    << std::endl
    << "has " << openCascades->size()
    << " open cascades. 1 is squarkToEwino, 2 is gluinoToEwino, 3 is"
    << " squarkToSquarkThenMore, 4 is squarkToGauginoThenMore, 5 is"
    << " gluinoOrElectroweakinoToSquarkThenMore. BR is for"
    << " ->getBrToEwino( &testSmParticleExclusion ), which"
    << " excludes top & antitop, \"acc\" is for ->getAcceptance("
    << " [particle], &testAcceptanceCutSet, 0, 0, 0, 1, 0 ).";
    for( std::vector< fullCascade* >::iterator
         cascadeIterator( openCascades->begin() );
         openCascades->end() > cascadeIterator;
         ++cascadeIterator )
    {
      std::cout
      << std::endl
      << *((*cascadeIterator)->getInitialSparticle()->get_name())
      << " => " << (*cascadeIterator)->getColorfulCascadeType()
      << " => "
      << *((*cascadeIterator)->getElectroweakinoAtEndOfScoloredness(
                                                               )->get_name())
      << ", BR = "
      << (*cascadeIterator)->getBrToEwino( &testSmParticleExclusion )
      << ", acc = "
      << (*cascadeIterator)->getAcceptance( true,
                                            &testAcceptanceCutSet,
                                            0,
                                            0,
                                            0,
                                            1,
                                            0 );
    }
    std::cout
    << std::endl
    << "remember that gluino -> (squark ->) chargino sums over both charges";
    std::cout << std::endl;

    std::cout
    << std::endl
    << "changing spectrum.";
    std::cout << std::endl;
    basicStuff->getSlha()->read_file( "SPS2_spectrum.out" );
    basicStuff->getReadier()->readyObserversForNewPoint();
    testFullCascadeSetFromFactory
    = testFullCascadeSetFactory.getFullCascadeSet(
                                                 testInputHandler->getSdownL(),
                                                   7 );

    std::cout
    << std::endl
    << "testFullCascadeSetFromFactory = " << testFullCascadeSetFromFactory;
    std::cout << std::endl;
    openCascades = testFullCascadeSetFromFactory->getOpenCascades();
    testSmParticleExclusion.clear();
    testSmParticleExclusion.push_back( CppSLHA::PDG_code::top );
    testSmParticleExclusion.push_back( -(CppSLHA::PDG_code::top) );
    std::cout
    << std::endl
    << "has " << openCascades->size()
    << " open cascades. 1 is squarkToEwino, 2 is gluinoToEwino, 3 is"
    << " squarkToSquarkThenMore, 4 is squarkToGauginoThenMore, 5 is"
    << " gluinoOrElectroweakinoToSquarkThenMore. BR is for"
    << " ->getBrToEwino( &testSmParticleExclusion ), which excludes top &"
    << " antitop, \"acc\" is for ->getAcceptance("
    << " [particle], &testAcceptanceCutSet, 0, 0, 0, 1, 0 ).";
    for( std::vector< fullCascade* >::iterator
         cascadeIterator( openCascades->begin() );
         openCascades->end() > cascadeIterator;
         ++cascadeIterator )
    {
      std::cout
      << std::endl
      << *((*cascadeIterator)->getInitialSparticle()->get_name())
      << " => " << (*cascadeIterator)->getColorfulCascadeType()
      << " => "
      << *((*cascadeIterator)->getElectroweakinoAtEndOfScoloredness(
                                                               )->get_name())
      << ", BR = "
      << (*cascadeIterator)->getBrToEwino( &testSmParticleExclusion )
      << ", acc = "
      << (*cascadeIterator)->getAcceptance( true,
                                            &testAcceptanceCutSet,
                                            0,
                                            0,
                                            0,
                                            1,
                                            0 );
    }
    std::cout
    << std::endl
    << "remember that squark -> gluino -> (squark ->) chargino sums over both"
    << " charges";
    std::cout << std::endl;
    testFullCascadeSetFromFactory
    = testFullCascadeSetFactory.getFullCascadeSet(
                                               testInputHandler->getGluino(),
                                                   7 );
    std::cout
    << std::endl
    << "testFullCascadeSetFromFactory = " << testFullCascadeSetFromFactory;
    std::cout << std::endl;
    openCascades = testFullCascadeSetFromFactory->getOpenCascades();
    std::cout
    << std::endl
    << "has " << openCascades->size()
    << " open cascades. 1 is squarkToEwino, 2 is gluinoToEwino, 3 is"
    << " squarkToSquarkThenMore, 4 is squarkToGauginoThenMore, 5 is"
    << " gluinoOrElectroweakinoToSquarkThenMore. BR is for"
    << " ->getBrToEwino( &testSmParticleExclusion ), which excludes top &"
    << " antitop, \"acc\" is for ->getAcceptance("
    << " [particle], &testAcceptanceCutSet, 0, 0, 0, 1, 0 ).";
    for( std::vector< fullCascade* >::iterator
         cascadeIterator( openCascades->begin() );
         openCascades->end() > cascadeIterator;
         ++cascadeIterator )
    {
      std::cout
      << std::endl
      << *((*cascadeIterator)->getInitialSparticle()->get_name())
      << " => " << (*cascadeIterator)->getColorfulCascadeType()
      << " => "
      << *((*cascadeIterator)->getElectroweakinoAtEndOfScoloredness(
                                                                 )->get_name())
      << ", BR = "
      << (*cascadeIterator)->getBrToEwino( &testSmParticleExclusion )
      << ", acc = "
      << (*cascadeIterator)->getAcceptance( true,
                                            &testAcceptanceCutSet,
                                            0,
                                            0,
                                            0,
                                            1,
                                            0 );
    }
    std::cout
    << std::endl
    << "remember that gluino -> (squark ->) chargino sums over both charges";
    std::cout << std::endl;
    std::cout << std::endl;
    basicStuff->getSlha()->read_file( "SPS1a_spectrum.out" );
    basicStuff->getReadier()->readyObserversForNewPoint();
  }


  protected:
    basicStuffTesting* basicStuff;
    inputTesting* inputTester;
    inputHandler const* const testInputHandler;
    acceptanceCutSet testAcceptanceCutSet;
  };

}  // end of LHC_FASER namespace

#endif /* FULL_CASCADE_TESTING_HPP_ */
