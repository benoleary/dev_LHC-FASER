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
      electroweakCascadeHandler
      testElectroweakCascadeHandler( testInputHandler );
      electroweakCascadeSet*
      testEwinoCascadeSet(
        testElectroweakCascadeHandler.getElectroweakCascadesForOneBeamEnergy( 7
                                                              )->getCascadeSet(
                                                 testInputHandler->getSdownL(),
                                      testInputHandler->getNeutralinoTwo() ) );
      sxFullCascade testSxFullCascade;
      testSxFullCascade.setProperties( testInputHandler,
                                       testInputHandler->getSdownL(),
                                       7,
                                       testEwinoCascadeSet );
      std::cout
      << std::endl
      << "testSxFullCascade set as sdown_L -> chi^0_2";
      std::cout << std::endl;
      if( testSxFullCascade.isOpen() )
      {
        std::cout
        << std::endl
        << "testSxFullCascade is open.";
        std::cout << std::endl;
      }
      else
      {
        std::cout
        << std::endl
        << "testSxFullCascade is not open?!";
        std::cout << std::endl;
      }
      std::list< int > testSmParticleExclusion;
      std::cout
      << std::endl
      << "BR with no excluded SM fermions = "
      << testSxFullCascade.getBrToEwino( &testSmParticleExclusion );
      testSmParticleExclusion.push_back( CppSLHA::PDG_code::down );
      std::cout
      << std::endl
      << "BR excluding downs (but not antidowns) in final state = "
      << testSxFullCascade.getBrToEwino( &testSmParticleExclusion );
      std::cout << std::endl;
      std::cout
      << std::endl
      << ".getAcceptance( [particle], &testAcceptanceCutSet, 0, 1, 1, 0, 0 )"
      << " = " << testSxFullCascade.getAcceptance( true,
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
      << testSxFullCascade.getAcceptance( false,
                                          &testAcceptanceCutSet,
                                          0,
                                          1,
                                          1,
                                          0,
                                          0 );
      std::cout
      << std::endl
      << "testSxFullCascade seems OK.";
      std::cout << std::endl;
      std::cout << std::endl;
      testEwinoCascadeSet
      = testElectroweakCascadeHandler.getElectroweakCascadesForOneBeamEnergy( 7
                                                              )->getCascadeSet(
                                                 testInputHandler->getSdownL(),
                                          testInputHandler->getCharginoOne() );
      testSxFullCascade.setProperties( testInputHandler,
                                       testInputHandler->getSdownL(),
                                       7,
                                       testEwinoCascadeSet );
      std::cout
      << std::endl
      << "testSxFullCascade set as sdown_L -> chi^-_1";
      std::cout << std::endl;
      if( testSxFullCascade.isOpen() )
      {
        std::cout
        << std::endl
        << "testSxFullCascade is open.";
        std::cout << std::endl;
      }
      else
      {
        std::cout
        << std::endl
        << "testSxFullCascade is not open?!";
        std::cout << std::endl;
      }
      testSmParticleExclusion.clear();
      std::cout
      << std::endl
      << "BR with no excluded SM fermions = "
      << testSxFullCascade.getBrToEwino( &testSmParticleExclusion );
      testSmParticleExclusion.push_back( CppSLHA::PDG_code::up );
      testSmParticleExclusion.push_back( -(CppSLHA::PDG_code::up) );
      std::cout
      << std::endl
      << "BR excluding ups & antiups in final state = "
      << testSxFullCascade.getBrToEwino( &testSmParticleExclusion );
      std::cout << std::endl;
      std::cout
      << std::endl
      << ".getAcceptance( [particle], &testAcceptanceCutSet, 0, 1, 0, 0, 0 )"
      << " = " << testSxFullCascade.getAcceptance( true,
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
      << testSxFullCascade.getAcceptance( false,
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
      << testSxFullCascade.getAcceptance( false,
                                          &testAcceptanceCutSet,
                                          0,
                                          0,
                                          1,
                                          0,
                                          0 );
      std::cout
      << std::endl
      << "testSxFullCascade again seems OK.";
      std::cout << std::endl;
      std::cout << std::endl;
      testEwinoCascadeSet
      = testElectroweakCascadeHandler.getElectroweakCascadesForOneBeamEnergy( 7
                                                              )->getCascadeSet(
                                                 testInputHandler->getSdownL(),
                                          testInputHandler->getCharginoTwo() );
      testSxFullCascade.setProperties( testInputHandler,
                                       testInputHandler->getSdownL(),
                                       7,
                                       testEwinoCascadeSet );
      std::cout
      << std::endl
      << "testSxFullCascade set as sdown_L -> chi^-_2";
      std::cout << std::endl;
      if( testSxFullCascade.isOpen() )
      {
        std::cout
        << std::endl
        << "testSxFullCascade is open.";
        std::cout << std::endl;
      }
      else
      {
        std::cout
        << std::endl
        << "testSxFullCascade is not open?!";
        std::cout << std::endl;
      }
      testSmParticleExclusion.clear();
      std::cout
      << std::endl
      << "BR with no excluded SM fermions = "
      << testSxFullCascade.getBrToEwino( &testSmParticleExclusion );
      testSmParticleExclusion.push_back( CppSLHA::PDG_code::up );
      testSmParticleExclusion.push_back( -(CppSLHA::PDG_code::up) );
      std::cout
      << std::endl
      << "BR excluding ups & antiups in final state = "
      << testSxFullCascade.getBrToEwino( &testSmParticleExclusion );
      std::cout << std::endl;
      std::cout
      << std::endl
      << ".getAcceptance( [particle], &testAcceptanceCutSet, 0, 0, 0, 0, 0 )"
      << " = " << testSxFullCascade.getAcceptance( true,
                                                   &testAcceptanceCutSet,
                                                   0,
                                                   0,
                                                   0,
                                                   0,
                                                   0 );
      std::cout
      << std::endl
      << "testSxFullCascade again seems OK.";
      std::cout << std::endl;
      std::cout << std::endl;


      testEwinoCascadeSet
      = testElectroweakCascadeHandler.getElectroweakCascadesForOneBeamEnergy( 7
                                                              )->getCascadeSet(
                                                 testInputHandler->getGluino(),
                                        testInputHandler->getNeutralinoTwo() );
      gxFullCascade testGxFullCascade;
      testGxFullCascade.setProperties( testInputHandler,
                                       7,
                                       testEwinoCascadeSet );
      std::cout
      << std::endl
      << "testGxFullCascade set as gluino -> chi^0_2";
      std::cout << std::endl;
      if( testGxFullCascade.isOpen() )
      {
        std::cout
        << std::endl
        << "testGxFullCascade is open.";
        std::cout << std::endl;
      }
      else
      {
        std::cout
        << std::endl
        << "testGxFullCascade is not open?!";
        std::cout << std::endl;
      }
      testSmParticleExclusion.clear();
      std::cout
      << std::endl
      << "BR with no excluded SM fermions = "
      << testGxFullCascade.getBrToEwino( &testSmParticleExclusion );
      testSmParticleExclusion.push_back( CppSLHA::PDG_code::down );
      std::cout
      << std::endl
      << "BR excluding downs (but not antidowns) in final state = "
      << testGxFullCascade.getBrToEwino( &testSmParticleExclusion );
      std::cout
      << std::endl
      << "testGxFullCascade is open, & gives non-zero acceptances, but the"
      << " BR is so low that it's snapped to zero.";
      std::cout << std::endl;
      std::cout
      << std::endl
      << ".getAcceptance( [particle], &testAcceptanceCutSet, 0, 1, 1, 0, 0 )"
      << " = " << testGxFullCascade.getAcceptance( true,
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
      << testGxFullCascade.getAcceptance( false,
                                          &testAcceptanceCutSet,
                                          0,
                                          1,
                                          1,
                                          0,
                                          0 );
      std::cout << std::endl;


      gjsxFullCascade testGjsxFullCascade;
      testGjsxFullCascade.setProperties( &testSxFullCascade );
      std::cout
      << std::endl
      << "testGjsxFullCascade set with testSxFullCascade (pointer is "
      << &testSxFullCascade << "), so is gluino => sdown_L => chargino_1";
      testSmParticleExclusion.clear();
      std::cout
      << std::endl
      << "BR with no excluded SM fermions = "
      << testGjsxFullCascade.getBrToEwino( &testSmParticleExclusion );
      testSmParticleExclusion.push_back( CppSLHA::PDG_code::down );
      std::cout
      << std::endl
      << "BR excluding downs (but not antidowns) in final state = "
      << testGjsxFullCascade.getBrToEwino( &testSmParticleExclusion );
      std::cout << std::endl;
      std::cout
      << std::endl
      << ".getAcceptance( [particle], &testAcceptanceCutSet, 0, 1, 0, 0, 0 )"
      << " = " << testGjsxFullCascade.getAcceptance( true,
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
      << testGjsxFullCascade.getAcceptance( false,
                                            &testAcceptanceCutSet,
                                            0,
                                            1,
                                            0,
                                            0,
                                            0 );
      std::cout
      << std::endl
      << "testGjsxFullCascade seems OK.";
      std::cout << std::endl;
      std::cout << std::endl;


      testEwinoCascadeSet
      = testElectroweakCascadeHandler.getElectroweakCascadesForOneBeamEnergy( 7
                                                              )->getCascadeSet(
                                                testInputHandler->getStopOne(),
                                          testInputHandler->getCharginoOne() );
      sxFullCascade stopOneCharginoOneCascade;
      stopOneCharginoOneCascade.setProperties( testInputHandler,
                                               testInputHandler->getStopOne(),
                                               7,
                                               testEwinoCascadeSet );
      electroweakCascadeSet*
      testBosonCascadeSet(
        testElectroweakCascadeHandler.getElectroweakCascadesForOneBeamEnergy( 7
                                                              )->getCascadeSet(
                                             testInputHandler->getSbottomTwo(),
                                                  testInputHandler->getWPlus(),
                                            testInputHandler->getStopOne() ) );
      sbsxFullCascade testSbsxFullCascade;
      testSbsxFullCascade.setProperties( testInputHandler->getSbottomTwo(),
                                         testBosonCascadeSet,
                                         &stopOneCharginoOneCascade );
      std::cout
      << std::endl
      << "testSbsxFullCascade set with sbottom_2 -> W^- stop_1 -> chargino_1";
      std::cout << std::endl;
      testSmParticleExclusion.clear();
      std::cout
      << std::endl
      << "BR with no excluded SM fermions = "
      << testSbsxFullCascade.getBrToEwino( &testSmParticleExclusion );
      testSmParticleExclusion.push_back( CppSLHA::PDG_code::down );
      std::cout
      << std::endl
      << "BR excluding downs (but not antidowns) in final state = "
      << testSbsxFullCascade.getBrToEwino( &testSmParticleExclusion );
      std::cout << std::endl;
      std::cout
      << std::endl
      << ".getAcceptance( [particle], &testAcceptanceCutSet, 0, 1, 1, 0, 0 )"
      << " = " << testSbsxFullCascade.getAcceptance( true,
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
      << testSbsxFullCascade.getAcceptance( false,
                                            &testAcceptanceCutSet,
                                            0,
                                            1,
                                            1,
                                            0,
                                            0 );
      std::cout
      << std::endl
      << "testSbsxFullCascade seems OK (~10% for W^- to pass e^-, ~3% for"
      << " chargino as shown by testSxFullCascade).";
      std::cout << std::endl;
      std::cout << std::endl;


      gjsbsxFullCascade testGjsbsxFullCascade;
      testGjsbsxFullCascade.setProperties( &testSbsxFullCascade );
      std::cout
      << std::endl
      << "testGjsbsxFullCascade set with testSbsxFullCascade";
      std::cout << std::endl;
      testSmParticleExclusion.clear();
      std::cout
      << std::endl
      << "BR with no excluded SM fermions = "
      << testGjsbsxFullCascade.getBrToEwino( &testSmParticleExclusion );
      testSmParticleExclusion.push_back( CppSLHA::PDG_code::down );
      std::cout
      << std::endl
      << "BR excluding downs (but not antidowns) in final state = "
      << testGjsbsxFullCascade.getBrToEwino( &testSmParticleExclusion );
      std::cout << std::endl;
      std::cout
      << std::endl
      << ".getAcceptance( [particle], &testAcceptanceCutSet, 0, 1, 1, 0, 0 )"
      << " = " << testGjsbsxFullCascade.getAcceptance( true,
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
      << testGjsbsxFullCascade.getAcceptance( false,
                                              &testAcceptanceCutSet,
                                              0,
                                              1,
                                              1,
                                              0,
                                              0 );
      std::cout
      << std::endl
      << "testGjsbsxFullCascade seems OK.";
      std::cout << std::endl;
      std::cout << std::endl;


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
      << "has " << openCascades->size() << " open cascades. 1 is sx, 3 is"
      << " sjgx. BR is for ->getBrToEwino( &testSmParticleExclusion ), which"
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
    << "has " << openCascades->size() << " open cascades. 2 is gx, 4 is gjsx,"
    << " 8 is gjsbsx."
    << " BR is for ->getBrToEwino( &testSmParticleExclusion ), which"
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
    << "has " << openCascades->size() << " open cascades. 1 is sx, 3 is"
    << " sjgx. BR is for ->getBrToEwino( &testSmParticleExclusion ), which"
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
    << "has " << openCascades->size() << " open cascades. 2 is gx, 4 is gjsx,"
    << " 8 is gjsbsx. BR is for ->getBrToEwino( &testSmParticleExclusion ),"
    << " which excludes top & antitop, \"acc\" is for ->getAcceptance("
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
