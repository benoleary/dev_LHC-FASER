/*
 * input_testing.hpp
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

#ifndef INPUT_TESTING_HPP_
#define INPUT_TESTING_HPP_

#include "basic_testing_stuff.hpp"
#include "../LHC-FASER_input_handling_stuff.hpp"

namespace LHC_FASER
{
  class inputTesting
  {
  public:
    inputTesting( basicStuffTesting* basicStuff ) :
      slhaData( basicStuff->getSlha() ),
      testReadier( basicStuff->getReadier() ),
      testInputHandler( slhaData,
                        slhaData->inspect_particle_spectrum(),
                        "./grids",
                        basicStuff->getReadier() )
    {
      testInputHandler.setVerbosity( true );
    }
    ~inputTesting()
    {
      // does nothing.
    }

    inputHandler const* const
    getInput()
    {
      return &testInputHandler;
    }

    void
    performTest()
    {
      testInputHandler.setVerbosity( true );
      std::cout
      << std::endl
      << "is verbose? " << testInputHandler.isVerbose()
      << std::endl
      << "set to be not verbose.";
      testInputHandler.setVerbosity( false );
      std::cout
      << std::endl
      << "is verbose? " << testInputHandler.isVerbose()
      << std::endl
      << "set to be verbose.";
      testInputHandler.setVerbosity( true );
      std::cout
      << std::endl
      << "is verbose? " << testInputHandler.isVerbose()
      << std::endl
      << "using NLO? " << testInputHandler.usingNlo()
      << std::endl
      << "set to not use NLO.";
      testInputHandler.setUsingNloFlag( false );
      std::cout
      << std::endl
      << "using NLO? " << testInputHandler.usingNlo()
      << std::endl
      << "set to use NLO.";
      testInputHandler.setUsingNloFlag( true );
      std::cout
      << std::endl
      << "using NLO? " << testInputHandler.usingNlo()
      << std::endl
      << std::endl
      << "pathToGrids = "
      << *(testInputHandler.getPathToGrids())
      << std::endl
      << "pathToCrossSectionGrids = "
      << *(testInputHandler.getPathToCrossSectionGrids())
      << std::endl
      << "(inspect: pathToCrossSectionGrids = "
      << *(testInputHandler.inspectPathToCrossSectionGrids()) << ")"
      << std::endl
      << "pathToKinematicsGrids = "
      << *(testInputHandler.getPathToKinematicsGrids())
      << std::endl
      << "(inspect: pathToKinematicsGrids = "
      << *(testInputHandler.inspectPathToKinematicsGrids()) << ")"
      << std::endl
      << std::endl
      << "g_EM = "
      << testInputHandler.getElectromagneticFineStructureConstant()
      << std::endl
      << std::endl
      << "SPS1a MASS:"
      << testInputHandler.getCppSlha()->get_MASS_as_string()
      << std::endl
      << std::endl
      << *(testInputHandler.getSpectrum()->get_PDG_coded_name(
                                                        CppSLHA::PDG_code::Z ))
      << " "
      << testInputHandler.getSpectrum()->get_PDG_coded_mass(
                                                         CppSLHA::PDG_code::Z )
      << " GeV"
      << std::endl
      << *(testInputHandler.getWPlus()->get_name()) << " "
      << testInputHandler.getWPlus()->get_mass() << " GeV"
      << std::endl
      << *(testInputHandler.getTop()->get_name()) << " "
      << testInputHandler.getTop()->get_mass() << " GeV"
      << std::endl
      << *(testInputHandler.getGluino()->get_name()) << " "
      << testInputHandler.getGluino()->get_mass() << " GeV"
      << std::endl
      << "averageSquarks4Mass = "
      << testInputHandler.getAverageSquarks4Mass()
      << std::endl
      << "scoloreds: " << std::endl;
      for( std::vector< particlePointer >::const_iterator
           setIterator( testInputHandler.getColoredSparticles()->begin() );
           testInputHandler.getColoredSparticles()->end() > setIterator;
           ++setIterator )
        {
          std::cout << *((*setIterator)->get_name()) << std::endl;
        }
      std::cout
      << std::endl
      << "charginos: " << std::endl;
      for( std::vector< particlePointer >::const_iterator
           setIterator( testInputHandler.getCharginos()->begin() );
           testInputHandler.getCharginos()->end() > setIterator;
           ++setIterator )
        {
          std::cout << *((*setIterator)->get_name()) << std::endl;
        }
      std::cout
      << std::endl
      << "neutralinos: " << std::endl;
      for( std::vector< particlePointer >::const_iterator
           setIterator( testInputHandler.getNeutralinos()->begin() );
           testInputHandler.getNeutralinos()->end() > setIterator;
           ++setIterator )
        {
          std::cout << *((*setIterator)->get_name()) << std::endl;
        }
      std::cout
      << std::endl
      << "unstable neutralinos: " << std::endl;
      for( std::vector< particlePointer >::const_iterator
           setIterator( testInputHandler.getUnstableNeutralinos()->begin() );
           testInputHandler.getUnstableNeutralinos()->end() > setIterator;
           ++setIterator )
        {
          std::cout << *((*setIterator)->get_name()) << std::endl;
        }
      std::cout
      << std::endl
      << "electroweakinos: " << std::endl;
      for( std::vector< particlePointer >::const_iterator
           setIterator( testInputHandler.getElectroweakinos()->begin() );
           testInputHandler.getElectroweakinos()->end() > setIterator;
           ++setIterator )
        {
          std::cout << *((*setIterator)->get_name()) << std::endl;
        }
      std::cout
      << std::endl
      << "unstable electroweakinos: " << std::endl;
      for( std::vector< particlePointer >::const_iterator
           setIterator(
                      testInputHandler.getUnstableElectroweakinos()->begin() );
           testInputHandler.getUnstableElectroweakinos()->end() > setIterator;
           ++setIterator )
        {
          std::cout << *((*setIterator)->get_name()) << std::endl;
        }
      std::cout
      << std::endl
      << "charged sleptons: " << std::endl;
      for( std::vector< particlePointer >::const_iterator
           setIterator( testInputHandler.getChargedSleptons()->begin() );
           testInputHandler.getChargedSleptons()->end() > setIterator;
           ++setIterator )
        {
          std::cout << *((*setIterator)->get_name()) << std::endl;
        }
      std::cout
      << std::endl
      << "sneutrinos: " << std::endl;
      for( std::vector< particlePointer >::const_iterator
           setIterator( testInputHandler.getSneutrinos()->begin() );
           testInputHandler.getSneutrinos()->end() > setIterator;
           ++setIterator )
        {
          std::cout << *((*setIterator)->get_name()) << std::endl;
        }
      std::cout
      << std::endl
      << "sdowns: " << std::endl;
      for( std::vector< particlePointer >::const_iterator
           setIterator( testInputHandler.getSdownTypes()->begin() );
           testInputHandler.getSdownTypes()->end() > setIterator;
           ++setIterator )
        {
          std::cout << *((*setIterator)->get_name()) << std::endl;
        }
      std::cout
      << std::endl
      << "sups: " << std::endl;
      for( std::vector< particlePointer >::const_iterator
           setIterator( testInputHandler.getSupTypes()->begin() );
           testInputHandler.getSupTypes()->end() > setIterator;
           ++setIterator )
        {
          std::cout << *((*setIterator)->get_name()) << std::endl;
        }
      std::cout
      << std::endl
      << "jets: " << std::endl;
      for( std::list< int >::const_iterator
           setIterator( testInputHandler.getJets()->begin() );
           testInputHandler.getJets()->end() != setIterator;
           ++setIterator )
        {
          std::cout
          << *(testInputHandler.getSpectrum()->get_PDG_coded_name(
                                                              (*setIterator) ))
          << std::endl;
        }
      std::cout
      << std::endl
      << "jets4: " << std::endl;
      for( std::list< int >::const_iterator
           setIterator( testInputHandler.getJets4()->begin() );
           testInputHandler.getJets4()->end() != setIterator;
           ++setIterator )
        {
          std::cout
          << *(testInputHandler.getSpectrum()->get_PDG_coded_name(
                                                              (*setIterator) ))
          << std::endl;
        }
      std::cout
      << std::endl
      << "jets5: " << std::endl;
      for( std::list< int >::const_iterator
           setIterator( testInputHandler.getJets5()->begin() );
           testInputHandler.getJets5()->end() != setIterator;
           ++setIterator )
        {
          std::cout
          << *(testInputHandler.getSpectrum()->get_PDG_coded_name(
                                                              (*setIterator) ))
          << std::endl;
        }
      std::cout
      << std::endl
      << "not in jets4: " << std::endl;
      for( std::list< int >::const_iterator
           setIterator( testInputHandler.getNotInJets4()->begin() );
           testInputHandler.getNotInJets4()->end() != setIterator;
           ++setIterator )
        {
          std::cout
          << *(testInputHandler.getSpectrum()->get_PDG_coded_name(
                                                              (*setIterator) ))
          << std::endl;
        }
      std::cout
      << std::endl
      << "not in jets5: " << std::endl;
      for( std::list< int >::const_iterator
           setIterator( testInputHandler.getNotInJets5()->begin() );
           testInputHandler.getNotInJets5()->end() != setIterator;
           ++setIterator )
        {
          std::cout
          << *(testInputHandler.getSpectrum()->get_PDG_coded_name(
                                                              (*setIterator) ))
          << std::endl;
        }
      std::cout
      << std::endl
      << "squarks: " << std::endl;
      for( std::vector< particlePointer >::const_iterator
           setIterator( testInputHandler.getSquarks()->begin() );
           testInputHandler.getSquarks()->end() > setIterator;
           ++setIterator )
        {
          std::cout << *((*setIterator)->get_name()) << std::endl;
        }
      std::cout
      << std::endl
      << "squarks4: " << std::endl;
      for( std::vector< particlePointer >::const_iterator
           setIterator( testInputHandler.getSquarks4()->begin() );
           testInputHandler.getSquarks4()->end() > setIterator;
           ++setIterator )
        {
          std::cout << *((*setIterator)->get_name()) << std::endl;
        }
      std::cout
      << std::endl
      << "squarks5: " << std::endl;
      for( std::vector< particlePointer >::const_iterator
           setIterator( testInputHandler.getSquarks5()->begin() );
           testInputHandler.getSquarks5()->end() > setIterator;
           ++setIterator )
        {
          std::cout << *((*setIterator)->get_name()) << std::endl;
        }
      std::cout
      << std::endl
      << "sbottom_1 is in squarks? "
      << testInputHandler.isIn(
                              testInputHandler.getSbottomOne()->get_PDG_code(),
                                testInputHandler.getSquarks() )
      << std::endl
      << "sbottom_1 is in squarks4? "
      << testInputHandler.isIn(
                              testInputHandler.getSbottomOne()->get_PDG_code(),
                                testInputHandler.getSquarks4() )
      << std::endl;
      std::cout
      << std::endl
      << std::endl
      << "sparticle production combinations: " << std::endl;
      for( std::vector< signedParticleShortcutPair* >::const_iterator
           setIterator( testInputHandler.getScoloredProductionCombinations(
                                                                  )->begin() );
           testInputHandler.getScoloredProductionCombinations()->end()
           > setIterator;
           ++setIterator )
        {
          if( (*setIterator)->firstIsNotAntiparticle() )
            {
              std::cout << *((*setIterator)->getFirstParticle()->get_name());
            }
          else
            {
              std::cout
              << *((*setIterator)->getFirstParticle()->get_antiname());
            }
          std::cout << " + ";
          if( (*setIterator)->secondIsNotAntiparticle() )
            {
              std::cout << *((*setIterator)->getSecondParticle()->get_name());
            }
          else
            {
              std::cout
              << *((*setIterator)->getSecondParticle()->get_antiname());
            }
          std::cout << std::endl;
        }

      effectiveSquarkMassHolder*
      sdownLEffectiveMass( testInputHandler.getOnShellEffectiveSquarkMass(
                                              testInputHandler.getSdownL() ) );
      effectiveSquarkMassHolder*
      sbottomOneEffectiveMass( testInputHandler.getOnShellEffectiveSquarkMass(
                                          testInputHandler.getSbottomOne() ) );
      effectiveSquarkMassHolder*
      sdownLEffectiveMass2( testInputHandler.getOnShellEffectiveSquarkMass(
                                              testInputHandler.getSdownL() ) );
      std::cout
      << std::endl
      << "sdown_L effective mass holder = "
      << sdownLEffectiveMass
      << ", mass = " << sdownLEffectiveMass->getEffectiveSquarkMass()
      << std::endl
      << "2nd sdown_L effective mass holder = "
      << sdownLEffectiveMass2
      << std::endl
      << "sbottom_1 effective mass holder = "
      << sbottomOneEffectiveMass
      << ", mass = " << sbottomOneEffectiveMass->getEffectiveSquarkMass()
      << std::endl;


      std::cout
      << std::endl
      << "sdown_L to neutralino_1 fermion handedness = "
      << testInputHandler.quarkOrLeptonLeftHandedness(
                                                    CppSLHA::PDG_code::sdown_L,
                                            CppSLHA::PDG_code::neutralino_one )
      << std::endl
      << "sdown_R to neutralino_1 fermion handedness = "
      << testInputHandler.quarkOrLeptonLeftHandedness(
                                                    CppSLHA::PDG_code::sdown_R,
                                            CppSLHA::PDG_code::neutralino_one )
      << std::endl
      << "sbottom_1 to neutralino_1 fermion handedness = "
      << testInputHandler.quarkOrLeptonLeftHandedness(
                                                CppSLHA::PDG_code::sbottom_one,
                                            CppSLHA::PDG_code::neutralino_one )
      << std::endl
      << "sbottom_2 to neutralino_1 fermion handedness = "
      << testInputHandler.quarkOrLeptonLeftHandedness(
                                                CppSLHA::PDG_code::sbottom_two,
                                            CppSLHA::PDG_code::neutralino_one )
      << std::endl;


      exclusiveBrCalculator*
      gluinoToSdownLBr( testInputHandler.getExclusiveBrCalculator(
                                                  testInputHandler.getGluino(),
                                                  testInputHandler.getSdownL(),
                                                                   true,
                                          testInputHandler.getNotInJets5() ) );
      exclusiveBrCalculator*
      gluinoToSdownRBr( testInputHandler.getExclusiveBrCalculator(
                                                  testInputHandler.getGluino(),
                                                  testInputHandler.getSdownR(),
                                                                   true,
                                          testInputHandler.getNotInJets5() ) );
      exclusiveBrCalculator*
      sdownLToGluinoBr( testInputHandler.getExclusiveBrCalculator(
                                                  testInputHandler.getSdownL(),
                                                  testInputHandler.getGluino(),
                                                                   true,
                                          testInputHandler.getNotInJets5() ) );
      exclusiveBrCalculator*
      sdownRToGluinoBr( testInputHandler.getExclusiveBrCalculator(
                                                  testInputHandler.getSdownR(),
                                                  testInputHandler.getGluino(),
                                                                   true,
                                          testInputHandler.getNotInJets5() ) );
      exclusiveBrCalculator*
      gluinoToSdownLBr2( testInputHandler.getExclusiveBrCalculator(
                                                  testInputHandler.getGluino(),
                                                  testInputHandler.getSdownL(),
                                                                    true,
                                          testInputHandler.getNotInJets5() ) );
      std::cout
      << std::endl
      << "exclusive BR checker: gluino -> sdown_L"
      << std::endl << "pointer = "
      << gluinoToSdownLBr
      << std::endl << *(gluinoToSdownLBr->getDecayer()->get_name())
      << " -> " << *(gluinoToSdownLBr->getProduct()->get_name())
      << " has BR " << gluinoToSdownLBr->getBr()
      << std::endl
      << "2nd gluino -> sdown_L"
      << std::endl << "pointer = "
      << gluinoToSdownLBr2
      << std::endl << *(gluinoToSdownLBr2->getDecayer()->get_name())
      << " -> " << *(gluinoToSdownLBr2->getProduct()->get_name())
      << " has BR " << gluinoToSdownLBr2->getBr()
      << std::endl
      << "exclusive BR checker: sdown_L -> gluino"
      << std::endl << "pointer = "
      << sdownLToGluinoBr
      << std::endl << *(sdownLToGluinoBr->getDecayer()->get_name())
      << " -> " << *(sdownLToGluinoBr->getProduct()->get_name())
      << " has BR " << sdownLToGluinoBr->getBr()
      << std::endl
      << "gluino -> sdown_R"
      << std::endl << "pointer = "
      << gluinoToSdownRBr
      << std::endl << *(gluinoToSdownRBr->getDecayer()->get_name())
      << " -> " << *(gluinoToSdownRBr->getProduct()->get_name())
      << " has BR " << gluinoToSdownRBr->getBr()
      << std::endl
      << "exclusive BR checker: sdown_R -> gluino"
      << std::endl << "pointer = "
      << sdownLToGluinoBr
      << std::endl << *(sdownRToGluinoBr->getDecayer()->get_name())
      << " -> " << *(sdownRToGluinoBr->getProduct()->get_name())
      << " has BR " << sdownRToGluinoBr->getBr()
      << std::endl
      << std::endl;


      slhaData->read_file( "SPS2_spectrum.out" );
      testReadier->readyObserversForNewPoint();
      std::cout
      << std::endl
      << "SPS2 MASS:"
      << testInputHandler.getCppSlha()->get_MASS_as_string()
      << std::endl
      << std::endl
      << *(testInputHandler.getGluino()->get_name()) << " "
      << testInputHandler.getGluino()->get_mass() << " GeV"
      << std::endl
      << std::endl
      << "charged sleptons: " << std::endl;
      for( std::vector< particlePointer >::const_iterator
           setIterator( testInputHandler.getChargedSleptons()->begin() );
           testInputHandler.getChargedSleptons()->end() > setIterator;
           ++setIterator )
        {
          std::cout << *((*setIterator)->get_name()) << std::endl;
        }
      std::cout
      << std::endl
      << "exclusive BR checker: gluino -> sdown_L"
      << std::endl << "pointer = "
      << gluinoToSdownLBr
      << std::endl << *(gluinoToSdownLBr->getDecayer()->get_name())
      << " -> " << *(gluinoToSdownLBr->getProduct()->get_name())
      << " has BR " << gluinoToSdownLBr->getBr()
      << std::endl
      << "2nd gluino -> sdown_L"
      << std::endl << "pointer = "
      << gluinoToSdownLBr2
      << std::endl << *(gluinoToSdownLBr2->getDecayer()->get_name())
      << " -> " << *(gluinoToSdownLBr2->getProduct()->get_name())
      << " has BR " << gluinoToSdownLBr2->getBr()
      << std::endl
      << "exclusive BR checker: sdown_L -> gluino"
      << std::endl << "pointer = "
      << sdownLToGluinoBr
      << std::endl << *(sdownLToGluinoBr->getDecayer()->get_name())
      << " -> " << *(sdownLToGluinoBr->getProduct()->get_name())
      << " has BR " << sdownLToGluinoBr->getBr()
      << std::endl
      << "gluino -> sdown_R"
      << std::endl << "pointer = "
      << gluinoToSdownRBr
      << std::endl << *(gluinoToSdownRBr->getDecayer()->get_name())
      << " -> " << *(gluinoToSdownRBr->getProduct()->get_name())
      << " has BR " << gluinoToSdownRBr->getBr()
      << std::endl
      << "exclusive BR checker: sdown_R -> gluino"
      << std::endl << "pointer = "
      << sdownLToGluinoBr
      << std::endl << *(sdownRToGluinoBr->getDecayer()->get_name())
      << " -> " << *(sdownRToGluinoBr->getProduct()->get_name())
      << " has BR " << sdownRToGluinoBr->getBr()
      << std::endl
      << std::endl;


      slhaData->read_file( "SPS1a_spectrum.out" );
      testReadier->readyObserversForNewPoint();
      std::cout
      << std::endl
      << "SPS1a MASS:"
      << testInputHandler.getCppSlha()->get_MASS_as_string()
      << std::endl
      << std::endl
      << *(testInputHandler.getGluino()->get_name()) << " "
      << testInputHandler.getGluino()->get_mass() << " GeV"
      << std::endl
      << std::endl
      << "charged sleptons: " << std::endl;
      for( std::vector< particlePointer >::const_iterator
           setIterator( testInputHandler.getChargedSleptons()->begin() );
           testInputHandler.getChargedSleptons()->end() > setIterator;
           ++setIterator )
        {
          std::cout << *((*setIterator)->get_name()) << std::endl;
        }
      std::cout
      << std::endl
      << "exclusive BR checker: gluino -> sdown_L"
      << std::endl << "pointer = "
      << gluinoToSdownLBr
      << std::endl << *(gluinoToSdownLBr->getDecayer()->get_name())
      << " -> " << *(gluinoToSdownLBr->getProduct()->get_name())
      << " has BR " << gluinoToSdownLBr->getBr()
      << std::endl
      << "2nd gluino -> sdown_L"
      << std::endl << "pointer = "
      << gluinoToSdownLBr2
      << std::endl << *(gluinoToSdownLBr2->getDecayer()->get_name())
      << " -> " << *(gluinoToSdownLBr2->getProduct()->get_name())
      << " has BR " << gluinoToSdownLBr2->getBr()
      << std::endl
      << "exclusive BR checker: sdown_L -> gluino"
      << std::endl << "pointer = "
      << sdownLToGluinoBr
      << std::endl << *(sdownLToGluinoBr->getDecayer()->get_name())
      << " -> " << *(sdownLToGluinoBr->getProduct()->get_name())
      << " has BR " << sdownLToGluinoBr->getBr()
      << std::endl
      << "gluino -> sdown_R"
      << std::endl << "pointer = "
      << gluinoToSdownRBr
      << std::endl << *(gluinoToSdownRBr->getDecayer()->get_name())
      << " -> " << *(gluinoToSdownRBr->getProduct()->get_name())
      << " has BR " << gluinoToSdownRBr->getBr()
      << std::endl
      << "exclusive BR checker: sdown_R -> gluino"
      << std::endl << "pointer = "
      << sdownLToGluinoBr
      << std::endl << *(sdownRToGluinoBr->getDecayer()->get_name())
      << " -> " << *(sdownRToGluinoBr->getProduct()->get_name())
      << " has BR " << sdownRToGluinoBr->getBr()
      << std::endl
      << std::endl;

      slhaData->read_file( "SPS2_spectrum.out" );
      testReadier->readyObserversForNewPoint();
      std::cout
      << std::endl
      << "exclusive BR checker: gluino -> sdown_L"
      << std::endl << "pointer = "
      << gluinoToSdownLBr
      << std::endl << *(gluinoToSdownLBr->getDecayer()->get_name())
      << " -> " << *(gluinoToSdownLBr->getProduct()->get_name())
      << " has BR " << gluinoToSdownLBr->getBr()
      << std::endl
      << "exclusive BR checker: sdown_L -> gluino"
      << std::endl << "pointer = "
      << sdownLToGluinoBr
      << std::endl << *(sdownLToGluinoBr->getDecayer()->get_name())
      << " -> " << *(sdownLToGluinoBr->getProduct()->get_name())
      << " has BR " << sdownLToGluinoBr->getBr()
      << std::endl
      << "gluino -> sdown_R"
      << std::endl << "pointer = "
      << gluinoToSdownRBr
      << std::endl << *(gluinoToSdownRBr->getDecayer()->get_name())
      << " -> " << *(gluinoToSdownRBr->getProduct()->get_name())
      << " has BR " << gluinoToSdownRBr->getBr()
      << std::endl
      << "exclusive BR checker: sdown_R -> gluino"
      << std::endl << "pointer = "
      << sdownLToGluinoBr
      << std::endl << *(sdownRToGluinoBr->getDecayer()->get_name())
      << " -> " << *(sdownRToGluinoBr->getProduct()->get_name())
      << " has BR " << sdownRToGluinoBr->getBr()
      << std::endl
      << std::endl;
      slhaData->read_file( "SPS1a_spectrum.out" );
      testReadier->readyObserversForNewPoint();
      std::cout
      << std::endl
      << "exclusive BR checker: gluino -> sdown_L"
      << std::endl << "pointer = "
      << gluinoToSdownLBr
      << std::endl << *(gluinoToSdownLBr->getDecayer()->get_name())
      << " -> " << *(gluinoToSdownLBr->getProduct()->get_name())
      << " has BR " << gluinoToSdownLBr->getBr()
      << std::endl
      << "exclusive BR checker: sdown_L -> gluino"
      << std::endl << "pointer = "
      << sdownLToGluinoBr
      << std::endl << *(sdownLToGluinoBr->getDecayer()->get_name())
      << " -> " << *(sdownLToGluinoBr->getProduct()->get_name())
      << " has BR " << sdownLToGluinoBr->getBr()
      << std::endl
      << "gluino -> sdown_R"
      << std::endl << "pointer = "
      << gluinoToSdownRBr
      << std::endl << *(gluinoToSdownRBr->getDecayer()->get_name())
      << " -> " << *(gluinoToSdownRBr->getProduct()->get_name())
      << " has BR " << gluinoToSdownRBr->getBr()
      << std::endl
      << "exclusive BR checker: sdown_R -> gluino"
      << std::endl << "pointer = "
      << sdownLToGluinoBr
      << std::endl << *(sdownRToGluinoBr->getDecayer()->get_name())
      << " -> " << *(sdownRToGluinoBr->getProduct()->get_name())
      << " has BR " << sdownRToGluinoBr->getBr()
      << std::endl
      << std::endl
      << "gluino pointer = " << testInputHandler.getGluino()
      << std::endl
      << "sdown_L pointer = " << testInputHandler.getSdownL()
      << std::endl
      << "sdown_R pointer = " << testInputHandler.getSdownR()
      << std::endl;

      LHC_FASER::readierForNewPoint*
      testReadier2( new LHC_FASER::readierForNewPoint() );
      exclusiveBrCalculator*
      sdownLCharginoOneMinus( new exclusiveBrCalculator(
                                                  testInputHandler.getSdownL(),
                                             testInputHandler.getCharginoOne(),
                                                         false,
                                              testInputHandler.getNotInJets4(),
                                                         testReadier2 ) );
      exclusiveBrCalculator*
      sdownLCharginoOneMinus2( new exclusiveBrCalculator(
                                                  testInputHandler.getSdownL(),
                                             testInputHandler.getCharginoOne(),
                                                          false,
                                              testInputHandler.getNotInJets4(),
                                                          testReadier2 ) );
      exclusiveBrCalculator*
      sdownLCharginoOnePlus( new exclusiveBrCalculator(
                                                  testInputHandler.getSdownL(),
                                             testInputHandler.getCharginoOne(),
                                                        true,
                                              testInputHandler.getNotInJets4(),
                                                              testReadier2 ) );
      std::cout
      << std::endl << "sdownLCharginoOneMinus = "
      << sdownLCharginoOneMinus << ", BR = " << sdownLCharginoOneMinus->getBr()
      << std::endl << "sdownLCharginoOneMinus2 = " << sdownLCharginoOneMinus2
      << ", BR = " << sdownLCharginoOneMinus2->getBr() << std::endl
      << "sdownLCharginoOnePlus BR = " << sdownLCharginoOnePlus->getBr();
      testReadier2->readyObserversForNewPoint();
      delete sdownLCharginoOnePlus;
      std::cout
      << std::endl << "sdownLCharginoOneMinus BR = "
      << sdownLCharginoOneMinus->getBr();
      testReadier2->readyObserversForNewPoint();
      delete testReadier2;
      delete sdownLCharginoOneMinus2;
      std::cout
      << std::endl << "sdownLCharginoOneMinus BR = "
      << sdownLCharginoOneMinus->getBr();
      std::cout
      << std::endl;
      delete sdownLCharginoOneMinus;
      testInputHandler.setVerbosity( true );
    }


  protected:
    CppSLHA::CppSLHA2* slhaData;
    readierForNewPoint* testReadier;
    inputHandler testInputHandler;
  };

}  // end of LHC_FASER namespace

#endif /* INPUT_TESTING_HPP_ */
