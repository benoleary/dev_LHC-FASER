/*
 * lhcFaser_testing.hpp
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

#ifndef LHCFASER_TESTING_HPP_
#define LHCFASER_TESTING_HPP_

#include "basic_testing_stuff.hpp"
#include "input_testing.hpp"
#include "../LHC-FASER.hpp"

namespace LHC_FASER
{
  class lhcFaserTesting
  {
  public:
    lhcFaserTesting( basicStuffTesting* basicStuff,
                     inputTesting* inputTester ) :
        basicStuff( basicStuff ),
        inputTester( inputTester ),
        testInputHandler( inputTester->getInput() )
    {
      // just an initialization list.
    }
    ~lhcFaserTesting()
    {
      // does nothing.
    }

    void
    performTest()
    {
      lhcFaser* testLhcFaserPointer( new lhcFaser( "SPS1a_spectrum.out" ) );
      testLhcFaserPointer->addSignal( "sigmaBreakdownTest" );
      testLhcFaserPointer->addSignal( "Atlas4jMET_0l_7TeV" );
      testLhcFaserPointer->addSignal( "badSignalName" );
      std::cout
      << std::endl
      << "added signals without obvious problem.";
      std::cout << std::endl;
      delete testLhcFaserPointer;
      std::cout
      << std::endl
      << "no obvious problems in lhcFaser constructor & destructor.";
      std::cout << std::endl;

      lhcFaser sps1aLhcFaser( "SPS1a_spectrum.out" );
      lhcFaser sps2LhcFaser( "SPS2_spectrum.out" );
      sps1aLhcFaser.addSignal( "sigmaBreakdownTest" );
      signalHandler*
      sps1aSigma( sps1aLhcFaser.getSignal( "sigmaBreakdownTest" ) );
      signalHandler*
      sps2Sigma( sps2LhcFaser.addSignal( "sigmaBreakdownTest" ) );
      signalHandler*
      sps1aAtlas4j0l( sps1aLhcFaser.addSignal( "Atlas4jMET_0l_7TeV" ) );
      sps2LhcFaser.addSignal( "Atlas4jMET_0l_7TeV" );
      signalHandler*
      sps2Atlas4j0l( sps2LhcFaser.getSignal( "Atlas4jMET_0l_7TeV" ) );
      sps1aLhcFaser.addSignal( "Atlas3jMET_1l_7TeV" );
      signalHandler*
      sps1aAtlas3j1l( sps1aLhcFaser.getSignal( "Atlas3jMET_1l_7TeV" ) );
      signalHandler*
      sps2Atlas3j1l( sps2LhcFaser.addSignal( "Atlas3jMET_1l_7TeV" ) );
      signalHandler*
      sps1aBad( sps1aLhcFaser.addSignal( "badSignalName" ) );
      sps2LhcFaser.addSignal( "badSignalName" );
      signalHandler*
      sps2Bad( sps2LhcFaser.getSignal( "badSignalName" ) );

      std::cout << std::endl;
      std::cout
      << std::endl
      << "all 4 test signals from test lhcFasers have pointers.";
      std::cout << std::endl;

      sps1aLhcFaser.updateForNewSlha();
      std::cout << std::endl;
      std::cout
      << std::endl
      << "sps1aLhcFaser.updateForNewSlha() successful.";
      std::cout << std::endl;
      sps2LhcFaser.updateForNewSlha();
      std::cout << std::endl;
      std::cout
      << std::endl
      << "sps2LhcFaser.updateForNewSlha() successful.";
      std::cout << std::endl;
      std::cout
      << std::endl << "after updating for new point, ";
      std::cout << std::endl
      << *(sps1aSigma->getName()) << " (sps1aSigma) has signal value "
      << sps1aSigma->getValue();
      std::cout << std::endl;
      std::cout
      << *(sps1aAtlas4j0l->getName()) << " (sps1aAtlas4j0l) has signal value "
      << sps1aAtlas4j0l->getValue();
      std::cout << std::endl;
      std::cout
      << *(sps1aAtlas3j1l->getName());
      std::cout << " (sps1aAtlas3j1l) has signal value ";
      std::cout << sps1aAtlas3j1l->getValue();
      std::cout << std::endl;
      std::cout
      << *(sps1aBad->getName()) << " (sps1aBad) has signal value "
      << sps1aBad->getValue();
      std::cout << std::endl
      << *(sps2Sigma->getName()) << " (sps2Sigma) has signal value "
      << sps2Sigma->getValue();
      std::cout << std::endl;
      std::cout
      << *(sps2Atlas4j0l->getName()) << " (sps2Atlas4j0l) has signal value "
      << sps2Atlas4j0l->getValue();
      std::cout << std::endl;
      std::cout
      << *(sps2Atlas3j1l->getName()) << " (sps2Atlas3j1l) has signal value "
      << sps2Atlas3j1l->getValue();
      std::cout << std::endl;
      std::cout
      << *(sps2Bad->getName()) << " (sps2Bad) has signal value "
      << sps2Bad->getValue();
      std::cout << std::endl;
    }


  protected:
    basicStuffTesting* basicStuff;
    inputTesting* inputTester;
    inputHandler const* const testInputHandler;
  };

}  // end of LHC_FASER namespace

#endif /* LHCFASER_TESTING_HPP_ */
