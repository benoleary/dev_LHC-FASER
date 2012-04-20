/*
 * LHC-FASER_for_JL.cpp
 *
 *  Created on: 20 Apr 2012
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

/* this program is to print out several rates for a given input file.
 * ID: scan ID corresponding the SLHA file
 * Nevent: number of produced events
 * ACCj2: CMS 2-jet  -  huh?
 * ACCj3: Atlas 3-jet  -  huh?
 * ACCj4: Atlas 4-jet  -  huh?
 * ACCl0: 0 lepton  -  huh?
 * ACCl1: 1 lepton  -  huh?
 * ACCl2: 2 lepton  -  huh?
 * ACCos: 2 os leptons  -  huh?
 * ACCss: 2 ss leptons  -  huh?
 * ACCossfsubst 2 ossf - osdf leptons  -  huh?
 * ACC1: Atlas4jMET0l
 * ACC2: Atlas4jjMET1l
 * ...(order as in 3.3 in your note)
 * ACC12: CMSsamesigndilepton
 */


// includes:

#include "LHC-FASER.hpp"


// the main program:

int main( int argumentCount,
          char* argumentStrings[] )
{
  if( 3 != argumentCount )
    // if the input was definitely wrong...
  {
    std::cout
    << std::endl
    << "error! this program requires the name of the SLHA spectrum file as"
    << " its 1st argument and the name of the file to write as output as its"
    << " 2nd argument!";
    std::cout << std::endl;  // let the user know the format.
  }
  else  // otherwise, do what is supposed to be done.
  {
    std::string spectrumFilename( argumentStrings[ 1 ] );
    std::string outputFilename( argumentStrings[ 2 ] );

    LHC_FASER::lhcFaser comparisonLhcFaser( spectrumFilename,
                                            "./testing/grids/",
                                            "fb" );
    comparisonLhcFaser.setVerbosity( true );
    std::vector< LHC_FASER::signalHandler* > signalSet;

    // ACCj2
    signalSet.push_back( comparisonLhcFaser.addSignal( "ACCj2" ) );
    // not implemented, will always give NaN.

    // ACCj3
    signalSet.push_back( comparisonLhcFaser.addSignal(
                                             "Atlas3jMET_noExtraCut_07TeV" ) );
    // is this just the 3jet+MET cut, with no cut on leptons?

    // ACCj4
    signalSet.push_back( comparisonLhcFaser.addSignal(
                                             "Atlas4jMET_noExtraCut_07TeV" ) );
    // is this just the 4jet+MET cut, with no cut on leptons?

    // ACCl0
    signalSet.push_back( comparisonLhcFaser.addSignal(
                                                  "noJetOrMETCut_0l_07TeV" ) );
    // is this just the cut of 0 leptons passing a 10GeV pT cut, with no cut on
    // jet or MET?

    // ACCl1
    signalSet.push_back( comparisonLhcFaser.addSignal(
                                                  "noJetOrMETCut_1l_07TeV" ) );
    // is this just the cut of 1 lepton passing a 10GeV pT cut, with no cut on
    // jet or MET?

    // ACCl2
    signalSet.push_back( comparisonLhcFaser.addSignal(
                                                  "noJetOrMETCut_2l_07TeV" ) );
    // is this just the cut of 2 leptons passing a 10GeV pT cut, with no cut on
    // jet or MET?

    // ACCos
    signalSet.push_back( comparisonLhcFaser.addSignal(
                                              "noJetOrMETCut_1lm1lp_07TeV" ) );
    // is this just the cut of 1 lepton of each charge (regardless of flavor)
    // passing a 10GeV pT cut, with no cut on jet or MET?

    // ACCss
    signalSet.push_back( comparisonLhcFaser.addSignal(
                                                "noJetOrMETCut_sssf_07TeV" ) );
    // is this just the cut of 1 lepton of each charge (regardless of flavor)
    // passing a 10GeV pT cut, with no cut on jet or MET?
    std::cout
    << std::endl
    << "ACCss: noJetOrMETCut_sssf_07TeV does _not_ count same-sign dileptons"
    << " when they are not the same flavor.";
    std::cout << std::endl;

    // ACCossfsubst
    signalSet.push_back( comparisonLhcFaser.addSignal(
                                                "noJetOrMETCut_ossf_07TeV" ) );
    // is this just the cut of 1 OSSF-OSDF pair of leptons passing a 20GeV pT
    // cut (other failing a 10GeV pT cut), with no cut on jet or MET?

    // ACC1
    signalSet.push_back(
         comparisonLhcFaser.addSignal( "Atlas4jMET_0l_07TeV_pTl20GeV10GeV" ) );

    // ACC2
    signalSet.push_back(
         comparisonLhcFaser.addSignal( "Atlas4jMET_1l_07TeV_pTl20GeV10GeV" ) );

    // ACC3
    signalSet.push_back(
         comparisonLhcFaser.addSignal( "Atlas4jMET_2l_07TeV_pTl20GeV10GeV" ) );

    // ACC4
    signalSet.push_back(
         comparisonLhcFaser.addSignal( "Atlas3jMET_0l_07TeV_pTl20GeV10GeV" ) );

    // ACC5
    signalSet.push_back(
         comparisonLhcFaser.addSignal( "Atlas3jMET_1l_07TeV_pTl20GeV10GeV" ) );

    // ACC6
    signalSet.push_back(
         comparisonLhcFaser.addSignal( "Atlas3jMET_2l_07TeV_pTl20GeV10GeV" ) );

    // ACC7
    signalSet.push_back(
    comparisonLhcFaser.addSignal( "noJetOrMETCut_sssf_07TeV_pTl20GeV20GeV" ) );

    // ACC8
    signalSet.push_back(
       comparisonLhcFaser.addSignal( "Atlas4jMET_ossf_07TeV_pTl20GeV10GeV" ) );

    // ACC9
    signalSet.push_back(
       comparisonLhcFaser.addSignal( "Atlas3jMET_ossf_07TeV_pTl20GeV10GeV" ) );

    // ACC10
    signalSet.push_back(
              comparisonLhcFaser.addSignal( "CMS2jAlphaT_noExtraCut_07TeV" ) );

    // ACC11
    signalSet.push_back(
      comparisonLhcFaser.addSignal( "CMS2jAlphaT_ossf_07TeV_pTl20GeV10GeV" ) );

    // ACC12
    signalSet.push_back(
      comparisonLhcFaser.addSignal( "CMS2jAlphaT_sssf_07TeV_pTl20GeV20GeV" ) );

    comparisonLhcFaser.updateForNewSlha();

    std::ofstream outputStream( outputFilename.c_str() );
    if( outputStream.is_open() )
    {
      for( std::vector< LHC_FASER::signalHandler* >::iterator
           signalIterator( signalSet.begin() );
           signalSet.end() > signalIterator;
           ++signalIterator )
      {
        outputStream << " " << (*signalIterator)->getValue();
      }
      outputStream << std::endl;
      outputStream.close();
    }
    else
    {
      std::cout
      << std::endl
      << "error! could not open \"" << outputFilename
      << "\"! no output written!";
      std::cout << std::endl;
    }
  }

    // this was a triumph! I'm making a note here:
  return EXIT_SUCCESS;

}
