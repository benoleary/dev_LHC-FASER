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
 * it appends the ranges of rates for the following signals as an extra SLHA
 * block at the end of the input file.
 * acceptance ID:
 * 1 CMS2jalphaT_noExtraCut_07TeV / noJetOrMETCut_noExtraCut_07TeV
 * 2 Atlas3jMET_noExtraCut_07TeV / noJetOrMETCut_noExtraCut_07TeV
 * 3 Atlas4jMET_noExtraCut_07TeV / noJetOrMETCut_noExtraCut_07TeV
 * 4 noJetOrMETCut_0l_07TeV / noJetOrMETCut_noExtraCut_07TeV
 * 5 noJetOrMETCut_1l_07TeV / noJetOrMETCut_noExtraCut_07TeV
 * 6 noJetOrMETCut_2l_07TeV / noJetOrMETCut_noExtraCut_07TeV
 * 7 noJetOrMETCut_sssf_07TeV / noJetOrMETCut_noExtraCut_07TeV
 * 8 noJetOrMETCut_ossf_07TeV / noJetOrMETCut_noExtraCut_07TeV
 * or the 14 TeV version, replacing "07TeV" with "14TeV" in all the signal
 * names above.
 */


// includes:

#include "LHC-FASER.hpp"


// the main program:

int main( int argumentCount,
          char* argumentStrings[] )
{
  if( 2 != argumentCount )
    // if the input was definitely wrong...
  {
    std::cout
    << std::endl
    << "error! this program requires the name of the SLHA spectrum file as"
    << " its argument!";
    std::cout << std::endl;  // let the user know the format.
    return EXIT_FAILURE;
  }

  std::string slhaFilename( argumentStrings[ 1 ] );
  LHC_FASER::lhcFaser comparisonLhcFaser( slhaFilename,
                                          "./grids/",
                                          "fb",
                                          false );
  comparisonLhcFaser.setVerbosity( true );
  std::vector< LHC_FASER::signalHandler* > sevenTevSignals;
  sevenTevSignals.push_back( comparisonLhcFaser.addSignal(
                                          "noJetOrMETCut_noExtraCut_07TeV" ) );
  sevenTevSignals.push_back(
                comparisonLhcFaser.addSignal( "CMSalphaT_noExtraCut_07TeV" ) );
  sevenTevSignals.push_back(
               comparisonLhcFaser.addSignal( "Atlas3jMET_noExtraCut_07TeV" ) );
  sevenTevSignals.push_back(
               comparisonLhcFaser.addSignal( "Atlas4jMET_noExtraCut_07TeV" ) );
  sevenTevSignals.push_back(
                    comparisonLhcFaser.addSignal( "noJetOrMETCut_0l_07TeV" ) );
  sevenTevSignals.push_back(
                    comparisonLhcFaser.addSignal( "noJetOrMETCut_1l_07TeV" ) );
  sevenTevSignals.push_back(
                    comparisonLhcFaser.addSignal( "noJetOrMETCut_2l_07TeV" ) );
  sevenTevSignals.push_back(
                  comparisonLhcFaser.addSignal( "noJetOrMETCut_sssf_07TeV" ) );
  sevenTevSignals.push_back(
                  comparisonLhcFaser.addSignal( "noJetOrMETCut_ossf_07TeV" ) );

  std::vector< LHC_FASER::signalHandler* > fourteenTevSignals;
  fourteenTevSignals.push_back( comparisonLhcFaser.addSignal(
                                          "noJetOrMETCut_noExtraCut_14TeV" ) );
  fourteenTevSignals.push_back(
                comparisonLhcFaser.addSignal( "CMSalphaT_noExtraCut_14TeV" ) );
  fourteenTevSignals.push_back(
               comparisonLhcFaser.addSignal( "Atlas3jMET_noExtraCut_14TeV" ) );
  fourteenTevSignals.push_back(
               comparisonLhcFaser.addSignal( "Atlas4jMET_noExtraCut_14TeV" ) );
  fourteenTevSignals.push_back(
                    comparisonLhcFaser.addSignal( "noJetOrMETCut_0l_14TeV" ) );
  fourteenTevSignals.push_back(
                    comparisonLhcFaser.addSignal( "noJetOrMETCut_1l_14TeV" ) );
  fourteenTevSignals.push_back(
                    comparisonLhcFaser.addSignal( "noJetOrMETCut_2l_14TeV" ) );
  fourteenTevSignals.push_back(
                  comparisonLhcFaser.addSignal( "noJetOrMETCut_sssf_14TeV" ) );
  fourteenTevSignals.push_back(
                  comparisonLhcFaser.addSignal( "noJetOrMETCut_ossf_14TeV" ) );


  /*
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
                comparisonLhcFaser.addSignal( "CMSalphaT_noExtraCut_07TeV" ) );

    // ACC11
    signalSet.push_back(
        comparisonLhcFaser.addSignal( "CMSalphaT_ossf_07TeV_pTl20GeV10GeV" ) );

    // ACC12
    signalSet.push_back(
        comparisonLhcFaser.addSignal( "CMSalphaT_sssf_07TeV_pTl20GeV20GeV" ) );


    // extra pure lepton acceptances with harder cuts:
    // ACCl0pT50
    signalSet.push_back( comparisonLhcFaser.addSignal(
                                    "noJetOrMETCut_0l_07TeV_pTl50GeV50GeV" ) );
    // just the cut of 0 leptons passing a 50GeV pT cut, with no cut on jet or
    // MET

    // ACCl1pT50
    signalSet.push_back( comparisonLhcFaser.addSignal(
                                    "noJetOrMETCut_1l_07TeV_pTl50GeV50GeV" ) );
    // just the cut of 1 lepton passing a 50GeV pT cut, with no cut on jet or
    // MET

    // ACCl2pT50
    signalSet.push_back( comparisonLhcFaser.addSignal(
                                    "noJetOrMETCut_2l_07TeV_pTl50GeV50GeV" ) );
    // is this just the cut of 2 leptons passing a 10GeV pT cut, with no cut on
    // jet or MET?

    // ACCospT50
    signalSet.push_back( comparisonLhcFaser.addSignal(
                                "noJetOrMETCut_1lm1lp_07TeV_pTl50GeV50GeV" ) );
    // just the cut of 1 lepton of each charge (regardless of flavor) passing a
    // 50GeV pT cut, with no cut on jet or MET.

    // ACCsspT50
    signalSet.push_back( comparisonLhcFaser.addSignal(
                                  "noJetOrMETCut_sssf_07TeV_pTl50GeV50GeV" ) );
    // just the cut of 1 lepton of each charge (regardless of flavor) passing a
    // 50GeV pT cut, with no cut on jet or MET.


    // the whole set again, but for 14 TeV:
    signalSet.push_back( comparisonLhcFaser.addSignal( "ACCj2" ) );
    signalSet.push_back( comparisonLhcFaser.addSignal(
                                             "Atlas3jMET_noExtraCut_14TeV" ) );
    signalSet.push_back( comparisonLhcFaser.addSignal(
                                             "Atlas4jMET_noExtraCut_14TeV" ) );
    LHC_FASER::signalHandler*
    fourteenTevNoLepton( comparisonLhcFaser.addSignal(
                                                  "noJetOrMETCut_0l_14TeV" ) );
    signalSet.push_back( fourteenTevNoLepton );
    LHC_FASER::signalHandler*
    fourteenTevOneLepton( comparisonLhcFaser.addSignal(
                                                  "noJetOrMETCut_1l_14TeV" ) );
    signalSet.push_back( fourteenTevOneLepton );
    LHC_FASER::signalHandler*
    fourteenTevTwoLeptons( comparisonLhcFaser.addSignal(
                                                  "noJetOrMETCut_2l_14TeV" ) );
    signalSet.push_back( fourteenTevTwoLeptons );
    signalSet.push_back( comparisonLhcFaser.addSignal(
                                              "noJetOrMETCut_1lm1lp_14TeV" ) );
    signalSet.push_back( comparisonLhcFaser.addSignal(
                                                "noJetOrMETCut_sssf_14TeV" ) );
    signalSet.push_back( comparisonLhcFaser.addSignal(
                                                "noJetOrMETCut_ossf_14TeV" ) );
    signalSet.push_back(
         comparisonLhcFaser.addSignal( "Atlas4jMET_0l_14TeV_pTl20GeV10GeV" ) );
    signalSet.push_back(
         comparisonLhcFaser.addSignal( "Atlas4jMET_1l_14TeV_pTl20GeV10GeV" ) );
    signalSet.push_back(
         comparisonLhcFaser.addSignal( "Atlas4jMET_2l_14TeV_pTl20GeV10GeV" ) );
    signalSet.push_back(
         comparisonLhcFaser.addSignal( "Atlas3jMET_0l_14TeV_pTl20GeV10GeV" ) );
    signalSet.push_back(
         comparisonLhcFaser.addSignal( "Atlas3jMET_1l_14TeV_pTl20GeV10GeV" ) );
    signalSet.push_back(
         comparisonLhcFaser.addSignal( "Atlas3jMET_2l_14TeV_pTl20GeV10GeV" ) );
    signalSet.push_back(
    comparisonLhcFaser.addSignal( "noJetOrMETCut_sssf_14TeV_pTl20GeV20GeV" ) );
    signalSet.push_back(
       comparisonLhcFaser.addSignal( "Atlas4jMET_ossf_14TeV_pTl20GeV10GeV" ) );
    signalSet.push_back(
       comparisonLhcFaser.addSignal( "Atlas3jMET_ossf_14TeV_pTl20GeV10GeV" ) );
    signalSet.push_back(
                comparisonLhcFaser.addSignal( "CMSalphaT_noExtraCut_14TeV" ) );
    signalSet.push_back(
        comparisonLhcFaser.addSignal( "CMSalphaT_ossf_14TeV_pTl20GeV10GeV" ) );
    signalSet.push_back(
        comparisonLhcFaser.addSignal( "CMSalphaT_sssf_14TeV_pTl20GeV20GeV" ) );
    signalSet.push_back( comparisonLhcFaser.addSignal(
                                    "noJetOrMETCut_0l_14TeV_pTl50GeV50GeV" ) );
    signalSet.push_back( comparisonLhcFaser.addSignal(
                                    "noJetOrMETCut_1l_14TeV_pTl50GeV50GeV" ) );
    signalSet.push_back( comparisonLhcFaser.addSignal(
                                    "noJetOrMETCut_2l_14TeV_pTl50GeV50GeV" ) );
    signalSet.push_back( comparisonLhcFaser.addSignal(
                                "noJetOrMETCut_1lm1lp_14TeV_pTl50GeV50GeV" ) );
    signalSet.push_back( comparisonLhcFaser.addSignal(
                                  "noJetOrMETCut_sssf_14TeV_pTl50GeV50GeV" ) );

    LHC_FASER::signalHandler*
    noCutsForSevenTev( comparisonLhcFaser.addSignal(
                                          "noJetOrMETCut_noExtraCut_07TeV" ) );
    signalSet.push_back( noCutsForSevenTev );
    LHC_FASER::signalHandler*
    noCutsForFourteenTev( comparisonLhcFaser.addSignal(
                                          "noJetOrMETCut_noExtraCut_14TeV" ) );
    signalSet.push_back( noCutsForFourteenTev );
*/


  comparisonLhcFaser.updateForNewSlha();

  LHC_FASER::signalHandler* normalizationSignal( sevenTevSignals[ 0 ] );
  double
  normalizationFactor( 1.0 / normalizationSignal->getValue() );
  LHC_FASER::signalHandler* currentSignal( NULL );

  std::ofstream outputStream( slhaFilename.c_str(),
                              std::ios::app );
  if( !(outputStream.good()) )
  {
    std::cout
    << std::endl
    << "error! could not open \"" << slhaFilename
    << "\"! no output written!";
    std::cout << std::endl;

    return EXIT_FAILURE;
  }

  std::cout
  << std::endl
  << "noJetOrMETCut_noExtraCut_..."
  << std::endl
  << "07TeV: " << sevenTevSignals[ 0 ]->getValue()
  << std::endl
  << "14TeV: " << fourteenTevSignals[ 0 ]->getValue();
  std::cout << std::endl;


  outputStream
  << std::endl << "BLOCK LHCFASERFORJL" << std::endl
  << "# sqrt(s)   ID    lower value        upper value" << std::endl;
  for( int whichSignal( 1 );
       sevenTevSignals.size() > whichSignal;
       ++whichSignal )
  {
    currentSignal = sevenTevSignals[ whichSignal ];
    outputStream
    << "       7     " << whichSignal << "   "
    << CppSLHA::CppSLHA_global::SLHA_double( currentSignal->getLowerValue()
                                             * normalizationFactor )
    << "   "
    << CppSLHA::CppSLHA_global::SLHA_double( currentSignal->getUpperValue()
                                             * normalizationFactor )
    << "   # " << *(currentSignal->getName()) << " / "
    << *(normalizationSignal->getName())
    << std::endl;
  }
  normalizationSignal = fourteenTevSignals[ 0 ];
  normalizationFactor = ( 1.0 / normalizationSignal->getValue() );
  for( int whichSignal( 1 );
       fourteenTevSignals.size() > whichSignal;
       ++whichSignal )
  {
    currentSignal = fourteenTevSignals[ whichSignal ];
    outputStream
    << "      14     " << whichSignal << "   "
    << CppSLHA::CppSLHA_global::SLHA_double( currentSignal->getLowerValue()
                                             * normalizationFactor )
    << "   "
    << CppSLHA::CppSLHA_global::SLHA_double( currentSignal->getUpperValue()
                                             * normalizationFactor )
    << "   # " << *(currentSignal->getName()) << " / "
    << *(normalizationSignal->getName())
    << std::endl;
  }
  outputStream << std::endl;
  outputStream.close();

  std::cout
  << std::endl
  << "finished!";
  std::cout << std::endl;

  // this was a triumph! I'm making a note here:
  return EXIT_SUCCESS;
}
