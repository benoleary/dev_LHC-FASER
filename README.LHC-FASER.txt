/*
 * README.LHC-FASER.txt
 *
 *  Created on: 31 Mar 2011
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

/* Usage Instructions!
 * 
 * Currently LHC-FASER is in alpha testing. The code performs the algorithm
 * correctly, as far as tests have shown; however, the accuracy of the
 * approximations used by the algorithm has not yet been studied.
 *
 * The intended usage of LHC-FASER is that an instance of the
 * LHC_FASER::lhcFaser class is created (there are several options which can be
 * set in the constructor, such as where to find the lookup grids and what unit
 * of cross-section to use for the calculated event rates), then that the
 * lhcFaser::addSignal( std::string const signalName ) function is called,
 * perhaps several times.
 * LHC_FASER::signalHandler pointers can be used to store the return value of
 * lhcFaser::addSignal( std::string const signalName ) so that the rate
 * calculated for this signal can be directly accessed. Otherwise, 
 * lhcFaser::getSignal( std::string const signalName ) has to be used to find
 * it again.
 * After this initial setup, each time that the input file changes,
 * lhcFaser::updateForNewSlha() should be called (or
 * updateForUpdatedSpectrum(), if the lhcFaser is getting its information from
 * an external CppSLHA object rather than reading in the SLHA file itself),
 * which flags all cached values as needing to be recalculated. If this
 * function is not called, no matter what happens to the SLHA file that
 * originally provided the data, the lhcFaser object will continue to return
 * cached values for the original input.
 * Once the data file has been (re)read in, the LHC_FASER::signalHandler
 * pointers can return the calculated event rates with
 * signalHandler::getValue(), which returns a double.
 *
 * example code would be:
 
// initial code stuff
 
#include "/path/to/LHC-FASER/LHC-FASER.hpp" // remember to include the header!
 
// more preamble
 
// code...
 
// somewhere in the main(...) function:
LHC_FASER::lhcFaser
exampleLhcFaser( "~/exampleWorkingDirectory/exampleSlhaFile.txt",
                 "~/exampleWorkingDirectory/exampleAlternativeGridDirectory/",
                 "pb" );
 /* exampleLhcFaser now reads in ~/exampleWorkingDirectory/exampleSlhaFile.txt
  * for its data, looks for its grids in
  * ~/exampleWorkingDirectory/exampleAlternativeGridDirectory/, and returns its
  * values in units of picobarns of cross-section.
  */
 LHC_FASER::signalHandler*
 exampleRateForAtlas4j0l( exampleLhcFaser.addSignal( "Atlas4jMET_0l_7TeV" ) );
 LHC_FASER::signalHandler*
 exampleRateForAtlas3j1l( exampleLhcFaser.addSignal( "Atlas3jMET_1l_7TeV" ) );
 /* exampleRateForAtlas4j0l is now a pointer to an LHC_FASER::signalHandler
  * which calculates the rate of SUSY events which pass the Atlas 4-jet, MET,
  * 0-lepton cuts based on the spectrum found in
  * ~/exampleWorkingDirectory/exampleSlhaFile.txt when
  * exampleRateForAtlas4j0l->getValue() is called.
  * likewise exampleRateForAtlas3j1l is for the Atlas 3-jet, MET, 1-lepton cuts
  * based on ~/exampleWorkingDirectory/exampleSlhaFile.txt
  */
  
  // more code
  
  bool notYetFinished( true );
  double exampleEventRateHolder;
  double exampleCombinedRate;
  for( int numberOfPointsStillToCheck( initialNumberOfPointsToCheck );
       0 < numberOfPointsStillToCheck;
      --numberOfPointsStillToCheck )
  {
    /* more code, which presumably modifies the contents of
     * ~/exampleWorkingDirectory/exampleSlhaFile.txt, possibly by re-writing
     * it, though it's not important how it changed, just that it changed.
     */
  
    exampleLhcFaser.updateForNewSlha();
    // exampleLhcFaser has now read in the new SLHA data.
  
    exampleEventRateHolder = exampleRateForAtlas4j0l->getValue();
  
    // code that does something based on the event rate now stored in
    // exampleEventRateHolder
  
    exampleEventRateHolder = exampleRateForAtlas3j1l->getValue();
  
    // code that does something based on the event rate now stored in
    // exampleEventRateHolder
  
    // more code
    
    /* oops, we forgot that we still need the Atlas 4-jet, MET, 0-lepton rate.
     * luckily exampleRateForAtlas4j0l->getValue() doesn't re-calculate the
     * value, taking up the CPU time again, since
     * exampleLhcFaser.updateForNewSlha() has not been called since last we
     * called exampleRateForAtlas4j0l->getValue(); it just returns the cached
     * value.
     */
     exampleCombinedRate
     = ( exampleEventRateHolder + exampleRateForAtlas4j0l->getValue() );
  
    // more code
    
  }
  
// more code
  
/* The example code above should demonstrate all the functions required for the
 * intended normal use.
 */

