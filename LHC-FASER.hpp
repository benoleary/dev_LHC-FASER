/*
 * LHC-FASER.hpp
 *
 *  Created on: 06 Oct 2010
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

#ifndef LHC_FASER_HPP_
#define LHC_FASER_HPP_

/* this is a C++ class for calculating event rates for signals of the Minimal
 * Supersymmetric Standard Model at the LHC. it calculates only the signal
 * event rate, assuming that the background is known.
 */

/* this version of the code was written mainly by Ben O'Leary ("BOL",
 * benjamin.oleary@gmail.com), based on a setup decided upon by Michael Kr√§mer
 * ("MK", mkraemer@physik.rwth-aachen.de), Jonas Lindert ("JL",
 * jonas.lindert@googlemail.com) and BOL. Quite a lot is a refactoring of
 * JL's code. the cross-section grids were created by JL, while the current
 * acceptance grids were written by Carsten Robens
 * ("CR", Carsten.Robens@rwth-aachen.de).
 */

/* the method of calculation is as laid out in the publication
 * "SUSY parameter determination at the LHC using cross sections and kinematic
 * edges" by Herbi Dreiner, Michael Krämer, Jonas Lindert and Ben O'Leary,
 * published in JHEP 1004 (2010) 109 available on the arXiv as number
 * 1003.2648 as well.
 * other publications *should* follow, I hope.
 */

// this particular file is the header file for the main class. there is an
// accompanying .cpp file with the functions.


// include the header file with useful global-ish stuff, as well as including
// other general headers such as iostream:
#include "LHC-FASER_global_stuff.hpp"
// include the classes for handling signals (hence also those for handling the
// data collections for the signals):
#include "LHC-FASER_signal_calculator_stuff.hpp"


namespace LHC_FASER
{
  /* this is a class for calculating MSSM signal event rates at the LHC.
   * it tries to be a Factory (creating appropriate subtypes of signals based
   * on input).
   */
  class lhcFaser
  {
  public:
    /* constructor with specified CppSLHA & user-defined path to grids,
     * whether returned event rates should be in nb, pb or fb& whether to use
     * LO or NLO:
     */
    lhcFaser( CppSLHA::CppSLHA0* const spectrumData,
              std::string const pathToGrids = "./grids/",
              std::string const crossSectionUnit = "fb",
              bool const usingNlo = true )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    /* constructor with specified SLHA file & user-defined path to grids,
     * whether returned event rates should be in nb, pb or fb& whether to use
     * LO or NLO:
     */
    lhcFaser( std::string const slhaFileName,
              std::string const pathToGrids = "./grids/",
              std::string const crossSectionUnit = "fb",
              bool const usingNlo = true )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~lhcFaser()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    signalHandler*
    addSignal( std::string const signalName )
    // this adds a new signal to the set of signals based on its name.
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    signalHandler*
    getSignal( std::string const signalName )
    // this returns the handler object for the requested signal name.
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    updateForUpdatedSpectrum()
    /* this assumes that the CppSLHA was updated & so sets each signal to be
     * recalculated next time its value is requested. I expect that I could do
     * this more elegantly with throwing exceptions, but I'll leave that for
     * another day.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    updateForNewSlha()
    // this reads in the CppSLHA's target file & recalculates all required
    // signals.
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    updateForNewSlha( std::string const slhaFileName )
    // this reads in the new file & recalculates all required signals.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    CppSLHA::CppSLHA0* spectrumData;
    // this holds the MSSM spectrum data.
    bool usingOwnCppSlha;
    // this notes whether the CppSLHA was supplied externally (in which case
    // the destructor does NOT delete the CppSLHA).
    std::string pathToGrids;
    // this is where the lookup tables live.
    std::string pathToJetPlusMetAcceptanceGrids;
    // this is where the lookup tables for jet+MET acceptances live.
    double crossSectionUnitFactor;
    // this is to allow for the user to specify event rates in fb, pb or nb.
    bool usingNlo;
    // this is to allow the user to use LO or NLO cross-sections.

    inputHandler* inputSource;
    // this keeps const pointers to useful objects together for ease of passing
    // around & for neater code.
    crossSectionHandler* crossSectionSource;
    // this holds the lookup tables for LHC colored sparticle pair production
    // cross-sections.
    jetPlusMetAcceptanceHandler* jetPlusMetAcceptanceSource;
    // this holds the kinematic data lookup tables for colored sparticle
    // production at the LHC.
    electroweakCascadeHandler* electroweakCascadeSource;
    fullCascadeSetFactory* fullCascadeSetSource;
    // this handles the cascade decays of colored sparticles.
    signalDefinitionSet* signalPreparationDefinitions;
    std::vector< signalHandler* > signalSet;
    // this tracks the various signals of LHC supersymmetric events & their
    // rates.

    signalShortcuts* shortcut;
    // this keeps const pointers to useful objects together for ease of passing
    // around & for neater code.
    readierForNewPoint* readier;

    void
    initialize( CppSLHA::CppSLHA0* const spectrumData,
                std::string const pathToGrids,
                std::string const crossSectionUnit,
                bool const usingNlo )
    // this is used by all the constructors to do most of the construction.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };  // end of lhcFaser class.




  // inline functions:


  inline signalHandler*
  lhcFaser::addSignal( std::string const signalName )
  // this adds a new signal to the set of signals based on its name.
  {
    signalSet.push_back( new signalHandler( signalName,
                                            crossSectionUnitFactor,
                                            signalPreparationDefinitions ) );
    return signalSet.back();
  }

  inline signalHandler*
  lhcFaser::getSignal( std::string const signalName )
  // this returns the handler object for the requested signal name.
  {
    signalHandler* returnPointer( NULL );
    for( std::vector< signalHandler* >::iterator
         signalIterator( signalSet.begin() );
         signalSet.end() > signalIterator;
         signalIterator++ )
      // look through all the signals...
      {
        if( 0 == (*signalIterator)->getName()->compare( signalName ) )
          // if we find the requested signal...
          {
            returnPointer = *signalIterator;
            signalIterator = signalSet.end();
            // stop looking.
          }
      }
    return returnPointer;
  }

  inline void
  lhcFaser::updateForUpdatedSpectrum()
  /* this assumes that the CppSLHA was updated & so sets each signal to be
   * recalculated next time its value is requested. I expect that I could do
   * this more elegantly with throwing exceptions, but I'll leave that for
   * another day.
   */
  {
    readier->readyObserversForNewPoint();
  }

  inline void
  lhcFaser::updateForNewSlha()
  // this reads in the CppSLHA's target file & recalculates all required
  // signals.
  {
    spectrumData->read_file();
    updateForUpdatedSpectrum();
  }

  inline void
  lhcFaser::updateForNewSlha( std::string const slhaFileName )
  // this reads in the new file & recalculates all required signals.
  {
    spectrumData->read_file( slhaFileName );
    updateForUpdatedSpectrum();
  }

}  // end of LHC_FASER namespace.

#endif /* LHC_FASER_HPP_ */
