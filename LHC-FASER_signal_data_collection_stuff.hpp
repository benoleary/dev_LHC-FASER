/*
 * LHC-FASER_signal_data_collection_stuff.hpp
 *
 *  Created on: 17 Dec 2010
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
 *      LHC-FASER_electroweak_cascade_stuff.hpp
 *      LHC-FASER_electroweak_cascade_stuff.cpp
 *      LHC-FASER_full_cascade_stuff.hpp
 *      LHC-FASER_full_cascade_stuff.cpp
 *      LHC-FASER_global_stuff.hpp
 *      LHC-FASER_global_stuff.cpp
 *      LHC-FASER_input_handling_stuff.hpp
 *      LHC-FASER_input_handling_stuff.cpp
 *      LHC-FASER_kinematics_stuff.hpp
 *      LHC-FASER_kinematics_stuff.cpp
 *      LHC-FASER_lepton_distributions.hpp
 *      LHC-FASER_lepton_distributions.cpp
 *      LHC-FASER_signal_calculator_stuff.hpp
 *      LHC-FASER_signal_calculator_stuff.cpp
 *      LHC-FASER_signal_data_collection_stuff.hpp
 *      LHC-FASER_signal_data_collection_stuff.cpp
 *      LHC-FASER_sparticle_decay_stuff.hpp
 *      LHC-FASER_sparticle_decay_stuff.cpp
 *      and README.LHC-FASER.txt which describes the package.
 *
 *      LHC-FASER also requires CppSLHA. It should be found in a subdirectory
 *      included with this package.
 *
 *      LHC-FASER also requires grids of lookup acceptanceValues. These should also be
 *      found in a subdirectory included with this package.
 */

#ifndef LHC_FASER_SIGNAL_DATA_COLLECTION_STUFF_HPP_
#define LHC_FASER_SIGNAL_DATA_COLLECTION_STUFF_HPP_

#include "LHC-FASER_input_handling_stuff.hpp"
#include "LHC-FASER_cross-section_stuff.hpp"
#include "LHC-FASER_jet_kinematics_stuff.hpp"
#include "LHC-FASER_full_cascade_stuff.hpp"

namespace LHC_FASER
{
  /* this class holds a pointer to an inputHandler instance along with
   * various pointers to cross-section, kinematics, & cascade handlers, for
   * ease of reference. it also keeps an instance of the sQCD_to_EWino_handler
   * class.
   */
  class signalShortcuts
  {
  public:
    signalShortcuts( inputHandler const* const inputShortcuts,
                     crossSectionHandler* const crossSections,
                     jetPlusMetAcceptanceHandler* const jetPlusMetAcceptances,
                     //leptonAcceptanceHandler* const leptonAcceptances,
                     fullCascadeSetFactory* const cascadeSets )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~signalShortcuts()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    inputHandler const*
    getInputShortcuts()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    crossSectionHandler*
    getCrossSections()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    jetPlusMetAcceptanceHandler*
    getJetPlusMetAcceptances()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    //leptonAcceptanceHandler*
    //getLeptonAcceptances()
    //const
    /* code after the classes in this .hpp file, or in the .cpp file. *///;
    fullCascadeSetFactory*
    getCascadeSets()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getUncertainty()
    /* this looks at the various mass differences in the hierarchy & then
     * guesses how uncertain the signal value is (as a multiplicative factor
     * > 1.0).
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    inputHandler const* const inputShortcuts;
    crossSectionHandler* const crossSections;
    jetPlusMetAcceptanceHandler* const jetPlusMetAcceptances;
    //leptonAcceptanceHandler* const leptonAcceptances;
    fullCascadeSetFactory* const cascadeSets;
  };  // end of signalShortcuts class.


  /* this holds together the various doubles & pointers that store information
   * that is particular to a signal. it doesn't hold all the information
   * required by every single signal type, & it does hold a little bit more
   * than is necessary for the least elaborate signal, but it should be a
   * convenient package of information. it builds on the acceptanceCutSet
   * class.
   */
  class signalDefinitionSet : public acceptanceCutSet
  {
  public:
    signalDefinitionSet( signalShortcuts* const shortcut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    signalDefinitionSet( signalDefinitionSet* const copySource )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    virtual
    ~signalDefinitionSet()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    signalShortcuts*
    getShortcuts()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setBeamEnergy( int const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    crossSectionTableSet*
    getCrossSections()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    jetAcceptanceTable*
    getJetPlusMetAcceptance()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setJetPlusMetAcceptance( jetAcceptanceTable* const inputTable )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setExcludedStandardModelProducts(
                                  std::list< int > const* const exclusionList )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    fullCascadeSetsForOneBeamEnergy*
    getCascadeSets()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    signalShortcuts* const shortcut;

    // each signal needs to look up specific tables, based on the beam energy:
    crossSectionTableSet* crossSections;
    jetAcceptanceTable* jetPlusMetAcceptance;
    // by default each signal has only a single jet+MET table to look up; if a
    // signal needs to look up more tables, it should store the extras itself.
    std::list< int > const* exclusionList;
    fullCascadeSetsForOneBeamEnergy* cascadeSets;
  };  // end of signal_definitions class.


  /* this holds each production channel for a combination of colored sparticles
   * as a set of pointers to the crossSectionTable (which holds a pointer to
   * the signedParticleShortcutPair) & the relevant fullCascadeSets for the
   * sparticles.
   */
  class productionChannelPointerSet
  {
  public:
    productionChannelPointerSet( signalDefinitionSet* const signalDefinitions,
                         signedParticleShortcutPair const* const scoloredPair )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~productionChannelPointerSet()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    signedParticleShortcutPair const*
    getScoloredPair()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getCrossSection()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    fullCascadeSet*
    getFirstCascadeSet()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    fullCascadeSet*
    getSecondCascadeSet()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    //signalShortcuts* const shortcut;
    crossSectionTable* crossSection;
    fullCascadeSet* firstCascadeSet;
    fullCascadeSet* secondCascadeSet;
  };  // end of productionChannelPointerSet class.





  // inline functions:


  inline inputHandler const*
  signalShortcuts::getInputShortcuts()
  const
  {
    return inputShortcuts;
  }

  inline crossSectionHandler*
  signalShortcuts::getCrossSections()
  const
  {
    return crossSections;
  }

  inline jetPlusMetAcceptanceHandler*
  signalShortcuts::getJetPlusMetAcceptances()
  const
  {
    return jetPlusMetAcceptances;
  }

  /*inline leptonAcceptanceHandler*
  signalShortcuts::getLeptonAcceptances()
  const
  {
    return leptonAcceptances;
  }*/

  inline fullCascadeSetFactory*
  signalShortcuts::getCascadeSets()
  const
  {
    return cascadeSets;
  }



  inline signalShortcuts*
  signalDefinitionSet::getShortcuts()
  {
    return shortcut;
  }

  inline void
  signalDefinitionSet::setBeamEnergy( int const inputValue )
  {
    beamEnergy = inputValue;
    crossSections = shortcut->getCrossSections()->getTableSet( inputValue );
    cascadeSets
    = shortcut->getCascadeSets()->getFullCascadeSetsForOneBeamEnergy(
                                                                  inputValue );
  }

  inline crossSectionTableSet*
  signalDefinitionSet::getCrossSections()
  {
    return crossSections;
  }

  inline jetAcceptanceTable*
  signalDefinitionSet::getJetPlusMetAcceptance()
  {
    return jetPlusMetAcceptance;
  }

  inline void
  signalDefinitionSet::setJetPlusMetAcceptance(
                                         jetAcceptanceTable* const inputTable )
  {
    jetPlusMetAcceptance = inputTable;
  }

  inline void
  signalDefinitionSet::setExcludedStandardModelProducts(
                                  std::list< int > const* const exclusionList )
  {
    this->exclusionList = exclusionList;
  }

  inline fullCascadeSetsForOneBeamEnergy*
  signalDefinitionSet::getCascadeSets()
  {
    return cascadeSets;
  }



  inline signedParticleShortcutPair const*
  productionChannelPointerSet::getScoloredPair()
  const
  {
    return crossSection->getPair();
  }

  inline double
  productionChannelPointerSet::getCrossSection()
  {
    return crossSection->getValue();
  }

  inline fullCascadeSet*
  productionChannelPointerSet::getFirstCascadeSet()
  {
    return firstCascadeSet;
  }

  inline fullCascadeSet*
  productionChannelPointerSet::getSecondCascadeSet()
  {
    return secondCascadeSet;
  }

}  // end of LHC_FASER namespace

#endif /* LHC_FASER_SIGNAL_DATA_COLLECTION_STUFF_HPP_ */
