/*
 * LHC-FASER_electroweak_cascade_collection_stuff.hpp
 *
 *  Created on: 26 Aug 2010
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

#ifndef LHC_FASER_ELECTROWEAK_CASCADE_COLLECTION_STUFF_HPP_
#define LHC_FASER_ELECTROWEAK_CASCADE_COLLECTION_STUFF_HPP_

#include "LHC-FASER_neutral_electroweak_cascade_stuff.hpp"
#include "LHC-FASER_charged_electroweak_cascade_stuff.hpp"

namespace LHC_FASER
{
  /* this is a class to hold a set of electroweakCascadeSource associated with
   * a single decaying electroweakino or EWSB boson (from the decay of a
   * single colored sparticle). the decaying colored sparticle is assumed to be
   * the version with a PDG code which is POSITIVE, rather than the
   * antiparticle with negative code. hence for example an instance of this
   * class would be for a positively-charged up squark which decays into a
   * negatively-charged down quark & a positively-charged chargino, or for
   * another example an instance could be for a negatively-charged down squark,
   * which decays into a positively-charged up quark & a negatively-charged
   * chargino (which has negative PDG code!). the charge-conjugate process,
   * which is assumed to have the same branching ratios & acceptances, though
   * for opposite-charge SM fermions, must be taken care of by whatever objects
   * are using this class.
   *
   * there are 2 different constructors: 1 for electroweakino cascades, the
   * other for EWSB boson cascades. I could use derived classes here, but I
   * think that it's not worth it.
   */
  class electroweakCascadeSet
  {
  public:
    electroweakCascadeSet( leptonAcceptanceParameterSet* const kinematicsTable,
                           particlePointer const coloredDecayer,
                           particlePointer const electroweakDecayer,
                           inputHandler const* const inputShortcut,
                effectiveSquarkMassHolder* effectiveSquarkMassPointer = NULL );
    // this is the version for electroweakinos.
    electroweakCascadeSet( leptonAcceptanceParameterSet* const kinematicsTable,
                           particlePointer const coloredDecayer,
                           particlePointer const electroweakDecayer,
                           particlePointer const lighterScolored,
                           inputHandler const* const inputShortcut,
                effectiveSquarkMassHolder* effectiveSquarkMassPointer = NULL );
    // this is the version for EWSB bosons.
    virtual
    ~electroweakCascadeSet();

    bool
    isEquivalent( particlePointer const coloredDecayer,
                  particlePointer const electroweakDecayer );
    // this returns true if the coloredDecayers & electroweakDecayers match.
    bool
    isEquivalent( particlePointer const coloredDecayer,
                  particlePointer const electroweakDecayer,
                  particlePointer const lighterScolored );
    // this returns true if the coloredDecayers, electroweakDecayers, &
    // lighterScoloreds match.
    particlePointer
    getColoredDecayer()
    const;
    particlePointer
    getElectroweakDecayer()
    const;
    particlePointer
    getLighterScolored()
    const;
    double
    getAcceptance( acceptanceCutSet const* const acceptanceCuts,
                   int const numberOfAdditionalJets,
                   int const numberOfNegativeElectrons,
                   int const numberOfPositiveElectrons,
                   int const numberOfNegativeMuons,
                   int const numberOfPositiveMuons );
    /* this does the common job of checking to see if the point has been
     * updated before calling the relevant protected virtual function, which
     * returns the set of values for acceptance * branching ratio FOR THE
     * CASCADES FROM electroweakDecayer ONWARDS (does NOT include the BR for
     * coloredDecayer to electroweakDecayer - this is dealt with by the
     * cutSpecifiedFullCascade, because it depends on what flavors of quark jet
     * are being considered for the signal).
     * this function returns 0.0 if the cascades never result in the requested
     * set of jets + leptons.
     * the number of jets is the number produced from the decay of the
     * decaying *electroweak* particle (neutralino, chargino or vector boson),
     * NOT from the decay of the colored sparticle to the decaying electroweak
     * particle.
     */
    double
    getOssfMinusOsdf( acceptanceCutSet const* const acceptanceCuts );
    /* this does the common job of checking to see if the point has been
     * updated before calling the relevant protected virtual function, which
     * returns the value for acceptance * branching ratio FOR THE
     * CASCADES FROM electroweakDecayer ONWARDS (does NOT include the BR for
     * coloredDecayer to electroweakDecayer - this is dealt with by the
     * cutSpecifiedFullCascade, because it depends on what flavors of quark jet
     * are being considered for the signal).
     * this function returns 0.0 if the cascades never result in an OSSF - OSDF
     * signal.
     */


  protected:
    leptonAcceptanceParameterSet* const kinematicsTable;
    inputHandler const* const inputShortcut;
    particlePointer const coloredDecayer;
    particlePointer const electroweakDecayer;
    particlePointer const lighterScolored;
    std::vector< electroweakCascade* > cascadeSet;
    electroweakCascade* currentCascade; // this is used for filling cascadeSet.
  };


  // this gives out pointers to electroweakCascadeSets at a fixed beam energy
  // based on the requested colored sparticle & electroweakino or vector boson.
  class electroweakCascadesForOneBeamEnergy
  {
  public:
    electroweakCascadesForOneBeamEnergy(
                               leptonAcceptanceHandler* const kinematicsSource,
                                         int const beamEnergy,
                                     inputHandler const* const inputShortcut );
    ~electroweakCascadesForOneBeamEnergy();

    int
    getBeamEnergy()
    const;
    electroweakCascadeSet*
    getCascadeSet( particlePointer const coloredDecayer,
                   particlePointer const electroweakino );
    /* this looks to see if it already has an electroweakCascadeSet
     * corresponding to the requested pairing, & if it does, it returns a
     * pointer to it, & if it doesn't, it constructs a new
     * electroweakCascadeSet & returns a pointer to that.
     */
    electroweakCascadeSet*
    getCascadeSet( particlePointer const coloredDecayer,
                   particlePointer const electroweakBoson,
                   particlePointer const lighterScolored );
    /* this looks to see if it already has an electroweakCascadeSet
     * corresponding to the requested particles, & if it does, it returns a
     * pointer to it, & if it doesn't, it constructs a new
     * electroweakCascadeSet & returns a pointer to that.
     */
    electroweakCascadeSet*
    getCascadeSet( particlePointer const coloredDecayer,
                   particlePointer const electroweakino,
                   effectiveSquarkMassHolder* const effectiveSquarkMass );
    /* this looks to see if it already has an electroweakCascadeSet
     * corresponding to the requested particles, & if it does, it returns a
     * pointer to it, & if it doesn't, it constructs a new
     * electroweakCascadeSet & returns a pointer to that.
     */
    electroweakCascadeSet*
    getCascadeSet( particlePointer const coloredDecayer,
                   particlePointer const electroweakBoson,
                   particlePointer const lighterScolored,
                   effectiveSquarkMassHolder* const effectiveSquarkMass );
    /* this looks to see if it already has an electroweakCascadeSet
     * corresponding to the requested particles, & if it does, it returns a
     * pointer to it, & if it doesn't, it constructs a new
     * electroweakCascadeSet & returns a pointer to that.
     */


  protected:
    leptonAcceptancesForOneBeamEnergy* const kinematicsTable;
    int const beamEnergy;
    inputHandler const* const inputShortcut;
    std::vector< electroweakCascadeSet* > electroweakinoCascadeSets;
    std::vector< electroweakCascadeSet* > bosonCascadeSets;
    std::vector< electroweakCascadeSet* > virtualSdownCascadeSets;
    std::vector< electroweakCascadeSet* > virtualTopCascadeSets;
    electroweakCascadeSet* currentCascadeSet;
    // this is used for filling electroweakinoCascadeSets.
  };


  // this class holds all the electroweakCascadesForOneBeamEnergy & passes out
  // pointers to them.
  class electroweakCascadeHandler
  {
  public:
    electroweakCascadeHandler( inputHandler const* const inputShortcut );
    ~electroweakCascadeHandler();

    electroweakCascadesForOneBeamEnergy*
    getElectroweakCascadesForOneBeamEnergy( int const beamEnergy );
    /* this looks to see if there is an existing
     * leptonAcceptancesForOneBeamEnergy with the requested values, & if not,
     * makes 1, & returns the pointer.
     */


  protected:
    leptonAcceptanceHandler kinematicsSource;
    inputHandler const* const inputShortcut;
    std::vector< electroweakCascadesForOneBeamEnergy* >
    cascadeSetAtBeamEnergies;
  };





  // inline functions:


  inline bool
  electroweakCascadeSet::isEquivalent( particlePointer const coloredDecayer,
                                     particlePointer const electroweakDecayer )
  // this returns true if the coloredDecayers & electroweakDecayers match.
  {
    if( ( coloredDecayer == this->coloredDecayer )
        &&
        ( electroweakDecayer == this->electroweakDecayer )
        &&
        ( NULL == lighterScolored ) )
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  inline bool
  electroweakCascadeSet::isEquivalent( particlePointer const coloredDecayer,
                                     particlePointer const electroweakDecayer,
                                       particlePointer const lighterScolored )
  // this returns true if the coloredDecayers, electroweakDecayers, &
  // lighterScoloreds match.
  {
    if( ( coloredDecayer == this->coloredDecayer )
        &&
        ( electroweakDecayer == this->electroweakDecayer )
        &&
        ( lighterScolored == this->lighterScolored ) )
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  inline particlePointer
  electroweakCascadeSet::getColoredDecayer()
  const
  {
    return coloredDecayer;
  }

  inline particlePointer
  electroweakCascadeSet::getElectroweakDecayer()
  const
  {
    return electroweakDecayer;
  }

  inline particlePointer
  electroweakCascadeSet::getLighterScolored()
  const
  {
    return lighterScolored;
  }

  inline double
  electroweakCascadeSet::getAcceptance(
                                  acceptanceCutSet const* const acceptanceCuts,
                                        int const numberOfAdditionalJets,
                                        int const numberOfNegativeElectrons,
                                        int const numberOfPositiveElectrons,
                                        int const numberOfNegativeMuons,
                                        int const numberOfPositiveMuons )
  // this adds up the acceptances from the stored electroweakCascadeSource.
  {
    double returnValue( 0.0 );
    for( std::vector< electroweakCascade* >::iterator
         cascadeIterator( cascadeSet.begin() );
         cascadeSet.end() > cascadeIterator;
         ++cascadeIterator )
    {
      returnValue += (*cascadeIterator)->getAcceptance( acceptanceCuts,
                                                        numberOfAdditionalJets,
                                                     numberOfNegativeElectrons,
                                                     numberOfPositiveElectrons,
                                                        numberOfNegativeMuons,
                                                       numberOfPositiveMuons );
    }
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "electroweakCascadeSet::getAcceptance( "
    << acceptanceCuts << ", " << numberOfAdditionalJets << ", "
    << numberOfNegativeElectrons << ", " << numberOfPositiveElectrons << ", "
    << numberOfNegativeMuons << ", " << numberOfPositiveMuons
    << " ) [ coloredDecayer = ";
    if( NULL != coloredDecayer )
    {
      std::cout << *(coloredDecayer->get_name());
    }
    else
    {
      std::cout << "NULL";
    }
    std::cout << ", electroweakDecayer = ";
    if( NULL != electroweakDecayer )
    {
      std::cout << *(electroweakDecayer->get_name());
    }
    else
    {
      std::cout << "NULL";
    }
    std::cout << " ] returning " << returnValue;
    std::cout << std::endl;**/

    return returnValue;
  }

  inline double
  electroweakCascadeSet::getOssfMinusOsdf(
                                 acceptanceCutSet const* const acceptanceCuts )
  /* this does the common job of checking to see if the point has been
   * updated before calling the relevant protected virtual function, which
   * returns the value for acceptance * branching ratio FOR THE
   * CASCADES FROM electroweakDecayer ONWARDS (does NOT include the BR for
   * coloredDecayer to electroweakDecayer - this is dealt with by the
   * cutSpecifiedFullCascade, because it depends on what flavors of quark jet
   * are being considered for the signal).
   * this function returns 0.0 if the cascadeSet never results in an
   * OSSF - OSDF signal.
   */
  {
    double returnValue( 0.0 );
    for( std::vector< electroweakCascade* >::iterator
         cascadeIterator( cascadeSet.begin() );
         cascadeSet.end() > cascadeIterator;
         ++cascadeIterator )
    {
      returnValue += (*cascadeIterator)->getOssfMinusOsdf( acceptanceCuts );
    }
    return returnValue;
  }



  inline int
  electroweakCascadesForOneBeamEnergy::getBeamEnergy()
  const
  {
    return beamEnergy;
  }

}

#endif /* LHC_FASER_ELECTROWEAK_CASCADE_COLLECTION_STUFF_HPP_ */
