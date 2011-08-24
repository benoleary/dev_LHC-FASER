/*
 * LHC-FASER_cross-section_stuff.hpp
 *
 *  Created on: 15 Nov 2010
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
 *      LHC-FASER also requires grids of lookup values. These should also be
 *      found in a subdirectory included with this package.
 */

#ifndef LHC_FASER_CROSS_SECTION_STUFF_HPP_
#define LHC_FASER_CROSS_SECTION_STUFF_HPP_

#include "LHC-FASER_global_stuff.hpp"
#include "LHC-FASER_input_handling_stuff.hpp"

namespace LHC_FASER
{
  /* this class reads in a file in the assumed format, stores it, & gives out
   * interpolated values. it was written with cross-sections in mind, with data
   * files in the format
   * squark_mass gluino_mass cross-section newline
   * or K-factors, in a similar format.
   */
  class squareGrid
  {
  public:
    squareGrid( std::string const* const gridFileLocation,
                std::string const* const gridName,
                squareGrid const* const scalingGrid )
    /* this constructor reads in a grid file, assumed to be in the format
     * x_coordinate y_coordinate value
     * in ascending order, y_coordinate varying first
     * (e.g.
     * 200.0 200.0 123.4
     * 200.0 210.0 567.8
     * ...
     * 210.0 200.0 876.5
     * 210.0 210.0 432.1
     * ...).
     * if provided with a non-NULL pointer to another squareGrid, it stores
     * the its values scaled by an interpolated value from this other
     * squareGrid.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    squareGrid( std::string baseGridFileLocation,
                std::string gridName,
                std::string scalingGridFileLocation )
    // this constructor uses the other constructor to make another squareGrid
    // with the scaling factors, which is then used to construct this instance.
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~squareGrid()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    std::string const*
    getName()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    valueAt( double const xCoordinate,
             double const yCoordinate )
    const
    /* this finds the grid square which the given point is in, & then uses
     * lhcFaserGlobal::squareBilinearInterpolation to get an interpolated
     * value.
     * N.B.: cross-section grids use
     * xCoordinate == SQUARK mass, yCoordinate == GLUINO mass!
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getLowestX()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getHighestX()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getLowestY()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getHighestY()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    std::string gridName;
    double gridStepSize;
    double lowestXCoordinate;
    double highestXCoordinate;
    double lowestYCoordinate;
    double highestYCoordinate;
    std::vector< std::vector< double >* > values;

    void
    readIn( std::string const* const gridFileLocation,
            squareGrid const* const scalingGrid )
    // this actually does most of the job of the constructors.
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };


  /* this class holds a squareGrid instance as a cross-section lookup table
   * for a particular colored sparticle combination, & returns the
   * interpolated value by looking up the masses of its sparticles.
   * this class is small enough to be kept entirely in this header file.
   * N.B.: currently, this class assumes that only mass-degenerate
   * squark-squark or squark-antisquark pairs are produced. this should be fine
   * at the accuracy we work at unless there are large splittings among the
   * sdown & sup mass eigenstates.
   */
  class crossSectionTable : public getsReadiedForNewPoint
  {

  public:

    crossSectionTable( squareGrid const* const lookupTable,
                       signedParticleShortcutPair const* const scoloredPair,
                       double const flavorFactor,
                       inputHandler const* const shortcut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~crossSectionTable()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    signedParticleShortcutPair const*
    getPair()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getValue()
    /* this checks to see if it needs to recalculate the value for this point,
     * & if so, it does, & if the point lies outside the grid, it takes the
     * value of the nearest point on the grid instead.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    signedParticleShortcutPair const* scoloredPair;
    CppSLHA::particle_property_set const* directlyProducedSquark;
    inputHandler const* const shortcut;
    double storedValue;
    squareGrid const* const lookupTable;
    double const flavorFactor;
    /* this is used when the flavor of squark for this table uses the grid for
     * a different flavor of squark, in the approximation that it has the same
     * cross-section scaled by some factor (e.g. the gluino+sup cross-section
     * is twice the gluino+sdown cross-section for the same squark mass).
     */
  };


  /* this class holds several squareGrid instances as cross-section lookup
   * tables for various colored sparticle combinations, & returns their
   * interpolated values for given signedParticleShortcutPair pointers.
   */
  class crossSectionTableSet
  {

  public:

    crossSectionTableSet( std::string const* const gridDirectory,
                          int const beamEnergy,
                          inputHandler const* const shortcut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~crossSectionTableSet()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    int
    getBeamEnergy()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    crossSectionTable*
    getTable( signedParticleShortcutPair const* const requestedChannel )
    // this returns the crossSectionTable for the requested pair.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    inputHandler const* const shortcut;
    std::string const gridDirectory;
    int const beamEnergy;
    std::vector< squareGrid* > grids;
    // this holds squareGrid instances with unique names.
    std::vector< crossSectionTable* > tables;
    // this holds crossSectionTable instances with appropriate
    // signedParticleShortcutPair pointers.

    double
    prepareGridName( std::string* const gridName,
                     signedParticleShortcutPair const* const scoloredPair )
    // this returns the flavor factor for the requested pair while putting the
    // squareGrid's associated string in gridName.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    squareGrid const*
    getGrid( std::string const* const gridName )
    // this either finds the appropriate squareGrid in grids or makes a new
    // instance, stores it, & returns its pointer.
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };


  class crossSectionHandler
  {
  public:
    crossSectionHandler( inputHandler const* const shortcut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~crossSectionHandler()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    crossSectionTableSet*
    getTableSet( int const beamEnergyInTev )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    crossSectionTable*
    getTable( int const beamEnergyInTev,
              signedParticleShortcutPair const* const requestedChannel )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    inputHandler const* const shortcut;
    std::vector< crossSectionTableSet* > tableSets;
    // this is for holding all the cross-sections which are looked up by the
    // signal handlers.
  };




  // inline functions:


  inline std::string const*
  squareGrid::getName()
  const
  {

    return &gridName;

  }

  inline double
  squareGrid::getLowestX()
  const
  {

    return lowestXCoordinate;

  }

  inline double
  squareGrid::getHighestX()
  const
  {

    return highestXCoordinate;

  }

  inline double
  squareGrid::getLowestY()
  const
  {

    return lowestYCoordinate;

  }

  inline double
  squareGrid::getHighestY()
  const
  {

    return highestYCoordinate;

  }



  inline signedParticleShortcutPair const*
  crossSectionTable::getPair()
  const
  {

    return scoloredPair;

  }



  inline int
  crossSectionTableSet::getBeamEnergy()
  const
  {

    return beamEnergy;

  }



  inline crossSectionTable*
  crossSectionHandler::getTable( int const beamEnergyInTev,
                  signedParticleShortcutPair const* const requestedChannel )
  {

    return getTableSet( beamEnergyInTev )->getTable( requestedChannel );

  }

}  // end of LHC_FASER namespace.

#endif /* LHC_FASER_CROSS_SECTION_STUFF_HPP_ */
