/*
 * LHC-FASER_kinematics_stuff.hpp
 *
 *  Created on: 15 Dec 2010
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

#ifndef LHC_FASER_JET_KINEMATICS_STUFF_HPP_
#define LHC_FASER_JET_KINEMATICS_STUFF_HPP_

#include "LHC-FASER_input_handling_stuff.hpp"
#include "LHC-FASER_full_cascade_stuff.hpp"


/* this file contains all the code relevant to obtaining the acceptance for
 * jet transverse momentum plus missing transverse momentum cuts.
 *
 * if the format of the grid files changes, a lot of the code here has to
 * change too. the following are format-dependent:
 * - all of the leptonAcceptanceParameterSet class.
 * - all of the acceptanceGrid class.
 * - all of the leptonAcceptanceGrid class, in particular
 *   -- leptonAcceptanceGrid::getValue()
 */

namespace LHC_FASER
{
  // this is an abstract base class for a set of functors for deciding what
  // squark mass to use for looking up the grid.
  class squarkMassForGridDecider
  {
  public:
    squarkMassForGridDecider( inputHandler const* const inputShortcut );
    virtual
    ~squarkMassForGridDecider();

    virtual double
    operator()( particlePointer const firstScolored,
                particlePointer const secondScolored )
    = 0;

  protected:
    inputHandler const* const inputShortcut;
  };


  /* this is a derived class returning the average mass of the non-stop,
   * non-sbottom squarks, or the gluino mass + 1.0 if it would be less than the
   * gluino mass.
   */
  class heavierThanGluinoSquarkMassForGrid : public squarkMassForGridDecider
  {
  public:
    heavierThanGluinoSquarkMassForGrid(
                                     inputHandler const* const inputShortcut );
    virtual
    ~heavierThanGluinoSquarkMassForGrid();

    virtual double
    operator()( particlePointer const firstScolored,
                particlePointer const secondScolored );

  //protected:
    // nothing
  };


  // this is a derived class returning the mass of the 1st provided squark.
  class firstMassForGrid : public squarkMassForGridDecider
  {
  public:
    firstMassForGrid( inputHandler const* const inputShortcut );
    virtual
    ~firstMassForGrid();

    virtual double
    operator()( particlePointer const firstScolored,
                particlePointer const secondScolored );

  //protected:
    // nothing
  };


  // this is a derived class returning the mass of the 2nd provided squark.
  class secondMassForGrid : public squarkMassForGridDecider
  {
  public:
    secondMassForGrid( inputHandler const* const inputShortcut );
    virtual
    ~secondMassForGrid();

    virtual double
    operator()( particlePointer const firstScolored,
                particlePointer const secondScolored );

  //protected:
    // nothing
  };


  // this is a derived class returning the average mass of the provided
  // squarks.
  class averageSquarkMassForGrid : public squarkMassForGridDecider
  {
  public:
    averageSquarkMassForGrid( inputHandler const* const inputShortcut );
    virtual
    ~averageSquarkMassForGrid();

    virtual double
    operator()( particlePointer const firstScolored,
                particlePointer const secondScolored );

  //protected:
    // nothing
  };




  /* this is a derived class to interpret an acceptanceGrid as acceptances for
   * lepton distributions. this class, as well as acceptanceGrid, has to be
   * changed if the format of the acceptance grids changes.
   */
  class jetAcceptanceGrid : public acceptanceGrid
  {
  public:
    jetAcceptanceGrid( std::string const* const gridFileLocation,
                       inputHandler const* const inputShortcut );
    ~jetAcceptanceGrid();

    double
    getAcceptance( double const squarkMass,
                   double const gluinoMass,
                   double firstElectroweakinoMass,
                   double secondElectroweakinoMass,
                   int acceptanceColumn );
    // this interpolates the grid to obtain values for the acceptances based on
    // the given masses.

  protected:
    static int const columnsNotInAcceptanceVector;
    // we want to skip the squark, gluino, & electroweakino masses, & we also
    // want to count from 1 rather than from 0.

    double lowerScoloredMass;
    double lighterElectroweakinoMassFraction;
    double heavierElectroweakinoMassFraction;
    bool shouldInterpolateOnElectroweakinos;
    double lighterLighterElectroweakinoPointValue;
    double lighterHeavierElectroweakinoPointValue;
    double heavierHeavierElectroweakinoPointValue;
    double heavierLighterElectroweakinoPointValue;
    double returnValue;

    double
    interpolateOnElectroweakinoAndAcceptanceIndices(
                                                 int const electroweakinoIndex,
                                                   int const acceptanceIndex );
    /* this checks to see if the appropriate entries in
     * lowerLeftVectorOfVectors, lowerRightVectorOfVectors,
     * upperRightVectorOfVectors, & upperLeftVectorOfVectors exist, & if so,
     * interpolates on the squark & gluino masses for these indices.
     */
  };


  typedef std::pair< jetAcceptanceGrid*,
                     squarkMassForGridDecider* > gridWithSquarkMass;
  // this class holds an acceptanceGrid with a string identifying the type of
  // jet+MET signal which the acceptancesPerCutSet are for.
  class jetAcceptanceTable
  {
  public:
    enum usedCascades
    {
      gx = 0,
      sx = 1,
      gsx = 2,
      sgx = 3,
      sizeOfEnum = 4 /* should be useful */
    };
    // this could be expanded, provided grids that account for multiple squark
    // flavors are provided.
    static usedCascades const lastEnumElement;
    static int const sizeOfUsedCascades;

    jetAcceptanceTable( std::string const* const gridFilesLocation,
                        std::string const* const jetCutName,
                        int const acceptanceColumn,
                        inputHandler const* const inputShortcut );
    ~jetAcceptanceTable();

    bool
    isRequested( std::string const* const jetCutName,
                 int const acceptanceColumn )
    const;
    double
    getAcceptance( signedParticleShortcutPair const* const initialPair,
                   fullCascade const* const firstCascade,
                   fullCascade const* const secondCascade );

  protected:
    std::string const jetCutName;
    int const acceptanceColumn;
    inputHandler const* const inputShortcut;
    jetAcceptanceGrid* gluinoGluinoGrid;
    jetAcceptanceGrid* squarkGluinoGrid;
    jetAcceptanceGrid* squarkAntisquarkGrid;
    jetAcceptanceGrid* squarkSquarkGrid;
    jetAcceptanceGrid* gridToUse;
    minimalAllocationVector< minimalAllocationVector< gridWithSquarkMass > >
    gridsMatrix;
    minimalAllocationVector< gridWithSquarkMass >*
    gridMatrixRow;
    gridWithSquarkMass* gridMatrixElement;
    heavierThanGluinoSquarkMassForGrid heavierThanGluinoSquarkMass;
    firstMassForGrid useFirstMass;
    secondMassForGrid useSecondMass;
    averageSquarkMassForGrid useAverageMass;
    usedCascades firstCascadeType;
    particlePointer firstCascadeSquark;
    // this is NULL if there is no squark in the cascade.
    usedCascades secondCascadeType;
    particlePointer secondCascadeSquark;
    // this is NULL if there is no squark in the cascade.

    int
    getIntForCascadeType( fullCascade const* const givenCascade,
                          particlePointer* const squarkFromCascade );
  };


  // this class holds all the jetAcceptanceTables for a given beam energy.
  class jetAcceptanceTablesForOneBeamEnergy
  {
  public:
    jetAcceptanceTablesForOneBeamEnergy(
                                       inputHandler const* const inputShortcut,
                                         int const beamEnergy,
                                std::string const* const gridFileSetLocation );
    ~jetAcceptanceTablesForOneBeamEnergy();

    jetAcceptanceTable*
    getTable( std::string const* const jetCutName,
              int const acceptanceColumn );
    // this returns the jetAcceptanceTable for the requested signal & column.
    int
    getBeamEnergy()
    const;

  protected:
    inputHandler const* const inputShortcut;
    int const beamEnergy;
    std::string gridFileSetLocation;
    std::vector< jetAcceptanceTable* > jetTables;
  };


  // this class holds all the jetAcceptanceTablesForOneBeamEnergys & passes out
  // pointers to requested jetAcceptanceTables.
  class jetPlusMetAcceptanceHandler
  {
  public:
    jetPlusMetAcceptanceHandler( inputHandler const* const inputShortcut,
                                std::string const* const gridFileSetLocation );
    ~jetPlusMetAcceptanceHandler();

    jetAcceptanceTable*
    getJetPlusMetAcceptanceTable( int const beamEnergy,
                                  std::string const* const jetCutName,
                                  int const acceptanceColumn );
    // this looks to see if there is an existing jetAcceptanceTable with the
    // requested acceptanceValues, & if not, makes 1, & returns the pointer.

  protected:
    inputHandler const* const inputShortcut;
    std::string const gridFileSetLocation;
    std::vector< jetAcceptanceTablesForOneBeamEnergy* > acceptanceTables;

    jetAcceptanceTablesForOneBeamEnergy*
    getJetAcceptanceTablesForOneBeamEnergy( int const beamEnergy );
    /* this looks to see if there is an existing
     * jetAcceptanceTablesForOneBeamEnergy with the requested beamEnergy, & if
     * not, makes 1, & returns the pointer.
     */
  };





  // inline functions:


  inline double
  heavierThanGluinoSquarkMassForGrid::operator()(
                                           particlePointer const firstScolored,
                                         particlePointer const secondScolored )
  {
    double returnValue( inputShortcut->getAverageSquarks4Mass() );
    if( inputShortcut->getGluinoMass() > returnValue )
    {
      returnValue = ( inputShortcut->getGluinoMass() + 1.0 );
    }
    return returnValue;
  }


  inline double
  firstMassForGrid::operator()( particlePointer const firstScolored,
                                particlePointer const secondScolored )
  {
    return firstScolored->get_absolute_mass();
  }


  inline double
  secondMassForGrid::operator()( particlePointer const firstScolored,
                                 particlePointer const secondScolored )
  {
    return secondScolored->get_absolute_mass();
  }


  inline double
  averageSquarkMassForGrid::operator()( particlePointer const firstScolored,
                                        particlePointer const secondScolored )
  {
    return ( 0.5 * ( firstScolored->get_absolute_mass()
                     + secondScolored->get_absolute_mass() ) );
  }



  inline bool
  jetAcceptanceTable::isRequested( std::string const* const jetCutName,
                                   int const acceptanceColumn )
  const
  {
    if( ( acceptanceColumn == this->acceptanceColumn )
        &&
        ( 0 == jetCutName->compare( this->jetCutName ) ) )
    {
      return true;
    }
    else
    {
      return false;
    }
  }



  inline int
  jetAcceptanceTablesForOneBeamEnergy::getBeamEnergy()
  const
  {
    return beamEnergy;
  }



  inline jetAcceptanceTable*
  jetPlusMetAcceptanceHandler::getJetPlusMetAcceptanceTable(
                                                          int const beamEnergy,
                                           std::string const* const jetCutName,
                                                   int const acceptanceColumn )
  // this looks to see if there is an existing jetAcceptanceTable with the
  // requested values, & if not, makes 1, & returns the pointer.
  {
    return
    getJetAcceptanceTablesForOneBeamEnergy( beamEnergy )->getTable( jetCutName,
                                                            acceptanceColumn );
  }

}  // end of LHC_FASER namespace.

#endif /* LHC_FASER_JET_KINEMATICS_STUFF_HPP_ */
