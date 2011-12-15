/*
 * LHC-FASER_lepton_kinematics_stuff.hpp
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

#ifndef LHC_FASER_LEPTON_KINEMATICS_STUFF_HPP_
#define LHC_FASER_LEPTON_KINEMATICS_STUFF_HPP_

#include "LHC-FASER_base_kinematics_stuff.hpp"
#include "LHC-FASER_derived_lepton_distributions.hpp"


/* this file contains the specific code relevant to obtaining the acceptance
 * for lepton transverse momentum and pseudorapidity cuts.
 *
 * if the format of the grid files changes, a lot of the code here has to
 * change too. the following are format-dependent:
 * - the leptonAcceptanceParameterSet class
 * - the leptonAcceptanceGrid class, in particular
 *   -- leptonAcceptanceGrid::getValue()
 * - the squarkBasedLeptonAcceptanceTable & gluinoBasedLeptonAcceptanceTable
 *   classes, in particular
 *   -- squarkBasedLeptonAcceptanceTable::getValue
 *   -- gluinoBasedLeptonAcceptanceTable::getValue
 *
 * other thoughts: actually, for lepton acceptances, it's probably just as
 * accurate to just take the grid points nearest the requested point rather
 * than bothering with a lot of interpolation. I might get around to
 * implementing that (without removing the interpolating versions though!).
 */

namespace LHC_FASER
{
  // this is a class derived from acceptanceCutSet, restricted to only look
  // at primaryLeptonCut.
  class singleLeptonCut : public acceptanceCutSet
  {
  public:
    singleLeptonCut();
    virtual
    ~singleLeptonCut();

    virtual bool
    isSameAcceptanceCutSet( acceptanceCutSet const* const comparisonCuts )
    const;

  //protected:
    // nothing
  };


  /* this is a derived class to interpret an acceptanceGrid as acceptances for
   * lepton distributions. this class, as well as acceptanceGrid, has to be
   * changed if the format of the acceptance grids changes.
   */
  class leptonAcceptanceGrid : public acceptanceGrid
  {
  public:
    leptonAcceptanceGrid( std::string const* const gridFileLocation,
                          inputHandler const* const inputShortcut );
    virtual
    ~leptonAcceptanceGrid();

    virtual void
    interpolateAcceptances( double squarkMass,
                            double gluinoMass,
                            double const electroweakinoMass,
                            double* const effectiveSquarkMass,
                            double* const pseudorapidityAcceptance,
                            std::vector< double >* const energyAcceptances )
    = 0;
    // this interpolates the grid to obtain values for the acceptances based on
    // the given masses.

  //protected:
    // nothing.
  };


  /* this is a derived class to interpret an acceptanceGrid as acceptances for
   * lepton distributions. this class, as well as acceptanceGrid, has to be
   * changed if the format of the acceptance grids changes.
   */
  class leptonAcceptanceFromSquarkGrid : public leptonAcceptanceGrid
  {
  public:
    leptonAcceptanceFromSquarkGrid( std::string const* const gridFileLocation,
                                    inputHandler const* const inputShortcut );
    virtual
    ~leptonAcceptanceFromSquarkGrid();

    virtual void
    interpolateAcceptances( double squarkMass,
                            double gluinoMass,
                            double const electroweakinoMass,
                            double* const effectiveSquarkMass,
                            double* const pseudorapidityAcceptance,
                            std::vector< double >* const energyAcceptances );
    // this interpolates the grid to obtain values for the acceptances based on
    // the given masses.

  protected:
    std::vector< double >* lowerLeftVector;
    std::vector< double >* lowerRightVector;
    std::vector< double >* upperRightVector;
    std::vector< double >* upperLeftVector;
  };


  /* this is a derived class to interpret an acceptanceGrid as acceptances for
   * lepton distributions. this class, as well as acceptanceGrid, has to be
   * changed if the format of the acceptance grids changes.
   */
  class leptonAcceptanceFromGluinoGrid : public leptonAcceptanceGrid
  {
  public:
    leptonAcceptanceFromGluinoGrid( std::string const* const gridFileLocation,
                                    inputHandler const* const inputShortcut );
    virtual
    ~leptonAcceptanceFromGluinoGrid();

    virtual void
    interpolateAcceptances( double squarkMass,
                            double gluinoMass,
                            double const electroweakinoMass,
                            double* const effectiveSquarkMass,
                            double* const pseudorapidityAcceptance,
                            std::vector< double >* const energyAcceptances );
    // this interpolates the grid to obtain values for the acceptances based on
    // the given masses.

  protected:
    unsigned int lowerElectroweakinoMassIndex;
    unsigned int upperElectroweakinoMassIndex;
    double electroweakinoMassFraction;
    bool shouldInterpolateOnElectroweakino;
    double otherElectroweakinoValue;
    std::vector< double >* foreLowerLeftVector;
    std::vector< double >* foreLowerRightVector;
    std::vector< double >* foreUpperRightVector;
    std::vector< double >* foreUpperLeftVector;
    std::vector< double >* rearLowerLeftVector;
    std::vector< double >* rearLowerRightVector;
    std::vector< double >* rearUpperRightVector;
    std::vector< double >* rearUpperLeftVector;
    // left-right is lighter-heavier squark mass,
    // lower-upper is lighter-heavier gluino mass,
    // fore-rear is lighter-heavier electroweakino mass.
  };


  /* this stores the binned lepton transverse momentum acceptances &
   * returns interpolated acceptances (given a transverse momentum cut to scale
   * to). it also stores the effective squark mass & the pseudorapidity cut
   * acceptance.
   */
  class leptonAcceptanceParameterSet : public getsReadiedForNewPoint,
                                       public effectiveSquarkMassHolder
  {
  public:
    static double const defaultBinSize;
    // by default, the grids provide acceptances at 2 GeV spacings for the
    // default cut on the lepton transverse momentum in the lab rest frame.
    static double const defaultTransverseMomentumCut;
    // by default, the grids provide acceptances assuming a 10 GeV cut on the
    // lepton transverse momentum in the lab rest frame.

    leptonAcceptanceParameterSet( inputHandler const* const inputShortcut,
                                  leptonAcceptanceGrid* const acceptanceGrid,
                                  particlePointer const coloredSparticle,
                                  particlePointer const ewinoPointer,
                                  double const binSize = defaultBinSize,
           double const transverseMomentumCut = defaultTransverseMomentumCut );
    ~leptonAcceptanceParameterSet();

    double
    getEffectiveSquarkMass();
    double
    acceptanceAt( double const givenEnergy,
                  double const givenCut );
    /* this checks to see if the acceptances need updating, then returns
     * calculateAcceptanceAt( givenEnergy,
     *                        givenCut ), which interpolates the acceptances in
     * acceptanceBins to the requested value, or returns
     * pseudorapidityAcceptance if it's lower, scaled to the given value for
     * the transverse momentum cut.
     */

    particlePointer
    getEwino()
    const;

  protected:
    inputHandler const* const inputShortcut;
    leptonAcceptanceGrid* const acceptanceGrid;
    particlePointer const coloredSparticle;
    particlePointer const ewinoPointer;
    double const binSize;
    double const transverseMomentumCut;
    double effectiveSquarkMass;
    std::vector< double > acceptanceBins;
    double pseudorapidityAcceptance;
    // these are used by acceptanceAt():
    double returnValue;
    double binFraction;
    int lowerBin;
    // these are used by calculate():
    int acceptanceCounter;
    double currentAcceptance;

    void
    resetValues();
    // this interpolates values from acceptanceTable to set up acceptanceBins
    // for the given colored sparticle.

    double
    calculateAcceptanceAt( double const givenEnergy,
                           double const givenCut );
    //const
    /* this interpolates the values in acceptanceBins to the requested value,
     * or returns pseudorapidityAcceptance if it's lower, scaled to the given
     * value for the transverse momentum cut.
     */
  };


  // this class holds all the leptonAcceptanceParameterSets for a given
  // coloredSparticle at a given beam energy.
  class leptonAcceptancesForOneScolored
  {
  public:
    leptonAcceptancesForOneScolored( inputHandler const* const inputShortcut,
                                     particlePointer const coloredSparticle,
                                    leptonAcceptanceGrid* const acceptanceGrid,
                                     double const binSize,
                                     double const transverseMomentumCut );
    ~leptonAcceptancesForOneScolored();

    leptonAcceptanceParameterSet*
    getParameterSet( particlePointer const ewinoPointer );
    // this returns the leptonAcceptanceParameterSet for the requested
    // electroweakino.
    particlePointer
    getScolored()
    const;

  protected:
    inputHandler const* const inputShortcut;
    particlePointer const coloredSparticle;
    leptonAcceptanceGrid* const acceptanceGrid;
    std::vector< leptonAcceptanceParameterSet* > parameterSets;
    double const binSize;
    // by default, the grids provide acceptances at 2 GeV spacings for the
    // default cut on the lepton transverse momentum in the lab rest frame.
    double const transverseMomentumCut;
    // by default, the grids provide acceptances assuming a 10 GeV cut on the
    // lepton transverse momentum in the lab rest frame.
  };


  // this class holds all the leptonAcceptancesForOneScoloreds for a given
  // beam energy.
  class leptonAcceptancesForOneBeamEnergy
  {
  public:
    leptonAcceptancesForOneBeamEnergy( inputHandler const* const inputShortcut,
                                       int const beamEnergy,
                                  std::string const* const gridFileSetLocation,
                                       double const binSize,
                                       double const transverseMomentumCut );
    ~leptonAcceptancesForOneBeamEnergy();

    leptonAcceptancesForOneScolored*
    getParameterSets( particlePointer const coloredSparticle );
    // this returns the leptonAcceptancesForOneScolored for the requested
    // colored sparticle.
    bool
    isRequested( int const beamEnergy,
                 double const binSize,
                 double const transverseMomentumCut )
    const;

  protected:
    inputHandler const* const inputShortcut;
    int const beamEnergy;
    std::string gridFileLocation;
    leptonAcceptanceGrid* acceptanceFromSquarkGrid;
    leptonAcceptanceGrid* acceptanceFromGluinoGrid;
    std::vector< leptonAcceptancesForOneScolored* > squarkAcceptanceSets;
    leptonAcceptancesForOneScolored* gluinoAcceptanceSet;
    double const binSize;
    // by default, the grids provide acceptances at 2 GeV spacings for the
    // default cut on the lepton transverse momentum in the lab rest frame.
    double const transverseMomentumCut;
    // by default, the grids provide acceptances assuming a 10 GeV cut on the
    // lepton transverse momentum in the lab rest frame.
  };


  // this class holds all the leptonAcceptancesForOneBeamEnergy & passes out
  // pointers to them & to requested leptonAcceptanceParameterSets.
  class leptonAcceptanceHandler
  {
  public:
    leptonAcceptanceHandler( inputHandler const* const inputShortcut,
                             std::string const* const gridFileSetLocation );
    ~leptonAcceptanceHandler();

    leptonAcceptanceParameterSet*
    getLeptonAcceptanceParameterSet( int const beamEnergy,
                                     particlePointer const coloredSparticle,
                                     particlePointer const ewinoPointer,
                                     double const binSize,
                                     double const transverseMomentumCut );
    /* this looks to see if there is an existing leptonAcceptanceParameterSet
     * with the requested acceptances, & if not, makes 1, & returns the
     * pointer.
     */
    leptonAcceptanceParameterSet*
    getLeptonAcceptanceParameterSet( int const beamEnergy,
                                     particlePointer const coloredSparticle,
                                     particlePointer const ewinoPointer );
    // this calls getLeptonAcceptanceParameterSet with the default binSize &
    // transverseMomentumCut.
    leptonAcceptancesForOneBeamEnergy*
    getLeptonAcceptancesForOneBeamEnergy( int const beamEnergy );
    /* this looks to see if there is an existing
     * leptonAcceptancesForOneBeamEnergy with the requested acceptances, & if
     * not, makes 1, & returns the pointer.
     */
    leptonAcceptancesForOneBeamEnergy*
    getLeptonAcceptancesForOneBeamEnergy( int const beamEnergy,
                                          double const binSize,
                                          double const transverseMomentumCut );
    /* this looks to see if there is an existing
     * leptonAcceptancesForOneBeamEnergy with the requested acceptances, & if
     * not, makes 1, & returns the pointer.
     */

  protected:
    inputHandler const* const inputShortcut;
    std::string const gridFileSetLocation;
    std::vector< leptonAcceptancesForOneBeamEnergy* > acceptanceTables;
  };





  // inline functions:


  inline bool
  singleLeptonCut::isSameAcceptanceCutSet(
                                 acceptanceCutSet const* const comparisonCuts )
  const
  {
    if( ( comparisonCuts->getBeamEnergy() == beamEnergy )
        &&
        ( comparisonCuts->getPrimaryLeptonCut() == primaryLeptonCut )  )
    {
      return true;
    }
    else
    {
      return false;
    }
  }



  inline double
  leptonAcceptanceParameterSet::getEffectiveSquarkMass()
  {
    if( needsToPrepareForThisPoint() )
    {
      resetValues();
    }
    return effectiveSquarkMass;
  }

  inline double
  leptonAcceptanceParameterSet::acceptanceAt( double const givenEnergy,
                                              double const givenCut )
  /* this checks to see if the acceptancesPerCutSet need updating, then returns
   * calculateAcceptanceAt( givenEnergy,
   *                        givenCut ), which interpolates the acceptances in
   * acceptanceBins to the requested value, or returns
   * pseudorapidityAcceptance if it's lower, scaled to the given value for
   * the transverse momentum cut.
   */
  {
    if( needsToPrepareForThisPoint() )
    {
      resetValues();
    }
    return calculateAcceptanceAt( givenEnergy,
                                  givenCut );
  }

  inline particlePointer
  leptonAcceptanceParameterSet::getEwino()
  const
  {
    return ewinoPointer;
  }

  inline void
  leptonAcceptanceParameterSet::resetValues()
  // this interpolates acceptanceValues from acceptanceTable to set up
  // acceptanceBins for the given colored sparticle.
  {
    acceptanceGrid->interpolateAcceptances(
                                         coloredSparticle->get_absolute_mass(),
                                            inputShortcut->getGluinoMass(),
                                            ewinoPointer->get_absolute_mass(),
                                            &effectiveSquarkMass,
                                            &pseudorapidityAcceptance,
                                            &acceptanceBins );
  }



  inline particlePointer
  leptonAcceptancesForOneScolored::getScolored()
  const
  {
    return coloredSparticle;
  }



  inline bool
  leptonAcceptancesForOneBeamEnergy::isRequested( int const beamEnergy,
                                                  double const binSize,
                                           double const transverseMomentumCut )
  const
  {
    if( ( beamEnergy == this->beamEnergy )
        &&
        ( binSize == this->binSize )
        &&
        ( transverseMomentumCut == this->transverseMomentumCut ) )
    {
      return true;
    }
    else
    {
      return false;
    }
  }



  inline leptonAcceptanceParameterSet*
  leptonAcceptanceHandler::getLeptonAcceptanceParameterSet(
                                                          int const beamEnergy,
                                        particlePointer const coloredSparticle,
                                            particlePointer const ewinoPointer,
                                                          double const binSize,
                                          double const transverseMomentumCut )
  // this looks to see if there is an existing leptonAcceptanceParameterSet
  // with the requested acceptances, & if not, makes 1, & returns the pointer.
  {
    return getLeptonAcceptancesForOneBeamEnergy( beamEnergy,
                                                 binSize,
                                                 transverseMomentumCut
                                                           )->getParameterSets(
                                                               coloredSparticle
                                                            )->getParameterSet(
                                                                ewinoPointer );
  }

  inline leptonAcceptanceParameterSet*
  leptonAcceptanceHandler::getLeptonAcceptanceParameterSet(
                                                          int const beamEnergy,
                                        particlePointer const coloredSparticle,
                                           particlePointer const ewinoPointer )
  // this calls getLeptonAcceptanceParameterSet with the default binSize &
  // transverseMomentumCut.
  {
    return getLeptonAcceptanceParameterSet( beamEnergy,
                                            coloredSparticle,
                                            ewinoPointer,
                                  leptonAcceptanceParameterSet::defaultBinSize,
                  leptonAcceptanceParameterSet::defaultTransverseMomentumCut );
  }

  inline leptonAcceptancesForOneBeamEnergy*
  leptonAcceptanceHandler::getLeptonAcceptancesForOneBeamEnergy(
                                                         int const beamEnergy )
  /* this looks to see if there is an existing
   * leptonAcceptancesForOneBeamEnergy with the requested values, & if not,
   * makes 1, & returns the pointer.
   */
  {
    return
    leptonAcceptanceHandler::getLeptonAcceptancesForOneBeamEnergy( beamEnergy,
                                  leptonAcceptanceParameterSet::defaultBinSize,
                  leptonAcceptanceParameterSet::defaultTransverseMomentumCut );
  }

}  // end of LHC_FASER namespace.

#endif /* LHC_FASER_LEPTON_KINEMATICS_STUFF_HPP_ */
