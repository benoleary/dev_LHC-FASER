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
 *   classes, inparticular
 *   -- squarkBasedLeptonAcceptanceTable::getValue
 *   -- gluinoBasedLeptonAcceptanceTable::getValue
 */

namespace LHC_FASER
{
  // this is a class derived from acceptanceCutSet, restricted to only look
  // at primaryLeptonCut.
  class singleLeptonCut : public acceptanceCutSet
  {
  public:
    singleLeptonCut()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~singleLeptonCut()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual bool
    isSameAcceptanceCutSet( acceptanceCutSet const* const comparisonCuts )
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  //protected:
    // nothing
  };


  /* this class holds acceptanceGrid instances for looking up lepton
   * acceptances & returns the interpolated values by looking up the masses of
   * its sparticles. this class, as well as acceptanceGrid, has to be changed
   * if the format of the acceptance grids changes.
   */
  class leptonAcceptanceGrid
  {
  public:
    leptonAcceptanceGrid( std::string const* const gridFileLocation,
                          input_handler const* const shortcut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~leptonAcceptanceGrid()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    getSquarkValue( CppSLHA::particle_property_set const* const squark,
                    CppSLHA::particle_property_set const* const ewino,
                    int const requestedColumn,
                    bool const lookingForPseudorapidity )
    const
    /* this interpolates the requested column based on the squark, gluino, &
     * electroweakino masses. it fudges some cases that were not properly done
     * in the single-quark-flavor approximation.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getGluinoValue( CppSLHA::particle_property_set const* const ewino,
                    int const requestedColumn,
                    bool const lookingForEffectiveSquarkMass,
                    bool const lookingForPseudorapidity )
    const
    /* this interpolates the requested column based on the squark, gluino, &
     * electroweakino masses. it fudges some cases that were not properly done
     * in the single-quark-flavor approximation.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    acceptanceGrid* gluinoTable;
    acceptanceGrid* squarkTable;
    input_handler const* const shortcut;
  };


  /* this class accesses a leptonAcceptanceGrid for a specific colored
   * sparticle & beam energy. it is an abstract base class, for separate
   * derived versions for the gluino & squarks.
   */
  class leptonAcceptanceTable
  {
  public:
    leptonAcceptanceTable( leptonAcceptanceGrid const* const lookupGrid,
                          CppSLHA::particle_property_set const* const scolored,
                           input_handler const* const shortcut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    virtual
    ~leptonAcceptanceTable()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual double
    getValue( CppSLHA::particle_property_set const* const ewino,
              int const requestedColumn )
    const
    // this gets the appropriate interpolated value from lookupGrid.
    = 0;

  protected:
    leptonAcceptanceGrid const* const lookupGrid;
    CppSLHA::particle_property_set const* const scolored;
    input_handler const* const shortcut;
  };

  // this class is the class derived from leptonAcceptanceTable specifically
  // for squarks.
  class squarkBasedLeptonAcceptanceTable : public leptonAcceptanceTable
  {
  public:
    squarkBasedLeptonAcceptanceTable(
                                  leptonAcceptanceGrid const* const lookupGrid,
                            CppSLHA::particle_property_set const* const squark,
                                      input_handler const* const shortcut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    virtual
    ~squarkBasedLeptonAcceptanceTable()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    getValue( CppSLHA::particle_property_set const* const ewino,
              int const requestedColumn )
    const
    // this gets the appropriate interpolated value from lookupGrid.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  //protected:
    // nothing.
  };

  // this class is the class derived from leptonAcceptanceTable specifically
  // for the gluino.
  class gluinoBasedLeptonAcceptanceTable : public leptonAcceptanceTable
  {
  public:
    gluinoBasedLeptonAcceptanceTable(
                                  leptonAcceptanceGrid const* const lookupGrid,
                                      input_handler const* const shortcut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    virtual
    ~gluinoBasedLeptonAcceptanceTable()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    getValue( CppSLHA::particle_property_set const* const ewino,
              int const requestedColumn )
    const
    // this gets the appropriate interpolated value from lookupGrid.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  //protected:
    // nothing.
  };


  /* this stores the binned lepton transverse momentum acceptance values &
   * returns interpolated values (given a transverse momentum cut to scale to).
   * it also stores the effective squark mass & the pseudorapidity cut
   * acceptance.
   */
  class leptonAcceptanceParameterSet : public readied_for_new_point,
                                              effectiveSquarkMassHolder
  {
  public:
    static double const defaultBinSize = 2.0;
    static double const defaultTransverseMomentumCut = 10.0;

    leptonAcceptanceParameterSet( input_handler* const shortcut,
                                  leptonAcceptanceTable const* acceptanceTable,
                          CppSLHA::particle_property_set const* const scolored,
                            CppSLHA::particle_property_set const* const ewino )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    leptonAcceptanceParameterSet( input_handler* const shortcut,
                                  leptonAcceptanceTable const* acceptanceTable,
                          CppSLHA::particle_property_set const* const scolored,
                             CppSLHA::particle_property_set const* const ewino,
                                  double const binSize,
                                  double const transverseMomentumCut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~leptonAcceptanceParameterSet()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    getEffectiveSquarkMass()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    acceptanceAt( double const givenEnergy,
                  double const givenCut )
    /* this checks to see if the acceptances need updating, then returns
     * calculateAcceptanceAt( givenEnergy,
     *                        givenCut ), which interpolates the values in
     * acceptanceBins to the requested value, or returns
     * pseudorapidityAcceptance if it's lower, scaled to the given value for
     * the transverse momentum cut.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    CppSLHA::particle_property_set const*
    getEwino()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    input_handler* const shortcut;
    leptonAcceptanceTable const* acceptanceTable;
    CppSLHA::particle_property_set const* const scolored;
    CppSLHA::particle_property_set const* const ewino;
    double const binSize;
    // by default, the grids provide acceptances at 2 GeV spacings for the
    // default cut on the lepton transverse momentum in the lab rest frame.
    double const transverseMomentumCut;
    // by default, the grids provide acceptances assuming a 10 GeV cut on the
    // lepton transverse momentum in the lab rest frame.
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
    resetValues()
    // this interpolates values from acceptanceTable to set up acceptanceBins
    // for the given colored sparticle.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    calculateAcceptanceAt( double const givenEnergy,
                           double const givenCut )
    /* this interpolates the values in acceptanceBins to the requested value,
     * or returns pseudorapidityAcceptance if it's lower, scaled to the given
     * value for the transverse momentum cut.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };


  // this class holds all the leptonAcceptanceParameterSets for a given
  // scolored at a given beam energy.
  class leptonAcceptancesForOneScolored
  {
  public:
    leptonAcceptancesForOneScolored( input_handler const* const shortcut,
                          CppSLHA::particle_property_set const* const scolored,
                                 leptonAcceptanceGrid const* acceptanceGrid,
                                     double const binSize,
                                     double const transverseMomentumCut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~leptonAcceptancesForOneScolored()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    leptonAcceptanceParameterSet*
    getParameterSet( CppSLHA::particle_property_set const* const ewino )
    // this returns the leptonAcceptanceParameterSet for the requested
    // electroweakino.
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    CppSLHA::particle_property_set const*
    getScolored()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    input_handler const* const shortcut;
    CppSLHA::particle_property_set const* const scolored;
    leptonAcceptanceTable const* acceptanceTable;
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
    leptonAcceptancesForOneBeamEnergy( input_handler const* const shortcut,
                                       int const beamEnergy,
                                  std::string const* const gridFileSetLocation,
                                       double const binSize,
                                       double const transverseMomentumCut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~leptonAcceptancesForOneBeamEnergy()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    leptonAcceptancesForOneScolored*
    getParameterSets( CppSLHA::particle_property_set const* const scolored )
    // this returns the leptonAcceptancesForOneScolored for the requested
    // colored sparticle.
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    bool
    isRequested( int const beamEnergy,
                 double const binSize,
                 double const transverseMomentumCut )
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    input_handler const* const shortcut;
    int const beamEnergy;
    leptonAcceptanceGrid acceptanceGrid;
    std::vector< leptonAcceptancesForOneScolored* > acceptanceSets;
    double const binSize;
    // by default, the grids provide acceptances at 2 GeV spacings for the
    // default cut on the lepton transverse momentum in the lab rest frame.
    double const transverseMomentumCut;
    // by default, the grids provide acceptances assuming a 10 GeV cut on the
    // lepton transverse momentum in the lab rest frame.
  };


  // this class holds all the leptonAcceptancesForOneBeamEnergy & passes out
  // pointers to requested leptonAcceptanceParameterSets.
  class leptonAcceptanceHandler
  {
  public:
    leptonAcceptanceHandler( input_handler const* const shortcut,
                             std::string const* const gridFileSetLocation )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~leptonAcceptanceHandler()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    leptonAcceptanceParameterSet*
    getLeptonAcceptanceParameterSet( int const beamEnergy,
                          CppSLHA::particle_property_set const* const scolored,
                             CppSLHA::particle_property_set const* const ewino,
                                     double const binSize,
                                     double const transverseMomentumCut )
    // this looks to see if there is an existing leptonAcceptanceParameterSet
    // with the requested values, & if not, makes 1, & returns the pointer.
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    leptonAcceptanceParameterSet*
    getLeptonAcceptanceParameterSet( int const beamEnergy,
                          CppSLHA::particle_property_set const* const scolored,
                            CppSLHA::particle_property_set const* const ewino )
    // this calls getLeptonAcceptanceParameterSet with the default binSize &
    // transverseMomentumCut.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    input_handler const* const shortcut;
    std::string const gridFileSetLocation;
    std::vector< leptonAcceptancesForOneBeamEnergy* > acceptanceTables;

    leptonAcceptancesForOneBeamEnergy*
    getleptonAcceptancesForOneBeamEnergy( int const beamEnergy,
                                          double const binSize,
                                          double const transverseMomentumCut )
    /* this looks to see if there is an existing
     * leptonAcceptancesForOneBeamEnergy with the requested values, & if not,
     * makes 1, & returns the pointer.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };



  // inline functions:

  inline bool
  singleLeptonCut::isSameAcceptanceCutSet(
                              acceptanceCutSet const* const comparison_cuts )
  const
  {

    if( ( comparison_cuts->getBeamEnergy() == beamEnergy )
        &&
        ( comparison_cuts->getPrimaryLeptonCut() == primaryLeptonCut )  )
    {

      return true;

    }
    else
    {

      return false;

    }

  }



  inline double
  squarkBasedLeptonAcceptanceTable::getValue(
                             CppSLHA::particle_property_set const* const ewino,
                                              int const requestedColumn )
  const
  // this gets the appropriate interpolated value from lookupGrid.
  {
    if( 0 == requestedColumn )
    {
      return scolored->get_absolute_mass();
    }
    else if( 1 == requestedColumn )
    {
      return lookupGrid->getSquarkValue( scolored,
                                         ewino,
                                         requestedColumn,
                                         true );
    }
    else
    {
      return lookupGrid->getSquarkValue( scolored,
                                         ewino,
                                         requestedColumn,
                                         false );
    }
  }



  inline double
  gluinoBasedLeptonAcceptanceTable::getValue(
                             CppSLHA::particle_property_set const* const ewino,
                                              int const requestedColumn,
                                      bool const lookingForEffectiveSquarkMass,
                                          bool const lookingForPseudorapidity )
  const
  // this gets the appropriate interpolated value from lookupGrid.
  {
    if( 0 == requestedColumn )
    {
      return lookupGrid->getGluinoValue( ewino,
                                         requestedColumn,
                                         true,
                                         false );
    }
    else if( 1 == requestedColumn )
    {
      return lookupGrid->getGluinoValue( ewino,
                                         requestedColumn,
                                         false,
                                         true );
    }
    else
    {
      return lookupGrid->getGluinoValue( ewino,
                                         requestedColumn,
                                         false,
                                         false );
    }
  }



  inline double
  leptonAcceptanceParameterSet::getEffectiveSquarkMass()
  {
    // debugging:
    /**std::cout
    << std::endl
    << "debugging: leptonAcceptanceParameterSet::getEffectiveSquarkMass()"
    << " called.";
    std::cout << std::endl;**/

    if( needs_to_prepare_for_this_point() )
    {
      resetValues();
    }
    return effectiveSquarkMass;
  }

  inline double
  leptonAcceptanceParameterSet::acceptanceAt( double const givenEnergy,
                                     double const givenCut )
  /* this checks to see if the acceptances need updating, then returns
   * calculateAcceptanceAt( givenEnergy,
   *                        givenCut ), which interpolates the values in
   * acceptance_bins to the requested value, or returns
   * pseudorapidityAcceptance if it's lower, scaled to the given value for
   * the transverse momentum cut.
   */
  {
    // debugging:
    /**std::cout
    << std::endl
    << "debugging: leptonAcceptanceParameterSet::acceptanceAt( " << givenEnergy
    << ", " << givenCut << " ) called.";
    std::cout << std::endl;**/

    if( needs_to_prepare_for_this_point() )
    {
      resetValues();
    }
    return calculateAcceptanceAt( givenEnergy,
                                  givenCut );
  }

  inline CppSLHA::particle_property_set const*
  leptonAcceptanceParameterSet::getEwino()
  const
  {
    return ewino;
  }



  inline CppSLHA::particle_property_set const*
  leptonAcceptancesForOneScolored::getScolored()
  const
  {
    return scolored;
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
                          CppSLHA::particle_property_set const* const scolored,
                             CppSLHA::particle_property_set const* const ewino,
                                                          double const binSize,
                                          double const transverseMomentumCut )
  // this looks to see if there is an existing leptonAcceptanceParameterSet
  // with the requested values, & if not, makes 1, & returns the pointer.
  {
    return getleptonAcceptancesForOneBeamEnergy( beamEnergy,
                                                 binSize,
                                                 transverseMomentumCut
                                                  )->getParameterSets( scolored
                                                   )->getParameterSet( ewino );
  }

  inline leptonAcceptanceParameterSet*
  leptonAcceptanceHandler::getLeptonAcceptanceParameterSet(
                                                          int const beamEnergy,
                           CppSLHA::particle_property_set const* const scolored,
                          CppSLHA::particle_property_set const* const ewino )
  // this calls getLeptonAcceptanceParameterSet with the default binSize &
  // transverseMomentumCut.
  {
    return getLeptonAcceptanceParameterSet( beamEnergy,
                                            scolored,
                                            ewino,
                                  leptonAcceptanceParameterSet::defaultBinSize,
                  leptonAcceptanceParameterSet::defaultTransverseMomentumCut );
  }
}  // end of LHC_FASER namespace.

#endif /* LHC_FASER_LEPTON_KINEMATICS_STUFF_HPP_ */
