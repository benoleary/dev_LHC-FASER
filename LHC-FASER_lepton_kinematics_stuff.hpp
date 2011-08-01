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
 * - all of the leptonAcceptanceParameterSet class
 * - all of the leptonAcceptanceGrid class, in particular
 *   -- leptonAcceptanceGrid::getValue()
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
                          CppSLHA::particle_property_set const* const scolored,
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

  protected:
    static double const defaultBinSize = 2.0;
    static double const defaultTransverseMomentumCut = 10.0;

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



  /* this class does stuff!
   */
  class leptonAcceptanceTableSet
  {

  public:

    leptonAcceptanceTableSet( int const given_LHC_energy_in_TeV,
                          std::string const* const given_grid_directory,
                          std::string const* const given_jet_subdirectory,
                          lepton_acceptance_table* const given_lepton_table,
                          input_handler const* const given_shortcuts )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~leptonAcceptanceTableSet()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    kinematics_table*
    get_table( int const given_acceptance_column )
    // this returns the kinematics_table for the requested column.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    bool
    is_requested( int const given_LHC_energy_in_TeV,
                  std::string const* const given_jet_subdirectory )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:

    input_handler const* const shortcut;

    int const LHC_energy_in_TeV;
    std::string const grid_directory;
    std::string jet_grid_subdirectory;

    lepton_acceptance_table* lepton_table;
    jet_acceptance_table* jet_table;
    std::vector< kinematics_table* > kinematics_tables;

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
}  // end of LHC_FASER namespace.

#endif /* LHC_FASER_LEPTON_KINEMATICS_STUFF_HPP_ */
