/*
 * LHC-FASER_base_kinematics_stuff.hpp
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

#ifndef LHC_FASER_BASE_KINEMATICS_STUFF_HPP_
#define LHC_FASER_BASE_KINEMATICS_STUFF_HPP_

#include "LHC-FASER_input_handling_stuff.hpp"


/* this file contains the common code relevant to obtaining the acceptance for
 * jet transverse momentum plus missing transverse momentum and lepton
 * transverse momentum and pseudorapidity cuts. the code specific to lepton or
 * jet acceptances are in LHC-FASER_lepton_kinematics_stuff.hpp &
 * LHC-FASER_jet_kinematics_stuff.hpp respectively (& the appropriate source
 * files).
 *
 * if the format of the grid files changes, a lot of the code here has to
 * change too. the following are format-dependent:
 * - the acceptanceGrid class, in particular
 *   -- acceptanceGrid::neutralinoIndependentValueAt
 */

namespace LHC_FASER
{
  // this stores the acceptance cut values which are associated with a
  // particular signal.
  class acceptanceCutSet
  {
  public:
    acceptanceCutSet()
    // this constructor sets the values to -1.0, which is the default "unset"
    // value.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    acceptanceCutSet( acceptanceCutSet* const copySource )
    // this constructor copies the values from a given acceptanceCutSet.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~acceptanceCutSet()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    static bool
    compareJetAndBothLeptonCuts( acceptanceCutSet* const firstPointer,
                                 acceptanceCutSet* const secondPointer )
    // this returns true if both acceptanceCutSets have the same jet cut & both
    // lepton cuts.
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    static bool
    compareJustBothLeptonCuts( acceptanceCutSet* const firstPointer,
                               acceptanceCutSet* const secondPointer )
    // this returns true if both acceptanceCutSets have the same
    // primaryLeptonCut & secondaryLeptonCut, regardless of the jet cut.
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    static bool
    compareJustPrimaryLeptonCut( acceptanceCutSet* const firstPointer,
                                 acceptanceCutSet* const secondPointer )
    // this returns true if both acceptanceCutSets have the same
    // primaryLeptonCut, regardless of the other values.
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    static bool
    justReturnTrue( acceptanceCutSet* const firstPointer,
                    acceptanceCutSet* const secondPointer )
    // this just returns true, for when the acceptance cut doesn't matter.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    becomeCopyOf( acceptanceCutSet* const copyPointer )
    // this copies the values from a given acceptanceCutSet.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    getBeamEnergy()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    // set_beam_energy is only defined in the signal_definitions derived class.
    double
    getPrimaryLeptonCut()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setPrimaryLeptonCut( double inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getSecondaryLeptonCut()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setSecondaryLeptonCut( double inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getLeptonCut()
    const
    // this just returns primary_lepton_cut. it is intended to be used when the
    // signal only has 1 lepton cut.
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setLeptonCut( double inputValue )
    // this sets both primary_lepton_cut & secondary_lepton_cut to input_value.
    // it is intended to be used when the signal only has 1 lepton cut.
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getJetCut()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setJetCut( double inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    std::list< int > const*
    getExcludedStandardModelProducts()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setExcludedStandardModelProducts( std::list< int > const* inputList )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    virtual bool
    isSameAcceptanceCutSet( acceptanceCutSet const* const comparisonCuts )
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    double beamEnergy;
    // the lepton transverse momentum cuts are easy enough to scale that we can
    // let the user specify them:
    double primaryLeptonCut;
    double secondaryLeptonCut;
    double jetCut;
    std::list< int > const* const excludedStandardModelProducts;
  };  // end of acceptanceCutSet class.


  /* this class reads in a file in the assumed format, stores it, & gives out
   * interpolated values. it was written with acceptances for leptons in mind,
   * to be adapted for jets plus missing transverse momentum as a special case,
   * with data files in the format
   * squarkMass gluinoMass lighterNeutralinoMass heavierNeutralinoMass
   * (continued) then either
   * 42 lepton acceptance values (effective squark mass, pseudorapidity cut
   * acceptance, then 40 bins for transverse momentum cut acceptances)
   * or
   * 7 acceptance values for different choices of which jets+MET combination to
   * use
   * then newline.
   */
  class acceptanceGrid
  {
  public:
    acceptanceGrid( std::string const* const gridFileLocation )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~acceptanceGrid()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    valueAt( double const squarkMass,
             double const gluinoMass,
             double const firstNeutralinoMass,
             double const secondNeutralinoMass,
             int const acceptanceElement,
             bool const heavyNeutralinoEdgeIsLighterScoloredMass,
             bool const heavyNeutralinoAreaIsConstant )
    const
    /* this finds the grid square which the given point is in, & then uses
     * LHC_FASER_global::square_bilinear_interpolation to get an interpolated
     * value, assuming that the heavy neutralino edge goes to 0.0 as the
     * heavier neutralino mass approaches the lighter scolored mass unless
     * heavy_neutralino_edge_is_lighter_scolored_mass is true, in which case
     * it interpolates to the lighter scolored mass, or unless
     * heavy_neutralino_area_is_constant is true, in which case it interpolates
     * to the same value as the grid points with the heaviest neutralino mass.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    valueAt( double const squarkMass,
             double const gluinoMass,
             double const degenerateNeutralinoMass,
             int const acceptanceElement,
             bool const heavyNeutralinoEdgeIsLighterScoloredMass,
             bool const heavyNeutralinoAreaIsConstant )
    const
    /* this finds the grid square which the given point is in, & then uses
     * LHC_FASER_global::square_bilinear_interpolation to get an interpolated
     * value, assuming that the heavy neutralino edge goes to 0.0 as the
     * heavier neutralino mass approaches the lighter scolored mass unless
     * heavyNeutralinoEdgeIsLighterScoloredMass is true, in which case it
     * interpolates to the lighter scolored mass, or unless
     * heavyNeutralinoAreaIsConstant is true, in which case it interpolates
     * to the same value as the grid points with the heaviest neutralino mass.
     * N.B.: this version is just to save a little calculation for the
     * lepton acceptance for a cascade because of the approximation that the
     * kinematics for the lepton acceptance of 1 cascade is independent of that
     * of the other cascade. the improvement is probably utterly negligible,
     * but I want to write this function anyway...
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    neutralinoIndependentValueAt( double const squarkMass,
                                  double const gluinoMass,
                                  int const acceptanceElement )
    const
    // this interpolates for acceptanceElement assuming that the neutralino
    // masses are on the light-light grid point.
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    int
    getNumberOfAcceptanceColumns()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getLowestSquarkMass()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getHighestSquarkMass()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getLowestGluinoMass()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getHighestGluinoMass()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    double scoloredMassStepSize;
    double lowestSquarkMass;
    double highestSquarkMass;
    double lowestGluinoMass;
    double highestGluinoMass;
    double lowNeutralinoMassRatio;
    double middleNeutralinoMassRatio;
    double highNeutralinoMassRatio;
    std::vector< std::vector< std::vector< std::vector< double >* >* >* >
    values;
    int acceptanceColumns;

    double
    vectorElementAt( double const squarkMass,
                     double const gluinoMass,
                     int const neutralinoElement,
                     int const acceptanceElement )
    const
    // this interpolates the values for the element of the vector for
    // neutralino masses requested on the squark & gluino masses.
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };






  // inline functions:

  inline void
  acceptanceCutSet::becomeCopyOf( acceptanceCutSet* const copyPointer )
  // this copies the values from a given acceptanceCutSet.
  {
    beamEnergy = copyPointer->getBeamEnergy();
    primaryLeptonCut = copyPointer->getPrimaryLeptonCut();
    secondaryLeptonCut = copyPointer->getSecondaryLeptonCut();
    jetCut = copyPointer->getJetCut();
  }

  inline bool
  acceptanceCutSet::compareJetAndBothLeptonCuts(
                                          acceptanceCutSet* const firstPointer,
                                        acceptanceCutSet* const secondPointer )
  // this returns true if both acceptanceCutSets have the same jet cut & both
  // lepton cuts.
  {
    if( firstPointer->isSameAcceptanceCutSet( secondPointer ) )
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  inline bool
  acceptanceCutSet::compareJustBothLeptonCuts(
                                          acceptanceCutSet* const firstPointer,
                                        acceptanceCutSet* const secondPointer )
  {
    if( ( firstPointer->getPrimaryLeptonCut()
          == secondPointer->getPrimaryLeptonCut() )
        &&
        ( firstPointer->getSecondaryLeptonCut()
          == secondPointer->getSecondaryLeptonCut() ) )
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  inline bool
  acceptanceCutSet::compareJustPrimaryLeptonCut(
                                          acceptanceCutSet* const firstPointer,
                                        acceptanceCutSet* const secondPointer )
  {
    if( firstPointer->getPrimaryLeptonCut()
        == secondPointer->getPrimaryLeptonCut() )
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  inline bool
  acceptanceCutSet::justReturnTrue( acceptanceCutSet* const firstPointer,
                                    acceptanceCutSet* const secondPointer )
  {
    return true;
  }

  inline double
  acceptanceCutSet::getBeamEnergy()
  const
  {

    return beamEnergy;

  }

  inline double
  acceptanceCutSet::getPrimaryLeptonCut()
  const
  {

    return primaryLeptonCut;

  }

  inline void
  acceptanceCutSet::setPrimaryLeptonCut( double inputValue )
  {

    primaryLeptonCut = inputValue;

  }

  inline double
  acceptanceCutSet::getSecondaryLeptonCut()
  const
  {

    return secondaryLeptonCut;

  }

  inline void
  acceptanceCutSet::setSecondaryLeptonCut( double inputValue )
  {

    secondaryLeptonCut = inputValue;

  }

  inline double
  acceptanceCutSet::getLeptonCut()
  const
  {

    return getPrimaryLeptonCut();

  }

  inline void
  acceptanceCutSet::setLeptonCut( double inputValue )
  {

    setPrimaryLeptonCut( inputValue );
    setSecondaryLeptonCut( inputValue );

  }

  inline double
  acceptanceCutSet::getJetCut()
  const
  {

    return jetCut;

  }

  inline void
  acceptanceCutSet::setJetCut( double inputValue )
  {

    jetCut = inputValue;

  }

  inline bool
  acceptanceCutSet::isSameAcceptanceCutSet(
                              acceptanceCutSet const* const comparisonCuts )
  const
  {

    if( ( comparisonCuts->getBeamEnergy() == beamEnergy )
        &&
        ( comparisonCuts->getJetCut() == jetCut )
        &&
        ( comparisonCuts->getPrimaryLeptonCut() == primaryLeptonCut )
        &&
        ( comparisonCuts->getSecondaryLeptonCut()
          == secondaryLeptonCut ) )
    {

      return true;

    }
    else
    {

      return false;

    }

  }



  inline int
  acceptanceGrid::getNumberOfAcceptanceColumns()
  const
  {
    return acceptanceColumns;
  }

  inline double
  acceptanceGrid::getLowestSquarkMass()
  const
  {
    return lowestSquarkMass;
  }

  inline double
  acceptanceGrid::getHighestSquarkMass()
  const
  {
    return  highestSquarkMass;
  }

  inline double
  acceptanceGrid::getLowestGluinoMass()
  const
  {
    return  lowestGluinoMass;
  }

  inline double
  acceptanceGrid::getHighestGluinoMass()
  const
  {
    return  highestGluinoMass;
  }

  inline double
  acceptanceGrid::neutralinoIndependentValueAt( double const squarkMass,
                                                double const gluinoMass,
                                                int const acceptanceElement )
  const
  // this interpolates for acceptanceElement assuming that the neutralino
  // masses are on the light-light grid point.
  {
    return vectorElementAt( squarkMass,
                            gluinoMass,
                            2,
                            acceptanceElement );
    // this is format-dependent, & corresponds to 2 being the index for
    // light-light neutralinos in the grid.
;  }
}  // end of LHC_FASER namespace.

#endif /* LHC_FASER_BASE_KINEMATICS_STUFF_HPP_ */
