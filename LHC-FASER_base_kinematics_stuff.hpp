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

#ifndef LHC_FASER_BASE_KINEMATICS_STUFF_HPP_
#define LHC_FASER_BASE_KINEMATICS_STUFF_HPP_

#include "LHC-FASER_input_handling_stuff.hpp"


/* this file contains the common code relevant to obtaining the acceptance for
 * jet transverse momentum plus missing transverse momentum and lepton
 * transverse momentum and pseudorapidity cuts. the code specific to lepton or
 * jet acceptancesPerCutSet are in LHC-FASER_lepton_kinematics_stuff.hpp &
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
    acceptanceCutSet();
    // this constructor sets the acceptances to -1.0, which is the default
    // "unset" value.
    acceptanceCutSet( acceptanceCutSet const* const copySource );
    // this constructor copies the acceptances from a given acceptanceCutSet.
    virtual
    ~acceptanceCutSet();

    static bool
    compareJetAndBothLeptonCuts( acceptanceCutSet const* const firstPointer,
                                 acceptanceCutSet const* const secondPointer );
    // this returns true if both acceptanceCutSets have the same jet cut & both
    // lepton cuts.
    static bool
    compareJustBothLeptonCuts( acceptanceCutSet const* const firstPointer,
                               acceptanceCutSet const* const secondPointer );
    // this returns true if both acceptanceCutSets have the same
    // primaryLeptonCut & secondaryLeptonCut, regardless of the jet cut.
    static bool
    compareJustPrimaryLeptonCut( acceptanceCutSet const* const firstPointer,
                                 acceptanceCutSet const* const secondPointer );
    // this returns true if both acceptanceCutSets have the same
    // primaryLeptonCut, regardless of the other acceptanceValues.
    static bool
    justReturnTrue( acceptanceCutSet const* const firstPointer,
                    acceptanceCutSet const* const secondPointer );
    // this just returns true, for when the acceptance cut doesn't matter.
    void
    becomeCopyOf( acceptanceCutSet const* const copyPointer );
    // this copies the acceptanceValues from a given acceptanceCutSet.
    int
    getBeamEnergy()
    const;
    // setBeamEnergy is only defined in the signalDefinitions derived class.
    bool
    leptonCutsAreEqual()
    const;
    double
    getPrimaryLeptonCut()
    const;
    void
    setPrimaryLeptonCut( double inputValue );
    double
    getSecondaryLeptonCut()
    const;
    void
    setSecondaryLeptonCut( double inputValue );
    double
    getLeptonCut()
    const;
    // this just returns primaryLeptonCut. it is intended to be used when the
    // signal only has 1 lepton cut.
    void
    setLeptonCut( double inputValue );
    // this sets both primaryLeptonCut & secondaryLeptonCut to inputValue.
    // it is intended to be used when the signal only has 1 lepton cut.
    double
    getJetCut()
    const;
    void
    setJetCut( double inputValue );
    std::list< int > const*
    getExcludedStandardModelProducts()
    const;
    void
    setExcludedStandardModelProducts( std::list< int > const* inputList );
    virtual bool
    isSameAcceptanceCutSet( acceptanceCutSet const* const comparisonCuts )
    const;

  protected:
    int beamEnergy;
    // the lepton transverse momentum cuts are easy enough to scale that we can
    // let the user specify them:
    double primaryLeptonCut;
    double secondaryLeptonCut;
    double jetCut;
    bool leptonCutsAreEqualFlag;
    std::list< int > const* excludedStandardModelProducts;

    void
    updateLeptonCutsAreEqualFlag();
  };


  /* this class reads in a file in the assumed format, stores it, & gives out
   * interpolated values. it was written with acceptances for leptons in mind,
   * to be adapted for jets plus missing transverse momentum as a special case,
   * with data files in the format
   * squarkMass gluinoMass lighterNeutralinoMass heavierNeutralinoMass
   * (continued) then either
   * 42 lepton acceptance values (effective squark mass, pseudorapidity cut
   * acceptance, then 40 bins for transverse momentum cut acceptancesPerCutSet)
   * or
   * 7 acceptance values for different choices of which jets+MET combination to
   * use
   * then newline.
   */
  class acceptanceGrid
  {
  public:
    acceptanceGrid( std::string const* const gridFileLocation,
                    inputHandler const* const inputShortcut );
    ~acceptanceGrid();

    double
    getLowestSquarkMass()
    const;
    double
    getHighestSquarkMass()
    const;
    double
    getLowestGluinoMass()
    const;
    double
    getHighestGluinoMass()
    const;

    double
    testFunction( double const squarkMass,
                  double const gluinoMass,
                  int const electroweakinoIndex,
                  int const acceptanceIndex );

  protected:
    static unsigned int const indexForLightLightNeutralinoPair;
    static unsigned int const indexForLightMediumNeutralinoPair;
    static unsigned int const indexForLightHeavyNeutralinoPair;
    static unsigned int const indexForMediumMediumNeutralinoPair;
    static unsigned int const indexForMediumHeavyNeutralinoPair;
    static unsigned int const indexForHeavyHeavyNeutralinoPair;

    inputHandler const* const inputShortcut;
    double scoloredMassStepSize;
    double lowestSquarkMass;
    double highestSquarkMass;
    double lowestGluinoMass;
    double highestGluinoMass;
    double lowElectroweakinoMassRatio;
    double mediumElectroweakinoMassRatio;
    double highElectroweakinoMassRatio;
    std::vector< std::vector< std::vector< std::vector< double >* >* >* >
    acceptanceValues;
    // these are pointers to vectors of pointers to vectors of acceptances
    // which are used to store the acceptances nearest the current point:
    std::vector< std::vector< double >* >* lowerLeftVectorOfVectors;
    std::vector< std::vector< double >* >* lowerRightVectorOfVectors;
    std::vector< std::vector< double >* >* upperRightVectorOfVectors;
    std::vector< std::vector< double >* >* upperLeftVectorOfVectors;
    double squarkMassFraction;
    // this is the fraction of a step in the squark mass direction that the
    // current point is from the grid point with squark mass just below.
    double gluinoMassFraction;
    // this is the fraction of a step in the gluino mass direction that the
    // current point is from the grid point with gluino mass just below.

    bool
    pointIsOnGrid( double const squarkMass,
                   double const gluinoMass );
    /* this looks to see if the squark & gluino masses are on the grid, & if
     * so, sets lowerLeftVectorOfVectors, lowerRightVectorOfVectors,
     * upperRightVectorOfVectors, upperLeftVectorOfVectors, squarkMassFraction,
     * & gluinoMassFraction appropriately, returning true. otherwise, it
     * returns false.
     */
  };





  // inline functions:


  inline void
  acceptanceCutSet::becomeCopyOf( acceptanceCutSet const* const copyPointer )
  // this copies the acceptanceValues from a given acceptanceCutSet.
  {
    beamEnergy = copyPointer->getBeamEnergy();
    primaryLeptonCut = copyPointer->getPrimaryLeptonCut();
    secondaryLeptonCut = copyPointer->getSecondaryLeptonCut();
    jetCut = copyPointer->getJetCut();
  }

  inline bool
  acceptanceCutSet::compareJetAndBothLeptonCuts(
                                    acceptanceCutSet const* const firstPointer,
                                  acceptanceCutSet const* const secondPointer )
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
                                    acceptanceCutSet const* const firstPointer,
                                  acceptanceCutSet const* const secondPointer )
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
                                    acceptanceCutSet const* const firstPointer,
                                  acceptanceCutSet const* const secondPointer )
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
  acceptanceCutSet::justReturnTrue( acceptanceCutSet const* const firstPointer,
                                  acceptanceCutSet const* const secondPointer )
  {
    return true;
  }

  inline int
  acceptanceCutSet::getBeamEnergy()
  const
  {
    return beamEnergy;
  }

  inline bool
  acceptanceCutSet::leptonCutsAreEqual()
  const
  {
    return leptonCutsAreEqualFlag;
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
    updateLeptonCutsAreEqualFlag();
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
    updateLeptonCutsAreEqualFlag();
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
    leptonCutsAreEqualFlag = true;
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

  inline std::list< int > const*
  acceptanceCutSet::getExcludedStandardModelProducts()
  const
  {
    return excludedStandardModelProducts;
  }

  inline void
  acceptanceCutSet::setExcludedStandardModelProducts(
                                            std::list< int > const* inputList )
  {
    excludedStandardModelProducts = inputList;
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

  inline void
  acceptanceCutSet::updateLeptonCutsAreEqualFlag()
  {
    if( primaryLeptonCut == secondaryLeptonCut )
    {
      leptonCutsAreEqualFlag = true;
    }
    else
    {
      leptonCutsAreEqualFlag = false;
    }
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
    return highestSquarkMass;
  }

  inline double
  acceptanceGrid::getLowestGluinoMass()
  const
  {
    return lowestGluinoMass;
  }

  inline double
  acceptanceGrid::getHighestGluinoMass()
  const
  {
    return highestGluinoMass;
  }

}  // end of LHC_FASER namespace.

#endif /* LHC_FASER_BASE_KINEMATICS_STUFF_HPP_ */
