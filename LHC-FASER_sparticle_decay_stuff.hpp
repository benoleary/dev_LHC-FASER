/*
 * LHC-FASER_sparticle_decay_stuff.hpp
 *
 *  Created on: 21 Mar 2011
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

#ifndef LHC_FASER_SPARTICLE_DECAY_STUFF_HPP_
#define LHC_FASER_SPARTICLE_DECAY_STUFF_HPP_

#include "CppSLHA/CppSLHA.hpp"
#include "LHC-FASER_global_stuff.hpp"

namespace LHC_FASER
{
  typedef CppSLHA::particle_property_set const* particlePointer;
  typedef std::vector< particlePointer > const* particleVectorPointer;
  // these are just for convenience.

  // this class holds a pair of particle_property_set along with a bool which
  // is true for particles, false for antiparticles.
  class signedParticleShortcutPair
  {
  public:
    signedParticleShortcutPair( particlePointer const firstParticlePointer,
                                bool const firstIsNotAntiparticleFlag,
                                particlePointer const secondParticlePointer,
                                bool const secondIsNotAntiparticleFlag );
    signedParticleShortcutPair(
                          signedParticleShortcutPair const* const copySource );
    ~signedParticleShortcutPair();

    particlePointer
    getFirstParticle()
    const;
    bool
    firstIsNotAntiparticle()
    const;
    particlePointer
    getSecondParticle()
    const;
    bool
    secondIsNotAntiparticle()
    const;
    std::string
    getAsString()
    const;
    int const
    getSignedFirstCode()
    const;
    int const
    getSignedSecondCode()
    const;

  protected:
    particlePointer firstParticlePointer;
    bool firstIsNotAntiparticleFlag;
    particlePointer secondParticlePointer;
    bool secondIsNotAntiparticleFlag;
  };


  /* this class calculates the branching ratio for the given decaying particle
   * into the given product particle for all channels *excluding* those which
   * also have any of the particles in the exclusion list.
   */
  class exclusiveBrCalculator : public getsReadiedForNewPoint
  {
  public:
    exclusiveBrCalculator( particlePointer const decayingParticle,
                           particlePointer const productParticle,
                           bool const productIsNotAntiparticle,
                           std::list< int > const* const exclusionList,
                           readierForNewPoint* const readierPointer );
    ~exclusiveBrCalculator();

    particlePointer
    getDecayer()
    const;
    particlePointer
    getProduct()
    const;
    std::list< int > const*
    getExcluded()
    const;
    double
    getBr();

  protected:
    particlePointer decayingParticle;
    particlePointer productParticle;
    std::list< int >* productCode;
    std::list< int > const* const exclusionList;
    double branchingRatioValue;

    double
    calculateBranchingRatio()
    const;
    // this calculates the branching ratio for this point.
  };


  // this class holds a set of exclusiveBrs & passes out pointers to them.
  class exclusiveBrHandler
  {
  public:
    exclusiveBrHandler( particleVectorPointer const sdownTypes,
                        particleVectorPointer const supTypes,
                        readierForNewPoint* const readierPointer );
    ~exclusiveBrHandler();

    exclusiveBrCalculator*
    getBrCalculator( particlePointer const decayingParticle,
                     particlePointer const productParticle,
                     bool const productIsNotAntiparticle,
                     std::list< int > const* const exclusionList );

  protected:
    std::vector< exclusiveBrCalculator* > exclusiveBrs;
    readierForNewPoint* readierPointer;
    std::vector< std::pair< int,
                            int >* > alwaysNeglectedDecays;

    exclusiveBrCalculator*
    addBrCalculator( particlePointer const decayingParticle,
                     particlePointer const productParticle,
                     bool const productIsNotAntiparticle,
                     std::list< int > const* const exclusionList );
  };





  // inline functions:


  inline particlePointer
  signedParticleShortcutPair::getFirstParticle()
  const
  {
    return firstParticlePointer;
  }

  inline bool
  signedParticleShortcutPair::firstIsNotAntiparticle()
  const
  {
    return firstIsNotAntiparticleFlag;
  }

  inline particlePointer
  signedParticleShortcutPair::getSecondParticle()
  const
  {
    return secondParticlePointer;
  }

  inline bool
  signedParticleShortcutPair::secondIsNotAntiparticle()
  const
  {
    return secondIsNotAntiparticleFlag;
  }

  inline std::string
  signedParticleShortcutPair::getAsString() const
  {
    std::string returnString( "{ " );
    returnString.append( *(firstParticlePointer->get_name_or_antiname(
                                               firstIsNotAntiparticleFlag )) );
    returnString.append( " + " );
    returnString.append( *(secondParticlePointer->get_name_or_antiname(
                                              secondIsNotAntiparticleFlag )) );
    returnString.append( " }" );
    return returnString;
  }

  inline int const
  signedParticleShortcutPair::getSignedFirstCode() const
  {
    if( firstIsNotAntiparticleFlag )
    {
      return firstParticlePointer->get_PDG_code();
    }
    else
    {
      return -(firstParticlePointer->get_PDG_code());
    }
  }

  inline int const
  signedParticleShortcutPair::getSignedSecondCode() const
  {
    if( secondIsNotAntiparticleFlag )
    {
      return secondParticlePointer->get_PDG_code();
    }
    else
    {
      return -(secondParticlePointer->get_PDG_code());
    }
  }



  inline particlePointer
  exclusiveBrCalculator::getDecayer()
  const
  {
    return decayingParticle;
  }

  inline particlePointer
  exclusiveBrCalculator::getProduct()
  const
  {
    return productParticle;
  }

  inline std::list< int > const*
  exclusiveBrCalculator::getExcluded()
  const
  {
    return exclusionList;
  }


  inline double
  exclusiveBrCalculator::getBr()
  {
    if( needsToPrepareForThisPoint() )
    {
      if( decayingParticle->get_absolute_mass()
          > productParticle->get_absolute_mass() )
        // we check 1st to see if the BR can be non-zero.
      {
        branchingRatioValue = calculateBranchingRatio();
      }
      else
      {
        branchingRatioValue = 0.0;
      }
    }
    return branchingRatioValue;
  }

  inline double
  exclusiveBrCalculator::calculateBranchingRatio()
  const
  // this calculates the branching ratio for this point.
  {
    return decayingParticle->inspect_direct_decay_handler(
                                )->get_branching_ratio_for_subset( productCode,
                                                               exclusionList );
  }



  inline exclusiveBrCalculator*
  exclusiveBrHandler::addBrCalculator( particlePointer const decayingParticle,
                                       particlePointer const productParticle,
                                       bool const productIsNotAntiparticle,
                                  std::list< int > const* const exclusionList )
  {
    exclusiveBrs.push_back( new exclusiveBrCalculator( decayingParticle,
                                                       productParticle,
                                                      productIsNotAntiparticle,
                                                       exclusionList,
                                                       readierPointer ) );
    return exclusiveBrs.back();
  }

}  // end of LHC_FASER namespace


#endif /* LHC_FASER_SPARTICLE_DECAY_STUFF_HPP_ */
