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
                                bool const secondIsNotAntiparticleFlag )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    signedParticleShortcutPair(
                           signedParticleShortcutPair const* const copySource )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~signedParticleShortcutPair()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particlePointer
    getFirstParticle()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    bool
    firstIsNotAntiparticle()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getSecondParticle()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    bool
    secondIsNotAntiparticle()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    particlePointer firstParticlePointer;
    bool firstIsNotAntiparticleFlag;
    particlePointer secondParticlePointer;
    bool secondIsNotAntiparticleFlag;
  };  // end of signedParticleShortcutPair class.


  // this class just checks whether the given particle can decay directly to
  // any of the particles in the given vector, or their antiparticles.
  class directDecayChecker : public getsReadiedForNewPoint
  {
  public:
    directDecayChecker( particlePointer const decayingParticle,
                        particleVectorPointer const productSet,
                        readierForNewPoint* const readier )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~directDecayChecker()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particlePointer
    getDecayer()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particleVectorPointer
    getProducts()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    bool
    hasDirectDecay()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    particlePointer decayingParticle;
    particleVectorPointer productSet;
    bool decayerDoesDecayToAtLeastOneProduct;

    bool
    checkForDecay()
    // this checks the particles to see if the decay exists for this point.
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };  // end of directDecayChecker class.


  // this class holds a set of directDecayCheckers & passes out pointers to
  // them.
  class decayCheckerHandler
  {
  public:
    decayCheckerHandler( readierForNewPoint* const readier )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~decayCheckerHandler()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    directDecayChecker*
    getDecayChecker( particlePointer const decayingParticle,
                     particleVectorPointer const productSet )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    std::vector< directDecayChecker* > decayCheckers;
    readierForNewPoint* readier;

    directDecayChecker*
    addDecayChecker( particlePointer const decayingParticle,
                     particleVectorPointer const productSet  )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };  // end of decayCheckerHandler class.


  /* this class calculates the branching ratio for the given decayingParticle into the
   * given productParticle for all channels *excluding* those which also have any of
   * the particles in the exclusion list.
   */
  class exclusiveBrCalculator : public getsReadiedForNewPoint
  {
  public:
    exclusiveBrCalculator( particlePointer const decayingParticle,
                           particlePointer const productParticle,
                           bool const productIsNotAntiparticle,
                           std::list< int > const* const exclusionList,
                           readierForNewPoint* const readier )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~exclusiveBrCalculator()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particlePointer
    getDecayer()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getProduct()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    std::list< int > const*
    getExcluded()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getBr()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    particlePointer decayingParticle;
    particlePointer productParticle;
    std::list< int >* productCode;
    std::list< int > const* const exclusionList;
    double branchingRatioValue;

    double
    calculateBranchingRatio()
    const
    // this calculates the branching ratio for this point.
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };  // end of exclusiveBrCalculator class.


  // this class holds a set of exclusiveBrs & passes out pointers to
  // them.
  class exclusiveBrHandler
  {
  public:
    exclusiveBrHandler( particleVectorPointer const sdowns,
                        particleVectorPointer const sups,
                        readierForNewPoint* const readier )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~exclusiveBrHandler()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    exclusiveBrCalculator*
    getBrCalculator( particlePointer const decayingParticle,
                     particlePointer const productParticle,
                     bool const productIsNotAntiparticle,
                     std::list< int > const* const exclusionList )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    std::vector< exclusiveBrCalculator* > exclusiveBrs;
    readierForNewPoint* readier;
    std::vector< std::pair< int,
                            int >* > alwaysNeglectedDecays;

    exclusiveBrCalculator*
    addBrCalculator( particlePointer const decayingParticle,
                     particlePointer const productParticle,
                     bool const productIsNotAntiparticle,
                     std::list< int > const* const exclusionList )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };  // end of exclusiveBrHandler class.





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



  inline particlePointer
  directDecayChecker::getDecayer()
  const
  {
    return decayingParticle;
  }

  inline particleVectorPointer
  directDecayChecker::getProducts()
  const
  {
    return productSet;
  }

  inline bool
  directDecayChecker::hasDirectDecay()
  {
    if( needsToPrepareForThisPoint() )
    {
      decayerDoesDecayToAtLeastOneProduct = checkForDecay();
      finishPreparingForThisPoint();
    }
    return decayerDoesDecayToAtLeastOneProduct;
  }

  inline bool
  directDecayChecker::checkForDecay()
  // this checks the particles to see if the decay exists for this point.
  {
    bool returnBool( false );
    // check each possible decay mode:
    for( std::vector< particlePointer >::const_iterator
         productIterator = productSet->begin();
         productSet->end() > productIterator;
         ++productIterator )
      {
        if( decayingParticle->get_absolute_mass()
            > (*productIterator)->get_absolute_mass() )
          // if the decayingParticle can decay into this productParticle...
          {
            // note that the decayingParticle decays to at least 1 of the productSet in
            // productSet with a 2-body decay:
            returnBool = true;
            productIterator = productSet->end();
            // stop searching.
          }
      }
    return returnBool;
  }



  inline directDecayChecker*
  decayCheckerHandler::addDecayChecker( particlePointer const decayingParticle,
                                       particleVectorPointer const productSet )
  {
    decayCheckers.push_back( new directDecayChecker( decayingParticle,
                                                     productSet,
                                                     readier ) );
    return decayCheckers.back();
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
      finishPreparingForThisPoint();
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
                                                       readier ) );
    return exclusiveBrs.back();
  }
}  // end of LHC_FASER namespace


#endif /* LHC_FASER_SPARTICLE_DECAY_STUFF_HPP_ */
