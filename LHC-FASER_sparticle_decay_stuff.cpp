/*
 * LHC-FASER_sparticle_decay_stuff.cpp
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
 *      LHC-FASER also requires grids of lookup acceptanceValues. These should also be
 *      found in a subdirectory included with this package.
 */

#include "LHC-FASER_sparticle_decay_stuff.hpp"

namespace LHC_FASER
{
  signedParticleShortcutPair::signedParticleShortcutPair(
                                    particlePointer const firstParticlePointer,
                                         bool const firstIsNotAntiparticleFlag,
                                   particlePointer const secondParticlePointer,
                                     bool const secondIsNotAntiparticleFlag ) :
    firstParticlePointer( firstParticlePointer ),
    firstIsNotAntiparticleFlag( firstIsNotAntiparticleFlag ),
    secondParticlePointer( secondParticlePointer ),
    secondIsNotAntiparticleFlag( secondIsNotAntiparticleFlag )
  {
    // just an initialization list.
  }

  signedParticleShortcutPair::signedParticleShortcutPair(
                         signedParticleShortcutPair const* const copySource ) :
    firstParticlePointer( copySource->getFirstParticle() ),
    firstIsNotAntiparticleFlag( copySource->firstIsNotAntiparticle() ),
    secondParticlePointer( copySource->getSecondParticle() ),
    secondIsNotAntiparticleFlag( copySource->secondIsNotAntiparticle() )
  {
    // just an initialization list.
  }

  signedParticleShortcutPair::~signedParticleShortcutPair()
  {
    // does nothing.
  }



  directDecayChecker::directDecayChecker(
                                        particlePointer const decayingParticle,
                                        particleVectorPointer const productSet,
                                          readierForNewPoint* const readier ) :
    getsReadiedForNewPoint( readier ),
    decayingParticle( decayingParticle ),
    productSet( productSet ),
    decayerDoesDecayToAtLeastOneProduct( false )
  {
    // just an initialization list.
  }

  directDecayChecker::~directDecayChecker()
  {
    // does nothing.
  }



  decayCheckerHandler::decayCheckerHandler(
                                          readierForNewPoint* const readier ) :
    readier( readier )
  {
    // just an initialization list.
  }

  decayCheckerHandler::~decayCheckerHandler()
  {
    for( std::vector< directDecayChecker* >::iterator
         deletionIterator = decayCheckers.begin();
         decayCheckers.end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
  }


  directDecayChecker*
  decayCheckerHandler::getDecayChecker( particlePointer const decayingParticle,
                                       particleVectorPointer const productSet )
  {
    directDecayChecker* returnPointer( NULL );
    for( std::vector< directDecayChecker* >::const_iterator
         checkerIterator = decayCheckers.begin();
         decayCheckers.end() > checkerIterator;
         ++checkerIterator )
    {
      if( ( decayingParticle == (*checkerIterator)->getDecayer() )
          &&
          ( productSet == (*checkerIterator)->getProducts() ) )
        // if we find a perfect match...
      {
        returnPointer = *checkerIterator;
        // return this exclusiveBrCalculator.
        checkerIterator = decayCheckers.end();
        // stop looking.
      }
    }
    if( NULL == returnPointer )
      // if we didn't find a pre-existing match...
    {
      returnPointer = addDecayChecker( decayingParticle,
                                       productSet );
      // add this calculator to the collection.
    }
    return returnPointer;
  }



  exclusiveBrCalculator::exclusiveBrCalculator(
                                        particlePointer const decayingParticle,
                                         particlePointer const productParticle,
                                           bool const productIsNotAntiparticle,
                                   std::list< int > const* const exclusionList,
                                          readierForNewPoint* const readier ) :
    getsReadiedForNewPoint( readier ),
    decayingParticle( decayingParticle ),
    productParticle( productParticle ),
    exclusionList( exclusionList ),
    branchingRatioValue( 0.0 )
  {
    if( productIsNotAntiparticle )
    {
      productCode = new std::list< int >( 1,
                                          productParticle->get_PDG_code() );
    }
    else
    {
      productCode = new std::list< int >( 1,
                                          -(productParticle->get_PDG_code()) );
    }
  }

  exclusiveBrCalculator::~exclusiveBrCalculator()
  {
    delete productCode;
  }



  exclusiveBrHandler::exclusiveBrHandler( particleVectorPointer const sdowns,
                                          particleVectorPointer const sups,
                                          readierForNewPoint* const readier ) :
    readier( readier )
  {
    for( std::vector< particlePointer >::const_iterator
         sdownIterator = sdowns->begin();
         sdowns->end() > sdownIterator;
         ++sdownIterator )
    {
      alwaysNeglectedDecays.push_back( new std::pair< int,
                                                      int >(
                                              (*sdownIterator)->get_PDG_code(),
                                           CppSLHA::PDG_code::chargino_one ) );
      alwaysNeglectedDecays.push_back( new std::pair< int,
                                                      int >(
                                              (*sdownIterator)->get_PDG_code(),
                                           CppSLHA::PDG_code::chargino_two ) );
    }
    for( std::vector< particlePointer >::const_iterator
         supIterator = sups->begin();
         sups->end() > supIterator;
         ++supIterator )
    {
      alwaysNeglectedDecays.push_back( new std::pair< int,
                                                      int >(
                                                (*supIterator)->get_PDG_code(),
                                        -(CppSLHA::PDG_code::chargino_one ) ) );
      alwaysNeglectedDecays.push_back( new std::pair< int,
                                                      int >(
                                                (*supIterator)->get_PDG_code(),
                                       -(CppSLHA::PDG_code::chargino_two ) ) );
    }
  }

  exclusiveBrHandler::~exclusiveBrHandler()
  {
    for( std::vector< exclusiveBrCalculator* >::iterator
         deletionIterator = exclusiveBrs.begin();
         exclusiveBrs.end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
    for( std::vector< std::pair< int,
                                 int >* >::iterator
         deletion_iterator = alwaysNeglectedDecays.begin();
         alwaysNeglectedDecays.end() > deletion_iterator;
         ++deletion_iterator )
    {
      delete *deletion_iterator;
    }
  }


  exclusiveBrCalculator*
  exclusiveBrHandler::getBrCalculator( particlePointer const decayingParticle,
                                       particlePointer const productParticle,
                                       bool const productIsNotAntiparticle,
                                  std::list< int > const* const exclusionList )
  {
    exclusiveBrCalculator* returnPointer( NULL );
    int productCode;
    if( productIsNotAntiparticle )
    {
      productCode = productParticle->get_PDG_code();
    }
    else
    {
      productCode = -(productParticle->get_PDG_code());
    }
    bool notAlwaysNeglected( true );
    for( std::vector< std::pair< int,
                                 int >* >::iterator
         neglectionIterator = alwaysNeglectedDecays.begin();
         alwaysNeglectedDecays.end() > neglectionIterator;
         ++neglectionIterator )
      // check whether this decay is always neglected because it cannot happen
      // often enough no matter what the scenario:
    {
      if( decayingParticle->get_PDG_code() == (*neglectionIterator)->first )
      {
        if( productCode == (*neglectionIterator)->second )
        {
          notAlwaysNeglected = false;
          neglectionIterator = alwaysNeglectedDecays.end();
          // stop looking.
        }
      }
    }
    if( notAlwaysNeglected )
    {
      for( std::vector< exclusiveBrCalculator* >::iterator
           brIterator = exclusiveBrs.begin();
           exclusiveBrs.end() > brIterator;
           ++brIterator )
      {
        if( ( decayingParticle == (*brIterator)->getDecayer() )
            &&
            ( productParticle == (*brIterator)->getProduct() )
            &&
            ( exclusionList == (*brIterator)->getExcluded() ) )
          // if we find a perfect match...
        {
          returnPointer = *brIterator;
          // return this exclusiveBrCalculator.
          brIterator = exclusiveBrs.end();
          // stop looking.
        }
      }
      if( NULL == returnPointer )
        // if we didn't find a pre-existing match...
      {
        returnPointer = addBrCalculator( decayingParticle,
                                         productParticle,
                                         productIsNotAntiparticle,
                                         exclusionList );
        // add this calculator to the collection.
      }
    }
    // end of if( notAlwaysNeglected ) - if notAlwaysNeglected == false, then
    // NULL is returned. we should always check whether we got a NULL!
    return returnPointer;
  }
}  // end of LHC_FASER namespace
