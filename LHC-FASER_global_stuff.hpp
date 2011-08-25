/*
 * LHC-FASER_global_stuff.hpp
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


#ifndef LHC_FASER_GLOBAL_STUFF_HPP_
#define LHC_FASER_GLOBAL_STUFF_HPP_

#include <stdlib.h>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <list>
#include <map>
#include "LHC-FASER_template_classes.hpp"


namespace LHC_FASER
{
  // this class has static methods for miscellaneous useful functions.
  class lhcFaserGlobal
  {
  public:
    static double
    unitLinearInterpolation( double const xFraction,
                             double const leftValue,
                             double const rightValue )
    // this does a linear interpolation from the values at the given 2
    // ends of a unit line segment for the given point in the segment.
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    static double
    squareBilinearInterpolation( double const xFraction,
                                 double const yFraction,
                                 double const lowerLeftValue,
                                 double const lowerRightValue,
                                 double const upperRightValue,
                                 double const upperLeftValue )
    // this does a bilinear interpolation from the values at the given 4
    // corners of a unit square for the given point in the unit square.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    static double const negligibleBr;
    // this is what we use to decide when we throw away cascades because the
    // acceptance times branching ratio dropped too low.
    static double const negligibleSigma;
    // this is what we use to decide when we throw away channels because the
    // cross-section times acceptance times branching ratio dropped too low.
  };  // end of lhcFaserGlobal class.


  /* this class does the job of a subject, observed by
   * getsReadiedForNewPoint objects. the
   * getsReadiedForNewPoint objects are responsible for registering &
   * de-registering themselves with their constructors & destructors.
   * this class sets the 1st bools of the pairs observed by its observers to
   * true when requested. it sets the 2nd bools to false when the
   * readierForNewPoint instance's destructor is called.
   */
  class readierForNewPoint
  {
  public:
    readierForNewPoint()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~readierForNewPoint()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    includeMe( bool* observerBoolForReadierExistence )
    // this registers the observer so that the readier's destructor can let
    // the observer know that the readier is retiring.
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    removeMe( bool* observerBoolForReadierExistence,
              bool* observerBoolForReadying )
    /* this removes the observer from the readier's list of observers which is
     * used when the readier's destructor is called to let the observers know
     * to stop asking to be reset.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    resetMeOnNextUpdate( bool* observerBoolForReadying )
    // this adds a pointer from the observer to the list of bools to be set to
    // true with the next call of readyObserversForNewPoint().
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    readyObserversForNewPoint()
    // this goes through the list of pointers to bools created by observers
    // calling resetMe() & sets all the bools to true.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    std::list< bool* > observerBoolsForReadying;
    std::list< bool* > observerBoolsForReadierExistence;
  };  // end of readierForNewPoint class.


  /* this class does the job of being an interface for objects which need to
   * observe a readierForNewPoint subject. it's not mimicking a Java-style
   * interface properly, since it is not an abstract class (no abstract
   * functions), so there is the possibility of the deadly diamond of death if
   * any class which inherits from this class also inherits from another class
   * which also defines readyObserversForNewPoint(). however, you have been warned.
   * I am not going to write any class like that, & if you write such a class,
   * it is entirely your own silly fault.  ;]
   * the constructor registers the object with the readierForNewPoint & the
   * destructor de-registers itself from the readierForNewPoint.
   */
  class getsReadiedForNewPoint
  {
  public:
    getsReadiedForNewPoint( readierForNewPoint* const readier )
    /* the constructor gives a pointer to a bool to given_readier which
     * has a value that given_readier's destructor changes. this object keeps
     * a pointer to given_readier so that it can de-register when its
     * destructor is called.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    virtual
    ~getsReadiedForNewPoint()
    // the destructor tells readier, if it still exists, to stop modifying its
    // bool pointers.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    bool needsToPrepare;
    bool readierStillExists;
    readierForNewPoint* const readier;

    bool
    needsToPrepareForThisPoint()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    finishPreparingForThisPoint()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };  // end of getsReadiedForNewPoint class.


  // this class is a wrapper to allow for a getsReadiedForNewPoint that is
  // publically manipulable.
  class publicGetsReadiedForNewPoint : public getsReadiedForNewPoint
  {
  public:
    publicGetsReadiedForNewPoint( readierForNewPoint* const readier )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    virtual
    ~publicGetsReadiedForNewPoint()
    // the destructor tells readier, if it still exists, to stop modifying its
    // bool pointers.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    bool
    publicNeedsToPrepareForThisPoint()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    publicFinishPreparingForThisPoint()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  //protected:
    // nothing.
  };




  // inline functions:

  inline double
  lhcFaserGlobal::unitLinearInterpolation( double const xFraction,
                                           double const leftValue,
                                           double const rightValue )
  /* this does a linear interpolation from the values at the given 2
   * ends of a unit line segment for the given point in the segment.
   */
  {
    return ( leftValue + ( rightValue - leftValue ) * xFraction );
  }



  inline void
  readierForNewPoint::includeMe( bool* observerBoolForReadierExistence )
  // this registers the observer so that the readier's destructor can let
  // the observer know that the readier is retiring.
  {
    observerBoolsForReadierExistence.push_back(
                                             observerBoolForReadierExistence );
  }

  inline void
  readierForNewPoint::resetMeOnNextUpdate( bool* observerBoolForReadying )
  // this adds a pointer from the observer to the list of bools to be set to
  // true with the next call of readyObserversForNewPoint().
  {
    observerBoolsForReadying.push_back( observerBoolForReadying );
  }

  inline void
  readierForNewPoint::readyObserversForNewPoint()
  // this goes through the list of pointers to bools created by observers
  // calling resetMeOnNextUpdate() & sets all the bools to true.
  {
    for( std::list< bool* >::iterator
         boolIterator = observerBoolsForReadying.begin();
         observerBoolsForReadying.end() != boolIterator;
         ++boolIterator )
      // go through the list of observed bools...
      {
        // set the observer's bool to true so that it knows that an update has
        // happened:
        *(*boolIterator) = true;
      }
    observerBoolsForReadying.clear();
    // the list needs to be reset afterwards. pointers to bools are added in
    // again by observers calling resetMeOnNextUpdate().
  }



  inline bool
  getsReadiedForNewPoint::needsToPrepareForThisPoint()
  const
  {
    return needsToPrepare;
  }

  inline void
  getsReadiedForNewPoint::finishPreparingForThisPoint()
  {
    needsToPrepare = false;
    readier->resetMeOnNextUpdate( &needsToPrepare );
  }



  inline bool
  publicGetsReadiedForNewPoint::publicNeedsToPrepareForThisPoint()
  const
  {
    return needsToPrepareForThisPoint();
  }

  inline void
  publicGetsReadiedForNewPoint::publicFinishPreparingForThisPoint()
  {
    finishPreparingForThisPoint();
  }
}  // end of LHC_FASER namespace.

#endif /* LHC_FASER_GLOBAL_STUFF_HPP_ */
