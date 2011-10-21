/*
 * LHC-FASER_template_classes.hpp
 *
 *  Created on: Jun 24, 2011
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
 *      LHC-FASER also requires grids of lookup acceptanceValues. These should also be
 *      found in a subdirectory included with this package.
 */

#ifndef LHC_FASER_TEMPLATE_CLASSES_HPP_
#define LHC_FASER_TEMPLATE_CLASSES_HPP_

#include <stdlib.h>
#include <math.h>
#include <utility>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <list>
#include <map>


namespace LHC_FASER
{
  /* this is a template class to hold a vector of pointers to objects. the main
   * use is when the objects are slow to construct but fast to reset.
   * it needs to be used with storedClass objects that have default
   * constructors that can be called with "new storedClass".
   * usage: the intended use is that addNewAtEnd() takes over the role of
   * memory allocation. instead of
   * storedClass* egPointer = new storedClass;
   * the allocation is done with
   * storedClass* egPointer = egMinimalAllocationVector.addNewAtEnd();
   * & then carry on as normal, except egMinimalAllocationVector will free up
   * the memory on its own, so you should *never* try to delete the pointer
   * yourself.
   * when the collection needs to be reset, call clearEntries(), but beware
   * that the next pointers returned by addNewAtEnd() will point to old objects
   * from before the calling of clearEntries()! you need to reset the objects
   * yourself.
   * I had made a version with a stored function pointer to allow the objects
   * to be reset as they are passed out, but in the end I felt that it couldn't
   * be flexible enough to be really useful. the functionality has been moved
   * into
   * addNewAtEnd( void (*makeAsDefaultConstructorWould)( storedClass* const ) )
   * for obtaining consistent new objects. (makeAsDefaultConstructorWould does
   * not _have_ to set the object's properties to be as they would be when
   * constructed by the default constructor, but this is the use intended.)
   */
  template< class storedClass >
  class minimalAllocationVector
  {
  public:
    minimalAllocationVector();
    ~minimalAllocationVector();

    storedClass*
    getPointer( int const requestedElement );
    // if an element out of range is asked for, this returns NULL.
    storedClass const*
    inspectPointer( int const requestedElement )
    const;
    // if an element out of range is asked for, this returns NULL.
    storedClass*
    getFront();
    // this returns the 1st element deemed current. if there is none, this
    // returns NULL.
    storedClass const*
    inspectFront()
    const;
    // this returns the 1st element deemed current. if there is none, this
    // returns NULL.
    storedClass*
    getBack();
    // this returns the last element deemed current. if there is none, this
    // returns NULL.
    storedClass const*
    inspectBack()
    const;
    // this returns the last element deemed current. if there is none, this
    // returns NULL.
    int
    getSize()
    const;
    // this returns the number of elements which are deemed current.
    int
    getLastIndex()
    const;
    // this returns the number for getPointer to access the last of elements
    // which are deemed current.
    std::vector< storedClass* >*
    getVector();
    // this returns a pointer to the std::vector holding pointers to the
    // storedClass instances which are deemed current.
    std::vector< storedClass* > const*
    inspectVector()
    const;
    // this returns a pointer to the std::vector holding pointers to the
    // storedClass instances which are deemed current.
    void
    clearEntries();
    // this notes that all current storedClass instances are now to be treated
    // as stored-but-not-current, & drops their pointers from currentPointers.
    storedClass*
    addNewAtEnd();
    /* this, if necessary, pushes back a new storedClass instance into
     * allConstructedPointers, & pushes back a pointer into currentPointers, &
     * returns the pointer at the end of currentPointers.
     */
    storedClass*
    addNewAtEnd( void (*makeAsDefaultConstructorWould)( storedClass* const ) );
    /* this pushes a pointer to a storedClass instance into currentPointers &
     * also returns that pointer. the returned pointer either was already in
     * allConstructedPointers, in which case it also gets reset in some way by
     * makeAsDefaultConstructorWould before being returned, or points to a new
     * storedClass instance which is constructed by this function with the
     * default storedClass constructor, & this pointer is also pushed back into
     * allConstructedPointers. N.B.: makeAsDefaultConstructorWould is NOT
     * applied to pointers which just made new storedClass instances - if you
     * want it applied to every pointer that comes out of addNewAtEnd(), you
     * might as well just use the argumentless addNewAtEnd() & run
     * makeAsDefaultConstructorWould on whatever comes out.
     * (makeAsDefaultConstructorWould does not _have_ to set the object's
     * properties to be as they would be when constructed by the default
     * constructor, but this is the use intended.)
     */

  protected:
    std::vector< storedClass* > allConstructedPointers;
    std::vector< storedClass* > currentPointers;
    void (*resetFunction)( storedClass* const );
  };


  /* this is a template class to hold a vector of pointers to objects, with
   * keys to find them. the main use is when the objects are slow to construct
   * but fast to reset.
   * it needs to be used with storedClass objects that can be constructed by
   * a single function which takes a single pointer to a keyClass object.
   */
  template< class storedClass,
            class keyClass >
  class minimalAllocationKeyedVector
  {
  public:
    minimalAllocationKeyedVector(
                             bool (*comparisonFunction)( keyClass const* const,
                                                       keyClass const* const ),
                             void (*constructionFunction)( storedClass** const,
                                                           keyClass** const,
                                                       keyClass const* const ),
                                  void (*resetFunction)( storedClass* const,
                                                         keyClass* const,
                                                     keyClass const* const ) );
    /* comparisonFunction should return true if the 2 arguments are "equal" in
     * whatever sense is required.
     * constructionFunction should assign a pointer to a newly-constructed
     * storedClass instance to the 1st argument, & a pointer to a
     * newly-constructed keyClass instance to the 2nd argument, with both
     * instances based on the 3rd argument: the constructed keyClass instance
     * should be such that comparisonFunction with it & the 3rd argument
     * returns true, & the constructed storedClass instance should be as it
     * should be to be associated with the newly-constructed keyClass instance;
     * i.e.
     * (*constructionFunction)( &egStoredClassPointer,
     *                          &egFirstKeyClassPointer,
     *                          egSecondKeyClassPointer ) should make a new
     * storedClass instance using egSecondKeyClassPointer, & assign
     * egStoredClassPointer to point at it, while also making a new
     * keyClass instance using egSecondKeyClassPointer, & assigning
     * egFirstKeyClassPointer to point at it (sorry about the mess with
     * pointers-to-pointers).
     * resetFunction should reset a storedClass instance & a keyClass instance
     * to be as they would be if they had been constructed with
     * constructionFunction, with the same arguments.
     */
    ~minimalAllocationKeyedVector();

    void
    setComparison( bool (*comparisonFunction)( keyClass const* const,
                                               keyClass const* const ) );
    storedClass*
    getPointer( keyClass const* const keyPointer );
    /* this returns a pointer to the storedClass instance that matches
     * keyPointer, if there is a stored pointer to such an instance, or it
     * constructs a new storedClass instance with keyPointer & returns the
     * pointer to that.
     */
    std::vector< storedClass* >*
    getVector();
    // this returns a pointer to the std::vector holding pointers to the
    // storedClass instances which are deemed current.
    void
    clearEntries();
    // this notes that all current storedClass instances are now to be treated
    // as stored-but-not-current, & drops their pointers from currentPointers.

  protected:
    std::vector< std::pair< storedClass*,
                            keyClass* >* > allConstructedPointers;
    std::vector< storedClass* > currentPointers;
    storedClass* returnPointer;
    keyClass* keyClassConstructionPointer;
    //unsigned int currentlyStoredNumber;
    unsigned int searchCounter;
    bool keyNotFound;
    bool (*comparisonFunction)( keyClass const* const,
                                keyClass const* const );
    void (*constructionFunction)( storedClass** const,
                                  keyClass** const,
                                  keyClass const* const );
    void (*resetFunction)( storedClass* const,
                           keyClass* const,
                           keyClass const* const );
  };





  // inline functions:


  template< class storedClass >
  inline
  minimalAllocationVector< storedClass >::minimalAllocationVector()
  {
    // just an empty initialization list.
  }

  template< class storedClass >
  inline
  minimalAllocationVector< storedClass >::~minimalAllocationVector()
  {
    for( unsigned int deletionCounter( 0 );
         allConstructedPointers.size() > deletionCounter;
         ++deletionCounter )
    {
      delete allConstructedPointers.at( deletionCounter );
    }
    // I cannot use an iterator here because of template complications. also,
    // I'm going from 0 upwards because of unsigned int complications.
  }


  template< class storedClass >
  inline storedClass*
  minimalAllocationVector< storedClass >::getPointer(
                                                   int const requestedElement )
  // if an element out of range is asked for, this returns NULL.
  {
    if( ( 0 <= requestedElement )
        &&
        ( currentPointers.size() > (unsigned int)requestedElement ) )
    {
      return currentPointers.at( requestedElement );
    }
    else
    {
      return NULL;
    }
  }

  template< class storedClass >
  inline storedClass const*
  minimalAllocationVector< storedClass >::inspectPointer(
                                                   int const requestedElement )
  const
  // if an element out of range is asked for, this returns NULL.
  {
    return getPointer( requestedElement );
  }

  template< class storedClass >
  inline storedClass*
  minimalAllocationVector< storedClass >::getFront()
  // this returns the 1st element deemed current. if there is none, this
  // returns NULL.
  {
    if( 0 < currentPointers.size() )
    {
      return currentPointers.front();
    }
    else
    {
      return NULL;
    }
  }

  template< class storedClass >
  inline storedClass const*
  minimalAllocationVector< storedClass >::inspectFront()
  const
  // this returns the 1st element deemed current. if there is none, this
  // returns NULL.
  {
    return getFront();
  }

  template< class storedClass >
  inline storedClass*
  minimalAllocationVector< storedClass >::getBack()
  // this returns the last element deemed current. if there is none, this
  // returns NULL.
  {
    if( 0 < currentPointers.size() )
    {
      return currentPointers.back();
    }
    else
    {
      return NULL;
    }
  }

  template< class storedClass >
  inline storedClass const*
  minimalAllocationVector< storedClass >::inspectBack()
  // this returns the 1st element deemed current. if there is none, this
  // returns NULL.
  const
  {
    return getBack();
  }

  template< class storedClass >
  inline int
  minimalAllocationVector< storedClass >::getSize()
  const
  // this returns the number of elements which are deemed current.
  {
    return (int)(currentPointers.size());
  }

  template< class storedClass >
  inline int
  minimalAllocationVector< storedClass >::getLastIndex()
  const
  // this returns the number for getPointer to access the last of elements
  // which are deemed current.
  {
    return ( getSize() - 1 );
  }

  template< class storedClass >
  inline std::vector< storedClass* >*
  minimalAllocationVector< storedClass >::getVector()
  // this returns a pointer to the std::vector holding pointers to the
  // storedClass instances which are deemed current.
  {
    return &currentPointers;
  }

  template< class storedClass >
  inline std::vector< storedClass* > const*
  minimalAllocationVector< storedClass >::inspectVector()
  const
  // this returns a pointer to the std::vector holding pointers to the
  // storedClass instances which are deemed current.
  {
    return &currentPointers;
  }

  template< class storedClass >
  inline void
  minimalAllocationVector< storedClass >::clearEntries()
  // this notes that all current storedClass instances are now to be treated
  // as stored-but-not-current, & drops their pointers from currentPointers.
  {
    currentPointers.clear();
  }

  template< class storedClass >
  inline storedClass*
  minimalAllocationVector< storedClass >::addNewAtEnd()
  /* this, if necessary, pushes back a new storedClass instance into
   * allConstructedPointers, & pushes back a pointer into currentPointers, &
   * returns the pointer at the end of currentPointers.
   */
  {
    while( currentPointers.size() >= allConstructedPointers.size() )
      // this should cover cases where currentPointers has been changed in size
      // by malevolent use of getVector...
    {
      allConstructedPointers.push_back( new storedClass );
    }
    currentPointers.push_back(
                         allConstructedPointers.at( currentPointers.size() ) );
    return currentPointers.back();
  }

  template< class storedClass >
  inline storedClass*
  minimalAllocationVector< storedClass >::addNewAtEnd(
                  void (*makeAsDefaultConstructorWould)( storedClass* const ) )
  /* this pushes a pointer to a storedClass instance into currentPointers &
   * also returns that pointer. the returned pointer either was already in
   * allConstructedPointers, in which case it also gets reset in some way by
   * makeAsDefaultConstructorWould before being returned, or points to a new
   * storedClass instance which is constructed by this function with the
   * default storedClass constructor, & this pointer is also pushed back into
   * allConstructedPointers. N.B.: makeAsDefaultConstructorWould is NOT
   * applied to pointers which just made new storedClass instances - if you
   * want it applied to every pointer that comes out of addNewAtEnd(), you
   * might as well just use the argumentless addNewAtEnd() & run
   * makeAsDefaultConstructorWould on whatever comes out.
   * (makeAsDefaultConstructorWould does not _have_ to set the object's
   * properties to be as they would be when constructed by the default
   * constructor, but this is the use intended.)
   */
  {
    if( currentPointers.size() < allConstructedPointers.size() )
    {
      currentPointers.push_back(
                         allConstructedPointers.at( currentPointers.size() ) );
      (*makeAsDefaultConstructorWould)( currentPointers.back() );
    }
    else
    {
      while( currentPointers.size() >= allConstructedPointers.size() )
      // this should cover cases where currentPointers has been changed in size
      // by malevolent use of getVector...
      {
        allConstructedPointers.push_back( new storedClass );
      }
      currentPointers.push_back(
                         allConstructedPointers.at( currentPointers.size() ) );
    }
    return currentPointers.back();
  }



  /* comparisonFunction should return true if the 2 arguments are "equal" in
   * whatever sense is required.
   * constructionFunction should assign a pointer to a newly-constructed
   * storedClass instance to the 1st argument, & a pointer to a
   * newly-constructed keyClass instance to the 2nd argument, with both
   * instances based on the 3rd argument: the constructed keyClass instance
   * should be such that comparisonFunction with it & the 3rd argument
   * returns true, & the constructed storedClass instance should be as it
   * should be to be associated with the newly-constructed keyClass instance;
   * i.e.
   * (*constructionFunction)( egStoredClassPointer,
   *                          egFirstKeyClassPointer,
   *                          egSecondKeyClassPointer ) should make a new
   * storedClass instance using egSecondKeyClassPointer, & assign
   * egStoredClassPointer to point at it, while also making a new
   * keyClass instance using egSecondKeyClassPointer, & assigning
   * egFirstKeyClassPointer to point at it.
   * resetFunction should reset a storedClass instance & a keyClass instance
   * to be as they would be if they had been constructed with
   * constructionFunction, with the same arguments.
   */
  template< class storedClass,
            class keyClass >
  inline
  minimalAllocationKeyedVector< storedClass,
                                keyClass >::minimalAllocationKeyedVector(
                             bool (*comparisonFunction)( keyClass const* const,
                                                       keyClass const* const ),
                             void (*constructionFunction)( storedClass** const,
                                                           keyClass** const,
                                                       keyClass const* const ),
                                     void (*resetFunction)( storedClass* const,
                                                            keyClass* const,
                                                    keyClass const* const ) ) :
  //currentlyStoredNumber( 0 ),
  comparisonFunction( comparisonFunction ),
  constructionFunction( constructionFunction ),
  resetFunction( resetFunction )
  {
    // just an initialization list.
  }

  template< class storedClass,
            class keyClass >
  inline
  minimalAllocationKeyedVector< storedClass,
                                keyClass >::~minimalAllocationKeyedVector()
  {
    {
      for( unsigned int deletionCounter( 0 );
           allConstructedPointers.size() > deletionCounter;
          ++deletionCounter )
      {
        delete allConstructedPointers.at( deletionCounter )->first;
        delete allConstructedPointers.at( deletionCounter )->second;
        delete allConstructedPointers.at( deletionCounter );
      }
      // I cannot use an iterator here because of template complications. also,
      // I'm going from 0 upwards because of unsigned int complications.
    }
  }


  template< class storedClass,
            class keyClass >
  inline void
  minimalAllocationKeyedVector< storedClass,
                                keyClass >::setComparison(
                             bool (*comparisonFunction)( keyClass const* const,
                                                      keyClass const* const ) )
  {
    this->comparisonFunction = comparisonFunction;
  }

  template< class storedClass,
            class keyClass >
  inline storedClass*
  minimalAllocationKeyedVector< storedClass,
                                keyClass >::getPointer(
                                             keyClass const* const keyPointer )
  /* this returns a pointer to the storedClass instance that matches
   * keyPointer, if there is a stored pointer to such an instance, or it
   * constructs a new storedClass instance with keyPointer & returns the
   * pointer to that.
   */
  {
    // set things up for the search:
    keyNotFound = true;
    searchCounter = 0;
    while( keyNotFound
           &&
           ( currentPointers.size() > searchCounter ) )
    {
      if( (*comparisonFunction)(
                            allConstructedPointers.at( searchCounter )->second,
                                 keyPointer ) )
      {
        returnPointer = allConstructedPointers.at( searchCounter )->first;
        keyNotFound = false;
      }
      ++searchCounter;
    }
    if( keyNotFound )
      // if we didn't find a match...
    {
      if( allConstructedPointers.size() > currentPointers.size() )
        // if we have stored-but-not-current storedClass instances...
      {
        (*resetFunction)(
                    allConstructedPointers.at( currentPointers.size() )->first,
                   allConstructedPointers.at( currentPointers.size() )->second,
                          keyPointer );
        // reset 1 of the stored-but-not-current storedClass instances with
        // keyPointer.
        returnPointer
        = allConstructedPointers.at( currentPointers.size() )->first;
        // note this pointer.
      }
      else
      {
        (*constructionFunction)( &returnPointer,
                                 &keyClassConstructionPointer,
                                 keyPointer );
        // construct a new storedClass instance & its associated keyClass
        // instance.
        allConstructedPointers.push_back(
                                     new std::pair< storedClass*,
                                                    keyClass* >( returnPointer,
                                               keyClassConstructionPointer ) );
        // hold the pointers in allConstructedPointers.
      }
      //++currentlyStoredNumber;
      currentPointers.push_back( returnPointer );
      // note that we have a new current storedClass instance.
    }
    return returnPointer;
  }

  template< class storedClass,
            class keyClass >
  inline std::vector< storedClass* >*
  minimalAllocationKeyedVector< storedClass,
                                keyClass >::getVector()
  // this returns a pointer to the std::vector holding pointers to the
  // storedClass instances which are deemed current.
  {
    return &currentPointers;
  }

  template< class storedClass,
            class keyClass >
  inline void
  minimalAllocationKeyedVector< storedClass,
                                keyClass >::clearEntries()
  // this notes that all current storedClass instances are now to be treated
  // as stored-but-not-current, & drops their pointers from currentPointers.
  {
    //currentlyStoredNumber = 0;
    currentPointers.clear();
  }

}  // end of LHC_FASER namespace.

#endif /* LHC_FASER_TEMPLATE_CLASSES_HPP_ */
