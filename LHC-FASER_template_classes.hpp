/*
 * LHC-FASER_template_classes.hpp
 *
 *  Created on: Jun 24, 2011
 *      Author: oleary
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

/* you might notice a change of code-writing style with this file. I (BOL) just
 * decided that I should be less idiosyncratic & adopt a conventional style, at
 * least to some extent...
 */

namespace LHC_FASER
{
  /* this is a template class to hold a vector of pointers to objects. the main
   * use is when the objects are slow to construct but fast to reset.
   * it needs to be used with storedClass objects that can be constructed by
   * a single function which takes no arguments.
   */
  template< class storedClass >
  class minimalAllocationVector
  {
  public:
    minimalAllocationVector( storedClass* const (*constructionFunction)() )
    /* code after the classes in this .hpp file. */;

    ~minimalAllocationVector()
    /* code after the classes in this .hpp file. */;

    storedClass*
    getPointer( int const requestedElement )
    // if an element out of range is asked for, this returns NULL.
    /* code after the classes in this .hpp file. */;

    storedClass*
    getFront()
    // this returns the 1st element deemed current. if there is none, this
    // returns NULL.
    /* code after the classes in this .hpp file. */;

    storedClass*
    getBack()
    // this returns the last element deemed current. if there is none, this
    // returns NULL.
    /* code after the classes in this .hpp file. */;

    int
    getSize()
    // this returns the number of elements which are deemed current.
    /* code after the classes in this .hpp file. */;

    std::vector< storedClass* >*
    getVector()
    // this returns a pointer to the std::vector holding pointers to the
    // storedClass instances which are deemed current.
    /* code after the classes in this .hpp file. */;

    void
    clear()
    // this notes that all current storedClass instances are now to be treated
    // as stored-but-not-current, & drops their pointers from currentPointers.
    /* code after the classes in this .hpp file. */;

    storedClass*
    addNewAtEnd()
    /* this, if necessary, pushes back a new storedClass instance into
     * allConstructedPointers, & pushes back a pointer into currentPointers, &
     * returns the pointer at the end ofcurrentPointers.
     */
    /* code after the classes in this .hpp file. */;


  protected:
    std::vector< storedClass* > allConstructedPointers;
    std::vector< storedClass* > currentPointers;
    storedClass* (*constructionFunction)();
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
    minimalAllocationKeyedVector( bool (*comparisonFunction)( keyClass* const,
                                                             keyClass* const ),
                             void (*constructionFunction)( storedClass** const,
                                                           keyClass** const,
                                                           keyClass* const ),
                                  void (*resetFunction)( storedClass* const,
                                                         keyClass* const,
                                                         keyClass* const ) )
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
    /* code after the classes in this .hpp file. */;

    ~minimalAllocationKeyedVector()
    /* code after the classes in this .hpp file. */;


    void
    setComparison( bool (*comparisonFunction)( keyClass* const,
                                               keyClass* const ) )
    /* code after the classes in this .hpp file. */;

    storedClass*
    getPointer( keyClass* const keyPointer )
    /* this returns a pointer to the storedClass instance that matches
     * keyPointer, if there is a stored pointer to such an instance, or it
     * constructs a new storedClass instance with keyPointer & returns the
     * pointer to that.
     */
    /* code after the classes in this .hpp file. */;

    std::vector< storedClass* >*
    getVector()
    // this returns a pointer to the std::vector holding pointers to the
    // storedClass instances which are deemed current.
    /* code after the classes in this .hpp file. */;

    void
    clear()
    // this notes that all current storedClass instances are now to be treated
    // as stored-but-not-current, & drops their pointers from currentPointers.
    /* code after the classes in this .hpp file. */;


  protected:
    std::vector< std::pair< storedClass*,
                            keyClass* >* > allConstructedPointers;
    std::vector< storedClass* > currentPointers;
    storedClass* returnPointer;
    keyClass* keyClassConstructionPointer;
    //unsigned int currentlyStoredNumber;
    unsigned int searchCounter;
    bool keyNotFound;
    bool (*comparisonFunction)( keyClass* const,
                                keyClass* const );
    void (*constructionFunction)( storedClass** const,
                                  keyClass** const,
                                  keyClass* const );
    void (*resetFunction)( storedClass* const,
                           keyClass* const,
                           keyClass* const );
  };





  // inline functions:

  template< class storedClass >
  inline
  minimalAllocationVector< storedClass >::minimalAllocationVector(
                               storedClass* const (*constructionFunction)() ) :
  constructionFunction( constructionFunction )
  {
    // just an initialization list.
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
  inline int
  minimalAllocationVector< storedClass >::getSize()
  // this returns the number of elements which are deemed current.
  {
    return (int)(currentPointers.size());
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
  inline void
  minimalAllocationVector< storedClass >::clear()
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
   * returns the pointer at the end ofcurrentPointers.
   */
  {
    while( currentPointers.size() >= allConstructedPointers.size() )
      // this should cover cases where currentPointers has been changed in size
      // by malevolent use of getVector...
    {
      allConstructedPointers.push_back( (*constructionFunction)() );
    }
    currentPointers.push_back(
                         allConstructedPointers.at( currentPointers.size() ) );
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
                                   bool (*comparisonFunction)( keyClass* const,
                                                             keyClass* const ),
                             void (*constructionFunction)( storedClass** const,
                                                           keyClass** const,
                                                           keyClass* const ),
                                     void (*resetFunction)( storedClass* const,
                                                            keyClass* const,
                                                          keyClass* const ) ) :
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
    for( unsigned int deletionCounter( 0 );
         allConstructedPointers.size() > deletionCounter;
         ++deletionCounter )
    {
      delete allConstructedPointers.at( deletionCounter )->first;
      delete allConstructedPointers.at( deletionCounter )->second;
      delete allConstructedPointers.at( deletionCounter );
    }
  }


  template< class storedClass,
            class keyClass >
  inline void
  minimalAllocationKeyedVector< storedClass,
                                keyClass >::setComparison(
                                   bool (*comparisonFunction)( keyClass* const,
                                                            keyClass* const ) )
  {
    this->comparisonFunction = comparisonFunction;
  }

  // this returns a pointer to the storedClass instance that matches
  // keyPointer, if there is a stored pointer to such an instance, or it
  // constructs a new storedClass instance with keyPointer & returns the
  // pointer to that.
  template< class storedClass,
            class keyClass >
  inline storedClass*
  minimalAllocationKeyedVector< storedClass,
                                keyClass >::getPointer(
                                                   keyClass* const keyPointer )
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

  // this returns a pointer to the std::vector holding pointers to the
  // storedClass instances which are deemed current.
  template< class storedClass,
            class keyClass >
  inline std::vector< storedClass* >*
  minimalAllocationKeyedVector< storedClass,
                                keyClass >::getVector()
  {
    return &currentPointers;
  }

  // this notes that all current storedClass instances are now to be treated
  // as stored-but-not-current, & drops their pointers from currentPointers.
  template< class storedClass,
            class keyClass >
  inline void
  minimalAllocationKeyedVector< storedClass,
                                keyClass >::clear()
  {
    //currentlyStoredNumber = 0;
    currentPointers.clear();
  }
}  // end of LHC_FASER namespace.

#endif /* LHC_FASER_TEMPLATE_CLASSES_HPP_ */
