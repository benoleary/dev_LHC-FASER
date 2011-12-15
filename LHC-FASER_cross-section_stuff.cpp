/*
 * LHC-FASER_cross-section_stuff.cpp
 *
 *  Created on: 15 Nov 2010
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

#include "LHC-FASER_cross-section_stuff.hpp"

namespace LHC_FASER
{
  squareGrid::squareGrid( std::string const* const gridFileLocation,
                          std::string const* const gridName,
                          squareGrid const* const scalingGrid,
                          inputHandler const* const inputShortcut ) :
    gridName( *gridName ),
    inputShortcut( inputShortcut )
  /* this constructor reads in a grid file, assumed to be in the format
   * x_coordinate y_coordinate value
   * in ascending order, y_coordinate varying first
   * (e.g.
   * 200.0 200.0 123.4
   * 200.0 210.0 567.8
   * ...
   * 210.0 200.0 876.5
   * 210.0 210.0 432.1
   * ...).
   * if provided with a non-NULL pointer to another squareGrid, it stores
   * the its acceptanceValues scaled by an interpolated value from this other
   * squareGrid.
   */
  {
    readIn( gridFileLocation,
            scalingGrid );

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: squareGrid() finished, value at ( 200.0, 200.0 ) = "
    << valueAt( 200.0,
                200.0 );
    std::cout << std::endl;**/
  }

  squareGrid::squareGrid( std::string baseGridFileLocation,
                          std::string gridName,
                          std::string scalingGridFileLocation,
                          inputHandler const* const inputShortcut ) :
    gridName( gridName ),
    inputShortcut( inputShortcut )
  // this constructor uses the other constructor to make another squareGrid
  // with the scaling factors, which is then used to construct this instance.
  {
    squareGrid* scalingGrid( NULL );
    if( !(scalingGridFileLocation.empty()) )
    {
      scalingGrid = new squareGrid( &scalingGridFileLocation,
                                    &gridName,
                                    NULL,
                                    inputShortcut );
    }
    readIn( &baseGridFileLocation,
            scalingGrid );
    if( NULL != scalingGrid )
    {
      delete scalingGrid;
    }
  }

  squareGrid::~squareGrid()
  {
    // delete all the vectors held as pointers by the vector of vectors:
    for( std::vector< std::vector< double >* >::iterator
         deletionIterator( gridValues.begin() );
         gridValues.end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
  }


  void
  squareGrid::readIn( std::string const* const gridFileLocation,
                      squareGrid const* const scalingGrid )
  // this actually does most of the job of the constructors.
  {
    // debugging:
    /**std::cout
    << std::endl
    << "debugging: squareGrid::readIn( "
    << *gridFileLocation << ", " << scalingGrid << " ) called";
    std::cout << std::endl;**/

    if( NULL != gridFileLocation )
    {
      // open the file:
      int const maxLengthOfGridFile( 100000 );
      CppSLHA::hash_commented_file_reader inputFileReader( gridFileLocation,
                                                           maxLengthOfGridFile,
                                                           false );
      /* a hundred thousand lines should be long enough for any grid file. if
       * readIn() is reading in more lines than that, it's probably because
       * of some file corruption.
       * also, if it's not enough, change it here.
       */
      std::string inputLineAsString;
      // prepare a string to parse the line.
      std::stringstream inputLineAsStream;
      // prepare a stringstream to parse the line.
      // the masses from the last line read & the currently-read line:
      double lastXCoordinate( 0.0 );
      double lastYCoordinate( 0.0 );
      double currentXCoordinate( 0.0 );
      double currentYCoordinate( 0.0 );
      double currentValue;
      // the value currently being read in.
      bool gridSizeStillUnknown( true );
      // we start out not knowing the step size for the grid (but we do
      // assume that it is a square grid).
      std::vector< double >* currentVector( NULL );
      // this holds all the acceptanceValues for a set of points with the same
      // x co-ordinate.
      unsigned int sizeOfLargestConstantXVector( 0 );
      // this holds the size of the largest vector that gets added to values.

      while( inputFileReader.read_line( &inputLineAsString ) )
        // the evaluation of the conditional reads in the next line.
      {
        // read in the line for interpretting:
        inputLineAsStream.clear();
        inputLineAsStream.str( inputLineAsString );
        // update what the masses for the last line were:
        lastYCoordinate = currentYCoordinate;
        lastXCoordinate = currentXCoordinate;
        inputLineAsStream >> currentXCoordinate
                          >> currentYCoordinate
                          >> currentValue;
        if( currentXCoordinate > lastXCoordinate )
          // if the x co-ordinate has changed, we need a new vector for the
          // acceptanceValues for varying y co-ordinates.
        {
          // 1st record the last set of points with the same
          // x co-ordinate:
          if( NULL != currentVector )
            // if we have a vector to add...
          {
            gridValues.push_back( currentVector );
            if( currentVector->size() > sizeOfLargestConstantXVector )
            {
              sizeOfLargestConstantXVector = currentVector->size();
            }
            // debugging:
            /**std::cout
            << std::endl
            << "pushed back a vector of size "
            << currentVector->size();**/
          }
          // prepare a new vector for the new x co-ordinate:
          currentVector = new std::vector< double >;
        }

        if( NULL != scalingGrid )
          // if we have to scale this value...
        {
          // debugging:
          /**std::cout
          << std::endl
          << "debugging: scaling point at ( "
          << currentXCoordinate << ", " << currentYCoordinate
          << " ) with value " << currentValue << " by factor "
          << scalingGrid->valueAt( currentXCoordinate,
                                   currentYCoordinate );
          std::cout << std::endl;**/

          currentValue *= scalingGrid->valueAt( currentXCoordinate,
                                                currentYCoordinate );
        }
        currentVector->push_back( currentValue );

        if( gridSizeStillUnknown
            &&
            ( 0.0 < lastXCoordinate )
            &&
            ( 0.0 < lastYCoordinate ) )
          // if we still need to work out the grid step size & we have read
          // in enough to do so...
        {
          // first we take note of the 1st point, which should have the
          // lowest co-ordinate acceptanceValues:
          lowestXCoordinate = lastXCoordinate;
          lowestYCoordinate = lastYCoordinate;
          // then we work out the grid step size:
          if( currentXCoordinate > lastXCoordinate )
          {
            gridStepSize = ( currentXCoordinate - lastXCoordinate );
            gridSizeStillUnknown = false;
          }
          else if( currentYCoordinate > lastYCoordinate )
          {
            gridStepSize = ( currentYCoordinate - lastYCoordinate );
            gridSizeStillUnknown = false;
          }
          // debugging:
          /**std::cout
          << std::endl
          << "gridStepSize = " << gridStepSize;**/
        }
      }  // end of while loop going over the lines of the file.

      // we still have to add the last vector of values for the last x
      // co-ordinate:
      gridValues.push_back( currentVector );
      if( currentVector->size() > sizeOfLargestConstantXVector )
      {
        sizeOfLargestConstantXVector = currentVector->size();
      }

      // debugging:
      /**std::cout
      << std::endl
      << "pushed back last vector of size "
      << current_vector->size();
      std::cout << std::endl;**/

      // the last co-ordinates should be the maximal acceptanceValues:
      highestXCoordinate = currentXCoordinate;
      highestYCoordinate = currentYCoordinate;

      // now we check that the grid is usable:
      for( std::vector< std::vector< double >* >::iterator
           xIterator( gridValues.begin() );
           gridValues.end() > xIterator;
           ++xIterator )
      {
        while( (*xIterator)->size() < sizeOfLargestConstantXVector )
        {
          (*xIterator)->push_back(
                                 CppSLHA::CppSLHA_global::really_wrong_value );
          std::cout
          << std::endl
          << "LHC-FASER::warning! A cross-section grid (\""
          << *gridFileLocation << "\") was malformed. Inserting "
          << CppSLHA::CppSLHA_global::really_wrong_value_string
          << " at each missing point. Interpolations near such a point should"
          << " also end up with this value.";
          std::cout << std::endl;
        }
      }

    }
    else
      // otherwise, a NULL pointer was given instead of a pointer to the
      // grid file name string.
    {
      std::cout
      << std::endl
      << "LHC-FASER::error! asked to open a file with a NULL pointer!";
      std::cout << std::endl;
    }
    // debugging:
    /**std::cout
    << std::endl
    << "debugging: squareGrid::readIn() finished: "
    << std::endl << "lowestXCoordinate = " << lowestXCoordinate
    << std::endl << "highestXCoordinate = " << highestXCoordinate
    << std::endl << "lowestYCoordinate = " << lowestYCoordinate
    << std::endl << "highestYCoordinate = " << highestYCoordinate
    << std::endl << "value at ( " << lowestXCoordinate << ", "
    << lowestYCoordinate << " ) is " << valueAt( lowestXCoordinate,
                                                 lowestYCoordinate )
    << std::endl << "value at ( "
    << ( 0.5 * ( lowestXCoordinate + highestXCoordinate ) + 0.2 ) << ", "
    << ( 0.5 * ( lowestYCoordinate + highestYCoordinate ) + 0.3 )
    << " ) is "
    << valueAt( ( 0.5 * ( lowestXCoordinate + highestXCoordinate ) + 0.2 ),
                ( 0.5 * ( lowestXCoordinate + highestXCoordinate ) + 0.3 ) );
    std::cout << std::endl;**/
  }

  double
  squareGrid::valueAt( double const xCoordinate,
                       double const yCoordinate )
  const
  /* this finds the grid square which the given point is in, & then uses
   * lhcFaserGlobal::squareBilinearInterpolation to get an interpolated
   * value. if the requested co-ordinates are outside the grid, or on the
   * edges of the grid corresponding to the maximum x or y acceptanceValues,
   * CppSLHA::CppSLHA_global::really_wrong_value is returned.
   * N.B.: cross-section grids use
   * x_coordinate == SQUARK mass, y_coordinate == GLUINO mass!
   */
  {
    if( ( 0.0 < gridStepSize )
        &&
        ( xCoordinate >= lowestXCoordinate )
        &&
        ( yCoordinate >= lowestYCoordinate ) )
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "gridValues.size() = " << gridValues.size()
      << std::endl
      << "gridValues.back()->size() = "
      << gridValues.back()->size();
      std::cout << std::endl;**/

      double xSteps( ( ( xCoordinate - lowestXCoordinate ) / gridStepSize ) );
      unsigned int lowerLeftX( (unsigned int)xSteps );
      double
      ySteps( ( ( yCoordinate - lowestYCoordinate ) / gridStepSize ) );
      unsigned int lowerLeftY( (unsigned int)ySteps );
      if( gridValues.size() > ( lowerLeftX + 1 ) )
        // if the x co-ordinate is less than its maximal grid value...
      {
        if( gridValues.at( lowerLeftX )->size() > ( lowerLeftY + 1 ) )
          // if the y co-ordinate is less than its maximal grid value...
        {
          return lhcFaserGlobal::squareBilinearInterpolation(
                                                       ( xSteps - lowerLeftX ),
                                                       ( ySteps - lowerLeftY ),
                                 gridValues.at( lowerLeftX )->at( lowerLeftY ),
                         gridValues.at( ( lowerLeftX + 1 ) )->at( lowerLeftY ),
                 gridValues.at( ( lowerLeftX + 1 ) )->at( ( lowerLeftY + 1 ) ),
                       gridValues.at( lowerLeftX )->at( ( lowerLeftY + 1 ) ) );
        }
        else if( ( gridValues.at( lowerLeftX )->size() == ( lowerLeftY + 1 ) )
                 &&
                 ( (double)lowerLeftY == ySteps ) )
          // otherwise, if it's on the maximal y edge...
        {
          return lhcFaserGlobal::unitLinearInterpolation(
                                                       ( xSteps - lowerLeftX ),
                                 gridValues.at( lowerLeftX )->at( lowerLeftY ),
                       gridValues.at( ( lowerLeftX + 1 ) )->at( lowerLeftY ) );
        }
        else
          // otherwise, it's off the grid:
        {
          //return CppSLHA::CppSLHA_global::really_wrong_value;
          // actually, let's do some extrapolation:
          if( inputShortcut->isVerbose() )
          {
            std::cout
            << std::endl
            << "LHC-FASER::warning! A cross-section for masses ("
            << xCoordinate << ", " << yCoordinate << ") outside the grid was"
            << " requested. Returning an exponential extrapolation.";
            std::cout << std::endl;
          }
          double
          lowerLeftValue( gridValues.at( lowerLeftX )->at( gridValues.at(
                                                  lowerLeftX )->size() - 2 ) );
          double upperLeftValue( gridValues.at( lowerLeftX )->back() );
          double
          lowerRightValue( gridValues.at( ( lowerLeftX + 1 ) )->at(
                                              gridValues.at( ( lowerLeftX + 1 )
                                                             )->size() - 2 ) );
          double
          upperRightValue( gridValues.at( ( lowerLeftX + 1 ) )->back() );
          double leftRatio( upperLeftValue / lowerLeftValue );
          if( 1.0 < leftRatio )
          {
            leftRatio = 1.0;
          }
          double rightRatio( upperRightValue / lowerRightValue );
          if( 1.0 < rightRatio )
          {
            rightRatio = 1.0;
          }
          double const yFraction( ySteps - lowerLeftY );
          // debugging:
          /**std::cout << std::endl << "debugging:"
          << std::endl
          << "lowerLeftX = " << lowerLeftX
          << "; off the grid in y. lowerLeftY = " << lowerLeftY
          << std::endl
          << "lowerLeftValue = " << lowerLeftValue
          << ", upperLeftValue = " << upperLeftValue
          << ", lowerRightValue = " << lowerRightValue
          << ", upperRightValue = " << upperRightValue
          << ", leftRatio = " << leftRatio
          << ", rightRatio = " << rightRatio;
          std::cout << std::endl;**/
          while( gridValues.at( lowerLeftX )->size() <= ( lowerLeftY + 1 ) )
          {
            lowerLeftValue *= leftRatio;
            upperLeftValue *= leftRatio;
            lowerRightValue *= rightRatio;
            upperRightValue *= rightRatio;
            --lowerLeftY;
            // note that we've extrapolated the grid a step.
          }
          /* at this point, the nearest grid square has been transformed into
           * the square of an exponential extrapolation of the grid for the
           * requested point.
           */
          // debugging:
          /**std::cout << std::endl << "debugging:"
          << std::endl
          << "after extrapolating, xFraction = " << ( xSteps - lowerLeftX )
          << ", yFraction = " << yFraction
          << ", lowerLeftValue = " << lowerLeftValue
          << ", upperLeftValue = " << upperLeftValue
          << ", lowerRightValue = " << lowerRightValue
          << ", upperRightValue = " << upperRightValue;
          std::cout << std::endl;**/
          return lhcFaserGlobal::squareBilinearInterpolation(
                                                       ( xSteps - lowerLeftX ),
                                                              yFraction,
                                                              lowerLeftValue,
                                                              lowerRightValue,
                                                              upperRightValue,
                                                              upperLeftValue );
        }
      }
      else if( ( gridValues.size() == ( lowerLeftX + 1 ) )
               &&
               ( (double)lowerLeftX == xSteps ) )
        // otherwise, if it's on the maximal x edge...
      {
        if( gridValues.back()->size() > ( lowerLeftY + 1 ) )
          // if the y co-ordinate is less than its maximal grid value...
        {
          return lhcFaserGlobal::unitLinearInterpolation(
                                                       ( ySteps - lowerLeftY ),
                                           gridValues.back()->at( lowerLeftY ),
                                 gridValues.back()->at( ( lowerLeftY + 1 ) ) );
        }
        else if( ( gridValues.back()->size() == ( lowerLeftY + 1 ) )
                 &&
                 ( (double)lowerLeftY == ySteps ) )
          // otherwise, if it's on the maximal x & y corner...
        {
          return gridValues.back()->back();
        }
        else
          // otherwise, it's off the grid:
        {
          //return CppSLHA::CppSLHA_global::really_wrong_value;
          // actually, let's do some extrapolation:
          if( inputShortcut->isVerbose() )
          {
            std::cout
            << std::endl
            << "LHC-FASER::warning! A cross-section for masses ("
            << xCoordinate << ", " << yCoordinate << ") outside the grid was"
            << " requested. Returning an exponential extrapolation.";
            std::cout << std::endl;
          }
          double
          lowerValue( gridValues.back()->at( gridValues.back()->size() - 2 ) );
          double upperValue( gridValues.back()->back() );
          double extrapolationRatio( upperValue / lowerValue );
          if( 1.0 < extrapolationRatio )
          {
            extrapolationRatio = 1.0;
          }
          double const yFraction( ySteps - lowerLeftY );
          // debugging:
          /**std::cout << std::endl << "debugging:"
          << std::endl
          << "lowerLeftX = " << lowerLeftX
          << "; off the grid in y. lowerLeftY = " << lowerLeftY
          << std::endl
          << "lowerValue = " << lowerValue
          << ", upperValue = " << upperValue
          << ", extrapolationRatio = " << extrapolationRatio;
          std::cout << std::endl;**/
          while( gridValues.back()->size() <= ( lowerLeftY + 1 ) )
          {
            lowerValue *= extrapolationRatio;
            upperValue *= extrapolationRatio;
            --lowerLeftY;
            // note that we've extrapolated the grid a step.
          }
          /* at this point, the nearest grid square has been transformed into
           * the square of an exponential extrapolation of the grid for the
           * requested point.
           */
          // debugging:
          /**std::cout << std::endl << "debugging:"
          << std::endl
          << "after extrapolating, yFraction = " << yFraction
          << ", lowerValue = " << lowerValue
          << ", upperValue = " << upperValue;
          std::cout << std::endl;**/
          return lhcFaserGlobal::unitLinearInterpolation( yFraction,
                                                          lowerValue,
                                                          upperValue );
        }
      }
      else
        // otherwise, it's off the grid:
      {
        //return CppSLHA::CppSLHA_global::really_wrong_value;
        // actually, let's do some extrapolation:
        if( inputShortcut->isVerbose() )
        {
          std::cout
          << std::endl
          << "LHC-FASER::warning! A cross-section for masses ("
          << xCoordinate << ", " << yCoordinate << ") outside the grid was"
          << " requested. Returning an exponential extrapolation.";
          std::cout << std::endl;
          // debugging:
          /**std::cout << std::endl << "debugging:"
          << std::endl
          << "lowerLeftX = " << lowerLeftX
          << ", lowerLeftY = " << lowerLeftY;
          std::cout << std::endl;**/
        }
        if( ( gridValues.back()->size() > lowerLeftY )
            &&
            ( (double)lowerLeftY == ySteps ) )
          // if we're in line with a row of grid points...
        {
          double leftValue( gridValues.at( gridValues.size()
                                           - 2 )->at( lowerLeftY ) );
          double rightValue( gridValues.back()->at( lowerLeftY ) );
          double extrapolationRatio( rightValue / leftValue );
          if( 1.0 < extrapolationRatio )
          {
            extrapolationRatio = 1.0;
          }
          double xFraction( xSteps - lowerLeftX );
          // debugging:
          /**std::cout << std::endl << "debugging:"
          << std::endl
          << "lowerLeftX = " << lowerLeftX
          << "; off the grid in y. lowerLeftY = " << lowerLeftY
          << std::endl
          << "leftValue = " << leftValue
          << ", rightValue = " << rightValue
          << ", extrapolationRatio = " << extrapolationRatio;
          std::cout << std::endl;**/
          while( gridValues.size() <= ( lowerLeftX + 1 ) )
          {
            leftValue *= extrapolationRatio;
            rightValue *= extrapolationRatio;
            --lowerLeftX;
            // note that we've extrapolated the grid a step.
          }
          /* at this point, the nearest grid square has been transformed into
           * the square of an exponential extrapolation of the grid for the
           * requested point.
           */
          // debugging:
          /**std::cout << std::endl << "debugging:"
          << std::endl
          << "after extrapolating, xFraction = " << xFraction
          << ", leftValue = " << leftValue
          << ", rightValue = " << rightValue;
          std::cout << std::endl;**/
          return lhcFaserGlobal::unitLinearInterpolation( xFraction,
                                                          leftValue,
                                                          rightValue );

        }
        else if( gridValues.back()->size() > ( lowerLeftY + 1 ) )
          // if the y co-ordinate is less than its maximal grid value...
        {
          double
          lowerLeftValue( gridValues.at( gridValues.size()
                                         - 2 )->at( lowerLeftY ) );
          double
          upperLeftValue( gridValues.at( gridValues.size() - 2 )->at(
                                                            lowerLeftY + 1 ) );
          double lowerRightValue( gridValues.back()->at( lowerLeftY ) );
          double upperRightValue( gridValues.back()->at( lowerLeftY + 1 ) );
          double lowerRatio( lowerRightValue / lowerLeftValue );
          if( 1.0 < lowerRatio )
          {
            lowerRatio = 1.0;
          }
          double upperRatio( upperRightValue / upperLeftValue );
          if( 1.0 < upperRatio )
          {
            upperRatio = 1.0;
          }
          double const xFraction( xSteps - lowerLeftX );
          // debugging:
          /**std::cout << std::endl << "debugging:"
          << std::endl
          << "lowerLeftY = " << lowerLeftY
          << "; off the grid in x. lowerLeftX = " << lowerLeftX
          << std::endl
          << "lowerLeftValue = " << lowerLeftValue
          << ", upperLeftValue = " << upperLeftValue
          << ", lowerRightValue = " << lowerRightValue
          << ", upperRightValue = " << upperRightValue
          << ", lowerRatio = " << lowerRatio
          << ", upperRatio = " << upperRatio;
          std::cout << std::endl;**/
          while( gridValues.size() <= ( lowerLeftX + 1 ) )
          {
            lowerLeftValue *= lowerRatio;
            lowerRightValue *= lowerRatio;
            upperLeftValue *= upperRatio;
            upperRightValue *= upperRatio;
            --lowerLeftX;
            // note that we've extrapolated the grid a step.
          }
          /* at this point, the nearest grid square has been transformed into
           * the square of an exponential extrapolation of the grid for the
           * requested point.
           */
          // debugging:
          /**std::cout << std::endl << "debugging:"
          << std::endl
          << "after extrapolating, xFraction = " << xFraction
          << ", yFraction = " << ( ySteps - lowerLeftY )
          << ", lowerLeftValue = " << lowerLeftValue
          << ", upperLeftValue = " << upperLeftValue
          << ", lowerRightValue = " << lowerRightValue
          << ", upperRightValue = " << upperRightValue;
          std::cout << std::endl;**/
          return lhcFaserGlobal::squareBilinearInterpolation( xFraction,
                                                       ( ySteps - lowerLeftY ),
                                                              lowerLeftValue,
                                                              lowerRightValue,
                                                              upperRightValue,
                                                              upperLeftValue );
        }
        else
          // otherwise both the x & y co-ordinates are outside the grid ranges
          // & we have to do more extraploation:
        {
          double
          lowerLeftValue( gridValues.at( gridValues.size() - 2 )->at(
                                           gridValues.at( gridValues.size() - 2
                                                             )->size() - 2 ) );
          double
          upperLeftValue( gridValues.at( gridValues.size() - 2 )->back() );
          double
          lowerRightValue( gridValues.back()->at( gridValues.back()->size()
                                                  - 2 ) );
          double upperRightValue( gridValues.back()->back() );
          double lowerRatio( lowerRightValue / lowerLeftValue );
          if( 1.0 < lowerRatio )
          {
            lowerRatio = 1.0;
          }
          double upperRatio( upperRightValue / upperLeftValue );
          if( 1.0 < upperRatio )
          {
            upperRatio = 1.0;
          }
          double leftRatio( upperLeftValue / lowerLeftValue );
          if( 1.0 < leftRatio )
          {
            leftRatio = 1.0;
          }
          double rightRatio( upperRightValue / lowerRightValue );
          if( 1.0 < rightRatio )
          {
            rightRatio = 1.0;
          }
          double const xFraction( xSteps - lowerLeftX );
          double const yFraction( ySteps - lowerLeftY );
          // debugging:
          /**std::cout << std::endl << "debugging:"
          << std::endl
          << "lowerLeftYY = " << lowerLeftY
          << ", lowerLeftX = " << lowerLeftX
          << "; off the grid in both x & y."
          << std::endl
          << "lowerLeftValue = " << lowerLeftValue
          << ", upperLeftValue = " << upperLeftValue
          << ", lowerRightValue = " << lowerRightValue
          << ", upperRightValue = " << upperRightValue
          << ", lowerRatio = " << lowerRatio
          << ", upperRatio = " << upperRatio
          << ", leftRatio = " << leftRatio
          << ", rightRatio = " << rightRatio;
          std::cout << std::endl;**/
          while( gridValues.back()->size() <= ( lowerLeftY + 1 ) )
          {
            lowerLeftValue *= leftRatio;
            upperLeftValue *= leftRatio;
            lowerRightValue *= rightRatio;
            upperRightValue *= rightRatio;
            --lowerLeftY;
            // note that we've extrapolated the grid a step.
          }
          while( gridValues.size() <= ( lowerLeftX + 1 ) )
          {
            lowerLeftValue *= lowerRatio;
            lowerRightValue *= lowerRatio;
            upperLeftValue *= upperRatio;
            upperRightValue *= upperRatio;
            --lowerLeftX;
            // note that we've extrapolated the grid a step.
          }
          /* at this point, the nearest grid square has been transformed into
           * the square of an exponential extrapolation of the grid for the
           * requested point.
           */
          // debugging:
          /**std::cout << std::endl << "debugging:"
          << std::endl
          << "after extrapolating, xFraction = " << xFraction
          << ", yFraction = " << yFraction
          << ", lowerLeftValue = " << lowerLeftValue
          << ", upperLeftValue = " << upperLeftValue
          << ", lowerRightValue = " << lowerRightValue
          << ", upperRightValue = " << upperRightValue;
          std::cout << std::endl;**/
          return lhcFaserGlobal::squareBilinearInterpolation( xFraction,
                                                              yFraction,
                                                              lowerLeftValue,
                                                              lowerRightValue,
                                                              upperRightValue,
                                                              upperLeftValue );
        }
      }
    }
    else
      // otherwise, it's off the grid (or the grid step size was not positive,
      // which is still a problem...):
    {
      return CppSLHA::CppSLHA_global::really_wrong_value;
      // we don't extrapolate in this case because we'd be extrapolating to
      // lower sparticle masses, & that's not particularly sensible.
    }
  }



  crossSectionTable::crossSectionTable( squareGrid const* const lookupTable,
                       signedParticleShortcutPair const* const scoloredPair,
                                        double const flavorFactor,
                                    inputHandler const* const inputShortcut ) :
    getsReadiedForNewPoint( inputShortcut->getReadier() ),
    inputShortcut( inputShortcut ),
    scoloredPair( scoloredPair ),
    lookupTable( lookupTable ),
    flavorFactor( flavorFactor )
  {
    // debugging:
    /**std::cout
    << std::endl
    << "debugging: crossSectionTable( "
    << given_lookup_table << ", { "
    << given_pair->getFirstParticle()->getName() << ": "
    << given_pair->is_first_particle() << "; "
    << given_pair->getSecondParticle()->getName() << ": "
    << given_pair->is_second_particle() << " }, "
    << given_flavor_factor << ", " << given_inputShortcuts << " ) called";
    std::cout << std::endl;**/

    // we take the 1st squark found to ask for its mass in looking up the
    // cross-section. if no squark was provided, we keep the pointer as NULL.
    if( inputShortcut->isIn( scoloredPair->getFirstParticle()->get_PDG_code(),
                         inputShortcut->getSquarks() ) )
    {
      directlyProducedSquark = scoloredPair->getFirstParticle();
    }
    else if( inputShortcut->isIn(
                            scoloredPair->getSecondParticle()->get_PDG_code(),
                              inputShortcut->getSquarks() ) )
    {
      directlyProducedSquark = scoloredPair->getSecondParticle();
    }
    else
    {
      directlyProducedSquark = NULL;
    }
  }

  crossSectionTable::~crossSectionTable()
  {
    // does nothing.
  }


  double
  crossSectionTable::getValue()
  /* this checks to see if it needs to recalculate the value for this point, &
   * if so, it does, & if the point lies outside the grid, it takes the value
   * of the nearest point on the grid instead.
   */
  {
    // debugging:
    /**std::cout
    << std::endl
    << "debugging: crossSectionTable::getValue() called."
    << " needsToPrepare = "
    << needsToPrepare << ", "
    << "flavorFactor = " << flavorFactor;
    std::cout << std::endl;**/

    if( needsToPrepareForThisPoint() )
    {
      if( 0.0 < flavorFactor )
      {
        double squarkMassToUse;
        if( NULL != directlyProducedSquark )
          // if we have a directly-produced squark...
        {
          squarkMassToUse = directlyProducedSquark->get_absolute_mass();
        }
        else
          // otherwise, we use the average squark mass...
        {
          squarkMassToUse = inputShortcut->getAverageSquarks4Mass();
        }
        if( squarkMassToUse < lookupTable->getLowestX() )
        {
          if( inputShortcut->isVerbose() )
          {
            std::cout
            << std::endl
            << "LHC-FASER::warning! cross-section requested for a squark"
            << " mass ( " << squarkMassToUse << " ) lower than the"
            << " lowest squark mass of the lookup table (";
            std::cout
            << lookupTable->getLowestX()
            << "); using this lowest mass instead.";
            std::cout << std::endl;
          }
          squarkMassToUse = lookupTable->getLowestX();
        }
        if( inputShortcut->getGluinoMass() < lookupTable->getLowestY() )
        {
          if( inputShortcut->isVerbose() )
          {
            std::cout
            << std::endl
            << "LHC-FASER::warning! cross-section requested for a gluino"
            << " mass ( " << inputShortcut->getGluinoMass()
            << " ) lower than the lowest gluino mass of the lookup table ( ";
            std::cout
            << lookupTable->getLowestY()
            << " ); using this lowest mass instead.";
            std::cout << std::endl;
          }
          storedValue = lookupTable->valueAt( squarkMassToUse,
                                              lookupTable->getLowestY() );
        }
        else
        {
          storedValue = lookupTable->valueAt( squarkMassToUse,
                                              inputShortcut->getGluinoMass() );
        }
        storedValue *= flavorFactor;
      }
      else
      {
        storedValue = 0.0;
      }
    }
    return storedValue;
  }



  crossSectionTableSet::crossSectionTableSet(
                                        std::string const* const gridDirectory,
                                              int const beamEnergy,
                                    inputHandler const* const inputShortcut ) :
    inputShortcut( inputShortcut ),
    gridDirectory( *gridDirectory ),
    beamEnergy( beamEnergy )
  {
    // just an initialization list.
  }

  crossSectionTableSet::~crossSectionTableSet()
  {
    for( std::vector< squareGrid* >::iterator
         deletionIterator( gridSet.begin() );
         gridSet.end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
    for( std::vector< crossSectionTable* >::iterator
         deletionIterator( tableSet.begin() );
         tableSet.end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
  }


  double
  crossSectionTableSet::prepareGridName( std::string* const gridName,
                         signedParticleShortcutPair const* const scoloredPair )
  // this returns the flavor factor for the requested pair while putting the
  // squareGrid's associated string in gridName.
  {
    /* here we make the approximations about how the cross-sections are in
     * squark flavors.
     *
     * gluino+gluino is assumed to depend on the average squark mass.
     *
     * gluino+squark is assumed to be non-zero only for sdown & sup,
     * & that the ratio of sup to sdown is 2:1 (twice as many ups in the
     * other proton for the gluon from its proton to fuse with as downs).
     *
     * squark+antisquark is assumed to be non-zero only for both squarks
     * being the same flavor. all the _L+_L/_R+_R (_one+_one/_two+_two) entries
     * are assumed to be equal for the different flavors. for squark+antisquark
     * pairs of different chiralities, only sdown & sup have non-zero entries &
     * sdown_L+antisdown_R:sdown_L+antisup_R:sup_L+antisdown_R:sup_L+antisup_L
     * is assumed to be 1:1:2:2 & likewise for _L <-> _R (twice as many ups in
     * the other proton for the antiquark from its proton to fuse with as
     * downs).
     *
     * squark+squark is assumed to be non-zero only for both sdown, both
     * sup, or sdown+sup, & only for both being the same chirality (LO sQCD t-
     * & u-channels are identically zero if the squarks are opposite
     * chirality).
     * the ratios for sdown+sdown : sdown+sup : sup+sup are assumed to be
     * 1 * 8/3 * 1/2 : 4 * 1 * 1 : 4 * 8/3 * 1/2
     * (ratios of getting the various valence quarks from both
     * protons * matrix-element-squared in the heavy gluino mass limit
     * (t-channel + u-channel + interference terms with 1/N_c suppression, * 2
     * for t*u + u*t) * phase space halving for identical final-state
     * particles) = 4:12:16 = 1:3:4
     */
    if( ( CppSLHA::PDG_code::gluino
          == scoloredPair->getFirstParticle()->get_PDG_code() )
        &&
        ( CppSLHA::PDG_code::gluino
          == scoloredPair->getSecondParticle()->get_PDG_code() ) )
    {
      gridName->assign( "gluino+gluino" );
      return 1.0;
    }
    else if( inputShortcut->isIn(
                              scoloredPair->getFirstParticle()->get_PDG_code(),
                                  inputShortcut->getSquarks() )
             &&
             ( CppSLHA::PDG_code::gluino
               == scoloredPair->getSecondParticle()->get_PDG_code() ) )
    {
      gridName->assign( "sdown_L+gluino" );
      if( CppSLHA::PDG_code::sdown_L
          == scoloredPair->getFirstParticle()->get_PDG_code() )
      {
        return 1.0;
      }
      else if( CppSLHA::PDG_code::sdown_R
               == scoloredPair->getFirstParticle()->get_PDG_code() )
      {
        return 1.0;
      }
      else if( CppSLHA::PDG_code::sup_L
               == scoloredPair->getFirstParticle()->get_PDG_code() )
      {
        return 2.0;
      }
      else if( CppSLHA::PDG_code::sup_R
               ==  scoloredPair->getFirstParticle()->get_PDG_code() )
      {
        return 2.0;
      }
      else
      {
        return 0.0;
      }
    }
    else if( inputShortcut->isIn(
                             scoloredPair->getSecondParticle()->get_PDG_code(),
                                  inputShortcut->getSquarks() )
             &&
             ( CppSLHA::PDG_code::gluino
               == scoloredPair->getFirstParticle()->get_PDG_code() ) )
    {
      gridName->assign( "sdown_L+gluino" );
      if( CppSLHA::PDG_code::sdown_L
          == scoloredPair->getSecondParticle()->get_PDG_code() )
      {
        return 1.0;
      }
      else if( CppSLHA::PDG_code::sdown_R
               == scoloredPair->getSecondParticle()->get_PDG_code() )
      {
        return 1.0;
      }
      else if( CppSLHA::PDG_code::sup_L
               == scoloredPair->getSecondParticle()->get_PDG_code() )
      {
        return 2.0;
      }
      else if( CppSLHA::PDG_code::sup_R
               ==  scoloredPair->getSecondParticle()->get_PDG_code() )
      {
        return 2.0;
      }
      else
      {
        return 0.0;
      }
    }
    else if( inputShortcut->isIn(
                              scoloredPair->getFirstParticle()->get_PDG_code(),
                                  inputShortcut->getSquarks() )
             &&
             inputShortcut->isIn(
                             scoloredPair->getSecondParticle()->get_PDG_code(),
                                  inputShortcut->getSquarks() ) )
    {
      if( ( scoloredPair->firstIsNotAntiparticle()
          &&
          !scoloredPair->secondIsNotAntiparticle() )
          ||
          ( !scoloredPair->firstIsNotAntiparticle()
              &&
              scoloredPair->secondIsNotAntiparticle() ) )
      {
        if( ( CppSLHA::PDG_code::sdown_L
              == scoloredPair->getFirstParticle()->get_PDG_code() )
            &&
            ( CppSLHA::PDG_code::sdown_L
              == scoloredPair->getSecondParticle()->get_PDG_code() ) )
        {
          gridName->assign( "sdown_L+antisdown_L" );
          return 1.0;
        }
        else if( ( CppSLHA::PDG_code::sdown_R
                   == scoloredPair->getFirstParticle()->get_PDG_code() )
                 &&
                 ( CppSLHA::PDG_code::sdown_R
                   == scoloredPair->getSecondParticle()->get_PDG_code() ) )
        {
          gridName->assign( "sdown_L+antisdown_L" );
          return 1.0;
        }
        else if( ( CppSLHA::PDG_code::sdown_L
                   == scoloredPair->getFirstParticle()->get_PDG_code() )
                 &&
                 ( CppSLHA::PDG_code::sdown_R
                   == scoloredPair->getSecondParticle()->get_PDG_code() ) )
        {
          gridName->assign( "sdown_L+antisdown_R" );
          return 1.0;
        }
        else if( ( CppSLHA::PDG_code::sdown_R
                   == scoloredPair->getFirstParticle()->get_PDG_code() )
                 &&
                 ( CppSLHA::PDG_code::sdown_L
                   == scoloredPair->getSecondParticle()->get_PDG_code() ) )
        {
          gridName->assign( "sdown_L+antisdown_R" );
          return 1.0;
        }
        else if( ( CppSLHA::PDG_code::sdown_L
                   == scoloredPair->getFirstParticle()->get_PDG_code() )
                 &&
                 ( CppSLHA::PDG_code::sup_R
                   == scoloredPair->getSecondParticle()->get_PDG_code() ) )
        {
          gridName->assign( "sdown_L+antisdown_R" );
          return 1.0;
        }
        else if( ( CppSLHA::PDG_code::sdown_R
                   == scoloredPair->getFirstParticle()->get_PDG_code() )
                 &&
                 ( CppSLHA::PDG_code::sup_L
                   == scoloredPair->getSecondParticle()->get_PDG_code() ) )
        {
          gridName->assign( "sdown_L+antisdown_R" );
          return 1.0;
        }
        else if( ( CppSLHA::PDG_code::sup_L
                   == scoloredPair->getFirstParticle()->get_PDG_code() )
                 &&
                 ( CppSLHA::PDG_code::sdown_R
                   == scoloredPair->getSecondParticle()->get_PDG_code() ) )
        {
          gridName->assign( "sdown_L+antisdown_R" );
          return 2.0;
        }
        else if( ( CppSLHA::PDG_code::sup_R
                   == scoloredPair->getFirstParticle()->get_PDG_code() )
                 &&
                 ( CppSLHA::PDG_code::sdown_L
                   == scoloredPair->getSecondParticle()->get_PDG_code() ) )
        {
          gridName->assign( "sdown_L+antisdown_R" );
          return 2.0;
        }
        else if( ( CppSLHA::PDG_code::sup_L
                   == scoloredPair->getFirstParticle()->get_PDG_code() )
                 &&
                 ( CppSLHA::PDG_code::sup_L
                   == scoloredPair->getSecondParticle()->get_PDG_code() ) )
        {
          gridName->assign( "sdown_L+antisdown_L" );
          return 1.0;
        }
        else if( ( CppSLHA::PDG_code::sup_R
                   == scoloredPair->getFirstParticle()->get_PDG_code() )
                 &&
                 ( CppSLHA::PDG_code::sup_R
                   == scoloredPair->getSecondParticle()->get_PDG_code() ) )
        {
          gridName->assign( "sdown_L+antisdown_L" );
          return 1.0;
        }
        else if( ( CppSLHA::PDG_code::sup_L
                   == scoloredPair->getFirstParticle()->get_PDG_code() )
                 &&
                 ( CppSLHA::PDG_code::sup_R
                   == scoloredPair->getSecondParticle()->get_PDG_code() ) )
        {
          gridName->assign( "sdown_L+antisdown_R" );
          return 2.0;
        }
        else if( ( CppSLHA::PDG_code::sup_R
                   == scoloredPair->getFirstParticle()->get_PDG_code() )
                 &&
                 ( CppSLHA::PDG_code::sup_L
                   == scoloredPair->getSecondParticle()->get_PDG_code() ) )
        {
          gridName->assign( "sdown_L+antisdown_R" );
          return 2.0;
        }
        else if( ( CppSLHA::PDG_code::sstrange_L
                   == scoloredPair->getFirstParticle()->get_PDG_code() )
                 &&
                 ( CppSLHA::PDG_code::sstrange_L
                   == scoloredPair->getSecondParticle()->get_PDG_code() ) )
        {
          gridName->assign( "sdown_L+antisdown_L" );
          return 1.0;
        }
        else if( ( CppSLHA::PDG_code::sstrange_R
                   == scoloredPair->getFirstParticle()->get_PDG_code() )
                 &&
                 ( CppSLHA::PDG_code::sstrange_R
                   == scoloredPair->getSecondParticle()->get_PDG_code() ) )
        {
          gridName->assign( "sdown_L+antisdown_L" );
          return 1.0;
        }
        else if( ( CppSLHA::PDG_code::scharm_L
                   == scoloredPair->getFirstParticle()->get_PDG_code() )
                 &&
                 ( CppSLHA::PDG_code::scharm_L
                   == scoloredPair->getSecondParticle()->get_PDG_code() ) )
        {
          gridName->assign( "sdown_L+antisdown_L" );
          return 1.0;
        }
        else if( ( CppSLHA::PDG_code::scharm_R
                   == scoloredPair->getFirstParticle()->get_PDG_code() )
                 &&
                 ( CppSLHA::PDG_code::scharm_R
                   == scoloredPair->getSecondParticle()->get_PDG_code() ) )
        {
          gridName->assign( "sdown_L+antisdown_L" );
          return 1.0;
        }
        else if( ( CppSLHA::PDG_code::sbottom_one
                   == scoloredPair->getFirstParticle()->get_PDG_code() )
                 &&
                 ( CppSLHA::PDG_code::sbottom_one
                   == scoloredPair->getSecondParticle()->get_PDG_code() ) )
        {
          gridName->assign( "sdown_L+antisdown_L" );
          return 1.0;
        }
        else if( ( CppSLHA::PDG_code::sbottom_two
                   == scoloredPair->getFirstParticle()->get_PDG_code() )
                 &&
                 ( CppSLHA::PDG_code::sbottom_two
                   == scoloredPair->getSecondParticle()->get_PDG_code() ) )
        {
          gridName->assign( "sdown_L+antisdown_L" );
          return 1.0;
        }
        else if( ( CppSLHA::PDG_code::stop_one
                   == scoloredPair->getFirstParticle()->get_PDG_code() )
                 &&
                 ( CppSLHA::PDG_code::stop_one
                   == scoloredPair->getSecondParticle()->get_PDG_code() ) )
        {
          gridName->assign( "sdown_L+antisdown_L" );
          return 1.0;
        }
        else if( ( CppSLHA::PDG_code::stop_two
                   == scoloredPair->getFirstParticle()->get_PDG_code() )
                 &&
                 ( CppSLHA::PDG_code::stop_two
                   == scoloredPair->getSecondParticle()->get_PDG_code() ) )
        {
          gridName->assign( "sdown_L+antisdown_L" );
          return 1.0;
        }
        else
        {
          return 0.0;
        }
      }
      else if( scoloredPair->firstIsNotAntiparticle()
               &&
               scoloredPair->secondIsNotAntiparticle() )
      {
        if( ( CppSLHA::PDG_code::sdown_L
              == scoloredPair->getFirstParticle()->get_PDG_code() )
            &&
            ( CppSLHA::PDG_code::sdown_L
              == scoloredPair->getSecondParticle()->get_PDG_code() ) )
        {
          gridName->assign( "sdown_L+sdown_L" );
          return 1.0;
        }
        else if( ( CppSLHA::PDG_code::sdown_R
                   == scoloredPair->getFirstParticle()->get_PDG_code() )
                 &&
                 ( CppSLHA::PDG_code::sdown_R
                   == scoloredPair->getSecondParticle()->get_PDG_code() ) )
        {
          gridName->assign( "sdown_L+sdown_L" );
          return 1.0;
        }
        else if( ( CppSLHA::PDG_code::sdown_L
                   == scoloredPair->getFirstParticle()->get_PDG_code() )
                 &&
                 ( CppSLHA::PDG_code::sup_L
                   == scoloredPair->getSecondParticle()->get_PDG_code() ) )
        {
          gridName->assign( "sdown_L+sdown_L" );
          return 3.0;
        }
        else if( ( CppSLHA::PDG_code::sdown_R
                   == scoloredPair->getFirstParticle()->get_PDG_code() )
                 &&
                 ( CppSLHA::PDG_code::sup_R
                   == scoloredPair->getSecondParticle()->get_PDG_code() ) )
        {
          gridName->assign( "sdown_L+sdown_L" );
          return 3.0;
        }
        else if( ( CppSLHA::PDG_code::sup_L
                   == scoloredPair->getFirstParticle()->get_PDG_code() )
                 &&
                 ( CppSLHA::PDG_code::sdown_L
                   == scoloredPair->getSecondParticle()->get_PDG_code() ) )
        {
          gridName->assign( "sdown_L+sdown_L" );
          return 3.0;
        }
        else if( ( CppSLHA::PDG_code::sup_R
                   == scoloredPair->getFirstParticle()->get_PDG_code() )
                 &&
                 ( CppSLHA::PDG_code::sdown_R
                   == scoloredPair->getSecondParticle()->get_PDG_code() ) )
        {
          gridName->assign( "sdown_L+sdown_L" );
          return 3.0;
        }
        else if( ( CppSLHA::PDG_code::sup_L
                   == scoloredPair->getFirstParticle()->get_PDG_code() )
                 &&
                 ( CppSLHA::PDG_code::sup_L
                   == scoloredPair->getSecondParticle()->get_PDG_code() ) )
        {
          gridName->assign( "sdown_L+sdown_L" );
          return 4.0;
        }
        else if( ( CppSLHA::PDG_code::sup_R
                   == scoloredPair->getFirstParticle()->get_PDG_code() )
                 &&
                 ( CppSLHA::PDG_code::sup_R
                   == scoloredPair->getSecondParticle()->get_PDG_code() ) )
        {
          gridName->assign( "sdown_L+sdown_L" );
          return 4.0;
        }
        else
        {
          return 0.0;
        }
      }
      else
      {
        return 0.0;
      }
    }
    else
    {
      return 0.0;
    }
  }


  squareGrid const*
  crossSectionTableSet::getGrid( std::string const* const gridName )
  // this either finds the appropriate squareGrid in grids or makes a new
  // instance, stores it, & returns its pointer.
  {
    // debugging:
    /**std::cout
    << std::endl
    << "debugging: crossSectionTableSet::getGrid( "
    << *gridName << " ) called";
    std::cout << std::endl;**/

    squareGrid* returnPointer( NULL );
    // this starts as NULL so that we know if it wasn't found among the
    // existing instances.
    for( std::vector< squareGrid* >::iterator gridIterator( gridSet.begin() );
         gridSet.end() > gridIterator;
         ++gridIterator )
    {
      if( 0 == (*gridIterator)->getName()->compare( *gridName ) )
        // if a table with this name already exists...
      {
        returnPointer = *gridIterator;
        // note the table.
        gridIterator = gridSet.end();
        // stop looking.
      }
    }

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: returnPointer = " << returnPointer;
    std::cout << std::endl;**/

    if( NULL == returnPointer )
      // if we didn't find an existing match, we make a new instance, store
      // it & return a pointer to it:
    {
      std::string loGridFileName( gridDirectory );
      loGridFileName.append( "/" );
      loGridFileName.append( *gridName );
      loGridFileName.append( "_LO.dat" );
      // debugging:
      /**std::cout
      << std::endl
      << "debugging: loGridFileName = " << loGridFileName;
      std::cout << std::endl;**/

      // we check for whether we should use NLO K-factors:
      squareGrid* nloGrid( NULL );
      if( inputShortcut->usingNlo() )
      {
        std::string nloGridFileName( gridDirectory );
        nloGridFileName.append( "/" );
        nloGridFileName.append( *gridName );
        nloGridFileName.append( "_K.dat" );

        // debugging:
        /**std::cout
        << std::endl
        << "debugging: nloGridFileName = " << nloGridFileName;
        std::cout << std::endl;**/

        nloGrid = new squareGrid( &nloGridFileName,
                                  gridName,
                                  NULL,
                                  inputShortcut );
        // the grid's name isn't important because it is only temporary.

        // debugging:
        /**std::cout
        << std::endl
        << "debugging: nloGrid = " << nloGrid
        << ", valueAt( 200.0, 200.0 ) = "
        << nloGrid->valueAt( 200.0,
                             200.0 );
        std::cout << std::endl;**/
      }
      // we make a new table:
      returnPointer = new squareGrid( &loGridFileName,
                                      gridName,
                                      nloGrid,
                                      inputShortcut );
      if( NULL != nloGrid )
        // if we needed a K-factor grid...
      {
        // we can let it go now:
        delete nloGrid;
      }
      // debugging:
      /**std::cout
        << std::endl
        << "debugging: return_pointer = " << return_pointer;
        std::cout << std::endl;**/

      gridSet.push_back( returnPointer );
    }
    // debugging:
    /**std::cout
    << std::endl
    << "debugging: returning = " << return_pointer;
    std::cout << std::endl;**/

    return returnPointer;
  }


  crossSectionTable*
  crossSectionTableSet::getTable(
                     signedParticleShortcutPair const* const requestedChannel )
  // this returns the crossSectionTable for the requested pair.
  {
    // debugging:
    /**std::cout
    << std::endl
    << "debugging: crossSectionTableSet::getTable( { "
    << *(requested_channel->getFirstParticle()->get_name_or_antiname(
                              requested_channel->firstIsNotAntiparticle() ))
    << "; "
    << *(requested_channel->getSecondParticle()->get_name_or_antiname(
                             requested_channel->secondIsNotAntiparticle() ))
    << " } ) called";
    std::cout << std::endl;**/

    crossSectionTable* returnPointer( NULL );
    // this starts as NULL so that we know if it wasn't found among the
    // existing instances.

    for( std::vector< crossSectionTable* >::iterator
         tableIterator( tableSet.begin() );
         tableSet.end() > tableIterator;
         ++tableIterator )
    {
      if( requestedChannel == (*tableIterator)->getPair() )
        // if we find the requested pair...
      {
        // debugging:
        /**std::cout
        << std::endl
        << "found existing table.";
        std::cout << std::endl;**/

        returnPointer = *tableIterator;
        // note the table to return.
        tableIterator = tableSet.end();
        // stop looking.
      }
    }
    // debugging:
    /**std::cout
    << std::endl
    << "debugging: returnPointer = " << returnPointer;
    std::cout << std::endl;**/

    if( NULL == returnPointer )
      // if we didn't find an existing match, we make a new instance, store
      // it & return a pointer to it:
    {
      // debugging:
      /**std::cout
      << std::endl
      << "did not find existing table.";
      std::cout << std::endl;**/

      // 1st we need to find out which squareGrid to use & what the flavor
      // factor is:
      std::string gridName;
      double flavorFactor( prepareGridName( &gridName,
                                            requestedChannel ) );

      // debugging:
      /**std::cout
      << std::endl
      << "debugging: gridName = " << gridName
      << std::endl
      << "debugging: flavorFactor = " << flavorFactor;
      std::cout << std::endl;**/
      if( 0.0 < flavorFactor )
        // if it is a non-zero channel as far as this code is concerned...
      {
        returnPointer = new crossSectionTable( getGrid( &gridName ),
                                               requestedChannel,
                                               flavorFactor,
                                               inputShortcut );
      }
      else
      {
        returnPointer = new crossSectionTable( NULL,
                                               requestedChannel,
                                               0.0,
                                               inputShortcut );
        /* a crossSectionTable initialized with a NULL squareGrid
         * pointer *should* be OK if its flavor factor is 0.0, so that it
         * never calls any of the squareGrid's member functions.
         */
      }
      tableSet.push_back( returnPointer );
    }
    return returnPointer;
  }




  crossSectionHandler::crossSectionHandler(
                                    inputHandler const* const inputShortcut ) :
    inputShortcut( inputShortcut )
  {
    // just an initialization list
  }

  crossSectionHandler::~crossSectionHandler()
  {
    for( std::vector< crossSectionTableSet* >::iterator
         deletionIterator( tableSets.begin() );
         tableSets.end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
  }


  crossSectionTableSet*
  crossSectionHandler::getTableSet( int const beamEnergyInTev )
  {
    crossSectionTableSet* returnPointer( NULL );
    // this starts as NULL so that we know if it wasn't found among the
    // existing instances.
    for( std::vector< crossSectionTableSet* >::iterator
         setIterator( tableSets.begin() );
         tableSets.end() > setIterator;
         ++setIterator )
    {
      if( beamEnergyInTev == (*setIterator)->getBeamEnergy() )
      {
        returnPointer = *setIterator;
        setIterator = tableSets.end();
        // stop looking.
      }
    }

    if( NULL == returnPointer )
      // if we didn't find an existing match, we make a new instance, store
      // it & return a pointer to it:
    {
      std::string
      energyName( *(inputShortcut->inspectPathToCrossSectionGrids()) );
      std::stringstream energyStream( "" );
      energyStream << "/" << beamEnergyInTev << "TeV";
      energyName.append( energyStream.str() );
      returnPointer = new crossSectionTableSet( &energyName,
                                                beamEnergyInTev,
                                                inputShortcut );
      tableSets.push_back( returnPointer );
    }
    return returnPointer;
  }

}  // end of LHC_FASER namespace.
