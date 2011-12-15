/*
 * LHC-FASER_kinematics_stuff.hpp
 *
 *  Created on: 02 Mar 2011
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

#include "LHC-FASER_jet_kinematics_stuff.hpp"

namespace LHC_FASER
{
  squarkMassForGridDecider::squarkMassForGridDecider(
                                    inputHandler const* const inputShortcut ) :
    inputShortcut( inputShortcut )
  {
    // just an initialization list.
  }

  squarkMassForGridDecider::~squarkMassForGridDecider()
  {
    // does nothing.
  }



  heavierThanGluinoSquarkMassForGrid::heavierThanGluinoSquarkMassForGrid(
                                    inputHandler const* const inputShortcut ) :
    squarkMassForGridDecider( inputShortcut )
  {
    // just an initialization list.
  }

  heavierThanGluinoSquarkMassForGrid::~heavierThanGluinoSquarkMassForGrid()
  {
    // does nothing.
  }



  firstMassForGrid::firstMassForGrid(
                                    inputHandler const* const inputShortcut ) :
    squarkMassForGridDecider( inputShortcut )
  {
    // just an initialization list.
  }

  firstMassForGrid::~firstMassForGrid()
  {
    // does nothing.
  }



  secondMassForGrid::secondMassForGrid(
                                    inputHandler const* const inputShortcut ) :
    squarkMassForGridDecider( inputShortcut )
  {
    // just an initialization list.
  }

  secondMassForGrid::~secondMassForGrid()
  {
    // does nothing.
  }



  averageSquarkMassForGrid::averageSquarkMassForGrid(
                                    inputHandler const* const inputShortcut ) :
    squarkMassForGridDecider( inputShortcut )
  {
    // just an initialization list.
  }

  averageSquarkMassForGrid::~averageSquarkMassForGrid()
  {
    // does nothing.
  }



  int const jetAcceptanceGrid::columnsNotInAcceptanceVector( 5 );
  // we want to skip the squark, gluino, & electroweakino masses, & we also
  // want to count from 1 rather than from 0.

  jetAcceptanceGrid::jetAcceptanceGrid(
                                     std::string const* const gridFileLocation,
                                    inputHandler const* const inputShortcut ) :
      acceptanceGrid( gridFileLocation,
                      inputShortcut )
  {
    // just an initialization list.
  }

  jetAcceptanceGrid::~jetAcceptanceGrid()
  {
    // does nothing.
  }


  double
  jetAcceptanceGrid::interpolateOnElectroweakinoAndAcceptanceIndices(
                                                 int const electroweakinoIndex,
                                                    int const acceptanceIndex )
  /* this checks to see if the appropriate entries in lowerLeftVectorOfVectors,
   * lowerRightVectorOfVectors, upperRightVectorOfVectors, &
   * upperLeftVectorOfVectors exist, & if so, interpolates on the squark &
   * gluino masses for these indices.
   */
  {
    if( ( 0 <= electroweakinoIndex )
        &&
        ( 0 <= acceptanceIndex )
        &&
        ( lowerLeftVectorOfVectors->size()
          > (unsigned int)electroweakinoIndex )
        &&
        ( lowerLeftVectorOfVectors->at( electroweakinoIndex )->size()
          > (unsigned int)acceptanceIndex )
        &&
        ( lowerRightVectorOfVectors->size()
          > (unsigned int)electroweakinoIndex )
          &&
        ( lowerRightVectorOfVectors->at( electroweakinoIndex )->size()
          > (unsigned int)acceptanceIndex )
        &&
        ( upperRightVectorOfVectors->size()
          > (unsigned int)electroweakinoIndex )
          &&
        ( upperRightVectorOfVectors->at( electroweakinoIndex )->size()
          > (unsigned int)acceptanceIndex )
        &&
        ( upperLeftVectorOfVectors->size()
          > (unsigned int)electroweakinoIndex )
        &&
        ( upperLeftVectorOfVectors->at( electroweakinoIndex )->size()
          > (unsigned int)acceptanceIndex ) )
    {
      return lhcFaserGlobal::squareBilinearInterpolation( squarkMassFraction,
                                                          gluinoMassFraction,
                                                  lowerLeftVectorOfVectors->at(
                                                     electroweakinoIndex )->at(
                                                             acceptanceIndex ),
                                                 lowerRightVectorOfVectors->at(
                                                     electroweakinoIndex )->at(
                                                             acceptanceIndex ),
                                                 upperRightVectorOfVectors->at(
                                                     electroweakinoIndex )->at(
                                                             acceptanceIndex ),
                                                  upperLeftVectorOfVectors->at(
                                                     electroweakinoIndex )->at(
                                                           acceptanceIndex ) );
    }
    else
    {
      shouldInterpolateOnElectroweakinos = false;
      // if the lower value was uncalculable, we shouldn't bother interpolating
      // with an upper value, even if it is calculable.
      return CppSLHA::CppSLHA_global::really_wrong_value;
    }
  }

  double
  jetAcceptanceGrid::getAcceptance( double const squarkMass,
                                    double const gluinoMass,
                                    double firstElectroweakinoMass,
                                    double secondElectroweakinoMass,
                                    int acceptanceColumn )
  // this interpolates the grid to obtain values for the acceptances based on
  // the given masses.
  {
    // 1st we ensure that firstElectroweakinoMass <= secondElectroweakinoMass:
    if( firstElectroweakinoMass > secondElectroweakinoMass )
    {
      lowerScoloredMass = firstElectroweakinoMass;
      firstElectroweakinoMass = secondElectroweakinoMass;
      secondElectroweakinoMass = lowerScoloredMass;
    }
    // then we find out which scolored mass is lower:
    if( squarkMass > gluinoMass )
    {
      lowerScoloredMass = gluinoMass;
    }
    else
    {
      lowerScoloredMass = squarkMass;
    }
    if( ( secondElectroweakinoMass < lowerScoloredMass )
        &&
        pointIsOnGrid( squarkMass,
                       gluinoMass ) )
      // if the point is on the grid...
    {
      acceptanceColumn -= columnsNotInAcceptanceVector;
      // now acceptanceColumn is the index of the element in the vector that
      // we want to find.
      if( secondElectroweakinoMass
          <= ( lowerScoloredMass * lowElectroweakinoMassRatio ) )
        // if both electroweakinos are very light, we use just the lightest
        // electroweakino grid point (interpolated on the scolored masses):
      {
        returnValue = interpolateOnElectroweakinoAndAcceptanceIndices(
                                              indexForLightLightNeutralinoPair,
                                                            acceptanceColumn );
        // this sets returnValue appropriately.

        // debugging:
        /**std::cout << std::endl << "debugging:"
        << std::endl
        << "below-left of light-light, returning " << returnValue;
        std::cout << std::endl;**/
      }
      else
        // otherwise, we need to interpolate...
      {
        shouldInterpolateOnElectroweakinos = true;
        /* shouldInterpolateOnElectroweakinos is set to true _before_ calling
         * interpolateOnElectroweakinoAndAcceptanceIndices because that
         * function might set it to false because it was unable to calculate
         * a lower value.
         */

        if( firstElectroweakinoMass
            <= ( lowerScoloredMass * lowElectroweakinoMassRatio ) )
          // if we're going to do a linear interpolation...
        {
          // debugging:
          /**std::cout << std::endl << "debugging:"
          << std::endl
          << "left of light lighter line.";
          std::cout << std::endl;**/

          if( secondElectroweakinoMass
              <= ( lowerScoloredMass * mediumElectroweakinoMassRatio ) )
            /* at this point we are sure that
             * secondElectroweakinoMass
             * > ( lowerScoloredMass * lowElectroweakinoMassRatio )
             */
          {
            // debugging:
            /**std::cout << std::endl << "debugging:"
            << std::endl
            << "below medium heavier line.";
            std::cout << std::endl;**/

            lighterLighterElectroweakinoPointValue
            = interpolateOnElectroweakinoAndAcceptanceIndices(
                                              indexForLightLightNeutralinoPair,
                                                            acceptanceColumn );
            if( shouldInterpolateOnElectroweakinos )
              // if that interpolation was successful...
            {
              lighterHeavierElectroweakinoPointValue
              = interpolateOnElectroweakinoAndAcceptanceIndices(
                                             indexForLightMediumNeutralinoPair,
                                                            acceptanceColumn );
              if( shouldInterpolateOnElectroweakinos )
                // if that interpolation was successful...
              {
                heavierElectroweakinoMassFraction
                = ( ( ( secondElectroweakinoMass / lowerScoloredMass )
                      - lowElectroweakinoMassRatio )
                    / ( mediumElectroweakinoMassRatio
                        - lowElectroweakinoMassRatio ) );
              }
            }
          }
          else if( secondElectroweakinoMass
                   <= ( lowerScoloredMass * highElectroweakinoMassRatio ) )
              /* at this point we are sure that
               * secondElectroweakinoMass
               * > ( lowerScoloredMass * mediumElectroweakinoMassRatio )
               */
          {
            // debugging:
            /**std::cout << std::endl << "debugging:"
            << std::endl
            << "below heavy heavier line.";
            std::cout << std::endl;**/
            lighterLighterElectroweakinoPointValue
            = interpolateOnElectroweakinoAndAcceptanceIndices(
                                             indexForLightMediumNeutralinoPair,
                                                            acceptanceColumn );
            if( shouldInterpolateOnElectroweakinos )
              // if that interpolation was successful...
            {
              lighterHeavierElectroweakinoPointValue
              = interpolateOnElectroweakinoAndAcceptanceIndices(
                                              indexForLightHeavyNeutralinoPair,
                                                            acceptanceColumn );
              if( shouldInterpolateOnElectroweakinos )
                // if that interpolation was successful...
              {
                heavierElectroweakinoMassFraction
                = ( ( ( secondElectroweakinoMass / lowerScoloredMass )
                      - mediumElectroweakinoMassRatio )
                   / ( highElectroweakinoMassRatio
                       - mediumElectroweakinoMassRatio ) );
              }
            }
          }
          else
          {
            // debugging:
            /**std::cout << std::endl << "debugging:"
            << std::endl
            << "above heavy heavier line.";
            std::cout << std::endl;**/

            lighterLighterElectroweakinoPointValue
            = interpolateOnElectroweakinoAndAcceptanceIndices(
                                              indexForLightHeavyNeutralinoPair,
                                                            acceptanceColumn );
            if( shouldInterpolateOnElectroweakinos )
              // if that interpolation was successful...
            {
              lighterHeavierElectroweakinoPointValue = 0.0;
              heavierElectroweakinoMassFraction
              = ( ( ( secondElectroweakinoMass / lowerScoloredMass )
                    - highElectroweakinoMassRatio )
                  / ( 1.0 - highElectroweakinoMassRatio ) );
            }
          }
          if( shouldInterpolateOnElectroweakinos )
            // if the interpolations so far were successful...
          {
            returnValue = lhcFaserGlobal::unitLinearInterpolation(
                                             heavierElectroweakinoMassFraction,
                                        lighterLighterElectroweakinoPointValue,
                                      lighterHeavierElectroweakinoPointValue );

            // debugging:
            /**std::cout << std::endl << "debugging:"
            << std::endl
            << "linear interpolation, from "
            << lighterLighterElectroweakinoPointValue << " to "
            << lighterHeavierElectroweakinoPointValue << " with fraction "
            << heavierElectroweakinoMassFraction
            << ", returning " << returnValue;
            std::cout << std::endl;**/
          }
          else
          {
            returnValue = CppSLHA::CppSLHA_global::really_wrong_value;
          }
        }
        else
          // otherwise we're going to do a bilinear interpolation:
        {
          // debugging:
          /**std::cout << std::endl << "debugging:"
          << std::endl
          << "right of light lighter line.";
          std::cout << std::endl;**/

          if( secondElectroweakinoMass
              <= ( lowerScoloredMass * mediumElectroweakinoMassRatio ) )
            /* at this point we are sure that
             * secondElectroweakinoMass
             * > ( lowerScoloredMass * lowElectroweakinoMassRatio ) & that
             * firstElectroweakinoMass
             *  > ( lowerScoloredMass * lowElectroweakinoMassRatio ) &
             *  <= ( lowerScoloredMass * mediumElectroweakinoMassRatio ) too.
             */
          {
            // debugging:
            /**std::cout << std::endl << "debugging:"
            << std::endl
            << "below medium heavier line.";
            std::cout << std::endl;**/

            lighterLighterElectroweakinoPointValue
            = interpolateOnElectroweakinoAndAcceptanceIndices(
                                              indexForLightLightNeutralinoPair,
                                                            acceptanceColumn );
            if( shouldInterpolateOnElectroweakinos )
              // if that interpolation was successful...
            {
              lighterHeavierElectroweakinoPointValue
              = interpolateOnElectroweakinoAndAcceptanceIndices(
                                             indexForLightMediumNeutralinoPair,
                                                            acceptanceColumn );
              if( shouldInterpolateOnElectroweakinos )
                // if that interpolation was successful...
              {
                heavierHeavierElectroweakinoPointValue
                = interpolateOnElectroweakinoAndAcceptanceIndices(
                                            indexForMediumMediumNeutralinoPair,
                                                            acceptanceColumn );
                if( shouldInterpolateOnElectroweakinos )
                  // if that interpolation was successful...
                {
                  heavierLighterElectroweakinoPointValue
                  = ( lighterLighterElectroweakinoPointValue
                      + heavierHeavierElectroweakinoPointValue
                      - lighterHeavierElectroweakinoPointValue );
                  lighterElectroweakinoMassFraction
                  = ( ( ( firstElectroweakinoMass / lowerScoloredMass )
                        - lowElectroweakinoMassRatio )
                      / ( mediumElectroweakinoMassRatio
                          - lowElectroweakinoMassRatio ) );
                  heavierElectroweakinoMassFraction
                  = ( ( ( secondElectroweakinoMass / lowerScoloredMass )
                        - lowElectroweakinoMassRatio )
                      / ( mediumElectroweakinoMassRatio
                          - lowElectroweakinoMassRatio ) );
                }
              }
            }
            /* end of checking for successful interpolations. at this point,
             * shouldInterpolateOnElectroweakinos will be true unless 1 of the
             * interpolations went wrong, & all the values for the final
             * bilinear interpolation will also be appropriately set.
             */
          }
          else if( secondElectroweakinoMass
                   <= ( lowerScoloredMass * highElectroweakinoMassRatio ) )
            /* at this point we are sure that
             * secondElectroweakinoMass
             * > ( lowerScoloredMass * mediumElectroweakinoMassRatio ) & that
             * firstElectroweakinoMass
             *  > ( lowerScoloredMass * lowElectroweakinoMassRatio ) too.
             */
          {
            // debugging:
            /**std::cout << std::endl << "debugging:"
            << std::endl
            << "below heavy heavier line.";
            std::cout << std::endl;**/

            if( firstElectroweakinoMass
                <= ( lowerScoloredMass * mediumElectroweakinoMassRatio ) )
              // if the heavier electroweakino is between the middle & high
              // masses & the lighter is between low & middle...
            {
              // debugging:
              /**std::cout << std::endl << "debugging:"
              << std::endl
              << "left of medium lighter line.";
              std::cout << std::endl;**/

              lighterLighterElectroweakinoPointValue
              = interpolateOnElectroweakinoAndAcceptanceIndices(
                                             indexForLightMediumNeutralinoPair,
                                                            acceptanceColumn );
              if( shouldInterpolateOnElectroweakinos )
                // if that interpolation was successful...
              {
                lighterHeavierElectroweakinoPointValue
                = interpolateOnElectroweakinoAndAcceptanceIndices(
                                              indexForLightHeavyNeutralinoPair,
                                                            acceptanceColumn );
                if( shouldInterpolateOnElectroweakinos )
                  // if that interpolation was successful...
                {
                  heavierHeavierElectroweakinoPointValue
                  = interpolateOnElectroweakinoAndAcceptanceIndices(
                                             indexForMediumHeavyNeutralinoPair,
                                                            acceptanceColumn );
                  if( shouldInterpolateOnElectroweakinos )
                    // if that interpolation was successful...
                  {
                    heavierLighterElectroweakinoPointValue
                    = interpolateOnElectroweakinoAndAcceptanceIndices(
                                            indexForMediumMediumNeutralinoPair,
                                                            acceptanceColumn );
                    if( shouldInterpolateOnElectroweakinos )
                      // if that interpolation was successful...
                    {
                      lighterElectroweakinoMassFraction
                      = ( ( ( firstElectroweakinoMass / lowerScoloredMass )
                            - lowElectroweakinoMassRatio )
                          / ( mediumElectroweakinoMassRatio
                              - lowElectroweakinoMassRatio ) );
                      heavierElectroweakinoMassFraction
                      = ( ( ( secondElectroweakinoMass / lowerScoloredMass )
                            - mediumElectroweakinoMassRatio )
                          / ( highElectroweakinoMassRatio
                              - mediumElectroweakinoMassRatio ) );
                    }
                  }
                }
              }
              /* end of checking for successful interpolations. at this point,
               * shouldInterpolateOnElectroweakinos will be true unless 1 of the
               * interpolations went wrong, & all the values for the final
               * bilinear interpolation will also be appropriately set.
               */
            }
            else
              // otherwise the heavier electroweakino is between the middle &
              // high masses & the lighter is between middle & high too...
            {
              // debugging:
              /**std::cout << std::endl << "debugging:"
              << std::endl
              << "right of medium lighter line.";
              std::cout << std::endl;**/

              lighterLighterElectroweakinoPointValue
              = interpolateOnElectroweakinoAndAcceptanceIndices(
                                            indexForMediumMediumNeutralinoPair,
                                                            acceptanceColumn );
              if( shouldInterpolateOnElectroweakinos )
                // if that interpolation was successful...
              {
                lighterHeavierElectroweakinoPointValue
                = interpolateOnElectroweakinoAndAcceptanceIndices(
                                             indexForMediumHeavyNeutralinoPair,
                                                            acceptanceColumn );
                if( shouldInterpolateOnElectroweakinos )
                  // if that interpolation was successful...
                {
                  heavierHeavierElectroweakinoPointValue
                  = interpolateOnElectroweakinoAndAcceptanceIndices(
                                              indexForHeavyHeavyNeutralinoPair,
                                                            acceptanceColumn );
                  if( shouldInterpolateOnElectroweakinos )
                    // if that interpolation was successful...
                  {
                    heavierLighterElectroweakinoPointValue
                    = ( lighterLighterElectroweakinoPointValue
                        + heavierHeavierElectroweakinoPointValue
                        - lighterHeavierElectroweakinoPointValue );
                    lighterElectroweakinoMassFraction
                    = ( ( ( firstElectroweakinoMass / lowerScoloredMass )
                          - mediumElectroweakinoMassRatio )
                        / ( highElectroweakinoMassRatio
                            - mediumElectroweakinoMassRatio ) );
                    heavierElectroweakinoMassFraction
                    = ( ( ( secondElectroweakinoMass / lowerScoloredMass )
                          - mediumElectroweakinoMassRatio )
                        / ( highElectroweakinoMassRatio
                            - mediumElectroweakinoMassRatio ) );
                  }
                }
              }
              /* end of checking for successful interpolations. at this point,
               * shouldInterpolateOnElectroweakinos will be true unless 1 of the
               * interpolations went wrong, & all the values for the final
               * bilinear interpolation will also be appropriately set.
               */
            }
          }
          else
            // otherwise we're in the region where the heavier electroweakino
            // is heavier than the heavier grid value...
          {
            // debugging:
            /**std::cout << std::endl << "debugging:"
            << std::endl
            << "above heavy heavier line.";
            std::cout << std::endl;**/

            if( firstElectroweakinoMass
                <= ( lowerScoloredMass * mediumElectroweakinoMassRatio ) )
              // if the lighter electroweakino is between the low & middle
              // masses...
            {
              // debugging:
              /**std::cout << std::endl << "debugging:"
              << std::endl
              << "left of medium lighter line.";
              std::cout << std::endl;**/

              lighterLighterElectroweakinoPointValue
              = interpolateOnElectroweakinoAndAcceptanceIndices(
                                              indexForLightHeavyNeutralinoPair,
                                                            acceptanceColumn );
              if( shouldInterpolateOnElectroweakinos )
                // if that interpolation was successful...
              {
                heavierLighterElectroweakinoPointValue
                = interpolateOnElectroweakinoAndAcceptanceIndices(
                                             indexForMediumHeavyNeutralinoPair,
                                                            acceptanceColumn );
              if( shouldInterpolateOnElectroweakinos )
                // if that interpolation was successful...
              {
                lighterHeavierElectroweakinoPointValue = 0.0;
                heavierHeavierElectroweakinoPointValue = 0.0;
                lighterElectroweakinoMassFraction
                = ( ( ( firstElectroweakinoMass / lowerScoloredMass )
                      - lowElectroweakinoMassRatio )
                    / ( mediumElectroweakinoMassRatio
                        - lowElectroweakinoMassRatio ) );
                heavierElectroweakinoMassFraction
                = ( ( ( secondElectroweakinoMass / lowerScoloredMass )
                      - highElectroweakinoMassRatio )
                    / ( 1.0 - highElectroweakinoMassRatio ) );
                }
              }
              /* end of checking for successful interpolations. at this point,
               * shouldInterpolateOnElectroweakinos will be true unless 1 of the
               * interpolations went wrong, & all the values for the final
               * bilinear interpolation will also be appropriately set.
               */
            }
            else if( firstElectroweakinoMass
                     <= ( lowerScoloredMass * highElectroweakinoMassRatio ) )
              // if the lighter electroweakino is between the middle & high
              // masses...
            {
              // debugging:
              /**std::cout << std::endl << "debugging:"
              << std::endl
              << "left of heavy lighter line.";
              std::cout << std::endl;**/

              lighterLighterElectroweakinoPointValue
              = interpolateOnElectroweakinoAndAcceptanceIndices(
                                             indexForMediumHeavyNeutralinoPair,
                                                            acceptanceColumn );
              if( shouldInterpolateOnElectroweakinos )
                // if that interpolation was successful...
              {
                heavierLighterElectroweakinoPointValue
                = interpolateOnElectroweakinoAndAcceptanceIndices(
                                              indexForHeavyHeavyNeutralinoPair,
                                                            acceptanceColumn );
                if( shouldInterpolateOnElectroweakinos )
                  // if that interpolation was successful...
                {
                  lighterHeavierElectroweakinoPointValue = 0.0;
                  heavierHeavierElectroweakinoPointValue = 0.0;
                  lighterElectroweakinoMassFraction
                  = ( ( ( firstElectroweakinoMass / lowerScoloredMass )
                        - mediumElectroweakinoMassRatio )
                      / ( highElectroweakinoMassRatio
                          - mediumElectroweakinoMassRatio ) );
                  heavierElectroweakinoMassFraction
                  = ( ( ( secondElectroweakinoMass / lowerScoloredMass )
                        - highElectroweakinoMassRatio )
                      / ( 1.0 - highElectroweakinoMassRatio ) );
                }
              }
              /* end of checking for successful interpolations. at this point,
               * shouldInterpolateOnElectroweakinos will be true unless 1 of the
               * interpolations went wrong, & all the values for the final
               * bilinear interpolation will also be appropriately set.
               */
            }
            else
              // otherwise both are in the region over the heavier grid value:
            {
              // debugging:
              /**std::cout << std::endl << "debugging:"
              << std::endl
              << "right of heavy lighter line.";
              std::cout << std::endl;**/

              lighterLighterElectroweakinoPointValue
              = interpolateOnElectroweakinoAndAcceptanceIndices(
                                              indexForHeavyHeavyNeutralinoPair,
                                                            acceptanceColumn );
              if( shouldInterpolateOnElectroweakinos )
                // if that interpolation was successful...
              {
                heavierLighterElectroweakinoPointValue = 0.0;
                lighterHeavierElectroweakinoPointValue = 0.0;
                heavierHeavierElectroweakinoPointValue = 0.0;
                lighterElectroweakinoMassFraction
                = ( ( ( firstElectroweakinoMass / lowerScoloredMass )
                      - highElectroweakinoMassRatio )
                    / ( 1.0 - highElectroweakinoMassRatio ) );
                heavierElectroweakinoMassFraction
                = ( ( ( secondElectroweakinoMass / lowerScoloredMass )
                      - highElectroweakinoMassRatio )
                    / ( 1.0 - highElectroweakinoMassRatio ) );
              }
              /* end of checking for successful interpolations. at this point,
               * shouldInterpolateOnElectroweakinos will be true unless 1 of the
               * interpolations went wrong, & all the values for the final
               * bilinear interpolation will also be appropriately set.
               */
            }
          }
          if( shouldInterpolateOnElectroweakinos )
            // if the interpolations so far were successful...
          {
            returnValue = lhcFaserGlobal::squareBilinearInterpolation(
                                             lighterElectroweakinoMassFraction,
                                             heavierElectroweakinoMassFraction,
                                        lighterLighterElectroweakinoPointValue,
                                        heavierLighterElectroweakinoPointValue,
                                        heavierHeavierElectroweakinoPointValue,
                                      lighterHeavierElectroweakinoPointValue );

            // debugging:
            /**std::cout << std::endl << "debugging:"
            << std::endl
            << "bilinear interpolation:" << std::endl
            << "lighterElectroweakinoMassFraction = "
            << lighterElectroweakinoMassFraction << std::endl
            << "heavierElectroweakinoMassFraction = "
            << heavierElectroweakinoMassFraction << std::endl
            << "lighterLighterElectroweakinoPointValue = "
            << lighterLighterElectroweakinoPointValue << std::endl
            << "heavierLighterElectroweakinoPointValue = "
            << heavierLighterElectroweakinoPointValue << std::endl
            << "heavierHeavierElectroweakinoPointValue = "
            << heavierHeavierElectroweakinoPointValue << std::endl
            << "lighterHeavierElectroweakinoPointValue = "
            << lighterHeavierElectroweakinoPointValue << std::endl
            << "returning " << returnValue;
            std::cout << std::endl;**/
          }
          else
          {
            returnValue = CppSLHA::CppSLHA_global::really_wrong_value;
          }
        }
      }  // end of interpolations.
    }  // end of if the point is on the grid.
    else
    {
      returnValue = CppSLHA::CppSLHA_global::really_wrong_value;
    }
    return returnValue;
  }



  jetAcceptanceTable::usedCascades const
  jetAcceptanceTable::lastEnumElement( jetAcceptanceTable::sizeOfEnum );
  int const
  jetAcceptanceTable::sizeOfUsedCascades(
                                    (int)jetAcceptanceTable::lastEnumElement );

  jetAcceptanceTable::jetAcceptanceTable(
                                    std::string const* const gridFilesLocation,
                                          std::string const* const jetCutName,
                                          int const acceptanceColumn,
                                    inputHandler const* const inputShortcut ) :
    jetCutName( *jetCutName ),
    acceptanceColumn( acceptanceColumn ),
    inputShortcut( inputShortcut ),
    heavierThanGluinoSquarkMass( inputShortcut ),
    useFirstMass( inputShortcut ),
    useSecondMass( inputShortcut ),
    useAverageMass( inputShortcut )
  {
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "pointers for" << std::endl
    << "heavierThanGluinoSquarkMass = " << &heavierThanGluinoSquarkMass
    << std::endl
    << "useFirstMass = " << &useFirstMass
    << std::endl
    << "useSecondMass = " << &useSecondMass
    << std::endl
    << "useAverageMass = " << &useAverageMass
    << std::endl;
    std::cout << std::endl;**/

    std::string gridFileBaseName( *gridFilesLocation );
    gridFileBaseName.append( "/" );
    gridFileBaseName.append( *jetCutName );
    std::string gridFileName( gridFileBaseName );
    gridFileName.append( "/gluino+gluino_acceptance.dat" );
    gluinoGluinoGrid = new jetAcceptanceGrid( &gridFileName,
                                              inputShortcut );
    gridFileName.assign( gridFileBaseName );
    gridFileName.append( "/squark+gluino_acceptance.dat" );
    squarkGluinoGrid = new jetAcceptanceGrid( &gridFileName,
                                              inputShortcut );
    gridFileName.assign( gridFileBaseName );
    gridFileName.append( "/squark+antisquark_acceptance.dat" );
    squarkAntisquarkGrid = new jetAcceptanceGrid( &gridFileName,
                                                  inputShortcut );
    gridFileName.assign( gridFileBaseName );
    gridFileName.append( "/squark+squark_acceptance.dat" );
    squarkSquarkGrid = new jetAcceptanceGrid( &gridFileName,
                                              inputShortcut );
    for( int vectorAdder( sizeOfUsedCascades );
         0 < vectorAdder;
         --vectorAdder )
    {
      gridMatrixRow = gridsMatrix.addNewAtEnd();
      for( int pairAdder( sizeOfUsedCascades );
           0 < pairAdder;
           --pairAdder )
      {
        gridMatrixRow->addNewAtEnd();
      }
    }

    // now we specify the various cases:
    gridMatrixRow = gridsMatrix.getPointer( (int)directGluino );
    gridMatrixElement = gridMatrixRow->getPointer( (int)directGluino );
    gridMatrixElement->first = gluinoGluinoGrid;
    gridMatrixElement->second = &heavierThanGluinoSquarkMass;
    gridMatrixElement = gridMatrixRow->getPointer( (int)directSquark );
    gridMatrixElement->first = squarkGluinoGrid;
    gridMatrixElement->second = &useSecondMass;
    gridMatrixElement = gridMatrixRow->getPointer( (int)gluinoToSquark );
    gridMatrixElement->first = gluinoGluinoGrid;
    gridMatrixElement->second = &useSecondMass;
    gridMatrixElement = gridMatrixRow->getPointer( (int)squarkToGluino );
    gridMatrixElement->first = squarkGluinoGrid;
    gridMatrixElement->second = &useSecondMass;

    gridMatrixRow = gridsMatrix.getPointer( (int)directSquark );
    gridMatrixElement = gridMatrixRow->getPointer( (int)directGluino );
    gridMatrixElement->first = squarkGluinoGrid;
    gridMatrixElement->second = &useFirstMass;
    gridMatrixElement = gridMatrixRow->getPointer( (int)directSquark );
    gridMatrixElement->first = squarkAntisquarkGrid;
    gridMatrixElement->second = &useAverageMass;
    gridMatrixElement = gridMatrixRow->getPointer( (int)gluinoToSquark );
    gridMatrixElement->first = squarkGluinoGrid;
    gridMatrixElement->second = &useAverageMass;
    gridMatrixElement = gridMatrixRow->getPointer( (int)squarkToGluino );
    gridMatrixElement->first = squarkAntisquarkGrid;
    gridMatrixElement->second = &useFirstMass;

    gridMatrixRow = gridsMatrix.getPointer( (int)gluinoToSquark );
    gridMatrixElement = gridMatrixRow->getPointer( (int)directGluino );
    gridMatrixElement->first = gluinoGluinoGrid;
    gridMatrixElement->second = &useFirstMass;
    gridMatrixElement = gridMatrixRow->getPointer( (int)directSquark );
    gridMatrixElement->first = squarkGluinoGrid;
    gridMatrixElement->second = &useAverageMass;
    gridMatrixElement = gridMatrixRow->getPointer( (int)gluinoToSquark );
    gridMatrixElement->first = gluinoGluinoGrid;
    gridMatrixElement->second = &useAverageMass;
    gridMatrixElement = gridMatrixRow->getPointer( (int)squarkToGluino );
    gridMatrixElement->first = gluinoGluinoGrid;
    gridMatrixElement->second = &useFirstMass;

    gridMatrixRow = gridsMatrix.getPointer( (int)squarkToGluino );
    gridMatrixElement = gridMatrixRow->getPointer( (int)directGluino );
    gridMatrixElement->first = squarkGluinoGrid;
    gridMatrixElement->second = &useFirstMass;
    gridMatrixElement = gridMatrixRow->getPointer( (int)directSquark );
    gridMatrixElement->first = squarkAntisquarkGrid;
    gridMatrixElement->second = &useSecondMass;
    gridMatrixElement = gridMatrixRow->getPointer( (int)gluinoToSquark );
    gridMatrixElement->first = gluinoGluinoGrid;
    gridMatrixElement->second = &useSecondMass;
    gridMatrixElement = gridMatrixRow->getPointer( (int)squarkToGluino );
    gridMatrixElement->first = squarkAntisquarkGrid;
    gridMatrixElement->second = &useAverageMass;

    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "gluinoGluinoGrid = " << gluinoGluinoGrid
    << ", squarkGluinoGrid = " << squarkGluinoGrid
    << ", squarkAntisquarkGrid = " << squarkAntisquarkGrid
    << ", squarkSquarkGrid = " << squarkSquarkGrid;
    std::cout << std::endl;**/
  }

  jetAcceptanceTable::~jetAcceptanceTable()
  {
    delete gluinoGluinoGrid;
    delete squarkGluinoGrid;
    delete squarkAntisquarkGrid;
    delete squarkSquarkGrid;
  }

  int
  jetAcceptanceTable::getIntForCascadeType(
                                         fullCascade const* const givenCascade,
                                     particlePointer* const squarkFromCascade )
  {
    usedCascades typeToCountAs;
    if( fullCascade::gluinoToEwino == givenCascade->getColorfulCascadeType() )
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "counts as gluinoToEwino";
      std::cout << std::endl;**/

      typeToCountAs = directGluino;
      *squarkFromCascade = NULL;
    }
    else if( fullCascade::squarkToEwino
             == givenCascade->getColorfulCascadeType() )
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "counts as squarkToEwino";
      std::cout << std::endl;**/

      typeToCountAs = directSquark;
      *squarkFromCascade = givenCascade->getInitialSparticle();
    }
    else if( fullCascade::squarkToSquarkThenMore
             == givenCascade->getColorfulCascadeType() )
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "counts as squarkToEwino";
      std::cout << std::endl;**/

      typeToCountAs = directSquark;
      *squarkFromCascade = givenCascade->getCascadeDefiner()->at( 1 )->first;
    }
    else if( fullCascade::squarkToGauginoThenMore
             == givenCascade->getColorfulCascadeType() )
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "counts as squarkToGluino";
      std::cout << std::endl;**/

      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "givenCascade->getCascadeDefiner()->size() = "
      << givenCascade->getCascadeDefiner()->size()
      << ", givenCascade = "
      << *(givenCascade->getInitialSparticle()->get_name());
      for( std::vector< fullCascade::particleWithInt*
                                                >::const_reverse_iterator
           cascadeParticleIterator(
                                 givenCascade->getCascadeDefiner()->rbegin() );
           givenCascade->getCascadeDefiner()->rend()
           > cascadeParticleIterator;
           ++cascadeParticleIterator )
      {
        std::cout
        << " => " << *((*cascadeParticleIterator)->first->get_name());
      }
      std::cout << std::endl;**/

      typeToCountAs = squarkToGluino;
      *squarkFromCascade = givenCascade->getInitialSparticle();
    }
    else
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "counts as gluinoToSquark";
      std::cout << std::endl;**/

      typeToCountAs = gluinoToSquark;
      *squarkFromCascade = givenCascade->getCascadeDefiner()->at( 1 )->first;
    }
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "returning " << (int)typeToCountAs;
    std::cout << std::endl;**/

    return (int)typeToCountAs;
  }

  double
  jetAcceptanceTable::getAcceptance(
                        signedParticleShortcutPair const* const initialPair,
                                     fullCascade const* const firstCascade,
                                     fullCascade const* const secondCascade )
  {
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "jetAcceptanceTable::getAcceptance( "
    << *(initialPair->getFirstParticle()->get_name_or_antiname(
                                       initialPair->firstIsNotAntiparticle() ))
    << "+"
    << *(initialPair->getSecondParticle()->get_name_or_antiname(
                                      initialPair->secondIsNotAntiparticle() ))
    << ", ";
    for( std::vector< fullCascade::particleWithInt* >::const_reverse_iterator
         cascadeParticleIterator(
                                 firstCascade->getCascadeDefiner()->rbegin() );
        firstCascade->getCascadeDefiner()->rend()
        > cascadeParticleIterator;
        ++cascadeParticleIterator )
    {
      std::cout << " => " << *((*cascadeParticleIterator)->first->get_name());
    }
    std::cout << ", ";
    for( std::vector< fullCascade::particleWithInt* >::const_reverse_iterator
         cascadeParticleIterator(
                                secondCascade->getCascadeDefiner()->rbegin() );
        secondCascade->getCascadeDefiner()->rend()
        > cascadeParticleIterator;
        ++cascadeParticleIterator )
    {
      std::cout << " => " << *((*cascadeParticleIterator)->first->get_name());
    }
    std::cout << " ) called.";
    std::cout << std::endl;**/

    gridMatrixRow = gridsMatrix.getPointer( getIntForCascadeType( firstCascade,
                                                       &firstCascadeSquark ) );
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "gridMatrixRow = " << gridMatrixRow;
    std::cout << std::endl;**/

    gridMatrixElement = gridMatrixRow->getPointer( getIntForCascadeType(
                                                                 secondCascade,
                                                      &secondCascadeSquark ) );
    /* gridMatrixElement is now a pair of pointers, 1 to the appropriate grid,
     * the other to the appropriate functor for obtaining the squark mass to
     * use. now we check to see if we need to redirect to squarkSquarkGrid:
     */
    gridToUse = gridMatrixElement->first;
    if( ( squarkAntisquarkGrid == gridToUse )
        &&
        ( initialPair->firstIsNotAntiparticle()
          == initialPair->secondIsNotAntiparticle() ) )
    {
      gridToUse = squarkSquarkGrid;
    }

    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "gridToUse = " << gridToUse
    << ", grid squark mass decider = " << gridMatrixElement->second
    << std::endl
    << "about to call"
    << " gridToUse->getAcceptance( (*(gridMatrixElement->second))( ";
    if( NULL == firstCascadeSquark )
    {
      std::cout << "NULL";
    }
    else
    {
      std::cout << *(firstCascadeSquark->get_name());
    }
    std::cout << ", ";
    if( NULL == secondCascadeSquark )
    {
      std::cout << "NULL";
    }
    else
    {
      std::cout << *(secondCascadeSquark->get_name());
    }
    std::cout << " ) = " << (*(gridMatrixElement->second))( firstCascadeSquark,
                                                          secondCascadeSquark )
    << ", " << inputShortcut->getGluinoMass() << ", "
    << firstCascade->getCascadeDefiner()->front()->first->get_absolute_mass()
    << ", "
    << secondCascade->getCascadeDefiner()->front()->first->get_absolute_mass()
    << ", " << acceptanceColumn << " ).";
    std::cout << std::endl;**/

    // now we ensure that we'll be asking for results from the grid:
    double squarkMass( (*(gridMatrixElement->second))( firstCascadeSquark,
                                                       secondCascadeSquark ) );
    if( gridToUse->getHighestSquarkMass() < squarkMass )
    {
      squarkMass = gridToUse->getHighestSquarkMass();
    }
    double gluinoMass( inputShortcut->getGluinoMass() );
    if( gridToUse->getHighestGluinoMass() < gluinoMass )
    {
      gluinoMass = gridToUse->getHighestGluinoMass();
    }
    double lighterScoloredMass( squarkMass );
    if( gluinoMass < lighterScoloredMass )
    {
      lighterScoloredMass = gluinoMass;
    }
    double firstEwinoMass( firstCascade->getCascadeDefiner()->front(
                                                )->first->get_absolute_mass());
    if( lighterScoloredMass < firstEwinoMass )
    {
      firstEwinoMass
      *= ( lighterScoloredMass / firstCascadeSquark->get_absolute_mass() );
    }
    double secondEwinoMass( secondCascade->getCascadeDefiner()->front(
                                               )->first->get_absolute_mass() );
    if( lighterScoloredMass < secondEwinoMass )
    {
      secondEwinoMass
      *= ( lighterScoloredMass / secondCascadeSquark->get_absolute_mass() );
    }

    // now we just use the general valueAt(...) function of acceptanceGrid:
    return
    gridToUse->getAcceptance( squarkMass,
                              gluinoMass,
                              firstEwinoMass,
                              secondEwinoMass,
                              acceptanceColumn );
  }



  jetAcceptanceTablesForOneBeamEnergy::jetAcceptanceTablesForOneBeamEnergy(
                                       inputHandler const* const inputShortcut,
                                                          int const beamEnergy,
                               std::string const* const gridFileSetLocation ) :
    inputShortcut( inputShortcut ),
    beamEnergy( beamEnergy ),
    gridFileSetLocation( *gridFileSetLocation )
  {
    std::stringstream energyStream( "" );
    energyStream << "/" << beamEnergy << "TeV/jets/";
    this->gridFileSetLocation.append( energyStream.str() );
  }

  jetAcceptanceTablesForOneBeamEnergy::~jetAcceptanceTablesForOneBeamEnergy()
  {
    for( std::vector< jetAcceptanceTable* >::iterator
         deletionIterator( jetTables.begin() );
         jetTables.end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
  }

  jetAcceptanceTable*
  jetAcceptanceTablesForOneBeamEnergy::getTable(
                                           std::string const* const jetCutName,
                                                 int const acceptanceColumn )
  // this returns the jetAcceptanceTable for the requested signal & column.
  {
    jetAcceptanceTable* returnPointer( NULL );
    // we look to see if we already have a jetAcceptanceTable for this jet+MET
    // cut set:
    for( std::vector< jetAcceptanceTable* >::iterator
         searchIterator( jetTables.begin() );
         jetTables.end() > searchIterator;
         ++searchIterator )
    {
      if((*searchIterator)->isRequested( jetCutName,
                                         acceptanceColumn ) )
      {
        returnPointer = *searchIterator;
        searchIterator = jetTables.end();
      }
    }
    if( NULL == returnPointer )
      // if we do not already have a jetAcceptanceTable for this jet+MET cut
      // set, we make a new instance:
    {
      returnPointer = new jetAcceptanceTable( &gridFileSetLocation,
                                              jetCutName,
                                              acceptanceColumn,
                                              inputShortcut );
      jetTables.push_back( returnPointer );
    }
    return returnPointer;
  }



  jetPlusMetAcceptanceHandler::jetPlusMetAcceptanceHandler(
                                    inputHandler const* const inputShortcut ) :
    inputShortcut( inputShortcut ),
    gridFileSetLocation( *(inputShortcut->inspectPathToKinematicsGrids()) )
  {
    // just an initialization list.
  }

  jetPlusMetAcceptanceHandler::~jetPlusMetAcceptanceHandler()
  {
    for( std::vector< jetAcceptanceTablesForOneBeamEnergy* >::iterator
         deletionIterator( acceptanceTables.begin() );
         acceptanceTables.end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
  }

  jetAcceptanceTablesForOneBeamEnergy*
  jetPlusMetAcceptanceHandler::getJetAcceptanceTablesForOneBeamEnergy(
                                                         int const beamEnergy )
  /* this looks to see if there is an existing
   * jetAcceptanceTablesForOneBeamEnergy with the requested beamEnergy, & if
   * not, makes 1, & returns the pointer.
   */
  {
    jetAcceptanceTablesForOneBeamEnergy* returnPointer( NULL );
    // we look to see if we already have a leptonAcceptancesForOneBeamEnergy
    // for these acceptanceValues:
    for( std::vector< jetAcceptanceTablesForOneBeamEnergy* >::iterator
         searchIterator( acceptanceTables.begin() );
         acceptanceTables.end() > searchIterator;
         ++searchIterator )
    {
      if( beamEnergy == (*searchIterator)->getBeamEnergy() )
      {
        returnPointer = *searchIterator;
        searchIterator = acceptanceTables.end();
      }
    }
    if( NULL == returnPointer )
      // if we do not already have a leptonAcceptanceParameterSet for this
      // electroweakino, we make a new instance:
    {
      returnPointer = new jetAcceptanceTablesForOneBeamEnergy( inputShortcut,
                                                               beamEnergy,
                                                        &gridFileSetLocation );
      acceptanceTables.push_back( returnPointer );
    }
    return returnPointer;
  }

}  // end of LHC_FASER namespace.

