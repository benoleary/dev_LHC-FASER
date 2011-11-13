/*
 * LHC-FASER_lepton_kinematics_stuff.cpp
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
 *      LHC-FASER also requires grids of lookup acceptanceValues. These should also be
 *      found in a subdirectory included with this package.
 */

#include "LHC-FASER_lepton_kinematics_stuff.hpp"

namespace LHC_FASER
{
  singleLeptonCut::singleLeptonCut() :
    acceptanceCutSet()
  {
    // just an initialization list.
  }

  singleLeptonCut::~singleLeptonCut()
  {
    // does nothing.
  }



  leptonAcceptanceGrid::leptonAcceptanceGrid(
                                     std::string const* const gridFileLocation,
                                    inputHandler const* const inputShortcut ) :
    acceptanceGrid( gridFileLocation,
                    inputShortcut )
  {
    // just an initialization list.
  }

  leptonAcceptanceGrid::~leptonAcceptanceGrid()
  {
    // does nothing.
  }



  leptonAcceptanceFromSquarkGrid::leptonAcceptanceFromSquarkGrid(
                                     std::string const* const gridFileLocation,
                                    inputHandler const* const inputShortcut ) :
      leptonAcceptanceGrid( gridFileLocation,
                            inputShortcut ),
      lowerLeftVector( NULL ),
      lowerRightVector( NULL ),
      upperRightVector( NULL ),
      upperLeftVector( NULL )
  {
    // just an initialization list.
  }

  leptonAcceptanceFromSquarkGrid::~leptonAcceptanceFromSquarkGrid()
  {
    // does nothing.
  }


  void
  leptonAcceptanceFromSquarkGrid::interpolateAcceptances( double squarkMass,
                                                          double gluinoMass,
                                               double const electroweakinoMass,
                                             double* const effectiveSquarkMass,
                                        double* const pseudorapidityAcceptance,
                               std::vector< double >* const energyAcceptances )
  // this interpolates the grid to obtain values for the acceptancesPerCutSet based on
  // the given masses.
  {
    if( gluinoMass < squarkMass)
    {
      gluinoMass = ( squarkMass + 1.0 );
    }
    energyAcceptances->clear();
    if( pointIsOnGrid( squarkMass,
                       gluinoMass ) )
    {
      if( ( indexForLightLightNeutralinoPair
            < lowerLeftVectorOfVectors->size() )
          &&
          ( indexForLightLightNeutralinoPair
            < lowerRightVectorOfVectors->size() )
          &&
          ( indexForLightLightNeutralinoPair
            < upperRightVectorOfVectors->size() )
          &&
          ( indexForLightLightNeutralinoPair
            < lowerRightVectorOfVectors->size() ) )
      {
        lowerLeftVector
        = lowerLeftVectorOfVectors->at( indexForLightLightNeutralinoPair );
        lowerRightVector
        = lowerRightVectorOfVectors->at( indexForLightLightNeutralinoPair );
        upperRightVector
        = upperRightVectorOfVectors->at( indexForLightLightNeutralinoPair );
        upperLeftVector
        = upperLeftVectorOfVectors->at( indexForLightLightNeutralinoPair );
        if( ( 3 <= lowerLeftVector->size() )
            &&
            ( lowerLeftVector->size() == lowerRightVector->size() )
            &&
            ( lowerLeftVector->size() == upperRightVector->size() )
            &&
            ( lowerLeftVector->size() == upperLeftVector->size() ) )
        {
          *effectiveSquarkMass = squarkMass;
          *pseudorapidityAcceptance
          = lhcFaserGlobal::squareBilinearInterpolation( squarkMassFraction,
                                                         gluinoMassFraction,
                                                      lowerLeftVector->at( 1 ),
                                                     lowerRightVector->at( 1 ),
                                                     upperRightVector->at( 1 ),
                                                    upperLeftVector->at( 1 ) );
          for( unsigned int acceptanceCounter( 2 );
               lowerLeftVector->size() > acceptanceCounter;
               ++acceptanceCounter )
          {
            energyAcceptances->push_back(
                                   lhcFaserGlobal::squareBilinearInterpolation(
                                                            squarkMassFraction,
                                                            gluinoMassFraction,
                                      lowerLeftVector->at( acceptanceCounter ),
                                     lowerRightVector->at( acceptanceCounter ),
                                     upperRightVector->at( acceptanceCounter ),
                                  upperLeftVector->at( acceptanceCounter ) ) );
          }
        }
      }
    }
    if( energyAcceptances->empty() )
      // if we didn't fill any entries in energyAcceptances, something was
      // wrong with trying to find grid points for interpolation.
    {
      *effectiveSquarkMass = CppSLHA::CppSLHA_global::really_wrong_value;
      *pseudorapidityAcceptance = CppSLHA::CppSLHA_global::really_wrong_value;
      energyAcceptances->push_back(
                                 CppSLHA::CppSLHA_global::really_wrong_value );
    }
  }



  leptonAcceptanceFromGluinoGrid::leptonAcceptanceFromGluinoGrid(
                                     std::string const* const gridFileLocation,
                                    inputHandler const* const inputShortcut ) :
      leptonAcceptanceGrid( gridFileLocation,
                            inputShortcut ),
      lowerElectroweakinoMassIndex( 0 ),
      upperElectroweakinoMassIndex( 0 ),
      electroweakinoMassFraction(
                                 CppSLHA::CppSLHA_global::really_wrong_value ),
      shouldInterpolateOnElectroweakino( false ),
      otherElectroweakinoValue( CppSLHA::CppSLHA_global::really_wrong_value ),
      foreLowerLeftVector( NULL ),
      foreLowerRightVector( NULL ),
      foreUpperRightVector( NULL ),
      foreUpperLeftVector( NULL ),
      rearLowerLeftVector( NULL ),
      rearLowerRightVector( NULL ),
      rearUpperRightVector( NULL ),
      rearUpperLeftVector( NULL )
  {
    // just an initialization list.
  }

  leptonAcceptanceFromGluinoGrid::~leptonAcceptanceFromGluinoGrid()
  {
    // does nothing.
  }


  void
  leptonAcceptanceFromGluinoGrid::interpolateAcceptances( double squarkMass,
                                                          double gluinoMass,
                                               double const electroweakinoMass,
                                             double* const effectiveSquarkMass,
                                        double* const pseudorapidityAcceptance,
                               std::vector< double >* const energyAcceptances )
  // this interpolates the grid to obtain values for the acceptancesPerCutSet based on
  // the given masses.
  {
    // 1st we fudge for gluinos not decaying to squarks even if they can at
    // this point:
    if( squarkMass < gluinoMass )
    {
      squarkMass = ( gluinoMass + 1.0 );
    }
    // then we ensure that we set up a fresh set of acceptancesPerCutSet:
    energyAcceptances->clear();
    if( ( electroweakinoMass < gluinoMass )
        &&
        pointIsOnGrid( squarkMass,
                       gluinoMass ) )
      // if the grid is valid for this point...
    {
      // either we approximate very light electroweakinos by the lightest grid
      // point...
      if( electroweakinoMass <= ( gluinoMass * lowElectroweakinoMassRatio ) )
      {
        lowerElectroweakinoMassIndex = indexForLightLightNeutralinoPair;
        shouldInterpolateOnElectroweakino = false;
      }
      // or we interpolate between light & medium electroweakino masses...
      else if( electroweakinoMass
               <= ( gluinoMass * mediumElectroweakinoMassRatio ) )
      {
        lowerElectroweakinoMassIndex = indexForLightLightNeutralinoPair;
        upperElectroweakinoMassIndex = indexForMediumMediumNeutralinoPair;
        electroweakinoMassFraction
        = ( ( ( electroweakinoMass / gluinoMass ) - lowElectroweakinoMassRatio )
            / ( mediumElectroweakinoMassRatio - lowElectroweakinoMassRatio ) );
        shouldInterpolateOnElectroweakino = true;
      }
      // or we interpolate between medium & heavy electroweakino masses...
      else if( electroweakinoMass
               <= ( gluinoMass * highElectroweakinoMassRatio ) )
      {
        lowerElectroweakinoMassIndex = indexForMediumMediumNeutralinoPair;
        upperElectroweakinoMassIndex = indexForHeavyHeavyNeutralinoPair;
        electroweakinoMassFraction
        = ( ( ( electroweakinoMass / gluinoMass ) - mediumElectroweakinoMassRatio )
            / ( highElectroweakinoMassRatio - mediumElectroweakinoMassRatio ) );
        shouldInterpolateOnElectroweakino = true;
      }
      else
        /* or we approximate very heavy electroweakinos by the heaviest grid
         * point (with effective squark mass interpolating to 0.0 at the
         * electroweakino mass = gluino mass edge of the grid)...
         */
      {
        lowerElectroweakinoMassIndex = indexForHeavyHeavyNeutralinoPair;
        electroweakinoMassFraction
        = ( ( ( electroweakinoMass / gluinoMass ) - highElectroweakinoMassRatio )
            / ( 1.0 - highElectroweakinoMassRatio ) );
        shouldInterpolateOnElectroweakino = false;
      }
      // now we check that the required grid point exists:
      if( ( lowerElectroweakinoMassIndex
            < lowerLeftVectorOfVectors->size() )
          &&
          ( lowerElectroweakinoMassIndex
            < lowerRightVectorOfVectors->size() )
          &&
          ( lowerElectroweakinoMassIndex
            < upperRightVectorOfVectors->size() )
          &&
          ( lowerElectroweakinoMassIndex
            < lowerRightVectorOfVectors->size() ) )
      {
        // if the grid point exists, we set up the vectors of acceptancesPerCutSet:
        foreLowerLeftVector
        = lowerLeftVectorOfVectors->at( lowerElectroweakinoMassIndex );
        foreLowerRightVector
        = lowerRightVectorOfVectors->at( lowerElectroweakinoMassIndex );
        foreUpperRightVector
        = upperRightVectorOfVectors->at( lowerElectroweakinoMassIndex );
        foreUpperLeftVector
        = upperLeftVectorOfVectors->at( lowerElectroweakinoMassIndex );
        if( ( 2 < foreLowerLeftVector->size() )
            &&
            ( foreLowerLeftVector->size() == foreLowerRightVector->size() )
            &&
            ( foreLowerLeftVector->size() == foreUpperRightVector->size() )
            &&
            ( foreLowerLeftVector->size() == foreUpperLeftVector->size() ) )
          // if the acceptance vectors have enough entries...
        {
          // we set up all the acceptancesPerCutSet assuming that they don't have to be
          // interpolated, & later we interpolate if necessary:
          *effectiveSquarkMass
          = lhcFaserGlobal::squareBilinearInterpolation( squarkMassFraction,
                                                         gluinoMassFraction,
                                                  foreLowerLeftVector->front(),
                                                 foreLowerRightVector->front(),
                                                 foreUpperRightVector->front(),
                                                foreUpperLeftVector->front() );
          if( electroweakinoMass
              > ( gluinoMass * highElectroweakinoMassRatio ) )
            // if we need to interpolate the effective squark mass because the
            // electroweakino mass is so high...
          {
            *effectiveSquarkMass = lhcFaserGlobal::unitLinearInterpolation(
                                                    electroweakinoMassFraction,
                                                          *effectiveSquarkMass,
                                                                         0.0 );
          }
          *pseudorapidityAcceptance
          = lhcFaserGlobal::squareBilinearInterpolation( squarkMassFraction,
                                                         gluinoMassFraction,
                                                  foreLowerLeftVector->at( 1 ),
                                                 foreLowerRightVector->at( 1 ),
                                                 foreUpperRightVector->at( 1 ),
                                                foreUpperLeftVector->at( 1 ) );
          if( shouldInterpolateOnElectroweakino
              &&
              ( upperElectroweakinoMassIndex
                < lowerLeftVectorOfVectors->size() )
              &&
              ( upperElectroweakinoMassIndex
                < lowerRightVectorOfVectors->size() )
              &&
              ( upperElectroweakinoMassIndex
                < upperRightVectorOfVectors->size() )
              &&
              ( upperElectroweakinoMassIndex
                < lowerRightVectorOfVectors->size() ) )
            /* if we need to interpolate the values,
             * shouldInterpolateOnElectroweakinos will have been set to true &
             * upperElectroweakinoMassIndex will have been set.
             */
          {
            rearLowerLeftVector
            = lowerLeftVectorOfVectors->at( upperElectroweakinoMassIndex );
            rearLowerRightVector
            = lowerRightVectorOfVectors->at( upperElectroweakinoMassIndex );
            rearUpperRightVector
            = upperRightVectorOfVectors->at( upperElectroweakinoMassIndex );
            rearUpperLeftVector
            = upperLeftVectorOfVectors->at( upperElectroweakinoMassIndex );
            if( ( foreLowerLeftVector->size() == rearLowerLeftVector->size() )
                &&
                ( foreLowerLeftVector->size() == rearLowerRightVector->size() )
                &&
                ( foreLowerLeftVector->size() == rearUpperRightVector->size() )
                &&
                ( foreLowerLeftVector->size()
                  == rearUpperLeftVector->size() ) )
            {
              otherElectroweakinoValue
              = lhcFaserGlobal::squareBilinearInterpolation(
                                                            squarkMassFraction,
                                                            gluinoMassFraction,
                                                  rearLowerLeftVector->front(),
                                                 rearLowerRightVector->front(),
                                                 rearUpperRightVector->front(),
                                                rearUpperLeftVector->front() );
              *effectiveSquarkMass = lhcFaserGlobal::unitLinearInterpolation(
                                                    electroweakinoMassFraction,
                                                          *effectiveSquarkMass,
                                                    otherElectroweakinoValue );
              otherElectroweakinoValue
              = lhcFaserGlobal::squareBilinearInterpolation(
                                                            squarkMassFraction,
                                                            gluinoMassFraction,
                                                  rearLowerLeftVector->at( 1 ),
                                                 rearLowerRightVector->at( 1 ),
                                                 rearUpperRightVector->at( 1 ),
                                                rearUpperLeftVector->at( 1 ) );
              *pseudorapidityAcceptance
              = lhcFaserGlobal::unitLinearInterpolation(
                                                    electroweakinoMassFraction,
                                                     *pseudorapidityAcceptance,
                                                    otherElectroweakinoValue );
              for( unsigned int acceptanceCounter( 2 );
                   foreLowerLeftVector->size() > acceptanceCounter;
                   ++acceptanceCounter )
              {
                energyAcceptances->push_back(
                                       lhcFaserGlobal::unitLinearInterpolation(
                                                    electroweakinoMassFraction,
                                   lhcFaserGlobal::squareBilinearInterpolation(
                                                            squarkMassFraction,
                                                            gluinoMassFraction,
                                  foreLowerLeftVector->at( acceptanceCounter ),
                                 foreLowerRightVector->at( acceptanceCounter ),
                                 foreUpperRightVector->at( acceptanceCounter ),
                                foreUpperLeftVector->at( acceptanceCounter ) ),
                                   lhcFaserGlobal::squareBilinearInterpolation(
                                                            squarkMassFraction,
                                                            gluinoMassFraction,
                                  rearLowerLeftVector->at( acceptanceCounter ),
                                 rearLowerRightVector->at( acceptanceCounter ),
                                 rearUpperRightVector->at( acceptanceCounter ),
                            rearUpperLeftVector->at( acceptanceCounter ) ) ) );
              }  // end of for loop filling energy acceptancesPerCutSet.
            }
            // end of if all the rear vectors have the same size as the fore
            // vectors.
          }
          else
            // otherwise, we don't interpolate on the electroweakino mass, but
            // just use the values from interpolation on the scolored masses:
          {
            for( unsigned int acceptanceCounter( 2 );
                 foreLowerLeftVector->size() > acceptanceCounter;
                ++acceptanceCounter )
            {
              energyAcceptances->push_back(
                                   lhcFaserGlobal::squareBilinearInterpolation(
                                                            squarkMassFraction,
                                                            gluinoMassFraction,
                                  foreLowerLeftVector->at( acceptanceCounter ),
                                 foreLowerRightVector->at( acceptanceCounter ),
                                 foreUpperRightVector->at( acceptanceCounter ),
                              foreUpperLeftVector->at( acceptanceCounter ) ) );
            }  // end of for loop filling energy acceptancesPerCutSet.
          }  // end of whether we should interpolate on the electroweakino.
        }  // end of if the acceptance vectors were the correct size.
      }
      // end of if vectors were large enough for the required electroweakino
      // indices.
    }  // end of if electroweakinoMass < gluinoMass && pointIsOnGrid(...).
    if( energyAcceptances->empty() )
      // if we didn't fill any entries in energyAcceptances, something was
      // wrong with trying to find grid points for interpolation.
    {
      *effectiveSquarkMass = CppSLHA::CppSLHA_global::really_wrong_value;
      *pseudorapidityAcceptance = CppSLHA::CppSLHA_global::really_wrong_value;
      energyAcceptances->push_back(
                                 CppSLHA::CppSLHA_global::really_wrong_value );
    }
  }



  double const leptonAcceptanceParameterSet::defaultBinSize( 2.0 );
  double const
  leptonAcceptanceParameterSet::defaultTransverseMomentumCut( 10.0 );

  leptonAcceptanceParameterSet::leptonAcceptanceParameterSet(
                                       inputHandler const* const inputShortcut,
                                    leptonAcceptanceGrid* const acceptanceGrid,
                                                particlePointer const scolored,
                                                   particlePointer const ewino,
                                                          double const binSize,
                                         double const transverseMomentumCut ) :
    getsReadiedForNewPoint( inputShortcut->getReadier() ),
    effectiveSquarkMassHolder(),
    inputShortcut( inputShortcut ),
    acceptanceGrid( acceptanceGrid ),
    scolored( scolored ),
    ewino( ewino ),
    binSize( binSize ),
    transverseMomentumCut( transverseMomentumCut ),
    effectiveSquarkMass( CppSLHA::CppSLHA_global::really_wrong_value ),
    acceptanceBins(),
    pseudorapidityAcceptance( CppSLHA::CppSLHA_global::really_wrong_value ),
    returnValue( CppSLHA::CppSLHA_global::really_wrong_value ),
    binFraction( CppSLHA::CppSLHA_global::really_wrong_value ),
    lowerBin( (int)CppSLHA::CppSLHA_global::really_wrong_value ),
    acceptanceCounter( (int)CppSLHA::CppSLHA_global::really_wrong_value ),
    currentAcceptance( CppSLHA::CppSLHA_global::really_wrong_value )
  {
    acceptanceBins.push_back( CppSLHA::CppSLHA_global::really_wrong_value );
    // this is to ensure that there is at least 1 entry for the transverse
    // momentum cut.
  }

  leptonAcceptanceParameterSet::~leptonAcceptanceParameterSet()
  {
    // does nothing.
  }


  double
  leptonAcceptanceParameterSet::calculateAcceptanceAt(
                                                      double const givenEnergy,
                                                       double const givenCut )
  //const
  /* this interpolates the acceptanceValues in acceptanceBins to the requested
   * value, or returns pseudorapidityAcceptance if it's lower, scaled to the
   * given value for the transverse momentum cut.
   */
  {
    // debugging:
    /**std::cout
    << std::endl
    << "debugging: leptonAcceptanceParameterSet::acceptanceAt( "
    << energy << ", " << cut << " ) called";
    std::cout << std::endl;**/

    returnValue = 0.0;
    if( ( 0.0 < givenEnergy )
        &&
        ( 0.0 < givenCut ) )
    {
      binFraction = ( givenEnergy * ( transverseMomentumCut
                                      / ( givenCut * binSize ) ) );
      // now binFraction is givenEnergy scaled to the given momentum cut,
      // in units of binSize. importantly, binFraction > 0.0 still.
      lowerBin = (int)binFraction;
      /* now lowerBin is 1 more than the lower bin for givenEnergy because
       * "bin -1" is the acceptance for 0.0 GeV, which is 0.0 by assumption.
       * "bin 0" is for 2.0 GeV by default, "bin 1" for 4.0 GeV, & so on.
       */
      binFraction -= (double)lowerBin;
      // now binFraction is the fraction of the bin along from lowerBin
      // that givenEnergy corresponds to.
      --lowerBin;
      // now lowerBin corresponds to the correct bin for givenEnergy.

      // debugging:
      /**std::cout
      << std::endl
      << "binFraction = " << binFraction
      << ", lowerBin = " << lowerBin
      << ", acceptanceBins.size() = " << acceptanceBins.size();
      std::cout << std::endl;**/

      if( -1 == lowerBin )
        // if we're in the region 0.0 GeV to binSize GeV...
      {
        // debugging:
        /**std::cout
        << std::endl
        << "( -1 == lowerBin ) is true, binFraction = "
        << binFraction << ", acceptanceBins.front() = "
        << acceptanceBins.front();
        std::cout << std::endl;**/

        returnValue = ( givenEnergy * acceptanceBins.front() );
        // we assume that the acceptance at 0.0 GeV is 0.0, so the linear
        // interpolation is simple.
      }
      else if( ( 0 <= lowerBin )
               &&
               ( acceptanceBins.size() > (unsigned int)( lowerBin + 1 ) ) )
        // if we're in a region between 2 bin entries...
      {
        // debugging:
        /**std::cout
        << std::endl
        << "( ( ( 0 <= lowerBin ) &&"
        << " ( acceptanceBins.size() > ( lowerBin + 1 ) ) ) is true,"
        << " acceptanceBins.at( lowerBin ) = "
        << acceptanceBins.at( lowerBin )
        << ", binFraction = "
        << binFraction
        << ", acceptanceBins.at( ( lowerBin + 1 ) ) = "
        << acceptanceBins.at( ( lowerBin + 1 ) );
        std::cout << std::endl;**/

        returnValue = ( acceptanceBins.at( lowerBin )
                        + binFraction * ( acceptanceBins.at( ( lowerBin + 1 ) )
                                          - acceptanceBins.at( lowerBin ) ) );
        // we return a linear interpolation.
      }
      else if( acceptanceBins.size() <= (unsigned int)( lowerBin + 1 ) )
        // if we're in a region beyond the bins...
      {
        // debugging:
        /**std::cout
        << std::endl
        << "( acceptanceBins.size() <= ( lowerBin + 1 ) ) is true,"
        << " acceptanceBins.back() = " << acceptanceBins.back();
        std::cout << std::endl;**/

        returnValue = acceptanceBins.back();
        // we assume the acceptance is dominated by the pseudorapidity cut
        // by this point, so is constant out to any higher energy.
      }
      // this is a sanity check that I have decided is unnecessary:
      /*if( pseudorapidityAcceptance < returnValue )
      {
        returnValue = pseudorapidityAcceptance;
      }*/
    }
    return returnValue;
  }



  leptonAcceptancesForOneScolored::leptonAcceptancesForOneScolored(
                                       inputHandler const* const inputShortcut,
                                                particlePointer const scolored,
                                    leptonAcceptanceGrid* const acceptanceGrid,
                                                          double const binSize,
                                         double const transverseMomentumCut ) :
    inputShortcut( inputShortcut ),
    scolored( scolored ),
    acceptanceGrid( acceptanceGrid ),
    parameterSets(),
    binSize( binSize ),
    transverseMomentumCut( transverseMomentumCut )
  {
    // just an initialization list.
  }

  leptonAcceptancesForOneScolored::~leptonAcceptancesForOneScolored()
  {
    for( std::vector< leptonAcceptanceParameterSet* >::iterator
         deletionIterator( parameterSets.begin() );
         parameterSets.end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
  }


  leptonAcceptanceParameterSet*
  leptonAcceptancesForOneScolored::getParameterSet(
                                                  particlePointer const ewino )
  // this returns the leptonAcceptanceParameterSet for the requested
  // electroweakino.
  {
    leptonAcceptanceParameterSet* returnPointer( NULL );
    // we look to see if we already have a leptonAcceptanceParameterSet for
    // this electroweakino:
    for( std::vector< leptonAcceptanceParameterSet* >::iterator
         searchIterator( parameterSets.begin() );
         parameterSets.end() > searchIterator;
         ++searchIterator )
    {
      if( ewino == (*searchIterator)->getEwino() )
      {
        returnPointer = *searchIterator;
        searchIterator = parameterSets.end();
      }
    }
    if( NULL == returnPointer )
      // if we do not already have a leptonAcceptanceParameterSet for this
      // electroweakino, we make a new instance:
    {
      returnPointer = new leptonAcceptanceParameterSet( inputShortcut,
                                                        acceptanceGrid,
                                                        scolored,
                                                        ewino,
                                                        binSize,
                                                       transverseMomentumCut );
      parameterSets.push_back( returnPointer );
    }
    return returnPointer;
  }



  leptonAcceptancesForOneBeamEnergy::leptonAcceptancesForOneBeamEnergy(
                                       inputHandler const* const inputShortcut,
                                                          int const beamEnergy,
                                  std::string const* const gridFileSetLocation,
                                                          double const binSize,
                                         double const transverseMomentumCut ) :
    inputShortcut( inputShortcut ),
    beamEnergy( beamEnergy ),
    gridFileLocation( *gridFileSetLocation ),
    acceptanceFromSquarkGrid( NULL ),
    acceptanceFromGluinoGrid( NULL ),
    squarkAcceptanceSets(),
    gluinoAcceptanceSet( NULL ),
    binSize( binSize ),
    transverseMomentumCut( transverseMomentumCut )
  {
    std::stringstream energyStream( "" );
    energyStream << "/" << beamEnergy << "TeV/leptons/";
    gridFileLocation.append( energyStream.str() );
  }

  leptonAcceptancesForOneBeamEnergy::~leptonAcceptancesForOneBeamEnergy()
  {
    for( std::vector< leptonAcceptancesForOneScolored* >::iterator
         deletionIterator( squarkAcceptanceSets.begin() );
         squarkAcceptanceSets.end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
    delete gluinoAcceptanceSet;
    delete acceptanceFromSquarkGrid;
    delete acceptanceFromGluinoGrid;
  }


  leptonAcceptancesForOneScolored*
  leptonAcceptancesForOneBeamEnergy::getParameterSets(
                                               particlePointer const scolored )
  // this returns the leptonAcceptancesForOneScolored for the requested
  // colored sparticle.
  {
    leptonAcceptancesForOneScolored* returnPointer( NULL );
    // we look to see if we already have a leptonAcceptancesForOneScolored for
    // this colored sparticle:
    if( inputShortcut->getGluino() == scolored )
    {
      if( NULL == gluinoAcceptanceSet )
      {
        std::string gridFileName( gridFileLocation );
        gridFileName.append( "/gluino+gluino_acceptance.dat" );
        acceptanceFromGluinoGrid
        = new leptonAcceptanceFromGluinoGrid( &gridFileName,
                                              inputShortcut );
        gluinoAcceptanceSet
        = new leptonAcceptancesForOneScolored( inputShortcut,
                                               scolored,
                                               acceptanceFromGluinoGrid,
                                               binSize,
                                               transverseMomentumCut );
      }
      returnPointer = gluinoAcceptanceSet;
    }
    else
    {
      for( std::vector< leptonAcceptancesForOneScolored* >::iterator
           searchIterator( squarkAcceptanceSets.begin() );
           squarkAcceptanceSets.end() > searchIterator;
           ++searchIterator )
      {
        if( scolored == (*searchIterator)->getScolored() )
        {
          returnPointer = *searchIterator;
          searchIterator = squarkAcceptanceSets.end();
        }
      }
      if( NULL == returnPointer )
        // if we do not already have a leptonAcceptanceParameterSet for this
        // squark, we make a new instance:
      {
        if( NULL == acceptanceFromSquarkGrid )
        {
          std::string gridFileName( gridFileLocation );
          gridFileName.append( "/squark+antisquark_acceptance.dat" );
          acceptanceFromSquarkGrid
          = new leptonAcceptanceFromSquarkGrid( &gridFileName,
                                                inputShortcut );
        }
        returnPointer = new leptonAcceptancesForOneScolored( inputShortcut,
                                                             scolored,
                                                      acceptanceFromSquarkGrid,
                                                             binSize,
                                                       transverseMomentumCut );
        squarkAcceptanceSets.push_back( returnPointer );
      }
    }
    return returnPointer;
  }



  leptonAcceptanceHandler::leptonAcceptanceHandler(
                                       inputHandler const* const inputShortcut,
                               std::string const* const gridFileSetLocation ) :
    inputShortcut( inputShortcut ),
    gridFileSetLocation( *gridFileSetLocation ),
    acceptanceTables()
  {
    // just an initialization list.
  }

  leptonAcceptanceHandler::~leptonAcceptanceHandler()
  {
    for( std::vector< leptonAcceptancesForOneBeamEnergy* >::iterator
         deletionIterator( acceptanceTables.begin() );
         acceptanceTables.end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
  }


  leptonAcceptancesForOneBeamEnergy*
  leptonAcceptanceHandler::getLeptonAcceptancesForOneBeamEnergy(
                                                          int const beamEnergy,
                                                          double const binSize,
                                           double const transverseMomentumCut )
  /* this looks to see if there is an existing
   * leptonAcceptancesForOneBeamEnergy with the requested acceptances, & if
   * not, makes 1, & returns the pointer.
   */
  {
    leptonAcceptancesForOneBeamEnergy* returnPointer( NULL );
    // we look to see if we already have a leptonAcceptancesForOneBeamEnergy
    // for these acceptances:
    for( std::vector< leptonAcceptancesForOneBeamEnergy* >::iterator
         searchIterator( acceptanceTables.begin() );
         acceptanceTables.end() > searchIterator;
         ++searchIterator )
    {
      if( (*searchIterator)->isRequested( beamEnergy,
                                          binSize,
                                          transverseMomentumCut ) )
      {
        returnPointer = *searchIterator;
        searchIterator = acceptanceTables.end();
      }
    }
    if( NULL == returnPointer )
      // if we do not already have a leptonAcceptancesForOneBeamEnergy for the
      // requested acceptancesPerCutSet, we make a new instance:
    {
      returnPointer = new leptonAcceptancesForOneBeamEnergy( inputShortcut,
                                                             beamEnergy,
                                                          &gridFileSetLocation,
                                                             binSize,
                                                       transverseMomentumCut );
      acceptanceTables.push_back( returnPointer );
    }
    return returnPointer;
  }

}  // end of LHC_FASER namespace.

