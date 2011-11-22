/*
 * LHC-FASER_full_cascade_stuff.cpp
 *
 *  Created on: 30 Mar 2011
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

#include "LHC-FASER_full_cascade_stuff.hpp"

namespace LHC_FASER
{
  int const
  fullCascade::numberOfSmFermionsFromElectroweakDecaysPerPartialCascade( 2 );
  int const
  fullCascade::numberOfSmFermionsFromElectroweakDecaysPerFullCascade( 4 );

  fullCascade::fullCascade( colorfulCascadeType const typeOfColorfulCascade,
                            int const firstDecayBodyNumber,
                            inputHandler const* const inputShortcut,
                            particlePointer const initialSparticle,
                            double const beamEnergy ) :
    inputShortcut( inputShortcut ),
    initialSparticle( initialSparticle ),
    beamEnergy( beamEnergy ),
    typeOfColorfulCascade( typeOfColorfulCascade ),
    firstDecayBodyNumber( firstDecayBodyNumber ),
    subcascadePointer( NULL ),
    cascadeDefiner(),
    cascadeSegment( NULL ),
    soughtDecayProductList( 1,
                            CppSLHA::CppSLHA_global::really_wrong_value ),
    cachedBranchingRatio( CppSLHA::CppSLHA_global::really_wrong_value ),
    branchingRatioNeedsToBeRecalculated( true )
  {
    // just an initialization list.
  }

  fullCascade::~fullCascade()
  {
    // does nothing.
  }


  double
  fullCascade::specifiedJetsOneOssfMinusOsdfPair(
                                        acceptanceCutSet* const acceptanceCuts,
                                             int const numberOfAdditionalJets )
  /* this calls
   * ewinoCascades->getOssfMinusOsdf( scoloredIsNotAntiparticle,
   *                                 acceptanceCuts ),
   * & if bosonCascades is not NULL, it calls
   * ewinoCascades->getAcceptance( [ { 0, 1, 2 } jets
   *                                 + 0 leptons for acceptanceCuts ] )
   * & the same for bosonCascades, as well as
   * ->getOssfMinusOsdf( acceptanceCuts ).
   */
  {
    // we should really check that ewinoCascade is not NULL, but this function
    // should never be called on a fullCascade which has a NULL ewinoCascade.
    if( NULL == bosonCascades )
    {
      /* if there is no possibility of extra jets from a vector decay, we can
       * only get a single OSSF-OSDF pair from an electroweakino decay, which
       * cannot produce additional jets.
       */
      if( 0 == numberOfAdditionalJets )
      {
        return ewinoCascades->getOssfMinusOsdf( acceptanceCuts );
      }
      else
      {
        return 0.0;
      }
    }
    else
    {
      return ( ewinoCascades->getOssfMinusOsdf( acceptanceCuts )
               * bosonCascades->getAcceptance( acceptanceCuts,
                                               numberOfAdditionalJets,
                                               0,
                                               0,
                                               0,
                                               0 )
               + bosonCascades->getOssfMinusOsdf( acceptanceCuts )
                 * ewinoCascades->getAcceptance( acceptanceCuts,
                                                 numberOfAdditionalJets,
                                                 0,
                                                 0,
                                                 0,
                                                 0 ) );
    }
  }

  double
  fullCascade::unspecifiedJetsSpecifiedOssfMinusOsdfPairs(
                                        acceptanceCutSet* const acceptanceCuts,
                                                int const numberOfLeptonPairs )
  /* this only calls ewinoCascades->getOssfMinusOsdf( acceptanceCuts ),
   * & if bosonCascades is not NULL,
   * bosonCascades->getOssfMinusOsdf( acceptanceCuts ).
   */
  {
    if( 0 == numberOfLeptonPairs )
    {
      return unspecifiedJetsSpecifiedChargeSummedLeptons( acceptanceCuts,
                                                          0 );
    }
    else
    {
      // we should check that ewinoCascade is not NULL, but this function
      // should never be called on a fullCascade which has a NULL ewinoCascade.
      if( NULL == bosonCascades )
      {
        /* if there is no possibility of extra jets from a vector decay, we can
         * only get a single OSSF-OSDF pair from an electroweakino decay, which
         * cannot produce additional jets.
         */
        if( 1 == numberOfLeptonPairs )
        {
          return ewinoCascades->getOssfMinusOsdf( acceptanceCuts );
        }
        else
          // if there's no vector cascade, the only possibilities for non-zero
          // results are 0 pairs or 1 pair.
        {
          return 0.0;
        }
      }
      else
      {
        if( 2 == numberOfLeptonPairs )
        {
          return ( ewinoCascades->getOssfMinusOsdf( acceptanceCuts )
                   * bosonCascades->getOssfMinusOsdf( acceptanceCuts ) );
        }
        else if( 1 == numberOfLeptonPairs )
        {
        return ( ewinoCascades->getOssfMinusOsdf( acceptanceCuts )
                 * ( bosonCascades->getAcceptance( acceptanceCuts,
                                                   0,
                                                   0,
                                                   0,
                                                   0,
                                                   0 )
                     + bosonCascades->getAcceptance( acceptanceCuts,
                                                     1,
                                                     0,
                                                     0,
                                                     0,
                                                     0 )
                     + bosonCascades->getAcceptance( acceptanceCuts,
                                                     2,
                                                     0,
                                                     0,
                                                     0,
                                                     0 ) )
                 + bosonCascades->getOssfMinusOsdf( acceptanceCuts )
                   * ( ewinoCascades->getAcceptance( acceptanceCuts,
                                                     0,
                                                     0,
                                                     0,
                                                     0,
                                                     0 )
                       + ewinoCascades->getAcceptance( acceptanceCuts,
                                                       1,
                                                       0,
                                                       0,
                                                       0,
                                                       0 )
                       + ewinoCascades->getAcceptance( acceptanceCuts,
                                                       2,
                                                       0,
                                                       0,
                                                       0,
                                                       0 ) ) );
        }
        else
          /* the case of 0 pairs has already been checked for, so by this
           * point, we've discovered that we were asked for a number of pairs
           * that is not 0, 1, or 2. thus we return 0.0 because we assume that
           * any cascades which could result in so many leptons are negligible.
           */
        {
          return 0.0;
        }
      }  // end of if vectorCascade is NULL or not.
    }  // end of if numberOfLeptonPairs == 0 or not.
  }

  double
  fullCascade::specifiedJetsSpecifiedChargeSummedLeptons(
                                        acceptanceCutSet* const acceptanceCuts,
                                              int const numberOfAdditionalJets,
                                                    int const numberOfLeptons )
  {
    if( ( 0 > numberOfAdditionalJets )
        ||
        ( 0 > numberOfLeptons )
        ||
        ( numberOfSmFermionsFromElectroweakDecaysPerFullCascade
          < ( numberOfLeptons + numberOfAdditionalJets ) ) )
    {
      return 0.0;
    }
    else
    {
      double returnDouble( 0.0 );
      for( int numberOfNegativeElectrons( numberOfLeptons );
           0 <= numberOfNegativeElectrons;
           --numberOfNegativeElectrons )
      {
        for( int numberOfPositiveElectrons( ( numberOfLeptons
                                              - numberOfNegativeElectrons ) );
             0 <= numberOfPositiveElectrons;
             --numberOfPositiveElectrons )
        {
          for( int numberOfNegativeMuons( ( numberOfLeptons
                                           - numberOfNegativeElectrons
                                           - numberOfPositiveElectrons ) );
               0 <= numberOfNegativeMuons;
               --numberOfNegativeMuons )
          {
            returnDouble += getAcceptance( true
                                /* doesn't matter for charge-summed leptons */,
                                           acceptanceCuts,
                                           numberOfAdditionalJets,
                                           numberOfNegativeElectrons,
                                           numberOfPositiveElectrons,
                                           numberOfNegativeMuons,
                                           ( numberOfLeptons
                                             - numberOfNegativeElectrons
                                             - numberOfPositiveElectrons
                                             - numberOfNegativeMuons ) );
          }  // end of loop over negative muons.
        }  // end of loop over positive electrons.
      }  // end of loop over negative electrons.
      return returnDouble;
    }  // end of if numberOfLeptons was in the allowed range.
  }

  double
  fullCascade::unspecifiedJetsSpecifiedChargeSummedLeptons(
                                        acceptanceCutSet* const acceptanceCuts,
                                                    int const numberOfLeptons )
  {
    if( ( 0 > numberOfLeptons )
        ||
        ( numberOfSmFermionsFromElectroweakDecaysPerFullCascade
          < numberOfLeptons ) )
    {
      return 0.0;
    }
    else
    {
      double returnDouble( 0.0 );
      for( int numberOfJets(
                        ( numberOfSmFermionsFromElectroweakDecaysPerFullCascade
                          - numberOfLeptons ) );
           0 <= numberOfJets;
           --numberOfJets )
      {
        returnDouble
        += specifiedJetsSpecifiedChargeSummedLeptons( acceptanceCuts,
                                                      numberOfJets,
                                                      numberOfLeptons );
      }  // end of loop over jets.
      return returnDouble;
    }  // end of if numberOfLeptons was in the allowed range.
  }



  namespace fullCascadeType
  {
    sdownDirectlyToElectroweak::sdownDirectlyToElectroweak(
                                       inputHandler const* const inputShortcut,
                                           particlePointer const initialSquark,
                                                       double const beamEnergy,
                                 electroweakCascadeSet* const ewinoCascades ) :
        fullCascade( squarkToEwino,
                     2,
                     inputShortcut,
                     initialSquark,
                     beamEnergy ),
        ewinoFlipsCharge( false ),
        ewinoCascades( ewinoCascades )
    {
      cascadeSegment = cascadeDefiner.addNewAtEnd();
      cascadeSegment->first = ewinoCascade->getElectroweakDecayer();
      // squarkToEwino means only 1 decay to be recorded.
      cascadeSegment->second = firstDecayBodyNumber;
      // squarkToEwino also means that the decay is 2-body.
      if( inputShortcut->isIn( cascadeSegment->first->get_PDG_code(),
                               inputShortcut->getCharginos() ) )
        // if we have to worry about which sign of PDG code to use...
      {
        ewinoFlipsCharge = true;
        soughtDecayProductList.front()
        = -(cascadeSegment->first->get_PDG_code());
      }
      else
      {
        ewinoFlipsCharge = false;
        soughtDecayProductList.front() = cascadeSegment->first->get_PDG_code();
      }
    }

    sdownDirectlyToElectroweak::~sdownDirectlyToElectroweak()
    {
      // does nothing.
    }



    supDirectlyToElectroweak::supDirectlyToElectroweak(
                                       inputHandler const* const inputShortcut,
                                           particlePointer const initialSquark,
                                                       double const beamEnergy,
                              electroweakCascadeSet* const directEwinoCascades,
                               electroweakCascadeSet* const ewinoWithWCascades,
                                 electroweakCascadeSet* const bosonCascades ) :
        fullCascade( squarkToEwino,
                     2,
                     inputShortcut,
                     initialSquark,
                     beamEnergy ),
        directEwinoCascades( directEwinoCascades ),
        ewinoWithWCascades( ewinoWithWCascades ),
        bosonCascades( bosonCascades ),
        wBoson( inputShortcut->getWPlus() ),
        ewinoMass( CppSLHA::CppSLHA_global::really_wrong_value ),
        shouldUseDecaysWithW( false ),
        directFraction( CppSLHA::CppSLHA_global::really_wrong_value ),
        wFraction( CppSLHA::CppSLHA_global::really_wrong_value ),
        decayProductListIncludingW()
    {
      cascadeSegment = cascadeDefiner.addNewAtEnd();
      cascadeSegment->first = ewinoCascade->getElectroweakDecayer();
      // squarkToEwino means only 1 decay to be recorded.
      cascadeSegment->second = firstDecayBodyNumber;
      // squarkToEwino also means that the decay is 2-body.
      soughtDecayProductList.front() = cascadeSegment->first->get_PDG_code();
      decayProductListIncludingW.push_back( soughtDecayProductList.front() );
      decayProductListIncludingW.push_back(
                 ewinoWithWCascades->getElectroweakDecayer()->get_PDG_code() );
    }

    supDirectlyToElectroweak::~supDirectlyToElectroweak()
    {
      // does nothing.
    }


    double
    supDirectlyToElectroweak::getAcceptance(
                                  bool const initialSparticleIsNotAntiparticle,
                                        acceptanceCutSet* const acceptanceCuts,
                                             int const numberOfAdditionalJets,
                                           int const numberOfNegativeElectrons,
                                           int const numberOfPositiveElectrons,
                                             int const numberOfNegativeMuons,
                                             int const numberOfPositiveMuons )
    /* this calls the appropriate functions on directEwinoCascades to build the
     * required acceptance or ewinoWithWCascades & bosonCascades if the decay
     * with a boson is not negligible, taking into account whether the charges
     * should be swapped if scoloredIsNotAntiparticle is false.
     */
    {
      if( !initialSparticleIsNotAntiparticle )
        // if the electroweakino is an antiparticle, we swap the charges.
      {
        int swappingInt( numberOfPositiveElectrons );
        numberOfPositiveElectrons = numberOfNegativeElectrons;
        numberOfNegativeElectrons = swappingInt;
        swappingInt = numberOfPositiveMuons;
        numberOfPositiveMuons = numberOfNegativeMuons;
        numberOfNegativeMuons = swappingInt;
      }
      if( shouldUseDecaysWithW )
      {
        return getCombinedAcceptance( acceptanceCuts,
                                      numberOfAdditionalJets,
                                      numberOfNegativeElectrons,
                                      numberOfPositiveElectrons,
                                      numberOfNegativeMuons,
                                      numberOfPositiveMuons );
      }
      else
      {
        return directEwinoCascades->getAcceptance( acceptanceCuts,
                                                   numberOfAdditionalJets,
                                                   numberOfNegativeElectrons,
                                                   numberOfPositiveElectrons,
                                                   numberOfNegativeMuons,
                                                   numberOfPositiveMuons );
      }
    }

    double
    supDirectlyToElectroweak::getCombinedAcceptance(
                                        acceptanceCutSet* const acceptanceCuts,
                                              int const numberOfAdditionalJets,
                                           int const numberOfNegativeElectrons,
                                           int const numberOfPositiveElectrons,
                                               int const numberOfNegativeMuons,
                                              int const numberOfPositiveMuons )
    {
      if( ( 0 > numberOfAdditionalJets )
          ||
          ( 0 > numberOfNegativeElectrons )
          ||
          ( 0 > numberOfPositiveElectrons )
          ||
          ( 0 > numberOfNegativeMuons )
          ||
          ( 0 > numberOfPositiveMuons )
          ||
          ( numberOfSmFermionsFromElectroweakDecaysPerFullCascade
            < ( numberOfAdditionalJets
                + numberOfNegativeElectrons + numberOfPositiveElectrons
                + numberOfNegativeMuons + numberOfPositiveMuons ) ) )
      {
        return 0.0;
      }
      else
      {
        double returnDouble( 0.0 );
        for( int ewinoJets( numberOfAdditionalJets );
             0 <= ewinoJets;
             --ewinoJets )
        {
          for( int ewinoNegativeElectrons( numberOfNegativeElectrons );
               0 <= ewinoNegativeElectrons;
               --ewinoNegativeElectrons )
          {
            for( int ewinoPositiveElectrons( numberOfPositiveElectrons );
                 0 <= ewinoPositiveElectrons;
                 --ewinoPositiveElectrons )
            {
              for( int ewinoNegativeMuons( numberOfNegativeMuons );
                   0 <= ewinoNegativeMuons;
                   --ewinoNegativeMuons )
              {
                for( int ewinoPositiveMuons( numberOfPositiveMuons );
                     0 <= ewinoPositiveMuons;
                     --ewinoPositiveMuons )
                {
                  returnDouble
                  += ( ewinoWithWCascades->getAcceptance( acceptanceCuts,
                                                          ewinoJets,
                                                        ewinoNegativeElectrons,
                                                        ewinoPositiveElectrons,
                                                          ewinoNegativeMuons,
                                                          ewinoPositiveMuons )
                       * bosonCascades->getAcceptance( acceptanceCuts,
                                        ( numberOfAdditionalJets - ewinoJets ),
                        ( numberOfNegativeElectrons - ewinoNegativeElectrons ),
                        ( numberOfPositiveElectrons - ewinoPositiveElectrons ),
                                ( numberOfNegativeMuons - ewinoNegativeMuons ),
                            ( numberOfPositiveMuons - ewinoPositiveMuons ) ) );
                }  // end of loop over positive electrons.
              }  // end of loop over negative muons.
            }  // end of loop over positive electrons.
          }  // end of loop over negative electrons.
        }  // end of loop over jets.
        return ( wFraction * returnDouble
                 + directFraction
                   * directEwinoCascades->getAcceptance( acceptanceCuts,
                                                        numberOfAdditionalJets,
                                                     numberOfNegativeElectrons,
                                                     numberOfPositiveElectrons,
                                                         numberOfNegativeMuons,
                                                     numberOfPositiveMuons ) );
      }  // end of if number of SM fermions was in the allowed range.
    }



    gluinoDirectlyToElectroweak::gluinoDirectlyToElectroweak() :
        fullCascade( gluinoToEwino,
                     3,
                     inputShortcut,
                     inputShortcut->getGluino(),
                     beamEnergy ),
        ewinoCascades( ewinoCascades )
    {
      cascadeSegment = cascadeDefiner.addNewAtEnd();
      cascadeSegment->first = ewinoCascade->getElectroweakDecayer();
      // gluinoToEwino means only 1 decay to be recorded.
      cascadeSegment->second = firstDecayBodyNumber;
      // gluinoToEwino also means that the decay is 3-body.
      soughtDecayProductList.front() = cascadeSegment->first->get_PDG_code();
      decayProductListIncludingW.push_back( soughtDecayProductList.front() );
      decayProductListIncludingW.push_back(
                 ewinoWithWCascades->getElectroweakDecayer()->get_PDG_code() );
    }

    gluinoDirectlyToElectroweak::~gluinoDirectlyToElectroweak()
    {
      // does nothing.
    }


    double
    gluinoDirectlyToElectroweak::getAcceptance(
                                  bool const initialSparticleIsNotAntiparticle,
                                        acceptanceCutSet* const acceptanceCuts,
                                              int const numberOfAdditionalJets,
                                                 int numberOfNegativeElectrons,
                                                 int numberOfPositiveElectrons,
                                                     int numberOfNegativeMuons,
                                                    int numberOfPositiveMuons )
    /* this calls the appropriate functions on ewinoCascade, summing over the 2
     * charge versions if the electroweakino is a chargino. it ignores the
     * value of scoloredIsNotAntiparticle.
     */
    {
      if( ewinoCascades->getElectroweakDecayer()->counts_as_self_conjugate() )
        // if the gluino has only 1 charge version of this decay...
      {
        return ewinoCascades->getAcceptance( acceptanceCuts,
                                             numberOfAdditionalJets,
                                             numberOfNegativeElectrons,
                                             numberOfPositiveElectrons,
                                             numberOfNegativeMuons,
                                             numberOfPositiveMuons );
        // the summing over charge-conjugate decays is already included in the
        // BR to the electroweakino.
      }
      else
        // otherwise if the gluino has 2 charge versions of this decay...
      {
        return ( 0.5 * ( ewinoCascades->getAcceptance( acceptanceCuts,
                                                       numberOfAdditionalJets,
                                                       numberOfNegativeElectrons,
                                                       numberOfPositiveElectrons,
                                                       numberOfNegativeMuons,
                                                       numberOfPositiveMuons )
                         + ewinoCascades->getAcceptance( acceptanceCuts,
                                                        numberOfAdditionalJets,
                                                     numberOfPositiveElectrons,
                                                     numberOfNegativeElectrons,
                                                         numberOfPositiveMuons,
                                                   numberOfNegativeMuons ) ) );
        // the summing over charge-conjugate decays is already included in the
        // BR to the electroweakino.
      }
    }



    squarkByBosonToCompound::squarkByBosonToCompound() :
        fullCascade( squarkToSquarkThenMore,
                     2 ),
        bosonCascades( NULL ),
        bosonFlipsCharge(false )
    // squarkToSquarkThenMore also means that the initial decay is 2-body.
    {
      soughtDecayProductList.push_back(
                                 CppSLHA::CppSLHA_global::really_wrong_value );
      // a cascade beginning with a squark decaying to another squark plus a
      // boson needs to look for 2 specific decay products.
    }

    squarkByBosonToCompound::~squarkByBosonToCompound()
    {
      // does nothing.
    }


    double
    squarkByBosonToCompound::getAcceptance(
                                  bool const initialSparticleIsNotAntiparticle,
                                        acceptanceCutSet* const acceptanceCuts,
                                            int const numberOfAdditionalJets,
                                            int numberOfNegativeElectrons,
                                            int numberOfPositiveElectrons,
                                            int numberOfNegativeMuons,
                                            int numberOfPositiveMuons )
    /* this calls the appropriate functions on subcascadePointer &
     * bosonCascades to build the required acceptance, taking into account
     * whether the charges should be swapped if scoloredIsNotAntiparticle is
     * false.
     */
    {
      if( ( 0 > numberOfAdditionalJets )
          ||
          ( 0 > numberOfNegativeElectrons )
          ||
          ( 0 > numberOfPositiveElectrons )
          ||
          ( 0 > numberOfNegativeMuons )
          ||
          ( 0 > numberOfPositiveMuons )
          ||
          ( numberOfSmFermionsFromElectroweakDecaysPerFullCascade
            < ( numberOfAdditionalJets
                + numberOfNegativeElectrons + numberOfPositiveElectrons
                + numberOfNegativeMuons + numberOfPositiveMuons ) ) )
      {
        return 0.0;
      }
      else
      {
        if( ( initialSparticleIsNotAntiparticle
              &&
              bosonFlipsCharge )
            ||
            ( !initialSparticleIsNotAntiparticle
              &&
              !bosonFlipsCharge ) )
          // if the electroweakino is an antiparticle, we swap the charges.
        {
          int swappingInt( numberOfPositiveElectrons );
          numberOfPositiveElectrons = numberOfNegativeElectrons;
          numberOfNegativeElectrons = swappingInt;
          swappingInt = numberOfPositiveMuons;
          numberOfPositiveMuons = numberOfNegativeMuons;
          numberOfNegativeMuons = swappingInt;
        }
        double returnDouble( 0.0 );
        for( int ewinoJets( numberOfAdditionalJets );
             0 <= ewinoJets;
             --ewinoJets )
        {
          for( int ewinoNegativeElectrons( numberOfNegativeElectrons );
               0 <= ewinoNegativeElectrons;
               --ewinoNegativeElectrons )
          {
            for( int ewinoPositiveElectrons( numberOfPositiveElectrons );
                 0 <= ewinoPositiveElectrons;
                 --ewinoPositiveElectrons )
            {
              for( int ewinoNegativeMuons( numberOfNegativeMuons );
                   0 <= ewinoNegativeMuons;
                   --ewinoNegativeMuons )
              {
                for( int ewinoPositiveMuons( numberOfPositiveMuons );
                     0 <= ewinoPositiveMuons;
                     --ewinoPositiveMuons )
                {
                returnDouble
                += ( bosonCascades->getAcceptance( acceptanceCuts,
                                        ( numberOfAdditionalJets - ewinoJets ),
                        ( numberOfNegativeElectrons - ewinoNegativeElectrons ),
                        ( numberOfPositiveElectrons - ewinoPositiveElectrons ),
                                ( numberOfNegativeMuons - ewinoNegativeMuons ),
                               ( numberOfPositiveMuons - ewinoPositiveMuons ) )
                     * subcascadePointer->getAcceptance(
                                             initialSparticleIsNotAntiparticle,
                                                         acceptanceCuts,
                                                         ewinoJets,
                                                        ewinoNegativeElectrons,
                                                        ewinoPositiveElectrons,
                                                         ewinoNegativeMuons,
                                                        ewinoPositiveMuons ) );
                // in the used conventions, an incoming squark decays to a
                // boson plus a squark, hence the same bool is used for both.
                }  // end of loop over positive electrons.
              }  // end of loop over negative muons.
            }  // end of loop over positive electrons.
          }  // end of loop over negative electrons.
        }  // end of loop over jets.
        return returnDouble;
      }  // end of if number of SM fermions was in the allowed range.
    }



    sdownByJetToCompound::sdownByJetToCompound() :
        fullCascade( sjm,
                     2 )
        // sjm also means that the initial decay is 2-body.
    {
      // just an initialization list.
    }

    sdownByJetToCompound::~sdownByJetToCompound()
    {
      // does nothing.
    }



    supByJetToCompound::supByJetToCompound() :
        fullCascade( sjm,
                     2 ),
        bosonCascades( NULL ),
        wBoson( NULL ),
        ewinoMass( CppSLHA::CppSLHA_global::really_wrong_value ),
        shouldUseDecaysWithW( false ),
        directFraction( CppSLHA::CppSLHA_global::really_wrong_value ),
        wFraction( CppSLHA::CppSLHA_global::really_wrong_value ),
        decayProductListIncludingW( 2,
                                  CppSLHA::CppSLHA_global::really_wrong_value )
        // sjm also means that the initial decay is 2-body.
    {
      // just an initialization list.
    }

    supByJetToCompound::~supByJetToCompound()
    {
      // does nothing.
    }


    double
    supByJetToCompound::getAcceptance(
                                  bool const initialSparticleIsNotAntiparticle,
                                       acceptanceCutSet* const acceptanceCuts,
                                       int const numberOfAdditionalJets,
                                       int numberOfNegativeElectrons,
                                       int numberOfPositiveElectrons,
                                       int numberOfNegativeMuons,
                                       int numberOfPositiveMuons )
    /* this calls the appropriate functions on subcascadePointer to build the
     * required acceptance, & also with bosonCascades if the decay with a boson
     * is not negligible, taking into account whether the charges should be
     * swapped if scoloredIsNotAntiparticle is false.
     */
    {
      if( !initialSparticleIsNotAntiparticle )
        // if the electroweakino is an antiparticle, we swap the charges.
      {
        int swappingInt( numberOfPositiveElectrons );
        numberOfPositiveElectrons = numberOfNegativeElectrons;
        numberOfNegativeElectrons = swappingInt;
        swappingInt = numberOfPositiveMuons;
        numberOfPositiveMuons = numberOfNegativeMuons;
        numberOfNegativeMuons = swappingInt;
      }
      if( shouldUseDecaysWithW )
      {
        return getCombinedAcceptance( acceptanceCuts,
                                      numberOfAdditionalJets,
                                      numberOfNegativeElectrons,
                                      numberOfPositiveElectrons,
                                      numberOfNegativeMuons,
                                      numberOfPositiveMuons );
      }
      else
      {
        return
        subcascadePointer->getAcceptance( initialSparticleIsNotAntiparticle,
                                          acceptanceCuts,
                                          numberOfAdditionalJets,
                                          numberOfNegativeElectrons,
                                          numberOfPositiveElectrons,
                                          numberOfNegativeMuons,
                                          numberOfPositiveMuons );
      }
    }

    double
    supByJetToCompound::getCombinedAcceptance(
                                        acceptanceCutSet* const acceptanceCuts,
                                              int const numberOfAdditionalJets,
                                           int const numberOfNegativeElectrons,
                                           int const numberOfPositiveElectrons,
                                               int const numberOfNegativeMuons,
                                              int const numberOfPositiveMuons )
    {
      if( ( 0 > numberOfAdditionalJets )
          ||
          ( 0 > numberOfNegativeElectrons )
          ||
          ( 0 > numberOfPositiveElectrons )
          ||
          ( 0 > numberOfNegativeMuons )
          ||
          ( 0 > numberOfPositiveMuons )
          ||
          ( numberOfSmFermionsFromElectroweakDecaysPerFullCascade
            < ( numberOfAdditionalJets
                + numberOfNegativeElectrons + numberOfPositiveElectrons
                + numberOfNegativeMuons + numberOfPositiveMuons ) ) )
      {
        return 0.0;
      }
      else
      {
        double returnDouble( 0.0 );
        for( int ewinoJets( numberOfAdditionalJets );
             0 <= ewinoJets;
             --ewinoJets )
        {
          for( int ewinoNegativeElectrons( numberOfNegativeElectrons );
               0 <= ewinoNegativeElectrons;
               --ewinoNegativeElectrons )
          {
            for( int ewinoPositiveElectrons( numberOfPositiveElectrons );
                 0 <= ewinoPositiveElectrons;
                 --ewinoPositiveElectrons )
            {
              for( int ewinoNegativeMuons( numberOfNegativeMuons );
                   0 <= ewinoNegativeMuons;
                   --ewinoNegativeMuons )
              {
                for( int ewinoPositiveMuons( numberOfPositiveMuons );
                     0 <= ewinoPositiveMuons;
                     --ewinoPositiveMuons )
                {
                  returnDouble
                  += ( subcascadePointer->getAcceptance( true,
                                                         acceptanceCuts,
                                                         ewinoJets,
                                                        ewinoNegativeElectrons,
                                                        ewinoPositiveElectrons,
                                                         ewinoNegativeMuons,
                                                         ewinoPositiveMuons )
                       * bosonCascades->getAcceptance( acceptanceCuts,
                                        ( numberOfAdditionalJets - ewinoJets ),
                        ( numberOfNegativeElectrons - ewinoNegativeElectrons ),
                        ( numberOfPositiveElectrons - ewinoPositiveElectrons ),
                                ( numberOfNegativeMuons - ewinoNegativeMuons ),
                            ( numberOfPositiveMuons - ewinoPositiveMuons ) ) );
                }  // end of loop over positive electrons.
              }  // end of loop over negative muons.
            }  // end of loop over positive electrons.
          }  // end of loop over negative electrons.
        }  // end of loop over jets.
        return ( wFraction * returnDouble
                 + directFraction
                   * directEwinoCascades->getAcceptance( acceptanceCuts,
                                                        numberOfAdditionalJets,
                                                     numberOfNegativeElectrons,
                                                     numberOfPositiveElectrons,
                                                         numberOfNegativeMuons,
                                                     numberOfPositiveMuons ) );
      }  // end of if number of SM fermions was in the allowed range.
    }



    gluinoOrNeutralinoToCompound::gluinoOrNeutralinoToCompound() :
        fullCascade( gjm,
                     2 ),
        bosonCascades( NULL ),
        wBoson( NULL ),
        ewinoMass( CppSLHA::CppSLHA_global::really_wrong_value ),
        shouldUseDecaysWithW( false ),
        directFraction( CppSLHA::CppSLHA_global::really_wrong_value ),
        wFraction( CppSLHA::CppSLHA_global::really_wrong_value ),
        decayProductListIncludingW( 2,
                                  CppSLHA::CppSLHA_global::really_wrong_value )
        // gjm also means that the initial decay is 2-body.
    {
      // just an initialization list.
    }

    gluinoOrNeutralinoToCompound::~gluinoOrNeutralinoToCompound()
    {
      // does nothing.
    }


    double
    gluinoOrNeutralinoToCompound::getAcceptance(
                                  bool const initialSparticleIsNotAntiparticle,
                                        acceptanceCutSet* const acceptanceCuts,
                                              int const numberOfAdditionalJets,
                                                 int numberOfNegativeElectrons,
                                                 int numberOfPositiveElectrons,
                                                 int numberOfNegativeMuons,
                                                 int numberOfPositiveMuons )
    /* this calls the appropriate functions on subcascadePointer to build the
     * required acceptance, & also with bosonCascades if the decay with a boson
     * is not negligible, taking into account whether the charges should be
     * swapped if scoloredIsNotAntiparticle is false.
     */
    {
      if( shouldUseDecaysWithW )
      {
        return getCombinedAcceptance( acceptanceCuts,
                                      numberOfAdditionalJets,
                                      numberOfNegativeElectrons,
                                      numberOfPositiveElectrons,
                                      numberOfNegativeMuons,
                                      numberOfPositiveMuons );
      }
      else
      {
        return
        ( 0.5 * ( subcascadePointer->getAcceptance(
                                             initialSparticleIsNotAntiparticle,
                                                    acceptanceCuts,
                                                    numberOfAdditionalJets,
                                                    numberOfNegativeElectrons,
                                                    numberOfPositiveElectrons,
                                                    numberOfNegativeMuons,
                                                    numberOfPositiveMuons )
                  + subcascadePointer->getAcceptance(
                                             initialSparticleIsNotAntiparticle,
                                                      acceptanceCuts,
                                                      numberOfAdditionalJets,
                                                     numberOfNegativeElectrons,
                                                     numberOfPositiveElectrons,
                                                      numberOfNegativeMuons,
                                                   numberOfPositiveMuons ) ) );
      }
    }

    double
    gluinoOrNeutralinoToCompound::getCombinedAcceptance(
                                        acceptanceCutSet* const acceptanceCuts,
                                              int const numberOfAdditionalJets,
                                           int const numberOfNegativeElectrons,
                                           int const numberOfPositiveElectrons,
                                               int const numberOfNegativeMuons,
                                              int const numberOfPositiveMuons )
    {
      if( ( 0 > numberOfAdditionalJets )
          ||
          ( 0 > numberOfNegativeElectrons )
          ||
          ( 0 > numberOfPositiveElectrons )
          ||
          ( 0 > numberOfNegativeMuons )
          ||
          ( 0 > numberOfPositiveMuons )
          ||
          ( numberOfSmFermionsFromElectroweakDecaysPerFullCascade
            < ( numberOfAdditionalJets
                + numberOfNegativeElectrons + numberOfPositiveElectrons
                + numberOfNegativeMuons + numberOfPositiveMuons ) ) )
      {
        return 0.0;
      }
      else
      {
        double returnDouble( 0.0 );
        for( int ewinoJets( numberOfAdditionalJets );
             0 <= ewinoJets;
             --ewinoJets )
        {
          for( int ewinoNegativeElectrons( numberOfNegativeElectrons );
               0 <= ewinoNegativeElectrons;
               --ewinoNegativeElectrons )
          {
            for( int ewinoPositiveElectrons( numberOfPositiveElectrons );
                 0 <= ewinoPositiveElectrons;
                 --ewinoPositiveElectrons )
            {
              for( int ewinoNegativeMuons( numberOfNegativeMuons );
                   0 <= ewinoNegativeMuons;
                   --ewinoNegativeMuons )
              {
                for( int ewinoPositiveMuons( numberOfPositiveMuons );
                     0 <= ewinoPositiveMuons;
                     --ewinoPositiveMuons )
                {
                  /* we average over the configurations with opposite charges
                   * (could use false for initialSparticleIsNotAntiparticle
                   * for subcascadePointer & appropriate negative/positive
                   * numbers, but this way seems simpler):
                   */
                  returnDouble
                  += ( 0.5 * ( subcascadePointer->getAcceptance( true,
                                                                acceptanceCuts,
                                                                ewinoJets,
                                                        ewinoNegativeElectrons,
                                                        ewinoPositiveElectrons,
                                                            ewinoNegativeMuons,
                                                           ewinoPositiveMuons )
                               * bosonCascades->getAcceptance( acceptanceCuts,
                                        ( numberOfAdditionalJets - ewinoJets ),
                        ( numberOfNegativeElectrons - ewinoNegativeElectrons ),
                        ( numberOfPositiveElectrons - ewinoPositiveElectrons ),
                                ( numberOfNegativeMuons - ewinoNegativeMuons ),
                               ( numberOfPositiveMuons - ewinoPositiveMuons ) )
                               + subcascadePointer->getAcceptance( true,
                                                                acceptanceCuts,
                                                                   ewinoJets,
                                                        ewinoPositiveElectrons,
                                                        ewinoNegativeElectrons,
                                                            ewinoPositiveMuons,
                                                           ewinoNegativeMuons )
                                 * bosonCascades->getAcceptance(
                                                                acceptanceCuts,
                                        ( numberOfAdditionalJets - ewinoJets ),
                        ( numberOfPositiveElectrons - ewinoPositiveElectrons ),
                        ( numberOfNegativeElectrons - ewinoNegativeElectrons ),
                                ( numberOfPositiveMuons - ewinoPositiveMuons ),
                          ( numberOfNegativeMuons - ewinoNegativeMuons ) ) ) );
                }  // end of loop over positive electrons.
              }  // end of loop over negative muons.
            }  // end of loop over positive electrons.
          }  // end of loop over negative electrons.
        }  // end of loop over jets.
        return ( wFraction * returnDouble
                 + directFraction * 0.5
                   * ( subcascadePointer->getAcceptance( true,
                                                         acceptanceCuts,
                                                        numberOfAdditionalJets,
                                                     numberOfNegativeElectrons,
                                                     numberOfPositiveElectrons,
                                                         numberOfNegativeMuons,
                                                        numberOfPositiveMuons )
                       + subcascadePointer->getAcceptance( true,
                                                           acceptanceCuts,
                                                        numberOfAdditionalJets,
                                                     numberOfPositiveElectrons,
                                                     numberOfNegativeElectrons,
                                                         numberOfPositiveMuons,
                                                   numberOfNegativeMuons ) ) );
      }  // end of if number of SM fermions was in the allowed range.
    }



    charginoToCompound::charginoToCompound() :
        fullCascade( gjm,
                     2 ),
        bosonCascades( NULL ),
        wBoson( NULL ),
        ewinoMass( CppSLHA::CppSLHA_global::really_wrong_value ),
        shouldUseDecaysWithW( false ),
        directFraction( CppSLHA::CppSLHA_global::really_wrong_value ),
        wFraction( CppSLHA::CppSLHA_global::really_wrong_value ),
        decayProductListIncludingW( 2,
                                 CppSLHA::CppSLHA_global::really_wrong_value ),
        squarkWithWFraction( CppSLHA::CppSLHA_global::really_wrong_value ),
        antisquarkWithWFraction( CppSLHA::CppSLHA_global::really_wrong_value ),
        decayingToSupType( false )
        // gjm also means that the initial decay is 2-body.
    {
      // just an initialization list.
    }

    charginoToCompound::~charginoToCompound()
    {
      // does nothing.
    }


    double
    charginoToCompound::getAcceptance(
                                  bool const initialSparticleIsNotAntiparticle,
                                       acceptanceCutSet* const acceptanceCuts,
                                       int const numberOfAdditionalJets,
                                       int numberOfNegativeElectrons,
                                       int numberOfPositiveElectrons,
                                       int numberOfNegativeMuons,
                                       int numberOfPositiveMuons )
    /* this calls the appropriate functions on subcascadePointer to build the
     * required acceptance, & also with bosonCascades if the decay with a boson
     * is not negligible, taking into account whether the charges should be
     * swapped if scoloredIsNotAntiparticle is false.
     */
    {
      if( shouldUseDecaysWithW )
      {
        if( initialSparticleIsNotAntiparticle )
        {
          return getCombinedAcceptance( acceptanceCuts,
                                        numberOfAdditionalJets,
                                        numberOfNegativeElectrons,
                                        numberOfPositiveElectrons,
                                        numberOfNegativeMuons,
                                        numberOfPositiveMuons );
        }
        else
        {
          return getCombinedAcceptance( acceptanceCuts,
                                        numberOfAdditionalJets,
                                        numberOfPositiveElectrons,
                                        numberOfNegativeElectrons,
                                        numberOfPositiveMuons,
                                        numberOfNegativeMuons );
        }
      }
      else
      {
        bool askForSquarkNotAntisquark( !decayingToSupType );
        if( initialSparticleIsNotAntiparticle )
        {
          askForSquarkNotAntisquark = decayingToSupType;
        }
        return subcascadePointer->getAcceptance( askForSquarkNotAntisquark,
                                                 acceptanceCuts,
                                                 numberOfAdditionalJets,
                                                 numberOfNegativeElectrons,
                                                 numberOfPositiveElectrons,
                                                 numberOfNegativeMuons,
                                                 numberOfPositiveMuons );
      }
    }

    double
    charginoToCompound::getCombinedAcceptance(
                                        acceptanceCutSet* const acceptanceCuts,
                                              int const numberOfAdditionalJets,
                                           int const numberOfNegativeElectrons,
                                           int const numberOfPositiveElectrons,
                                               int const numberOfNegativeMuons,
                                              int const numberOfPositiveMuons )
    // if this is called, it is the decay of a positive chargino.
    {
      if( ( 0 > numberOfAdditionalJets )
          ||
          ( 0 > numberOfNegativeElectrons )
          ||
          ( 0 > numberOfPositiveElectrons )
          ||
          ( 0 > numberOfNegativeMuons )
          ||
          ( 0 > numberOfPositiveMuons )
          ||
          ( numberOfSmFermionsFromElectroweakDecaysPerFullCascade
            < ( numberOfAdditionalJets
                + numberOfNegativeElectrons + numberOfPositiveElectrons
                + numberOfNegativeMuons + numberOfPositiveMuons ) ) )
      {
        return 0.0;
      }
      else
      {
        double returnDouble( 0.0 );
        for( int ewinoJets( numberOfAdditionalJets );
             0 <= ewinoJets;
             --ewinoJets )
        {
          for( int ewinoNegativeElectrons( numberOfNegativeElectrons );
               0 <= ewinoNegativeElectrons;
               --ewinoNegativeElectrons )
          {
            for( int ewinoPositiveElectrons( numberOfPositiveElectrons );
                 0 <= ewinoPositiveElectrons;
                 --ewinoPositiveElectrons )
            {
              for( int ewinoNegativeMuons( numberOfNegativeMuons );
                   0 <= ewinoNegativeMuons;
                   --ewinoNegativeMuons )
              {
                for( int ewinoPositiveMuons( numberOfPositiveMuons );
                     0 <= ewinoPositiveMuons;
                     --ewinoPositiveMuons )
                {
                  returnDouble
                  += ( ( squarkWithWFraction
                         * subcascadePointer->getAcceptance( false,
                                                             acceptanceCuts,
                                                             ewinoJets,
                                                        ewinoNegativeElectrons,
                                                        ewinoPositiveElectrons,
                                                            ewinoNegativeMuons,
                                                           ewinoPositiveMuons )
                         + antisquarkWithWFraction
                           * subcascadePointer->getAcceptance( false,
                                                               acceptanceCuts,
                                                               ewinoJets,
                                                        ewinoNegativeElectrons,
                                                        ewinoPositiveElectrons,
                                                            ewinoNegativeMuons,
                                                         ewinoPositiveMuons ) )
                       * bosonCascades->getAcceptance( acceptanceCuts,
                                        ( numberOfAdditionalJets - ewinoJets ),
                        ( numberOfNegativeElectrons - ewinoNegativeElectrons ),
                        ( numberOfPositiveElectrons - ewinoPositiveElectrons ),
                                ( numberOfNegativeMuons - ewinoNegativeMuons ),
                            ( numberOfPositiveMuons - ewinoPositiveMuons ) ) );
                }  // end of loop over positive electrons.
              }  // end of loop over negative muons.
            }  // end of loop over positive electrons.
          }  // end of loop over negative electrons.
        }  // end of loop over jets.
        return ( returnDouble
                 + directFraction
                   * subcascadePointer->getAcceptance( decayingToSupType,
                                                       acceptanceCuts,
                                                       numberOfAdditionalJets,
                                                     numberOfNegativeElectrons,
                                                     numberOfPositiveElectrons,
                                                         numberOfNegativeMuons,
                                                     numberOfPositiveMuons ) );
      }  // end of if number of SM fermions was in the allowed range.
    }

  }  // end of fullCascadeType namespace



  fullCascadeSet::fullCascadeSet( inputHandler const* const inputShortcut,
                                  particlePointer const initialScolored,
           electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource,
                                  double const beamEnergy,
                                  fullCascadeSet* const gluinoFullCascade ) :
    getsReadiedForNewPoint( inputShortcut->getReadier() ),
    inputShortcut( inputShortcut ),
    initialSparticle( initialScolored ),
    electroweakCascadeSource( electroweakCascadeSource ),
    openCascades(),
    gluinoFullCascade( gluinoFullCascade ),
    orderedCascadeSets( NULL ),
    potentialSubcascades( NULL ),
    subcascadeBranchingRatio( CppSLHA::CppSLHA_global::really_wrong_value ),
    beamEnergy( beamEnergy ),
    singleSpecifiedDecayProductList( 1,
                                 CppSLHA::CppSLHA_global::really_wrong_value ),
    singleSpecifiedDecayProduct( CppSLHA::CppSLHA_global::really_wrong_value ),
    potentialDecayProducts( NULL )
  {
    // just an initialization list.
  }

  fullCascadeSet::~fullCascadeSet()
  {
    // does nothing.
  }

  void
  fullCascadeSet::setUpCascades()
  // this clears openCascades, then sets it to be filled with all open
  // fullCascades that initialSparticle has for this point.
  {
    openCascades.clear();
    initialSparticleMass = initialSparticle->get_absolute_mass();
    buildLongerCascades();
  }



  fullCascadeSetOrderer::fullCascadeSetOrderer(
                                       inputHandler const* const inputShortcut,
                                    fullCascadeSet* const gluinoFullCascade ) :
      getsReadiedForNewPoint( inputShortcut->getReadier() ),
      inputShortcut( inputShortcut ),
      gluinoFullCascade( gluinoFullCascade ),
      sdownTypeColoredCascades(),
      supTypeColoredCascades(),
      allNeutralinoColoredCascades(),
      allCharginoColoredCascades(),
      relevantNeutralinoColoredCascades(),
      relevantCharginoColoredCascades(),
      lightestSquarkMass( CppSLHA::CppSLHA_global::really_wrong_value )
  {
    // just an initialization list.
  }

  fullCascadeSetOrderer::~fullCascadeSetOrderer()
  {
    // does nothing.
  }


  void
  fullCascadeSetOrderer::orderCascades()
  {
    sdownTypeColoredCascades.sort( &(fullCascadeSet::massOrdered) );
    lightestSquarkMass
    = sdownTypeColoredCascades.front()->getInitialSparticle(
                                                        )->get_absolute_mass();
    supTypeColoredCascades.sort( &(fullCascadeSet::massOrdered) );
    if( supTypeColoredCascades.front()->getInitialSparticle(
                                                         )->get_absolute_mass()
        < lightestSquarkMass )
    {
      lightestSquarkMass
      = supTypeColoredCascades.front()->getInitialSparticle(
                                                        )->get_absolute_mass();
    }
    allNeutralinoColoredCascades.sort( &(fullCascadeSet::massOrdered) );
    relevantNeutralinoColoredCascades.clear();
    for( std::list< fullCascadeSet* >::reverse_iterator
         setIterator( allNeutralinoColoredCascades.rbegin() );
         allNeutralinoColoredCascades.rend() != setIterator;
         ++setIterator )
    {
      if( (*setIterator)->getInitialSparticle()->get_absolute_mass()
          > lightestSquarkMass )
      {
        relevantNeutralinoColoredCascades.push_front( *setIterator );
      }
    }
    allCharginoColoredCascades.sort( &(fullCascadeSet::massOrdered) );
    relevantCharginoColoredCascades.clear();
    for( std::list< fullCascadeSet* >::reverse_iterator
         setIterator( allCharginoColoredCascades.rbegin() );
         allCharginoColoredCascades.rend() != setIterator;
         ++setIterator )
    {
      if( (*setIterator)->getInitialSparticle()->get_absolute_mass()
          > lightestSquarkMass )
      {
        relevantCharginoColoredCascades.push_front( *setIterator );
      }
    }
  }


  namespace fullCascadeSetType
  {
    sdownTypeSet::sdownTypeSet( inputHandler const* const inputShortcut,
           electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource,
                                particlePointer const initialScolored,
                                fullCascadeSetOrderer* const setOrderer,
                                double const beamEnergy,
                                fullCascadeSet* const gluinoFullCascade ) :
    fullCascadeSet( inputShortcut,
                    initialScolored,
                    electroweakCascadeSource,
                    beamEnergy ),
    setOrderer( setOrderer ),
    directToEwinoCascades(),
    compoundByBosonCascades(),
    compoundByJetCascades()
    {
      // we have to set up the cascades directly to electroweakinos now:
      for( std::vector< particlePointer >::const_iterator
           ewinoIterator( inputShortcut->getElectroweakinos()->begin() );
           inputShortcut->getElectroweakinos()->end() > ewinoIterator;
           ++ewinoIterator )
      {
        directToEwinoCascades.push_back(
                               new fullCascadeType::sdownDirectlyToElectroweak(
                                                                 inputShortcut,
                                                               initialScolored,
                                                                    beamEnergy,
                                       electroweakCascadeSource->getCascadeSet(
                                                               initialScolored,
                                                          *ewinoIterator ) ) );
      }
    }

    sdownTypeSet::~sdownTypeSet()
    {
      for( std::vector< fullCascadeType::sdownDirectlyToElectroweak*
                                                                    >::iterator
           deletionIterator( directToEwinoCascades.begin() );
           directToEwinoCascades.end() > deletionIterator;
           ++deletionIterator )
      {
        delete *deletionIterator;
      }
    }


    void
    sdownTypeSet::buildLongerCascades()
    {
      // 1st we clear the compound cascades:
      compoundByBosonCascades.clearEntries();
      compoundByJetCascades.clearEntries();

      // next we check to see if we should add compound cascades from the
      // gluino:
      if( inputShortcut->getGluinoMass() < initialSparticleMass )
      {
        singleSpecifiedDecayProductList.front() = CppSLHA::PDG_code::gluino;
        subcascadeBranchingRatio
        = initialSparticle->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                            &singleSpecifiedDecayProductList );
        if( lhcFaserGlobal::negligibleBr < subcascadeBranchingRatio )
          // if the BR to the next stage of the cascade is not negligible...
        {
          potentialSubcascades = gluinoFullCascade->getOpenCascades();
          for( std::vector< fullCascade* >::iterator
               cascadeIterator( potentialSubcascades->begin() );
               potentialSubcascades->end() > cascadeIterator;
               ++cascadeIterator )
          {
            if( lhcFaserGlobal::negligibleBr
                < ( subcascadeBranchingRatio
                    * (*cascadeIterator)->getBrToEwino() ) )
            {
              // we add each cascade with an overall BR that is not negligible:
              openCascades.push_back(
                            compoundByJetCascades.addNewAtEnd()->setProperties(
                                                              initialSparticle,
                                                          *cascadeIterator ) );
            }
          }
        }
      }

      // now we look at compound cascades from neutralinos:
      orderedCascadeSets = setOrderer->getNeutralinoColoredCascades();
      std::list< fullCascadeSet* >::iterator
      setIterator( orderedCascadeSets->begin() );
      while( ( setIterator != orderedCascadeSets->end() )
             &&
             ( initialSparticle->get_absolute_mass()
               > (*setIterator)->getInitialSparticle()->get_absolute_mass() ) )
        /* we go through the list of neutralinos heavier than the lightest
         * squark until we get to either the end of the list or a neutralino
         * heavier than this down-type squark:
         */
      {
        singleSpecifiedDecayProductList.front()
        = (*setIterator)->getInitialSparticle()->get_PDG_code();
        subcascadeBranchingRatio
        = initialSparticle->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                            &singleSpecifiedDecayProductList );
        if( lhcFaserGlobal::negligibleBr < subcascadeBranchingRatio )
          // if the BR to the next stage of the cascade is not negligible...
        {
          potentialSubcascades = (*setIterator)->getOpenCascades();
          for( std::vector< fullCascade* >::iterator
               cascadeIterator( potentialSubcascades->begin() );
               potentialSubcascades->end() > cascadeIterator;
               ++cascadeIterator )
          {
            if( lhcFaserGlobal::negligibleBr
                < ( subcascadeBranchingRatio
                    * (*cascadeIterator)->getBrToEwino() ) )
            {
              // we add each cascade with an overall BR that is not negligible:
              openCascades.push_back(
                            compoundByJetCascades.addNewAtEnd()->setProperties(
                                                              initialSparticle,
                                                          *cascadeIterator ) );
            }
          }
        }
        // now we move on to the next potential neutralino cascade:
        ++setIterator;
      }

      // now we look at compound cascades from charginos:
      orderedCascadeSets = setOrderer->getCharginoColoredCascades();
      setIterator = orderedCascadeSets->begin();
      while( ( setIterator != orderedCascadeSets->end() )
             &&
             ( initialSparticle->get_absolute_mass()
               > (*setIterator)->getInitialSparticle()->get_absolute_mass() ) )
        /* we go through the list of charginos heavier than the lightest squark
         * until we get to either the end of the list or a chargino heavier
         * than this down-type squark:
         */
      {
        singleSpecifiedDecayProductList.front()
        = -((*setIterator)->getInitialSparticle()->get_PDG_code());
        // down-type squarks decay to negatively-charged charginos.
        subcascadeBranchingRatio
        = initialSparticle->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                            &singleSpecifiedDecayProductList );
        if( lhcFaserGlobal::negligibleBr < subcascadeBranchingRatio )
          // if the BR to the next stage of the cascade is not negligible...
        {
          potentialSubcascades = (*setIterator)->getOpenCascades();
          for( std::vector< fullCascade* >::iterator
               cascadeIterator( potentialSubcascades->begin() );
               potentialSubcascades->end() > cascadeIterator;
               ++cascadeIterator )
          {
            if( lhcFaserGlobal::negligibleBr
                < ( subcascadeBranchingRatio
                    * (*cascadeIterator)->getBrToEwino() ) )
            {
              // we add each cascade with an overall BR that is not negligible:
              openCascades.push_back(
                            compoundByJetCascades.addNewAtEnd()->setProperties(
                                                              initialSparticle,
                                                          *cascadeIterator ) );
            }
          }
        }
        // now we move on to the next potential chargino cascade:
        ++setIterator;
      }

      // now we look at compound cascades from down-type squarks:
      orderedCascadeSets = setOrderer->getSdownTypeCascades();
      setIterator = orderedCascadeSets->begin();
      while( ( setIterator != orderedCascadeSets->end() )
             &&
             ( initialSparticle->get_absolute_mass()
               > (*setIterator)->getInitialSparticle()->get_absolute_mass() ) )
        /* we go through the list of down-type squarks until we get to either
         * the end of the list or run out of down-type squarks lighter than
         * this down-type squark:
         */
      {
        singleSpecifiedDecayProduct
        = (*setIterator)->getInitialSparticle()->get_PDG_code();
        potentialDecayProducts
        = inputShortcut->getNeutralEwsbBosonsAndMassiveVectorBosons();
        for( std::vector< particlePointer >::const_iterator
             bosonIterator( potentialDecayProducts->begin() );
             potentialDecayProducts->end() > bosonIterator;
             ++bosonIterator )
        {
          subcascadeBranchingRatio
          = initialSparticle->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                                   singleSpecifiedDecayProduct,
                                            (*bosonIterator)->get_PDG_code() );
          if( lhcFaserGlobal::negligibleBr < subcascadeBranchingRatio )
            // if the BR to the next stage of the cascade is not negligible...
          {
            potentialSubcascades = (*setIterator)->getOpenCascades();
            for( std::vector< fullCascade* >::iterator
                 cascadeIterator( potentialSubcascades->begin() );
                 potentialSubcascades->end() > cascadeIterator;
                 ++cascadeIterator )
            {
              if( lhcFaserGlobal::negligibleBr
                  < ( subcascadeBranchingRatio
                      * (*cascadeIterator)->getBrToEwino() ) )
              {
                // we add each cascade with an overall BR that is not
                // negligible:
                openCascades.push_back(
                          squarkByBosonToCompound.addNewAtEnd()->setProperties(
                                                              initialSparticle,
                                       electroweakCascadeSource->getCascadeSet(
                                                              initialSparticle,
                                                                *bosonIterator,
                                       (*setIterator)->getInitialSparticle() ),
                                                          *cascadeIterator ) );
              }
            }
          }
        }
        // now we move on to the next potential sdown-type cascade:
        ++setIterator;
      }

      // now we look at compound cascades from up-type squarks:
      orderedCascadeSets = setOrderer->getSupTypeCascades();
      setIterator = orderedCascadeSets->begin();
      while( ( setIterator != orderedCascadeSets->end() )
             &&
             ( initialSparticle->get_absolute_mass()
               > (*setIterator)->getInitialSparticle()->get_absolute_mass() ) )
        /* we go through the list of up-type squarks until we get to either the
         * end of the list or run out of up-type squarks lighter than this
         * down-type squark:
         */
      {
        singleSpecifiedDecayProduct
        = (*setIterator)->getInitialSparticle()->get_PDG_code();
        potentialDecayProducts
        = inputShortcut->getChargedEwsbBosonsAndMassiveVectorBosons();
        for( std::vector< particlePointer >::const_iterator
             bosonIterator( potentialDecayProducts->begin() );
             potentialDecayProducts->end() > bosonIterator;
             ++bosonIterator )
        {
          subcascadeBranchingRatio
          = initialSparticle->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                                   singleSpecifiedDecayProduct,
                                         -((*bosonIterator)->get_PDG_code()) );
          // sdown-types decay to negatively-charged EW bosons.
          if( lhcFaserGlobal::negligibleBr < subcascadeBranchingRatio )
            // if the BR to the next stage of the cascade is not negligible...
          {
            potentialSubcascades = (*setIterator)->getOpenCascades();
            for( std::vector< fullCascade* >::iterator
                 cascadeIterator( potentialSubcascades->begin() );
                 potentialSubcascades->end() > cascadeIterator;
                 ++cascadeIterator )
            {
              if( lhcFaserGlobal::negligibleBr
                  < ( subcascadeBranchingRatio
                      * (*cascadeIterator)->getBrToEwino() ) )
              {
                // we add each cascade with an overall BR that is not
                // negligible:
                openCascades.push_back(
                          squarkByBosonToCompound.addNewAtEnd()->setProperties(
                                                              initialSparticle,
                                       electroweakCascadeSource->getCascadeSet(
                                                              initialSparticle,
                                                                *bosonIterator,
                                       (*setIterator)->getInitialSparticle() ),
                                                          *cascadeIterator ) );
              }
            }
          }
        }
        // now we move on to the next potential sup-type cascade:
        ++setIterator;
      }
    }



    supTypeSet::supTypeSet( inputHandler const* const inputShortcut,
           electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource,
                            particlePointer const initialScolored,
                            fullCascadeSetOrderer* const setOrderer,
                            double const beamEnergy,
                            fullCascadeSet* const gluinoFullCascade ) :
    fullCascadeSet( inputShortcut,
                    initialScolored,
                    electroweakCascadeSource,
                    beamEnergy ),
    setOrderer( setOrderer ),
    directToEwinoCascades(),
    compoundByBosonCascades(),
    compoundByJetCascades(),
    twoSpecifiedDecayProductsList( 2,
                                 CppSLHA::CppSLHA_global::really_wrong_value ),
    appropriateSdownForWDecay( NULL ),
    effectiveSdownMass( NULL ),
    effectiveSupMass( NULL )
    {
      if( inputShortcut->getSupL() == initialScolored )
      {
        appropriateSdownForWDecay = inputShortcut->getSdownL();
      }
      else if( inputShortcut->getSupR() == initialScolored )
      {
        appropriateSdownForWDecay = inputShortcut->getSdownR();
      }
      else if( inputShortcut->getScharmL() == initialScolored )
      {
        appropriateSdownForWDecay = inputShortcut->getSstrangeL();
      }
      else if( inputShortcut->getScharmR() == initialScolored )
      {
        appropriateSdownForWDecay = inputShortcut->getSstrangeR();
      }
      else if( inputShortcut->getStopTwo() == initialScolored )
      {
        appropriateSdownForWDecay = inputShortcut->getSbottomTwo();
      }
      else
      {
        appropriateSdownForWDecay = inputShortcut->getSbottomOne();
      }
      // we have to set up the cascades directly to electroweakinos now:
      for( std::vector< particlePointer >::const_iterator
           ewinoIterator( inputShortcut->getElectroweakinos()->begin() );
           inputShortcut->getElectroweakinos()->end() > ewinoIterator;
           ++ewinoIterator )
      {
        effectiveSdownMass
        = inputShortcut->getSquarkMinusBosonEffectiveMass( initialScolored,
                                                     inputShortcut->getWPlus(),
                                                           *ewinoIterator );
        effectiveSupMass
        = inputShortcut->getSquarkPlusBosonEffectiveMass(
                                                     appropriateSdownForWDecay,
                                                     inputShortcut->getWPlus(),
                                                          *ewinoIterator );
        directToEwinoCascades.push_back(
                                new fullCascadeType::supDirectlyToElectroweak(
                                                                 inputShortcut,
                                                               initialScolored,
                                                                    beamEnergy,
                                       electroweakCascadeSource->getCascadeSet(
                                                               initialScolored,
                                                              *ewinoIterator ),
                                       electroweakCascadeSource->getCascadeSet(
                                                              appropriateSdown,
                                                                *ewinoIterator,
                                                          effectiveSdownMass ),
                                       electroweakCascadeSource->getCascadeSet(
                                                               initialScolored,
                                                                *ewinoIterator,
                                                        effectiveSupMass ) ) );
      }
    }

    supTypeSet::~supTypeSet()
    {
      for( std::vector< fullCascadeType::supDirectlyToElectroweak* >::iterator
           deletionIterator( directToEwinoCascades.begin() );
           directToEwinoCascades.end() > deletionIterator;
           ++deletionIterator )
      {
        delete *deletionIterator;
      }
    }


    void
    supTypeSet::buildLongerCascades()
    {
      // 1st we clear the compound cascades:
      compoundByBosonCascades.clearEntries();
      compoundByJetCascades.clearEntries();

      // next we check to see if we should add compound cascades from the
      // gluino:
      if( inputShortcut->getGluinoMass() < initialSparticleMass )
      {
        singleSpecifiedDecayProductList.front() = CppSLHA::PDG_code::gluino;
        subcascadeBranchingRatio
        = initialSparticle->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                            &singleSpecifiedDecayProductList );
        if( lhcFaserGlobal::negligibleBr < subcascadeBranchingRatio )
          // if the BR to the next stage of the cascade is not negligible...
        {
          potentialSubcascades = gluinoFullCascade->getOpenCascades();
          for( std::vector< fullCascade* >::iterator
               cascadeIterator( potentialSubcascades->begin() );
               potentialSubcascades->end() > cascadeIterator;
               ++cascadeIterator )
          {
            if( lhcFaserGlobal::negligibleBr
                < ( subcascadeBranchingRatio
                    * (*cascadeIterator)->getBrToEwino() ) )
            {
              // we add each cascade with an overall BR that is not negligible:
              openCascades.push_back(
                            compoundByJetCascades.addNewAtEnd()->setProperties(
                                                              initialSparticle,
                                                          *cascadeIterator ) );
            }
          }
        }
      }

      // now we look at compound cascades from neutralinos:
      orderedCascadeSets = setOrderer->getNeutralinoColoredCascades();
      std::list< fullCascadeSet* >::iterator
      setIterator( orderedCascadeSets->begin() );
      while( ( setIterator != orderedCascadeSets->end() )
             &&
             ( initialSparticle->get_absolute_mass()
               > (*setIterator)->getInitialSparticle()->get_absolute_mass() ) )
        /* we go through the list of neutralinos heavier than the lightest
         * squark until we get to either the end of the list or a neutralino
         * heavier than this down-type squark:
         */
      {
        singleSpecifiedDecayProductList.front()
        = (*setIterator)->getInitialSparticle()->get_PDG_code();
        subcascadeBranchingRatio
        = initialSparticle->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                            &singleSpecifiedDecayProductList );
        if( lhcFaserGlobal::negligibleBr < subcascadeBranchingRatio )
          // if the BR to the next stage of the cascade is not negligible...
        {
          potentialSubcascades = (*setIterator)->getOpenCascades();
          for( std::vector< fullCascade* >::iterator
               cascadeIterator( potentialSubcascades->begin() );
               potentialSubcascades->end() > cascadeIterator;
               ++cascadeIterator )
          {
            if( lhcFaserGlobal::negligibleBr
                < ( subcascadeBranchingRatio
                    * (*cascadeIterator)->getBrToEwino() ) )
            {
              // we add each cascade with an overall BR that is not negligible:
              openCascades.push_back(
                            compoundByJetCascades.addNewAtEnd()->setProperties(
                                                              initialSparticle,
                                                          *cascadeIterator ) );
            }
          }
        }
        // now we move on to the next potential neutralino cascade:
        ++setIterator;
      }

      // now we look at compound cascades from charginos:
      orderedCascadeSets = setOrderer->getCharginoColoredCascades();
      setIterator = orderedCascadeSets->begin();
      while( ( setIterator != orderedCascadeSets->end() )
             &&
             ( initialSparticle->get_absolute_mass()
               > (*setIterator)->getInitialSparticle()->get_absolute_mass() ) )
        /* we go through the list of charginos heavier than the lightest squark
         * until we get to either the end of the list or a chargino heavier
         * than this down-type squark:
         */
      {
        singleSpecifiedDecayProductList.front()
        = -((*setIterator)->getInitialSparticle()->get_PDG_code());
        // down-type squarks decay to negatively-charged charginos.
        subcascadeBranchingRatio
        = initialSparticle->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                            &singleSpecifiedDecayProductList );
        if( lhcFaserGlobal::negligibleBr < subcascadeBranchingRatio )
          // if the BR to the next stage of the cascade is not negligible...
        {
          potentialSubcascades = (*setIterator)->getOpenCascades();
          for( std::vector< fullCascade* >::iterator
               cascadeIterator( potentialSubcascades->begin() );
               potentialSubcascades->end() > cascadeIterator;
               ++cascadeIterator )
          {
            if( lhcFaserGlobal::negligibleBr
                < ( subcascadeBranchingRatio
                    * (*cascadeIterator)->getBrToEwino() ) )
            {
              // we add each cascade with an overall BR that is not negligible:
              openCascades.push_back(
                            compoundByJetCascades.addNewAtEnd()->setProperties(
                                                              initialSparticle,
                                                          *cascadeIterator ) );
            }
          }
        }
        // now we move on to the next potential chargino cascade:
        ++setIterator;
      }

      // now we look at compound cascades from down-type squarks:
      orderedCascadeSets = setOrderer->getSdownTypeCascades();
      setIterator = orderedCascadeSets->begin();
      while( ( setIterator != orderedCascadeSets->end() )
             &&
             ( initialSparticle->get_absolute_mass()
               > (*setIterator)->getInitialSparticle()->get_absolute_mass() ) )
        /* we go through the list of down-type squarks until we get to either
         * the end of the list or run out of down-type squarks lighter than
         * this down-type squark:
         */
      {
        singleSpecifiedDecayProduct
        = (*setIterator)->getInitialSparticle()->get_PDG_code();
        potentialDecayProducts
        = inputShortcut->getNeutralEwsbBosonsAndMassiveVectorBosons();
        for( std::vector< particlePointer >::const_iterator
             bosonIterator( potentialDecayProducts->begin() );
             potentialDecayProducts->end() > bosonIterator;
             ++bosonIterator )
        {
          subcascadeBranchingRatio
          = initialSparticle->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                                   singleSpecifiedDecayProduct,
                                            (*bosonIterator)->get_PDG_code() );
          if( lhcFaserGlobal::negligibleBr < subcascadeBranchingRatio )
            // if the BR to the next stage of the cascade is not negligible...
          {
            potentialSubcascades = (*setIterator)->getOpenCascades();
            for( std::vector< fullCascade* >::iterator
                 cascadeIterator( potentialSubcascades->begin() );
                 potentialSubcascades->end() > cascadeIterator;
                 ++cascadeIterator )
            {
              if( lhcFaserGlobal::negligibleBr
                  < ( subcascadeBranchingRatio
                      * (*cascadeIterator)->getBrToEwino() ) )
              {
                // we add each cascade with an overall BR that is not
                // negligible:
                openCascades.push_back(
                          squarkByBosonToCompound.addNewAtEnd()->setProperties(
                                                              initialSparticle,
                                       electroweakCascadeSource->getCascadeSet(
                                                              initialSparticle,
                                                                *bosonIterator,
                                       (*setIterator)->getInitialSparticle() ),
                                                          *cascadeIterator ) );
              }
            }
          }
        }
        // now we move on to the next potential sdown-type cascade:
        ++setIterator;
      }

      // now we look at compound cascades from up-type squarks:
      orderedCascadeSets = setOrderer->getSupTypeCascades();
      setIterator = orderedCascadeSets->begin();
      while( ( setIterator != orderedCascadeSets->end() )
             &&
             ( initialSparticle->get_absolute_mass()
               > (*setIterator)->getInitialSparticle()->get_absolute_mass() ) )
        /* we go through the list of up-type squarks until we get to either the
         * end of the list or run out of up-type squarks lighter than this
         * down-type squark:
         */
      {
        singleSpecifiedDecayProduct
        = (*setIterator)->getInitialSparticle()->get_PDG_code();
        potentialDecayProducts
        = inputShortcut->getChargedEwsbBosonsAndMassiveVectorBosons();
        for( std::vector< particlePointer >::const_iterator
             bosonIterator( potentialDecayProducts->begin() );
             potentialDecayProducts->end() > bosonIterator;
             ++bosonIterator )
        {
          subcascadeBranchingRatio
          = initialSparticle->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                                   singleSpecifiedDecayProduct,
                                         -((*bosonIterator)->get_PDG_code()) );
          // sdown-types decay to negatively-charged EW bosons.
          if( lhcFaserGlobal::negligibleBr < subcascadeBranchingRatio )
            // if the BR to the next stage of the cascade is not negligible...
          {
            potentialSubcascades = (*setIterator)->getOpenCascades();
            for( std::vector< fullCascade* >::iterator
                 cascadeIterator( potentialSubcascades->begin() );
                 potentialSubcascades->end() > cascadeIterator;
                 ++cascadeIterator )
            {
              if( lhcFaserGlobal::negligibleBr
                  < ( subcascadeBranchingRatio
                      * (*cascadeIterator)->getBrToEwino() ) )
              {
                // we add each cascade with an overall BR that is not
                // negligible:
                openCascades.push_back(
                          squarkByBosonToCompound.addNewAtEnd()->setProperties(
                                                              initialSparticle,
                                       electroweakCascadeSource->getCascadeSet(
                                                              initialSparticle,
                                                                *bosonIterator,
                                       (*setIterator)->getInitialSparticle() ),
                                                          *cascadeIterator ) );
              }
            }
          }
        }
        // now we move on to the next potential sup-type cascade:
        ++setIterator;
      }
    }


  gluinoSet::gluinoSet(
                                       inputHandler const* const inputShortcut,
           electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource,
                           readiableCascadeSetList* const squarkCascadeSetList,
                                              double const beamEnergy ) :
    fullCascadeSet( inputShortcut,
                    inputShortcut->getGluino(),
                    electroweakCascadeSource,
                    squarkCascadeSetList,
                    beamEnergy )
  {
    // we have to set up gxCascades now:
    for( std::vector< particlePointer >::const_iterator
         ewinoIterator( inputShortcut->getElectroweakinos()->begin() );
         inputShortcut->getElectroweakinos()->end() > ewinoIterator;
         ++ewinoIterator )
    {
      gxCascades.addNewAtEnd()->setProperties( inputShortcut,
                                               beamEnergy,
                                       electroweakCascadeSource->getCascadeSet(
                                                               initialSparticle,
                                                            *ewinoIterator ) );
    }
  }

  gluinoSet::~gluinoSet()
  {
    // does nothing.
  }


  void
  gluinoSet::buildLongerCascades()
  {
    gjsxCascades.clearEntries();
    gbsxCascades.clearEntries();
    gjsbsxCascades.clearEntries();
    // now we look for decays into squarks (squarkCascadeSetList has been
    // ordered for this point already in fullCascadeSet::setUpCascades()):
    for( std::list< fullCascadeSet* >::iterator
         cascadeIterator( squarkCascadeSetList->first->begin() );
         ( ( squarkCascadeSetList->first->end() != cascadeIterator )
           &&
           ( (*cascadeIterator)->getInitialSparticle()->get_absolute_mass()
             < initialSparticle->get_absolute_mass() ) );
         ++cascadeIterator )
    {
      // each of the squarks looked at in this loop are lighter than the
      // gluino, so these cascades should be built:
      for( std::vector< sxFullCascade* >::const_iterator
           sxIterator( (*cascadeIterator)->getSxCascades()->begin() );
           (*cascadeIterator)->getSxCascades()->end() > sxIterator;
           ++sxIterator )
      {
        gjsxCascades.addNewAtEnd()->setProperties( *sxIterator );
        addOpenCascade( gjsxCascades.getBack() );
      }
      for( std::vector< sbsxFullCascade* >::const_iterator
           svsxIterator( (*cascadeIterator)->getSbsxCascades()->begin() );
           (*cascadeIterator)->getSbsxCascades()->end() > svsxIterator;
           ++svsxIterator )
      {
        gjsbsxCascades.addNewAtEnd()->setProperties( *svsxIterator );
        addOpenCascade( gjsbsxCascades.getBack() );
      }

      // now we check for gvs:
      for( std::vector< particlePointer >::const_iterator
           ewIterator( inputShortcut->getEwsbBosonsAndMassiveVectorBosons(
                                                                  )->begin() );
           inputShortcut->getEwsbBosonsAndMassiveVectorBosons()->end()
           > ewIterator;
           ++ewIterator )
      {
        soughtPositivePdgCodeList.clear();
        soughtPositivePdgCodeList.push_back(
                    (*cascadeIterator)->getInitialSparticle()->get_PDG_code() );
        soughtPositivePdgCodeList.push_back( (*ewIterator)->get_PDG_code() );
        soughtNegativePdgCodeList.clear();
        soughtNegativePdgCodeList.push_back(
                    (*cascadeIterator)->getInitialSparticle()->get_PDG_code() );
        soughtNegativePdgCodeList.push_back(
                                            -((*ewIterator)->get_PDG_code()) );
        if( ( initialSparticle->get_absolute_mass()
              > ( (*cascadeIterator)->getInitialSparticle()->get_absolute_mass()
                  + (*ewIterator)->get_absolute_mass() ) )
            &&
            ( ( lhcFaserGlobal::negligibleBr
                < initialSparticle->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                                 &soughtPositivePdgCodeList ) )
              || ( lhcFaserGlobal::negligibleBr
                               < initialSparticle->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                             &soughtNegativePdgCodeList ) ) ) )
        {
          for( std::vector< sxFullCascade* >::const_iterator
               sxIterator( (*cascadeIterator)->getSxCascades()->begin() );
               (*cascadeIterator)->getSxCascades()->end() > sxIterator;
               ++sxIterator )
          {
            gbsxCascades.addNewAtEnd()->setProperties(
                                       electroweakCascadeSource->getCascadeSet(
                                                               initialSparticle,
                                                                 *ewIterator,
                                    (*cascadeIterator)->getInitialSparticle() ),
                                                       *sxIterator );
            addOpenCascade( gbsxCascades.getBack() );
          }
        }
      }  // end of loop over if gvs is open
    }  // end of loop over squarkFullCascadeSets
  }



  }  // end of fullCascadeSetType namespace

  fullCascadeSetsForOneBeamEnergy::fullCascadeSetsForOneBeamEnergy(
                                       inputHandler const* const inputShortcut,
           electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource,
                                                    double const beamEnergy ) :
    inputShortcut( inputShortcut ),
    electroweakCascadeSource( electroweakCascadeSource ),
    beamEnergy( beamEnergy ),
    squarkCascadeSetList( new std::list< fullCascadeSet* >,
              new publicGetsReadiedForNewPoint( inputShortcut->getReadier() ) )
  {
    gluinoCascadeSet = new gluinoSet( inputShortcut,
                                                 electroweakCascadeSource,
                                                 &squarkCascadeSetList,
                                                 beamEnergy );
    for( std::vector< particlePointer >::const_iterator
         squarkIterator( inputShortcut->getSquarks()->begin() );
         inputShortcut->getSquarks()->end() > squarkIterator;
         ++squarkIterator )
    {
      squarkCascadeSets.push_back( new sdownTypeSet( inputShortcut,
                                                      electroweakCascadeSource,
                                                             *squarkIterator,
                                                         &squarkCascadeSetList,
                                                             gluinoCascadeSet,
                                                             beamEnergy ) );
      squarkCascadeSetList.first->push_back( squarkCascadeSets.back() );
    }
  }

  fullCascadeSetsForOneBeamEnergy::~fullCascadeSetsForOneBeamEnergy()
  {
    for( std::vector< sdownTypeSet* >::iterator
         deletionIterator( squarkCascadeSets.begin() );
         squarkCascadeSets.end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
    delete gluinoCascadeSet;
    delete squarkCascadeSetList.first;
    delete squarkCascadeSetList.second;
  }


  fullCascadeSet*
  fullCascadeSetsForOneBeamEnergy::getFullCascadeSet(
                                        particlePointer const initialScolored )
  /* this returns the fullCascadeSet for the requested colored sparticle, or
   * NULL if we were asked for a sparticle that is not the gluino or in
   * inputShortcut->getSquarks().
   */
  {
    if( CppSLHA::PDG_code::gluino == initialScolored->get_PDG_code() )
    {
      return gluinoCascadeSet;
    }
    else
    {
      fullCascadeSet* returnPointer( NULL );
      // we look to see if we already have a fullCascadeSet for this squark:
      for( std::vector< sdownTypeSet* >::iterator
           searchIterator( squarkCascadeSets.begin() );
           squarkCascadeSets.end() > searchIterator;
           ++searchIterator )
      {
        if( initialScolored == (*searchIterator)->getInitialSparticle() )
        {
          returnPointer = *searchIterator;
          searchIterator = squarkCascadeSets.end();
        }
      }
      return returnPointer;
    }
  }



  fullCascadeSetFactory::fullCascadeSetFactory(
                                       inputHandler const* const inputShortcut,
                  electroweakCascadeHandler* const electroweakCascadeSource ) :
      inputShortcut( inputShortcut ),
      electroweakCascadeSource( electroweakCascadeSource )
  {
    // just an initialization list.
  }


  fullCascadeSetFactory::~fullCascadeSetFactory()
  {
    for( std::vector< fullCascadeSetsForOneBeamEnergy* >::iterator
         deletionIterator( cascadeSetsPerEnergy.begin() );
         cascadeSetsPerEnergy.end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
  }



  fullCascadeSetsForOneBeamEnergy*
  fullCascadeSetFactory::getFullCascadeSetsForOneBeamEnergy(
                                                         int const beamEnergy )
  /* this looks to see if there is an existing
   * fullCascadeSetsForOneBeamEnergy with the requested values, & if not,
   * makes 1, & returns the pointer.
   */
  {
    fullCascadeSetsForOneBeamEnergy* returnPointer( NULL );
    // we look to see if we already have a fullCascadeSetsForOneBeamEnergy for
    // these acceptanceValues:
    for( std::vector< fullCascadeSetsForOneBeamEnergy* >::iterator
         searchIterator( cascadeSetsPerEnergy.begin() );
         cascadeSetsPerEnergy.end() > searchIterator;
         ++searchIterator )
    {
      if( beamEnergy == (*searchIterator)->getBeamEnergy() )
      {
        returnPointer = *searchIterator;
        searchIterator = cascadeSetsPerEnergy.end();
      }
    }
    if( NULL == returnPointer )
      // if we do not already have a fullCascadeSetsForOneBeamEnergy for this
      // beam energy, we make a new instance:
    {
      returnPointer
      = new fullCascadeSetsForOneBeamEnergy( inputShortcut,
              electroweakCascadeSource->getElectroweakCascadesForOneBeamEnergy(
                                                                  beamEnergy ),
                                             beamEnergy );
      cascadeSetsPerEnergy.push_back( returnPointer );
    }
    return returnPointer;
  }


}  // end of LHC_FASER namespace.
