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
    squarkDirectlyToElectroweak::squarkDirectlyToElectroweak(
                                       inputHandler const* const inputShortcut,
                                           particlePointer const initialSquark,
                                                       double const beamEnergy,
                           electroweakCascadeSet* const directEwinoCascades ) :
        fullCascade( squarkToEwino,
                     2,
                     inputShortcut,
                     initialSquark,
                     beamEnergy ),
        directEwinoCascades( directEwinoCascades )
    {
      cascadeSegment = cascadeDefiner.addNewAtEnd();
      cascadeSegment->first = directEwinoCascades->getElectroweakDecayer();
      // squarkToEwino means only 1 decay to be recorded.
      cascadeSegment->second = firstDecayBodyNumber;
      // squarkToEwino also means that the decay is 2-body.
    }

    squarkDirectlyToElectroweak::~squarkDirectlyToElectroweak()
    {
      // does nothing.
    }



    namespace squarkDirectlyToElectroweakType
    {
      sdownType::sdownType( inputHandler const* const inputShortcut,
                            particlePointer const initialSquark,
                            double const beamEnergy,
                            electroweakCascadeSet* const directEwinoCascades ) :
          squarkDirectlyToElectroweak( inputShortcut,
                                       initialSquark,
                                       beamEnergy,
                                       directEwinoCascades ),
          ewinoFlipsCharge( false )
      {
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
          soughtDecayProductList.front()
          = cascadeSegment->first->get_PDG_code();
        }
      }

      sdownType::~sdownType()
      {
        // does nothing.
      }

      double
      sdownType::getAcceptance( bool const initialSparticleIsNotAntiparticle,
                                acceptanceCutSet* const acceptanceCuts,
                                int const numberOfAdditionalJets,
                                int numberOfNegativeElectrons,
                                int numberOfPositiveElectrons,
                                int numberOfNegativeMuons,
                                int numberOfPositiveMuons )
      /* this calls the appropriate functions on ewinoCascades to build the
       * required acceptance, taking into account whether the charges should be
       * swapped if scoloredIsNotAntiparticle is false.
       */
      {
        if( ( initialSparticleIsNotAntiparticle
              &&
              ewinoFlipsCharge )
            ||
            ( !initialSparticleIsNotAntiparticle
              &&
              !ewinoFlipsCharge ) )
          // if the electroweakino is an antiparticle, we swap the charges.
        {
          int swappingInt( numberOfPositiveElectrons );
          numberOfPositiveElectrons = numberOfNegativeElectrons;
          numberOfNegativeElectrons = swappingInt;
          swappingInt = numberOfPositiveMuons;
          numberOfPositiveMuons = numberOfNegativeMuons;
          numberOfNegativeMuons = swappingInt;
        }
        return ewinoCascades->getAcceptance( acceptanceCuts,
                                             numberOfAdditionalJets,
                                             numberOfNegativeElectrons,
                                             numberOfPositiveElectrons,
                                             numberOfNegativeMuons,
                                             numberOfPositiveMuons );
      }



      supType::supType( inputHandler const* const inputShortcut,
                        particlePointer const initialSquark,
                        double const beamEnergy,
                        electroweakCascadeSet* const directEwinoCascades,
                        electroweakCascadeSet* const ewinoWithWCascades,
                        electroweakCascadeSet* const bosonCascades ) :
          squarkDirectlyToElectroweak( inputShortcut,
                                       initialSquark,
                                       beamEnergy,
                                       directEwinoCascades ),
          ewinoWithWCascades( ewinoWithWCascades ),
          bosonCascades( bosonCascades ),
          wBoson( inputShortcut->getWPlus() ),
          ewinoMass( CppSLHA::CppSLHA_global::really_wrong_value ),
          shouldUseDecaysWithW( false ),
          directFraction( CppSLHA::CppSLHA_global::really_wrong_value ),
          wFraction( CppSLHA::CppSLHA_global::really_wrong_value ),
          decayProductListIncludingW()
      {
        soughtDecayProductList.front() = cascadeSegment->first->get_PDG_code();
        decayProductListIncludingW.push_back( soughtDecayProductList.front() );
        decayProductListIncludingW.push_back(
                 ewinoWithWCascades->getElectroweakDecayer()->get_PDG_code() );
      }

      supType::~supType()
      {
        // does nothing.
      }


      bool
      supType::isOpen()
      /* this returns true if the squark is heavy enough to decay into the
       * electroweakino, false otherwise. it also sorts out whether it should
       * be using the decays involving a W boson.
       */
      {
        ewinoMass
        = directEwinoCascades->getElectroweakDecayer()->get_absolute_mass();
        if( ( initialSparticle->get_absolute_mass()
              > ewinoCascades->getElectroweakDecayer()->get_absolute_mass() )
            &&
            ( lhcFaserGlobal::negligibleBr < getTotalBrToEwino() ) )
        {
          if( initialSparticle->get_absolute_mass()
              > ( ewinoMass + wBoson->get_absolute_mass() ) )
          {
            shouldUseDecaysWithW = decayWithWIsNotNegligible();
          }
          else
          {
            shouldUseDecaysWithW = false;
          }
          return true;
        }
        else
        {
          return false;
        }
      }

      double
      supType::getAcceptance( bool const initialSparticleIsNotAntiparticle,
                              acceptanceCutSet* const acceptanceCuts,
                              int const numberOfAdditionalJets,
                              int const numberOfNegativeElectrons,
                              int const numberOfPositiveElectrons,
                              int const numberOfNegativeMuons,
                              int const numberOfPositiveMuons )
      /* this calls the appropriate functions on directEwinoCascades to build
       * the required acceptance or ewinoWithWCascades & bosonCascades if the
       * decay with a boson is not negligible, taking into account whether the
       * charges should be swapped if scoloredIsNotAntiparticle is false.
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
      supType::getCombinedAcceptance( acceptanceCutSet* const acceptanceCuts,
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

    }  // end of squarkDirectlyToElectroweakType namespace

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



    squarkByJetToCompound::squarkByJetToCompound() :
        fullCascade( sjm,
                     2 )
    // sjm also means that the initial decay is 2-body.
    {
      // just an initialization list.
    }

    squarkByJetToCompound::~squarkByJetToCompound()
    {
      // does nothing.
    }



    namespace squarkByJetToCompoundType
    {
      sdownByJetToCompound::sdownByJetToCompound() :
          squarkByJetToCompound()
      {
        // just an initialization list.
      }

      sdownByJetToCompound::~sdownByJetToCompound()
      {
        // does nothing.
      }



      supByJetToCompound::supByJetToCompound() :
          squarkByJetToCompound(),
          bosonCascades( NULL ),
          wBoson( NULL ),
          ewinoMass( CppSLHA::CppSLHA_global::really_wrong_value ),
          shouldUseDecaysWithW( false ),
          directFraction( CppSLHA::CppSLHA_global::really_wrong_value ),
          wFraction( CppSLHA::CppSLHA_global::really_wrong_value ),
          decayProductListIncludingW( 2,
                                  CppSLHA::CppSLHA_global::really_wrong_value )
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
       * required acceptance, & also with bosonCascades if the decay with a
       * boson is not negligible, taking into account whether the charges
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

      bool
      supByJetToCompound::decayWithWIsNotNegligible()
      // this returns true if the decay involving a W boson is not negligible,
      // also setting up the relevant fractions, false otherwise.
      {
        wFraction = ( initialSparticle->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                                  &decayProductListIncludingW )
                      / initialSparticle->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                                   &soughtDecayProductList ) );
        if( lhcFaserGlobal::negligibleBr < wFraction )
        {
          directFraction = ( 1.0 - wFraction );
          bosonCascades
          = electroweakCascadeSource->getCascadeSet( initialSparticle,
                                                     wBoson,
                                      subcascadePointer->getInitialSparticle(),
                               inputShortcut->getSquarkMinusBosonEffectiveMass(
                                                              initialSparticle,
                                                                        wBoson,
                                  subcascadePointer->getInitialSparticle() ) );
          return true;
        }
        else
        {
          return false;
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

    }  // end of squarkByJetToCompoundType namespace



    gluinoOrElectroweakinoToCompound::gluinoOrElectroweakinoToCompound() :
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

    gluinoOrElectroweakinoToCompound::~gluinoOrElectroweakinoToCompound()
    {
      // does nothing.
    }



    namespace gluinoOrElectroweakinoToCompoundType
    {
      gluinoOrNeutralinoSet::gluinoOrNeutralinoSet() :
          gluinoOrElectroweakinoToCompound()
      {
        // just an initialization list.
      }

      gluinoOrNeutralinoSet::~gluinoOrNeutralinoSet()
      {
        // does nothing.
      }


      double
      gluinoOrNeutralinoSet::getAcceptance(
                                  bool const initialSparticleIsNotAntiparticle,
                                        acceptanceCutSet* const acceptanceCuts,
                                            int const numberOfAdditionalJets,
                                            int numberOfNegativeElectrons,
                                            int numberOfPositiveElectrons,
                                            int numberOfNegativeMuons,
                                            int numberOfPositiveMuons )
      /* this calls the appropriate functions on subcascadePointer to build the
       * required acceptance, & also with bosonCascades if the decay with a
       * boson is not negligible, taking into account whether the charges
       * should be swapped if scoloredIsNotAntiparticle is false.
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

      bool
      gluinoOrNeutralinoSet::decayWithWIsNotNegligible()
      // this returns true if the decay involving a W boson is not negligible,
      // also setting up the relevant fractions, false otherwise.
      {
        wFraction = ( initialSparticle->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                                  &decayProductListIncludingW )
                      / initialSparticle->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                                   &soughtDecayProductList ) );
        if( lhcFaserGlobal::negligibleBr < wFraction )
        {
          directFraction = ( 1.0 - wFraction );
          bosonCascades
          = electroweakCascadeSource->getCascadeSet( initialSparticle,
                                                     wBoson,
                                      subcascadePointer->getInitialSparticle(),
                               inputShortcut->getSquarkMinusBosonEffectiveMass(
                                                              initialSparticle,
                                                                        wBoson,
                                  subcascadePointer->getInitialSparticle() ) );
          return true;
        }
        else
        {
          return false;
        }
      }

      double
      gluinoOrNeutralinoSet::getCombinedAcceptance(
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
                                 * bosonCascades->getAcceptance(
                                                                acceptanceCuts,
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



      charginoSet::charginoSet() :
          gluinoOrElectroweakinoToCompound()
      {
        // just an initialization list.
      }

      charginoSet::~charginoSet()
      {
        // does nothing.
      }


      fullCascade*
      charginoSet::setProperties( particlePointer const initialSparticle,
                                  fullCascade* const subcascadePointer,
          electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource )
      {
        this->initialSparticle = initialSparticle;
        this->electroweakCascadeSource = electroweakCascadeSource;
        buildOn( subcascadePointer );
        resetCachedBranchingRatio();
        decayProductListIncludingW.front()
        = subcascadePointer->getInitialSparticle()->get_PDG_code();
        decayProductListIncludingW.back()
        = bosonCascades->getElectroweakDecayer()->get_PDG_code();
        // positive charginos decay to positive W bosons.
        if( inputShortcut->isIn( decayProductListIncludingW.front(),
                                 inputShortcut->getSdownTypes() ) )
        {
          decayingToSupType = false;
        }
        else
        {
          decayProductListIncludingW.front()
          = subcascadePointer->getInitialSparticle()->get_PDG_code();
          decayingToSupType = true;
        }
        shouldUseDecaysWithW = decayWithWIsNotNegligible();
        return this;
      }

      double
      charginoSet::getAcceptance( bool const initialSparticleIsNotAntiparticle,
                                  acceptanceCutSet* const acceptanceCuts,
                                  int const numberOfAdditionalJets,
                                  int numberOfNegativeElectrons,
                                  int numberOfPositiveElectrons,
                                  int numberOfNegativeMuons,
                                  int numberOfPositiveMuons )
      /* this calls the appropriate functions on subcascadePointer to build the
       * required acceptance, & also with bosonCascades if the decay with a
       * boson is not negligible, taking into account whether the charges
       * should be swapped if scoloredIsNotAntiparticle is false.
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

      bool
      charginoSet::decayWithWIsNotNegligible()
      // this returns true if the decay involving a W boson is not negligible,
      // also setting up the relevant fractions, false otherwise.
      {
        // 1st we look for chargino -> squark + antiquark + W:
        wFraction = ( initialSparticle->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                                  &decayProductListIncludingW )
                      / initialSparticle->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                                   &soughtDecayProductList ) );
        squarkWithWFraction = wFraction;
        // now we look for chargino -> antisquark + quark + W:
        decayProductListIncludingW.front()
        = -(decayProductListIncludingW.front());
        wFraction += ( initialSparticle->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                                  &decayProductListIncludingW )
                       / initialSparticle->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                                   &soughtDecayProductList ) );
        if( lhcFaserGlobal::negligibleBr < wFraction )
        {
          directFraction = ( 1.0 - wFraction );
          antisquarkWithWFraction = ( wFraction - squarkWithWFraction );
          bosonCascades
          = electroweakCascadeSource->getCascadeSet( initialSparticle,
                                                     wBoson,
                                      subcascadePointer->getInitialSparticle(),
                               inputShortcut->getSquarkMinusBosonEffectiveMass(
                                                              initialSparticle,
                                                                        wBoson,
                                  subcascadePointer->getInitialSparticle() ) );
          return true;
        }
        else
        {
          return false;
        }
      }

      double
      charginoSet::getCombinedAcceptance(
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

    }  // end of gluinoOrElectroweakinoType namespace

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
    setIterator(),
    potentialSubcascades( NULL ),
    subcascadeBranchingRatio( CppSLHA::CppSLHA_global::really_wrong_value ),
    beamEnergy( beamEnergy ),
    singleSpecifiedDecayProductList( 1,
                                 CppSLHA::CppSLHA_global::really_wrong_value ),
    singleSpecifiedDecayProduct( CppSLHA::CppSLHA_global::really_wrong_value ),
    potentialDecayProducts( NULL ),
    potentialSubcascadeSet( NULL )
  {
    // just an initialization list.
  }

  fullCascadeSet::~fullCascadeSet()
  {
    // does nothing.
  }



  fullCascadeSetOrderer::fullCascadeSetOrderer(
                                    inputHandler const* const inputShortcut ) :
      getsReadiedForNewPoint( inputShortcut->getReadier() ),
      inputShortcut( inputShortcut ),
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
    squarkSet::squarkSet( inputHandler const* const inputShortcut,
           electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource,
                          particlePointer const initialScolored,
                          fullCascadeSetOrderer* const setOrderer,
                          fullCascadeSet* const gluinoFullCascade,
                          double const beamEnergy,
                          bool const ewinoCodeIsAlwaysPositive,
                          bool const bosonCodeIsAlwaysPositive ) :
        fullCascadeSet( inputShortcut,
                        initialScolored,
                        electroweakCascadeSource,
                        beamEnergy ),
        setOrderer( setOrderer ),
        gluinoFullCascade( gluinoFullCascade ),
        directToEwinoCascades(),
        compoundByBosonCascades(),
        compoundByJetCascades(),
        ewinoCodeIsAlwaysPositive( ewinoCodeIsAlwaysPositive ),
        bosonCodeIsAlwaysPositive( bosonCodeIsAlwaysPositive )
    {
      // just an initialization list.
    }

    squarkSet::~squarkSet()
    {
      for( std::vector< fullCascadeType::squarkDirectlyToElectroweak*
                                                                    >::iterator
           deletionIterator( directToEwinoCascades.begin() );
           directToEwinoCascades.end() > deletionIterator;
           ++deletionIterator )
      {
        delete *deletionIterator;
      }
    }


    void
    squarkSet::setUpCascades()
    {
      // 1st we add the direct cascades:
      for( std::vector< fullCascadeType::squarkDirectlyToElectroweak*
                                                                    >::iterator
           cascadeIterator( directToEwinoCascades.begin() );
           directToEwinoCascades.end() > cascadeIterator;
           ++cascadeIterator )
      {
        if( (*cascadeIterator)->isOpen() )
        {
          openCascades.push_back( *cascadeIterator );
        }
      }

      // now we clear the compound cascades:
      compoundByBosonCascades.clearEntries();
      compoundByJetCascades.clearEntries();

      // next we check to see if we should add compound cascades from the
      // gluino:
      if( inputShortcut->getGluinoMass() < initialSparticleMass )
      {
        singleSpecifiedDecayProductList.front() = CppSLHA::PDG_code::gluino;
        potentialSubcascadeSet = gluinoFullCascade;
        addByJetCascade();
      }

      // now we look at compound cascades from neutralinos:
      orderedCascadeSets = setOrderer->getNeutralinoColoredCascades();
      addEwinoCompoundCascade( true );

      // now we look at compound cascades from charginos:
      orderedCascadeSets = setOrderer->getCharginoColoredCascades();
      addEwinoCompoundCascade( ewinoCodeIsAlwaysPositive );
      // down-type squarks decay to negatively-charged charginos, up-type to
      // positively-charged.

      // now we look at compound cascades from squarks:
      buildSquarkCompoundCascades();
      // now we look at compound cascades from down-type squarks:
      orderedCascadeSets = setOrderer->getSdownTypeCascades();
      potentialDecayProducts
      = inputShortcut->getNeutralEwsbBosonsAndMassiveVectorBosons();
      addSquarkCompoundCascade( true );

      // now we look at compound cascades from up-type squarks:
      orderedCascadeSets = setOrderer->getSupTypeCascades();
      potentialDecayProducts
      = inputShortcut->getChargedEwsbBosonsAndMassiveVectorBosons();
      addSquarkCompoundCascade( false );
    }

    void
    squarkSet::addByJetCascade()
    // this relies on singleSpecifiedDecayProductList.front() &
    // potentialSubcascadeSet having already been set correctly.
    {
      subcascadeBranchingRatio
      = initialSparticle->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                            &singleSpecifiedDecayProductList );
      if( lhcFaserGlobal::negligibleBr < subcascadeBranchingRatio )
        // if the BR to the next stage of the cascade is not negligible...
      {
        potentialSubcascades = potentialSubcascadeSet->getOpenCascades();
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
                                                              *cascadeIterator,
                                                  electroweakCascadeSource ) );
          }
        }
      }
    }

    void
    squarkSet::addSquarkCompoundCascade( bool positiveBosonCode )
    // this relies on orderedCascadeSets & potentialDecayProducts having
    // already been set correctly.
    {
      setIterator = orderedCascadeSets->begin();
      int bosonCode;
      while( ( setIterator != orderedCascadeSets->end() )
             &&
             ( initialSparticleMass
               > (*setIterator)->getInitialSparticle()->get_absolute_mass() ) )
        // we go through the list of squarks until we get to either the end of
        // the list or run out of squarks lighter than this down-type squark:
      {
        singleSpecifiedDecayProduct
        = (*setIterator)->getInitialSparticle()->get_PDG_code();
        for( std::vector< particlePointer >::const_iterator
             bosonIterator( potentialDecayProducts->begin() );
             potentialDecayProducts->end() > bosonIterator;
             ++bosonIterator )
        {
          if( positiveBosonCode )
          {
            bosonCode = (*bosonIterator)->get_PDG_code();
          }
          else
          {
            bosonCode = -((*bosonIterator)->get_PDG_code());
          }
          subcascadeBranchingRatio
          = initialSparticle->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                                   singleSpecifiedDecayProduct,
                                                            bosonCode );
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
        // now we move on to the next potential squark cascade:
        ++setIterator;
      }
    }



    namespace squarkSetType
    {
      sdownType::sdownType( inputHandler const* const inputShortcut,
           electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource,
                            particlePointer const initialScolored,
                            fullCascadeSetOrderer* const setOrderer,
                            fullCascadeSet* const gluinoFullCascade,
                            double const beamEnergy ) :
          squarkSet( inputShortcut,
                     electroweakCascadeSource,
                     initialScolored,
                     setOrderer,
                     gluinoFullCascade,
                     beamEnergy,
                     false,
                     false )
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

      sdownType::~sdownType()
      {
        // does nothing.
      }



      supType::supType( inputHandler const* const inputShortcut,
           electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource,
                        particlePointer const initialScolored,
                        fullCascadeSetOrderer* const setOrderer,
                        fullCascadeSet* const gluinoFullCascade,
                        double const beamEnergy ) :
          squarkSet( inputShortcut,
                     electroweakCascadeSource,
                     initialScolored,
                     setOrderer,
                     gluinoFullCascade,
                     beamEnergy ),
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

      supType::~supType()
      {
        // does nothing.
      }

      void
      supType::findOpenDirectCascades()
      // this puts all open direct cascades into openCascades.
      {
        for( std::vector< fullCascadeType::sdownDirectlyToElectroweak*
                                                                    >::iterator
             cascadeIterator( directToEwinoCascades.begin() );
             directToEwinoCascades.end() > cascadeIterator;
             ++cascadeIterator )
        {
          if( (*cascadeIterator)->isOpen() )
          {
            openCascades.push_back( *cascadeIterator );
          }
        }
      }

    }  // end of squarkSetType namespace



    gluinoOrElectroweakinoSet::gluinoOrElectroweakinoSet(
                                       inputHandler const* const inputShortcut,
           electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource,
                                        particlePointer const initialSparticle,
                                       fullCascadeSetOrderer* const setOrderer,
                                                  double const beamEnergy ) :
    fullCascadeSet( inputShortcut,
                    initialSparticle,
                    electroweakCascadeSource,
                    beamEnergy ),
    setOrderer( setOrderer )
    {
      // just an initialization list.
    }

    gluinoOrElectroweakinoSet::~gluinoOrElectroweakinoSet()
    {
      // does nothing.
    }


    void
    gluinoOrElectroweakinoSet::buildLongerCascades()
    {
      // 1st we clear the compound cascades:
      compoundCascades.clearEntries();

      // now we look at compound cascades from squarks:
      orderedCascadeSets = setOrderer->getSdownTypeCascades();
      buildSquarkCompoundCascades();
      orderedCascadeSets = setOrderer->getSupTypeCascades();
      buildSquarkCompoundCascades();
    }

    void
    gluinoOrElectroweakinoSet::buildSquarkCompoundCascades()
    // this does the job of finding the right squark subcascades.
    {
      setIterator = orderedCascadeSets->begin();
      while( ( setIterator != orderedCascadeSets->end() )
             &&
             ( initialSparticleMass
               > (*setIterator)->getInitialSparticle()->get_absolute_mass() ) )
        // we go through the list of squarks until we get to either the end of
        // the list or run out of squarks lighter than the gluino:
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
                                 compoundCascades.addNewAtEnd()->setProperties(
                                                              initialSparticle,
                                                              *cascadeIterator,
                                                  electroweakCascadeSource ) );
            }
          }
        }
        // now we move on to the next potential squark cascade:
        ++setIterator;
      }
    }



    namespace gluinoOrElectroweakinoSetType
    {
      gluinoSet::gluinoSet( inputHandler const* const inputShortcut,
           electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource,
                            fullCascadeSetOrderer* const setOrderer,
                            double const beamEnergy ) :
          gluinoOrNeutralinoSet( inputShortcut,
                                 inputShortcut->getGluino(),
                                 electroweakCascadeSource,
                                 beamEnergy ),
      directToEwinoCascades()
      {
        // we have to set up the cascades directly to electroweakinos now:
        for( std::vector< particlePointer >::const_iterator
             ewinoIterator( inputShortcut->getElectroweakinos()->begin() );
             inputShortcut->getElectroweakinos()->end() > ewinoIterator;
             ++ewinoIterator )
        {
          directToEwinoCascades.push_back(
                              new fullCascadeType::gluinoDirectlyToElectroweak(
                                                                 inputShortcut,
                                                                    beamEnergy,
                                       electroweakCascadeSource->getCascadeSet(
                                                              initialSparticle,
                                                          *ewinoIterator ) ) );
        }
      }

      gluinoSet::~gluinoSet()
      {
        for( std::vector< fullCascadeType::gluinoDirectlyToElectroweak*
                                                                    >::iterator
             deletionIterator( directToEwinoCascades.begin() );
             directToEwinoCascades.end() > deletionIterator;
             ++deletionIterator )
        {
          delete *deletionIterator;
        }
      }

    }  // end of gluinoOrNeutralinoSetType namespace

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
      squarkCascadeSets.push_back( new sdownType( inputShortcut,
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
    for( std::vector< sdownType* >::iterator
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
      for( std::vector< sdownType* >::iterator
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
