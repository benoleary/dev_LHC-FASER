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
 *      LHC-FASER also requires grids of lookup acceptanceValues. These should also be
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
                            double const chargeConjugateSumFactor ) :
    inputShortcut( NULL ),
    initialScolored( NULL ),
    beamEnergy( CppSLHA::CppSLHA_global::really_wrong_value ),
    typeOfColorfulCascade( typeOfColorfulCascade ),
    firstDecayBodyNumber( firstDecayBodyNumber ),
    chargeConjugateSumFactor( chargeConjugateSumFactor ),
    ewinoFlipsCharge( false ),
    bosonFlipsCharge( false ),
    subcascadePointer( NULL ),
    cascadeDefiner(),
    cascadeSegment( NULL ),
    ewinoCascades( NULL ),
    bosonCascades( NULL ),
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

  double
  fullCascade::getCombinedAcceptance( bool const bosonChargeNotFlipped,
                                      bool const ewinoChargeNotFlipped,
                                      acceptanceCutSet* const acceptanceCuts,
                                      int const numberOfAdditionalJets,
                                      int const numberOfNegativeElectrons,
                                      int const numberOfPositiveElectrons,
                                      int const numberOfNegativeMuons,
                                      int const numberOfPositiveMuons )
  /* this combines the acceptances from ewinoCascade & vectorCascade. Whether
   * the relevant scoloreds are particles or antiparticles should be decided
   * by the derived class using this function.
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
      double returnDouble( 0.0 );
      double ewinoAcceptance;
      double bosonAcceptance;
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
                if( ewinoChargeNotFlipped )
                {
                  ewinoAcceptance
                  = ewinoCascades->getAcceptance( acceptanceCuts,
                                                  ewinoJets,
                                                  ewinoNegativeElectrons,
                                                  ewinoPositiveElectrons,
                                                  ewinoNegativeMuons,
                                                  ewinoPositiveMuons );
                }
                else
                {
                  ewinoAcceptance
                  = ewinoCascades->getAcceptance( acceptanceCuts,
                                                  ewinoJets,
                                                  ewinoPositiveElectrons,
                                                  ewinoNegativeElectrons,
                                                  ewinoPositiveMuons,
                                                  ewinoNegativeMuons );
                }
                if( bosonChargeNotFlipped )
                {
                  bosonAcceptance
                  = bosonCascades->getAcceptance( acceptanceCuts,
                                        ( numberOfAdditionalJets - ewinoJets ),
                        ( numberOfNegativeElectrons - ewinoNegativeElectrons ),
                        ( numberOfPositiveElectrons - ewinoPositiveElectrons ),
                                ( numberOfNegativeMuons - ewinoNegativeMuons ),
                              ( numberOfPositiveMuons - ewinoPositiveMuons ) );
                }
                else
                {
                  bosonAcceptance
                  = bosonCascades->getAcceptance( acceptanceCuts,
                                        ( numberOfAdditionalJets - ewinoJets ),
                        ( numberOfPositiveElectrons - ewinoPositiveElectrons ),
                        ( numberOfNegativeElectrons - ewinoNegativeElectrons ),
                                ( numberOfPositiveMuons - ewinoPositiveMuons ),
                              ( numberOfNegativeMuons - ewinoNegativeMuons ) );
                }
              returnDouble += ( ewinoAcceptance * bosonAcceptance );
              }  // end of loop over positive electrons.
            }  // end of loop over negative muons.
          }  // end of loop over positive electrons.
        }  // end of loop over negative electrons.
      }  // end of loop over jets.
      return returnDouble;
    }  // end of if numberOfLeptons was in the allowed range.
  }



  sxFullCascade::sxFullCascade() :
    fullCascade( sx,
                 2,
                 1.0 )
  {
    // just an initialization list.
  }

  sxFullCascade::~sxFullCascade()
  {
    // does nothing.
  }


  void
  sxFullCascade::setProperties( inputHandler const* const inputShortcut,
                                particlePointer const initialSquark,
                                double const beamEnergy,
                                electroweakCascadeSet* const ewinoCascades )
  {
    this->inputShortcut = inputShortcut;
    this->initialScolored = initialSquark;
    this->beamEnergy = beamEnergy;
    this->ewinoCascades = ewinoCascades;
    int ewinoPdgCode( ewinoCascade->getElectroweakDecayer()->get_PDG_code() );
    if( ( inputShortcut->isIn( ewinoPdgCode,
                               inputShortcut->getCharginos() ) )
        &&
        ( inputShortcut->isIn( initialSquark->get_PDG_code(),
                               inputShortcut->getSdownTypes() ) ) )
      // if we have to worry about which sign of PDG code to use...
    {
      ewinoFlipsCharge = true;
      soughtDecayProductList.front() = -(ewinoPdgCode);
    }
    else
    {
      ewinoFlipsCharge = false;
      soughtDecayProductList.front() = ewinoPdgCode;
    }
    resetCachedBranchingRatio();
    cascadeDefiner.clearEntries();
    // reset cascadeDefiner.
    cascadeSegment = cascadeDefiner.addNewAtEnd();
    cascadeSegment->first = ewinoCascade->getElectroweakDecayer();
    // sx means only 1 decay to be recorded.
    cascadeSegment->second = firstDecayBodyNumber;
    // sx also means that the decay is 2-body.
  }



  gxFullCascade::gxFullCascade() :
      fullCascade( gx,
                   3,
                   2.0 )
  {
    // just an initialization list.
  }

  gxFullCascade::~gxFullCascade()
  {
    // does nothing.
  }

  double
  gxFullCascade::getAcceptance( bool const scoloredIsNotAntiparticle,
                                acceptanceCutSet* const acceptanceCuts,
                                int const numberOfAdditionalJets,
                                int const numberOfNegativeElectrons,
                                int const numberOfPositiveElectrons,
                                int const numberOfNegativeMuons,
                                int const numberOfPositiveMuons )
  /* this calls the appropriate functions on ewinoCascade, summing over the 2
   * charge versions if the electroweakino is a chargino. it ignores the value
   * of scoloredIsNotAntiparticle.
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
      // the summing over charge-conjugate decays is already included in the BR
      // to the electroweakino.
    }
  }


  supxFullCascade::supxFullCascade() :
      fullCascade( sx,
                   3,
                   1.0 )
  {
    // just an initialization list.
  }

  supxFullCascade::~sbxFullCascade()
  {
    // does nothing.
  }


  void
  supxFullCascade::setProperties( inputHandler const* const inputShortcut,
                                  particlePointer const initialSquark,
                                  double const beamEnergy,
                              electroweakCascadeSet* const directEwinoCascades,
                               electroweakCascadeSet* const ewinoWithWCascades,
                                  electroweakCascadeSet* const bosonCascades )
  {
    this->inputShortcut = inputShortcut;
    this->initialScolored = initialSquark;
    this->beamEnergy = beamEnergy;
    this->directEwinoCascades = directEwinoCascades;
    this->ewinoWithWCascades = ewinoWithWCascades;
    int ewinoPdgCode( ewinoCascade->getElectroweakDecayer()->get_PDG_code() );
    if( ( inputShortcut->isIn( ewinoPdgCode,
                               inputShortcut->getCharginos() ) )
        &&
        ( inputShortcut->isIn( initialSquark->get_PDG_code(),
                               inputShortcut->getSdownTypes() ) ) )
      // if we have to worry about which sign of PDG code to use...
    {
      ewinoFlipsCharge = true;
      soughtDecayProductList.front() = -(ewinoPdgCode);
    }
    else
    {
      ewinoFlipsCharge = false;
      soughtDecayProductList.front() = ewinoPdgCode;
    }
    resetCachedBranchingRatio();
    cascadeDefiner.clearEntries();
    // reset cascadeDefiner.
    cascadeSegment = cascadeDefiner.addNewAtEnd();
    cascadeSegment->first = ewinoCascade->getElectroweakDecayer();
    // sx means only 1 decay to be recorded.
    cascadeSegment->second = firstDecayBodyNumber;
    // sx also means that the decay is 2-body.
  }



  sjgxFullCascade::sjgxFullCascade() :
    fullCascade( sjgx,
                 2,
                 1.0 )
  // sjgx also means that the initial decay is 2-body.
  {
    // just an initialization list.
  }

  sjgxFullCascade::~sjgxFullCascade()
  {
    // does nothing.
  }



  gjsxFullCascade::gjsxFullCascade() :
    fullCascade( gjsx,
                 2,
                 2.0 )
  // gjsx also means that the initial decay is 2-body.
  {
    // just an initialization list.
  }

  gjsxFullCascade::~gjsxFullCascade()
  {
    // does nothing.
  }



  sjgjsxFullCascade::sjgjsxFullCascade() :
    fullCascade( sjgjsx,
                 2,
                 1.0 )
  // sjgjsx also means that the initial decay is 2-body.
  {
    // just an initialization list.
  }

  sjgjsxFullCascade::~sjgjsxFullCascade()
  {
    // does nothing.
  }



  sbsxFullCascade::sbsxFullCascade() :
    fullCascade( sbsx,
                 2,
                 1.0 )
  // sbsx also means that the initial decay is 2-body.
  {
    soughtDecayProductList.push_back(
                                 CppSLHA::CppSLHA_global::really_wrong_value );
    // a cascade beginning with sv needs to look for 2 specific decay products.
  }

  sbsxFullCascade::~sbsxFullCascade()
  {
    // does nothing.
  }



  gbsxFullCascade::gbsxFullCascade() :
    fullCascade( gbsx,
                 3,
                 2.0 )
  // gbsx also means that the initial decay is 3-body.
  {
    soughtDecayProductList.push_back(
                                 CppSLHA::CppSLHA_global::really_wrong_value );
    // a cascade beginning with gv needs to look for 2 specific decay products.
  }

  gbsxFullCascade::~gbsxFullCascade()
  {
    // does nothing.
  }



  gjsbsxFullCascade::gjsbsxFullCascade() :
    fullCascade( gjsbsx,
                 2,
                 2.0 )
  // gjsbsx also means that the initial decay is 2-body.
  {
    // just an initialization list.
  }

  gjsbsxFullCascade::~gjsbsxFullCascade()
  {
    // does nothing.
  }



  sbgxFullCascade::sbgxFullCascade() :
    fullCascade( sbgx,
                 3,
                 1.0 )
  // sbgx also means that the initial decay is 3-body.
  {
    // just an initialization list.
  }

  sbgxFullCascade::~sbgxFullCascade()
  {
    // does nothing.
  }



  sbsjgxFullCascade::sbsjgxFullCascade() :
    fullCascade( sbsjgx,
                 2,
                 1.0 )
  // sbsjgx also means that the initial decay is 2-body.
  {
    // just an initialization list.
  }

  sbsjgxFullCascade::~sbsjgxFullCascade()
  {
    // does nothing.
  }



  sbsjgjsxFullCascade::sbsjgjsxFullCascade() :
    fullCascade( sbsjgjsx,
                 2,
                 1.0 )
  // sbsjgjsx also means that the initial decay is 2-body.
  {
    // just an initialization list.
  }

  sbsjgjsxFullCascade::~sbsjgjsxFullCascade()
  {
    // does nothing.
  }



  sbgjsxFullCascade::sbgjsxFullCascade() :
    fullCascade( sbgjsx,
                 3,
                 1.0 )
  // sbgjsx also means that the initial decay is 3-body.
  {
    // just an initialization list.
  }

  sbgjsxFullCascade::~sbgjsxFullCascade()
  {
    // does nothing.
  }



  sjgbsxFullCascade::sjgbsxFullCascade() :
    fullCascade( sjgbsx,
                 2,
                 1.0 )
  // sjgbsx also means that the initial decay is 2-body.
  {
    // just an initialization list.
  }

  sjgbsxFullCascade::~sjgbsxFullCascade()
  {
    // does nothing.
  }



  sjgjsbsxFullCascade::sjgjsbsxFullCascade() :
    fullCascade( sjgjsbsx,
                 2,
                 1.0 )
  // sjgjsbsx also means that the initial decay is 2-body.
  {
    // just an initialization list.
  }

  sjgjsbsxFullCascade::~sjgjsbsxFullCascade()
  {
    // does nothing.
  }



  fullCascadeSet::fullCascadeSet( inputHandler const* const inputShortcut,
                                  particlePointer const initialScolored,
           electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource,
                           readiableCascadeSetList* const squarkCascadeSetList,
                                  double const beamEnergy ) :
    getsReadiedForNewPoint( inputShortcut->getReadier() ),
    inputShortcut( inputShortcut ),
    initialScolored( initialScolored ),
    electroweakCascadeSource( electroweakCascadeSource ),
    squarkCascadeSetList( squarkCascadeSetList ),
    //squarkCascadeSetListNotYetOrdered( true ),
    beamEnergy( beamEnergy )
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
  // fullCascades that initialScolored has for this point.
  {
    openCascades.clear();
    openSxCascades.clear();
    openGxCascades.clear();
    for( std::vector< sxFullCascade* >::iterator
         sxIterator( sxCascades.getVector()->begin() );
         sxCascades.getVector()->end() > sxIterator;
         ++sxIterator )
    {
      if( (*sxIterator)->isOpen() )
      {
        addOpenCascade( *sxIterator );
        openSxCascades.push_back( *sxIterator );
      }
    }
    for( std::vector< gxFullCascade* >::iterator
         gxIterator( gxCascades.getVector()->begin() );
         gxCascades.getVector()->end() > gxIterator;
         ++gxIterator )
    {
      if( (*gxIterator)->isOpen() )
      {
        addOpenCascade( *gxIterator );
        openGxCascades.push_back( *gxIterator );
      }
    }
    if( squarkCascadeSetList->second->publicNeedsToPrepareForThisPoint() )
    {
      squarkCascadeSetList->first->sort( &massOrdered );
    }
    // now squarkCascadeSetList has been ordered such that the cascades of the
    // lighter squarks come before those of the heavier squarks.
    buildLongerCascades();
  }



  squarkFullCascadeSet::squarkFullCascadeSet(
                                       inputHandler const* const inputShortcut,
           electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource,
                                         particlePointer const initialScolored,
                           readiableCascadeSetList* const squarkCascadeSetList,
                                        fullCascadeSet* const gluinoCascadeSet,
                                              double const beamEnergy ) :
    fullCascadeSet( inputShortcut,
                    initialScolored,
                    electroweakCascadeSource,
                    squarkCascadeSetList,
                    beamEnergy ),
    gluinoCascadeSet( gluinoCascadeSet )
  {
    // we have to set up sxCascades now:
    for( std::vector< particlePointer >::const_iterator
         ewinoIterator( inputShortcut->getElectroweakinos()->begin() );
         inputShortcut->getElectroweakinos()->end() > ewinoIterator;
         ++ewinoIterator )
    {
      sxCascades.addNewAtEnd()->setProperties( inputShortcut,
                                               initialScolored,
                                               beamEnergy,
                                     electroweakCascadeSource->getCascadeSet(
                                                               initialScolored,
                                                            *ewinoIterator ) );
    }
  }

  squarkFullCascadeSet::~squarkFullCascadeSet()
  {
    // does nothing.
  }


  void
  squarkFullCascadeSet::buildLongerCascades()
  {
    sjgxCascades.clearEntries();
    sjgjsxCascades.clearEntries();
    sjgbsxCascades.clearEntries();
    sjgjsbsxCascades.clearEntries();
    sbgxCascades.clearEntries();
    sbgjsxCascades.clearEntries();
    if( initialScolored->get_absolute_mass() > inputShortcut->getGluinoMass() )
    {
      for( std::vector< gxFullCascade* >::const_iterator
           gxIterator( gluinoCascadeSet->getGxCascades()->begin() );
           gluinoCascadeSet->getGxCascades()->end() > gxIterator;
           ++gxIterator )
      {
        sjgxCascades.addNewAtEnd()->setProperties( initialScolored,
                                                   *gxIterator );
        addOpenCascade( sjgxCascades.getBack() );
      }
      for( std::vector< gjsxFullCascade* >::const_iterator
           gjsxIterator( gluinoCascadeSet->getGjsxCascades()->begin() );
           gluinoCascadeSet->getGjsxCascades()->end()
           > gjsxIterator;
           ++gjsxIterator )
      {
        sjgjsxCascades.addNewAtEnd()->setProperties( initialScolored,
                                                     *gjsxIterator );
        addOpenCascade( sjgjsxCascades.getBack() );
      }
      for( std::vector< gbsxFullCascade* >::const_iterator
           gvsxIterator( gluinoCascadeSet->getGbsxCascades()->begin() );
           gluinoCascadeSet->getGbsxCascades()->end() > gvsxIterator;
           ++gvsxIterator )
      {
        sjgbsxCascades.addNewAtEnd()->setProperties( initialScolored,
                                                     *gvsxIterator );
        addOpenCascade( sjgbsxCascades.getBack() );
      }
      for( std::vector< gjsbsxFullCascade* >::const_iterator
           gjsvsxIterator( gluinoCascadeSet->getGjsbsxCascades()->begin() );
           gluinoCascadeSet->getGjsbsxCascades()->end()
           > gjsvsxIterator;
           ++gjsvsxIterator )
      {
        sjgjsbsxCascades.addNewAtEnd()->setProperties( initialScolored,
                                                       *gjsvsxIterator );
        addOpenCascade( sjgjsbsxCascades.getBack() );
      }
      for( std::vector< particlePointer >::const_iterator
           ewIterator( inputShortcut->getEwsbBosonsAndMassiveVectorBosons(
                                                                  )->begin() );
           inputShortcut->getEwsbBosonsAndMassiveVectorBosons()->end()
           > ewIterator;
           ++ewIterator )
      {
        soughtPositivePdgCodeList.clear();
        soughtPositivePdgCodeList.push_back( (*ewIterator)->get_PDG_code() );
        soughtPositivePdgCodeList.push_back( CppSLHA::PDG_code::gluino );
        soughtNegativePdgCodeList.clear();
        soughtNegativePdgCodeList.push_back(
                                            -((*ewIterator)->get_PDG_code()) );
        soughtNegativePdgCodeList.push_back( CppSLHA::PDG_code::gluino );
        if( ( initialScolored->get_absolute_mass()
              > ( inputShortcut->getGluinoMass()
                  + (*ewIterator)->get_absolute_mass() ) )
            &&
            ( ( lhcFaserGlobal::negligibleBr
                < initialScolored->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                                 &soughtPositivePdgCodeList ) )
              || ( lhcFaserGlobal::negligibleBr
                               < initialScolored->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                             &soughtNegativePdgCodeList ) ) ) )
        {
          for( std::vector< gxFullCascade* >::const_iterator
               gxIterator( gluinoCascadeSet->getGxCascades()->begin() );
               gluinoCascadeSet->getGxCascades()->end() > gxIterator;
               ++gxIterator )
          {
            sbgxCascades.addNewAtEnd()->setProperties( initialScolored,
                                     electroweakCascadeSource->getCascadeSet(
                                                               initialScolored,
                                                                 *ewIterator,
                                                  inputShortcut->getGluino() ),
                                                       *gxIterator );
            addOpenCascade( sbgxCascades.getBack() );
          }
          for( std::vector< gjsxFullCascade* >::const_iterator
               gjsxIterator( gluinoCascadeSet->getGjsxCascades()->begin() );
               gluinoCascadeSet->getGjsxCascades()->end() > gjsxIterator;
               ++gjsxIterator )
          {
            sbgjsxCascades.addNewAtEnd()->setProperties( initialScolored,
                                       electroweakCascadeSource->getCascadeSet(
                                                               initialScolored,
                                                                 *ewIterator ),
                                                         *gjsxIterator );
            addOpenCascade( sbgjsxCascades.getBack() );
          }
        }
      }  // end of loop over if svg is open
    }  // end of if squark mass is greater than gluino mass

    // now we look for decays into lighter squarks (squarkCascadeSetList has
    // been ordered for this point already in fullCascadeSet::setUpCascades()):
    sbsxCascades.clearEntries();
    sbsjgxCascades.clearEntries();
    sbsjgjsxCascades.clearEntries();
    for( std::list< fullCascadeSet* >::iterator
         cascadeIterator( squarkCascadeSetList->first->begin() );
         (*cascadeIterator)->getInitialScolored()->get_absolute_mass()
         < initialScolored->get_absolute_mass();
         // this conditional should prevent the iterator falling off the end
         // of the list provided that the list is as it should be.
         ++cascadeIterator )
    {
      for( std::vector< particlePointer >::const_iterator
           ewIterator( inputShortcut->getEwsbBosonsAndMassiveVectorBosons(
                                                                  )->begin() );
           inputShortcut->getEwsbBosonsAndMassiveVectorBosons()->end()
           > ewIterator;
           ++ewIterator )
      {
        soughtPositivePdgCodeList.clear();
        soughtPositivePdgCodeList.push_back(
                    (*cascadeIterator)->getInitialScolored()->get_PDG_code() );
        soughtPositivePdgCodeList.push_back( (*ewIterator)->get_PDG_code() );
        soughtNegativePdgCodeList.clear();
        soughtNegativePdgCodeList.push_back(
                    (*cascadeIterator)->getInitialScolored()->get_PDG_code() );
        soughtNegativePdgCodeList.push_back(
                                            -((*ewIterator)->get_PDG_code()) );
        if( ( initialScolored->get_absolute_mass()
              > ( (*cascadeIterator)->getInitialScolored()->get_absolute_mass()
                  + (*ewIterator)->get_absolute_mass() ) )
            &&
            ( ( lhcFaserGlobal::negligibleBr
                < initialScolored->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                                 &soughtPositivePdgCodeList ) )
              || ( lhcFaserGlobal::negligibleBr
                               < initialScolored->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                             &soughtNegativePdgCodeList ) ) ) )
        {
          for( std::vector< sxFullCascade* >::const_iterator
               sxIterator( (*cascadeIterator)->getSxCascades()->begin() );
               (*cascadeIterator)->getSxCascades()->end() > sxIterator;
               ++sxIterator )
          {
            sbsxCascades.addNewAtEnd()->setProperties( initialScolored,
                                       electroweakCascadeSource->getCascadeSet(
                                                               initialScolored,
                                                                   *ewIterator,
                                    (*cascadeIterator)->getInitialScolored() ),
                                                       *sxIterator );
            addOpenCascade( sbsxCascades.getBack() );
          }
          for( std::vector< sjgxFullCascade* >::const_iterator
               sjgxIterator( (*cascadeIterator)->getSjgxCascades()->begin() );
               (*cascadeIterator)->getSjgxCascades()->end() > sjgxIterator;
               ++sjgxIterator )
          {
            sbsjgxCascades.addNewAtEnd()->setProperties( initialScolored,
                                       electroweakCascadeSource->getCascadeSet(
                                                               initialScolored,
                                                                   *ewIterator,
                                    (*cascadeIterator)->getInitialScolored() ),
                                                         *sjgxIterator );
            addOpenCascade( sbsjgxCascades.getBack() );
          }
          for( std::vector< sjgjsxFullCascade* >::const_iterator
               sjgjsxIterator( (*cascadeIterator)->getSjgjsxCascades(
                                                                  )->begin() );
               (*cascadeIterator)->getSjgjsxCascades()->end() > sjgjsxIterator;
               ++sjgjsxIterator )
          {
            sbsjgjsxCascades.addNewAtEnd()->setProperties( initialScolored,
                                       electroweakCascadeSource->getCascadeSet(
                                                               initialScolored,
                                                                   *ewIterator,
                                    (*cascadeIterator)->getInitialScolored() ),
                                                           *sjgjsxIterator );
            addOpenCascade( sbsjgjsxCascades.getBack() );
          }
        }
      }  // end of loop over if svs is open
    }
  }



  gluinoFullCascadeSet::gluinoFullCascadeSet(
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
                                                               initialScolored,
                                                            *ewinoIterator ) );
    }
  }

  gluinoFullCascadeSet::~gluinoFullCascadeSet()
  {
    // does nothing.
  }


  void
  gluinoFullCascadeSet::buildLongerCascades()
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
           ( (*cascadeIterator)->getInitialScolored()->get_absolute_mass()
             < initialScolored->get_absolute_mass() ) );
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
                    (*cascadeIterator)->getInitialScolored()->get_PDG_code() );
        soughtPositivePdgCodeList.push_back( (*ewIterator)->get_PDG_code() );
        soughtNegativePdgCodeList.clear();
        soughtNegativePdgCodeList.push_back(
                    (*cascadeIterator)->getInitialScolored()->get_PDG_code() );
        soughtNegativePdgCodeList.push_back(
                                            -((*ewIterator)->get_PDG_code()) );
        if( ( initialScolored->get_absolute_mass()
              > ( (*cascadeIterator)->getInitialScolored()->get_absolute_mass()
                  + (*ewIterator)->get_absolute_mass() ) )
            &&
            ( ( lhcFaserGlobal::negligibleBr
                < initialScolored->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                                 &soughtPositivePdgCodeList ) )
              || ( lhcFaserGlobal::negligibleBr
                               < initialScolored->inspect_direct_decay_handler(
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
                                                               initialScolored,
                                                                 *ewIterator,
                                    (*cascadeIterator)->getInitialScolored() ),
                                                       *sxIterator );
            addOpenCascade( gbsxCascades.getBack() );
          }
        }
      }  // end of loop over if gvs is open
    }  // end of loop over squarkFullCascadeSets
  }



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
    gluinoCascadeSet = new gluinoFullCascadeSet( inputShortcut,
                                                 electroweakCascadeSource,
                                                 &squarkCascadeSetList,
                                                 beamEnergy );
    for( std::vector< particlePointer >::const_iterator
         squarkIterator( inputShortcut->getSquarks()->begin() );
         inputShortcut->getSquarks()->end() > squarkIterator;
         ++squarkIterator )
    {
      squarkCascadeSets.push_back( new squarkFullCascadeSet( inputShortcut,
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
    for( std::vector< squarkFullCascadeSet* >::iterator
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
      for( std::vector< squarkFullCascadeSet* >::iterator
           searchIterator( squarkCascadeSets.begin() );
           squarkCascadeSets.end() > searchIterator;
           ++searchIterator )
      {
        if( initialScolored == (*searchIterator)->getInitialScolored() )
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
