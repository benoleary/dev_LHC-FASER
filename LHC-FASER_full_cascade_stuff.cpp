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
    shortcut( NULL ),
    initialScolored( NULL ),
    beamEnergy( CppSLHA::CppSLHA_global::really_wrong_value ),
    typeOfColorfulCascade( typeOfColorfulCascade ),
    firstDecayBodyNumber( firstDecayBodyNumber ),
    chargeConjugateSumFactor( chargeConjugateSumFactor ),
    subcascade( NULL ),
    cascadeDefiner(),
    ewinoCascades( NULL ),
    vectorCascades( NULL ),
    soughtDecayProductList( 1,
                            CppSLHA::CppSLHA_global::really_wrong_value ),
    cachedBranchingRatio( CppSLHA::CppSLHA_global::really_wrong_value ),
    branchingRatioNeedsToBeReCalculated( true )
  {
    // just an initialization list.
  }

  fullCascade::~fullCascade()
  {
    // does nothing.
  }


  double
  fullCascade::specifiedJetsOneOssfMinusOsdfPair( acceptanceCutSet* const cuts,
                                             int const numberOfAdditionalJets )
  /* this calls
   * ewinoCascade->getOssfMinusOsdf( scoloredIsNotAntiparticle,
   *                                 cuts ),
   * & if vectorCascade is not NULL, it calls
   * ewinoCascade->getAcceptance( [ { 0, 1, 2 } jets + 0 leptons for cuts ] )
   * & the same for vectorCascade, as well as ->getOssfMinusOsdf( cuts ).
   */
  {
    // we should really check that ewinoCascade is not NULL, but this function
    // should never be called on a fullCascade which has a NULL ewinoCascade.
    if( NULL == vectorCascades )
    {
      /* if there is no possibility of extra jets from a vector decay, we can
       * only get a single OSSF-OSDF pair from an electroweakino decay, which
       * cannot produce additional jets.
       */
      if( 0 == numberOfAdditionalJets )
      {
        return ewinoCascades->getOssfMinusOsdf( cuts );
      }
      else
      {
        return 0.0;
      }
    }
    else
    {
      return ( ewinoCascades->getOssfMinusOsdf( cuts )
               * vectorCascades->getAcceptance( cuts,
                                                numberOfAdditionalJets,
                                                0,
                                                0,
                                                0,
                                                0 )
               + vectorCascades->getOssfMinusOsdf( cuts )
                 * ewinoCascades->getAcceptance( cuts,
                                                 numberOfAdditionalJets,
                                                 0,
                                                 0,
                                                 0,
                                                 0 ) );
    }
  }

  double
  fullCascade::unspecifiedJetsSpecifiedOssfMinusOsdfPairs(
                                                  acceptanceCutSet* const cuts,
                                                int const numberOfLeptonPairs )
  // this only calls ewinoCascade->getOssfMinusOsdf( cuts ),
  // & if vectorCascade is not NULL, vectorCascade->getOssfMinusOsdf( cuts ).
  {
    if( 0 == numberOfLeptonPairs )
    {
      return unspecifiedJetsSpecifiedChargeSummedLeptons( cuts,
                                                          0 );
    }
    else
    {
      // we should check that ewinoCascade is not NULL, but this function
      // should never be called on a fullCascade which has a NULL ewinoCascade.
      if( NULL == vectorCascades )
      {
        /* if there is no possibility of extra jets from a vector decay, we can
         * only get a single OSSF-OSDF pair from an electroweakino decay, which
         * cannot produce additional jets.
         */
        if( 1 == numberOfLeptonPairs )
        {
          return ewinoCascades->getOssfMinusOsdf( cuts );
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
          return ( ewinoCascades->getOssfMinusOsdf( cuts )
                   * vectorCascades->getOssfMinusOsdf( cuts ) );
        }
        else if( 1 == numberOfLeptonPairs )
        {
        return ( ewinoCascades->getOssfMinusOsdf( cuts )
                 * ( vectorCascades->getAcceptance( cuts,
                                                    0,
                                                    0,
                                                    0,
                                                    0,
                                                    0 )
                     + vectorCascades->getAcceptance( cuts,
                                                      1,
                                                      0,
                                                      0,
                                                      0,
                                                      0 )
                     + vectorCascades->getAcceptance( cuts,
                                                      2,
                                                      0,
                                                      0,
                                                      0,
                                                      0 ) )
                 + vectorCascades->getOssfMinusOsdf( cuts )
                   * ( ewinoCascades->getAcceptance( cuts,
                                                     0,
                                                     0,
                                                     0,
                                                     0,
                                                     0 )
                       + ewinoCascades->getAcceptance( cuts,
                                                       1,
                                                       0,
                                                       0,
                                                       0,
                                                       0 )
                       + ewinoCascades->getAcceptance( cuts,
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
                                                  acceptanceCutSet* const cuts,
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
                                           cuts,
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
                                                  acceptanceCutSet* const cuts,
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
        returnDouble += specifiedJetsSpecifiedChargeSummedLeptons( cuts,
                                                                  numberOfJets,
                                                             numberOfLeptons );
      }  // end of loop over jets.
      return returnDouble;
    }  // end of if numberOfLeptons was in the allowed range.
  }

  double
  fullCascade::getCombinedAcceptance(
                                    bool const vectorScoloredIsNotAntiparticle,
                                     bool const ewinoScoloredIsNotAntiparticle,
                                      acceptanceCutSet* const cuts,
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
      double vectorAcceptance;
      for( int ewinoJets = numberOfAdditionalJets;
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
                if( ewinoScoloredIsNotAntiparticle )
                {
                  ewinoAcceptance = ewinoCascades->getAcceptance( cuts,
                                                                 ewinoJets,
                                                        ewinoNegativeElectrons,
                                                        ewinoPositiveElectrons,
                                                            ewinoNegativeMuons,
                                                          ewinoPositiveMuons );
                }
                else
                {
                  ewinoAcceptance = ewinoCascades->getAcceptance( cuts,
                                                                 ewinoJets,
                                                        ewinoPositiveElectrons,
                                                        ewinoNegativeElectrons,
                                                            ewinoPositiveMuons,
                                                          ewinoNegativeMuons );
                }
                if( vectorScoloredIsNotAntiparticle )
                {
                  vectorAcceptance = vectorCascades->getAcceptance( cuts,
                                        ( numberOfAdditionalJets - ewinoJets ),
                        ( numberOfNegativeElectrons - ewinoNegativeElectrons ),
                        ( numberOfPositiveElectrons - ewinoPositiveElectrons ),
                                ( numberOfNegativeMuons - ewinoNegativeMuons ),
                              ( numberOfPositiveMuons - ewinoPositiveMuons ) );
                }
                else
                {
                  vectorAcceptance = vectorCascades->getAcceptance( cuts,
                                        ( numberOfAdditionalJets - ewinoJets ),
                        ( numberOfPositiveElectrons - ewinoPositiveElectrons ),
                        ( numberOfNegativeElectrons - ewinoNegativeElectrons ),
                                ( numberOfPositiveMuons - ewinoPositiveMuons ),
                              ( numberOfNegativeMuons - ewinoNegativeMuons ) );
                }
              returnDouble += ( ewinoAcceptance * vectorAcceptance );
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
                                acceptanceCutSet* const cuts,
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
      return ewinoCascades->getAcceptance( cuts,
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
      return ( 0.5 * ( ewinoCascades->getAcceptance( cuts,
                                                     numberOfAdditionalJets,
                                                     numberOfNegativeElectrons,
                                                     numberOfPositiveElectrons,
                                                     numberOfNegativeMuons,
                                                     numberOfPositiveMuons )
                       + ewinoCascades->getAcceptance( cuts,
                                                       numberOfAdditionalJets,
                                                     numberOfPositiveElectrons,
                                                     numberOfNegativeElectrons,
                                                       numberOfPositiveMuons,
                                                   numberOfNegativeMuons ) ) );
      // the summing over charge-conjugate decays is already included in the BR
      // to the electroweakino.
    }
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



  svsxFullCascade::svsxFullCascade() :
    fullCascade( svsx,
                 2,
                 1.0 )
  // svsx also means that the initial decay is 2-body.
  {
    soughtDecayProductList.push_back(
                                 CppSLHA::CppSLHA_global::really_wrong_value );
    // a cascade beginning with sv needs to look for 2 specific decay products.
  }

  svsxFullCascade::~svsxFullCascade()
  {
    // does nothing.
  }



  gvsxFullCascade::gvsxFullCascade() :
    fullCascade( gvsx,
                 3,
                 2.0 )
  // gvsx also means that the initial decay is 3-body.
  {
    soughtDecayProductList.push_back(
                                 CppSLHA::CppSLHA_global::really_wrong_value );
    // a cascade beginning with gv needs to look for 2 specific decay products.
  }

  gvsxFullCascade::~gvsxFullCascade()
  {
    // does nothing.
  }



  gjsvsxFullCascade::gjsvsxFullCascade() :
    fullCascade( gjsvsx,
                 2,
                 2.0 )
  // gjsvsx also means that the initial decay is 2-body.
  {
    // just an initialization list.
  }

  gjsvsxFullCascade::~gjsvsxFullCascade()
  {
    // does nothing.
  }



  svgxFullCascade::svgxFullCascade() :
    fullCascade( svgx,
                 3,
                 1.0 )
  // svgx also means that the initial decay is 3-body.
  {
    // just an initialization list.
  }

  svgxFullCascade::~svgxFullCascade()
  {
    // does nothing.
  }



  svsjgxFullCascade::svsjgxFullCascade() :
    fullCascade( svsjgx,
                 2,
                 1.0 )
  // svsjgx also means that the initial decay is 2-body.
  {
    // just an initialization list.
  }

  svsjgxFullCascade::~svsjgxFullCascade()
  {
    // does nothing.
  }



  svsjgjsxFullCascade::svsjgjsxFullCascade() :
    fullCascade( svsjgjsx,
                 2,
                 1.0 )
  // svsjgjsx also means that the initial decay is 2-body.
  {
    // just an initialization list.
  }

  svsjgjsxFullCascade::~svsjgjsxFullCascade()
  {
    // does nothing.
  }



  svgjsxFullCascade::svgjsxFullCascade() :
    fullCascade( svgjsx,
                 3,
                 1.0 )
  // svgjsx also means that the initial decay is 3-body.
  {
    // just an initialization list.
  }

  svgjsxFullCascade::~svgjsxFullCascade()
  {
    // does nothing.
  }



  sjgvsxFullCascade::sjgvsxFullCascade() :
    fullCascade( sjgvsx,
                 2,
                 1.0 )
  // sjgvsx also means that the initial decay is 2-body.
  {
    // just an initialization list.
  }

  sjgvsxFullCascade::~sjgvsxFullCascade()
  {
    // does nothing.
  }



  sjgjsvsxFullCascade::sjgjsvsxFullCascade() :
    fullCascade( sjgjsvsx,
                 2,
                 1.0 )
  // sjgjsvsx also means that the initial decay is 2-body.
  {
    // just an initialization list.
  }

  sjgjsvsxFullCascade::~sjgjsvsxFullCascade()
  {
    // does nothing.
  }



  fullCascadeSet::fullCascadeSet( inputHandler const* const shortcut,
                                  particlePointer const initialScolored,
           electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource,
                           readiableCascadeSetList* const squarkCascadeSetList,
                                  double const beamEnergy ) :
    getsReadiedForNewPoint( shortcut->getReadier() ),
    shortcut( shortcut ),
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
      squarkCascadeSetList->second->publicFinishPreparingForThisPoint();
    }
    // now squarkCascadeSetList has been ordered such that the cascades of the
    // lighter squarks come before those of the heavier squarks.
    buildLongerCascades();
  }



  squarkFullCascadeSet::squarkFullCascadeSet(
                                            inputHandler const* const shortcut,
           electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource,
                                         particlePointer const initialScolored,
                           readiableCascadeSetList* const squarkCascadeSetList,
                                        fullCascadeSet* const gluinoCascadeSet,
                                              double const beamEnergy ) :
    fullCascadeSet( shortcut,
                    initialScolored,
                    electroweakCascadeSource,
                    squarkCascadeSetList,
                    beamEnergy ),
    gluinoCascadeSet( gluinoCascadeSet )
  {
    // we have to set up sxCascades now:
    for( std::vector< particlePointer >::const_iterator
         ewinoIterator( shortcut->getElectroweakinos()->begin() );
         shortcut->getUnstableElectroweakinos()->end() > ewinoIterator;
         ++ewinoIterator )
    {
      sxCascades.addNewAtEnd()->setProperties( shortcut,
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
    sjgvsxCascades.clearEntries();
    sjgjsvsxCascades.clearEntries();
    svgxCascades.clearEntries();
    svgjsxCascades.clearEntries();
    if( initialScolored->get_absolute_mass() > shortcut->getGluinoMass() )
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
      for( std::vector< gvsxFullCascade* >::const_iterator
           gvsxIterator( gluinoCascadeSet->getGvsxCascades()->begin() );
           gluinoCascadeSet->getGvsxCascades()->end()
           > gvsxIterator;
           ++gvsxIterator )
      {
        sjgvsxCascades.addNewAtEnd()->setProperties( initialScolored,
                                                     *gvsxIterator );
        addOpenCascade( sjgvsxCascades.getBack() );
      }
      for( std::vector< gjsvsxFullCascade* >::const_iterator
           gjsvsxIterator( gluinoCascadeSet->getGjsvsxCascades()->begin() );
           gluinoCascadeSet->getGjsvsxCascades()->end()
           > gjsvsxIterator;
           ++gjsvsxIterator )
      {
        sjgjsvsxCascades.addNewAtEnd()->setProperties( initialScolored,
                                                       *gjsvsxIterator );
        addOpenCascade( sjgjsvsxCascades.getBack() );
      }
      for( std::vector< particlePointer >::const_iterator
           ewIterator( shortcut->getEwsbBosonsAndMassiveVectorBosons(
                                                                  )->begin() );
           shortcut->getEwsbBosonsAndMassiveVectorBosons()->end() > ewIterator;
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
              > ( shortcut->getGluinoMass()
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
               gluinoCascadeSet->getGxCascades()->end()
               > gxIterator;
               ++gxIterator )
          {
            svgxCascades.addNewAtEnd()->setProperties( initialScolored,
                                     electroweakCascadeSource->getCascadeSet(
                                                               initialScolored,
                                                                 *ewIterator,
                                                       shortcut->getGluino() ),
                                                       *gxIterator );
            addOpenCascade( svgxCascades.getBack() );
          }
          for( std::vector< gjsxFullCascade* >::const_iterator
               gjsxIterator( gluinoCascadeSet->getGjsxCascades()->begin() );
               gluinoCascadeSet->getGjsxCascades()->end()
               > gjsxIterator;
               ++gjsxIterator )
          {
            svgjsxCascades.addNewAtEnd()->setProperties( initialScolored,
                                       electroweakCascadeSource->getCascadeSet(
                                                               initialScolored,
                                                                 *ewIterator ),
                                                         *gjsxIterator );
            addOpenCascade( svgjsxCascades.getBack() );
          }
        }
      }  // end of loop over if svg is open
    }  // end of if squark mass is greater than gluino mass

    // now we look for decays into lighter squarks (squarkCascadeSetList has
    // been ordered for this point already in fullCascadeSet::setUpCascades()):
    svsxCascades.clearEntries();
    svsjgxCascades.clearEntries();
    svsjgjsxCascades.clearEntries();
    for( std::list< fullCascadeSet* >::iterator
         cascadeIterator( squarkCascadeSetList->first->begin() );
         (*cascadeIterator)->getInitialScolored()->get_absolute_mass()
         < initialScolored->get_absolute_mass();
         // this conditional should prevent the iterator falling off the end
         // of the list provided that the list is as it should be.
         ++cascadeIterator )
    {
      for( std::vector< particlePointer >::const_iterator
           ewIterator( shortcut->getEwsbBosonsAndMassiveVectorBosons(
                                                                  )->begin() );
           shortcut->getEwsbBosonsAndMassiveVectorBosons()->end() > ewIterator;
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
               (*cascadeIterator)->getSxCascades()->end()
               > sxIterator;
               ++sxIterator )
          {
            svsxCascades.addNewAtEnd()->setProperties( initialScolored,
                                       electroweakCascadeSource->getCascadeSet(
                                                               initialScolored,
                                                                   *ewIterator,
                                    (*cascadeIterator)->getInitialScolored() ),
                                                       *sxIterator );
            addOpenCascade( svsxCascades.getBack() );
          }
          for( std::vector< sjgxFullCascade* >::const_iterator
               sjgxIterator( (*cascadeIterator)->getSjgxCascades()->begin() );
               (*cascadeIterator)->getSjgxCascades()->end()
               > sjgxIterator;
               ++sjgxIterator )
          {
            svsjgxCascades.addNewAtEnd()->setProperties( initialScolored,
                                       electroweakCascadeSource->getCascadeSet(
                                                               initialScolored,
                                                                   *ewIterator,
                                    (*cascadeIterator)->getInitialScolored() ),
                                                         *sjgxIterator );
            addOpenCascade( svsjgxCascades.getBack() );
          }
          for( std::vector< sjgjsxFullCascade* >::const_iterator
               sjgjsxIterator( (*cascadeIterator)->getSjgjsxCascades(
                                                                  )->begin() );
               (*cascadeIterator)->getSjgjsxCascades()->end()
               > sjgjsxIterator;
               ++sjgjsxIterator )
          {
            svsjgjsxCascades.addNewAtEnd()->setProperties( initialScolored,
                                       electroweakCascadeSource->getCascadeSet(
                                                               initialScolored,
                                                                   *ewIterator,
                                    (*cascadeIterator)->getInitialScolored() ),
                                                           *sjgjsxIterator );
            addOpenCascade( svsjgjsxCascades.getBack() );
          }
        }
      }  // end of loop over if svs is open
    }
  }



  gluinoFullCascadeSet::gluinoFullCascadeSet(
                                            inputHandler const* const shortcut,
                     electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource,
                           readiableCascadeSetList* const squarkCascadeSetList,
                                              double const beamEnergy ) :
    fullCascadeSet( shortcut,
                    shortcut->getGluino(),
                    electroweakCascadeSource,
                    squarkCascadeSetList,
                    beamEnergy )
  {
    // we have to set up gxCascades now:
    for( std::vector< particlePointer >::const_iterator
         ewinoIterator( shortcut->getElectroweakinos()->begin() );
         shortcut->getUnstableElectroweakinos()->end() > ewinoIterator;
         ++ewinoIterator )
    {
      gxCascades.addNewAtEnd()->setProperties( shortcut,
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
    gvsxCascades.clearEntries();
    gjsvsxCascades.clearEntries();
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
           (*cascadeIterator)->getSxCascades()->end()
           > sxIterator;
           ++sxIterator )
      {
        gjsxCascades.addNewAtEnd()->setProperties( *sxIterator );
        addOpenCascade( gjsxCascades.getBack() );
      }
      for( std::vector< svsxFullCascade* >::const_iterator
           svsxIterator( (*cascadeIterator)->getSvsxCascades()->begin() );
           (*cascadeIterator)->getSvsxCascades()->end()
           > svsxIterator;
           ++svsxIterator )
      {
        gjsvsxCascades.addNewAtEnd()->setProperties( *svsxIterator );
        addOpenCascade( gjsvsxCascades.getBack() );
      }

      // now we check for gvs:
      for( std::vector< particlePointer >::const_iterator
           ewIterator( shortcut->getEwsbBosonsAndMassiveVectorBosons(
                                                                  )->begin() );
           shortcut->getEwsbBosonsAndMassiveVectorBosons()->end() > ewIterator;
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
            gvsxCascades.addNewAtEnd()->setProperties(
                                       electroweakCascadeSource->getCascadeSet(
                                                               initialScolored,
                                                                 *ewIterator,
                                    (*cascadeIterator)->getInitialScolored() ),
                                                       *sxIterator );
            addOpenCascade( gvsxCascades.getBack() );
          }
        }
      }  // end of loop over if gvs is open
    }  // end of loop over squarkFullCascadeSets
  }



  fullCascadeSetsForOneBeamEnergy::fullCascadeSetsForOneBeamEnergy(
                                            inputHandler const* const shortcut,
           electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource,
                                                    double const beamEnergy ) :
    shortcut( shortcut ),
    electroweakCascadeSource( electroweakCascadeSource ),
    beamEnergy( beamEnergy ),
    squarkCascadeSetList( new std::list< fullCascadeSet* >,
                   new publicGetsReadiedForNewPoint( shortcut->getReadier() ) )
  {
    gluinoCascadeSet = new gluinoFullCascadeSet( shortcut,
                                                 electroweakCascadeSource,
                                                 &squarkCascadeSetList,
                                                 beamEnergy );
    for( std::vector< particlePointer >::const_iterator
         squarkIterator( shortcut->getSquarks()->begin() );
         shortcut->getSquarks()->end() > squarkIterator;
         ++squarkIterator )
    {
      squarkCascadeSets.push_back( new squarkFullCascadeSet( shortcut,
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
   * shortcut->getSquarks().
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
                                            inputHandler const* const shortcut,
        electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource ) :
    shortcut( shortcut ),
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
   * fullCascadeSetsForOneBeamEnergy with the requested acceptanceValues, & if not,
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
      returnPointer = new fullCascadeSetsForOneBeamEnergy( shortcut,
                                                      electroweakCascadeSource,
                                                           beamEnergy );
      cascadeSetsPerEnergy.push_back( returnPointer );
    }
    return returnPointer;
  }


}  // end of LHC_FASER namespace.
