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

#include "LHC-FASER_full_cascade_stuff.hpp"

namespace LHC_FASER
{
  int const fullCascade::maximumJetsFromEwCascade( 2 );
  int const fullCascade::maximumLeptonsOfEachTypeFromEwCascade( 1 );

  fullCascade::fullCascade( colorfulCascadeType const typeOfColorfulCascade,
                            int const firstDecayBodyNumber,
                            int const maximumSmFermionsFromElectroweakCascades,
                            inputHandler const* const inputShortcut,
                            particlePointer const initialSparticle,
                            double const beamEnergy ) :
    inputShortcut( inputShortcut ),
    initialSparticle( initialSparticle ),
    beamEnergy( beamEnergy ),
    typeOfColorfulCascade( typeOfColorfulCascade ),
    maximumSmFermionsFromElectroweakCascades(
                                     maximumSmFermionsFromElectroweakCascades ),
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
  fullCascade::specifiedJetsSpecifiedChargeSummedLeptons(
                                  acceptanceCutSet const* const acceptanceCuts,
                                              int const numberOfAdditionalJets,
                                                    int const numberOfLeptons )
  {
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "fullCascade::specifiedJetsSpecifiedChargeSummedLeptons( "
    << acceptanceCuts << ", " << numberOfAdditionalJets << ", "
    << numberOfLeptons << " ) called.";
    std::cout << std::endl;**/

    if( ( 0 > numberOfAdditionalJets )
        ||
        ( 0 > numberOfLeptons )
        ||
        ( maximumSmFermionsFromElectroweakCascades
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

      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "fullCascade::specifiedJetsSpecifiedChargeSummedLeptons( "
      << acceptanceCuts << ", " << numberOfAdditionalJets << ", "
      << numberOfLeptons << " ) returning " << returnDouble;
      std::cout << std::endl;**/

      return returnDouble;
    }  // end of if numberOfLeptons was in the allowed range.
  }

  double
  fullCascade::unspecifiedJetsSpecifiedChargeSummedLeptons(
                                  acceptanceCutSet const* const acceptanceCuts,
                                                    int const numberOfLeptons )
  {
    if( ( 0 > numberOfLeptons )
        ||
        ( maximumSmFermionsFromElectroweakCascades
          < numberOfLeptons ) )
    {
      return 0.0;
    }
    else
    {
      double returnDouble( 0.0 );
      for( int numberOfJets( ( maximumSmFermionsFromElectroweakCascades
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
                              electroweakCascadeSet* const directEwinoCascades,
                         int const maximumSmFermionsFromElectroweakCascades ) :
        fullCascade( squarkToEwino,
                     2,
                     maximumSmFermionsFromElectroweakCascades,
                     inputShortcut,
                     initialSquark,
                     beamEnergy ),
        getsReadiedForNewPoint( inputShortcut->getReadier() ),
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
                                       directEwinoCascades,
                                       2 ),
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
                                acceptanceCutSet const* const acceptanceCuts,
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
        return directEwinoCascades->getAcceptance( acceptanceCuts,
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
                                       directEwinoCascades,
                                       4 ),
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
                      bosonCascades->getElectroweakDecayer()->get_PDG_code() );
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
        if( ( initialSparticle->get_absolute_mass() > ewinoMass )
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
                              acceptanceCutSet const* const acceptanceCuts,
                              int const numberOfAdditionalJets,
                              int numberOfNegativeElectrons,
                              int numberOfPositiveElectrons,
                              int numberOfNegativeMuons,
                              int numberOfPositiveMuons )
      /* this calls the appropriate functions on directEwinoCascades to build
       * the required acceptance or ewinoWithWCascades & bosonCascades if the
       * decay with a boson is not negligible, taking into account whether the
       * charges should be swapped if scoloredIsNotAntiparticle is false.
       */
      {
        // debugging:
        /**std::cout << std::endl << "debugging:"
        << std::endl
        << "fullCascadeType::squarkDirectlyToElectroweakType::supType::"
        << "getAcceptance( " << initialSparticleIsNotAntiparticle
        << ", " << acceptanceCuts << ", " << numberOfAdditionalJets << ", "
        << numberOfNegativeElectrons << ", " << numberOfPositiveElectrons
        << ", " << numberOfNegativeMuons << ", " << numberOfPositiveMuons
        << " ) called. shouldUseDecaysWithW = " << shouldUseDecaysWithW;
        std::cout << std::endl;**/

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
      supType::specifiedJetsOneOssfMinusOsdfPair(
                                  acceptanceCutSet const* const acceptanceCuts,
                                             int const numberOfAdditionalJets )
      // this should add up all combinations of jets while looking only for
      // one OSSF-OSDF pair.
      {
        /* if there is no possibility of extra jets from a boson decay, we can
         * only get a single OSSF-OSDF pair from an electroweakino decay, which
         * cannot produce additional jets.
         */
        if( shouldUseDecaysWithW )
        {
          double
          returnValue( wFraction
                       * bosonCascades->getAcceptance( acceptanceCuts,
                                                       numberOfAdditionalJets,
                                                       0,
                                                       0,
                                                       0,
                                                       0 ) );
          if( 0 == numberOfAdditionalJets )
          {
            returnValue += directFraction;
          }
          return ( returnValue
                   * directEwinoCascades->getOssfMinusOsdf( acceptanceCuts ) );

        }
        else if( 0 == numberOfAdditionalJets )
        {
          return directEwinoCascades->getOssfMinusOsdf( acceptanceCuts );
        }
        else
        {
          return 0.0;
        }
      }

      double
      supType::getCombinedAcceptance(
                                  acceptanceCutSet const* const acceptanceCuts,
                                      int const numberOfAdditionalJets,
                                      int const numberOfNegativeElectrons,
                                      int const numberOfPositiveElectrons,
                                      int const numberOfNegativeMuons,
                                      int const numberOfPositiveMuons )
      {
        // debugging:
        /**std::cout << std::endl << "debugging:"
        << std::endl
        << "fullCascadeType::squarkDirectlyToElectroweakType::supType::"
        << "getCombinedAcceptance( " << acceptanceCuts << ", "
        << numberOfAdditionalJets << ", " << numberOfNegativeElectrons << ", "
        << numberOfPositiveElectrons << ", " << numberOfNegativeMuons << ", "
        << numberOfPositiveMuons << " ) called."
        << " maximumSmFermionsFromElectroweakCascades = "
        << maximumSmFermionsFromElectroweakCascades;
        std::cout << std::endl;**/
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
            ( maximumSmFermionsFromElectroweakCascades
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
                    // debugging:
                    /**std::cout << std::endl << "debugging:"
                    << std::endl
                    << "ewinoWithWCascades->getAcceptance( "
                    << acceptanceCuts << ", " << ewinoJets << ", "
                    << ewinoNegativeElectrons << ", " << ewinoPositiveElectrons
                    << ", " << ewinoNegativeMuons << ", " << ewinoPositiveMuons
                    << " ) = "
                    << ewinoWithWCascades->getAcceptance( acceptanceCuts,
                                                          ewinoJets,
                                                        ewinoNegativeElectrons,
                                                        ewinoPositiveElectrons,
                                                          ewinoNegativeMuons,
                                                         ewinoPositiveMuons )
                    << std::endl
                    << "bosonCascades->getAcceptance( "
                    << acceptanceCuts << ", "
                    << ( numberOfAdditionalJets - ewinoJets ) << ", "
                    << ( numberOfNegativeElectrons - ewinoNegativeElectrons )
                    << ", "
                    << ( numberOfPositiveElectrons - ewinoPositiveElectrons )
                    << ", " << ( numberOfNegativeMuons - ewinoNegativeMuons )
                    << ", " << ( numberOfPositiveMuons - ewinoPositiveMuons )
                    << " ) = "
                    << bosonCascades->getAcceptance( acceptanceCuts,
                                        ( numberOfAdditionalJets - ewinoJets ),
                        ( numberOfNegativeElectrons - ewinoNegativeElectrons ),
                        ( numberOfPositiveElectrons - ewinoPositiveElectrons ),
                                ( numberOfNegativeMuons - ewinoNegativeMuons ),
                              ( numberOfPositiveMuons - ewinoPositiveMuons ) );
                    std::cout << std::endl;**/

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



    gluinoDirectlyToElectroweak::gluinoDirectlyToElectroweak(
                                       inputHandler const* const inputShortcut,
                                                       double const beamEnergy,
                                 electroweakCascadeSet* const ewinoCascades ) :
        fullCascade( gluinoToEwino,
                     3,
                     2,
                     inputShortcut,
                     inputShortcut->getGluino(),
                     beamEnergy ),
        getsReadiedForNewPoint( inputShortcut->getReadier() ),
        ewinoCascades( ewinoCascades )
    {
      cascadeSegment = cascadeDefiner.addNewAtEnd();
      cascadeSegment->first = ewinoCascades->getElectroweakDecayer();
      // gluinoToEwino means only 1 decay to be recorded.
      cascadeSegment->second = firstDecayBodyNumber;
      // gluinoToEwino also means that the decay is 3-body.
      soughtDecayProductList.front() = cascadeSegment->first->get_PDG_code();
    }

    gluinoDirectlyToElectroweak::~gluinoDirectlyToElectroweak()
    {
      // does nothing.
    }


    double
    gluinoDirectlyToElectroweak::getAcceptance(
                                  bool const initialSparticleIsNotAntiparticle,
                                  acceptanceCutSet const* const acceptanceCuts,
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


    fullCascade*
    squarkByBosonToCompound::setProperties(
                                           particlePointer const initialSquark,
                                    electroweakCascadeSet* const bosonCascades,
                                         fullCascade* const subcascadePointer )
    {
      initialSparticle = initialSquark;
      buildOn( subcascadePointer );
      maximumSmFermionsFromElectroweakCascades
      = ( subcascadePointer->getMaximumSmFermionsFromElectroweakCascades()
          + 2 );
      resetCachedBranchingRatio();
      this->bosonCascades = bosonCascades;
      if( ( CppSLHA::PDG_code::W_plus
            == bosonCascades->getElectroweakDecayer()->get_PDG_code() )
          &&
          ( inputShortcut->isIn( initialSquark->get_PDG_code(),
                                 inputShortcut->getSdownTypes() ) ) )
        // if we have to worry about which sign of PDG code to use...
      {
        bosonFlipsCharge = true;
        soughtDecayProductList.back() = -(CppSLHA::PDG_code::W_plus);
      }
      else
      {
        bosonFlipsCharge = false;
        soughtDecayProductList.back()
        = bosonCascades->getElectroweakDecayer()->get_PDG_code();
      }
      return this;
    }

    double
    squarkByBosonToCompound::getAcceptance(
                                  bool const initialSparticleIsNotAntiparticle,
                                  acceptanceCutSet const* const acceptanceCuts,
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
          ( maximumSmFermionsFromElectroweakCascades
            < ( numberOfAdditionalJets
                + numberOfNegativeElectrons + numberOfPositiveElectrons
                + numberOfNegativeMuons + numberOfPositiveMuons ) ) )
      {
        return 0.0;
      }
      else
      {
        /* we have to rearrange the number of charges of lepton sought based on
         * bosonCascades assuming a positively-charged W boson.
         * sup decaying:
         * initialSparticleIsNotAntiparticle = true,
         * bosonFlipsCharge = false
         * -> charges left alone, no problem for subsequent sdown.
         * antisup decaying:
         * initialSparticleIsNotAntiparticle = false,
         * bosonFlipsCharge = false
         * -> charges flipped, no problem for subsequent sdown
         * (e.g. 2 positive muons sought from antisup & subsequent antisdown ->
         * 2 negative muons sought from W^+ & subsequent sdown).
         * sdown decaying:
         * initialSparticleIsNotAntiparticle = true,
         * bosonFlipsCharge = true
         * -> charges flipped, subsequent sup must flip charges again
         * (e.g. 2 positive muons sought from antisdown & subsequent antisup ->
         * 2 negative muons sought from W^+ & subsequent antisup ->
         * positive muons again sought from sup).
         * antisup decaying:
         * initialSparticleIsNotAntiparticle = false,
         * bosonFlipsCharge = true
         * -> charges left alone, subsequent sup must flip charges
         * (e.g. 2 positive muons sought from antisdown & subsequent antisup ->
         * 2 positive muons sought from W^+ & subsequent antisup -> negative
         * muons sought from sup).
         * (above is fine for sdown decaying by Z or spin-0, since those cases
         * have bosonFlipsCharge = false.)
         */
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
        for( int
            ewinoJets( lhcFaserGlobal::smallerInteger( numberOfAdditionalJets,
                                                  maximumJetsFromEwCascade ) );
             0 <= ewinoJets;
             --ewinoJets )
        {
          for( int ewinoNegativeElectrons( lhcFaserGlobal::smallerInteger(
                                                     numberOfNegativeElectrons,
                                     maximumLeptonsOfEachTypeFromEwCascade ) );
               0 <= ewinoNegativeElectrons;
               --ewinoNegativeElectrons )
          {
            for( int ewinoPositiveElectrons( lhcFaserGlobal::smallerInteger(
                                                     numberOfPositiveElectrons,
                                     maximumLeptonsOfEachTypeFromEwCascade ) );
                 0 <= ewinoPositiveElectrons;
                 --ewinoPositiveElectrons )
            {
              for( int ewinoNegativeMuons( lhcFaserGlobal::smallerInteger(
                                                         numberOfNegativeMuons,
                                     maximumLeptonsOfEachTypeFromEwCascade ) );
                   0 <= ewinoNegativeMuons;
                   --ewinoNegativeMuons )
              {
                for( int ewinoPositiveMuons( lhcFaserGlobal::smallerInteger(
                                                         numberOfPositiveMuons,
                                     maximumLeptonsOfEachTypeFromEwCascade ) );
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
                       * subcascadePointer->getAcceptance( !bosonFlipsCharge,
                                                           acceptanceCuts,
                                                           ewinoJets,
                                                        ewinoNegativeElectrons,
                                                        ewinoPositiveElectrons,
                                                           ewinoNegativeMuons,
                                                        ewinoPositiveMuons ) );
                  /* in the used conventions, an incoming squark decays to a
                   * boson plus a squark, hence the same bool is used for both.
                   * also, initialSparticleIsNotAntiparticle was used to swap
                   * the numbers of charged leptons sought, & now should be
                   * replaced by the opposite of bosonFlipsCharge, for reasons
                   * outlined above.
                   */
                }  // end of loop over positive electrons.
              }  // end of loop over negative muons.
            }  // end of loop over positive electrons.
          }  // end of loop over negative electrons.
        }  // end of loop over jets.
        return returnDouble;
      }  // end of if number of SM fermions was in the allowed range.
    }

    double
    squarkByBosonToCompound::unspecifiedJetsSpecifiedOssfMinusOsdfPairs(
                                  acceptanceCutSet const* const acceptanceCuts,
                                                int const numberOfLeptonPairs )
    // this should add up all appropriate OSSF-OSDF pair combinations.
    {
      if( 0 <= numberOfLeptonPairs )
      {
        return 1.0;
      }
      else
      {
        return
        ( ( bosonCascades->getOssfMinusOsdf( acceptanceCuts )
            * subcascadePointer->unspecifiedJetsSpecifiedOssfMinusOsdfPairs(
                                                                acceptanceCuts,
                                                ( numberOfLeptonPairs - 1 ) ) )
          + ( ( bosonCascades->getAcceptance( acceptanceCuts,
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
               * subcascadePointer->unspecifiedJetsSpecifiedOssfMinusOsdfPairs(
                                                                acceptanceCuts,
                                                     numberOfLeptonPairs ) ) );
      }
    }

    double
    squarkByBosonToCompound::specifiedJetsOneOssfMinusOsdfPair(
                                  acceptanceCutSet const* const acceptanceCuts,
                                             int const numberOfAdditionalJets )
    // this should add up all combinations of jets while looking only for
    // one OSSF-OSDF pair.
    {
      return
      ( ( bosonCascades->getOssfMinusOsdf( acceptanceCuts )
          * subcascadePointer->specifiedJetsSpecifiedChargeSummedLeptons(
                                                                acceptanceCuts,
                                                        numberOfAdditionalJets,
                                                                          0 ) )
          + ( bosonCascades->getAcceptance( acceptanceCuts,
                                            2,
                                            0,
                                            0,
                                            0,
                                            0 )
              * subcascadePointer->specifiedJetsOneOssfMinusOsdfPair(
                                                                acceptanceCuts,
                                             ( numberOfAdditionalJets - 2 ) ) )
              + ( bosonCascades->getAcceptance( acceptanceCuts,
                                                1,
                                                0,
                                                0,
                                                0,
                                                0 )
                  * subcascadePointer->specifiedJetsOneOssfMinusOsdfPair(
                                                                acceptanceCuts,
                                             ( numberOfAdditionalJets - 1 ) ) )
              + ( bosonCascades->getAcceptance( acceptanceCuts,
                                                0,
                                                0,
                                                0,
                                                0,
                                                0 )
                  * subcascadePointer->specifiedJetsOneOssfMinusOsdfPair(
                                                                acceptanceCuts,
                                                  numberOfAdditionalJets ) ) );
    }



    squarkByJetToCompound::squarkByJetToCompound() :
        fullCascade( squarkToGauginoThenMore,
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
      sdownType::sdownType() :
          squarkByJetToCompound()
      {
        // just an initialization list.
      }

      sdownType::~sdownType()
      {
        // does nothing.
      }



      supType::supType() :
          squarkByJetToCompound(),
          bosonCascades( NULL ),
          wBoson( NULL ),
          ewinoMass( CppSLHA::CppSLHA_global::really_wrong_value ),
          shouldUseDecaysWithW( false ),
          directFraction( CppSLHA::CppSLHA_global::really_wrong_value ),
          wFraction( CppSLHA::CppSLHA_global::really_wrong_value ),
          decayProductListIncludingW( 2,
                                      CppSLHA::PDG_code::W_plus )
      {
        // just an initialization list.
      }

      supType::~supType()
      {
        // does nothing.
      }


      double
      supType::getAcceptance( bool const initialSparticleIsNotAntiparticle,
                              acceptanceCutSet const* const acceptanceCuts,
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
        // debugging:
        /**std::cout << std::endl << "debugging:"
        << std::endl
        << "fullCascadeType::squarkByJetToCompound::supType::getAcceptance( "
        << initialSparticleIsNotAntiparticle << ", " << acceptanceCuts << ", "
        << numberOfAdditionalJets << ", " << numberOfNegativeElectrons << ", "
        << numberOfPositiveElectrons << ", " << numberOfNegativeMuons << ", "
        << numberOfPositiveMuons << " ) called. shouldUseDecaysWithW = "
        << shouldUseDecaysWithW;
        std::cout << std::endl;**/

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
          return subcascadePointer->getAcceptance( true,
                                                   acceptanceCuts,
                                                   numberOfAdditionalJets,
                                                   numberOfNegativeElectrons,
                                                   numberOfPositiveElectrons,
                                                   numberOfNegativeMuons,
                                                   numberOfPositiveMuons );
          // initialSparticleIsNotAntiparticle was used to swap the numbers of
          // charged leptons sought, & now should be replaced by true.
        }
      }

      double
      supType::specifiedJetsOneOssfMinusOsdfPair(
                                  acceptanceCutSet const* const acceptanceCuts,
                                             int const numberOfAdditionalJets )
      // this should add up all combinations of jets while looking only for
      // one OSSF-OSDF pair.
      {
        if( 0 <= numberOfAdditionalJets )
        {
          /* if there is no possibility of extra jets from a boson decay, we
           * can only get a single OSSF-OSDF pair from an electroweakino decay,
           * which cannot produce additional jets.
           */
          if( shouldUseDecaysWithW )
          {
            double
            wContribution( bosonCascades->getAcceptance( acceptanceCuts,
                                                         0,
                                                         0,
                                                         0,
                                                         0,
                                                         0 )
                        * subcascadePointer->specifiedJetsOneOssfMinusOsdfPair(
                                                                acceptanceCuts,
                                                    numberOfAdditionalJets ) );
            // first the contribution of no jets from the boson is worked out,
            // then any other contributions.
            for( int jetCounter( lhcFaserGlobal::smallerInteger(
                                                        numberOfAdditionalJets,
                                                  maximumJetsFromEwCascade ) );
                 0 < jetCounter;
                 --jetCounter )
            {
              wContribution
              += ( bosonCascades->getAcceptance( acceptanceCuts,
                                                 jetCounter,
                                                 0,
                                                 0,
                                                 0,
                                                 0 )
                   * subcascadePointer->specifiedJetsOneOssfMinusOsdfPair(
                                                                acceptanceCuts,
                                   ( numberOfAdditionalJets - jetCounter ) ) );
            }
            return ( ( wFraction * wContribution )
                     + ( directFraction
                        * subcascadePointer->specifiedJetsOneOssfMinusOsdfPair(
                                                                acceptanceCuts,
                                                  numberOfAdditionalJets ) ) );
          }
          else
          {
            return subcascadePointer->specifiedJetsOneOssfMinusOsdfPair(
                                                                acceptanceCuts,
                                                      numberOfAdditionalJets );
          }
        }
        else
        {
          return 0.0;
        }
      }

      bool
      supType::decayWithWIsNotNegligible()
      // this returns true if the decay involving a W boson is not negligible,
      // also setting up the relevant fractions, false otherwise.
      {
        // debugging:
        /**std::cout << std::endl << "debugging:"
        << std::endl
        << "fullCascadeType::squarkByJetToCompound::supType::"
        << "decayWithWIsNotNegligible() called. decayProductListIncludingW = ";
        for( std::list< int >::iterator
             codeIterator( decayProductListIncludingW.begin() );
             decayProductListIncludingW.end() != codeIterator;
             ++codeIterator )
        {
          std::cout << " " << *codeIterator;
        }
        std::cout
        << std::endl
        << "soughtDecayProductList = ";
        for( std::list< int >::iterator
             codeIterator( soughtDecayProductList.begin() );
             soughtDecayProductList.end() != codeIterator;
             ++codeIterator )
        {
          std::cout << " " << *codeIterator;
        }
        std::cout << std::endl;
        std::cout
        << std::endl
        << "numerator = "
        << initialSparticle->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                                  &decayProductListIncludingW )
        << std::endl
        << "denominator = "
        << initialSparticle->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                                     &soughtDecayProductList );
        std::cout << std::endl;**/

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
      supType::getCombinedAcceptance(
                                  acceptanceCutSet const* const acceptanceCuts,
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
            ( maximumSmFermionsFromElectroweakCascades
              < ( numberOfAdditionalJets
                  + numberOfNegativeElectrons + numberOfPositiveElectrons
                  + numberOfNegativeMuons + numberOfPositiveMuons ) ) )
        {
          return 0.0;
        }
        else
        {
          double wContribution( 0.0 );
          for( int ewinoJets( lhcFaserGlobal::smallerInteger(
                                                        numberOfAdditionalJets,
                                                  maximumJetsFromEwCascade ) );
               0 <= ewinoJets;
               --ewinoJets )
          {
            for( int ewinoNegativeElectrons( lhcFaserGlobal::smallerInteger(
                                                     numberOfNegativeElectrons,
                                     maximumLeptonsOfEachTypeFromEwCascade ) );
                 0 <= ewinoNegativeElectrons;
                 --ewinoNegativeElectrons )
            {
              for( int ewinoPositiveElectrons( lhcFaserGlobal::smallerInteger(
                                                     numberOfPositiveElectrons,
                                     maximumLeptonsOfEachTypeFromEwCascade ) );
                   0 <= ewinoPositiveElectrons;
                   --ewinoPositiveElectrons )
              {
                for( int ewinoNegativeMuons( lhcFaserGlobal::smallerInteger(
                                                         numberOfNegativeMuons,
                                     maximumLeptonsOfEachTypeFromEwCascade ) );
                     0 <= ewinoNegativeMuons;
                     --ewinoNegativeMuons )
                {
                  for( int ewinoPositiveMuons( lhcFaserGlobal::smallerInteger(
                                                         numberOfPositiveMuons,
                                     maximumLeptonsOfEachTypeFromEwCascade ) );
                       0 <= ewinoPositiveMuons;
                       --ewinoPositiveMuons )
                  {
                    wContribution
                    += ( bosonCascades->getAcceptance( acceptanceCuts,
                                                       ewinoJets,
                                                       ewinoNegativeElectrons,
                                                       ewinoPositiveElectrons,
                                                       ewinoNegativeMuons,
                                                       ewinoPositiveMuons )
                         * subcascadePointer->getAcceptance( true,
                                                             acceptanceCuts,
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
          return ( ( wFraction * wContribution )
                   + ( directFraction
                       * subcascadePointer->getAcceptance( true,
                                                           acceptanceCuts,
                                                        numberOfAdditionalJets,
                                                     numberOfNegativeElectrons,
                                                     numberOfPositiveElectrons,
                                                         numberOfNegativeMuons,
                                                   numberOfPositiveMuons ) ) );
        }  // end of if number of SM fermions was in the allowed range.
      }

    }  // end of squarkByJetToCompoundType namespace



    gluinoOrElectroweakinoToCompound::gluinoOrElectroweakinoToCompound() :
        fullCascade( gluinoOrElectroweakinoToSquarkThenMore,
                     2 ),
        bosonCascades( NULL ),
        wBoson( NULL ),
        ewinoMass( CppSLHA::CppSLHA_global::really_wrong_value ),
        shouldUseDecaysWithW( false ),
        directFraction( CppSLHA::CppSLHA_global::really_wrong_value ),
        wFraction( CppSLHA::CppSLHA_global::really_wrong_value ),
        decayProductListIncludingW( 2,
                                    CppSLHA::PDG_code::W_plus )
        // gjm also means that the initial decay is 2-body.
    {
      // just an initialization list.
    }

    gluinoOrElectroweakinoToCompound::~gluinoOrElectroweakinoToCompound()
    {
      // does nothing.
    }


    double
    gluinoOrElectroweakinoToCompound::specifiedJetsOneOssfMinusOsdfPair(
                                  acceptanceCutSet const* const acceptanceCuts,
                                             int const numberOfAdditionalJets )
    // this should add up all combinations of jets while looking only for
    // one OSSF-OSDF pair.
    {
      if( 0 <= numberOfAdditionalJets )
      {
        /* if there is no possibility of extra jets from a boson decay, we
         * can only get a single OSSF-OSDF pair from an electroweakino decay,
         * which cannot produce additional jets.
         */
        if( shouldUseDecaysWithW )
        {
          double
          wContribution( bosonCascades->getAcceptance( acceptanceCuts,
                                                       0,
                                                       0,
                                                       0,
                                                       0,
                                                       0 )
                        * subcascadePointer->specifiedJetsOneOssfMinusOsdfPair(
                                                                acceptanceCuts,
                                                    numberOfAdditionalJets ) );
          // first the contribution of no jets from the boson is worked out,
          // then any contributions
          for( int jetCounter( lhcFaserGlobal::smallerInteger(
                                                        numberOfAdditionalJets,
                                                  maximumJetsFromEwCascade ) );
               0 < jetCounter;
               --jetCounter )
          {
            wContribution
            += ( bosonCascades->getAcceptance( acceptanceCuts,
                                               jetCounter,
                                               0,
                                               0,
                                               0,
                                               0 )
                 * subcascadePointer->specifiedJetsOneOssfMinusOsdfPair(
                                                                acceptanceCuts,
                                   ( numberOfAdditionalJets - jetCounter ) ) );
          }
          return ( ( wFraction * wContribution )
                   + ( directFraction
                       * subcascadePointer->specifiedJetsOneOssfMinusOsdfPair(
                                                                acceptanceCuts,
                                                  numberOfAdditionalJets ) ) );
        }
        else
        {
          return subcascadePointer->specifiedJetsOneOssfMinusOsdfPair(
                                                                acceptanceCuts,
                                                      numberOfAdditionalJets );
        }
      }
      else
      {
        return 0.0;
      }
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


      fullCascade*
      gluinoOrNeutralinoSet::setProperties(
                                        particlePointer const initialSparticle,
                                          fullCascade* const subcascadePointer,
          electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource )
      {
        // debugging:
        /**std::cout << std::endl << "debugging:"
        << std::endl
        << "gluinoOrNeutralinoSet::setProperties( "
        << *(initialSparticle->get_name()) << ", " << subcascadePointer;
        if( NULL != subcascadePointer )
        {
          std::cout << " = [ " << subcascadePointer->getAsString() << " ]";
        }
        std::cout << ", " << electroweakCascadeSource << " ) called.";
        std::cout << std::endl;**/

        this->initialSparticle = initialSparticle;
        this->electroweakCascadeSource = electroweakCascadeSource;
        buildOn( subcascadePointer );
        maximumSmFermionsFromElectroweakCascades
        = ( subcascadePointer->getMaximumSmFermionsFromElectroweakCascades()
            + 2 );
        resetCachedBranchingRatio();
        decayProductListIncludingW.front()
        = subcascadePointer->getInitialSparticle()->get_PDG_code();
        // only sup-types are going to be produced with virtual antitops by
        // gluinos or neutralinos:
        if( inputShortcut->isIn( decayProductListIncludingW.front(),
                                 inputShortcut->getSupTypes() ) )
        {
          shouldUseDecaysWithW = decayWithWIsNotNegligible();
        }
        else
        {
          shouldUseDecaysWithW = false;
        }
        return this;
      }

      double
      gluinoOrNeutralinoSet::getAcceptance(
                                  bool const initialSparticleIsNotAntiparticle,
                                  acceptanceCutSet const* const acceptanceCuts,
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
          // debugging:
          /**std::cout << std::endl << "debugging:"
          << std::endl
          << "gluinoOrNeutralinoSet::getAcceptance( "
          << initialSparticleIsNotAntiparticle << ", " << acceptanceCuts
          << ", " << numberOfAdditionalJets << ", "
          << numberOfNegativeElectrons << ", " << numberOfPositiveElectrons
          << ", " << numberOfNegativeMuons << ", " << numberOfPositiveMuons
          << " ) has subcascadePointer->getAcceptance( "
          << initialSparticleIsNotAntiparticle << " ... ) = "
          << subcascadePointer->getAcceptance(
                                            initialSparticleIsNotAntiparticle,
                                               acceptanceCuts,
                                               numberOfAdditionalJets,
                                               numberOfNegativeElectrons,
                                               numberOfPositiveElectrons,
                                               numberOfNegativeMuons,
                                               numberOfPositiveMuons );
          std::cout << std::endl;**/

          return
          ( 0.5 * ( subcascadePointer->getAcceptance( true,
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
        }
      }

      bool
      gluinoOrNeutralinoSet::decayWithWIsNotNegligible()
      // this returns true if the decay involving a W boson is not negligible,
      // also setting up the relevant fractions, false otherwise.
      {
        // debugging:
        /**std::cout << std::endl << "debugging:"
        << std::endl
        << "fullCascadeType::gluinoOrNeutralinoSet::decayWithWIsNotNegligible("
        << ") called. decayProductListIncludingW = ";
        std::cout << std::endl << "wtf?" << std::endl;
        std::cout << std::endl;
        for( std::list< int >::iterator
             codeIterator( decayProductListIncludingW.begin() );
             decayProductListIncludingW.end() != codeIterator;
             ++codeIterator )
        {
          std::cout << " " << *codeIterator;
        }
        std::cout
        << std::endl
        << "soughtDecayProductList = ";
        for( std::list< int >::iterator
             codeIterator( soughtDecayProductList.begin() );
             soughtDecayProductList.end() != codeIterator;
             ++codeIterator )
        {
          std::cout << " " << *codeIterator;
        }
        std::cout << std::endl;
        std::cout
        << std::endl
        << "numerator = "
        << initialSparticle->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                                  &decayProductListIncludingW )
        << std::endl
        << "denominator = "
        << initialSparticle->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                                     &soughtDecayProductList );
        std::cout << std::endl;**/

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
                                  acceptanceCutSet const* const acceptanceCuts,
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
            ( maximumSmFermionsFromElectroweakCascades
              < ( numberOfAdditionalJets
                  + numberOfNegativeElectrons + numberOfPositiveElectrons
                  + numberOfNegativeMuons + numberOfPositiveMuons ) ) )
        {
          return 0.0;
        }
        else
        {
          double returnDouble( 0.0 );
          for( int ewinoJets( lhcFaserGlobal::smallerInteger(
                                                        numberOfAdditionalJets,
                                                  maximumJetsFromEwCascade ) );
               0 <= ewinoJets;
               --ewinoJets )
          {
            for( int ewinoNegativeElectrons( lhcFaserGlobal::smallerInteger(
                                                     numberOfNegativeElectrons,
                                     maximumLeptonsOfEachTypeFromEwCascade ) );
                 0 <= ewinoNegativeElectrons;
                 --ewinoNegativeElectrons )
            {
              for( int ewinoPositiveElectrons( lhcFaserGlobal::smallerInteger(
                                                     numberOfPositiveElectrons,
                                     maximumLeptonsOfEachTypeFromEwCascade ) );
                   0 <= ewinoPositiveElectrons;
                   --ewinoPositiveElectrons )
              {
                for( int ewinoNegativeMuons( lhcFaserGlobal::smallerInteger(
                                                         numberOfNegativeMuons,
                                     maximumLeptonsOfEachTypeFromEwCascade ) );
                     0 <= ewinoNegativeMuons;
                     --ewinoNegativeMuons )
                {
                  for( int ewinoPositiveMuons( lhcFaserGlobal::smallerInteger(
                                                         numberOfPositiveMuons,
                                     maximumLeptonsOfEachTypeFromEwCascade ) );
                       0 <= ewinoPositiveMuons;
                       --ewinoPositiveMuons )
                  {
                    /* we average over the configurations with opposite charges
                     * (could use false for initialSparticleIsNotAntiparticle
                     * for subcascadePointer & appropriate negative/positive
                     * numbers, but this way seems simpler):
                     */
                    returnDouble
                    += ( 0.5 * ( bosonCascades->getAcceptance( acceptanceCuts,
                                                               ewinoJets,
                                                        ewinoNegativeElectrons,
                                                        ewinoPositiveElectrons,
                                                            ewinoNegativeMuons,
                                                           ewinoPositiveMuons )
                                 * subcascadePointer->getAcceptance( true,
                                                                acceptanceCuts,
                                        ( numberOfAdditionalJets - ewinoJets ),
                        ( numberOfNegativeElectrons - ewinoNegativeElectrons ),
                        ( numberOfPositiveElectrons - ewinoPositiveElectrons ),
                                ( numberOfNegativeMuons - ewinoNegativeMuons ),
                               ( numberOfPositiveMuons - ewinoPositiveMuons ) )
                                 + bosonCascades->getAcceptance( acceptanceCuts,
                                                                     ewinoJets,
                                                        ewinoPositiveElectrons,
                                                        ewinoNegativeElectrons,
                                                            ewinoPositiveMuons,
                                                           ewinoNegativeMuons )
                                   * subcascadePointer->getAcceptance( true,
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
          return ( ( wFraction * returnDouble )
                   + ( directFraction * 0.5
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
                                                 numberOfNegativeMuons ) ) ) );
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
        maximumSmFermionsFromElectroweakCascades
        = ( subcascadePointer->getMaximumSmFermionsFromElectroweakCascades()
            + 2 );
        resetCachedBranchingRatio();
        decayProductListIncludingW.front()
        = subcascadePointer->getInitialSparticle()->get_PDG_code();
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
                                  acceptanceCutSet const* const acceptanceCuts,
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
                                  acceptanceCutSet const* const acceptanceCuts,
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
            ( maximumSmFermionsFromElectroweakCascades
              < ( numberOfAdditionalJets
                  + numberOfNegativeElectrons + numberOfPositiveElectrons
                  + numberOfNegativeMuons + numberOfPositiveMuons ) ) )
        {
          return 0.0;
        }
        else
        {
          double returnDouble( 0.0 );
          for( int ewinoJets( lhcFaserGlobal::smallerInteger(
                                                        numberOfAdditionalJets,
                                                  maximumJetsFromEwCascade ) );
               0 <= ewinoJets;
               --ewinoJets )
          {
            for( int ewinoNegativeElectrons( lhcFaserGlobal::smallerInteger(
                                                     numberOfNegativeElectrons,
                                     maximumLeptonsOfEachTypeFromEwCascade ) );
                 0 <= ewinoNegativeElectrons;
                 --ewinoNegativeElectrons )
            {
              for( int ewinoPositiveElectrons( lhcFaserGlobal::smallerInteger(
                                                     numberOfPositiveElectrons,
                                     maximumLeptonsOfEachTypeFromEwCascade ) );
                   0 <= ewinoPositiveElectrons;
                   --ewinoPositiveElectrons )
              {
                for( int ewinoNegativeMuons( lhcFaserGlobal::smallerInteger(
                                                         numberOfNegativeMuons,
                                     maximumLeptonsOfEachTypeFromEwCascade ) );
                     0 <= ewinoNegativeMuons;
                     --ewinoNegativeMuons )
                {
                  for( int ewinoPositiveMuons( lhcFaserGlobal::smallerInteger(
                                                         numberOfPositiveMuons,
                                     maximumLeptonsOfEachTypeFromEwCascade ) );
                       0 <= ewinoPositiveMuons;
                       --ewinoPositiveMuons )
                  {
                    returnDouble
                    += ( ( squarkWithWFraction
                           * subcascadePointer->getAcceptance( false,
                                                               acceptanceCuts,
                                                               ewinoJets,
                        ( numberOfNegativeElectrons - ewinoNegativeElectrons ),
                        ( numberOfPositiveElectrons - ewinoPositiveElectrons ),
                                ( numberOfNegativeMuons - ewinoNegativeMuons ),
                               ( numberOfPositiveMuons - ewinoPositiveMuons ) )
                           + antisquarkWithWFraction
                             * subcascadePointer->getAcceptance( false,
                                                                acceptanceCuts,
                                                                 ewinoJets,
                        ( numberOfNegativeElectrons - ewinoNegativeElectrons ),
                        ( numberOfPositiveElectrons - ewinoPositiveElectrons ),
                                ( numberOfNegativeMuons - ewinoNegativeMuons ),
                             ( numberOfPositiveMuons - ewinoPositiveMuons ) ) )
                         * bosonCascades->getAcceptance( acceptanceCuts,
                                        ( numberOfAdditionalJets - ewinoJets ),
                                                        ewinoNegativeElectrons,
                                                        ewinoPositiveElectrons,
                                                        ewinoNegativeMuons,
                                                        ewinoPositiveMuons ) );
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
                                  double const beamEnergy ) :
    getsReadiedForNewPoint( inputShortcut->getReadier() ),
    inputShortcut( inputShortcut ),
    initialSparticle( initialScolored ),
    electroweakCascadeSource( electroweakCascadeSource ),
    openCascades(),
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
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "fullCascadeSetOrderer::orderCascades() called."
    << " sdownTypeColoredCascades.size() = "
    << sdownTypeColoredCascades.size()
    << ", supTypeColoredCascades.size() = "
    << supTypeColoredCascades.size()
    << ", allNeutralinoColoredCascades.size() = "
    << allNeutralinoColoredCascades.size()
    << ", relevantNeutralinoColoredCascades.size() = "
    << relevantNeutralinoColoredCascades.size()
    << ", allCharginoColoredCascades.size() = "
    << allCharginoColoredCascades.size()
    << ", relevantCharginoColoredCascades.size() = "
    << relevantCharginoColoredCascades.size();
    std::cout << std::endl;**/

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
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "fullCascadeSetOrderer::orderCascades() finished. lightestSquarkMass = "
    << lightestSquarkMass << std::endl << "sdownTypeColoredCascades =";
    for( std::list< fullCascadeSet* >::iterator
         setIterator( sdownTypeColoredCascades.begin() );
         sdownTypeColoredCascades.end() != setIterator;
         ++setIterator )
    {
      std::cout << " " << *((*setIterator)->getInitialSparticle()->get_name());
    }
    std::cout
    << std::endl << "supTypeColoredCascades =";
    for( std::list< fullCascadeSet* >::iterator
         setIterator( supTypeColoredCascades.begin() );
         supTypeColoredCascades.end() != setIterator;
         ++setIterator )
    {
      std::cout << " " << *((*setIterator)->getInitialSparticle()->get_name());
    }
    std::cout
    << std::endl << "allNeutralinoColoredCascades =";
    for( std::list< fullCascadeSet* >::iterator
         setIterator( allNeutralinoColoredCascades.begin() );
         allNeutralinoColoredCascades.end() != setIterator;
         ++setIterator )
    {
      std::cout << " " << *((*setIterator)->getInitialSparticle()->get_name());
    }
    std::cout
    << std::endl << "relevantNeutralinoColoredCascades =";
    for( std::list< fullCascadeSet* >::iterator
         setIterator( relevantNeutralinoColoredCascades.begin() );
         relevantNeutralinoColoredCascades.end() != setIterator;
         ++setIterator )
    {
      std::cout << " " << *((*setIterator)->getInitialSparticle()->get_name());
    }
    std::cout
    << std::endl << "allCharginoColoredCascades =";
    for( std::list< fullCascadeSet* >::iterator
         setIterator( allCharginoColoredCascades.begin() );
         allCharginoColoredCascades.end() != setIterator;
         ++setIterator )
    {
      std::cout << " " << *((*setIterator)->getInitialSparticle()->get_name());
    }
    std::cout
    << std::endl << "relevantCharginoColoredCascades =";
    for( std::list< fullCascadeSet* >::iterator
         setIterator( relevantCharginoColoredCascades.begin() );
         relevantCharginoColoredCascades.end() != setIterator;
         ++setIterator )
    {
      std::cout << " " << *((*setIterator)->getInitialSparticle()->get_name());
    }
    std::cout << std::endl;**/
  }



  namespace fullCascadeSetType
  {
    squarkSet::squarkSet( inputHandler const* const inputShortcut,
                          particlePointer const initialScolored,
           electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource,
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
        ewinoCodeIsAlwaysPositive( ewinoCodeIsAlwaysPositive ),
        bosonCodeIsAlwaysPositive( bosonCodeIsAlwaysPositive )
    {
      // just an initialization list.
    }

    squarkSet::~squarkSet()
    {
      for( std::vector< squarkDirectlyToElectroweakCascade* >::iterator
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
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "squarkSet::setUpCascades() for " << *(initialSparticle->get_name())
      << " called.";
      std::cout << std::endl;**/

      // 1st we add the direct cascades:
      for( std::vector< squarkDirectlyToElectroweakCascade* >::iterator
           cascadeIterator( directToEwinoCascades.begin() );
           directToEwinoCascades.end() > cascadeIterator;
           ++cascadeIterator )
      {
        // debugging:
        /**std::cout << std::endl << "debugging:"
        << std::endl
        << "squarkSet::setUpCascades() about to ask if "
        << (*cascadeIterator)->getAsString() << " is open";
        std::cout << std::endl;**/

        if( (*cascadeIterator)->isOpen() )
        {
          // debugging:
          /**std::cout << std::endl << "debugging:"
          << std::endl
          << "it was open.";
          std::cout << std::endl;**/

          openCascades.push_back( *cascadeIterator );
        }
      }

      // now we clear the compound cascades:
      compoundByBosonCascades.clearEntries();
      clearCompoundByJetCascades();

      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "squarkSet::setUpCascades() about to ask for gluino open cascades.";
      std::cout << std::endl;**/

      // next we check to see if we should add compound cascades from the
      // gluino:
      if( inputShortcut->getGluinoMass() < initialSparticleMass )
      {
        singleSpecifiedDecayProductList.front() = CppSLHA::PDG_code::gluino;
        potentialSubcascadeSet = gluinoFullCascade;
        addByJetCascade();
      }

      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "squarkSet::setUpCascades() about to ask for neutralino open"
      << " cascades.";
      std::cout << std::endl;**/

      // now we look at compound cascades from neutralinos:
      orderedCascadeSets = setOrderer->getNeutralinoColoredCascades();
      addEwinoCompoundCascade( true );

      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "squarkSet::setUpCascades() about to ask for chargino open"
      << " cascades.";
      std::cout << std::endl;**/

      // now we look at compound cascades from charginos:
      orderedCascadeSets = setOrderer->getCharginoColoredCascades();
      addEwinoCompoundCascade( ewinoCodeIsAlwaysPositive );
      // down-type squarks decay to negatively-charged charginos, up-type to
      // positively-charged.

      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "squarkSet::setUpCascades() about to ask for squark open cascades.";
      std::cout << std::endl;**/

      // now we look at compound cascades from squarks:
      buildSquarkCompoundCascades();
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
                  * (*cascadeIterator)->getTotalBrToEwino() ) )
          {
            // we add each cascade with an overall BR that is not negligible:
            openCascades.push_back( getNewCompoundByJetCascade(
                                                          *cascadeIterator ) );
          }
        }
      }
    }

    void
    squarkSet::addSquarkCompoundCascade( bool positiveBosonCode )
    // this relies on orderedCascadeSets & potentialDecayProducts having
    // already been set correctly.
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "squarkSet::addSquarkCompoundCascade( " << positiveBosonCode << " )";
      std::cout << std::endl;**/

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
                      * (*cascadeIterator)->getTotalBrToEwino() ) )
              {
                // we add each cascade with an overall BR that is not
                // negligible:
                openCascades.push_back(
                          compoundByBosonCascades.addNewAtEnd()->setProperties(
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
                            particlePointer const initialScolored,
           electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource,
                            fullCascadeSetOrderer* const setOrderer,
                            fullCascadeSet* const gluinoFullCascade,
                            double const beamEnergy ) :
          squarkSet( inputShortcut,
                     initialScolored,
                     electroweakCascadeSource,
                     setOrderer,
                     gluinoFullCascade,
                     beamEnergy,
                     false,
                     false ),
          compoundByJetCascades()
      {
        // we have to set up the cascades directly to electroweakinos now:
        for( std::vector< particlePointer >::const_iterator
             ewinoIterator( inputShortcut->getElectroweakinos()->begin() );
             inputShortcut->getElectroweakinos()->end() > ewinoIterator;
             ++ewinoIterator )
        {
          directToEwinoCascades.push_back(
               new fullCascadeType::squarkDirectlyToElectroweakType::sdownType(
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
                        particlePointer const initialScolored,
           electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource,
                        fullCascadeSetOrderer* const setOrderer,
                        fullCascadeSet* const gluinoFullCascade,
                        double const beamEnergy ) :
          squarkSet( inputShortcut,
                     initialScolored,
                     electroweakCascadeSource,
                     setOrderer,
                     gluinoFullCascade,
                     beamEnergy ),
          compoundByJetCascades(),
          twoSpecifiedDecayProductsList( 2,
                                 CppSLHA::CppSLHA_global::really_wrong_value ),
          appropriateSdownForWDecay( NULL )
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
          directToEwinoCascades.push_back(
                 new fullCascadeType::squarkDirectlyToElectroweakType::supType(
                                                                 inputShortcut,
                                                               initialScolored,
                                                                    beamEnergy,
                                       electroweakCascadeSource->getCascadeSet(
                                                               initialScolored,
                                                              *ewinoIterator ),
                                       electroweakCascadeSource->getCascadeSet(
                                                     appropriateSdownForWDecay,
                                                                *ewinoIterator,
                               inputShortcut->getSquarkMinusBosonEffectiveMass(
                                                               initialScolored,
                                                     inputShortcut->getWPlus(),
                                                            *ewinoIterator ) ),
                                       electroweakCascadeSource->getCascadeSet(
                                                               initialScolored,
                                                                *ewinoIterator,
                                inputShortcut->getSquarkPlusBosonEffectiveMass(
                                                     appropriateSdownForWDecay,
                                                     inputShortcut->getWPlus(),
                                                        *ewinoIterator ) ) ) );
        }
      }

      supType::~supType()
      {
        // does nothing.
      }

    }  // end of squarkSetType namespace



    gluinoOrElectroweakinoSet::gluinoOrElectroweakinoSet(
                                       inputHandler const* const inputShortcut,
                                        particlePointer const initialSparticle,
           electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource,
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
    gluinoOrElectroweakinoSet::buildSquarkCompoundCascades()
    // this does the job of finding the right squark subcascades.
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "gluinoOrElectroweakinoSet::buildSquarkCompoundCascades() called.";
      std::cout << std::endl;**/

      setIterator = orderedCascadeSets->begin();
      while( ( setIterator != orderedCascadeSets->end() )
             &&
             ( initialSparticleMass
               > (*setIterator)->getInitialSparticle()->get_absolute_mass() ) )
        // we go through the list of squarks until we get to either the end of
        // the list or run out of squarks lighter than this particle:
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
            // debugging:
            /**std::cout << std::endl << "debugging:"
            << std::endl
            << "about to ask " << (*cascadeIterator)->getAsString()
            << " what its total BR is.";
            std::cout << std::endl;**/

            if( lhcFaserGlobal::negligibleBr
                < ( subcascadeBranchingRatio
                    * (*cascadeIterator)->getTotalBrToEwino() ) )
            {
              // debugging:
              /**std::cout << std::endl << "debugging:"
              << std::endl
              << "subcascadeBranchingRatio = " << subcascadeBranchingRatio
              << ", ->getTotalBrToEwino() = "
              << (*cascadeIterator)->getTotalBrToEwino();
              std::cout << std::endl;**/

              // we add each cascade with an overall BR that is not negligible:
              openCascades.push_back( getNewCompoundCascade(
                                                          *cascadeIterator ) );

              // debugging:
              /**std::cout << std::endl << "debugging:"
              << std::endl
              << "added compound cascade.";
              std::cout << std::endl;**/
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
          gluinoOrElectroweakinoSet( inputShortcut,
                                     inputShortcut->getGluino(),
                                     electroweakCascadeSource,
                                     setOrderer,
                                     beamEnergy ),
      directToEwinoCascades(),
      compoundCascades()
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
        for( std::vector< gluinoDirectlyToElectroweakCascade* >::iterator
             deletionIterator( directToEwinoCascades.begin() );
             directToEwinoCascades.end() > deletionIterator;
             ++deletionIterator )
        {
          delete *deletionIterator;
        }
      }

      void
      gluinoSet::setUpCascades()
      {
        // 1st we add the direct cascades:
        for( std::vector< gluinoDirectlyToElectroweakCascade* >::iterator
             cascadeIterator( directToEwinoCascades.begin() );
             directToEwinoCascades.end() > cascadeIterator;
             ++cascadeIterator )
        {
          // debugging:
          /**std::cout << std::endl << "debugging:"
          << std::endl
          << "gluinoOrElectroweakinoSetType::gluinoSet::setUpCascades() about"
          << " to ask if " << (*cascadeIterator)->getAsString() << " is open";
          std::cout << std::endl;**/

          if( (*cascadeIterator)->isOpen() )
          {
            // debugging:
            /**std::cout << std::endl << "debugging:"
            << std::endl
            << "it was open.";
            std::cout << std::endl;**/

            openCascades.push_back( *cascadeIterator );
          }
        }

        // then we clear the compound cascades:
        compoundCascades.clearEntries();

        // now we look at compound cascades from squarks:
        orderedCascadeSets = setOrderer->getSdownTypeCascades();
        buildSquarkCompoundCascades();
        // debugging:
        /**std::cout << std::endl << "debugging:"
        << std::endl
        << "finished building gluino -> sdown-types.";
        std::cout << std::endl;**/

        orderedCascadeSets = setOrderer->getSupTypeCascades();
        buildSquarkCompoundCascades();

        // debugging:
        /**std::cout << std::endl << "debugging:"
        << std::endl
        << "finished building gluino -> sup-types.";
        std::cout << std::endl;**/
      }



      neutralinoSet::neutralinoSet( inputHandler const* const inputShortcut,
                                    particlePointer const initialSparticle,
           electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource,
                                    fullCascadeSetOrderer* const setOrderer,
                                    double const beamEnergy ) :
          gluinoOrElectroweakinoSet( inputShortcut,
                                     initialSparticle,
                                     electroweakCascadeSource,
                                     setOrderer,
                                     beamEnergy ),
      compoundCascades()
      {
        // just an initialization list.
      }

      neutralinoSet::~neutralinoSet()
      {
        // does nothing.
      }



      charginoSet::charginoSet( inputHandler const* const inputShortcut,
                                particlePointer const initialSparticle,
           electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource,
                                fullCascadeSetOrderer* const setOrderer,
                                double const beamEnergy ) :
          gluinoOrElectroweakinoSet( inputShortcut,
                                     initialSparticle,
                                     electroweakCascadeSource,
                                     setOrderer,
                                     beamEnergy ),
      compoundCascades()
      {
        // just an initialization list.
      }

      charginoSet::~charginoSet()
      {
        // does nothing.
      }

    }  // end of gluinoOrElectroweakinoSetType namespace

  }  // end of fullCascadeSetType namespace



  fullCascadeSetsForOneBeamEnergy::fullCascadeSetsForOneBeamEnergy(
                                       inputHandler const* const inputShortcut,
           electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource,
                                                    double const beamEnergy ) :
      inputShortcut( inputShortcut ),
      electroweakCascadeSource( electroweakCascadeSource ),
      beamEnergy( beamEnergy ),
      cascadeOrderer( inputShortcut ),
      gluinoCascadeSet( inputShortcut,
                        electroweakCascadeSource,
                        &cascadeOrderer,
                        beamEnergy ),
      cascadeSetList( cascadeOrderer.getSdownTypeCascadesWithoutOrdering() )
  {
    for( std::vector< particlePointer >::const_iterator
         squarkIterator( inputShortcut->getSdownTypes()->begin() );
         inputShortcut->getSdownTypes()->end() > squarkIterator;
         ++squarkIterator )
    {
      cascadeSetList->push_back(
               new fullCascadeSetType::squarkSetType::sdownType( inputShortcut,
                                                               *squarkIterator,
                                                      electroweakCascadeSource,
                                                               &cascadeOrderer,
                                                             &gluinoCascadeSet,
                                                                beamEnergy ) );
    }
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "number of fullCascadeSetType::squarkSetType::sdownTypes = "
    << cascadeSetList->size();
    std::cout << std::endl;**/
    cascadeSetList = cascadeOrderer.getSupTypeCascadesWithoutOrdering();
    for( std::vector< particlePointer >::const_iterator
         squarkIterator( inputShortcut->getSupTypes()->begin() );
         inputShortcut->getSupTypes()->end() > squarkIterator;
         ++squarkIterator )
    {
      cascadeSetList->push_back(
               new fullCascadeSetType::squarkSetType::supType( inputShortcut,
                                                               *squarkIterator,
                                                      electroweakCascadeSource,
                                                               &cascadeOrderer,
                                                             &gluinoCascadeSet,
                                                                beamEnergy ) );
    }
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "number of fullCascadeSetType::squarkSetType::supTypes = "
    << cascadeSetList->size();
    std::cout << std::endl;**/
    cascadeSetList
    = cascadeOrderer.getNeutralinoColoredCascadesWithoutOrdering();
    for( std::vector< particlePointer >::const_iterator
         sparticleIterator( inputShortcut->getUnstableNeutralinos()->begin() );
         inputShortcut->getUnstableNeutralinos()->end() > sparticleIterator;
         ++sparticleIterator )
    {
      cascadeSetList->push_back(
          new fullCascadeSetType::gluinoOrElectroweakinoSetType::neutralinoSet(
                                                                 inputShortcut,
                                                            *sparticleIterator,
                                                      electroweakCascadeSource,
                                                               &cascadeOrderer,
                                                                beamEnergy ) );
    }
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "number of"
    << " fullCascadeSetType::gluinoOrElectroweakinoSetType::neutralinoSets = "
    << cascadeSetList->size();
    std::cout << std::endl;**/
    cascadeSetList
    = cascadeOrderer.getCharginoColoredCascadesWithoutOrdering();
    for( std::vector< particlePointer >::const_iterator
         sparticleIterator( inputShortcut->getCharginos()->begin() );
         inputShortcut->getCharginos()->end() > sparticleIterator;
         ++sparticleIterator )
    {
      cascadeSetList->push_back(
            new fullCascadeSetType::gluinoOrElectroweakinoSetType::charginoSet(
                                                                 inputShortcut,
                                                            *sparticleIterator,
                                                      electroweakCascadeSource,
                                                               &cascadeOrderer,
                                                                beamEnergy ) );
    }
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "number of"
    << " fullCascadeSetType::gluinoOrElectroweakinoSetType::charginoSets = "
    << cascadeSetList->size();
    std::cout << std::endl;**/

    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "fullCascadeSetsForOneBeamEnergy::fullCascadeSetsForOneBeamEnergy( "
    << inputShortcut << ", " << electroweakCascadeSource << ", " << beamEnergy
    << " ) finished.";
    std::cout << std::endl;**/
  }

  fullCascadeSetsForOneBeamEnergy::~fullCascadeSetsForOneBeamEnergy()
  {
    // does nothing.
  }


  fullCascadeSet*
  fullCascadeSetsForOneBeamEnergy::getFullCascadeSet(
                                        particlePointer const initialScolored )
  /* this returns the fullCascadeSet for the requested colored sparticle, or
   * NULL if we were asked for a sparticle that is not the gluino or in
   * inputShortcut->getSquarks() or ->getElectroweakinos() (though not set up
   * for using compound decays of electroweakinos directly (used just by the
   * cascades themselves for building compound cascades, with the intention
   * of only cascades beginning with colored sparticles being requested,
   * there seems to be no problem in letting them be given out).
   */
  {
    if( CppSLHA::PDG_code::gluino == initialScolored->get_PDG_code() )
    {
      return &gluinoCascadeSet;
    }
    else
    {
      bool notFoundYet( true );
      fullCascadeSet* returnPointer( NULL );
      // we look to see if we already have a fullCascadeSet for this sparticle:
      cascadeSetList = cascadeOrderer.getSdownTypeCascades();
      for( std::list< fullCascadeSet* >::iterator
           setIterator( cascadeSetList->begin() );
           ( notFoundYet
             &&
             ( cascadeSetList->end() != setIterator ) );
           ++setIterator )
      {
        if( initialScolored == (*setIterator)->getInitialSparticle() )
        {
          returnPointer = *setIterator;
          notFoundYet = false;
        }
      }
      if( notFoundYet )
      {
        cascadeSetList = cascadeOrderer.getSupTypeCascades();
        for( std::list< fullCascadeSet* >::iterator
             setIterator( cascadeSetList->begin() );
             ( notFoundYet
               &&
               ( cascadeSetList->end() != setIterator ) );
             ++setIterator )
        {
          if( initialScolored == (*setIterator)->getInitialSparticle() )
          {
            returnPointer = *setIterator;
            notFoundYet = false;
          }
        }
      }
      if( notFoundYet )
      {
        cascadeSetList = cascadeOrderer.getNeutralinoColoredCascades();
        for( std::list< fullCascadeSet* >::iterator
             setIterator( cascadeSetList->begin() );
             ( notFoundYet
               &&
               ( cascadeSetList->end() != setIterator ) );
             ++setIterator )
        {
          if( initialScolored == (*setIterator)->getInitialSparticle() )
          {
            returnPointer = *setIterator;
            notFoundYet = false;
          }
        }
      }
      if( notFoundYet )
      {
        cascadeSetList = cascadeOrderer.getCharginoColoredCascades();
        for( std::list< fullCascadeSet* >::iterator
             setIterator( cascadeSetList->begin() );
             ( notFoundYet
               &&
               ( cascadeSetList->end() != setIterator ) );
             ++setIterator )
        {
          if( initialScolored == (*setIterator)->getInitialSparticle() )
          {
            returnPointer = *setIterator;
            notFoundYet = false;
          }
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
