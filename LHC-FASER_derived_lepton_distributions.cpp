/*
 * LHC-FASER_derived_lepton_distributions.cpp
 *
 *  Created on: 05 Apr 5 2011
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
 *      LHC-FASER_base_lepton_distribution_stuff.hpp
 *      LHC-FASER_base_lepton_distribution_stuff.cpp
 *      LHC-FASER_derived_lepton_distributions.hpp
 *      LHC-FASER_derived_lepton_distributions.cpp
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

#include "LHC-FASER_derived_lepton_distributions.hpp"

namespace LHC_FASER
{
  flatNearMuonPlusAntimuon::flatNearMuonPlusAntimuon(
                                      readierForNewPoint* const readierPointer,
                                   CppSLHA::CppSLHA0 const* const spectrumData,
                     CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                  CppSLHA::particle_property_set const* const thirdParticle ) :
    leptonEnergyDistribution( readierPointer,
                              spectrumData,
                              firstParticle,
                              effectiveSquarkMass,
                              secondParticle,
                              thirdParticle,
                              NULL ),
    minToMaxConst( minToMaxSegment.addTerm( 0,
                                            0,
                                            NULL,
                                            1.0 ) )
  {
    segments.push_back( &minToMaxSegment );
  }

  flatNearMuonPlusAntimuon::~flatNearMuonPlusAntimuon()
  {
    // does nothing.
  }



  sameChiralityNearMuon::sameChiralityNearMuon(
                                      readierForNewPoint* const readierPointer,
                                   CppSLHA::CppSLHA0 const* const spectrumData,
                     CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                  CppSLHA::particle_property_set const* const thirdParticle ) :
    leptonEnergyDistribution( readierPointer,
                              spectrumData,
                              firstParticle,
                              effectiveSquarkMass,
                              secondParticle,
                              thirdParticle,
                              NULL ),
    minToMaxConst( minToMaxSegment.addTerm( 0,
                                            0 ) ),
    minToMaxLin( minToMaxSegment.addTerm( 1,
                                          0 ) )
  {
    segments.push_back( &minToMaxSegment );
  }

  sameChiralityNearMuon::~sameChiralityNearMuon()
  {

    // does nothing.

  }



  oppositeChiralityNearMuon::oppositeChiralityNearMuon(
                                      readierForNewPoint* const readierPointer,
                                   CppSLHA::CppSLHA0 const* const spectrumData,
                     CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                  CppSLHA::particle_property_set const* const thirdParticle ) :
    leptonEnergyDistribution( readierPointer,
                              spectrumData,
                              firstParticle,
                              effectiveSquarkMass,
                              secondParticle,
                              thirdParticle,
                              NULL ),
    minToMaxConst( minToMaxSegment.addTerm( 0,
                                            0 ) ),
    minToMaxLin( minToMaxSegment.addTerm( 1,
                                          0 ) )
  {
    segments.push_back( &minToMaxSegment );
  }

  oppositeChiralityNearMuon::~oppositeChiralityNearMuon()
  {

    // does nothing.

  }



  flatFarMuonPlusAntimuon::flatFarMuonPlusAntimuon(
                                      readierForNewPoint* const readierPointer,
                                   CppSLHA::CppSLHA0 const* const spectrumData,
                     CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                     CppSLHA::particle_property_set const* const thirdParticle,
                 CppSLHA::particle_property_set const* const fourthParticle ) :
    leptonEnergyDistribution( readierPointer,
                              spectrumData,
                              firstParticle,
                              effectiveSquarkMass,
                              secondParticle,
                              thirdParticle,
                              fourthParticle ),
    minToLkConst( minToLkSegment.addTerm( 0,
                                          0 ) ),
    minToLkLog( minToLkSegment.addTerm( 0,
                                        1 ) ),
    lkToHkConst( lkToHkSegment.addTerm( 0,
                                        0 ) ),
    hkToMaxConst( hkToMaxSegment.addTerm( 0,
                                          0 ) ),
    hkToMaxLog( hkToMaxSegment.addTerm( 0,
                                        1 ) )
  {
    segments.push_back( &minToLkSegment );
    segments.push_back( &lkToHkSegment );
    segments.push_back( &hkToMaxSegment );
  }

  flatFarMuonPlusAntimuon::~flatFarMuonPlusAntimuon()
  {

    // does nothing.

  }


  void
  flatFarMuonPlusAntimuon::calculateCoefficients()
  {
    double mlSq( thirdMass * thirdMass );
    double mXSq( fourthMass * fourthMass );

    productionFrameEnergy = ( ( mlSq - mXSq ) / ( 2.0 * thirdMass ) );
    minimumEnergy = ( ( thirdMass / firstMass ) * productionFrameEnergy );
    // there are 2 kinematic cases for determining the kink energies in terms
    // of the cascade masses:
    if( ( firstMass * thirdMass ) > ( secondMass * secondMass ) )
    {
      Elk = ( ( ( secondMass / firstMass ) * ( secondMass / thirdMass ) )
              * productionFrameEnergy );
      Ehk = ( ( ( firstMass / secondMass ) * ( thirdMass / secondMass ) )
              * productionFrameEnergy );
    }
    else
    {
      Elk = ( ( ( firstMass / secondMass ) * ( thirdMass / secondMass ) )
              * productionFrameEnergy );
      Ehk = ( ( ( secondMass / firstMass ) * ( secondMass / thirdMass ) )
              * productionFrameEnergy );
    }
    maximumEnergy = ( ( firstMass / thirdMass ) * productionFrameEnergy );

    // now we set up the coefficients of various terms:

    // MIN to LK segment:
    minToLkSegment.setSegmentRange( minimumEnergy,
                                    Elk );
    minToLkConst->setCoefficient( -log( minimumEnergy ) );
    minToLkLog->setCoefficient( 1.0 );

    // LK to HK segment:
    lkToHkSegment.setSegmentRange( Elk,
                                   Ehk );
    lkToHkConst->setCoefficient( log( ( maximumEnergy / Ehk ) ) );

    // HK to MAX segment:
    hkToMaxSegment.setSegmentRange( Ehk,
                                    maximumEnergy );
    hkToMaxConst->setCoefficient( log( maximumEnergy ) );
    hkToMaxLog->setCoefficient( -1.0 );

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: flat_far_muon_plus_antimuon_distribution:";
    std::cout << std::endl << "normalization = " << normalization;
    std::cout << std::endl << "minimumEnergy = " << minimumEnergy;
    std::cout << std::endl << "Elk = " << Elk;
    std::cout << std::endl << "Ehk = " << Ehk;
    std::cout << std::endl << "maximumEnergy = " << maximumEnergy;
    std::cout << std::endl;**/
  }



  sameChiralityFarMuon::sameChiralityFarMuon(
                                      readierForNewPoint* const readierPointer,
                                   CppSLHA::CppSLHA0 const* const spectrumData,
                     CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                     CppSLHA::particle_property_set const* const thirdParticle,
                 CppSLHA::particle_property_set const* const fourthParticle ) :
    leptonEnergyDistribution( readierPointer,
                              spectrumData,
                              firstParticle,
                              effectiveSquarkMass,
                              secondParticle,
                              thirdParticle,
                              fourthParticle ),
    minToLkInv( minToLkSegment.addTerm( -1,
                                        0 ) ),
    minToLkConst( minToLkSegment.addTerm( 0,
                                          0 ) ),
    minToLkLog( minToLkSegment.addTerm( 0,
                                        1 ) ),
    minToLkLin( minToLkSegment.addTerm( 1,
                                        0 ) ),
    lkToHkInv( lkToHkSegment.addTerm( -1,
                                      0 ) ),
    lkToHkConst( lkToHkSegment.addTerm( 0,
                                        0 ) ),
    hkToMaxInv( hkToMaxSegment.addTerm( -1,
                                        0 ) ),
    hkToMaxConst( hkToMaxSegment.addTerm( 0,
                                          0 ) ),
    hkToMaxLog( hkToMaxSegment.addTerm( 0,
                                        1 ) ),
    hkToMaxLin( hkToMaxSegment.addTerm( 1,
                                        0 ) )
  {
    segments.push_back( &minToLkSegment );
    segments.push_back( &lkToHkSegment );
    segments.push_back( &hkToMaxSegment );
  }

  sameChiralityFarMuon::~sameChiralityFarMuon()
  {

    // does nothing.

  }


  void
  sameChiralityFarMuon::calculateCoefficients()
  {
    double mq( firstMass );
    double mqSq( mq * mq );
    double mNSq( secondMass * secondMass );
    double mlSq( thirdMass * thirdMass );
    double mXSq( fourthMass * fourthMass );

    productionFrameEnergy = ( ( mlSq - mXSq ) / ( 2.0 * thirdMass ) );
    minimumEnergy
    = ( ( thirdMass / firstMass ) * productionFrameEnergy );

    // there are 2 kinematic cases for determining the kink energies in terms
    // of the cascade masses.
    if( ( firstMass * thirdMass ) > ( secondMass * secondMass ) )
    {
      Elk = ( ( ( secondMass / firstMass ) * ( secondMass / thirdMass ) )
              * productionFrameEnergy );
      Ehk = ( ( ( firstMass / secondMass ) * ( thirdMass / secondMass ) )
              * productionFrameEnergy );
    }
    else
    {
      Elk = ( ( ( firstMass / secondMass ) * ( thirdMass / secondMass ) )
              * productionFrameEnergy );
      Ehk = ( ( ( secondMass / firstMass ) * ( secondMass / thirdMass ) )
              * productionFrameEnergy );
    }
    maximumEnergy = ( ( firstMass / thirdMass ) * productionFrameEnergy );

    // now we set up the coefficients of various terms:

    // MIN to LK segment:
    minToLkSegment.setSegmentRange( minimumEnergy,
                                    Elk );

    minToLkLin->setCoefficient( ( 4.0 * mlSq * mNSq * mq )
                                / ( ( mlSq - mNSq ) * ( mNSq + mqSq )
                                    * ( mlSq - mXSq ) ) );

    minToLkLog->setCoefficient( ( 2.0 * mNSq * ( mqSq - mlSq ) )
                                / ( ( mNSq - mlSq ) * ( mNSq + mqSq ) ) );

    minToLkConst->setCoefficient(
          ( 2.0 * mNSq
            * ( 2.0 * mNSq - mqSq + mlSq
              + ( mqSq - mlSq ) * log( ( ( 2.0 * mq ) / ( mlSq - mXSq ) ) ) ) )
        / ( ( mNSq - mlSq ) * ( mNSq + mqSq ) ) );

    minToLkInv->setCoefficient( ( mNSq * ( 2.0 * mNSq - mqSq )
                                  * ( mlSq - mXSq ) )
                                / ( mq * ( mlSq - mNSq ) * ( mNSq + mqSq ) ) );

    // LK to HK segment:
    lkToHkSegment.setSegmentRange( Elk,
                                   Ehk );
    // there are 2 kinematic cases for determining the kink energies in terms
    // of the cascade masses.
    if( ( firstMass * thirdMass ) < ( secondMass * secondMass ) )
      // if the assumption that the mass spectrum was SPS1a-like was
      // incorrect...
    {
      lkToHkConst->setCoefficient(
              ( 2.0 * ( ( mlSq + mNSq ) * ( mNSq - mqSq )
                      + mNSq * ( mlSq - mqSq ) * log( ( mNSq / mqSq ) ) ) )
              / ( ( mNSq - mlSq ) * ( mNSq + mqSq ) ) );

      lkToHkInv->setCoefficient(
                            ( 2.0 * mNSq * ( mNSq - mqSq ) * ( mlSq - mXSq ) )
                            / ( mq * ( mlSq - mNSq ) * ( mNSq + mqSq ) ) );
    }
    else
    {

      lkToHkConst->setCoefficient(
                ( 2.0 * ( ( mqSq - mNSq ) * ( mlSq - mNSq )
                          + mNSq * ( mlSq - mqSq ) * log( ( mlSq / mNSq ) ) ) )
              / ( ( mNSq - mlSq ) * ( mNSq + mqSq ) ) );

      lkToHkInv->setCoefficient( 0.0 );
    }

    // HK to MAX segment:
    hkToMaxSegment.setSegmentRange( Ehk,
                                    maximumEnergy );

    hkToMaxLin->setCoefficient( -minToLkLin->getCoefficient() );

    hkToMaxLog->setCoefficient( -minToLkLog->getCoefficient() );

    hkToMaxConst->setCoefficient(
     ( 2.0 * mNSq * ( mlSq + mqSq
                      + ( mqSq - mlSq ) * log( ( 2.0 * mlSq )
                                               / ( ( mlSq - mXSq ) * mq ) ) ) )
          / ( ( mlSq - mNSq ) * ( mNSq + mqSq ) ) );

    hkToMaxInv->setCoefficient( ( mNSq * mq * ( mlSq - mXSq ) )
                                / ( ( mNSq - mlSq ) * ( mNSq + mqSq ) ) );

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: lepton_distribution_same_chirality_far_muon:";
    std::cout << std::endl << "normalization = " << normalization;
    std::cout << std::endl << "minimum_energy = " << minimum_energy;
    std::cout << std::endl << "Elk = " << Elk;
    std::cout << std::endl << "Ehk = " << Ehk;
    std::cout << std::endl << "maximum_energy = " << maximum_energy;
    std::cout << std::endl << "minToLkLin = " << minToLkLin;
    std::cout << std::endl << "minToLkLog = " << minToLkLog;
    std::cout << std::endl << "minToLkConst = " << minToLkConst;
    std::cout << std::endl << "minToLkInv = " << minToLkInv;
    std::cout << std::endl << "lkToHkConst = " << lkToHkConst;
    std::cout << std::endl << "lkToHkInv = " << lkToHkInv;
    std::cout << std::endl;**/
  }



  oppositeChiralityFarMuon::oppositeChiralityFarMuon(
                                      readierForNewPoint* const readierPointer,
                                   CppSLHA::CppSLHA0 const* const spectrumData,
                     CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                     CppSLHA::particle_property_set const* const thirdParticle,
                 CppSLHA::particle_property_set const* const fourthParticle ) :
    leptonEnergyDistribution( readierPointer,
                              spectrumData,
                              firstParticle,
                              effectiveSquarkMass,
                              secondParticle,
                              thirdParticle,
                              fourthParticle ),
    minToLkInv( minToLkSegment.addTerm( -1,
                                        0 ) ),
    minToLkConst( minToLkSegment.addTerm( 0,
                                          0 ) ),
    minToLkLog( minToLkSegment.addTerm( 0,
                                        1 ) ),
    minToLkLin( minToLkSegment.addTerm( 1,
                                        0 ) ),
    lkToHkInv( lkToHkSegment.addTerm( -1,
                                      0 ) ),
    lkToHkConst( lkToHkSegment.addTerm( 0,
                                        0 ) ),
    hkToMaxInv( hkToMaxSegment.addTerm( -1,
                                        0 ) ),
    hkToMaxConst( hkToMaxSegment.addTerm( 0,
                                          0 ) ),
    hkToMaxLog( hkToMaxSegment.addTerm( 0,
                                        1 ) ),
    hkToMaxLin( hkToMaxSegment.addTerm( 1,
                                        0 ) )
  {
    segments.push_back( &minToLkSegment );
    segments.push_back( &lkToHkSegment );
    segments.push_back( &hkToMaxSegment );
  }

  oppositeChiralityFarMuon::~oppositeChiralityFarMuon()
  {

    // does nothing.

  }


  void
  oppositeChiralityFarMuon::calculateCoefficients()
  {
    double mq( firstMass );
    double mqSq( mq * mq );
    double mNSq( secondMass * secondMass );
    double mlSq( thirdMass * thirdMass );
    double mXSq( fourthMass * fourthMass );

    productionFrameEnergy = ( ( mlSq - mXSq ) / ( 2.0 * thirdMass ) );
    minimumEnergy = ( ( thirdMass / firstMass ) * productionFrameEnergy );
    // there are 2 kinematic cases for determining the kink energies in terms
    // of the cascade masses.
    if( ( firstMass * thirdMass ) > ( secondMass * secondMass ) )
    {
      Elk = ( ( ( secondMass / firstMass ) * ( secondMass / thirdMass ) )
              * productionFrameEnergy );
      Ehk = ( ( ( firstMass / secondMass ) * ( thirdMass / secondMass ) )
              * productionFrameEnergy );
    }
    else
    {
      Elk = ( ( ( firstMass / secondMass ) * ( thirdMass / secondMass ) )
              * productionFrameEnergy );
      Ehk = ( ( ( secondMass / firstMass ) * ( secondMass / thirdMass ) )
              * productionFrameEnergy );
    }
    maximumEnergy = ( ( firstMass / thirdMass ) * productionFrameEnergy );


    // now we set up the coefficients of various terms:

    // MIN to LK segment:
    minToLkSegment.setSegmentRange( minimumEnergy,
                                    Elk );

    minToLkLin->setCoefficient(
        ( ( 4.0 * mlSq * mNSq * mq )
          / ( ( mNSq - mlSq ) * ( mNSq + mqSq ) * ( mlSq - mXSq ) ) ) );

    minToLkLog->setCoefficient(
        ( ( 2.0 * ( mNSq * mNSq - mlSq * mqSq ) )
          / ( ( mNSq - mlSq ) * ( mNSq + mqSq ) ) ) );

    minToLkConst->setCoefficient(
        ( ( 2.0 * ( mNSq * ( 2.0 * mNSq - mqSq + mlSq )
                  + ( mlSq * mqSq - mNSq * mNSq )
                    * log( ( ( 2.0 * mq ) / ( mlSq - mXSq ) ) ) ) )
          / ( ( mlSq - mNSq ) * ( mNSq + mqSq ) ) ) );

    minToLkInv->setCoefficient(
        ( ( mNSq * ( ( 2.0 * mNSq - mqSq ) * ( mlSq - mXSq ) ) )
          / ( mq * ( mNSq - mlSq ) * ( mNSq + mqSq ) ) ) );


    // LK to HK segment:
    lkToHkSegment.setSegmentRange( Elk,
                                   Ehk );

    // there are 2 kinematic cases for determining the kink energies in terms
    // of the cascade masses.
    if( ( firstMass * thirdMass ) < ( secondMass * secondMass ) )
      // if the assumption that the mass spectrum was SPS1a-like was
      // incorrect...
    {
      lkToHkConst->setCoefficient(
            ( ( 2.0 * ( ( mlSq + mNSq ) * ( mNSq - mqSq )
                        + ( mNSq * mNSq - mlSq * mqSq )
                          * log( ( mNSq / mqSq ) ) ) )
              / ( ( mlSq - mNSq ) * ( mNSq + mqSq ) ) ) );

      lkToHkInv->setCoefficient(
            ( ( 2.0 * mNSq * ( mNSq - mqSq ) * ( mlSq - mXSq ) )
              / ( mq * ( mNSq - mlSq ) * ( mNSq + mqSq ) ) ) );
    }
    else
    {
      lkToHkConst->setCoefficient(
            ( ( 2.0 * ( ( mqSq - mNSq ) * ( mlSq - mNSq )
                        + ( mNSq * mNSq - mlSq * mqSq )
                          * log( ( mlSq / mNSq ) ) ) )
              / ( ( mlSq - mNSq ) * ( mNSq + mqSq ) ) ) );

      lkToHkInv->setCoefficient( 0.0 );
    }


    // HK to MAX segment:
    hkToMaxSegment.setSegmentRange( Ehk,
                                    maximumEnergy );

    hkToMaxLin->setCoefficient( -minToLkLin->getCoefficient() );

    hkToMaxLog->setCoefficient( -minToLkLog->getCoefficient() );

    hkToMaxConst->setCoefficient(
        ( ( 2.0 * ( ( mNSq * mNSq - mlSq * mqSq )
                    * log( ( ( 2.0 * mlSq ) / ( ( mlSq - mXSq ) * mq ) ) )
                    - mNSq * ( mqSq + mlSq ) ) )
          / ( ( mlSq - mNSq ) * ( mNSq + mqSq ) ) ) );

    hkToMaxInv->setCoefficient(
        ( ( mNSq * mq * ( mlSq - mXSq ) )
          / ( ( mlSq - mNSq ) * ( mNSq + mqSq ) ) ) );

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: lepton_distribution_opposite_chirality_far_muon:";
    std::cout << std::endl << "normalization = " << normalization;
    std::cout << std::endl << "minimum_energy = " << minimum_energy;
    std::cout << std::endl << "Elk = " << Elk;
    std::cout << std::endl << "Ehk = " << Ehk;
    std::cout << std::endl << "maximum_energy = " << maximum_energy;
    std::cout << std::endl << "lkToHkConst = " << lkToHkConst;
    std::cout << std::endl
    << "log( ( mlSq / mNSq ) ) = " << log( ( mlSq / mNSq ) );
    std::cout << std::endl
    << "( 2.0 / ( ( mlSq - mNSq ) * ( mNSq + mqSq ) ) ) = "
    << ( 2.0 / ( ( mlSq - mNSq ) * ( mNSq + mqSq ) ) );
    std::cout << std::endl
    << "( ( mqSq - mNSq ) * ( mlSq - mNSq ) "
    << "+ ( mNSq * mNSq - mlSq * mqSq ) * log( ( mlSq / mNSq ) ) ) = "
    << ( ( mqSq - mNSq ) * ( mlSq - mNSq )
         + ( mNSq * mNSq - mlSq * mqSq ) * log( ( mlSq / mNSq ) ) );
    std::cout << std::endl;**/
  }



  HiggsMuonPlusAntimuon::HiggsMuonPlusAntimuon(
                                    readierForNewPoint* const readierPointer,
                                 CppSLHA::CppSLHA0 const* const spectrumData,
              CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const secondParticle,
              CppSLHA::particle_property_set const* const thirdParticle,
          CppSLHA::particle_property_set const* const fourthParticle ) :
    leptonEnergyDistribution( readierPointer,
                                spectrumData,
                                firstParticle,
                                effectiveSquarkMass,
                                secondParticle,
                                thirdParticle,
                                fourthParticle ),
    minToLkConst( minToLkSegment.addTerm( 0,
                                          0 ) ),
    minToLkLog( minToLkSegment.addTerm( 0,
                                        1 ) ),
    lkToHkConst( lkToHkSegment.addTerm( 0,
                                        0 ) ),
    hkToMaxConst( hkToMaxSegment.addTerm( 0,
                                          0 ) ),
    hkToMaxLog( hkToMaxSegment.addTerm( 0,
                                        1 ) )
  {
    segments.push_back( &minToLkSegment );
    segments.push_back( &lkToHkSegment );
    segments.push_back( &hkToMaxSegment );
  }

  HiggsMuonPlusAntimuon::~HiggsMuonPlusAntimuon()
  {
    // does nothing.
  }


  void
  HiggsMuonPlusAntimuon::calculateCoefficients()
  {
    double mq( firstMass );
    double mN( secondMass );
    double mNSq( mN * mN );
    double mH( thirdMass );
    double mHSq( mH * mH );
    double mX( fourthMass );
    double mXSq( mX * mX );

    productionFrameEnergy = ( 0.5 * mH );

    double EHN( ( mNSq + mHSq - mXSq ) / ( 2.0 * mN ) );
    double gammaNH( EHN / mH );
    double betaNH( sqrt( ( 1.0 - ( 1.0 / ( gammaNH * gammaNH ) ) ) ) );
    double EmaxN( gammaNH * ( 1.0 + betaNH ) * productionFrameEnergy );
    double EminN( mHSq / ( 4.0 * EmaxN ) );
    // using ( EmaxN * EminN ) = ( mHSq / 4 )

    minimumEnergy = ( ( mN / mq ) * EminN );

    // there are 2 kinematic cases for determining the kink energies in terms
    // of the cascade masses:
    // 1st make an assumption:
    Elk = ( ( mN / mq ) * EmaxN );
    Ehk = ( ( mq / mN ) * EminN );
    // then check:
    if( Elk > Ehk )
    {
      Elk = Ehk;
      Ehk = ( ( mN / mq ) * EmaxN );
    }
    maximumEnergy = ( ( mq / mN ) * EmaxN );


    // now we set up the coefficients of various terms:

    // MIN to LK segment:
    minToLkSegment.setSegmentRange( minimumEnergy,
                                    Elk );
    minToLkConst->setCoefficient( -log( minimumEnergy ) );
    minToLkLog->setCoefficient( 1.0 );

    // LK to HK segment:
    lkToHkSegment.setSegmentRange( Elk,
                                   Ehk );
    lkToHkConst->setCoefficient( log( ( maximumEnergy / Ehk ) ) );

    // HK to MAX segment:
    hkToMaxSegment.setSegmentRange( Ehk,
                                    maximumEnergy );
    hkToMaxConst->setCoefficient( log( maximumEnergy ) );
    hkToMaxLog->setCoefficient( -1.0 );

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: lepton_distribution_Higgs_muon_plus_antimuon_distribution:";
    std::cout << std::endl << "normalization = " << normalization;
    std::cout << std::endl << "minimum_energy = " << minimum_energy;
    std::cout << std::endl << "Elk = " << Elk;
    std::cout << std::endl << "Ehk = " << Ehk;
    std::cout << std::endl << "maximum_energy = " << maximum_energy;
    std::cout << std::endl << "minToLkLog = " << minToLkLog;
    std::cout << std::endl << "minToLkConst = " << minToLkConst;
    std::cout << std::endl << "lkToHkConst = " << lkToHkConst;
    std::cout << std::endl;**/
  }



  zHandedMuon::zHandedMuon( readierForNewPoint* const readier,
                            CppSLHA::CppSLHA0 const* const spectrumData,
              CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const secondParticle,
              CppSLHA::particle_property_set const* const thirdParticle,
          CppSLHA::particle_property_set const* const fourthParticle,
                            bool const negativeMuonIsSameHandednessAsJet,
                            bool const shouldSumOverHandedness ) :
    leptonEnergyDistribution( readier,
                              spectrumData,
                              firstParticle,
                              effectiveSquarkMass,
                              secondParticle,
                              thirdParticle,
                              fourthParticle ),
    negativeMuonIsSameHandednessAsJet( negativeMuonIsSameHandednessAsJet ),
    shouldNotSumOverHandedness( !shouldSumOverHandedness ),
    minToLkInv( minToLkSegment.addTerm( -1,
                                        0 ) ),
    minToLkConst( minToLkSegment.addTerm( 0,
                                          0 ) ),
    minToLkLog( minToLkSegment.addTerm( 0,
                                        1 ) ),
    minToLkLin( minToLkSegment.addTerm( 1,
                                        0 ) ),
    minToLkLinlog( minToLkSegment.addTerm( 1,
                                           1 ) ),
    minToLkSq( minToLkSegment.addTerm( 2,
                                       0 ) ),
    lkToHkConst( lkToHkSegment.addTerm( 0,
                                        0 ) ),
    lkToHkLin( lkToHkSegment.addTerm( 1,
                                      0 ) ),
    lkToHkSq( lkToHkSegment.addTerm( 2,
                                     0 ) ),
    hkToMaxInv( hkToMaxSegment.addTerm( -1,
                                        0 ) ),
    hkToMaxConst( hkToMaxSegment.addTerm( 0,
                                          0 ) ),
    hkToMaxLog( hkToMaxSegment.addTerm( 0,
                                        1 ) ),
    hkToMaxLin( hkToMaxSegment.addTerm( 1,
                                        0 ) ),
    hkToMaxLinlog( hkToMaxSegment.addTerm( 1,
                                           1 ) ),
    hkToMaxSq( hkToMaxSegment.addTerm( 2,
                                       0 ) )
  {
    segments.push_back( &minToLkSegment );
    segments.push_back( &lkToHkSegment );
    segments.push_back( &hkToMaxSegment );
  }

  zHandedMuon::~zHandedMuon()
  {
    // does nothing.
  }


  void
  zHandedMuon::calculateCoefficients()
  {
    /* it doesn't matter if the neutralino mass is negative, since only the
     * absolute square appears, & the couplings are only affected by it as an
     * overall phase.  I (BOL) assume no CP violation, so the neutralino mixing
     * matrices are strictly real. just as a quirk of the way I did my
     * workings, axialCouplingFactor needs to be +1 when the Z couples as a
     * vector to the neutralinos, & -1 when it couples as an axial vector.
     */
    if( ( ( 0.0 < secondParticle->get_mass() )
          &&
          ( 0.0 > fourthParticle->get_mass() ) )
        ||
        ( ( 0.0 < secondParticle->get_mass() )
          &&
          ( 0.0 > fourthParticle->get_mass() ) ) )
      // if the relative phases end up such that the product of couplings is
      // purely imaginary...
    {
      couplesAsVectorNotAxial = false;
    }
    else
    {
      couplesAsVectorNotAxial = true;
    }
    mQ = firstMass;  // this is just for compact equations.
    mQsq = ( mQ * mQ );
    mQcu = ( mQ * mQsq );
    mC = secondMass;  // this is just for compact equations.
    mCsq = ( mC * mC );
    mVsq = ( thirdMass * thirdMass );
    mCmX = ( mC * fourthMass );
    mXsq = ( fourthMass * fourthMass );

    EV = ( ( mCsq + mVsq - mXsq ) / ( 2.0 * mC ) );
    ElMin = ( 0.5 * EV * ( 1.0 - sqrt( 1 - ( mVsq / ( EV * EV ) ) ) ) );
    minimumEnergy = ( ElMin * ( mC / mQ ) );
    Elk = ( ElMin * ( mQ / mC ) );
    Ehk = ( mVsq / ( 4.0 * Elk ) );
    maximumEnergy = ( mVsq / ( 4.0 * minimumEnergy ) );
    if( Ehk < Elk )
      /* if we're in the regime where the kinematics mean that there is no
       * range of energies which can be reached for all squark-rest-frame polar
       * muon momentum angles...
       */
    {
      // we swap the energies (borrowing currentLkToHkLinCoefficient for a
      // moment) & note that we had to swap them:
      currentLkToHkLinCoefficient = Elk;
      Elk = Ehk;
      Ehk = currentLkToHkLinCoefficient;
      cosinesLimitedByEnergy = true;
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "cosinesLimitedByEnergy = true";
      std::cout << std::endl;**/
    }
    else
    {
      cosinesLimitedByEnergy = false;
      // for comparison with the Mathematica notebook: underlap is Ehk > Elk
      // => if( cosinesLimitedByEnergy ){ [use overlap] }

      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "cosinesLimitedByEnergy = false";
      std::cout << std::endl;**/
    }
    mQCsqDiff = ( mQsq - mCsq );
    mQCsqSum = ( mQsq + mCsq );
    EightmCEVElMin = ( 8.0 * mC * EV * ElMin );
    lnElMin = log( ElMin );
    ElMinsq = ( ElMin * ElMin );
    lnmQC = log( ( mQ / mC ) );
    lnVsqOverFourElMinsq = log( ( mVsq / ( 4.0 * ElMinsq ) ) );
    EightmCXElMinmQsq = ( 8.0 * mC * fourthMass * ElMin * mQsq );
    lnEmin = log( minimumEnergy );
    lnElk = log( Elk );
    lnEhk = log( Ehk );
    lnEmax = log( maximumEnergy );

    minToLkSegment.setSegmentRange( minimumEnergy,
                                     Elk );
    lkToHkSegment.setSegmentRange( Elk,
                                    Ehk );
    hkToMaxSegment.setSegmentRange( Ehk,
                                     maximumEnergy );

    /* symmetric in lepton-antilepton =>
     *     ( Vv + jL * Aa )^2 * VvPlusJjAaAllSqSymDistribution
     *     - 2 * Jj * Vv * Aa * TwiceVvAaSymDistribution
     *     + ( Vv^( 2 ) - Aa^( 2 ) ) * VvSqMinusAaSqSymDistribution
     * antisymmetric in lepton-antilepton =>
     *     ( Vv + jL * Aa )^2 * VvPlusJjAaAllSqAntiDistribution
     *     + ( Vv^( 2 ) + Aa^( 2 ) ) * VvSqPlusAaSqAntiDistribution
     *     + ( Vv^( 2 ) - Aa^( 2 ) ) * VvSqMinusAaSqAntiDistribution
     * summing over handedness =>
     * symmetric
     * lepton same handedness as jet =>
     * symmetric + antisymmetric
     * lepton opposite handedness to jet =>
     * symmetric - antisymmetric
     *
     * for Z bosons, either Vv is 1 & Aa is 0 or Aa is 1 & Vv is 0.
     * jL & Jj are each either +1 or -1.
     */


    /* ends up in the order:
     * low (MINtoLK),
     * over (LKtoHK, true == cosinesLimitedByEnergy),
     * under (LKtoHK, false == cosinesLimitedByEnergy),
     * high (HKtoMAX)
     */


    // VvPlusJjAaAllSqSymDistribution:

    // MIN to LK:

    currentMinToLkConstCoefficient
    = ( ElMin * mQsq * mQCsqDiff
        * ( EightmCEVElMin * ( lnElMin - lnmQC ) * mVsq
            - 4.0 * ElMin * mC * ( EightmCEVElMin - 4.0 * ElMinsq * mC ) ) );
    // right!

    currentMinToLkLogCoefficient
    = ( -EightmCEVElMin * ElMin * mQsq * mQCsqDiff * mVsq );
    // right!

    currentMinToLkLinCoefficient
    = ( 4.0 * EightmCEVElMin * ElMin * mQcu * mQCsqDiff );
    // right!

    currentMinToLkSqCoefficient
    = ( -16.0 * ElMinsq * mQsq * mQsq * mQCsqDiff );
    // right!


    // LK to HK:

    if( cosinesLimitedByEnergy )
    {
      currentLkToHkConstCoefficient
      = ( mQsq * mQCsqDiff
          * ( EightmCEVElMin
              * ( mC * ( mVsq - 4.0 * ElMinsq )
                  - ElMin * lnVsqOverFourElMinsq * mVsq )
              + mCsq * ( 16.0 * ElMinsq * ElMinsq - mVsq * mVsq ) ) );
      // dunno...
      currentLkToHkLinCoefficient = 0.0;
      // dunno...
      currentLkToHkSqCoefficient = 0.0;
      // dunno...
    }
    else
    {
      currentLkToHkConstCoefficient
      = ( -2.0 * EightmCEVElMin * ElMin * lnmQC * mQsq * mQCsqDiff * mVsq );
      // right!

      currentLkToHkLinCoefficient
      = ( 4.0 * EightmCEVElMin * ElMin * mQ * mQCsqDiff * mQCsqDiff );
      // right!

      currentLkToHkSqCoefficient
      = ( -16.0 * ElMinsq * mQCsqDiff * mQCsqDiff * mQCsqSum );
      // right!
    }


    // HK to MAX:

    currentHkToMaxConstCoefficient
    = ( mQsq * mQCsqDiff * mVsq
        * ( EightmCEVElMin
            * ( mC - ElMin * ( lnmQC + lnVsqOverFourElMinsq + lnElMin ) )
            - mCsq * mVsq ) );
    // right!

    currentHkToMaxLogCoefficient
    = ( EightmCEVElMin * ElMin * mQsq * mQCsqDiff * mVsq );
    // right!

    currentHkToMaxLinCoefficient
    = ( -4.0 * EightmCEVElMin * ElMin * mCsq * mQ * mQCsqDiff );
    // right!

    currentHkToMaxSqCoefficient = ( 16.0 * ElMinsq * mCsq * mCsq * mQCsqDiff );
    // right!


    minToLkInv->setCoefficient( 0.0 );
    minToLkConst->setCoefficient( currentMinToLkConstCoefficient );
    minToLkLog->setCoefficient( currentMinToLkLogCoefficient );
    minToLkLin->setCoefficient( currentMinToLkLinCoefficient );
    minToLkLinlog->setCoefficient( 0.0 );
    minToLkSq->setCoefficient( currentMinToLkSqCoefficient );

    lkToHkConst->setCoefficient( currentLkToHkConstCoefficient );
    lkToHkLin->setCoefficient( currentLkToHkLinCoefficient );
    lkToHkSq->setCoefficient( currentLkToHkSqCoefficient );

    hkToMaxInv->setCoefficient( 0.0 );
    hkToMaxConst->setCoefficient( currentHkToMaxConstCoefficient );
    hkToMaxLog->setCoefficient( currentHkToMaxLogCoefficient );
    hkToMaxLin->setCoefficient( currentHkToMaxLinCoefficient );
    hkToMaxLinlog->setCoefficient( 0.0 );
    hkToMaxSq->setCoefficient( currentHkToMaxSqCoefficient );

    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "VvPlusJjAaAllSqSymDistribution:"
    << std::endl << "MINtoLK: "
    << std::endl << minToLkInv->getCoefficient() << "/E"
    << std::endl << minToLkConst->getCoefficient()
    << std::endl << minToLkLog->getCoefficient() << "L"
    << std::endl << minToLkLin->getCoefficient() << "E"
    << std::endl << minToLkLinlog->getCoefficient() << "EL"
    << std::endl << minToLkSq->getCoefficient() << "E^2"
    << std::endl << "LKtoHK: "
    << std::endl << lkToHkConst->getCoefficient()
    << std::endl << lkToHkLin->getCoefficient() << "E"
    << std::endl << lkToHkSq->getCoefficient() << "E^2"
    << std::endl << "HKtoMAX: "
    << std::endl << hkToMaxInv->getCoefficient() << "/E"
    << std::endl << hkToMaxConst->getCoefficient()
    << std::endl << hkToMaxLog->getCoefficient() << "L"
    << std::endl << hkToMaxLin->getCoefficient() << "E"
    << std::endl << hkToMaxLinlog->getCoefficient() << "EL"
    << std::endl << hkToMaxSq->getCoefficient() << "E^2";
    std::cout << std::endl;**/


    // VvSqMinusAaSqSymDistribution:

    // MIN to LK:

    // debugging:
    //currentMinToLkInvCoefficient = 0.0;
    // right!

    currentMinToLkConstCoefficient
    = ( EightmCXElMinmQsq * ElMin * ( lnmQC - lnElMin ) * mQCsqDiff * mVsq );
    // right!

    currentMinToLkLogCoefficient
    = ( EightmCXElMinmQsq * ElMin * mQCsqDiff * mVsq );
    // right!

    // debugging:
    //currentMinToLkLinCoefficient = 0.0;
    // right!

    // debugging:
    //currentMinToLkLinLogCoefficient = 0.0;
    // right!

    // debugging:
    //currentMinToLkSqCoefficient = 0.0;
    // right!


    // LK to HK:

    if( cosinesLimitedByEnergy )
    {
      currentLkToHkConstCoefficient = ( EightmCXElMinmQsq * ElMin * mQCsqDiff
                                        * mVsq * lnVsqOverFourElMinsq );
      // dunno...
    }
    else
    {
      currentLkToHkConstCoefficient
      = ( 2.0 * EightmCXElMinmQsq * ElMin * lnmQC * mQCsqDiff * mVsq );
      // right!
    }

    // debugging:
    //currentLkToHkLinCoefficient = 0.0;
    // right!

    // debugging:
    //currentLkToHkSqCoefficient = 0.0;
    // right!


    // HK to MAX:

    // debugging:
    //currentHkToMaxInvCoefficient = 0.0;
    // right!

    currentHkToMaxConstCoefficient
    = ( EightmCXElMinmQsq * ElMin * mQCsqDiff * mVsq
        * ( lnmQC + lnVsqOverFourElMinsq + lnElMin ) );
    // right!

    currentHkToMaxLogCoefficient = -currentMinToLkLogCoefficient;
    // right!

    // debugging:
    //currentHkToMaxLinCoefficient = 0.0;
    // right!

    // debugging:
    //currentHkToMaxLinLogCoefficient = 0.0;
    // right!

    // debugging:
    //currentHkToMaxSqCoefficient = 0.0;
    // right!


    if( couplesAsVectorNotAxial )
    {
      minToLkConst->addToCoefficient( currentMinToLkConstCoefficient );
      minToLkLog->addToCoefficient( currentMinToLkLogCoefficient );

      lkToHkConst->addToCoefficient( currentLkToHkConstCoefficient );

      hkToMaxConst->addToCoefficient( currentHkToMaxConstCoefficient );
      hkToMaxLog->addToCoefficient( currentHkToMaxLogCoefficient );
    }
    else
    {
      minToLkConst->addToCoefficient( -currentMinToLkConstCoefficient );
      minToLkLog->addToCoefficient( -currentMinToLkLogCoefficient );

      lkToHkConst->addToCoefficient( -currentLkToHkConstCoefficient );

      hkToMaxConst->addToCoefficient( -currentHkToMaxConstCoefficient );
      hkToMaxLog->addToCoefficient( -currentHkToMaxLogCoefficient );
    }

    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "VvSqMinusAaSqSymDistribution:"
    << std::endl << "MINtoLK: "
    << std::endl << currentMinToLkInvCoefficient << "/E"
    << std::endl << currentMinToLkConstCoefficient
    << std::endl << currentMinToLkLogCoefficient << "L"
    << std::endl << currentMinToLkLinCoefficient << "E"
    << std::endl << currentMinToLkLinLogCoefficient << "EL"
    << std::endl << currentMinToLkSqCoefficient << "E^2"
    << std::endl << "LKtoHK: "
    << std::endl << currentLkToHkConstCoefficient
    << std::endl << currentLkToHkLinCoefficient << "E"
    << std::endl << currentLkToHkSqCoefficient << "E^2"
    << std::endl << "HKtoMAX: "
    << std::endl << currentHkToMaxInvCoefficient << "/E"
    << std::endl << currentHkToMaxConstCoefficient
    << std::endl << currentHkToMaxLogCoefficient << "L"
    << std::endl << currentHkToMaxLinCoefficient << "E"
    << std::endl << currentHkToMaxLinLogCoefficient << "EL"
    << std::endl << currentHkToMaxSqCoefficient << "E^2";**/


    if( shouldNotSumOverHandedness )
    {
      // VvPlusJjAaAllSqAntiDistribution:

      // MIN to LK:

      // debugging:
      //currentMinToLkInvCoefficient = 0.0;
      // right!

      currentMinToLkConstCoefficient
      =  ( ElMin * ( EightmCEVElMin * ( lnElMin - lnmQC )
                     - 16.0 * ElMinsq * mC ) * mQsq * mQCsqDiff * mVsq );
      // right!

      currentMinToLkLogCoefficient
      = ( -EightmCEVElMin * ElMin * mQsq * mQCsqDiff * mVsq );
      // right!

      currentMinToLkLinCoefficient
      = ( 16.0 * ElMinsq * mQcu * mQCsqDiff * mVsq );
      // right!

      // debugging:
      //currentMinToLkLinLogCoefficient = 0.0;
      // right!

      // debugging:
      //currentMinToLkSqCoefficient = 0.0;
      // right!


      // LK to HK:

      if( cosinesLimitedByEnergy )
      {
        currentLkToHkConstCoefficient
        = ( ElMin * mQsq * mQCsqDiff * mVsq
            * ( 4.0 * mC * ( mVsq - 4.0 * ElMinsq )
                - EightmCEVElMin * lnVsqOverFourElMinsq ) );
        // dunno...
        currentLkToHkLinCoefficient = 0.0;
        // dunno...
      }
      else
      {
        currentLkToHkConstCoefficient
        = ( -2.0 * EightmCEVElMin * ElMin * lnmQC * mQsq * mQCsqDiff * mVsq );
        // right!
        currentLkToHkLinCoefficient
        = ( 16.0 * ElMinsq * mQ * mQCsqDiff * mQCsqDiff * mVsq );
        // right!
      }

      // debugging:
      //currentLkToHkSqCoefficient = 0.0;


      // HK to MAX:

      //currentHkToMaxInvCoefficient = 0.0;
      // right!

      currentHkToMaxConstCoefficient
      = ( ElMin * mQsq * mQCsqDiff * mVsq
          * ( 4.0 * mC * mVsq
              - EightmCEVElMin
                * ( lnmQC + lnVsqOverFourElMinsq + lnElMin ) ) );
      // right!

      currentHkToMaxLogCoefficient
      = ( EightmCEVElMin * ElMin * mQsq * mQCsqDiff * mVsq );
      // right!

      currentHkToMaxLinCoefficient
      = ( -16.0 * ElMinsq * mCsq * mQ * mQCsqDiff * mVsq );
      // right!

      // debugging:
      //currentHkToMaxLinLogCoefficient = 0.0;
      // right!

      // debugging:
      //currentHkToMaxSqCoefficient = 0.0;
      // right!


      minToLkConst->addToCoefficient( currentMinToLkConstCoefficient );
      minToLkLog->addToCoefficient( currentMinToLkLogCoefficient );
      minToLkLin->addToCoefficient( currentMinToLkLinCoefficient );

      lkToHkConst->addToCoefficient( currentLkToHkConstCoefficient );
      lkToHkLin->addToCoefficient( currentLkToHkLinCoefficient );

      hkToMaxConst->addToCoefficient( currentHkToMaxConstCoefficient );
      hkToMaxLog->addToCoefficient( currentHkToMaxLogCoefficient );
      hkToMaxLin->addToCoefficient( currentHkToMaxLinCoefficient );

      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "VvPlusJjAaAllSqAntiDistribution:"
      << std::endl << "MINtoLK: "
      << std::endl << currentMinToLkInvCoefficient << "/E"
      << std::endl << currentMinToLkConstCoefficient
      << std::endl << currentMinToLkLogCoefficient << "L"
      << std::endl << currentMinToLkLinCoefficient << "E"
      << std::endl << currentMinToLkLinLogCoefficient << "EL"
      << std::endl << currentMinToLkSqCoefficient << "E^2"
      << std::endl << "LKtoHK: "
      << std::endl << currentLkToHkConstCoefficient
      << std::endl << currentLkToHkLinCoefficient << "E"
      << std::endl << currentLkToHkSqCoefficient << "E^2"
      << std::endl << "HKtoMAX: "
      << std::endl << currentHkToMaxInvCoefficient << "/E"
      << std::endl << currentHkToMaxConstCoefficient
      << std::endl << currentHkToMaxLogCoefficient << "L"
      << std::endl << currentHkToMaxLinCoefficient << "E"
      << std::endl << currentHkToMaxLinLogCoefficient << "EL"
      << std::endl << currentHkToMaxSqCoefficient << "E^2";**/


      // VvSqPlusAaSqAntiDistribution:

      // MIN to LK:

      currentMinToLkInvCoefficient
      = ( 4.0 * ElMinsq * mCsq * mQcu * mVsq * mVsq );
      // right!

      currentMinToLkConstCoefficient
      = ( ElMin * mQsq * mVsq
          * ( 32.0 * ElMinsq * mC * mQsq
              + 2.0 * EightmCEVElMin * ( mCsq + ( lnmQC - lnElMin ) * mQsq )
              + 8.0 * ElMin * mCsq
                * ( 2.0 * mCsq + ( lnmQC - lnElMin ) * mQCsqSum )
              - 4.0 * mC * mQCsqSum * mVsq ) );
      // right!

      currentMinToLkLogCoefficient
      = ( 2.0 * ElMin * mQsq * mVsq
              * ( EightmCEVElMin * mQsq + 4.0 * ElMin * mCsq * mQCsqSum ) );
      // right!

      currentMinToLkLinCoefficient
      = ( 2.0 * mQcu * mVsq
          * ( 2.0 * mCsq * mVsq
              - EightmCEVElMin * mC
              - 8.0 * ElMin
                * ( mC * mCsq
                    + 2.0 * ElMin
                      * ( ( lnElMin - lnmQC ) * mCsq + mQsq ) ) ) );
      // right!

      currentMinToLkLinLogCoefficient
      = ( 32.0 * ElMinsq * mCsq * mQcu * mVsq );
      // right!

      // debugging:
      //currentMinToLkSqCoefficient = 0.0;
      // right!


      // LK to HK:

      if( cosinesLimitedByEnergy )
      {
        currentLkToHkConstCoefficient
        = ( 2.0 * ElMin * mQsq * mVsq
            * ( lnVsqOverFourElMinsq
                * ( EightmCEVElMin * mQsq + 4.0 * ElMin * mCsq * mQCsqSum )
                + 2.0 * mC * ( mCsq + 3.0 * mQsq )
                  * ( 4.0 * ElMinsq - mVsq ) ) );
        // dunno...
        currentLkToHkLinCoefficient
        = ( 2.0 * mC * mQcu
            * ( 4.0 * ElMinsq
                * ( EightmCEVElMin
                    + 8.0 * ElMin * mC * ( mC - ElMin ) )
                - mVsq * ( EightmCEVElMin
                           + 8.0 * ElMin * mC
                             * ( mC - 2.0 * ElMin * lnVsqOverFourElMinsq ) )
                + 2.0 * mC * mVsq * mVsq ) );
        // dunno...
      }
      else
      {
        currentLkToHkConstCoefficient
           = ( 2.0 * ElMin * mQsq * mVsq
               * ( 2.0 * mQsq * lnmQC
                   * ( EightmCEVElMin + 4.0 * ElMin * mCsq )
                   - mQCsqDiff * ( 8.0 * ElMin * mCsq + EightmCEVElMin )
                   + 8.0 * ElMin * mCsq * mCsq * lnmQC ) );
        // right!
        currentLkToHkLinCoefficient = ( 32.0 * ElMinsq * mQcu * mVsq
                                        * ( 2.0 * mCsq * lnmQC - mQCsqDiff ) );
        // right!
      }

      // debugging:
      //currentLkToHkSqCoefficient = 0.0;
      // right!


      // HK to MAX:

      currentHkToMaxInvCoefficient = -currentMinToLkInvCoefficient;
      // depends...

      currentHkToMaxConstCoefficient
      = ( 2.0 * ElMin * mQsq * mVsq
          * ( 4.0 * ElMin * mCsq
              * ( ( lnmQC + lnVsqOverFourElMinsq + lnElMin ) * mQCsqSum
                  - 2.0 * mQsq )
              + 8.0 * ElMinsq * mC * mQCsqSum
              + mQsq * ( EightmCEVElMin
                         * ( lnmQC + lnVsqOverFourElMinsq + lnElMin - 1.0 )
                         - 4.0 * mC * mVsq ) ) );
      // right!

      currentHkToMaxLogCoefficient
      = ( -2.0 * ElMin * mQsq * mVsq
          * ( EightmCEVElMin * mQsq + 4.0 * ElMin * mCsq * mQCsqSum ) );
      // right!

      currentHkToMaxLinCoefficient
      = ( 8.0 * ElMinsq * mC * mQcu
          * ( EightmCEVElMin + 8.0 * ElMin * mC * ( mC - ElMin )
              + 4.0 * mC * mVsq
                * ( 1.0 + lnmQC + lnVsqOverFourElMinsq + lnElMin ) ) );
      // right!

      currentHkToMaxLinLogCoefficient = -currentMinToLkLinLogCoefficient;
      // right!

      // debugging:
      //currentHkToMaxSqCoefficient = 0.0;
      // right!


      minToLkInv->addToCoefficient( currentMinToLkInvCoefficient );
      minToLkConst->addToCoefficient( currentMinToLkConstCoefficient );
      minToLkLog->addToCoefficient( currentMinToLkLogCoefficient );
      minToLkLin->addToCoefficient( currentMinToLkLinCoefficient );
      minToLkLinlog->addToCoefficient( currentMinToLkLinLogCoefficient );

      lkToHkConst->addToCoefficient( currentLkToHkConstCoefficient );
      lkToHkLin->addToCoefficient( currentLkToHkLinCoefficient );

      hkToMaxInv->addToCoefficient( currentHkToMaxInvCoefficient );
      hkToMaxConst->addToCoefficient( currentHkToMaxConstCoefficient );
      hkToMaxLog->addToCoefficient( currentHkToMaxLogCoefficient );
      hkToMaxLin->addToCoefficient( currentHkToMaxLinCoefficient );
      hkToMaxLinlog->addToCoefficient( currentHkToMaxLinLogCoefficient );

      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "VvSqPlusAaSqAntiDistribution:"
      << std::endl << "MINtoLK: "
      << std::endl << currentMinToLkInvCoefficient << "/E"
      << std::endl << currentMinToLkConstCoefficient
      << std::endl << currentMinToLkLogCoefficient << "L"
      << std::endl << currentMinToLkLinCoefficient << "E"
      << std::endl << currentMinToLkLinLogCoefficient << "EL"
      << std::endl << currentMinToLkSqCoefficient << "E^2"
      << std::endl << "LKtoHK: "
      << std::endl << currentLkToHkConstCoefficient
      << std::endl << currentLkToHkLinCoefficient << "E"
      << std::endl << currentLkToHkSqCoefficient << "E^2"
      << std::endl << "HKtoMAX: "
      << std::endl << currentHkToMaxInvCoefficient << "/E"
      << std::endl << currentHkToMaxConstCoefficient
      << std::endl << currentHkToMaxLogCoefficient << "L"
      << std::endl << currentHkToMaxLinCoefficient << "E"
      << std::endl << currentHkToMaxLinLogCoefficient << "EL"
      << std::endl << currentHkToMaxSqCoefficient << "E^2";**/


      // VvSqMinusAaSqAntiDistribution

      // MIN to LK:

      // debugging:
      //currentMinToLkInvCoefficient = 0.0;
      // right!

      currentMinToLkConstCoefficient
      = ( EightmCXElMinmQsq * ElMin * mVsq
          * ( ( 1.0 - lnElMin + lnmQC ) * mQCsqSum - mQCsqDiff ) );
      // right!

      currentMinToLkLogCoefficient
      = ( EightmCXElMinmQsq * ElMin * mQCsqSum * mVsq );
      // right!

      currentMinToLkLinCoefficient
      = ( -2.0 * EightmCXElMinmQsq * mC * mQ * mVsq );
      // right!

      // debugging:
      //currentMinToLkLinLogCoefficient = 0.0;
      // right!

      // debugging:
      //currentMinToLkSqCoefficient = 0.0;
      // right!


      // LK to HK:

      if( cosinesLimitedByEnergy )
      {
        currentLkToHkConstCoefficient = ( EightmCXElMinmQsq * ElMin * mQCsqSum
                                          * mVsq * lnVsqOverFourElMinsq );
        // dunno...
        currentLkToHkLinCoefficient
        = ( 2.0 * EightmCXElMinmQsq * mC * mQ * ( 4.0 * ElMinsq - mVsq ) );
        // dunno...
      }
      else
      {
        currentLkToHkConstCoefficient
        = ( 2.0 * EightmCXElMinmQsq * ElMin * mVsq
            * ( lnmQC * mQCsqSum - mQCsqDiff ) );
        // right!
        currentLkToHkLinCoefficient = 0.0;
        // right!
      }

      // debugging:
      //currentLkToHkSqCoefficient = 0.0;
      // right!


      // HK to MAX:

      // debugging:
      //currentHkToMaxInvCoefficient = 0.0;
      // right!

      currentHkToMaxConstCoefficient
      = ( EightmCXElMinmQsq * ElMin * mVsq
          * ( ( lnVsqOverFourElMinsq + lnElMin + lnmQC - 1.0 ) * mQCsqSum
              - mQCsqDiff ) );
      // right!

      currentHkToMaxLogCoefficient = -currentMinToLkLogCoefficient;
      // right!

      currentHkToMaxLinCoefficient
      = ( 8.0 * EightmCXElMinmQsq * mC * mQ * ElMinsq );
      // right!

      // debugging:
      //currentHkToMaxLinLogCoefficient = 0.0;
      // right!

      // debugging:
      //currentHkToMaxSqCoefficient = 0.0;
      // right!


      if( couplesAsVectorNotAxial )
      {
        minToLkConst->addToCoefficient( currentMinToLkConstCoefficient );
        minToLkLog->addToCoefficient( currentMinToLkLogCoefficient );
        minToLkLin->addToCoefficient( currentMinToLkLinCoefficient );

        lkToHkConst->addToCoefficient( currentLkToHkConstCoefficient );
        lkToHkLin->addToCoefficient( currentLkToHkLinCoefficient );

        hkToMaxConst->addToCoefficient( currentHkToMaxConstCoefficient );
        hkToMaxLog->addToCoefficient( currentHkToMaxLogCoefficient );
        hkToMaxLin->addToCoefficient( currentHkToMaxLinCoefficient );
      }
      else
      {
        minToLkConst->addToCoefficient( -currentMinToLkConstCoefficient );
        minToLkLog->addToCoefficient( -currentMinToLkLogCoefficient );
        minToLkLin->addToCoefficient( -currentMinToLkLinCoefficient );

        lkToHkConst->addToCoefficient( -currentLkToHkConstCoefficient );
        lkToHkLin->addToCoefficient( -currentLkToHkLinCoefficient );

        hkToMaxConst->addToCoefficient( -currentHkToMaxConstCoefficient );
        hkToMaxLog->addToCoefficient( -currentHkToMaxLogCoefficient );
        hkToMaxLin->addToCoefficient( -currentHkToMaxLinCoefficient );
      }

      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "VvSqMinusAaSqAntiDistribution:"
      << std::endl << "MINtoLK: "
      << std::endl << currentMinToLkInvCoefficient << "/E"
      << std::endl << currentMinToLkConstCoefficient
      << std::endl << currentMinToLkLogCoefficient << "L"
      << std::endl << currentMinToLkLinCoefficient << "E"
      << std::endl << currentMinToLkLinLogCoefficient << "EL"
      << std::endl << currentMinToLkSqCoefficient << "E^2"
      << std::endl << "LKtoHK: "
      << std::endl << currentLkToHkConstCoefficient
      << std::endl << currentLkToHkLinCoefficient << "E"
      << std::endl << currentLkToHkSqCoefficient << "E^2"
      << std::endl << "HKtoMAX: "
      << std::endl << currentHkToMaxInvCoefficient << "/E"
      << std::endl << currentHkToMaxConstCoefficient
      << std::endl << currentHkToMaxLogCoefficient << "L"
      << std::endl << currentHkToMaxLinCoefficient << "E"
      << std::endl << currentHkToMaxLinLogCoefficient << "EL"
      << std::endl << currentHkToMaxSqCoefficient << "E^2";**/
    }

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: Z_handed_muon:calculateCoefficients() produced "
    << std::endl << "firstMass = " << firstMass
    << std::endl << "secondMass = " << secondMass
    << std::endl << "cosinesLimitedByEnergy = " << cosinesLimitedByEnergy
    << std::endl << "(true = " << true << ")"
    << std::endl << "minimumEnergy = " << minimumEnergy
    << std::endl << "maximumEnergy = " << maximumEnergy
    << std::endl << "minToLkInv = " << minToLkInv->getCoefficient()
    << std::endl << "minToLkConst = " << minToLkConst->getCoefficient()
    << std::endl << "minToLkLog = " << minToLkLog->getCoefficient()
    << std::endl << "minToLkLin = " << minToLkLin->getCoefficient()
    << std::endl << "minToLkLinlog = " << minToLkLinlog->getCoefficient()
    << std::endl << "minToLkSq = " << minToLkSq->getCoefficient()
    << std::endl << "lkToHkConst = " << lkToHkConst->getCoefficient()
    << std::endl << "lkToHkLin = " << lkToHkLin->getCoefficient()
    << std::endl << "lkToHkSq = " << lkToHkSq->getCoefficient()
    << std::endl << "hkToMaxInv = " << hkToMaxInv->getCoefficient()
    << std::endl << "hkToMaxConst = " << hkToMaxConst->getCoefficient()
    << std::endl << "hkToMaxLog = " << hkToMaxLog->getCoefficient()
    << std::endl << "hkToMaxLin = " << hkToMaxLin->getCoefficient()
    << std::endl << "hkToMaxLinlog = " << hkToMaxLinlog->getCoefficient()
    << std::endl << "hkToMaxSq = " << hkToMaxSq->getCoefficient();
    std::cout << std::endl;**/
  }



  wMinusHandedMuon::wMinusHandedMuon( readierForNewPoint* const readierPointer,
                                   CppSLHA::CppSLHA0 const* const spectrumData,
                     CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                     CppSLHA::particle_property_set const* const thirdParticle,
                    CppSLHA::particle_property_set const* const fourthParticle,
                                       bool const jetIsLeftHanded ) :
    leptonEnergyDistribution( readierPointer,
                              spectrumData,
                              firstParticle,
                              effectiveSquarkMass,
                              secondParticle,
                              thirdParticle,
                              fourthParticle ),
    jetIsLeftHanded( jetIsLeftHanded ),
    minToLkInv( minToLkSegment.addTerm( -1,
                                        0 ) ),
    minToLkConst( minToLkSegment.addTerm( 0,
                                          0 ) ),
    minToLkLog( minToLkSegment.addTerm( 0,
                                        1 ) ),
    minToLkLin( minToLkSegment.addTerm( 1,
                                        0 ) ),
    minToLkLinlog( minToLkSegment.addTerm( 1,
                                           1 ) ),
    minToLkSq( minToLkSegment.addTerm( 2,
                                       0 ) ),
    lkToHkConst( lkToHkSegment.addTerm( 0,
                                        0 ) ),
    lkToHkLin( lkToHkSegment.addTerm( 1,
                                      0 ) ),
    lkToHkSq( lkToHkSegment.addTerm( 2,
                                     0 ) ),
    hkToMaxInv( hkToMaxSegment.addTerm( -1,
                                        0 ) ),
    hkToMaxConst( hkToMaxSegment.addTerm( 0,
                                          0 ) ),
    hkToMaxLog( hkToMaxSegment.addTerm( 0,
                                        1 ) ),
    hkToMaxLin( hkToMaxSegment.addTerm( 1,
                                        0 ) ),
    hkToMaxLinlog( hkToMaxSegment.addTerm( 1,
                                           1 ) ),
    hkToMaxSq( hkToMaxSegment.addTerm( 2,
                                       0 ) )
  {
    NMIX = spectrum->inspect_BLOCK( "NMIX" );
    UMIX = spectrum->inspect_BLOCK( "UMIX" );
    VMIX = spectrum->inspect_BLOCK( "VMIX" );
    if( CppSLHA::PDG_code::chargino_one == secondParticle->get_PDG_code() )
    {
      whichChargino = 1;
    }
    else if( CppSLHA::PDG_code::chargino_two
             == secondParticle->get_PDG_code() )
    {
      whichChargino = 2;
    }
    else
    {
      whichChargino = -1;
      // this will cause UMIX or VMIX to return NaN, which should show the
      // result up as obviously wrong.
      std::cout
      << std::endl
      << "LHC-FASER::error! W_minus_plus_up_L_type_muon expected to be given"
      << " a chargino as its 2nd particle, but did not receive a particle"
      << " with code +-" << CppSLHA::PDG_code::chargino_one << " or +-"
      << CppSLHA::PDG_code::chargino_two << ", rather, it received "
      << secondParticle->get_name() << " (code: "
      << secondParticle->get_PDG_code() << ")";
      std::cout << std::endl;
    }
    if( CppSLHA::PDG_code::neutralino_one == fourthParticle->get_PDG_code() )
    {
      whichNeutralino = 1;
    }
    else if( CppSLHA::PDG_code::neutralino_two
             == fourthParticle->get_PDG_code() )
    {
      whichNeutralino = 2;
    }
    else if( CppSLHA::PDG_code::neutralino_three
             == fourthParticle->get_PDG_code() )
    {
      whichNeutralino = 3;
    }
    else if( CppSLHA::PDG_code::neutralino_four
             == fourthParticle->get_PDG_code() )
    {
      whichNeutralino = 4;
    }
    else
    {
      whichNeutralino = -1;
      // this will cause NMIX to return NaN, which should show the result up
      // as obviously wrong.
      std::cout
      << std::endl
      << "LHC-FASER::error! W_minus_plus_up_L_type_muon expected to be given"
      << " a neutralino as its 4th particle, but did not receive a particle"
      << " with code " << CppSLHA::PDG_code::neutralino_one << ", "
      << CppSLHA::PDG_code::neutralino_two << ", "
      << CppSLHA::PDG_code::neutralino_three << ", or "
      << CppSLHA::PDG_code::neutralino_four << ", rather, it received "
      << fourthParticle->get_name() << " (code: "
      << fourthParticle->get_PDG_code() << ")";
      std::cout << std::endl;
    }
    segments.push_back( &minToLkSegment );
    segments.push_back( &lkToHkSegment );
    segments.push_back( &hkToMaxSegment );
  }

  wMinusHandedMuon::~wMinusHandedMuon()
  {
    // does nothing.
  }


  void
  wMinusHandedMuon::calculateCoefficients()
  {
    /* it doesn't matter if the neutralino mass is negative, since only the
     * absolute square appears, & the couplings are only affected by it as an
     * overall phase.  I (BOL) assume no CP violation, so the neutralino mixing
     * matrices are strictly real.
     */
    Vv = ( NMIX->get_entry( whichNeutralino,
                            2 ) * ( UMIX->get_entry( whichChargino,
                                                     1 )
                                    + VMIX->get_entry( whichChargino,
                                                       1 ) )
           + M_SQRT1_2 * ( NMIX->get_entry( whichNeutralino,
                                            3 )
                           * UMIX->get_entry( whichChargino,
                                              2 )
                           - NMIX->get_entry( whichNeutralino,
                                              4 )
                             * VMIX->get_entry( whichChargino,
                                                2 ) ) );
    Aa = ( NMIX->get_entry( whichNeutralino,
                            2 ) * ( VMIX->get_entry( whichChargino,
                                                     1 )
                                    - UMIX->get_entry( whichChargino,
                                                       1 ) )
           - M_SQRT1_2 * ( NMIX->get_entry( whichNeutralino,
                                            3 )
                           * UMIX->get_entry( whichChargino,
                                              2 )
                           + NMIX->get_entry( whichNeutralino,
                                              4 )
                             * VMIX->get_entry( whichChargino,
                                                2 ) ) );
    if( 0.0 > fourthParticle->get_mass() )
      // if the relative phases end up such that the product of couplings is
      // purely imaginary...
    {
      // we swap the couplings (borrowing currentLkToHkLinCoefficient for a
      // moment):
      currentLkToHkLinCoefficient = Vv;
      Vv = Aa;
      Aa = currentLkToHkLinCoefficient;
    }
    Vvsq = ( Vv * Vv );
    Aasq = ( Aa * Aa );
    VvSqMinusAaSq = ( Vvsq - Aasq );

    if( jetIsLeftHanded )
    {
      MinusTwiceVvJjAa = ( -2.0 * Vv * Aa );
      JjVvPlusJjAaAllSq = ( Vvsq - MinusTwiceVvJjAa + Aasq );
      JjVvSqMinusAaSq = VvSqMinusAaSq;
      JjVvSqPlusAaSq = ( Vvsq + Aasq );
    }
    else
    {
      MinusTwiceVvJjAa = ( 2.0 * Vv * Aa );
      JjVvPlusJjAaAllSq = ( MinusTwiceVvJjAa - Vvsq - Aasq );
      JjVvSqMinusAaSq = -VvSqMinusAaSq;
      JjVvSqPlusAaSq = ( -Vvsq - Aasq );
    }
    VvPlusJjAaAllSq = ( Vvsq - MinusTwiceVvJjAa + Aasq );

    mQ = firstMass;
    mQsq = ( mQ * mQ );
    mC = secondMass;
    mCsq = ( mC * mC );
    mVsq = ( thirdMass * thirdMass );
    mCmX = ( mC * fourthMass );
    mXsq = ( fourthMass * fourthMass );

    EV = ( ( mCsq + mVsq - mXsq ) / ( 2.0 * mC ) );
    ElMin = ( 0.5 * EV * ( 1.0 - sqrt( 1 - ( mVsq / ( EV * EV ) ) ) ) );
    minimumEnergy = ( ElMin * ( mC / mQ ) );
    Elk = ( ElMin * ( mQ / mC ) );
    Ehk = ( mVsq / ( 4.0 * Elk ) );
    maximumEnergy = ( mVsq / ( 4.0 * minimumEnergy ) );
    if( Ehk > Elk )
      /* if we're in the regime where the kinematics mean that there is no
       * range of energies which can be reached for all squark-rest-frame polar
       * muon momentum angles...
       */
    {
      // we swap the energies (borrowing currentLkToHkLinCoefficient for a
      // moment) & note that we had to swap them:
      currentLkToHkLinCoefficient = Elk;
      Elk = Ehk;
      Ehk = currentLkToHkLinCoefficient;
      cosinesLimitedByEnergy = true;
    }
    else
    {
      cosinesLimitedByEnergy = false;
    }
    mQCsqDiff = ( mQsq - mCsq );
    mQCsqDiffsq = ( mQCsqDiff * mQCsqDiff );
    mQCsqSum = ( mQsq + mCsq );
    EightmCEVElMin = ( 8.0 * mC * EV * ElMin );
    lnElMin = log( ElMin );
    ElMinsq = ( ElMin * ElMin );
    lnmQC = log( ( mQ / mC ) );
    mQcu = ( mQ * mQsq );
    mVqu = ( mVsq * mVsq );
    lnVsqOverFourElMinsq = log( ( mVsq / ( 4.0 * ElMinsq ) ) );
    EightmCXElMinmQsq = ( 8.0 * mC * fourthMass * ElMin * mQsq );
    lnEmin = log( minimumEnergy );
    lnElk = log( Elk );
    lnEhk = log( Ehk );
    lnEmax = log( maximumEnergy );

    minToLkSegment.setSegmentRange( minimumEnergy,
                                    Elk );
    lkToHkSegment.setSegmentRange( Elk,
                                    Ehk );
    hkToMaxSegment.setSegmentRange( Ehk,
                                    maximumEnergy );

    /* symmetric in lepton-antilepton =>
     *     ( Vv + jL * Aa )^2 * VvPlusJjAaAllSqSymDistribution
     *     - 2 * Jj * Vv * Aa * TwiceVvAaSymDistribution
     *     + ( Vv^( 2 ) - Aa^( 2 ) ) * VvSqMinusAaSqSymDistribution
     * antisymmetric in lepton-antilepton =>
     *     ( Vv + jL * Aa )^2 * VvPlusJjAaAllSqAntiDistribution
     *     + ( Vv^( 2 ) + Aa^( 2 ) ) * VvSqPlusAaSqAntiDistribution
     *     + ( Vv^( 2 ) - Aa^( 2 ) ) * VvSqMinusAaSqAntiDistribution
     * summing over handedness =>
     * symmetric
     * lepton same handedness as jet =>
     * symmetric + antisymmetric
     * lepton opposite handedness to jet =>
     * symmetric - antisymmetric
     */


    /* ends up in the order:
     * low (MINtoLK),
     * over (LKtoHK, true == cosinesLimitedByEnergy),
     * under (LKtoHK, false == cosinesLimitedByEnergy),
     * high (HKtoMAX)
     */


    // VvPlusJjAaAllSqSymDistribution:

    // MIN to LK:
    currentMinToLkConstCoefficient
    = ( ElMin * mQsq * mQCsqDiff
        * ( EightmCEVElMin * ( lnElMin - lnmQC ) * mVsq
            - 4.0 * ElMin * mC * ( EightmCEVElMin - 4.0 * ElMinsq * mC ) ) );
    // right!
    currentMinToLkLogCoefficient
    = ( -EightmCEVElMin * ElMin * mQsq * mQCsqDiff * mVsq );
    // right!
    currentMinToLkLinCoefficient
    = ( 4.0 * EightmCEVElMin * ElMin * mQcu * mQCsqDiff );
    // right!
    currentMinToLkSqCoefficient
    = ( -16.0 * ElMinsq * mQsq * mQsq * mQCsqDiff );
    // right!

    // LK to HK:
    if( cosinesLimitedByEnergy )
    {
      currentLkToHkConstCoefficient
      = ( mQsq * mQCsqDiff
          * ( EightmCEVElMin
              * ( mC * ( mVsq - 4.0 * ElMinsq )
                  - ElMin * lnVsqOverFourElMinsq * mVsq )
              + mCsq * ( 16.0 * ElMinsq * ElMinsq - mVsq * mVsq ) ) );
      // dunno...
      currentLkToHkLinCoefficient = 0.0;
      // dunno...
      currentLkToHkSqCoefficient = 0.0;
      // dunno...
    }
    else
    {
      currentLkToHkConstCoefficient
      = ( -2.0 * EightmCEVElMin * ElMin * lnmQC * mQsq * mQCsqDiff * mVsq );
      // right!
      currentLkToHkLinCoefficient
      = ( 4.0 * EightmCEVElMin * ElMin * mQ * mQCsqDiff * mQCsqDiff );
      // right!
      currentLkToHkSqCoefficient
      = ( -16.0 * ElMinsq * mQCsqDiff * mQCsqDiff * mQCsqSum );
      // right!
    }

    // HK to MAX:
    currentHkToMaxConstCoefficient
    = ( mQsq * mQCsqDiff * mVsq
        * ( EightmCEVElMin
            * ( mC - ElMin * ( lnmQC + lnVsqOverFourElMinsq + lnElMin ) )
            - mCsq * mVsq ) );
    // right!
    currentHkToMaxLogCoefficient
    = ( EightmCEVElMin * ElMin * mQsq * mQCsqDiff * mVsq );
    // right!
    currentHkToMaxLinCoefficient
    = ( -4.0 * EightmCEVElMin * ElMin * mCsq * mQ * mQCsqDiff );
    // right!
    currentHkToMaxSqCoefficient = ( 16.0 * ElMinsq * mCsq * mCsq * mQCsqDiff );
    // right!

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "VvPlusJjAaAllSq = " << VvPlusJjAaAllSq;
    std::cout << std::endl;/**/

    minToLkInv->setCoefficient( 0.0 );
    minToLkConst->setCoefficient( VvPlusJjAaAllSq
                                  * currentMinToLkConstCoefficient );
    minToLkLog->setCoefficient( VvPlusJjAaAllSq
                                * currentMinToLkLogCoefficient );
    minToLkLin->setCoefficient( VvPlusJjAaAllSq
                                * currentMinToLkLinCoefficient );
    minToLkLinlog->setCoefficient( 0.0 );
    minToLkSq->setCoefficient( VvPlusJjAaAllSq
                               * currentMinToLkSqCoefficient );

    lkToHkConst->setCoefficient( VvPlusJjAaAllSq
                                 * currentLkToHkConstCoefficient );
    lkToHkLin->setCoefficient( VvPlusJjAaAllSq
                               * currentLkToHkLinCoefficient );
    lkToHkSq->setCoefficient( VvPlusJjAaAllSq
                              * currentLkToHkSqCoefficient );

    hkToMaxInv->setCoefficient( 0.0 );
    hkToMaxConst->setCoefficient( VvPlusJjAaAllSq
                                  * currentHkToMaxConstCoefficient );
    hkToMaxLog->setCoefficient( VvPlusJjAaAllSq
                                * currentHkToMaxLogCoefficient );
    hkToMaxLin->setCoefficient( VvPlusJjAaAllSq
                                * currentHkToMaxLinCoefficient );
    hkToMaxLinlog->setCoefficient( 0.0 );
    hkToMaxSq->setCoefficient( VvPlusJjAaAllSq
                               * currentHkToMaxSqCoefficient );

    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "VvPlusJjAaAllSqSymDistribution:"
    << std::endl << "MINtoLK: "
    << std::endl << minToLkInv->getCoefficient() << "/E"
    << std::endl << minToLkConst->getCoefficient()
    << std::endl << minToLkLog->getCoefficient() << "L"
    << std::endl << minToLkLin->getCoefficient() << "E"
    << std::endl << minToLkLinlog->getCoefficient() << "EL"
    << std::endl << minToLkSq->getCoefficient() << "E^2"
    << std::endl << "LKtoHK: "
    << std::endl << lkToHkConst->getCoefficient()
    << std::endl << lkToHkLin->getCoefficient() << "E"
    << std::endl << lkToHkSq->getCoefficient() << "E^2"
    << std::endl << "HKtoMAX: "
    << std::endl << hkToMaxInv->getCoefficient() << "/E"
    << std::endl << hkToMaxConst->getCoefficient()
    << std::endl << hkToMaxLog->getCoefficient() << "L"
    << std::endl << hkToMaxLin->getCoefficient() << "E"
    << std::endl << hkToMaxLinlog->getCoefficient() << "EL"
    << std::endl << hkToMaxSq->getCoefficient() << "E^2";
    std::cout << std::endl;**/



    // TwiceVvAaSymDistribution:

    // MIN to LK:
    currentMinToLkInvCoefficient
    = ( -4.0 * ElMinsq * mCsq * mQcu * mVsq * mVsq );
    // right!
    currentMinToLkConstCoefficient
    = ( ElMin * mQsq
        * ( 2.0 * EightmCEVElMin
                * ( mVsq * ( mQsq * ( lnElMin - lnmQC ) - mCsq )
                    - 4.0 * ElMin * mC * mQsq )
            + 4.0 * mC
              * ( 2.0 * ElMin * mC
                  * ( 4.0 * ElMinsq * mQsq
                      + mVsq * ( mQCsqSum * ( lnElMin - lnmQC )
                                 - 2.0 * mCsq ) )
                  + mQCsqSum * mVsq * mVsq ) ) );
    // right!
    currentMinToLkLogCoefficient
    = ( -2.0 * ElMin * mQsq * mVsq
        * ( EightmCEVElMin * mQsq + 4.0 * ElMin * mCsq * mQCsqSum ) );
    // right!
    currentMinToLkLinCoefficient
    = ( 2.0 * mQcu
        * ( 4.0 * EightmCEVElMin * ElMin
            * ( ( lnElMin - lnmQC ) * mCsq + mQsq )
            + EightmCEVElMin * mC * mVsq
            - 2.0 * mCsq * ( 4.0 * ElMin * mC * ( 4.0 * ElMinsq - mVsq )
                             + mVsq * mVsq) ) );
    // right!
    currentMinToLkLinLogCoefficient
    = ( -8.0 * EightmCEVElMin * ElMin * mCsq * mQcu );
    // right!
    currentMinToLkSqCoefficient
    = ( -32.0 * ElMinsq * mQsq * mQsq * ( mQCsqDiff - mCsq ) );
    // right!

    // LK to HK:
    if( cosinesLimitedByEnergy )
    {
      currentLkToHkConstCoefficient
      = ( 2.0 * mQsq
          * ( EightmCEVElMin * mQsq
              * ( mC * ( mVsq - 4.0 * ElMinsq )
                  - ElMin * mVsq * lnVsqOverFourElMinsq )
              + mC * ( 16.0 * ElMinsq * ElMinsq * mC * mQsq
                       - 4.0 * ElMinsq * mQCsqSum * mVsq
                         * ( 2.0 * ElMin + lnVsqOverFourElMinsq * mC )
                       + ( 2.0 * ElMin * mQCsqSum
                           - mC * mQsq ) * mVsq * mVsq ) ) );
      // dunno...
      currentLkToHkLinCoefficient
      = ( 2.0 * mC * mQcu
          * ( EightmCEVElMin
              * ( mVsq
                  - 4.0 * ElMin * ( lnVsqOverFourElMinsq * mC + ElMin ) )
              + 2.0 * mC * ( 8.0 * ElMin
                             * ( 2.0 * ElMinsq * ( ElMin - 2.0 * mC )
                                 + mC * mVsq )
                             - mVsq * mVsq ) ) );
      // dunno...
      currentLkToHkSqCoefficient = 0.0;
      // dunno...
    }
    else
    {
      currentLkToHkConstCoefficient
      = ( 2.0 * ElMin * mQsq * mVsq
          * ( mQCsqDiff
              * ( EightmCEVElMin + 8.0 * ElMin * mCsq )
              - 2.0 * lnmQC * ( 4.0 * mQCsqSum * mCsq * ElMin
                                + mQsq * EightmCEVElMin ) ) );
      // right!
      currentLkToHkLinCoefficient
      = ( 8.0 * EightmCEVElMin * ElMin * mQcu
          * ( mQCsqDiff - 2.0 * lnmQC * mCsq ) );
      // right!
      currentLkToHkSqCoefficient
      = ( -32.0 * ElMinsq * mQsq * mQCsqDiff * mQCsqDiff );
      // right!
    }

    // HK to MAX:
    currentHkToMaxInvCoefficient
    = ( 4.0 * ElMinsq * mCsq * mQcu * mVsq * mVsq );
    // right!
    currentHkToMaxConstCoefficient
    = ( 2.0 * mQsq * mVsq
        * ( EightmCEVElMin * ElMin * mQsq
            * ( 1.0 - lnmQC - lnVsqOverFourElMinsq - lnElMin )
            - 4.0 * ElMinsq * mCsq
              * ( ( lnmQC + lnVsqOverFourElMinsq + lnElMin ) * mQCsqSum
                  - 2.0 * mQsq )
            - 8.0 * ElMin * ElMinsq * mC * mQCsqSum
            + mC * mQsq * ( EightmCEVElMin - mC * mVsq ) ) );
    // right!
    currentHkToMaxLogCoefficient
    = ( 2.0 * ElMin * mQsq * mVsq
        * ( EightmCEVElMin * mQsq + 4.0 * ElMin * mCsq * mQCsqSum ) );
    // right!
    currentHkToMaxLinCoefficient
    = ( 8.0  * ElMin * mC * mQcu
        * ( 2.0 * mC * ( 4.0 * ElMinsq * ( ElMin - mC ) + mC * mVsq )
            - EightmCEVElMin
            * ( ElMin
                + ( 1.0 + lnmQC + lnElMin + lnVsqOverFourElMinsq ) * mC ) ) );
    // right!
    currentHkToMaxLinLogCoefficient
    = ( 8.0 * EightmCEVElMin * ElMin * mCsq * mQcu );
    // right!
    currentHkToMaxSqCoefficient = ( -32.0 * ElMinsq * mCsq * mCsq * mQsq );
    // right!

    //debugging:
    MinusTwiceVvJjAa = 0.0;
    minToLkInv->addToCoefficient( MinusTwiceVvJjAa
                                  * currentMinToLkInvCoefficient );
    minToLkConst->addToCoefficient( MinusTwiceVvJjAa
                                    * currentMinToLkConstCoefficient );
    minToLkLog->addToCoefficient( MinusTwiceVvJjAa
                                  * currentMinToLkLogCoefficient );
    minToLkLin->addToCoefficient( MinusTwiceVvJjAa
                                  * currentMinToLkLinCoefficient );
    minToLkLinlog->addToCoefficient( MinusTwiceVvJjAa
                                     * currentMinToLkLinLogCoefficient );
    minToLkSq->addToCoefficient( MinusTwiceVvJjAa
                                 * currentMinToLkSqCoefficient );

    lkToHkConst->addToCoefficient( MinusTwiceVvJjAa
                                   * currentLkToHkConstCoefficient );
    lkToHkLin->addToCoefficient( MinusTwiceVvJjAa
                                 * currentLkToHkLinCoefficient );
    lkToHkSq->addToCoefficient( MinusTwiceVvJjAa
                                * currentLkToHkSqCoefficient );

    hkToMaxInv->addToCoefficient( MinusTwiceVvJjAa
                                  * currentHkToMaxInvCoefficient );
    hkToMaxConst->addToCoefficient( MinusTwiceVvJjAa
                                    * currentHkToMaxConstCoefficient );
    hkToMaxLog->addToCoefficient( MinusTwiceVvJjAa
                                  * currentHkToMaxLogCoefficient );
    hkToMaxLin->addToCoefficient( MinusTwiceVvJjAa
                                  * currentHkToMaxLinCoefficient );
    hkToMaxLinlog->addToCoefficient( MinusTwiceVvJjAa
                                     * currentHkToMaxLinLogCoefficient );
    hkToMaxSq->addToCoefficient( MinusTwiceVvJjAa
                                 * currentHkToMaxSqCoefficient );

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "TwiceVvAaSymDistribution:"
    << std::endl << "MINtoLK: "
    << std::endl << currentMinToLkInvCoefficient << "/E"
    << std::endl << currentMinToLkConstCoefficient
    << std::endl << currentMinToLkLogCoefficient << "L"
    << std::endl << currentMinToLkLinCoefficient << "E"
    << std::endl << currentMinToLkLinLogCoefficient << "EL"
    << std::endl << currentMinToLkSqCoefficient << "E^2"
    << std::endl << "LKtoHK: "
    << std::endl << currentLkToHkConstCoefficient
    << std::endl << currentLkToHkLinCoefficient << "E"
    << std::endl << currentLkToHkSqCoefficient << "E^2"
    << std::endl << "HKtoMAX: "
    << std::endl << currentHkToMaxInvCoefficient << "/E"
    << std::endl << currentHkToMaxConstCoefficient
    << std::endl << currentHkToMaxLogCoefficient << "L"
    << std::endl << currentHkToMaxLinCoefficient << "E"
    << std::endl << currentHkToMaxLinLogCoefficient << "EL"
    << std::endl << currentHkToMaxSqCoefficient << "E^2";/**/


    // VvSqMinusAaSqSymDistribution:

    // MIN to LK:
    // debugging:
    //currentMinToLkInvCoefficient = 0.0;
    // right!
    currentMinToLkConstCoefficient
    = ( EightmCXElMinmQsq * ElMin * ( lnmQC - lnElMin ) * mQCsqDiff * mVsq );
    // right!
    currentMinToLkLogCoefficient
    = ( EightmCXElMinmQsq * ElMin * mQCsqDiff * mVsq );
    // right!
    // debugging:
    //currentMinToLkLinCoefficient = 0.0;
    // right!
    // debugging:
    //currentMinToLkLinLogCoefficient = 0.0;
    // right!
    // debugging:
    //currentMinToLkSqCoefficient = 0.0;
    // right!

    // LK to HK:
    if( cosinesLimitedByEnergy )
    {
      currentLkToHkConstCoefficient = ( EightmCXElMinmQsq * ElMin * mQCsqDiff
                                        * mVsq * lnVsqOverFourElMinsq );
      // dunno...
    }
    else
    {
      currentLkToHkConstCoefficient
      = ( 2.0 * EightmCXElMinmQsq * ElMin * lnmQC * mQCsqDiff * mVsq );
      // right!
    }
    // debugging:
    //currentLkToHkLinCoefficient = 0.0;
    // right!
    // debugging:
    //currentLkToHkSqCoefficient = 0.0;
    // right!

    // HK to MAX:
    // debugging:
    //currentHkToMaxInvCoefficient = 0.0;
    // right!
    currentHkToMaxConstCoefficient
    = ( EightmCXElMinmQsq * ElMin * mQCsqDiff * mVsq
        * ( lnmQC + lnVsqOverFourElMinsq + lnElMin ) );
    // right!
    currentHkToMaxLogCoefficient = -currentMinToLkLogCoefficient;
    // right!
    // debugging:
    //currentHkToMaxLinCoefficient = 0.0;
    // right!
    // debugging:
    //currentHkToMaxLinLogCoefficient = 0.0;
    // right!
    // debugging:
    //currentHkToMaxSqCoefficient = 0.0;
    // right!

    //debugging:
    VvSqMinusAaSq = 0.0;
    minToLkConst->addToCoefficient( VvSqMinusAaSq
                                    * currentMinToLkConstCoefficient );
    minToLkLog->addToCoefficient( VvSqMinusAaSq
                                  * currentMinToLkLogCoefficient );

    lkToHkConst->addToCoefficient( VvSqMinusAaSq
                                   * currentLkToHkConstCoefficient );

    hkToMaxConst->addToCoefficient( VvSqMinusAaSq
                                    * currentHkToMaxConstCoefficient );
    hkToMaxLog->addToCoefficient( VvSqMinusAaSq
                                  * currentHkToMaxLogCoefficient );

    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "VvSqMinusAaSqSymDistribution:"
    << std::endl << "MINtoLK: "
    << std::endl << currentMinToLkInvCoefficient << "/E"
    << std::endl << currentMinToLkConstCoefficient
    << std::endl << currentMinToLkLogCoefficient << "L"
    << std::endl << currentMinToLkLinCoefficient << "E"
    << std::endl << currentMinToLkLinLogCoefficient << "EL"
    << std::endl << currentMinToLkSqCoefficient << "E^2"
    << std::endl << "LKtoHK: "
    << std::endl << currentLkToHkConstCoefficient
    << std::endl << currentLkToHkLinCoefficient << "E"
    << std::endl << currentLkToHkSqCoefficient << "E^2"
    << std::endl << "HKtoMAX: "
    << std::endl << currentHkToMaxInvCoefficient << "/E"
    << std::endl << currentHkToMaxConstCoefficient
    << std::endl << currentHkToMaxLogCoefficient << "L"
    << std::endl << currentHkToMaxLinCoefficient << "E"
    << std::endl << currentHkToMaxLinLogCoefficient << "EL"
    << std::endl << currentHkToMaxSqCoefficient << "E^2";**/


    // VvPlusJjAaAllSqAntiDistribution:

    // MIN to LK:
    // debugging:
    //currentMinToLkInvCoefficient = 0.0;
    // right!
    currentMinToLkConstCoefficient
    =  ( ElMin * ( EightmCEVElMin * ( lnElMin - lnmQC )
                   - 16.0 * ElMinsq * mC ) * mQsq * mQCsqDiff * mVsq );
    // right!
    currentMinToLkLogCoefficient
    = ( -EightmCEVElMin * ElMin * mQsq * mQCsqDiff * mVsq );
    // right!
    currentMinToLkLinCoefficient
    = ( 16.0 * ElMinsq * mQcu * mQCsqDiff * mVsq );
    // right!
    // debugging:
    //currentMinToLkLinLogCoefficient = 0.0;
    // right!
    // debugging:
    //currentMinToLkSqCoefficient = 0.0;
    // right!

    // LK to HK:
    if( cosinesLimitedByEnergy )
    {
      currentLkToHkConstCoefficient
      = ( ElMin * mQsq * mQCsqDiff * mVsq
          * ( 4.0 * mC * ( mVsq - 4.0 * ElMinsq )
          - EightmCEVElMin * lnVsqOverFourElMinsq ) );
      // dunno...
      currentLkToHkLinCoefficient = 0.0;
      // dunno...
    }
    else
    {
      currentLkToHkConstCoefficient
      = ( -2.0 * EightmCEVElMin * ElMin * lnmQC * mQsq * mQCsqDiff * mVsq );
      // right!
      currentLkToHkLinCoefficient
      = ( 16.0 * ElMinsq * mQ * mQCsqDiff * mQCsqDiff * mVsq );
      // right!
    }
    // debugging:
    //currentLkToHkSqCoefficient = 0.0;

    // HK to MAX:
    //currentHkToMaxInvCoefficient = 0.0;
    // right!
    currentHkToMaxConstCoefficient
    = ( ElMin * mQsq * mQCsqDiff * mVsq
        * ( 4.0 * mC * mVsq - EightmCEVElMin
                              * ( lnmQC + lnVsqOverFourElMinsq + lnElMin ) ) );
    // right!
    currentHkToMaxLogCoefficient
    = ( EightmCEVElMin * ElMin * mQsq * mQCsqDiff * mVsq );
    // right!
    currentHkToMaxLinCoefficient
    = ( -16.0 * ElMinsq * mCsq * mQ * mQCsqDiff * mVsq );
    // right!
    // debugging:
    //currentHkToMaxLinLogCoefficient = 0.0;
    // right!
    // debugging:
    //currentHkToMaxSqCoefficient = 0.0;
    // right!

    //debugging:
    JjVvPlusJjAaAllSq = 0.0;
    minToLkConst->addToCoefficient( JjVvPlusJjAaAllSq
                                    * currentMinToLkConstCoefficient );
    minToLkLog->addToCoefficient( JjVvPlusJjAaAllSq
                                  * currentMinToLkLogCoefficient );
    minToLkLin->addToCoefficient( JjVvPlusJjAaAllSq
                                  * currentMinToLkLinCoefficient );

    lkToHkConst->addToCoefficient( JjVvPlusJjAaAllSq
                                   * currentLkToHkConstCoefficient );
    lkToHkLin->addToCoefficient( JjVvPlusJjAaAllSq
                                 * currentLkToHkLinCoefficient );

    hkToMaxConst->addToCoefficient( JjVvPlusJjAaAllSq
                                    * currentHkToMaxConstCoefficient );
    hkToMaxLog->addToCoefficient( JjVvPlusJjAaAllSq
                                  * currentHkToMaxLogCoefficient );
    hkToMaxLin->addToCoefficient( JjVvPlusJjAaAllSq
                                  * currentHkToMaxLinCoefficient );

    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "VvPlusJjAaAllSqAntiDistribution:"
    << std::endl << "MINtoLK: "
    << std::endl << currentMinToLkInvCoefficient << "/E"
    << std::endl << currentMinToLkConstCoefficient
    << std::endl << currentMinToLkLogCoefficient << "L"
    << std::endl << currentMinToLkLinCoefficient << "E"
    << std::endl << currentMinToLkLinLogCoefficient << "EL"
    << std::endl << currentMinToLkSqCoefficient << "E^2"
    << std::endl << "LKtoHK: "
    << std::endl << currentLkToHkConstCoefficient
    << std::endl << currentLkToHkLinCoefficient << "E"
    << std::endl << currentLkToHkSqCoefficient << "E^2"
    << std::endl << "HKtoMAX: "
    << std::endl << currentHkToMaxInvCoefficient << "/E"
    << std::endl << currentHkToMaxConstCoefficient
    << std::endl << currentHkToMaxLogCoefficient << "L"
    << std::endl << currentHkToMaxLinCoefficient << "E"
    << std::endl << currentHkToMaxLinLogCoefficient << "EL"
    << std::endl << currentHkToMaxSqCoefficient << "E^2";**/


    // VvSqPlusAaSqAntiDistribution:

    // MIN to LK:
    currentMinToLkInvCoefficient
    = ( 4.0 * ElMinsq * mCsq * mQcu * mVsq * mVsq );
    // right!
    currentMinToLkConstCoefficient
    = ( ElMin * mQsq * mVsq
        * ( 32.0 * ElMinsq * mC * mQsq
            + 2.0 * EightmCEVElMin * ( mCsq + ( lnmQC - lnElMin ) * mQsq )
            + 8.0 * ElMin * mCsq
                  * ( 2.0 * mCsq + ( lnmQC - lnElMin ) * mQCsqSum )
            - 4.0 * mC * mQCsqSum * mVsq ) );
    // right!
    currentMinToLkLogCoefficient
    = ( 2.0 * ElMin * mQsq * mVsq
            * ( EightmCEVElMin * mQsq + 4.0 * ElMin * mCsq * mQCsqSum ) );
    // right!
    currentMinToLkLinCoefficient
    = ( 2.0 * mQcu * mVsq
        * ( 2.0 * mCsq * mVsq
            - EightmCEVElMin * mC
            - 8.0 * ElMin
              * ( mC * mCsq + 2.0 * ElMin
                              * ( ( lnElMin - lnmQC ) * mCsq + mQsq ) ) ) );
    // right!
    currentMinToLkLinLogCoefficient = ( 32.0 * ElMinsq * mCsq * mQcu * mVsq );
    // right!
    // debugging:
    //currentMinToLkSqCoefficient = 0.0;
    // right!

    // LK to HK:
    if( cosinesLimitedByEnergy )
    {
      currentLkToHkConstCoefficient
      = ( 2.0 * ElMin * mQsq * mVsq
          * ( lnVsqOverFourElMinsq
              * ( EightmCEVElMin * mQsq + 4.0 * ElMin * mCsq * mQCsqSum )
              + 2.0 * mC * ( mCsq + 3.0 * mQsq )
                * ( 4.0 * ElMinsq - mVsq ) ) );
      // dunno...
      currentLkToHkLinCoefficient
      = ( 2.0 * mC * mQcu
          * ( 4.0 * ElMinsq
              * ( EightmCEVElMin
                  + 8.0 * ElMin * mC * ( mC - ElMin ) )
              - mVsq * ( EightmCEVElMin
                         + 8.0 * ElMin * mC
                           * ( mC - 2.0 * ElMin * lnVsqOverFourElMinsq ) )
              + 2.0 * mC * mVsq * mVsq ) );
      // dunno...
    }
    else
    {
      currentLkToHkConstCoefficient
      = ( 2.0 * ElMin * mQsq * mVsq
          * ( 2.0 * mQsq * lnmQC
              * ( EightmCEVElMin + 4.0 * ElMin * mCsq )
              - mQCsqDiff * ( 8.0 * ElMin * mCsq + EightmCEVElMin )
              + 8.0 * ElMin * mCsq * mCsq * lnmQC ) );
      // right!
      currentLkToHkLinCoefficient = ( 32.0 * ElMinsq * mQcu * mVsq
                                      * ( 2.0 * mCsq * lnmQC - mQCsqDiff ) );
      // right!
    }
    // debugging:
    //currentLkToHkSqCoefficient = 0.0;
    // right!

    // HK to MAX:
    currentHkToMaxInvCoefficient = -currentMinToLkInvCoefficient;
    // depends...
    currentHkToMaxConstCoefficient
    = ( 2.0 * ElMin * mQsq * mVsq
        * ( 4.0 * ElMin * mCsq
            * ( ( lnmQC + lnVsqOverFourElMinsq + lnElMin ) * mQCsqSum
                  - 2.0 * mQsq )
            + 8.0 * ElMinsq * mC * mQCsqSum
            + mQsq * ( EightmCEVElMin
                       * ( lnmQC + lnVsqOverFourElMinsq + lnElMin - 1.0 )
                       - 4.0 * mC * mVsq ) ) );
    // right!
    currentHkToMaxLogCoefficient
    = ( -2.0 * ElMin * mQsq * mVsq
        * ( EightmCEVElMin * mQsq + 4.0 * ElMin * mCsq * mQCsqSum ) );
    // right!
    currentHkToMaxLinCoefficient
    = ( 8.0 * ElMinsq * mC * mQcu
        * ( EightmCEVElMin + 8.0 * ElMin * mC * ( mC - ElMin )
            + 4.0 * mC * mVsq
              * ( 1.0 + lnmQC + lnVsqOverFourElMinsq + lnElMin ) ) );
    // right!
    currentHkToMaxLinLogCoefficient = -currentMinToLkLinLogCoefficient;
    // right!
    // debugging:
    //currentHkToMaxSqCoefficient = 0.0;
    // right!

    //debugging:
    JjVvSqPlusAaSq = 0.0;
    minToLkInv->addToCoefficient( JjVvSqPlusAaSq
                                  * currentMinToLkInvCoefficient );
    minToLkConst->addToCoefficient( JjVvSqPlusAaSq
                                    * currentMinToLkConstCoefficient );
    minToLkLog->addToCoefficient( JjVvSqPlusAaSq
                                  * currentMinToLkLogCoefficient );
    minToLkLin->addToCoefficient( JjVvSqPlusAaSq
                                  * currentMinToLkLinCoefficient );
    minToLkLinlog->addToCoefficient( JjVvSqPlusAaSq
                                     * currentMinToLkLinLogCoefficient );

    lkToHkConst->addToCoefficient( JjVvSqPlusAaSq
                                   * currentLkToHkConstCoefficient );
    lkToHkLin->addToCoefficient( JjVvSqPlusAaSq
                                 * currentLkToHkLinCoefficient );

    hkToMaxInv->addToCoefficient( JjVvSqPlusAaSq
                                  * currentHkToMaxInvCoefficient );
    hkToMaxConst->addToCoefficient( JjVvSqPlusAaSq
                                    * currentHkToMaxConstCoefficient );
    hkToMaxLog->addToCoefficient( JjVvSqPlusAaSq
                                  * currentHkToMaxLogCoefficient );
    hkToMaxLin->addToCoefficient( JjVvSqPlusAaSq
                                  * currentHkToMaxLinCoefficient );
    hkToMaxLinlog->addToCoefficient( JjVvSqPlusAaSq
                                     * currentHkToMaxLinLogCoefficient );

    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "VvSqPlusAaSqAntiDistribution:"
    << std::endl << "MINtoLK: "
    << std::endl << currentMinToLkInvCoefficient << "/E"
    << std::endl << currentMinToLkConstCoefficient
    << std::endl << currentMinToLkLogCoefficient << "L"
    << std::endl << currentMinToLkLinCoefficient << "E"
    << std::endl << currentMinToLkLinLogCoefficient << "EL"
    << std::endl << currentMinToLkSqCoefficient << "E^2"
    << std::endl << "LKtoHK: "
    << std::endl << currentLkToHkConstCoefficient
    << std::endl << currentLkToHkLinCoefficient << "E"
    << std::endl << currentLkToHkSqCoefficient << "E^2"
    << std::endl << "HKtoMAX: "
    << std::endl << currentHkToMaxInvCoefficient << "/E"
    << std::endl << currentHkToMaxConstCoefficient
    << std::endl << currentHkToMaxLogCoefficient << "L"
    << std::endl << currentHkToMaxLinCoefficient << "E"
    << std::endl << currentHkToMaxLinLogCoefficient << "EL"
    << std::endl << currentHkToMaxSqCoefficient << "E^2";**/


    // VvSqMinusAaSqAntiDistribution

    // MIN to LK:
    // debugging:
    //currentMinToLkInvCoefficient = 0.0;
    // right!
    currentMinToLkConstCoefficient
    = ( EightmCXElMinmQsq * ElMin * mVsq
        * ( ( 1.0 - lnElMin + lnmQC ) * mQCsqSum - mQCsqDiff ) );
    // right!
    currentMinToLkLogCoefficient
    = ( EightmCXElMinmQsq * ElMin * mQCsqSum * mVsq );
    // right!
    currentMinToLkLinCoefficient
    = ( -2.0 * EightmCXElMinmQsq * mC * mQ * mVsq );
    // right!
    // debugging:
    //currentMinToLkLinLogCoefficient = 0.0;
    // right!
    // debugging:
    //currentMinToLkSqCoefficient = 0.0;
    // right!

    // LK to HK:
    if( cosinesLimitedByEnergy )
    {
      currentLkToHkConstCoefficient = ( EightmCXElMinmQsq * ElMin * mQCsqSum
                                        * mVsq * lnVsqOverFourElMinsq );
      // dunno...
      currentLkToHkLinCoefficient
      = ( 2.0 * EightmCXElMinmQsq * mC * mQ * ( 4.0 * ElMinsq - mVsq ) );
      // dunno...
    }
    else
    {
      currentLkToHkConstCoefficient = ( 2.0 * EightmCXElMinmQsq * ElMin * mVsq
                                        * ( lnmQC * mQCsqSum - mQCsqDiff ) );
      // right!
      currentLkToHkLinCoefficient = 0.0;
      // right!
    }
    // debugging:
    //currentLkToHkSqCoefficient = 0.0;
    // right!

    // HK to MAX:
    // debugging:
    //currentHkToMaxInvCoefficient = 0.0;
    // right!
    currentHkToMaxConstCoefficient
    = ( EightmCXElMinmQsq * ElMin * mVsq
        * ( ( lnVsqOverFourElMinsq + lnElMin + lnmQC - 1.0 ) * mQCsqSum
            - mQCsqDiff ) );
    // right!
    currentHkToMaxLogCoefficient = -currentMinToLkLogCoefficient;
    // right!
    currentHkToMaxLinCoefficient
    = ( 8.0 * EightmCXElMinmQsq * mC * mQ * ElMinsq );
    // right!
    // debugging:
    //currentHkToMaxLinLogCoefficient = 0.0;
    // right!
    // debugging:
    //currentHkToMaxSqCoefficient = 0.0;
    // right!

    //debugging:
    JjVvSqMinusAaSq = 0.0;
    minToLkConst->addToCoefficient( JjVvSqMinusAaSq
                                    * currentMinToLkConstCoefficient );
    minToLkLog->addToCoefficient( JjVvSqMinusAaSq
                                  * currentMinToLkLogCoefficient );
    minToLkLin->addToCoefficient( JjVvSqMinusAaSq
                                  * currentMinToLkLinCoefficient );

    lkToHkConst->addToCoefficient( JjVvSqMinusAaSq
                                   * currentLkToHkConstCoefficient );
    lkToHkLin->addToCoefficient( JjVvSqMinusAaSq
                                 * currentLkToHkLinCoefficient );

    hkToMaxConst->addToCoefficient( JjVvSqMinusAaSq
                                    * currentHkToMaxConstCoefficient );
    hkToMaxLog->addToCoefficient( JjVvSqMinusAaSq
                                  * currentHkToMaxLogCoefficient );
    hkToMaxLin->addToCoefficient( JjVvSqMinusAaSq
                                  * currentHkToMaxLinCoefficient );

    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "VvSqMinusAaSqAntiDistribution:"
    << std::endl << "MINtoLK: "
    << std::endl << currentMinToLkInvCoefficient << "/E"
    << std::endl << currentMinToLkConstCoefficient
    << std::endl << currentMinToLkLogCoefficient << "L"
    << std::endl << currentMinToLkLinCoefficient << "E"
    << std::endl << currentMinToLkLinLogCoefficient << "EL"
    << std::endl << currentMinToLkSqCoefficient << "E^2"
    << std::endl << "LKtoHK: "
    << std::endl << currentLkToHkConstCoefficient
    << std::endl << currentLkToHkLinCoefficient << "E"
    << std::endl << currentLkToHkSqCoefficient << "E^2"
    << std::endl << "HKtoMAX: "
    << std::endl << currentHkToMaxInvCoefficient << "/E"
    << std::endl << currentHkToMaxConstCoefficient
    << std::endl << currentHkToMaxLogCoefficient << "L"
    << std::endl << currentHkToMaxLinCoefficient << "E"
    << std::endl << currentHkToMaxLinLogCoefficient << "EL"
    << std::endl << currentHkToMaxSqCoefficient << "E^2";**/

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: wMinusHandedMuon::calculateCoefficients() produced "
    << std::endl << "firstMass = " << firstMass
    << std::endl << "secondMass = " << secondMass
    << std::endl << "cosinesLimitedByEnergy = " << cosinesLimitedByEnergy
    << std::endl << "(true = " << true << ")"
    << std::endl << "minimumEnergy = " << minimumEnergy
    << std::endl << "maximumEnergy = " << maximumEnergy
    << std::endl << "minToLkInv = " << minToLkInv->getCoefficient()
    << std::endl << "minToLkConst = " << minToLkConst->getCoefficient()
    << std::endl << "minToLkLog = " << minToLkLog->getCoefficient()
    << std::endl << "minToLkLin = " << minToLkLin->getCoefficient()
    << std::endl << "minToLkLinlog = " << minToLkLinlog->getCoefficient()
    << std::endl << "minToLkSq = " << minToLkSq->getCoefficient()
    << std::endl << "lkToHkConst = " << lkToHkConst->getCoefficient()
    << std::endl << "lkToHkLin = " << lkToHkLin->getCoefficient()
    << std::endl << "lkToHkSq = " << lkToHkSq->getCoefficient()
    << std::endl << "hkToMaxInv = " << hkToMaxInv->getCoefficient()
    << std::endl << "hkToMaxConst = " << hkToMaxConst->getCoefficient()
    << std::endl << "hkToMaxLog = " << hkToMaxLog->getCoefficient()
    << std::endl << "hkToMaxLin = " << hkToMaxLin->getCoefficient()
    << std::endl << "hkToMaxLinlog = " << hkToMaxLinlog->getCoefficient()
    << std::endl << "hkToMaxSq = " << hkToMaxSq->getCoefficient();
    std::cout << std::endl;**/
  }



  vectorFromSquarkToMuon::vectorFromSquarkToMuon(
                                      readierForNewPoint* const readierPointer,
                                   CppSLHA::CppSLHA0 const* const spectrumData,
                     CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                  CppSLHA::particle_property_set const* const thirdParticle ) :
    leptonEnergyDistribution( readierPointer,
                              spectrumData,
                              firstParticle,
                              effectiveSquarkMass,
                              secondParticle,
                              thirdParticle,
                              NULL ),
    minToMaxConst( minToMaxSegment.addTerm( 0,
                                            0 ) ),
    minToMaxLin( minToMaxSegment.addTerm( 1,
                                          0 ) ),
    minToMaxSq( minToMaxSegment.addTerm( 2,
                                         0 ) )
  {
    segments.push_back( &minToMaxSegment );
  }

  vectorFromSquarkToMuon::~vectorFromSquarkToMuon()
  {
    // does nothing.
  }



  zDirectJet::zDirectJet( readierForNewPoint* const readierPointer,
                          CppSLHA::CppSLHA0 const* const spectrumData,
                     CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                     CppSLHA::particle_property_set const* const thirdParticle,
                 CppSLHA::particle_property_set const* const fourthParticle ) :
    zHandedMuon( readierPointer,
                 spectrumData,
                 firstParticle,
                 effectiveSquarkMass,
                 secondParticle,
                 thirdParticle,
                 fourthParticle,
                 true,
                 true )
  {
    // just an initialization of the base class object.
  }

  zDirectJet::~zDirectJet()
  {
    // does nothing.
  }



  wMinusDirectJet::wMinusDirectJet( readierForNewPoint* const readierPointer,
                                   CppSLHA::CppSLHA0 const* const spectrumData,
                     CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                     CppSLHA::particle_property_set const* const thirdParticle,
                 CppSLHA::particle_property_set const* const fourthParticle ) :
    zHandedMuon( readierPointer,
                 spectrumData,
                 firstParticle,
                 effectiveSquarkMass,
                 secondParticle,
                 thirdParticle,
                 fourthParticle,
                 true,
                 true )
  {
    // just an initialization of the base class object.
  }

  wMinusDirectJet::~wMinusDirectJet()
  {
    // does nothing.
  }



  negativelyChargedHiggsMuon::negativelyChargedHiggsMuon(
                                      readierForNewPoint* const readierPointer,
                                   CppSLHA::CppSLHA0 const* const spectrumData,
                     CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                     CppSLHA::particle_property_set const* const thirdParticle,
                 CppSLHA::particle_property_set const* const fourthParticle ) :
    HiggsMuonPlusAntimuon( readierPointer,
                           spectrumData,
                           firstParticle,
                           effectiveSquarkMass,
                           secondParticle,
                           thirdParticle,
                           fourthParticle )
  {
    // just an initialization of the base class object.
  }

  negativelyChargedHiggsMuon::~negativelyChargedHiggsMuon()
  {
    // does nothing.
  }



  negativelyChargedHiggsJet::negativelyChargedHiggsJet(
                                      readierForNewPoint* const readierPointer,
                                   CppSLHA::CppSLHA0 const* const spectrumData,
                     CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                     CppSLHA::particle_property_set const* const thirdParticle,
                 CppSLHA::particle_property_set const* const fourthParticle ) :
    HiggsMuonPlusAntimuon( readierPointer,
                           spectrumData,
                           firstParticle,
                           effectiveSquarkMass,
                           secondParticle,
                           thirdParticle,
                           fourthParticle )
  {
    // just an initialization of the base class object.
  }

  negativelyChargedHiggsJet::~negativelyChargedHiggsJet()
  {
    // does nothing.
  }



  neutralinoThreeBodyDecay::neutralinoThreeBodyDecay(
                                      readierForNewPoint* const readierPointer,
                                   CppSLHA::CppSLHA0 const* const spectrumData,
                     CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                     CppSLHA::particle_property_set const* const thirdParticle,
                      CppSLHA::particle_property_set const* const leftSfermion,
                  CppSLHA::particle_property_set const* const rightSfermion ) :
    leptonEnergyDistribution( readierPointer,
                                spectrumData,
                                firstParticle,
                                effectiveSquarkMass,
                                secondParticle,
                                thirdParticle,
                                NULL ),
    leftSfermion( leftSfermion ),
    rightSfermion( rightSfermion ),
    minToMaxConst( minToMaxSegment.addTerm( 0,
                                            0,
                                            NULL,
                                            1.0 ) )
  {
    segments.push_back( &minToMaxSegment );
  }

  neutralinoThreeBodyDecay::~neutralinoThreeBodyDecay()
  {
    // does nothing.
  }


  void
  neutralinoThreeBodyDecay::calculateCoefficients()
  {
    minimumEnergy = 0.0;
    maximumEnergy = ( secondMass - thirdMass );

    // now we set up the coefficients of various terms:

    // MIN to MAX segment:
    minToMaxSegment.setSegmentRange( minimumEnergy,
                                     maximumEnergy );
  }



  charginoThreeBodyDecay::charginoThreeBodyDecay(
                                      readierForNewPoint* const readierPointer,
                                   CppSLHA::CppSLHA0 const* const spectrumData,
                     CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                     CppSLHA::particle_property_set const* const thirdParticle,
             CppSLHA::particle_property_set const* const leftUpIsospinSfermion,
            CppSLHA::particle_property_set const* const rightUpIsospinSfermion,
           CppSLHA::particle_property_set const* const leftDownIsospinSfermion,
       CppSLHA::particle_property_set const* const rightDownIsospinSfermion ) :
    leptonEnergyDistribution( readierPointer,
                              spectrumData,
                              firstParticle,
                              effectiveSquarkMass,
                              secondParticle,
                              thirdParticle,
                              NULL ),
    leftUpIsospinSfermion( leftUpIsospinSfermion ),
    rightUpIsospinSfermion( rightUpIsospinSfermion ),
    leftDownIsospinSfermion( leftDownIsospinSfermion ),
    rightDownIsospinSfermion( rightDownIsospinSfermion ),
    minToMaxConst( minToMaxSegment.addTerm( 0,
                                            0,
                                            NULL,
                                            1.0 ) )
  {
    segments.push_back( &minToMaxSegment );
  }

  charginoThreeBodyDecay::~charginoThreeBodyDecay()
  {
    // does nothing.
  }


  void
  charginoThreeBodyDecay::calculateCoefficients()
  {
    minimumEnergy = 0.0;
    maximumEnergy = ( secondMass - thirdMass );

    // now we set up the coefficients of various terms:

    // MIN to MAX segment:
    minToMaxSegment.setSegmentRange( minimumEnergy,
                                     maximumEnergy );
  }

}  // end of LHC_FASER namespace.

