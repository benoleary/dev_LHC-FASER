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
 *      LHC-FASER also requires grids of lookup acceptanceValues. These should also be
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
    double mQSq( mq * mq );
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
                                / ( ( mlSq - mNSq ) * ( mNSq + mQSq )
                                    * ( mlSq - mXSq ) ) );

    minToLkLog->setCoefficient( ( 2.0 * mNSq * ( mQSq - mlSq ) )
                                / ( ( mNSq - mlSq ) * ( mNSq + mQSq ) ) );

    minToLkConst->setCoefficient(
          ( 2.0 * mNSq
            * ( 2.0 * mNSq - mQSq + mlSq
              + ( mQSq - mlSq ) * log( ( ( 2.0 * mq ) / ( mlSq - mXSq ) ) ) ) )
        / ( ( mNSq - mlSq ) * ( mNSq + mQSq ) ) );

    minToLkInv->setCoefficient( ( mNSq * ( 2.0 * mNSq - mQSq )
                                  * ( mlSq - mXSq ) )
                                / ( mq * ( mlSq - mNSq ) * ( mNSq + mQSq ) ) );

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
              ( 2.0 * ( ( mlSq + mNSq ) * ( mNSq - mQSq )
                      + mNSq * ( mlSq - mQSq ) * log( ( mNSq / mQSq ) ) ) )
              / ( ( mNSq - mlSq ) * ( mNSq + mQSq ) ) );

      lkToHkInv->setCoefficient(
                            ( 2.0 * mNSq * ( mNSq - mQSq ) * ( mlSq - mXSq ) )
                            / ( mq * ( mlSq - mNSq ) * ( mNSq + mQSq ) ) );
    }
    else
    {

      lkToHkConst->setCoefficient(
                ( 2.0 * ( ( mQSq - mNSq ) * ( mlSq - mNSq )
                          + mNSq * ( mlSq - mQSq ) * log( ( mlSq / mNSq ) ) ) )
              / ( ( mNSq - mlSq ) * ( mNSq + mQSq ) ) );

      lkToHkInv->setCoefficient( 0.0 );
    }

    // HK to MAX segment:
    hkToMaxSegment.setSegmentRange( Ehk,
                                    maximumEnergy );

    hkToMaxLin->setCoefficient( -minToLkLin->getCoefficient() );

    hkToMaxLog->setCoefficient( -minToLkLog->getCoefficient() );

    hkToMaxConst->setCoefficient(
     ( 2.0 * mNSq * ( mlSq + mQSq
                      + ( mQSq - mlSq ) * log( ( 2.0 * mlSq )
                                               / ( ( mlSq - mXSq ) * mq ) ) ) )
          / ( ( mlSq - mNSq ) * ( mNSq + mQSq ) ) );

    hkToMaxInv->setCoefficient( ( mNSq * mq * ( mlSq - mXSq ) )
                                / ( ( mNSq - mlSq ) * ( mNSq + mQSq ) ) );

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
    double mQSq( mq * mq );
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
          / ( ( mNSq - mlSq ) * ( mNSq + mQSq ) * ( mlSq - mXSq ) ) ) );

    minToLkLog->setCoefficient(
        ( ( 2.0 * ( mNSq * mNSq - mlSq * mQSq ) )
          / ( ( mNSq - mlSq ) * ( mNSq + mQSq ) ) ) );

    minToLkConst->setCoefficient(
        ( ( 2.0 * ( mNSq * ( 2.0 * mNSq - mQSq + mlSq )
                  + ( mlSq * mQSq - mNSq * mNSq )
                    * log( ( ( 2.0 * mq ) / ( mlSq - mXSq ) ) ) ) )
          / ( ( mlSq - mNSq ) * ( mNSq + mQSq ) ) ) );

    minToLkInv->setCoefficient(
        ( ( mNSq * ( ( 2.0 * mNSq - mQSq ) * ( mlSq - mXSq ) ) )
          / ( mq * ( mNSq - mlSq ) * ( mNSq + mQSq ) ) ) );


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
            ( ( 2.0 * ( ( mlSq + mNSq ) * ( mNSq - mQSq )
                        + ( mNSq * mNSq - mlSq * mQSq )
                          * log( ( mNSq / mQSq ) ) ) )
              / ( ( mlSq - mNSq ) * ( mNSq + mQSq ) ) ) );

      lkToHkInv->setCoefficient(
            ( ( 2.0 * mNSq * ( mNSq - mQSq ) * ( mlSq - mXSq ) )
              / ( mq * ( mNSq - mlSq ) * ( mNSq + mQSq ) ) ) );
    }
    else
    {
      lkToHkConst->setCoefficient(
            ( ( 2.0 * ( ( mQSq - mNSq ) * ( mlSq - mNSq )
                        + ( mNSq * mNSq - mlSq * mQSq )
                          * log( ( mlSq / mNSq ) ) ) )
              / ( ( mlSq - mNSq ) * ( mNSq + mQSq ) ) ) );

      lkToHkInv->setCoefficient( 0.0 );
    }


    // HK to MAX segment:
    hkToMaxSegment.setSegmentRange( Ehk,
                                    maximumEnergy );

    hkToMaxLin->setCoefficient( -minToLkLin->getCoefficient() );

    hkToMaxLog->setCoefficient( -minToLkLog->getCoefficient() );

    hkToMaxConst->setCoefficient(
        ( ( 2.0 * ( ( mNSq * mNSq - mlSq * mQSq )
                    * log( ( ( 2.0 * mlSq ) / ( ( mlSq - mXSq ) * mq ) ) )
                    - mNSq * ( mQSq + mlSq ) ) )
          / ( ( mlSq - mNSq ) * ( mNSq + mQSq ) ) ) );

    hkToMaxInv->setCoefficient(
        ( ( mNSq * mq * ( mlSq - mXSq ) )
          / ( ( mlSq - mNSq ) * ( mNSq + mQSq ) ) ) );

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
    << "( 2.0 / ( ( mlSq - mNSq ) * ( mNSq + mQSq ) ) ) = "
    << ( 2.0 / ( ( mlSq - mNSq ) * ( mNSq + mQSq ) ) );
    std::cout << std::endl
    << "( ( mQSq - mNSq ) * ( mlSq - mNSq ) "
    << "+ ( mNSq * mNSq - mlSq * mQSq ) * log( ( mlSq / mNSq ) ) ) = "
    << ( ( mQSq - mNSq ) * ( mlSq - mNSq )
         + ( mNSq * mNSq - mlSq * mQSq ) * log( ( mlSq / mNSq ) ) );
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


  weakVectorBosonHandedMuon::weakVectorBosonHandedMuon(
                                             readierForNewPoint* const readier,
                                   CppSLHA::CppSLHA0 const* const spectrumData,
                     CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                     CppSLHA::particle_property_set const* const thirdParticle,
                    CppSLHA::particle_property_set const* const fourthParticle,
                               bool const negativeMuonIsSameHandednessAsJet ) :
    leptonEnergyDistribution( readier,
                              spectrumData,
                              firstParticle,
                              effectiveSquarkMass,
                              secondParticle,
                              thirdParticle,
                              fourthParticle ),
    negativeMuonIsSameHandednessAsJet( negativeMuonIsSameHandednessAsJet ),
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

  weakVectorBosonHandedMuon::~weakVectorBosonHandedMuon()
  {
    // does nothing.
  }


  void
  weakVectorBosonHandedMuon::calculateEnergiesAndFactors()
  {
    mQ = firstMass;  // this is just for compact equations.
    //mQSq = ( mQ * mQ );
    //mQCu = ( mQ * mQSq );
    mC = secondMass;  // this is just for compact equations.
    mCSq = ( mC * mC );
    mCCu = ( mCSq * mC );
    mVBSq = ( thirdMass * thirdMass );
    mVBQu = ( mVBSq * mVBSq );
    mCX = ( mC * fourthMass );
    mXSq = ( fourthMass * fourthMass );

    gammaCQ = ( ( mQ * mQ + mCSq ) / ( 2.0 * mC * mQ ) );
    gammaCQSq = ( gammaCQ * gammaCQ );
    gammaCQCu = ( gammaCQSq * gammaCQ );
    gammaCQQu = ( gammaCQCu * gammaCQ );
    betaCQSq = ( 1.0 - ( 1.0 / gammaCQSq ) );
    betaCQ = sqrt( betaCQSq );
    oneMinusBetaCQ = ( 1.0 - betaCQ );
    oneMinusBetaCQSq = ( oneMinusBetaCQ * oneMinusBetaCQ );
    oneMinusBetaCQCu = ( oneMinusBetaCQSq * oneMinusBetaCQ );
    oneMinusBetaCQQu = ( oneMinusBetaCQCu * oneMinusBetaCQ );
    onePlusBetaCQ = ( 1.0 + betaCQ );
    lnOnePlusBetaCQOverOneMinusBetaCQ = log( onePlusBetaCQ / oneMinusBetaCQ );

    EVB = ( ( mCSq + mVBSq - mXSq ) / ( 2.0 * mC ) );
    ElMin = ( 0.5 * EVB * ( 1.0 - sqrt( 1.0 - ( mVBSq / ( EVB * EVB ) ) ) ) );
    lnElMin = log( ElMin );
    ElMinSq = ( ElMin * ElMin );
    ElMinCu = ( ElMinSq * ElMin );
    ElMinQu = ( ElMinCu * ElMin );
    minimumEnergy = ( ElMin * ( mC / mQ ) );
    Elk = ( ElMin * ( mQ / mC ) );
    Ehk = ( mVBSq / ( 4.0 * Elk ) );
    maximumEnergy = ( mVBSq / ( 4.0 * minimumEnergy ) );
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
    //EightmCEVElMin = ( 8.0 * mC * EVB * ElMin );
    //lnmQC = log( ( mQ / mC ) );
    lnmVBSqOverFourElMinSq = log( ( mVBSq / ( 4.0 * ElMinSq ) ) );
    //EightmCXElMinmQSq = ( 8.0 * mC * fourthMass * ElMin * mQSq );
    //lnEmin = log( minimumEnergy );
    //lnElk = log( Elk );
    //lnEhk = log( Ehk );
    //lnEmax = log( maximumEnergy );

    minToLkSegment.setSegmentRange( minimumEnergy,
                                    Elk );
    lkToHkSegment.setSegmentRange( Elk,
                                   Ehk );
    hkToMaxSegment.setSegmentRange( Ehk,
                                    maximumEnergy );
  }

  void
  weakVectorBosonHandedMuon::calculateVvPlusJjAaAllSqSymCoefficients()
  {
    /* ends up in the order:
     * low (MINtoLK),
     * over (LKtoHK, true == cosinesLimitedByEnergy),
     * under (LKtoHK, false == cosinesLimitedByEnergy),
     * high (HKtoMAX)
     */

    // MIN to LK:
    currentMinToLkInvCoefficient = 0.0;
    // right!
    currentMinToLkConstCoefficient
    = ( 4.0 * betaCQ * ElMinQu * gammaCQSq * mCSq * oneMinusBetaCQSq
        * ( 4.0 * ElMin * ( ElMin - 2.0 * EVB ) * mC
            + EVB * mVBSq
              * ( 2.0 * lnElMin - lnOnePlusBetaCQOverOneMinusBetaCQ ) ) );
    // right!
    currentMinToLkLogCoefficient
    = ( -8.0 * betaCQ * ElMinQu * EVB * gammaCQSq * mCSq * mVBSq
        * oneMinusBetaCQSq );
    // right!
    currentMinToLkLinCoefficient
    = ( 32.0 * betaCQ * ElMinQu * EVB * gammaCQ * mCCu * oneMinusBetaCQ );
    // right!
    currentMinToLkLinLogCoefficient = 0.0;
    // right!
    currentMinToLkSqCoefficient = ( -16.0 * betaCQ * ElMinQu * mCCu );
    // right!
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "betaCQ = " << betaCQ
    << std::endl
    << "ElMinQu = " << ElMinQu
    << std::endl
    << "mCCu = " << mCCu;
    std::cout << std::endl;**/

    // LK to HK:
    if( cosinesLimitedByEnergy )
    {
      currentLkToHkConstCoefficient
      = ( betaCQ * ElMinSq * gammaCQSq * mCSq * oneMinusBetaCQSq
          * ( 16.0 * ElMinCu * ( ElMin - 2.0 * EVB ) * mC
              - mC * mVBQu
              + 8.0 * ElMin * EVB * mVBSq
                * ( mC - ElMin * lnmVBSqOverFourElMinSq ) ) );
      // right!
      currentLkToHkLinCoefficient = 0.0;
      // right!
      currentLkToHkSqCoefficient = 0.0;
      // right!
    }
    else
    {
      currentLkToHkConstCoefficient
      = ( -8.0 * betaCQ * ElMinQu * EVB * gammaCQSq * mCSq * mVBSq
          * lnOnePlusBetaCQOverOneMinusBetaCQ * oneMinusBetaCQSq );
      // right!
      currentLkToHkLinCoefficient
      = ( 64.0 * betaCQSq * ElMinQu * EVB * gammaCQCu * mCCu
          * oneMinusBetaCQSq );
      // right!
      currentLkToHkSqCoefficient
      = ( -64.0 * betaCQSq * ElMinQu * gammaCQQu * mCCu * oneMinusBetaCQSq );
      // right!
    }

    // HK to MAX:
    currentHkToMaxInvCoefficient = 0.0;
    // right!
    currentHkToMaxConstCoefficient
    = ( betaCQ * ElMinSq * gammaCQSq * mCSq * mVBSq
        * ( 4.0 * ElMin * EVB
            * ( 2.0 * ( mC - ElMin * ( lnElMin + lnmVBSqOverFourElMinSq ) )
                - ElMin * lnOnePlusBetaCQOverOneMinusBetaCQ )
            - mC * mVBSq ) * oneMinusBetaCQSq );
    // right!
    currentHkToMaxLogCoefficient = -currentMinToLkLogCoefficient;
    // right!
    currentHkToMaxLinCoefficient
    = ( -currentMinToLkLinCoefficient * gammaCQSq * oneMinusBetaCQSq );
    //= ( -32.0 * betaCQ * ElMinQu * EVB * gammaCQCu * mCCu
    //* oneMinusBetaCQCu );
    // right!
    currentHkToMaxLinLogCoefficient = 0.0;
    // right!
    currentHkToMaxSqCoefficient
    = ( -currentMinToLkSqCoefficient * gammaCQQu * oneMinusBetaCQQu );
    //= ( 16.0 * betaCQ * ElMinQu * gammaCQQu * mCCu * oneMinusBetaCQQu );
    // right!
  }

  void
  weakVectorBosonHandedMuon::calculateTwiceVvAaSymCoefficients()
  {
    /* ends up in the order:
     * low (MINtoLK),
     * over (LKtoHK, true == cosinesLimitedByEnergy),
     * under (LKtoHK, false == cosinesLimitedByEnergy),
     * high (HKtoMAX)
     */

    // MIN to LK:
    currentMinToLkInvCoefficient
    = ( 2.0 * ElMinQu * gammaCQ * mCSq * mVBQu * oneMinusBetaCQSq );
    // right!
    currentMinToLkConstCoefficient
    = ( 4.0 * ElMinCu * gammaCQSq * mCSq * oneMinusBetaCQSq
        * ( ElMin * mVBSq
            * ( EVB * oneMinusBetaCQ
                * ( 2.0 + lnOnePlusBetaCQOverOneMinusBetaCQ - 2.0 * lnElMin )
                + mC * ( lnOnePlusBetaCQOverOneMinusBetaCQ
                         - 2.0 * lnElMin
                         + 2.0 * oneMinusBetaCQ ) )
            - mVBQu
            - 4.0 * ElMinSq * ( ElMin - 2.0 * EVB ) * mC * oneMinusBetaCQ ) );
    // right!
    currentMinToLkLogCoefficient
    = ( 8.0 * ElMinQu * gammaCQSq * mCSq * mVBSq
        * ( mC + EVB * oneMinusBetaCQ ) * oneMinusBetaCQSq );
    // right!
    currentMinToLkLinCoefficient
    = ( 2.0 * ElMinSq * gammaCQ * mCSq
        * ( 8.0 * ElMinSq * mC
            * ( 2.0 * ElMin
                + EVB * ( lnOnePlusBetaCQOverOneMinusBetaCQ
                          - 2.0 * ( 1.0 + lnElMin ) ) )
            + mVBQu
            - 4.0 * ElMin * ( EVB + mC ) * mVBSq ) * oneMinusBetaCQSq );
    // right!
    currentMinToLkLinLogCoefficient
    = ( 32.0 * ElMinQu * EVB * gammaCQ * mCCu * oneMinusBetaCQSq );
    // right!
    currentMinToLkSqCoefficient = ( -16.0 * ElMinQu * mCCu * oneMinusBetaCQ );
    // right!

    // LK to HK:
    if( cosinesLimitedByEnergy )
    {
      currentLkToHkConstCoefficient
      = ( ElMinSq * gammaCQSq * mCSq
          * ( ( mC - 4.0 * ElMin - betaCQ * mC ) * mVBQu
              - 16.0 * ElMinCu * ( ElMin - 2.0 * EVB ) * mC * oneMinusBetaCQ
              + 8.0 * ElMin * mVBSq
                * ( 2.0 * ElMinSq
                    + ElMin * lnmVBSqOverFourElMinSq
                      * ( oneMinusBetaCQ * EVB + mC )
                    - EVB * mC * oneMinusBetaCQ ) ) * oneMinusBetaCQSq );
      // right!
      currentLkToHkLinCoefficient
      = ( 2.0 * ElMinSq * gammaCQ * mCSq * oneMinusBetaCQSq
          * ( mVBQu
              - 16.0 * ElMinSq
                * ( ElMinSq - EVB * lnmVBSqOverFourElMinSq * mC
                    - ElMin * ( EVB + 2.0 * mC ) )
              - 4.0 * ElMin * ( EVB + 2.0 * mC ) * mVBSq ) );
      // right!
      currentLkToHkSqCoefficient = 0.0;
      // right!
    }
    else
    {
      currentLkToHkConstCoefficient
      = ( 8.0 * ElMinQu * gammaCQSq * mCSq * mVBSq * oneMinusBetaCQSq
          * ( lnOnePlusBetaCQOverOneMinusBetaCQ * ( EVB + mC )
              - betaCQ * ( EVB * ( 2.0 + lnOnePlusBetaCQOverOneMinusBetaCQ )
                           + 2.0 * mC ) ) );
      // dunno...
      currentLkToHkLinCoefficient
      = ( 32.0 * ElMinQu * EVB * gammaCQCu * mCCu * oneMinusBetaCQCu
          * ( onePlusBetaCQ * lnOnePlusBetaCQOverOneMinusBetaCQ
              - 2.0 * betaCQ ) );
      // dunno...
      currentLkToHkSqCoefficient
      = ( -64.0 * betaCQSq * ElMinQu * gammaCQQu * mCCu * oneMinusBetaCQCu );
      // dunno...
    }

    // HK to MAX:
    currentHkToMaxInvCoefficient = -currentMinToLkInvCoefficient;
    // right!
    currentHkToMaxConstCoefficient
    = ( ElMinSq * gammaCQSq * mCSq * mVBSq * oneMinusBetaCQSq
        * ( 16.0 * ElMinCu
            - 4.0 * ElMinSq
              * ( 2.0 * onePlusBetaCQ * ( EVB + mC )
                  - ( EVB * oneMinusBetaCQ + mC )
                    * ( lnOnePlusBetaCQOverOneMinusBetaCQ
                        + 2.0 * ( lnElMin + lnmVBSqOverFourElMinSq ) ) )
            + oneMinusBetaCQ * mC * ( mVBSq - 8.0 * ElMin * EVB ) ) );
    // right!
    currentHkToMaxLogCoefficient = -currentMinToLkLogCoefficient;
    // right!
    currentHkToMaxLinCoefficient
    = ( 8.0 * ElMinCu * gammaCQCu * mCSq * oneMinusBetaCQCu
        * ( 2.0 * ElMin * EVB
            * ( lnOnePlusBetaCQOverOneMinusBetaCQ
                + 2.0 * ( 1.0 + lnElMin + lnmVBSqOverFourElMinSq )
                + betaCQ
                  * ( lnOnePlusBetaCQOverOneMinusBetaCQ
                      + 2.0
                      * (-1.0 + lnElMin + lnmVBSqOverFourElMinSq ) ) ) * mC
                      - 4.0 * ElMinCu * onePlusBetaCQ
                      + 4.0 * ElMinSq * ( EVB + mC ) * onePlusBetaCQ
                      - mC * mVBSq * onePlusBetaCQ ) );
    // right!
    currentHkToMaxLinLogCoefficient = -currentMinToLkLinLogCoefficient;
    // right!
    currentHkToMaxSqCoefficient
    = ( 16.0 * ( 1.0 + 3.0 * betaCQ ) * ElMinQu * gammaCQQu * mCCu
        * oneMinusBetaCQQu );
    // right!
  }

  void
  weakVectorBosonHandedMuon::calculateVvSqMinusAaSqSymCoefficients()
  {
    /* ends up in the order:
     * low (MINtoLK),
     * over (LKtoHK, true == cosinesLimitedByEnergy),
     * under (LKtoHK, false == cosinesLimitedByEnergy),
     * high (HKtoMAX)
     */

    // MIN to LK:
    currentMinToLkInvCoefficient = 0.0;
    // right!
    currentMinToLkConstCoefficient
    = ( 4.0 * betaCQ * ElMinQu * gammaCQSq * mC * mCX * oneMinusBetaCQSq
         * mVBSq * ( lnOnePlusBetaCQOverOneMinusBetaCQ - 2.0 * lnElMin ) );
    // right!
    currentMinToLkLogCoefficient
    = ( 8.0 * betaCQ * ElMinQu * gammaCQSq * mC * mCX * mVBSq
        * oneMinusBetaCQSq );
    // right!
    currentMinToLkLinCoefficient = 0.0;
    // right!
    currentMinToLkLinLogCoefficient = 0.0;
    // right!
    currentMinToLkSqCoefficient = 0.0;
    // right!

    // LK to HK:
    if( cosinesLimitedByEnergy )
    {
      currentLkToHkConstCoefficient
      = ( 8.0 * betaCQ * ElMinQu * gammaCQSq * lnmVBSqOverFourElMinSq * mC
          * mCX * mVBSq * oneMinusBetaCQSq );
      // right!
    }
    else
    {
      currentLkToHkConstCoefficient
      = ( 8.0 * betaCQ * ElMinQu * gammaCQSq * mC * mCX * mVBSq
          * oneMinusBetaCQSq * lnOnePlusBetaCQOverOneMinusBetaCQ );
      // right!
    }
    currentLkToHkLinCoefficient = 0.0;
    // right!
    currentLkToHkSqCoefficient = 0.0;
    // right!

    // HK to MAX:
    currentHkToMaxInvCoefficient = 0.0;
    // right!
    currentHkToMaxConstCoefficient
    = ( 4.0 * betaCQ * ElMinQu * gammaCQSq * mC * mCX * mVBSq
        * ( lnOnePlusBetaCQOverOneMinusBetaCQ
            + 2.0 * ( lnElMin + lnmVBSqOverFourElMinSq ) )
        * oneMinusBetaCQSq );
    // right!
    currentHkToMaxLogCoefficient = -currentMinToLkLogCoefficient;
    // right!
    currentHkToMaxLinCoefficient = 0.0;
    // right!
    currentHkToMaxLinLogCoefficient = 0.0;
    // right!
    currentHkToMaxSqCoefficient = 0.0;
    // right!
  }

  void
  weakVectorBosonHandedMuon::calculateVvPlusJjAaAllSqAntiCoefficients()
  {
    /* ends up in the order:
     * low (MINtoLK),
     * over (LKtoHK, true == cosinesLimitedByEnergy),
     * under (LKtoHK, false == cosinesLimitedByEnergy),
     * high (HKtoMAX)
     */

    // MIN to LK:
    currentMinToLkInvCoefficient = 0.0;
    // right!
    currentMinToLkConstCoefficient
    = ( 4.0 * betaCQ * ElMinQu * gammaCQSq * mCSq * mVBSq * oneMinusBetaCQSq
        * ( EVB * ( 2.0 * lnElMin - lnOnePlusBetaCQOverOneMinusBetaCQ )
            - 4.0 * ElMin ) );
    // right!
    currentMinToLkLogCoefficient
    = ( -8.0 * betaCQ * ElMinQu * EVB * gammaCQSq * mCSq * mVBSq
        * oneMinusBetaCQSq );
    // right!
    currentMinToLkLinCoefficient
    = ( 16.0 * betaCQ * ElMinQu * gammaCQ * mCSq * mVBSq * oneMinusBetaCQ );
    // right!
    currentMinToLkLinLogCoefficient = 0.0;
    // right!
    currentMinToLkSqCoefficient = 0.0;
    // right!

    // LK to HK:
    if( cosinesLimitedByEnergy )
    {
      currentLkToHkConstCoefficient
      = ( 4.0 * betaCQ * ElMinCu * gammaCQSq * mCSq * mVBSq * oneMinusBetaCQSq
          * ( mVBSq
              - 4.0 * ElMinSq
              - 2.0 * ElMin * EVB * lnmVBSqOverFourElMinSq ) );
      // right!
      currentLkToHkLinCoefficient = 0.0;
      // right!
    }
    else
    {
      currentLkToHkConstCoefficient
      = ( -8.0 * betaCQ * ElMinQu * EVB * gammaCQSq * mVBSq * oneMinusBetaCQSq
          * lnOnePlusBetaCQOverOneMinusBetaCQ * mCSq );
      // right!
      currentLkToHkLinCoefficient
      = ( 32.0 * betaCQSq * ElMinQu * gammaCQCu * mCSq * mVBSq
          * oneMinusBetaCQSq );
      // right!
    }
    currentLkToHkSqCoefficient = 0.0;
    // right!

    // HK to MAX:
    currentHkToMaxInvCoefficient = 0.0;
    // right!
    currentHkToMaxConstCoefficient
    = ( 4.0 * betaCQ * ElMinCu * gammaCQSq * mCSq * mVBSq * oneMinusBetaCQSq
        * ( mVBSq
            - ElMin * EVB
            * ( lnOnePlusBetaCQOverOneMinusBetaCQ
                + 2.0 * ( lnElMin + lnmVBSqOverFourElMinSq ) ) ) );
    // right!
    currentHkToMaxLogCoefficient = -currentMinToLkLogCoefficient;
    // right!
    currentHkToMaxLinCoefficient
    = ( -currentMinToLkLinCoefficient * gammaCQSq * oneMinusBetaCQSq );
    // right!
    currentHkToMaxLinLogCoefficient = 0.0;
    // right!
    currentHkToMaxSqCoefficient = 0.0;
    // right!
  }

  void
  weakVectorBosonHandedMuon::calculateVvSqPlusAaSqAntiCoefficients()
  {
    /* ends up in the order:
     * low (MINtoLK),
     * over (LKtoHK, true == cosinesLimitedByEnergy),
     * under (LKtoHK, false == cosinesLimitedByEnergy),
     * high (HKtoMAX)
     */

    // MIN to LK:
    currentMinToLkInvCoefficient
    = ( -2.0 * ElMinQu * gammaCQ * mCSq * mVBQu * oneMinusBetaCQSq );
    // right!
    currentMinToLkConstCoefficient
    = ( 4.0 * ElMinCu * gammaCQSq * mCSq * mVBSq * oneMinusBetaCQSq
        * ( mVBSq
            - ElMin * oneMinusBetaCQ
              * ( 4.0 * ElMin
                  + EVB * ( 2.0
                            + lnOnePlusBetaCQOverOneMinusBetaCQ
                            - 2.0 * lnElMin ) )
            - ElMin * mC * ( lnOnePlusBetaCQOverOneMinusBetaCQ
                             + 2.0 * ( oneMinusBetaCQ - lnElMin ) ) ) );
    // right!
    currentMinToLkLogCoefficient
    = ( -8.0 * ElMinQu * gammaCQSq * mCSq * mVBSq
        * ( mC + EVB * oneMinusBetaCQ ) * oneMinusBetaCQSq );
    // right!
    currentMinToLkLinCoefficient
    = ( 2.0 * ElMinSq * gammaCQ * mCSq
        * ( 4.0 * ElMin
            * ( EVB
                + ElMin
                  * ( 2.0 + 2.0 * lnElMin - lnOnePlusBetaCQOverOneMinusBetaCQ )
                + mC )
            - mVBSq ) * mVBSq * oneMinusBetaCQSq );
    // right!
    currentMinToLkLinLogCoefficient
    = ( -16.0 * ElMinQu * gammaCQ * mCSq * mVBSq * oneMinusBetaCQSq );
    // right!
    currentMinToLkSqCoefficient = 0.0;
    // right!

    // LK to HK:
    if( cosinesLimitedByEnergy )
    {
      currentLkToHkConstCoefficient
      = ( 4.0 * ElMinCu * gammaCQSq * mCSq * mVBSq
          * ( ( 2.0 - betaCQ ) * ( mVBSq - 4.0 * ElMinSq )
              - 2.0 * ElMin * lnmVBSqOverFourElMinSq
              * ( mC + EVB * oneMinusBetaCQ ) ) * oneMinusBetaCQSq );
      // right!
      currentLkToHkLinCoefficient
      = ( 2.0 * ElMinSq * gammaCQ * mCSq
          * ( 16.0 * ElMinCu * ( ElMin - EVB - mC )
              - mVBQu
              + 4.0 * ElMin
                * ( EVB - 2.0 * ElMin * lnmVBSqOverFourElMinSq + mC ) * mVBSq )
                * oneMinusBetaCQSq );
      // right!
    }
    else
    {
      currentLkToHkConstCoefficient
      = ( 8.0 * ElMinQu * gammaCQSq * mCSq * mVBSq * oneMinusBetaCQSq
          * ( betaCQ * EVB * ( 2.0 + lnOnePlusBetaCQOverOneMinusBetaCQ )
              + 2.0 * betaCQ * mC
              - lnOnePlusBetaCQOverOneMinusBetaCQ * ( EVB + mC ) ) );
      // right!
      currentLkToHkLinCoefficient
      = ( 16 * ElMinQu * gammaCQCu * mCSq * mVBSq * oneMinusBetaCQCu
          * ( 2.0 * betaCQ
              - onePlusBetaCQ * lnOnePlusBetaCQOverOneMinusBetaCQ ) );
      // right!
    }
    currentLkToHkSqCoefficient = 0.0;
    // right!

    // HK to MAX:
    currentHkToMaxInvCoefficient = -currentMinToLkInvCoefficient;
    // right!
    currentHkToMaxConstCoefficient
    = ( 4.0 * ElMinCu * gammaCQSq * mCSq * mVBSq * oneMinusBetaCQSq
        * ( mVBSq
            - 4.0 * ElMinSq
            - ElMin * EVB
              * ( lnOnePlusBetaCQOverOneMinusBetaCQ
                  + 2.0 * ( lnElMin + lnmVBSqOverFourElMinSq - 1.0 )
                  - betaCQ
                    * ( lnOnePlusBetaCQOverOneMinusBetaCQ
                        + 2.0 * ( 1.0 + lnElMin + lnmVBSqOverFourElMinSq ) ) )
            + ElMin * mC
              * ( 2.0 * ( onePlusBetaCQ - lnElMin - lnmVBSqOverFourElMinSq )
                  - lnOnePlusBetaCQOverOneMinusBetaCQ )
            - betaCQ * mVBSq ) );
    // right!
    currentHkToMaxLogCoefficient = -currentMinToLkLogCoefficient;
    // right!
    currentHkToMaxLinCoefficient
    = ( 8.0 * ElMinQu * gammaCQCu * mCSq * oneMinusBetaCQCu
        * ( 4.0 * ElMin * onePlusBetaCQ * ( ElMin - EVB - mC )
            - mVBSq * ( onePlusBetaCQ
                        * ( lnOnePlusBetaCQOverOneMinusBetaCQ
                            + 2.0 * ( lnElMin + lnmVBSqOverFourElMinSq ) )
                        + 2.0 * oneMinusBetaCQ ) ) );
    // right!
    currentHkToMaxLinLogCoefficient = -currentMinToLkLinLogCoefficient;
    // right!
    currentHkToMaxSqCoefficient = 0.0;
    // right!
  }

  void
  weakVectorBosonHandedMuon::calculateVvSqMinusAaSqAntiCoefficients()
  {
    /* ends up in the order:
     * low (MINtoLK),
     * over (LKtoHK, true == cosinesLimitedByEnergy),
     * under (LKtoHK, false == cosinesLimitedByEnergy),
     * high (HKtoMAX)
     */

    // MIN to LK:
    currentMinToLkInvCoefficient = 0.0;
    // right!
    currentMinToLkConstCoefficient
    = ( 4.0 * ElMinQu * gammaCQSq * mC * mCX * mVBSq * oneMinusBetaCQSq
        * ( 2.0 * ( lnElMin - oneMinusBetaCQ )
            - lnOnePlusBetaCQOverOneMinusBetaCQ ) );
    // right!
    currentMinToLkLogCoefficient
    = ( -8.0 * ElMinQu * gammaCQSq * mC * mCX * mVBSq * oneMinusBetaCQSq );
    // right!
    currentMinToLkLinCoefficient
    = ( 8.0 * ElMinCu * gammaCQ * mC * mCX * mVBSq * oneMinusBetaCQSq );
    // right!
    currentMinToLkLinLogCoefficient = 0.0;
    // right!
    currentMinToLkSqCoefficient = 0.0;
    // right!

    // LK to HK:
    if( cosinesLimitedByEnergy )
    {
      currentLkToHkConstCoefficient
      = ( -8.0 * ElMinQu * gammaCQSq * lnmVBSqOverFourElMinSq * mC * mCX
          * mVBSq * oneMinusBetaCQSq );
      // right!
      currentLkToHkLinCoefficient
      = ( 8.0 * ElMinCu * gammaCQ * mC * mCX * ( mVBSq - 4.0 * ElMinSq )
          * oneMinusBetaCQSq );
      // right!
    }
    else
    {
      currentLkToHkConstCoefficient
      = ( 8.0 * ElMinQu * gammaCQSq * mC * mCX * mVBSq * oneMinusBetaCQSq
          * ( 2.0 * betaCQ - lnOnePlusBetaCQOverOneMinusBetaCQ ) );
      // right!
      currentLkToHkLinCoefficient = 0.0;
      // right!
    }
    currentLkToHkSqCoefficient = 0.0;
    // right!

    // HK to MAX:
    currentHkToMaxInvCoefficient = 0.0;
    // right!
    currentHkToMaxConstCoefficient
    = ( 4.0 * ElMinQu * gammaCQSq * mC * mCX * mVBSq * oneMinusBetaCQSq
        * ( 2.0 * ( onePlusBetaCQ - lnElMin - lnmVBSqOverFourElMinSq )
            - lnOnePlusBetaCQOverOneMinusBetaCQ ) );
    // right!
    currentHkToMaxLogCoefficient = -currentMinToLkLogCoefficient;
    // right!
    currentHkToMaxLinCoefficient
    = ( -32.0 * ElMinQu * ElMin * gammaCQ * mC * mCX * oneMinusBetaCQSq );
    // right!
    currentHkToMaxLinLogCoefficient = 0.0;
    // right!
    currentHkToMaxSqCoefficient = 0.0;
    // right!
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
    weakVectorBosonHandedMuon( readier,
                               spectrumData,
                               firstParticle,
                               effectiveSquarkMass,
                               secondParticle,
                               thirdParticle,
                               fourthParticle,
                               negativeMuonIsSameHandednessAsJet ),
    shouldNotSumOverHandedness( !shouldSumOverHandedness )
  {
    // just an initialization list.
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
      couplesAsAxialNotVector = true;
    }
    else
    {
      couplesAsAxialNotVector = false;
    }
    calculateEnergiesAndFactors();

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

    calculateVvPlusJjAaAllSqSymCoefficients();
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "VvPlusJjAaAllSqSymDistribution:"
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
    << std::endl << currentHkToMaxSqCoefficient << "E^2";
    std::cout << std::endl;**/
    setCurrentCoefficientsAsTotals();

    calculateVvSqMinusAaSqSymCoefficients();
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
    << std::endl << currentHkToMaxSqCoefficient << "E^2";
    std::cout << std::endl;**/
    if( couplesAsAxialNotVector )
    {
      flipSignsOfCurrentCoefficients();
    }
    addCurrentCoefficientsToTotals();

    if( shouldNotSumOverHandedness )
    {
      calculateVvPlusJjAaAllSqAntiCoefficients();
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
      << std::endl << currentHkToMaxSqCoefficient << "E^2";
      std::cout << std::endl;**/
      addCurrentCoefficientsToTotals();

      calculateVvSqPlusAaSqAntiCoefficients();
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
      << std::endl << currentHkToMaxSqCoefficient << "E^2";
      std::cout << std::endl;**/
      addCurrentCoefficientsToTotals();

      calculateVvSqMinusAaSqAntiCoefficients();
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
      << std::endl << currentHkToMaxSqCoefficient << "E^2";
      std::cout << std::endl;**/
      if( couplesAsAxialNotVector )
      {
        flipSignsOfCurrentCoefficients();
      }
      addCurrentCoefficientsToTotals();

    }
    // debugging:
    /**std::cout
    << std::endl
    << "debugging: zHandedMuon:calculateCoefficients() produced "
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
    weakVectorBosonHandedMuon( readierPointer,
                               spectrumData,
                               firstParticle,
                               effectiveSquarkMass,
                               secondParticle,
                               thirdParticle,
                               fourthParticle,
                               jetIsLeftHanded
                               /* W^- decays to left-handed mu^- */ )
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
      << "LHC-FASER::error! wMinusHandedMuon::wMinusHandedMuon expected to be"
      << " given a chargino as its 2nd particle, but did not receive a"
      << " particle with code +-" << CppSLHA::PDG_code::chargino_one
      << " or +-" << CppSLHA::PDG_code::chargino_two
      << ", rather, it received " << secondParticle->get_name() << " (code: "
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
                            2 ) * ( UMIX->get_entry( whichChargino,
                                                     1 )
                                    - VMIX->get_entry( whichChargino,
                                                       1 ) )
           + M_SQRT1_2 * ( NMIX->get_entry( whichNeutralino,
                                            3 )
                           * UMIX->get_entry( whichChargino,
                                              2 )
                           + NMIX->get_entry( whichNeutralino,
                                              4 )
                             * VMIX->get_entry( whichChargino,
                                                2 ) ) );
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "Vv = " << Vv << ", Aa = " << Aa;
    std::cout << std::endl;**/
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

    if( negativeMuonIsSameHandednessAsJet )
    {
      MinusTwiceVvJjAa = ( -2.0 * Vv * Aa );
      JjVvPlusJjAaAllSq = ( Vvsq - MinusTwiceVvJjAa + Aasq );
      JjVvSqMinusAaSq = VvSqMinusAaSq;
      JjVvSqPlusAaSq = ( Vvsq + Aasq );
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "jetIsLeftHanded == true";
      std::cout << std::endl;**/

    }
    else
    {
      MinusTwiceVvJjAa = ( 2.0 * Vv * Aa );
      JjVvPlusJjAaAllSq = ( MinusTwiceVvJjAa - Vvsq - Aasq );
      JjVvSqMinusAaSq = -VvSqMinusAaSq;
      JjVvSqPlusAaSq = ( -Vvsq - Aasq );
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "jetIsLeftHanded != true";
      std::cout << std::endl;**/
    }
    VvPlusJjAaAllSq = ( Vvsq - MinusTwiceVvJjAa + Aasq );

    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "VvPlusJjAaAllSq = " << VvPlusJjAaAllSq
    << std::endl
    << "MinusTwiceVvJjAa = " << MinusTwiceVvJjAa
    << std::endl
    << "VvSqMinusAaSq = " << VvSqMinusAaSq
    << std::endl
    << "JjVvPlusJjAaAllSq = " << JjVvPlusJjAaAllSq
    << std::endl
    << "JjVvSqPlusAaSq = " << JjVvSqPlusAaSq
    << std::endl
    << "JjVvSqMinusAaSq = " << JjVvSqMinusAaSq;
    std::cout << std::endl;**/

    calculateEnergiesAndFactors();

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

    calculateVvPlusJjAaAllSqSymCoefficients();
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "VvPlusJjAaAllSqSymDistribution:"
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
    << std::endl << currentHkToMaxSqCoefficient << "E^2";
    std::cout << std::endl;**/
    scaleCurrentCoefficients( VvPlusJjAaAllSq );
    setCurrentCoefficientsAsTotals();

    calculateTwiceVvAaSymCoefficients();
    // debugging:
    /**std::cout << std::endl << "debugging:"
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
    << std::endl << currentHkToMaxSqCoefficient << "E^2";
    std::cout << std::endl;**/
    scaleCurrentCoefficients( MinusTwiceVvJjAa );
    addCurrentCoefficientsToTotals();

    calculateVvSqMinusAaSqSymCoefficients();
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
    << std::endl << currentHkToMaxSqCoefficient << "E^2";
    std::cout << std::endl;**/
    scaleCurrentCoefficients( VvSqMinusAaSq );
    addCurrentCoefficientsToTotals();

    calculateVvPlusJjAaAllSqAntiCoefficients();
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
    << std::endl << currentHkToMaxSqCoefficient << "E^2";
    std::cout << std::endl;**/
    scaleCurrentCoefficients( JjVvPlusJjAaAllSq );
    addCurrentCoefficientsToTotals();

    calculateVvSqPlusAaSqAntiCoefficients();
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
    << std::endl << currentHkToMaxSqCoefficient << "E^2";
    std::cout << std::endl;**/
    scaleCurrentCoefficients( JjVvSqPlusAaSq );
    addCurrentCoefficientsToTotals();

    calculateVvSqMinusAaSqAntiCoefficients();
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
    << std::endl << currentHkToMaxSqCoefficient << "E^2";
    std::cout << std::endl;**/
    scaleCurrentCoefficients( JjVvSqMinusAaSq );
    addCurrentCoefficientsToTotals();

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

