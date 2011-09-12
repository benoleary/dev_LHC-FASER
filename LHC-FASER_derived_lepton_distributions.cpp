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

  flat_near_muon_plus_antimuon::flat_near_muon_plus_antimuon(
                                       readierForNewPoint* const given_readier,
                                 CppSLHA::CppSLHA0 const* const given_spectrum,
              CppSLHA::particle_property_set const* const given_first_particle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const given_second_particle,
           CppSLHA::particle_property_set const* const given_third_particle ) :
    leptonEnergyDistribution( given_readier,
                              given_spectrum,
                              given_first_particle,
                              effectiveSquarkMass,
                              given_second_particle,
                              given_third_particle,
                              NULL ),
    MINtoMAX_const( MINtoMAX_segment.addTerm( 0,
                                              0,
                                              NULL,
                                              1.0 ) )
  {
    segments.push_back( &MINtoMAX_segment );
  }

  flat_near_muon_plus_antimuon::~flat_near_muon_plus_antimuon()
  {
    // does nothing.
  }


  void
  flat_near_muon_plus_antimuon::calculateCoefficients()
  {
    productionFrameEnergy
    = ( ( ( secondMass * secondMass ) - ( thirdMass * thirdMass ) )
        / ( 2.0 * secondMass ) );
    minimumEnergy = ( ( secondMass / firstMass ) * productionFrameEnergy );
    maximumEnergy = ( ( firstMass / secondMass ) * productionFrameEnergy );


    // now we set up the coefficients of various terms.

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: flat_near_muon_plus_antimuon::calculateCoefficients has"
    << std::endl << "firstMass = " << firstMass
    << std::endl << "secondMass = " << secondMass
    << std::endl << "productionFrameEnergy = " << productionFrameEnergy
    << std::endl << "minimumEnergy = " << minimumEnergy
    << std::endl << "maximumEnergy = " << maximumEnergy;
    std::cout << std::endl;**/


    // finally we update the sets of terms, segment by segment:

    // MIN to MAX segment:
    MINtoMAX_segment.setSegmentRange( minimumEnergy,
                                      maximumEnergy );
  }



  same_chirality_near_muon::same_chirality_near_muon(
                                    readierForNewPoint* const given_readier,
                                 CppSLHA::CppSLHA0 const* const given_spectrum,
              CppSLHA::particle_property_set const* const given_first_particle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const given_second_particle,
           CppSLHA::particle_property_set const* const given_third_particle ) :
    leptonEnergyDistribution( given_readier,
                                given_spectrum,
                                given_first_particle,
                                effectiveSquarkMass,
                                given_second_particle,
                                given_third_particle,
                                NULL ),
    MINtoMAX_const( MINtoMAX_segment.addTerm( 0,
                                              0 ) ),
    MINtoMAX_lin( MINtoMAX_segment.addTerm( 1,
                                            0 ) )
  {
    segments.push_back( &MINtoMAX_segment );
  }

  same_chirality_near_muon::~same_chirality_near_muon()
  {

    // does nothing.

  }


  void
  same_chirality_near_muon::calculateCoefficients()
  {

    productionFrameEnergy
    = ( ( ( secondMass * secondMass ) - ( thirdMass * thirdMass ) )
        / ( 2.0 * secondMass ) );

    minimumEnergy
    = ( ( secondMass / firstMass ) * productionFrameEnergy );

    maximumEnergy
    = ( ( firstMass / secondMass ) * productionFrameEnergy );


    // now we set up the coefficients of various terms, & use them to determine
    // the normalization:

    // MIN to MAX segment:
    MINtoMAX_segment.setSegmentRange( minimumEnergy,
                                        maximumEnergy );

    MINtoMAX_const->setCoefficient( -minimumEnergy );

    MINtoMAX_lin->setCoefficient( 1.0 );
  }



  opposite_chirality_near_muon::opposite_chirality_near_muon(
                                    readierForNewPoint* const given_readier,
                                 CppSLHA::CppSLHA0 const* const given_spectrum,
              CppSLHA::particle_property_set const* const given_first_particle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const given_second_particle,
           CppSLHA::particle_property_set const* const given_third_particle ) :
    leptonEnergyDistribution( given_readier,
                                given_spectrum,
                                given_first_particle,
                                effectiveSquarkMass,
                                given_second_particle,
                                given_third_particle,
                                NULL ),
    MINtoMAX_const( MINtoMAX_segment.addTerm( 0,
                                              0 ) ),
    MINtoMAX_lin( MINtoMAX_segment.addTerm( 1,
                                            0 ) )
  {
    segments.push_back( &MINtoMAX_segment );
  }

  opposite_chirality_near_muon::~opposite_chirality_near_muon()
  {

    // does nothing.

  }


  void
  opposite_chirality_near_muon::calculateCoefficients()
  {

    productionFrameEnergy
    = ( ( ( secondMass * secondMass ) - ( thirdMass * thirdMass ) )
        / ( 2.0 * secondMass ) );

    minimumEnergy
    = ( ( secondMass / firstMass ) * productionFrameEnergy );

    maximumEnergy
    = ( ( firstMass / secondMass ) * productionFrameEnergy );


    // now we set up the coefficients of various terms, & use them to determine
    // the normalization:

    // MIN to MAX segment:
    MINtoMAX_segment.setSegmentRange( minimumEnergy,
                                        maximumEnergy );

    MINtoMAX_const->setCoefficient( maximumEnergy );

    MINtoMAX_lin->setCoefficient( -1.0 );
  }



  flat_far_muon_plus_antimuon::flat_far_muon_plus_antimuon(
                                    readierForNewPoint* const given_readier,
                                 CppSLHA::CppSLHA0 const* const given_spectrum,
              CppSLHA::particle_property_set const* const given_first_particle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const given_second_particle,
              CppSLHA::particle_property_set const* const given_third_particle,
          CppSLHA::particle_property_set const* const given_fourth_particle ) :
    leptonEnergyDistribution( given_readier,
                                given_spectrum,
                                given_first_particle,
                                effectiveSquarkMass,
                                given_second_particle,
                                given_third_particle,
                                given_fourth_particle ),
    MINtoLK_const( MINtoLK_segment.addTerm( 0,
                                            0 ) ),
    MINtoLK_log( MINtoLK_segment.addTerm( 0,
                                          1 ) ),
    LKtoHK_const( LKtoHK_segment.addTerm( 0,
                                          0 ) ),
    HKtoMAX_const( HKtoMAX_segment.addTerm( 0,
                                            0 ) ),
    HKtoMAX_log( HKtoMAX_segment.addTerm( 0,
                                          1 ) )
  {
    segments.push_back( &MINtoLK_segment );
    segments.push_back( &LKtoHK_segment );
    segments.push_back( &HKtoMAX_segment );
  }

  flat_far_muon_plus_antimuon::~flat_far_muon_plus_antimuon()
  {

    // does nothing.

  }


  void
  flat_far_muon_plus_antimuon::calculateCoefficients()
  {

    double mlSq = ( thirdMass * thirdMass );
    double mXSq = ( fourthMass * fourthMass );

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


    maximumEnergy
    = ( ( firstMass / thirdMass ) * productionFrameEnergy );


    // now we set up the coefficients of various terms, & use them to determine
    // the normalization:

    // MIN to LK segment:
    MINtoLK_segment.setSegmentRange( minimumEnergy,
                                       Elk );

    MINtoLK_const->setCoefficient( -log( minimumEnergy ) );

    MINtoLK_log->setCoefficient( 1.0 );


    // LK to HK segment:
    LKtoHK_segment.setSegmentRange( Elk,
                                      Ehk );

    LKtoHK_const->setCoefficient( log( ( maximumEnergy / Ehk ) ) );


    // HK to MAX segment:
    HKtoMAX_segment.setSegmentRange( Ehk,
                                       maximumEnergy );

    HKtoMAX_const->setCoefficient( log( maximumEnergy ) );

    HKtoMAX_log->setCoefficient( -1.0 );

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



  same_chirality_far_muon::same_chirality_far_muon(
                                    readierForNewPoint* const given_readier,
                                 CppSLHA::CppSLHA0 const* const given_spectrum,
              CppSLHA::particle_property_set const* const given_first_particle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const given_second_particle,
              CppSLHA::particle_property_set const* const given_third_particle,
          CppSLHA::particle_property_set const* const given_fourth_particle ) :
    leptonEnergyDistribution( given_readier,
                                given_spectrum,
                                given_first_particle,
                                effectiveSquarkMass,
                                given_second_particle,
                                given_third_particle,
                                given_fourth_particle ),
    MINtoLK_inv( MINtoLK_segment.addTerm( -1,
                                          0 ) ),
    MINtoLK_const( MINtoLK_segment.addTerm( 0,
                                            0 ) ),
    MINtoLK_log( MINtoLK_segment.addTerm( 0,
                                          1 ) ),
    MINtoLK_lin( MINtoLK_segment.addTerm( 1,
                                          0 ) ),
    LKtoHK_inv( LKtoHK_segment.addTerm( -1,
                                        0 ) ),
    LKtoHK_const( LKtoHK_segment.addTerm( 0,
                                          0 ) ),
    HKtoMAX_inv( HKtoMAX_segment.addTerm( -1,
                                          0 ) ),
    HKtoMAX_const( HKtoMAX_segment.addTerm( 0,
                                            0 ) ),
    HKtoMAX_log( HKtoMAX_segment.addTerm( 0,
                                          1 ) ),
    HKtoMAX_lin( HKtoMAX_segment.addTerm( 1,
                                          0 ) )
  {
    segments.push_back( &MINtoLK_segment );
    segments.push_back( &LKtoHK_segment );
    segments.push_back( &HKtoMAX_segment );
  }

  same_chirality_far_muon::~same_chirality_far_muon()
  {

    // does nothing.

  }


  void
  same_chirality_far_muon::calculateCoefficients()
  {

    double mq = firstMass;
    double mqSq = ( mq * mq );
    double mNSq = ( secondMass * secondMass );
    double mlSq = ( thirdMass * thirdMass );
    double mXSq = ( fourthMass * fourthMass );

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


    // now we set up the coefficients of various terms, & use them to determine
    // the normalization:

    // MIN to LK segment:
    MINtoLK_segment.setSegmentRange( minimumEnergy,
                                       Elk );

    MINtoLK_lin->setCoefficient(
        ( ( 4.0 * mlSq * mNSq * mq )
          / ( ( mlSq - mNSq ) * ( mNSq + mqSq ) * ( mlSq - mXSq ) ) ) );

    MINtoLK_log->setCoefficient(
        ( ( 2.0 * mNSq * ( mqSq - mlSq ) )
          / ( ( mNSq - mlSq ) * ( mNSq + mqSq ) ) ) );

    MINtoLK_const->setCoefficient(
        ( ( 2.0 * mNSq
            * ( 2.0 * mNSq - mqSq + mlSq
              + ( mqSq - mlSq ) * log( ( ( 2.0 * mq ) / ( mlSq - mXSq ) ) ) ) )
        / ( ( mNSq - mlSq ) * ( mNSq + mqSq ) ) ) );

    MINtoLK_inv->setCoefficient(
        ( ( mNSq * ( 2.0 * mNSq - mqSq ) * ( mlSq - mXSq ) )
          / ( mq * ( mlSq - mNSq ) * ( mNSq + mqSq ) ) ) );


    // LK to HK segment:
    LKtoHK_segment.setSegmentRange( Elk,
                                      Ehk );

    // there are 2 kinematic cases for determining the kink energies in terms
    // of the cascade masses.
    if( ( firstMass * thirdMass ) < ( secondMass * secondMass ) )
      // if the assumption that the mass spectrum was SPS1a-like was
      // incorrect...
      {

        LKtoHK_const->setCoefficient(
            ( ( 2.0 * ( ( mlSq + mNSq ) * ( mNSq - mqSq )
                      + mNSq * ( mlSq - mqSq ) * log( ( mNSq / mqSq ) ) ) )
              / ( ( mNSq - mlSq ) * ( mNSq + mqSq ) ) ) );

        LKtoHK_inv->setCoefficient(
            ( ( 2.0 * mNSq * ( mNSq - mqSq ) * ( mlSq - mXSq ) )
              / ( mq * ( mlSq - mNSq ) * ( mNSq + mqSq ) ) ) );

      }
    else
      {

        LKtoHK_const->setCoefficient(
            ( ( 2.0 * ( ( mqSq - mNSq ) * ( mlSq - mNSq )
                        + mNSq * ( mlSq - mqSq ) * log( ( mlSq / mNSq ) ) ) )
              / ( ( mNSq - mlSq ) * ( mNSq + mqSq ) ) ) );

        LKtoHK_inv->setCoefficient( 0.0 );

      }


    // HK to MAX segment:
    HKtoMAX_segment.setSegmentRange( Ehk,
                                       maximumEnergy );

    HKtoMAX_lin->setCoefficient( -MINtoLK_lin->getCoefficient() );

    HKtoMAX_log->setCoefficient( -MINtoLK_log->getCoefficient() );

    HKtoMAX_const->setCoefficient(
        ( ( 2.0 * mNSq * ( mlSq + mqSq
                           + ( mqSq - mlSq )
                             * log( ( ( 2.0 * mlSq )
                                      / ( ( mlSq - mXSq ) * mq ) ) ) ) )
          / ( ( mlSq - mNSq ) * ( mNSq + mqSq ) ) ) );

    HKtoMAX_inv->setCoefficient(
        ( ( mNSq * mq * ( mlSq - mXSq ) )
          / ( ( mNSq - mlSq ) * ( mNSq + mqSq ) ) ) );

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: lepton_distribution_same_chirality_far_muon:";
    std::cout << std::endl << "normalization = " << normalization;
    std::cout << std::endl << "minimum_energy = " << minimum_energy;
    std::cout << std::endl << "Elk = " << Elk;
    std::cout << std::endl << "Ehk = " << Ehk;
    std::cout << std::endl << "maximum_energy = " << maximum_energy;
    std::cout << std::endl << "MINtoLK_lin = " << MINtoLK_lin;
    std::cout << std::endl << "MINtoLK_log = " << MINtoLK_log;
    std::cout << std::endl << "MINtoLK_const = " << MINtoLK_const;
    std::cout << std::endl << "MINtoLK_inv = " << MINtoLK_inv;
    std::cout << std::endl << "LKtoHK_const = " << LKtoHK_const;
    std::cout << std::endl << "LKtoHK_inv = " << LKtoHK_inv;
    std::cout << std::endl;**/
  }



  opposite_chirality_far_muon::opposite_chirality_far_muon(
                                    readierForNewPoint* const given_readier,
                                 CppSLHA::CppSLHA0 const* const given_spectrum,
              CppSLHA::particle_property_set const* const given_first_particle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const given_second_particle,
              CppSLHA::particle_property_set const* const given_third_particle,
          CppSLHA::particle_property_set const* const given_fourth_particle ) :
    leptonEnergyDistribution( given_readier,
                                given_spectrum,
                                given_first_particle,
                                effectiveSquarkMass,
                                given_second_particle,
                                given_third_particle,
                                given_fourth_particle ),
    MINtoLK_inv( MINtoLK_segment.addTerm( -1,
                                          0 ) ),
    MINtoLK_const( MINtoLK_segment.addTerm( 0,
                                            0 ) ),
    MINtoLK_log( MINtoLK_segment.addTerm( 0,
                                          1 ) ),
    MINtoLK_lin( MINtoLK_segment.addTerm( 1,
                                          0 ) ),
    LKtoHK_inv( LKtoHK_segment.addTerm( -1,
                                        0 ) ),
    LKtoHK_const( LKtoHK_segment.addTerm( 0,
                         0 ) ),
    HKtoMAX_inv( HKtoMAX_segment.addTerm( -1,
                                          0 ) ),
    HKtoMAX_const( HKtoMAX_segment.addTerm( 0,
                                            0 ) ),
    HKtoMAX_log( HKtoMAX_segment.addTerm( 0,
                                          1 ) ),
    HKtoMAX_lin( HKtoMAX_segment.addTerm( 1,
                                          0 ) )
  {
    segments.push_back( &MINtoLK_segment );
    segments.push_back( &LKtoHK_segment );
    segments.push_back( &HKtoMAX_segment );
  }

  opposite_chirality_far_muon::~opposite_chirality_far_muon()
  {

    // does nothing.

  }


  void
  opposite_chirality_far_muon::calculateCoefficients()
  {

    double mq = firstMass;
    double mqSq = ( mq * mq );
    double mNSq = ( secondMass * secondMass );
    double mlSq = ( thirdMass * thirdMass );
    double mXSq = ( fourthMass * fourthMass );

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


    // now we set up the coefficients of various terms, & use them to determine
    // the normalization.

    // MIN to LK segment:
    MINtoLK_segment.setSegmentRange( minimumEnergy,
                                       Elk );

    MINtoLK_lin->setCoefficient(
        ( ( 4.0 * mlSq * mNSq * mq )
          / ( ( mNSq - mlSq ) * ( mNSq + mqSq ) * ( mlSq - mXSq ) ) ) );

    MINtoLK_log->setCoefficient(
        ( ( 2.0 * ( mNSq * mNSq - mlSq * mqSq ) )
          / ( ( mNSq - mlSq ) * ( mNSq + mqSq ) ) ) );

    MINtoLK_const->setCoefficient(
        ( ( 2.0 * ( mNSq * ( 2.0 * mNSq - mqSq + mlSq )
                  + ( mlSq * mqSq - mNSq * mNSq )
                    * log( ( ( 2.0 * mq ) / ( mlSq - mXSq ) ) ) ) )
          / ( ( mlSq - mNSq ) * ( mNSq + mqSq ) ) ) );

    MINtoLK_inv->setCoefficient(
        ( ( mNSq * ( ( 2.0 * mNSq - mqSq ) * ( mlSq - mXSq ) ) )
          / ( mq * ( mNSq - mlSq ) * ( mNSq + mqSq ) ) ) );


    // LK to HK segment:
    LKtoHK_segment.setSegmentRange( Elk,
                                      Ehk );

    // there are 2 kinematic cases for determining the kink energies in terms
    // of the cascade masses.
    if( ( firstMass * thirdMass ) < ( secondMass * secondMass ) )
      // if the assumption that the mass spectrum was SPS1a-like was
      // incorrect...
      {

        LKtoHK_const->setCoefficient(
            ( ( 2.0 * ( ( mlSq + mNSq ) * ( mNSq - mqSq )
                        + ( mNSq * mNSq - mlSq * mqSq )
                          * log( ( mNSq / mqSq ) ) ) )
              / ( ( mlSq - mNSq ) * ( mNSq + mqSq ) ) ) );

        LKtoHK_inv->setCoefficient(
            ( ( 2.0 * mNSq * ( mNSq - mqSq ) * ( mlSq - mXSq ) )
              / ( mq * ( mNSq - mlSq ) * ( mNSq + mqSq ) ) ) );

      }
    else
      {

        LKtoHK_const->setCoefficient(
            ( ( 2.0 * ( ( mqSq - mNSq ) * ( mlSq - mNSq )
                        + ( mNSq * mNSq - mlSq * mqSq )
                          * log( ( mlSq / mNSq ) ) ) )
              / ( ( mlSq - mNSq ) * ( mNSq + mqSq ) ) ) );

        LKtoHK_inv->setCoefficient( 0.0 );

      }


    // HK to MAX segment:
    HKtoMAX_segment.setSegmentRange( Ehk,
                                       maximumEnergy );

    HKtoMAX_lin->setCoefficient( -MINtoLK_lin->getCoefficient() );

    HKtoMAX_log->setCoefficient( -MINtoLK_log->getCoefficient() );

    HKtoMAX_const->setCoefficient(
        ( ( 2.0 * ( ( mNSq * mNSq - mlSq * mqSq )
                    * log( ( ( 2.0 * mlSq ) / ( ( mlSq - mXSq ) * mq ) ) )
                    - mNSq * ( mqSq + mlSq ) ) )
          / ( ( mlSq - mNSq ) * ( mNSq + mqSq ) ) ) );

    HKtoMAX_inv->setCoefficient(
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
    std::cout << std::endl << "LKtoHK_const = " << LKtoHK_const;
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



  Higgs_muon_plus_antimuon::Higgs_muon_plus_antimuon(
                                    readierForNewPoint* const given_readier,
                                 CppSLHA::CppSLHA0 const* const given_spectrum,
              CppSLHA::particle_property_set const* const given_first_particle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const given_second_particle,
              CppSLHA::particle_property_set const* const given_third_particle,
          CppSLHA::particle_property_set const* const given_fourth_particle ) :
    leptonEnergyDistribution( given_readier,
                                given_spectrum,
                                given_first_particle,
                                effectiveSquarkMass,
                                given_second_particle,
                                given_third_particle,
                                given_fourth_particle ),
    MINtoLK_const( MINtoLK_segment.addTerm( 0,
                                            0 ) ),
    MINtoLK_log( MINtoLK_segment.addTerm( 0,
                                          1 ) ),
    LKtoHK_const( LKtoHK_segment.addTerm( 0,
                                          0 ) ),
    HKtoMAX_const( HKtoMAX_segment.addTerm( 0,
                                            0 ) ),
    HKtoMAX_log( HKtoMAX_segment.addTerm( 0,
                                          1 ) )
  {
    segments.push_back( &MINtoLK_segment );
    segments.push_back( &LKtoHK_segment );
    segments.push_back( &HKtoMAX_segment );
  }

  Higgs_muon_plus_antimuon::~Higgs_muon_plus_antimuon()
  {

    // does nothing.

  }


  void
  Higgs_muon_plus_antimuon::calculateCoefficients()
  {

    double mq = firstMass;
    double mN = secondMass;
    double mNSq = ( mN * mN );
    double mH = thirdMass;
    double mHSq = ( mH * mH );
    double mX = fourthMass;
    double mXSq = ( mX * mX );

    productionFrameEnergy = ( 0.5 * mH );

    double EHN = ( ( mNSq + mHSq - mXSq ) / ( 2.0 * mN ) );
    double gammaNH = ( EHN / mH );
    double betaNH = sqrt( ( 1.0 - ( 1.0 / ( gammaNH * gammaNH ) ) ) );
    double EmaxN = ( gammaNH * ( 1.0 + betaNH ) * productionFrameEnergy );
    double EminN = ( mHSq / ( 4.0 * EmaxN ) );
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


    // now we set up the coefficients of various terms, & use them to determine
    // the normalization:

    // MIN to LK segment:
    MINtoLK_segment.setSegmentRange( minimumEnergy,
                                       Elk );

    MINtoLK_const->setCoefficient( -log( minimumEnergy ) );
    MINtoLK_log->setCoefficient( 1.0 );



    // LK to HK segment:
    LKtoHK_segment.setSegmentRange( Elk,
                                      Ehk );
    LKtoHK_const->setCoefficient( log( ( maximumEnergy / Ehk ) ) );


    // HK to MAX segment:
    HKtoMAX_segment.setSegmentRange( Ehk,
                                       maximumEnergy );

    HKtoMAX_const->setCoefficient( log( maximumEnergy ) );
    HKtoMAX_log->setCoefficient( -1.0 );

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: lepton_distribution_Higgs_muon_plus_antimuon_distribution:";
    std::cout << std::endl << "normalization = " << normalization;
    std::cout << std::endl << "minimum_energy = " << minimum_energy;
    std::cout << std::endl << "Elk = " << Elk;
    std::cout << std::endl << "Ehk = " << Ehk;
    std::cout << std::endl << "maximum_energy = " << maximum_energy;
    std::cout << std::endl << "MINtoLK_log = " << MINtoLK_log;
    std::cout << std::endl << "MINtoLK_const = " << MINtoLK_const;
    std::cout << std::endl << "LKtoHK_const = " << LKtoHK_const;
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
    MINtoLK_inv( MINtoLK_segment.addTerm( -1,
                                          0 ) ),
    MINtoLK_const( MINtoLK_segment.addTerm( 0,
                                            0 ) ),
    MINtoLK_log( MINtoLK_segment.addTerm( 0,
                                          1 ) ),
    MINtoLK_lin( MINtoLK_segment.addTerm( 1,
                                          0 ) ),
    MINtoLK_linlog( MINtoLK_segment.addTerm( 1,
                                             1 ) ),
    MINtoLK_sq( MINtoLK_segment.addTerm( 2,
                                         0 ) ),
    LKtoHK_const( LKtoHK_segment.addTerm( 0,
                                          0 ) ),
    LKtoHK_lin( LKtoHK_segment.addTerm( 1,
                                        0 ) ),
    LKtoHK_sq( LKtoHK_segment.addTerm( 2,
                                       0 ) ),
    HKtoMAX_inv( HKtoMAX_segment.addTerm( -1,
                                          0 ) ),
    HKtoMAX_const( HKtoMAX_segment.addTerm( 0,
                                            0 ) ),
    HKtoMAX_log( HKtoMAX_segment.addTerm( 0,
                                          1 ) ),
    HKtoMAX_lin( HKtoMAX_segment.addTerm( 1,
                                          0 ) ),
    HKtoMAX_linlog( HKtoMAX_segment.addTerm( 1,
                                             1 ) ),
    HKtoMAX_sq( HKtoMAX_segment.addTerm( 2,
                                         0 ) )
  {
    segments.push_back( &MINtoLK_segment );
    segments.push_back( &LKtoHK_segment );
    segments.push_back( &HKtoMAX_segment );
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
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "cosinesLimitedByEnergy = true";
      std::cout << std::endl;/**/
    }
    else
    {
      cosinesLimitedByEnergy = false;
      // for comparison with the Mathematica notebook: underlap is Ehk > Elk
      // => if( cosinesLimitedByEnergy ){ [use overlap] }

      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "cosinesLimitedByEnergy = false";
      std::cout << std::endl;/**/
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

    MINtoLK_segment.setSegmentRange( minimumEnergy,
                                     Elk );
    LKtoHK_segment.setSegmentRange( Elk,
                                    Ehk );
    HKtoMAX_segment.setSegmentRange( Ehk,
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


    MINtoLK_inv->setCoefficient( 0.0 );
    MINtoLK_const->setCoefficient( currentMinToLkConstCoefficient );
    MINtoLK_log->setCoefficient( currentMinToLkLogCoefficient );
    MINtoLK_lin->setCoefficient( currentMinToLkLinCoefficient );
    MINtoLK_linlog->setCoefficient( 0.0 );
    MINtoLK_sq->setCoefficient( currentMinToLkSqCoefficient );

    LKtoHK_const->setCoefficient( currentLkToHkConstCoefficient );
    LKtoHK_lin->setCoefficient( currentLkToHkLinCoefficient );
    LKtoHK_sq->setCoefficient( currentLkToHkSqCoefficient );

    HKtoMAX_inv->setCoefficient( 0.0 );
    HKtoMAX_const->setCoefficient( currentHkToMaxConstCoefficient );
    HKtoMAX_log->setCoefficient( currentHkToMaxLogCoefficient );
    HKtoMAX_lin->setCoefficient( currentHkToMaxLinCoefficient );
    HKtoMAX_linlog->setCoefficient( 0.0 );
    HKtoMAX_sq->setCoefficient( currentHkToMaxSqCoefficient );

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "VvPlusJjAaAllSqSymDistribution:"
    << std::endl << "MINtoLK: "
    << std::endl << MINtoLK_inv->getCoefficient() << "/E"
    << std::endl << MINtoLK_const->getCoefficient()
    << std::endl << MINtoLK_log->getCoefficient() << "L"
    << std::endl << MINtoLK_lin->getCoefficient() << "E"
    << std::endl << MINtoLK_linlog->getCoefficient() << "EL"
    << std::endl << MINtoLK_sq->getCoefficient() << "E^2"
    << std::endl << "LKtoHK: "
    << std::endl << LKtoHK_const->getCoefficient()
    << std::endl << LKtoHK_lin->getCoefficient() << "E"
    << std::endl << LKtoHK_sq->getCoefficient() << "E^2"
    << std::endl << "HKtoMAX: "
    << std::endl << HKtoMAX_inv->getCoefficient() << "/E"
    << std::endl << HKtoMAX_const->getCoefficient()
    << std::endl << HKtoMAX_log->getCoefficient() << "L"
    << std::endl << HKtoMAX_lin->getCoefficient() << "E"
    << std::endl << HKtoMAX_linlog->getCoefficient() << "EL"
    << std::endl << HKtoMAX_sq->getCoefficient() << "E^2";
    std::cout << std::endl;/**/


    // VvSqMinusAaSqSymDistribution:

    // MIN to LK:

    // debugging:
    currentMinToLkInvCoefficient = 0.0;
    // right!

    currentMinToLkConstCoefficient
    = ( EightmCXElMinmQsq * ElMin * ( lnmQC - lnElMin ) * mQCsqDiff * mVsq );
    // right!

    currentMinToLkLogCoefficient
    = ( EightmCXElMinmQsq * ElMin * mQCsqDiff * mVsq );
    // right!

    // debugging:
    currentMinToLkLinCoefficient = 0.0;
    // right!

    // debugging:
    currentMinToLkLinLogCoefficient = 0.0;
    // right!

    // debugging:
    currentMinToLkSqCoefficient = 0.0;
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
    currentLkToHkLinCoefficient = 0.0;
    // right!

    // debugging:
    currentLkToHkSqCoefficient = 0.0;
    // right!


    // HK to MAX:

    // debugging:
    currentHkToMaxInvCoefficient = 0.0;
    // right!

    currentHkToMaxConstCoefficient
    = ( EightmCXElMinmQsq * ElMin * mQCsqDiff * mVsq
        * ( lnmQC + lnVsqOverFourElMinsq + lnElMin ) );
    // right!

    currentHkToMaxLogCoefficient = -currentMinToLkLogCoefficient;
    // right!

    // debugging:
    currentHkToMaxLinCoefficient = 0.0;
    // right!

    // debugging:
    currentHkToMaxLinLogCoefficient = 0.0;
    // right!

    // debugging:
    currentHkToMaxSqCoefficient = 0.0;
    // right!


    if( couplesAsVectorNotAxial )
    {
      MINtoLK_const->addToCoefficient( currentMinToLkConstCoefficient );
      MINtoLK_log->addToCoefficient( currentMinToLkLogCoefficient );

      LKtoHK_const->addToCoefficient( currentLkToHkConstCoefficient );

      HKtoMAX_const->addToCoefficient( currentHkToMaxConstCoefficient );
      HKtoMAX_log->addToCoefficient( currentHkToMaxLogCoefficient );
    }
    else
    {
      MINtoLK_const->addToCoefficient( -currentMinToLkConstCoefficient );
      MINtoLK_log->addToCoefficient( -currentMinToLkLogCoefficient );

      LKtoHK_const->addToCoefficient( -currentLkToHkConstCoefficient );

      HKtoMAX_const->addToCoefficient( -currentHkToMaxConstCoefficient );
      HKtoMAX_log->addToCoefficient( -currentHkToMaxLogCoefficient );
    }

    // debugging:
    /**/std::cout << std::endl << "debugging:"
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
    << std::endl << currentHkToMaxSqCoefficient << "E^2";/**/


    if( shouldNotSumOverHandedness )
    {
      // VvPlusJjAaAllSqAntiDistribution:

      // MIN to LK:

      // debugging:
      currentMinToLkInvCoefficient = 0.0;
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
      currentMinToLkLinLogCoefficient = 0.0;
      // right!

      // debugging:
      currentMinToLkSqCoefficient = 0.0;
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
      currentLkToHkSqCoefficient = 0.0;


      // HK to MAX:

      currentHkToMaxInvCoefficient = 0.0;
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
      currentHkToMaxLinLogCoefficient = 0.0;
      // right!

      // debugging:
      currentHkToMaxSqCoefficient = 0.0;
      // right!


      MINtoLK_const->addToCoefficient( currentMinToLkConstCoefficient );
      MINtoLK_log->addToCoefficient( currentMinToLkLogCoefficient );
      MINtoLK_lin->addToCoefficient( currentMinToLkLinCoefficient );

      LKtoHK_const->addToCoefficient( currentLkToHkConstCoefficient );
      LKtoHK_lin->addToCoefficient( currentLkToHkLinCoefficient );

      HKtoMAX_const->addToCoefficient( currentHkToMaxConstCoefficient );
      HKtoMAX_log->addToCoefficient( currentHkToMaxLogCoefficient );
      HKtoMAX_lin->addToCoefficient( currentHkToMaxLinCoefficient );

      // debugging:
      /**/std::cout << std::endl << "debugging:"
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
      << std::endl << currentHkToMaxSqCoefficient << "E^2";/**/


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
      currentMinToLkSqCoefficient = 0.0;
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
      currentLkToHkSqCoefficient = 0.0;
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
      currentHkToMaxSqCoefficient = 0.0;
      // right!


      MINtoLK_inv->addToCoefficient( currentMinToLkInvCoefficient );
      MINtoLK_const->addToCoefficient( currentMinToLkConstCoefficient );
      MINtoLK_log->addToCoefficient( currentMinToLkLogCoefficient );
      MINtoLK_lin->addToCoefficient( currentMinToLkLinCoefficient );
      MINtoLK_linlog->addToCoefficient( currentMinToLkLinLogCoefficient );

      LKtoHK_const->addToCoefficient( currentLkToHkConstCoefficient );
      LKtoHK_lin->addToCoefficient( currentLkToHkLinCoefficient );

      HKtoMAX_inv->addToCoefficient( currentHkToMaxInvCoefficient );
      HKtoMAX_const->addToCoefficient( currentHkToMaxConstCoefficient );
      HKtoMAX_log->addToCoefficient( currentHkToMaxLogCoefficient );
      HKtoMAX_lin->addToCoefficient( currentHkToMaxLinCoefficient );
      HKtoMAX_linlog->addToCoefficient( currentHkToMaxLinLogCoefficient );

      // debugging:
      /**/std::cout << std::endl << "debugging:"
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
      << std::endl << currentHkToMaxSqCoefficient << "E^2";/**/


      // VvSqMinusAaSqAntiDistribution

      // MIN to LK:

      // debugging:
      currentMinToLkInvCoefficient = 0.0;
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
      currentMinToLkLinLogCoefficient = 0.0;
      // right!

      // debugging:
      currentMinToLkSqCoefficient = 0.0;
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
      currentLkToHkSqCoefficient = 0.0;
      // right!


      // HK to MAX:

      currentHkToMaxInvCoefficient = 0.0;
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
      currentHkToMaxLinLogCoefficient = 0.0;
      // right!

      // debugging:
      currentHkToMaxSqCoefficient = 0.0;
      // right!


      if( couplesAsVectorNotAxial )
      {
        MINtoLK_const->addToCoefficient( currentMinToLkConstCoefficient );
        MINtoLK_log->addToCoefficient( currentMinToLkLogCoefficient );
        MINtoLK_lin->addToCoefficient( currentMinToLkLinCoefficient );

        LKtoHK_const->addToCoefficient( currentLkToHkConstCoefficient );
        LKtoHK_lin->addToCoefficient( currentLkToHkLinCoefficient );

        HKtoMAX_const->addToCoefficient( currentHkToMaxConstCoefficient );
        HKtoMAX_log->addToCoefficient( currentHkToMaxLogCoefficient );
        HKtoMAX_lin->addToCoefficient( currentHkToMaxLinCoefficient );
      }
      else
      {
        MINtoLK_const->addToCoefficient( -currentMinToLkConstCoefficient );
        MINtoLK_log->addToCoefficient( -currentMinToLkLogCoefficient );
        MINtoLK_lin->addToCoefficient( -currentMinToLkLinCoefficient );

        LKtoHK_const->addToCoefficient( -currentLkToHkConstCoefficient );
        LKtoHK_lin->addToCoefficient( -currentLkToHkLinCoefficient );

        HKtoMAX_const->addToCoefficient( -currentHkToMaxConstCoefficient );
        HKtoMAX_log->addToCoefficient( -currentHkToMaxLogCoefficient );
        HKtoMAX_lin->addToCoefficient( -currentHkToMaxLinCoefficient );
      }

      // debugging:
      /**/std::cout << std::endl << "debugging:"
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
      << std::endl << currentHkToMaxSqCoefficient << "E^2";/**/
    }

    // debugging:
    std::cout
    << std::endl
    << "debugging: Z_handed_muon:calculateCoefficients() produced "
    << std::endl << "firstMass = " << firstMass
    << std::endl << "secondMass = " << secondMass
    << std::endl << "cosinesLimitedByEnergy = " << cosinesLimitedByEnergy
    << std::endl << "(true = " << true << ")"
    << std::endl << "minimumEnergy = " << minimumEnergy
    << std::endl << "maximumEnergy = " << maximumEnergy
    << std::endl << "MINtoLK_inv = " << MINtoLK_inv->getCoefficient()
    << std::endl << "MINtoLK_const = " << MINtoLK_const->getCoefficient()
    << std::endl << "MINtoLK_log = " << MINtoLK_log->getCoefficient()
    << std::endl << "MINtoLK_lin = " << MINtoLK_lin->getCoefficient()
    << std::endl << "MINtoLK_linlog = " << MINtoLK_linlog->getCoefficient()
    << std::endl << "MINtoLK_sq = " << MINtoLK_sq->getCoefficient()
    << std::endl << "LKtoHK_const = " << LKtoHK_const->getCoefficient()
    << std::endl << "LKtoHK_lin = " << LKtoHK_lin->getCoefficient()
    << std::endl << "LKtoHK_sq = " << LKtoHK_sq->getCoefficient()
    << std::endl << "HKtoMAX_inv = " << HKtoMAX_inv->getCoefficient()
    << std::endl << "HKtoMAX_const = " << HKtoMAX_const->getCoefficient()
    << std::endl << "HKtoMAX_log = " << HKtoMAX_log->getCoefficient()
    << std::endl << "HKtoMAX_lin = " << HKtoMAX_lin->getCoefficient()
    << std::endl << "HKtoMAX_linlog = " << HKtoMAX_linlog->getCoefficient()
    << std::endl << "HKtoMAX_sq = " << HKtoMAX_sq->getCoefficient();
    std::cout << std::endl;/**/
  }



  W_minus_handed_muon::W_minus_handed_muon(
                                    readierForNewPoint* const given_readier,
                                 CppSLHA::CppSLHA0 const* const given_spectrum,
              CppSLHA::particle_property_set const* const given_first_particle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const given_second_particle,
              CppSLHA::particle_property_set const* const given_third_particle,
          CppSLHA::particle_property_set const* const given_fourth_particle,
                                        bool const given_jet_is_left_handed ) :
    leptonEnergyDistribution( given_readier,
                                given_spectrum,
                                given_first_particle,
                                effectiveSquarkMass,
                                given_second_particle,
                                given_third_particle,
                                given_fourth_particle ),
    jet_is_left_handed( given_jet_is_left_handed ),
    MINtoLK_inv( MINtoLK_segment.addTerm( -1,
                                          0 ) ),
    MINtoLK_const( MINtoLK_segment.addTerm( 0,
                                            0 ) ),
    MINtoLK_log( MINtoLK_segment.addTerm( 0,
                                          1 ) ),
    MINtoLK_lin( MINtoLK_segment.addTerm( 1,
                                          0 ) ),
    MINtoLK_linlog( MINtoLK_segment.addTerm( 1,
                                             1 ) ),
    MINtoLK_sq( MINtoLK_segment.addTerm( 2,
                                         0 ) ),
    LKtoHK_const( LKtoHK_segment.addTerm( 0,
                                          0 ) ),
    LKtoHK_lin( LKtoHK_segment.addTerm( 1,
                                        0 ) ),
    LKtoHK_sq( LKtoHK_segment.addTerm( 2,
                                       0 ) ),
    HKtoMAX_inv( HKtoMAX_segment.addTerm( -1,
                                          0 ) ),
    HKtoMAX_const( HKtoMAX_segment.addTerm( 0,
                                            0 ) ),
    HKtoMAX_log( HKtoMAX_segment.addTerm( 0,
                                          1 ) ),
    HKtoMAX_lin( HKtoMAX_segment.addTerm( 1,
                                          0 ) ),
    HKtoMAX_linlog( HKtoMAX_segment.addTerm( 1,
                                             1 ) ),
    HKtoMAX_sq( HKtoMAX_segment.addTerm( 2,
                                         0 ) )
  {
    NMIX = spectrum->inspect_BLOCK( "NMIX" );
    UMIX = spectrum->inspect_BLOCK( "UMIX" );
    VMIX = spectrum->inspect_BLOCK( "VMIX" );
    if( CppSLHA::PDG_code::chargino_one
        == given_second_particle->get_PDG_code() )
    {
      which_chargino = 1;
    }
    else if( CppSLHA::PDG_code::chargino_two
        == given_second_particle->get_PDG_code() )
    {
      which_chargino = 2;
    }
    else
    {
      which_chargino = -1;
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
    if( CppSLHA::PDG_code::neutralino_one
        == given_fourth_particle->get_PDG_code() )
    {
      which_neutralino = 1;
    }
    else if( CppSLHA::PDG_code::neutralino_two
        == given_fourth_particle->get_PDG_code() )
    {
      which_neutralino = 2;
    }
    else if( CppSLHA::PDG_code::neutralino_three
        == given_fourth_particle->get_PDG_code() )
    {
      which_neutralino = 3;
    }
    else if( CppSLHA::PDG_code::neutralino_four
        == given_fourth_particle->get_PDG_code() )
    {
      which_neutralino = 4;
    }
    else
    {
      which_neutralino = -1;
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
    segments.push_back( &MINtoLK_segment );
    segments.push_back( &LKtoHK_segment );
    segments.push_back( &HKtoMAX_segment );
  }

  W_minus_handed_muon::~W_minus_handed_muon()
  {

    // does nothing.

  }


  void
  W_minus_handed_muon::calculateCoefficients()
  {

    /* it doesn't matter if the neutralino mass is negative, since only the
     * absolute square appears, & the couplings are only affected by it as an
     * overall phase.  I (BOL) assume no CP violation, so the neutralino mixing
     * matrices are strictly real.
     */

    Vv = ( NMIX->get_entry( which_neutralino,
                            2 ) * ( UMIX->get_entry( which_chargino,
                                                     1 )
                                    + VMIX->get_entry( which_chargino,
                                                       1 ) )
           + M_SQRT1_2 * ( NMIX->get_entry( which_neutralino,
                                            3 )
                           * UMIX->get_entry( which_chargino,
                                              2 )
                           - NMIX->get_entry( which_neutralino,
                                              4 )
                             * VMIX->get_entry( which_chargino,
                                                2 ) ) );
    Aa = ( NMIX->get_entry( which_neutralino,
                            2 ) * ( VMIX->get_entry( which_chargino,
                                                     1 )
                                    - UMIX->get_entry( which_chargino,
                                                       1 ) )
           - M_SQRT1_2 * ( NMIX->get_entry( which_neutralino,
                                            3 )
                           * UMIX->get_entry( which_chargino,
                                              2 )
                           + NMIX->get_entry( which_neutralino,
                                              4 )
                             * VMIX->get_entry( which_chargino,
                                                2 ) ) );

    if( 0.0 > fourthParticle->get_mass() )
      // if the relative phases end up such that the product of couplings is
      // purely imaginary...
      {

        // we swap the coupings (borrowing current_term_coefficient for a
        // moment):
        current_term_coefficient = Vv;
        Vv = Aa;
        Aa = current_term_coefficient;

      }

    Vvsq = ( Vv * Vv );
    Aasq = ( Aa * Aa );
    VvSqMinusAaSq = ( Vvsq - Aasq );

    if( jet_is_left_handed )
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

        // we swap the energies (borrowing current_term_coefficient for a
        // moment) & note that we had to swap them:
        current_term_coefficient = Elk;
        Elk = Ehk;
        Ehk = current_term_coefficient;
        energy_limited_cosines = true;

      }
    else
      {

        energy_limited_cosines = false;

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





    // TwiceVvAaSymDistribution:
    current_term_coefficient
    = ( -4.0 * ElMinsq * mCsq * mQcu * mVsq * mVsq );
    MINtoLK_inv->setCoefficient( current_term_coefficient );
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "TwiceVvAaSymDistribution:"
    << std::endl << "MINtoLK: "
    << std::endl << current_term_coefficient << "/E";

    current_term_coefficient
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
    MINtoLK_const->setCoefficient( current_term_coefficient );
    // debugging:
    /**/std::cout
    << std::endl << current_term_coefficient;

    current_term_coefficient
    = ( -2.0 * ElMin * mQsq * ( EightmCEVElMin * mQsq + 4.0 * mCsq * mQCsqSum )
        * mVsq );
    MINtoLK_log->setCoefficient( current_term_coefficient );
    // debugging:
    /**/std::cout
    << std::endl << current_term_coefficient << "L";

    current_term_coefficient
    = ( 2.0 * mQcu
        * ( 4.0 * EightmCEVElMin * ElMin
            * ( ( lnElMin - lnmQC ) * mCsq + mQsq )
            + EightmCEVElMin * mC * mVsq
            - 2.0 * mCsq * ( 4.0 * ElMin * mC * ( 4.0 * ElMinsq - mVsq )
                             + mVsq * mVsq) ) );
    MINtoLK_lin->setCoefficient( current_term_coefficient );
    // debugging:
    /**/std::cout
    << std::endl << current_term_coefficient << "E";

    current_term_coefficient
    = ( -8.0 * EightmCEVElMin * ElMin * mCsq * mQcu );
    MINtoLK_linlog->setCoefficient( 0.0 );
    // debugging:
    /**/std::cout
    << std::endl << current_term_coefficient << "EL";

    current_term_coefficient
    = ( -32.0 * ElMinsq * mQsq * mQsq * ( mQCsqDiff - mCsq ) );
    MINtoLK_sq->setCoefficient( current_term_coefficient );
    // debugging:
    /**/std::cout
    << std::endl << current_term_coefficient << "E^2";

    if( energy_limited_cosines )
    {
      current_term_coefficient
      = ( 2.0 * mQsq
          * ( EightmCEVElMin * mQsq
              * ( mC * ( mVsq - 4.0 * ElMinsq )
                  - ElMin * mVsq * lnVsqOverFourElMinsq )
              + mC * ( 16.0 * ElMinsq * ElMinsq * mC * mQsq
                       - 4.0 * ElMinsq * mQCsqSum * mVsq
                         * ( 2.0 * ElMin
                             + lnVsqOverFourElMinsq * mC )
                       - ( mC * mQsq
                           - 2.0 * ElMin * mQCsqSum ) * mVsq * mVsq ) ) );
      LKtoHK_const->setCoefficient( current_term_coefficient );
      // debugging:
      /**/std::cout
      << std::endl << "LKtoHK: "
      << std::endl << current_term_coefficient;
      current_term_coefficient
      = ( 2.0 * mC * mQcu
          * ( EightmCEVElMin
              * ( mVsq
                  - 4.0 * ElMin * ( lnVsqOverFourElMinsq * mC + ElMin ) )
              + 2.0 * mC * ( 8.0 * ElMin
                             * ( 2.0 * ElMinsq * ( ElMin - 2.0 * mC )
                                 + mC * mVsq )
                             - mVsq * mVsq ) ) );
      LKtoHK_lin->setCoefficient( current_term_coefficient );
      // debugging:
      /**/std::cout
      << std::endl << current_term_coefficient << "E";
      current_term_coefficient = 0.0;
      LKtoHK_sq->setCoefficient( 0.0 );
      // debugging:
      /**/std::cout
      << std::endl << current_term_coefficient << "E^2";
    }
    else
    {
      current_term_coefficient
      = ( 2.0 * ElMin * mQsq * mVsq
          * ( mQCsqDiff
              * ( EightmCEVElMin + 8.0 * ElMin * mCsq )
              - 8.0 * lnmQC * ( 4.0 * mQCsqSum * mCsq * ElMin
                                + mQsq * EightmCEVElMin ) ) );
      LKtoHK_const->setCoefficient( current_term_coefficient );
      // debugging:
      /**/std::cout
      << std::endl << "LKtoHK: "
      << std::endl << current_term_coefficient;

      current_term_coefficient
      = ( 8.0 * EightmCEVElMin * ElMin * mQcu
          * ( mQCsqDiff - 2.0 * lnmQC * mCsq ) );
      LKtoHK_lin->setCoefficient( current_term_coefficient );
      // debugging:
      /**/std::cout
      << std::endl << current_term_coefficient << "E";

      current_term_coefficient
      = ( -32.0 * ElMinsq * mQsq * mQCsqDiff * mQCsqDiff );
      LKtoHK_sq->setCoefficient( current_term_coefficient );
      // debugging:
      /**/std::cout
      << std::endl << current_term_coefficient << "E^2";
    }

    current_term_coefficient = ( 4.0 * ElMinsq * mCsq * mQcu * mVsq * mVsq );
    HKtoMAX_inv->setCoefficient( current_term_coefficient );
    // debugging:
    /**/std::cout
    << std::endl << "HKtoMAX: "
    << std::endl << current_term_coefficient << "/E";

    current_term_coefficient
    = ( 2.0 * mQsq * mVsq
        * ( EightmCEVElMin * ElMin * mQsq
            * ( 1.0 - lnmQC - lnVsqOverFourElMinsq - lnElMin )
            - 4.0 * ElMinsq * mCsq
              * ( ( lnmQC + lnVsqOverFourElMinsq + lnElMin ) * mQCsqSum
                  - 2.0 * mQsq )
            - 8.0 * ElMin * ElMinsq * mC * mQCsqSum
            + mC * mQsq * ( EightmCEVElMin - mC * mVsq ) ) );
    HKtoMAX_const->setCoefficient( current_term_coefficient );
    // debugging:
    /**/std::cout
    << std::endl << current_term_coefficient;

    current_term_coefficient
    = ( 2.0 * ElMin * mQsq * mVsq
        * ( EightmCEVElMin * mQsq + 4.0 * ElMin * mCsq * mQCsqSum ) );
    HKtoMAX_log->setCoefficient( current_term_coefficient );
    // debugging:
    /**/std::cout
    << std::endl << current_term_coefficient << "L";

    current_term_coefficient
    = ( 8.0  * ElMin * mC * mQcu
        * ( 2.0 * mC * ( 4.0 * ElMinsq * ( ElMin - mC ) + mC * mVsq )
            - EightmCEVElMin
            * ( ElMin + ( 1.0 + lnmQC + lnVsqOverFourElMinsq ) * mC ) ) );
    HKtoMAX_lin->setCoefficient( current_term_coefficient );
    // debugging:
    /**/std::cout
    << std::endl << current_term_coefficient << "E";

    current_term_coefficient = ( 8.0 * EightmCEVElMin * ElMin * mCsq * mQcu );
    HKtoMAX_linlog->setCoefficient( current_term_coefficient );
    // debugging:
    /**/std::cout
    << std::endl << current_term_coefficient << "EL";

    current_term_coefficient = ( -32.0 * ElMinsq * mCsq * mCsq * mQsq );
    HKtoMAX_sq->setCoefficient( current_term_coefficient );
    // debugging:
    /**/std::cout
    << std::endl << current_term_coefficient << "E^2";
    std::cout << std::endl;/**/





















    // MIN to LK:
    MINtoLK_segment.setSegmentRange( minimumEnergy,
                                       Elk );

    // inv:

    // VvPlusJjAaAllSq
    // 0.0

    // MinusTwiceVvJjAa
    current_term_coefficient
    = ( MinusTwiceVvJjAa * ( -4.0 * ElMinsq * mCsq * mQcu * mVqu ) );

    // VvSqMinusAaSq
    // 0.0

    // JjVvPlusJjAaAllSq
    // 0.0

    // JjVvSqPlusAaSq
    current_term_coefficient
    += ( JjVvSqPlusAaSq
         * ( 4.0 * ElMinsq * mCsq * mQcu * mVsq ) );

    // JjVvSqMinusAaSq
    // 0.0

    MINtoLK_inv->setCoefficient( current_term_coefficient );


    // const:

    // VvPlusJjAaAllSq
    current_term_coefficient
    = ( VvPlusJjAaAllSq
        * ( ElMin * mQsq * mQCsqDiff
            * ( EightmCEVElMin * ( lnElMin - lnmQC ) * mVsq
                - 4.0 * ElMin * mC
                  * ( EightmCEVElMin - 4.0 * ElMinsq * mC ) ) ) );

    // MinusTwiceVvJjAa
    current_term_coefficient
    += ( MinusTwiceVvJjAa
         * ( ElMin * mQsq
             * ( 2.0 * ( lnElMin - lnmQC ) * mVsq
                 * ( EightmCEVElMin * mQsq + 4.0 * ElMin * mCsq * mQCsqSum )
                 - 2.0 * EightmCEVElMin
                   * ( 4.0 * ElMin * mC * mQsq - mCsq * mVsq )
                 + 16.0 * ElMin * mCsq
                   * ( 2.0 * ElMinsq * mQsq - mCsq * mVsq )
                 + 4.0 * mC * mQCsqSum * mVqu ) ) );

    // VvSqMinusAaSq
    current_term_coefficient
    += ( VvPlusJjAaAllSq
         * ( EightmCXElMinmQsq * ElMin * ( lnmQC - lnElMin )
             * mQCsqDiff * mVsq ) );

    // JjVvPlusJjAaAllSq
    current_term_coefficient
    += ( JjVvPlusJjAaAllSq
         * ( ( ElMin * mQsq * mQCsqDiff * mVsq
             * ( EightmCEVElMin * ( lnElMin - lnmQC )
                 - 16.0 * ElMinsq * mC ) ) ) );

    // JjVvSqPlusAaSq
    current_term_coefficient
    += ( JjVvSqPlusAaSq
         * ( ElMin * mQsq * mVsq
             * ( 32.0 * ElMinsq * mC * mQsq
                 + ( EightmCEVElMin * mQsq
                     + 4.0 * ElMin * mCsq * mQCsqSum )
                   * 2.0 * ( lnmQC - lnElMin )
                 + 2.0 * mCsq * ( EightmCEVElMin + 4.0 * ElMin * mCsq )
                 - 4.0 * mC * mQCsqSum * mVsq ) ) );

    // JjVvSqMinusAaSq
    current_term_coefficient
    += ( JjVvSqMinusAaSq
         * ( EightmCXElMinmQsq * ElMin * mVsq
             * ( 2.0 * mCsq + ( lnmQC - lnElMin ) * mQCsqSum ) ) );

    MINtoLK_const->setCoefficient( current_term_coefficient );


    // log:

    // VvPlusJjAaAllSq
    current_term_coefficient
    = ( VvPlusJjAaAllSq
        * ( -EightmCEVElMin * ElMin * mQsq * mQCsqDiff * mVsq ) );

    // MinusTwiceVvJjAa
    current_term_coefficient
    += ( MinusTwiceVvJjAa
         * ( -2.0 * ElMin * mQsq * mVsq
             * ( EightmCEVElMin * mQsq + 4.0 * ElMin * mCsq * mQCsqSum ) ) );

    // VvSqMinusAaSq
    current_term_coefficient
    += ( VvSqMinusAaSq
         * ( EightmCXElMinmQsq * ElMin * mQCsqDiff * mVsq ) );

    // JjVvPlusJjAaAllSq
    current_term_coefficient
    += ( JjVvPlusJjAaAllSq
         * ( -EightmCEVElMin * ElMin * mQsq * mQCsqDiff * mVsq ) );

    // JjVvSqPlusAaSq
    current_term_coefficient
    += ( JjVvSqPlusAaSq
         * ( 2.0 * ElMin * mQsq * mVsq
             * ( EightmCEVElMin * mQsq
                 + 4.0 * ElMin * mCsq * mQCsqSum ) ) );

    // JjVvSqMinusAaSq
    current_term_coefficient
    += ( JjVvSqMinusAaSq
         * ( EightmCXElMinmQsq * ElMin * mQCsqSum * mVsq ) );

    MINtoLK_log->setCoefficient( current_term_coefficient );


    // lin:

    // VvPlusJjAaAllSq
    current_term_coefficient
    = ( VvPlusJjAaAllSq
        * ( 4.0 * EightmCEVElMin * ElMin * mQcu * mQCsqDiff ) );

    // MinusTwiceVvJjAa
    current_term_coefficient
    += ( MinusTwiceVvJjAa
         * ( 2.0 * mQcu
             * ( EightmCEVElMin
                 * ( 4.0 * ElMin * ( ( lnElMin - lnmQC ) * mCsq + mQsq )
                     + mC * mVsq )
                 - 2.0 * mCsq
                   * ( 4.0 * ElMin * mC * ( 4.0 * ElMinsq - mVsq )
                       + mVqu ) ) ) );

    // VvSqMinusAaSq
    // 0.0

    // JjVvPlusJjAaAllSq
    current_term_coefficient
    += ( JjVvPlusJjAaAllSq
         * ( 16.0 * ElMinsq * mQcu * mQCsqDiff * mVsq ) );

    // JjVvSqPlusAaSq
    current_term_coefficient
    += ( JjVvSqPlusAaSq
         * ( 2.0 * mQcu * mVsq
             * ( 2.0 * mCsq * mVsq
                 - EightmCEVElMin * mC
                 - 8.0 * ElMin * ( mCsq * mC
                                   + 2.0 * ElMin * ( ( lnElMin - lnmQC ) * mCsq
                                                     + mQsq ) ) ) ) );

    // JjVvSqMinusAaSq
    current_term_coefficient
    += ( JjVvSqMinusAaSq
         * ( -2.0 * EightmCXElMinmQsq * mC * mQ * mVsq ) );

    MINtoLK_lin->setCoefficient( current_term_coefficient );


    // linlog:

    // VvPlusJjAaAllSq
    // 0.0

    // MinusTwiceVvJjAa
    current_term_coefficient
    = ( MinusTwiceVvJjAa
        * ( -8.0 * EightmCEVElMin * ElMin * mCsq * mQcu ) );

    // VvSqMinusAaSq
    // 0.0

    // JjVvPlusJjAaAllSq
    // 0.0

    // JjVvSqPlusAaSq
    current_term_coefficient
    += ( JjVvSqPlusAaSq
         * ( 32.0 * ElMinsq * mCsq * mQcu * mVsq ) );

    // JjVvSqMinusAaSq
    // 0.0

    MINtoLK_linlog->setCoefficient( current_term_coefficient );


    // sq:

    // VvPlusJjAaAllSq
    current_term_coefficient
    = ( VvPlusJjAaAllSq
        * ( -16.0 * ElMinsq * mQsq * mQsq * mQCsqDiff ) );

    // MinusTwiceVvJjAa
    current_term_coefficient
    += ( MinusTwiceVvJjAa
         * ( 32.0 * ElMinsq * mQsq * mQsq * ( 2.0 * mCsq - mQsq ) ) );

    // VvSqMinusAaSq
    // 0.0

    // JjVvPlusJjAaAllSq
    // 0.0

    // JjVvSqPlusAaSq
    // 0.0

    // JjVvSqMinusAaSq
    // 0.0

    MINtoLK_sq->setCoefficient( current_term_coefficient );


    // LK to HK:
    LKtoHK_segment.setSegmentRange( Elk,
                                      Ehk );

    if( energy_limited_cosines )
      // if the cosine limits are limited by kinematics...
      {

        // const:

        // VvPlusJjAaAllSq
        current_term_coefficient
        = ( VvPlusJjAaAllSq
            * ( mQsq * mQCsqDiff
                * ( EightmCEVElMin * mVsq
                    * ( ElMin * ( lnElMin - lnVsqOverFourElMinsq ) + mC )
                    - 4.0 * EightmCEVElMin * ElMinsq * mC
                    + mCsq * ( 16.0 * ElMinsq * ElMinsq - mVqu ) ) ) );

        // MinusTwiceVvJjAa
        current_term_coefficient
        += ( MinusTwiceVvJjAa
             * ( 2.0 * mQsq
                 * ( ( lnElMin - lnVsqOverFourElMinsq ) * mVsq * ElMin
                     * ( EightmCEVElMin * mQsq
                         + 4.0 * ElMin * mCsq * mQCsqSum )
                     + EightmCEVElMin * mQsq * mC
                       * ( mVsq - 4.0 * ElMinsq )
                     + 8.0 * ElMinsq * ElMin * mC
                       * ( 2.0 * ElMin * mC * mQsq - mQCsqSum * mVsq )
                     + ( 2.0 * ElMin * mQCsqSum - mC * mQsq ) * mVqu ) ) );

        // VvSqMinusAaSq
        current_term_coefficient
        += ( VvSqMinusAaSq
             * ( EightmCXElMinmQsq * ElMin * mQCsqDiff * mVsq
                 * ( lnVsqOverFourElMinsq - lnElMin ) ) );

        // JjVvPlusJjAaAllSq
        current_term_coefficient
        += ( JjVvPlusJjAaAllSq
             * ( ElMin * mQsq * mQCsqDiff * mVsq
                 * ( EightmCEVElMin * ( lnElMin - lnVsqOverFourElMinsq )
                     + 4.0 * mC * ( mVsq - 4.0 * ElMinsq ) ) ) );

        // JjVvSqPlusAaSq
        current_term_coefficient
        += ( JjVvSqPlusAaSq
             * ( 2.0 * ElMin * mQsq * mVsq
                 * ( EightmCEVElMin * mQsq
                     * ( lnVsqOverFourElMinsq - lnElMin )
                     + 8.0 * ElMinsq * mC * ( mCsq + 3.0 * mQsq )
                     + 4.0 * ElMin * mCsq * mQCsqSum
                       * ( lnVsqOverFourElMinsq - lnElMin )
                     - 2.0 * mC * ( mCsq + 3.0 * mQsq ) * mVsq ) ) );

        // JjVvSqMinusAaSq
        current_term_coefficient
        += ( JjVvSqMinusAaSq
             * ( EightmCXElMinmQsq * ElMin * mQCsqSum * mVsq
                 * ( lnVsqOverFourElMinsq - lnElMin ) ) );

        LKtoHK_const->setCoefficient( current_term_coefficient );


        // lin:

        // VvPlusJjAaAllSq
        // 0.0

        // MinusTwiceVvJjAa
        current_term_coefficient
        = ( MinusTwiceVvJjAa
            * ( 2.0 * mC * mQcu
                * ( EightmCEVElMin
                    * ( mVsq
                        - 4.0 * ElMin
                          * ( ElMin
                              + ( lnVsqOverFourElMinsq - lnElMin ) * mC ) )
                    - 2.0 * mC
                      * ( 8.0 * ElMin * ( 2.0 * ElMinsq * ( 2.0 * mC - ElMin )
                                          - mC * mVsq )
                          + mVqu ) ) ) );

        // VvSqMinusAaSq
        // 0.0

        // JjVvPlusJjAaAllSq
        // 0.0

        // JjVvSqPlusAaSq
        current_term_coefficient
        += ( JjVvSqPlusAaSq
             * ( 2.0 * mC * mQcu
                 * ( 4.0 * ElMinsq
                     * ( EightmCEVElMin
                         + 8.0 * ElMin * mC * ( mC - ElMin ) )
                     - mVsq * ( EightmCEVElMin
                                + 8.0 * ElMin * mC
                                  * ( 2.0 * ElMin
                                      * ( lnElMin - lnVsqOverFourElMinsq )
                                      + mC ) )
                   - 2.0 * mC * mVqu ) ) );

        // JjVvSqMinusAaSq
        current_term_coefficient
        += ( JjVvSqMinusAaSq
             * ( 2.0 * EightmCXElMinmQsq * mC * mQ
                 * (4.0 * ElMinsq - mVsq ) ) );

        LKtoHK_lin->setCoefficient( current_term_coefficient );


        // sq:

        // VvPlusJjAaAllSq
        // 0.0

        // MinusTwiceVvJjAa
        // 0.0

        // VvSqMinusAaSq
        // 0.0

        // JjVvPlusJjAaAllSq
        // 0.0

        // JjVvSqPlusAaSq
        // 0.0

        // JjVvSqMinusAaSq
        // 0.0

        LKtoHK_sq->setCoefficient( 0.0 );

      }
    else
      // otherwise, if the cosine limits are not limited by kinematics...
      {

        // const:

        // VvPlusJjAaAllSq
        current_term_coefficient
        = ( VvPlusJjAaAllSq
            * ( -2.0 * EightmCEVElMin * ElMin * lnmQC * mQsq * mQCsqDiff
                * mVsq ) );

        // MinusTwiceVvJjAa
        current_term_coefficient
        += ( MinusTwiceVvJjAa
             * ( 2.0 * ElMin * mQsq * mVsq
                 * ( mQCsqDiff * ( EightmCEVElMin + 8.0 * mCsq * ElMin )
                     - 2.0 * lnmQC * ( EightmCEVElMin * mQsq
                                       + 4.0 * ElMin * mCsq * mQCsqSum ) ) ) );

        // VvSqMinusAaSq
        current_term_coefficient
        += ( VvSqMinusAaSq
             * ( 2.0 * EightmCXElMinmQsq * ElMin * lnmQC * mQCsqDiff
                 * mVsq ) );

        // JjVvPlusJjAaAllSq
        current_term_coefficient
        += ( JjVvPlusJjAaAllSq
             * ( -2.0 * EightmCEVElMin * ElMin * lnmQC * mQsq * mQCsqDiff
                 * mVsq ) );

        // JjVvSqPlusAaSq
        current_term_coefficient
        += ( JjVvSqPlusAaSq
             * ( 2.0 * ElMin * mQsq * mVsq
                 * ( 2.0 * lnmQC * ( EightmCEVElMin * mQsq
                                     + 4.0 * ElMin * mCsq * mQCsqSum )
                     - mQCsqDiff
                       * ( EightmCEVElMin + 8.0 * ElMin * mCsq ) ) ) );

        // JjVvSqMinusAaSq
        current_term_coefficient
        += ( JjVvSqMinusAaSq
             * ( 2.0 * EightmCXElMinmQsq * ElMin * mVsq
                 * ( lnmQC * mQCsqSum - mQCsqDiff ) ) );

        LKtoHK_const->setCoefficient( current_term_coefficient );


        // lin:

        // VvPlusJjAaAllSq
        current_term_coefficient
        = ( VvPlusJjAaAllSq
            * ( 4.0 * EightmCEVElMin * ElMin * mQ * mQCsqDiffsq ) );

        // MinusTwiceVvJjAa
        current_term_coefficient
        += ( MinusTwiceVvJjAa
             * ( 8.0 * EightmCEVElMin * ElMin * mQcu
                 * ( mQCsqDiff - 2.0 * mCsq * lnmQC ) ) );

        // VvSqMinusAaSq
        // 0.0

        // JjVvPlusJjAaAllSq
        current_term_coefficient
        += ( JjVvPlusJjAaAllSq
             * ( 16.0 * ElMinsq * mQ * mQCsqDiffsq * mVsq ) );

        // JjVvSqPlusAaSq
        current_term_coefficient
        += ( JjVvSqPlusAaSq
             * ( 32.0 * ElMinsq * mQcu * mVsq
                 * ( 2.0 * mCsq * lnmQC - mQCsqDiff ) ) );

        // JjVvSqMinusAaSq
        // 0.0

        LKtoHK_lin->setCoefficient( current_term_coefficient );


        // sq:

        // VvPlusJjAaAllSq
        current_term_coefficient
        = ( VvPlusJjAaAllSq
            * ( -16.0 * ElMinsq * mQCsqDiffsq * mQCsqSum ) );

        // MinusTwiceVvJjAa
        current_term_coefficient
        += ( MinusTwiceVvJjAa
             * ( -32.0 * ElMinsq * mQsq * mQCsqDiffsq ) );

        // VvSqMinusAaSq
        // 0.0

        // JjVvPlusJjAaAllSq
        // 0.0

        // JjVvSqPlusAaSq
        // 0.0

        // JjVvSqMinusAaSq
        // 0.0

        LKtoHK_sq->setCoefficient( current_term_coefficient );

      }


    // HK to MAX:
    HKtoMAX_segment.setSegmentRange( Ehk,
                                       maximumEnergy );

    // inv:

    // VvPlusJjAaAllSq
    // 0.0

    // MinusTwiceVvJjAa
    current_term_coefficient
    = ( MinusTwiceVvJjAa
        * ( 4.0 * ElMinsq * mCsq * mQcu * mVqu ) );

    // VvSqMinusAaSq
    // 0.0

    // JjVvPlusJjAaAllSq
    // 0.0

    // JjVvSqPlusAaSq
    current_term_coefficient
    += ( JjVvSqPlusAaSq
         * ( -4.0 * ElMinsq * mCsq * mQcu * mVqu ));

    // JjVvSqMinusAaSq
    // 0.0

    HKtoMAX_inv->setCoefficient( current_term_coefficient );


    // const:

    // VvPlusJjAaAllSq
    current_term_coefficient
    = ( VvPlusJjAaAllSq
        * ( mQsq * mQCsqDiff * mVsq
            * ( EightmCEVElMin
                * ( mC - ElMin * ( lnmQC + lnVsqOverFourElMinsq ) )
                - mCsq * mVsq ) ) ) ;

    // MinusTwiceVvJjAa
    current_term_coefficient
    += ( MinusTwiceVvJjAa
         * ( 2.0 * mQsq * mVsq
             * ( EightmCEVElMin * mQsq * ( ElMin + mC )
                 - ( lnmQC + lnVsqOverFourElMinsq )
                   * ( mQsq * EightmCEVElMin * ElMin
                       + 4.0 * ElMinsq * mCsq * mQCsqSum )
                 + 8.0 * ElMinsq * mC
                   * ( mC * mQsq - ElMin * mQCsqSum )
                 - mCsq * mQsq * mVsq ) ) );

    // VvSqMinusAaSq
    current_term_coefficient
    += ( VvSqMinusAaSq
         * ( EightmCXElMinmQsq * ElMin * mQCsqDiff * mVsq
             * ( lnmQC + lnVsqOverFourElMinsq ) ) );

    // JjVvPlusJjAaAllSq
    current_term_coefficient
    += ( JjVvPlusJjAaAllSq
         * ( ElMin * mQsq * mQCsqDiff * mVsq
             * ( 4.0 * mC * mVsq
                 - EightmCEVElMin * ( lnmQC + lnVsqOverFourElMinsq ) ) ) );

    // JjVvSqPlusAaSq
    current_term_coefficient
    += ( JjVvSqPlusAaSq
         * ( 2.0 * ElMin * mQsq * mVsq
             * ( ( lnmQC + lnVsqOverFourElMinsq )
                 * ( 4.0 * ElMin * mCsq * mQCsqSum + mQsq * EightmCEVElMin )
                 + 8.0 * ElMin * mC * ( ElMin * mQCsqSum - mC * mQsq )
                 - EightmCEVElMin * mQsq
                 - 4.0 * mQsq * mC * mVsq ) ) );

    // JjVvSqMinusAaSq
    current_term_coefficient
    += ( JjVvSqMinusAaSq
         * ( EightmCXElMinmQsq * ElMin * mVsq
             * ( ( lnmQC + lnVsqOverFourElMinsq ) * mQCsqSum
                 - 2.0 * mQsq ) ) );

    HKtoMAX_const->setCoefficient( current_term_coefficient );


    // log:

    // VvPlusJjAaAllSq
    current_term_coefficient
    = ( VvPlusJjAaAllSq
        * ( EightmCEVElMin * ElMin * mQsq * mQCsqDiff * mVsq ) );

    // MinusTwiceVvJjAa
    current_term_coefficient
    += ( MinusTwiceVvJjAa
         * ( 2.0 * ElMin * mQsq * mVsq
             * ( EightmCEVElMin * mQsq + 4.0 * ElMin * mCsq * mQCsqSum ) ) );

    // VvSqMinusAaSq
    current_term_coefficient
    += ( VvSqMinusAaSq
         * ( -EightmCXElMinmQsq * ElMin * mQCsqDiff * mVsq ) );

    // JjVvPlusJjAaAllSq
    current_term_coefficient
    += ( JjVvPlusJjAaAllSq
         * ( EightmCEVElMin * ElMin * mQsq * mQCsqDiff * mVsq ) );

    // JjVvSqPlusAaSq
    current_term_coefficient
    += ( JjVvSqPlusAaSq
         * ( -2.0 * ElMin * mQsq * mVsq
             * ( EightmCEVElMin * mQsq
                 + 4.0 * ElMin * mCsq * mQCsqSum ) ) );

    // JjVvSqMinusAaSq
    current_term_coefficient
    += ( JjVvSqMinusAaSq
         * ( -EightmCXElMinmQsq * ElMin * mQCsqSum * mVsq ) );

    HKtoMAX_log->setCoefficient( current_term_coefficient );


    // lin:

    // VvPlusJjAaAllSq
    current_term_coefficient
    = ( VvPlusJjAaAllSq
        * ( -4.0 * EightmCEVElMin * ElMin * mCsq * mQ * mQCsqDiff ) );

    // MinusTwiceVvJjAa
    current_term_coefficient
    += ( MinusTwiceVvJjAa
         * ( 8.0 * ElMin * mC * mQcu
             * ( 2.0 * mC * ( 4.0 * ElMinsq * ( ElMin - mC ) + mC * mVsq )
                 - EightmCEVElMin
                   * ( ElMin
                       + ( 1.0 + lnmQC + lnVsqOverFourElMinsq ) * mC ) ) ) );

    // VvSqMinusAaSq
    // 0.0

    // JjVvPlusJjAaAllSq
    current_term_coefficient
    += ( JjVvPlusJjAaAllSq
         * ( -16.0 * ElMinsq * mCsq * mQ * mQCsqDiff * mVsq ) );

    // JjVvSqPlusAaSq
    current_term_coefficient
    += ( JjVvSqPlusAaSq
         * ( 8.0 * ElMinsq * mC * mQcu
             * ( EightmCEVElMin + 8.0 * ElMin * mC * ( mC - ElMin )
                 + 4.0 * mC * mVsq
                   * ( 1.0 + lnmQC + lnVsqOverFourElMinsq ) ) ) );

    // JjVvSqMinusAaSq
    current_term_coefficient
    += ( JjVvSqMinusAaSq
         * ( 8.0 * EightmCXElMinmQsq * ElMinsq * mC * mQ ) );

    HKtoMAX_lin->setCoefficient( current_term_coefficient );


    // linlog:

    // VvPlusJjAaAllSq
    // 0.0

    // MinusTwiceVvJjAa
    current_term_coefficient
    = ( MinusTwiceVvJjAa
        * ( 8.0 * EightmCEVElMin * ElMin * mCsq * mQcu ) );

    // VvSqMinusAaSq
    // 0.0

    // JjVvPlusJjAaAllSq
    // 0.0

    // JjVvSqPlusAaSq
    current_term_coefficient
    += ( JjVvSqPlusAaSq
         * ( -32.0 * ElMinsq * mCsq * mQcu * mVsq ) );

    // JjVvSqMinusAaSq
    // 0.0

    HKtoMAX_linlog->setCoefficient( current_term_coefficient );


    // sq:

    // VvPlusJjAaAllSq
    current_term_coefficient
    = ( VvPlusJjAaAllSq
        * ( 16.0 * ElMinsq * mCsq * mCsq * mQCsqDiff ) );

    // MinusTwiceVvJjAa
    current_term_coefficient
    += ( MinusTwiceVvJjAa
         * ( -32.0 * ElMinsq * mCsq * mCsq * mQsq ) );

    // VvSqMinusAaSq
    // 0.0

    // JjVvPlusJjAaAllSq
    // 0.0

    // JjVvSqPlusAaSq
    // 0.0

    // JjVvSqMinusAaSq
    // 0.0

    HKtoMAX_sq->setCoefficient( current_term_coefficient );

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: lepton_distribution_same_chirality_far_muon:";
    std::cout << std::endl << "normalization = " << normalization;
    std::cout << std::endl << "minimum_energy = " << minimum_energy;
    std::cout << std::endl << "Elk = " << Elk;
    std::cout << std::endl << "Ehk = " << Ehk;
    std::cout << std::endl << "maximum_energy = " << maximum_energy;
    std::cout << std::endl << "MINtoLK_lin = " << MINtoLK_lin;
    std::cout << std::endl << "MINtoLK_log = " << MINtoLK_log;
    std::cout << std::endl << "MINtoLK_const = " << MINtoLK_const;
    std::cout << std::endl << "MINtoLK_inv = " << MINtoLK_inv;
    std::cout << std::endl << "LKtoHK_const = " << LKtoHK_const;
    std::cout << std::endl << "LKtoHK_inv = " << LKtoHK_inv;
    std::cout << std::endl;**/
  }



  vector_from_squark_to_muon::vector_from_squark_to_muon(
                                       readierForNewPoint* const given_readier,
                                 CppSLHA::CppSLHA0 const* const given_spectrum,
              CppSLHA::particle_property_set const* const given_first_particle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const given_second_particle,
           CppSLHA::particle_property_set const* const given_third_particle ) :
    leptonEnergyDistribution( given_readier,
                              given_spectrum,
                              given_first_particle,
                              effectiveSquarkMass,
                              given_second_particle,
                              given_third_particle,
                              NULL ),
    MINtoMAX_const( MINtoMAX_segment.addTerm( 0,
                                              0 ) ),
    MINtoMAX_lin( MINtoMAX_segment.addTerm( 1,
                                            0 ) ),
    MINtoMAX_sq( MINtoMAX_segment.addTerm( 2,
                                           0 ) )
  {
    segments.push_back( &MINtoMAX_segment );
  }

  vector_from_squark_to_muon::~vector_from_squark_to_muon()
  {

    // does nothing.

  }

  void
  vector_from_squark_to_muon::calculateCoefficients()
  {
    mQhSq = ( firstMass * firstMass );
    mQlSq = ( secondMass * secondMass );
    mVSq = ( thirdMass * thirdMass );
    gammaV = ( ( mQhSq - mQlSq + mVSq ) / ( 2.0 * firstMass * thirdMass ) );
    betaV = sqrt( ( 1.0 - ( 1.0 / ( gammaV * gammaV ) ) ) );
    productionFrameEnergy = ( 0.5 * thirdMass );
    minimumEnergy = ( gammaV * ( 1.0 - betaV ) * productionFrameEnergy );
    maximumEnergy = ( gammaV * ( 1.0 + betaV ) * productionFrameEnergy );
    energyDifference = ( maximumEnergy - minimumEnergy );

    // now we set up the coefficients of various terms:
    // MIN to MAX segment:
    MINtoMAX_segment.setSegmentRange( minimumEnergy,
                                      maximumEnergy );
    MINtoMAX_const->setCoefficient( -minimumEnergy * maximumEnergy );
    MINtoMAX_lin->setCoefficient( minimumEnergy + maximumEnergy );
    MINtoMAX_sq->setCoefficient( -1.0 );

    //normalization
    //= ( ( energyDifference * energyDifference * energyDifference ) / 6.0 );
  }



  Z_direct_jet::Z_direct_jet( readierForNewPoint* const given_readier,
                              CppSLHA::CppSLHA0 const* const given_spectrum,
              CppSLHA::particle_property_set const* const given_first_particle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const given_second_particle,
              CppSLHA::particle_property_set const* const given_third_particle,
          CppSLHA::particle_property_set const* const given_fourth_particle ) :
    zHandedMuon( given_readier,
                   given_spectrum,
                   given_first_particle,
                   effectiveSquarkMass,
                   given_second_particle,
                   given_third_particle,
                   given_fourth_particle,
                   true,
                   true )
  {

    // just an initialization of the base class object.

  }

  Z_direct_jet::~Z_direct_jet()
  {

    // does nothing.

  }



  W_minus_direct_jet::W_minus_direct_jet(
                                    readierForNewPoint* const given_readier,
                                 CppSLHA::CppSLHA0 const* const given_spectrum,
              CppSLHA::particle_property_set const* const given_first_particle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const given_second_particle,
              CppSLHA::particle_property_set const* const given_third_particle,
          CppSLHA::particle_property_set const* const given_fourth_particle ) :
    zHandedMuon( given_readier,
                   given_spectrum,
                   given_first_particle,
                   effectiveSquarkMass,
                   given_second_particle,
                   given_third_particle,
                   given_fourth_particle,
                   true,
                   true )
  {

    // just an initialization of the base class object.

  }

  W_minus_direct_jet::~W_minus_direct_jet()
  {

    // does nothing.

  }



  negatively_charged_Higgs_muon::negatively_charged_Higgs_muon(
                                    readierForNewPoint* const given_readier,
                                 CppSLHA::CppSLHA0 const* const given_spectrum,
              CppSLHA::particle_property_set const* const given_first_particle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const given_second_particle,
              CppSLHA::particle_property_set const* const given_third_particle,
          CppSLHA::particle_property_set const* const given_fourth_particle ) :
    Higgs_muon_plus_antimuon( given_readier,
                              given_spectrum,
                              given_first_particle,
                              effectiveSquarkMass,
                              given_second_particle,
                              given_third_particle,
                              given_fourth_particle )
  {

    // just an initialization of the base class object.

  }

  negatively_charged_Higgs_muon::~negatively_charged_Higgs_muon()
  {

    // does nothing.

  }



  negatively_charged_Higgs_jet::negatively_charged_Higgs_jet(
                                    readierForNewPoint* const given_readier,
                                 CppSLHA::CppSLHA0 const* const given_spectrum,
              CppSLHA::particle_property_set const* const given_first_particle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const given_second_particle,
              CppSLHA::particle_property_set const* const given_third_particle,
          CppSLHA::particle_property_set const* const given_fourth_particle ) :
    Higgs_muon_plus_antimuon( given_readier,
                              given_spectrum,
                              given_first_particle,
                              effectiveSquarkMass,
                              given_second_particle,
                              given_third_particle,
                              given_fourth_particle )
  {

    // just an initialization of the base class object.

  }

  negatively_charged_Higgs_jet::~negatively_charged_Higgs_jet()
  {

    // does nothing.

  }



  neutralino_three_body_decay::neutralino_three_body_decay(
                                    readierForNewPoint* const given_readier,
                                 CppSLHA::CppSLHA0 const* const given_spectrum,
              CppSLHA::particle_property_set const* const given_first_particle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const given_second_particle,
              CppSLHA::particle_property_set const* const given_third_particle,
               CppSLHA::particle_property_set const* const given_left_sfermion,
           CppSLHA::particle_property_set const* const given_right_sfermion ) :
    leptonEnergyDistribution( given_readier,
                                given_spectrum,
                                given_first_particle,
                                effectiveSquarkMass,
                                given_second_particle,
                                given_third_particle,
                                NULL ),
    left_sfermion( given_left_sfermion ),
    right_sfermion( given_right_sfermion ),
    MINtoMAX_const( MINtoMAX_segment.addTerm( 0,
                                              0,
                                              NULL,
                                              1.0 ) )
  {
    segments.push_back( &MINtoMAX_segment );
  }

  neutralino_three_body_decay::~neutralino_three_body_decay()
  {

    // does nothing.

  }


  void
  neutralino_three_body_decay::calculateCoefficients()
  {

    minimumEnergy = 0.0;

    maximumEnergy = ( secondMass - thirdMass );


    // now we set up the coefficients of various terms.

    // MIN to MAX segment:
    MINtoMAX_segment.setSegmentRange( minimumEnergy,
                                        maximumEnergy );
  }



  chargino_three_body_decay::chargino_three_body_decay(
                                    readierForNewPoint* const given_readier,
                                 CppSLHA::CppSLHA0 const* const given_spectrum,
              CppSLHA::particle_property_set const* const given_first_particle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const given_second_particle,
              CppSLHA::particle_property_set const* const given_third_particle,
             CppSLHA::particle_property_set const* const leftUpIsospinSfermion,
            CppSLHA::particle_property_set const* const rightUpIsospinSfermion,
           CppSLHA::particle_property_set const* const leftDownIsospinSfermion,
       CppSLHA::particle_property_set const* const rightDownIsospinSfermion ) :
    leptonEnergyDistribution( given_readier,
                              given_spectrum,
                              given_first_particle,
                              effectiveSquarkMass,
                              given_second_particle,
                              given_third_particle,
                              NULL ),
    leftUpIsospinSfermion( leftUpIsospinSfermion ),
    rightUpIsospinSfermion( rightUpIsospinSfermion ),
    leftDownIsospinSfermion( leftDownIsospinSfermion ),
    rightDownIsospinSfermion( rightDownIsospinSfermion ),
    MINtoMAX_const( MINtoMAX_segment.addTerm( 0,
                                              0,
                                              NULL,
                                              1.0 ) )
  {
    segments.push_back( &MINtoMAX_segment );
  }

  chargino_three_body_decay::~chargino_three_body_decay()
  {

    // does nothing.

  }


  void
  chargino_three_body_decay::calculateCoefficients()
  {

    minimumEnergy = 0.0;

    maximumEnergy = ( secondMass - thirdMass );


    // now we set up the coefficients of various terms, & use them to determine
    // the normalization:

    // MIN to MAX segment:
    MINtoMAX_segment.setSegmentRange( minimumEnergy,
                                        maximumEnergy );
  }

}  // end of LHC_FASER namespace.

