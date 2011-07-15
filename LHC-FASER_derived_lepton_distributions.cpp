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
                                    readier_for_new_point* const given_readier,
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
    MINtoMAX_const( 0,
                    0,
                    CppSLHA::CppSLHA_global::really_wrong_value )
  {

    MINtoMAX_segment.add_term( &MINtoMAX_const );

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

    minimumEnergy
    = ( ( secondMass / firstMass ) * productionFrameEnergy );

    maximumEnergy
    = ( ( firstMass / secondMass ) * productionFrameEnergy );


    // now we set up the coefficients of various terms, & use them to determine
    // the normalization:

    // well, in this case, the distribution is a constant in the only segment's
    // range, so I won't bother with writing it.

    normalization = ( maximumEnergy - minimumEnergy );


    // debugging:
    /**std::cout
    << std::endl
    << "debugging: flat_near_muon_plus_antimuon constructor has"
    << std::endl << "first_mass = " << first_mass
    << std::endl << "second_mass = " << second_mass
    << std::endl << "production_frame_energy = " << production_frame_energy
    << std::endl << "minimum_energy = " << minimum_energy
    << std::endl << "maximum_energy = " << maximum_energy
    << std::endl << "normalization = " << normalization;
    std::cout << std::endl;**/


    // finally we update the sets of terms, segment by segment:

    // MIN to MAX segment:
    MINtoMAX_segment.set_segment_range( minimumEnergy,
                                        maximumEnergy );

    // constant term:
    MINtoMAX_const.set_coefficient( ( 1.0 / normalization ) );

  }



  same_chirality_near_muon::same_chirality_near_muon(
                                    readier_for_new_point* const given_readier,
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
    MINtoMAX_const( 0,
                    0,
                    CppSLHA::CppSLHA_global::really_wrong_value ),
    MINtoMAX_lin( 1,
                  0,
                  CppSLHA::CppSLHA_global::really_wrong_value )
  {

    MINtoMAX_segment.add_term( &MINtoMAX_const );
    MINtoMAX_segment.add_term( &MINtoMAX_lin );

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
    MINtoMAX_segment.set_segment_range( minimumEnergy,
                                        maximumEnergy );

    MINtoMAX_const.set_coefficient( -minimumEnergy );

    MINtoMAX_lin.set_coefficient( 1.0 );


    normalization = ( 0.5 * ( maximumEnergy - minimumEnergy )
                          * ( maximumEnergy - minimumEnergy ) );
    normalizeCoefficients();

  }



  opposite_chirality_near_muon::opposite_chirality_near_muon(
                                    readier_for_new_point* const given_readier,
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
    MINtoMAX_const( 0,
                    0,
                    CppSLHA::CppSLHA_global::really_wrong_value ),
    MINtoMAX_lin( 1,
                  0,
                  CppSLHA::CppSLHA_global::really_wrong_value )
  {

    MINtoMAX_segment.add_term( &MINtoMAX_const );
    MINtoMAX_segment.add_term( &MINtoMAX_lin );

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

    normalization = ( maximumEnergy - minimumEnergy );


    // now we set up the coefficients of various terms, & use them to determine
    // the normalization:

    // MIN to MAX segment:
    MINtoMAX_segment.set_segment_range( minimumEnergy,
                                        maximumEnergy );

    MINtoMAX_const.set_coefficient( maximumEnergy );

    MINtoMAX_lin.set_coefficient( -1.0 );


    normalization = ( 0.5 * ( maximumEnergy - minimumEnergy )
                          * ( maximumEnergy - minimumEnergy ) );
    normalizeCoefficients();

  }



  flat_far_muon_plus_antimuon::flat_far_muon_plus_antimuon(
                                    readier_for_new_point* const given_readier,
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
    MINtoLK_const( 0,
                   0,
                   CppSLHA::CppSLHA_global::really_wrong_value ),
    MINtoLK_log( 0,
                 1,
                 CppSLHA::CppSLHA_global::really_wrong_value ),
    LKtoHK_const( 0,
                  0,
                  CppSLHA::CppSLHA_global::really_wrong_value ),
    HKtoMAX_const( 0,
                   0,
                   CppSLHA::CppSLHA_global::really_wrong_value ),
    HKtoMAX_log( 0,
                 1,
                 CppSLHA::CppSLHA_global::really_wrong_value )
  {

    MINtoLK_segment.add_term( &MINtoLK_const );
    MINtoLK_segment.add_term( &MINtoLK_log );

    segments.push_back( &MINtoLK_segment );


    LKtoHK_segment.add_term( &LKtoHK_const );

    segments.push_back( &LKtoHK_segment );


    HKtoMAX_segment.add_term( &HKtoMAX_const );
    HKtoMAX_segment.add_term( &HKtoMAX_log );

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
    if( ( firstMass / secondMass ) > ( secondMass / thirdMass ) )
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
    MINtoLK_segment.set_segment_range( minimumEnergy,
                                       Elk );

    MINtoLK_const.set_coefficient( -log( minimumEnergy ) );

    MINtoLK_log.set_coefficient( 1.0 );


    // LK to HK segment:
    LKtoHK_segment.set_segment_range( Elk,
                                      Ehk );

    LKtoHK_const.set_coefficient( log( ( maximumEnergy / Ehk ) ) );


    // HK to MAX segment:
    HKtoMAX_segment.set_segment_range( Elk,
                                       Ehk );

    HKtoMAX_const.set_coefficient( log( maximumEnergy ) );

    HKtoMAX_log.set_coefficient( -1.0 );


    normalization
    = ( MINtoLK_const.get_coefficient() * ( Elk - minimumEnergy )
        + MINtoLK_log.get_coefficient()
          * ( Elk * ( log( Elk ) - 1.0 )
              - minimumEnergy * ( log( minimumEnergy ) - 1.0 ) )
        + LKtoHK_const.get_coefficient() * ( Ehk - Elk )
        + HKtoMAX_const.get_coefficient() * ( maximumEnergy - Ehk )
        + HKtoMAX_log.get_coefficient()
          * ( maximumEnergy * ( log( maximumEnergy ) - 1.0 )
              - Ehk * ( log( Ehk ) - 1.0 ) ) );
    normalizeCoefficients();

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: flat_far_muon_plus_antimuon_distribution:";
    std::cout << std::endl << "normalization = " << normalization;
    std::cout << std::endl << "minimum_energy = " << minimum_energy;
    std::cout << std::endl << "Elk = " << Elk;
    std::cout << std::endl << "Ehk = " << Ehk;
    std::cout << std::endl << "maximum_energy = " << maximum_energy;
    std::cout << std::endl;**/

  }



  same_chirality_far_muon::same_chirality_far_muon(
                                    readier_for_new_point* const given_readier,
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
    MINtoLK_inv( -1,
                 0,
                 CppSLHA::CppSLHA_global::really_wrong_value ),
    MINtoLK_const( 0,
                   0,
                   CppSLHA::CppSLHA_global::really_wrong_value ),
    MINtoLK_log( 0,
                 1,
                 CppSLHA::CppSLHA_global::really_wrong_value ),
    MINtoLK_lin( 1,
                 0,
                 CppSLHA::CppSLHA_global::really_wrong_value ),
    LKtoHK_inv( -1,
                0,
                CppSLHA::CppSLHA_global::really_wrong_value ),
    LKtoHK_const( 0,
                  0,
                  CppSLHA::CppSLHA_global::really_wrong_value ),
    HKtoMAX_inv( -1,
                 0,
                 CppSLHA::CppSLHA_global::really_wrong_value ),
    HKtoMAX_const( 0,
                   0,
                   CppSLHA::CppSLHA_global::really_wrong_value ),
    HKtoMAX_log( 0,
                 1,
                 CppSLHA::CppSLHA_global::really_wrong_value ),
    HKtoMAX_lin( 1,
                 0,
                 CppSLHA::CppSLHA_global::really_wrong_value )
  {

    MINtoLK_segment.add_term( &MINtoLK_inv );
    MINtoLK_segment.add_term( &MINtoLK_const );
    MINtoLK_segment.add_term( &MINtoLK_log );
    MINtoLK_segment.add_term( &MINtoLK_lin );

    segments.push_back( &MINtoLK_segment );


    LKtoHK_segment.add_term( &LKtoHK_inv );
    LKtoHK_segment.add_term( &LKtoHK_const );

    segments.push_back( &LKtoHK_segment );


    HKtoMAX_segment.add_term( &HKtoMAX_inv );
    HKtoMAX_segment.add_term( &HKtoMAX_const );
    HKtoMAX_segment.add_term( &HKtoMAX_log );
    HKtoMAX_segment.add_term( &HKtoMAX_lin );

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
    if( ( firstMass / secondMass ) > ( secondMass / thirdMass ) )
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
    MINtoLK_segment.set_segment_range( minimumEnergy,
                                       Elk );

    MINtoLK_lin.set_coefficient(
        ( ( 4.0 * mlSq * mNSq * mq )
          / ( ( mlSq - mNSq ) * ( mNSq + mqSq ) * ( mlSq - mXSq ) ) ) );

    MINtoLK_log.set_coefficient(
        ( ( 2.0 * mNSq * ( mqSq - mlSq ) )
          / ( ( mNSq - mlSq ) * ( mNSq + mqSq ) ) ) );

    MINtoLK_const.set_coefficient(
        ( ( 2.0 * mNSq
            * ( 2.0 * mNSq - mqSq + mlSq
              + ( mqSq - mlSq ) * log( ( ( 2.0 * mq ) / ( mlSq - mXSq ) ) ) ) )
        / ( ( mNSq - mlSq ) * ( mNSq + mqSq ) ) ) );

    MINtoLK_inv.set_coefficient(
        ( ( mNSq * ( 2.0 * mNSq - mqSq ) * ( mlSq - mXSq ) )
          / ( mq * ( mlSq - mNSq ) * ( mNSq + mqSq ) ) ) );


    // LK to HK segment:
    LKtoHK_segment.set_segment_range( Elk,
                                      Ehk );

    // there are 2 kinematic cases for determining the kink energies in terms
    // of the cascade masses.
    if( ( firstMass / secondMass ) < ( secondMass / thirdMass ) )
      // if the assumption that the mass spectrum was SPS1a-like was
      // incorrect...
      {

        LKtoHK_const.set_coefficient(
            ( ( 2.0 * ( ( mlSq + mNSq ) * ( mNSq - mqSq )
                      + mNSq * ( mlSq - mqSq ) * log( ( mNSq / mqSq ) ) ) )
              / ( ( mNSq - mlSq ) * ( mNSq + mqSq ) ) ) );

        LKtoHK_inv.set_coefficient(
            ( ( 2.0 * mNSq * ( mNSq - mqSq ) * ( mlSq - mXSq ) )
              / ( mq * ( mlSq - mNSq ) * ( mNSq + mqSq ) ) ) );

      }
    else
      {

        LKtoHK_const.set_coefficient(
            ( ( 2.0 * ( ( mqSq - mNSq ) * ( mlSq - mNSq )
                        + mNSq * ( mlSq - mqSq ) * log( ( mlSq / mNSq ) ) ) )
              / ( ( mNSq - mlSq ) * ( mNSq + mqSq ) ) ) );

        LKtoHK_inv.set_coefficient( 0.0 );

      }


    // HK to MAX segment:
    HKtoMAX_segment.set_segment_range( Ehk,
                                       maximumEnergy );

    HKtoMAX_lin.set_coefficient( -MINtoLK_lin.get_coefficient() );

    HKtoMAX_log.set_coefficient( -MINtoLK_log.get_coefficient() );

    HKtoMAX_const.set_coefficient(
        ( ( 2.0 * mNSq * ( mlSq + mqSq
                           + ( mqSq - mlSq )
                             * log( ( ( 2.0 * mlSq )
                                      / ( ( mlSq - mXSq ) * mq ) ) ) ) )
          / ( ( mlSq - mNSq ) * ( mNSq + mqSq ) ) ) );

    HKtoMAX_inv.set_coefficient(
        ( ( mNSq * mq * ( mlSq - mXSq ) )
          / ( ( mNSq - mlSq ) * ( mNSq + mqSq ) ) ) );


    normalization
    = ( MINtoLK_lin.get_coefficient()
        * 0.5 * ( Elk * Elk - minimumEnergy * minimumEnergy )
        + MINtoLK_log.get_coefficient()
          * ( Elk * ( log( Elk ) - 1.0 )
              - minimumEnergy * ( log( minimumEnergy ) - 1.0 ) )
        + MINtoLK_const.get_coefficient() * ( Elk - minimumEnergy )
        + MINtoLK_inv.get_coefficient() * log( ( Elk / minimumEnergy ) )
        + LKtoHK_const.get_coefficient() * ( Ehk - Elk )
        + LKtoHK_inv.get_coefficient() * log( ( Ehk / Elk ) )
        + HKtoMAX_lin.get_coefficient()
          * 0.5 * ( maximumEnergy * maximumEnergy - Ehk * Ehk )
        + HKtoMAX_log.get_coefficient()
          * ( maximumEnergy * ( log( maximumEnergy ) - 1.0 )
              - Ehk * ( log( Ehk ) - 1.0 ) )
        + HKtoMAX_const.get_coefficient() * ( maximumEnergy - Ehk )
        + HKtoMAX_inv.get_coefficient() * log( ( maximumEnergy / Ehk ) ) );
    normalizeCoefficients();


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
                                    readier_for_new_point* const given_readier,
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
    MINtoLK_inv( -1,
                 0,
                 CppSLHA::CppSLHA_global::really_wrong_value ),
    MINtoLK_const( 0,
                   0,
                   CppSLHA::CppSLHA_global::really_wrong_value ),
    MINtoLK_log( 0,
                 1,
                 CppSLHA::CppSLHA_global::really_wrong_value ),
    MINtoLK_lin( 1,
                 0,
                 CppSLHA::CppSLHA_global::really_wrong_value ),
    LKtoHK_inv( -1,
                0,
                CppSLHA::CppSLHA_global::really_wrong_value ),
    LKtoHK_const( 0,
                  0,
                  CppSLHA::CppSLHA_global::really_wrong_value ),
    HKtoMAX_inv( -1,
                 0,
                 CppSLHA::CppSLHA_global::really_wrong_value ),
    HKtoMAX_const( 0,
                   0,
                   CppSLHA::CppSLHA_global::really_wrong_value ),
    HKtoMAX_log( 0,
                 1,
                 CppSLHA::CppSLHA_global::really_wrong_value ),
    HKtoMAX_lin( 1,
                 0,
                 CppSLHA::CppSLHA_global::really_wrong_value )
  {

    MINtoLK_segment.add_term( &MINtoLK_inv );
    MINtoLK_segment.add_term( &MINtoLK_const );
    MINtoLK_segment.add_term( &MINtoLK_log );
    MINtoLK_segment.add_term( &MINtoLK_lin );

    segments.push_back( &MINtoLK_segment );


    LKtoHK_segment.add_term( &LKtoHK_inv );
    LKtoHK_segment.add_term( &LKtoHK_const );

    segments.push_back( &LKtoHK_segment );


    HKtoMAX_segment.add_term( &HKtoMAX_inv );
    HKtoMAX_segment.add_term( &HKtoMAX_const );
    HKtoMAX_segment.add_term( &HKtoMAX_log );
    HKtoMAX_segment.add_term( &HKtoMAX_lin );

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
    if( ( firstMass / secondMass ) > ( secondMass / thirdMass ) )
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
    MINtoLK_segment.set_segment_range( minimumEnergy,
                                       Elk );

    MINtoLK_lin.set_coefficient(
        ( ( 4.0 * mlSq * mNSq * mq )
          / ( ( mNSq - mlSq ) * ( mNSq + mqSq ) * ( mlSq - mXSq ) ) ) );

    MINtoLK_log.set_coefficient(
        ( ( 2.0 * ( mNSq * mNSq - mlSq * mqSq ) )
          / ( ( mNSq - mlSq ) * ( mNSq + mqSq ) ) ) );

    MINtoLK_const.set_coefficient(
        ( ( 2.0 * ( mNSq * ( 2.0 * mNSq - mqSq + mlSq )
                  + ( mlSq * mqSq - mNSq * mNSq )
                    * log( ( ( 2.0 * mq ) / ( mlSq - mXSq ) ) ) ) )
          / ( ( mlSq - mNSq ) * ( mNSq + mqSq ) ) ) );

    MINtoLK_inv.set_coefficient(
        ( ( mNSq * ( ( 2.0 * mNSq - mqSq ) * ( mlSq - mXSq ) ) )
          / ( mq * ( mNSq - mlSq ) * ( mNSq + mqSq ) ) ) );


    // LK to HK segment:
    LKtoHK_segment.set_segment_range( Elk,
                                      Ehk );

    // there are 2 kinematic cases for determining the kink energies in terms of the cascade masses.
    if( ( firstMass / secondMass ) < ( secondMass / thirdMass ) )
      // if the assumption that the mass spectrum was SPS1a-like was
      // incorrect...
      {

        LKtoHK_const.set_coefficient(
            ( ( 2.0 * ( ( mlSq + mNSq ) * ( mNSq - mqSq )
                        + ( mNSq * mNSq - mlSq * mqSq )
                          * log( ( mNSq / mqSq ) ) ) )
              / ( ( mlSq - mNSq ) * ( mNSq + mqSq ) ) ) );

        LKtoHK_inv.set_coefficient(
            ( ( 2.0 * mNSq * ( mNSq - mqSq ) * ( mlSq - mXSq ) )
              / ( mq * ( mNSq - mlSq ) * ( mNSq + mqSq ) ) ) );

      }
    else
      {

        LKtoHK_const.set_coefficient(
            ( ( 2.0 * ( ( mqSq - mNSq ) * ( mlSq - mNSq )
                        + ( mNSq * mNSq - mlSq * mqSq )
                          * log( ( mlSq / mNSq ) ) ) )
              / ( ( mlSq - mNSq ) * ( mNSq + mqSq ) ) ) );

        LKtoHK_inv.set_coefficient( 0.0 );

      }


    // HK to MAX segment:
    HKtoMAX_segment.set_segment_range( Ehk,
                                       maximumEnergy );

    HKtoMAX_lin.set_coefficient( -MINtoLK_lin.get_coefficient() );

    HKtoMAX_log.set_coefficient( -MINtoLK_log.get_coefficient() );

    HKtoMAX_const.set_coefficient(
        ( ( 2.0 * ( ( mNSq * mNSq - mlSq * mqSq )
                    * log( ( ( 2.0 * mlSq ) / ( ( mlSq - mXSq ) * mq ) ) )
                    - mNSq * ( mqSq + mlSq ) ) )
          / ( ( mlSq - mNSq ) * ( mNSq + mqSq ) ) ) );

    HKtoMAX_inv.set_coefficient(
        ( ( mNSq * mq * ( mlSq - mXSq ) )
          / ( ( mlSq - mNSq ) * ( mNSq + mqSq ) ) ) );


    normalization
    = ( MINtoLK_lin.get_coefficient()
        * 0.5 * ( Elk * Elk - minimumEnergy * minimumEnergy )
        + MINtoLK_log.get_coefficient()
          * ( Elk * ( log( Elk ) - 1.0 )
              - minimumEnergy * ( log( minimumEnergy ) - 1.0 ) )
        + MINtoLK_const.get_coefficient() * ( Elk - minimumEnergy )
        + MINtoLK_inv.get_coefficient() * log( ( Elk / minimumEnergy ) )
        + LKtoHK_const.get_coefficient() * ( Ehk - Elk )
        + LKtoHK_inv.get_coefficient() * log( ( Ehk / Elk ) )
        + HKtoMAX_lin.get_coefficient()
          * 0.5 * ( maximumEnergy * maximumEnergy - Ehk * Ehk )
        + HKtoMAX_log.get_coefficient()
          * ( maximumEnergy * ( log( maximumEnergy ) - 1.0 )
              - Ehk * ( log( Ehk ) - 1.0 ) )
        + HKtoMAX_const.get_coefficient() * ( maximumEnergy - Ehk )
        + HKtoMAX_inv.get_coefficient() * log( ( maximumEnergy / Ehk ) ) );
    normalizeCoefficients();


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
                                    readier_for_new_point* const given_readier,
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
    MINtoLK_const( 0,
                   0,
                   CppSLHA::CppSLHA_global::really_wrong_value ),
    MINtoLK_log( 0,
                 1,
                 CppSLHA::CppSLHA_global::really_wrong_value ),
    LKtoHK_const( 0,
                  0,
                  CppSLHA::CppSLHA_global::really_wrong_value ),
    HKtoMAX_const( 0,
                   0,
                   CppSLHA::CppSLHA_global::really_wrong_value ),
    HKtoMAX_log( 0,
                 1,
                 CppSLHA::CppSLHA_global::really_wrong_value )
  {

    MINtoLK_segment.add_term( &MINtoLK_const );
    MINtoLK_segment.add_term( &MINtoLK_log );

    segments.push_back( &MINtoLK_segment );


    LKtoHK_segment.add_term( &LKtoHK_const );

    segments.push_back( &LKtoHK_segment );


    HKtoMAX_segment.add_term( &HKtoMAX_const );
    HKtoMAX_segment.add_term( &HKtoMAX_log );

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
    MINtoLK_segment.set_segment_range( minimumEnergy,
                                       Elk );

    MINtoLK_const.set_coefficient( -log( minimumEnergy ) );
    MINtoLK_log.set_coefficient( 1.0 );



    // LK to HK segment:
    LKtoHK_segment.set_segment_range( Elk,
                                      Ehk );
    LKtoHK_const.set_coefficient( log( ( maximumEnergy / Ehk ) ) );


    // HK to MAX segment:
    HKtoMAX_segment.set_segment_range( Ehk,
                                       maximumEnergy );

    HKtoMAX_const.set_coefficient( log( maximumEnergy ) );
    HKtoMAX_log.set_coefficient( -1.0 );

    normalization
    = ( MINtoLK_const.get_coefficient() * ( Elk - minimumEnergy )
        + MINtoLK_log.get_coefficient()
          * ( Elk * ( log( Elk ) - 1.0 )
              - minimumEnergy * ( log( minimumEnergy ) - 1.0 ) )
        + LKtoHK_const.get_coefficient() * ( Ehk - Elk )
        + HKtoMAX_const.get_coefficient() * ( maximumEnergy - Ehk )
        + HKtoMAX_log.get_coefficient()
          * ( maximumEnergy * ( log( maximumEnergy ) - 1.0 )
              - Ehk * ( log( Ehk ) - 1.0 ) ) );
    normalizeCoefficients();

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



  Z_handed_muon::Z_handed_muon( readier_for_new_point* const given_readier,
                                CppSLHA::CppSLHA0 const* const given_spectrum,
              CppSLHA::particle_property_set const* const given_first_particle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const given_second_particle,
              CppSLHA::particle_property_set const* const given_third_particle,
          CppSLHA::particle_property_set const* const given_fourth_particle,
                                bool const negative_muon_is_same_handedness,
                                bool const should_sum_over_handedness ) :
    leptonEnergyDistribution( given_readier,
                                given_spectrum,
                                given_first_particle,
                                effectiveSquarkMass,
                                given_second_particle,
                                given_third_particle,
                                given_fourth_particle ),
    MINtoLK_inv( -1,
                 0,
                 CppSLHA::CppSLHA_global::really_wrong_value ),
    MINtoLK_const( 0,
                   0,
                   CppSLHA::CppSLHA_global::really_wrong_value ),
    MINtoLK_log( 0,
                 1,
                 CppSLHA::CppSLHA_global::really_wrong_value ),
    MINtoLK_lin( 1,
                 0,
                 CppSLHA::CppSLHA_global::really_wrong_value ),
    MINtoLK_linlog( 1,
                    1,
                    CppSLHA::CppSLHA_global::really_wrong_value ),
    MINtoLK_sq( 2,
                0,
                CppSLHA::CppSLHA_global::really_wrong_value ),
    LKtoHK_const( 0,
                  0,
                  CppSLHA::CppSLHA_global::really_wrong_value ),
    LKtoHK_lin( 1,
                0,
                CppSLHA::CppSLHA_global::really_wrong_value ),
    LKtoHK_sq( 2,
               0,
               CppSLHA::CppSLHA_global::really_wrong_value ),
    HKtoMAX_inv( -1,
                 0,
                 CppSLHA::CppSLHA_global::really_wrong_value ),
    HKtoMAX_const( 0,
                   0,
                   CppSLHA::CppSLHA_global::really_wrong_value ),
    HKtoMAX_log( 0,
                 1,
                 CppSLHA::CppSLHA_global::really_wrong_value ),
    HKtoMAX_lin( 1,
                 0,
                 CppSLHA::CppSLHA_global::really_wrong_value ),
    HKtoMAX_linlog( 1,
                    1,
                    CppSLHA::CppSLHA_global::really_wrong_value ),
    HKtoMAX_sq( 2,
                0,
                CppSLHA::CppSLHA_global::really_wrong_value )
  {

    if( should_sum_over_handedness )
      {

        same_handedness_factor = 0.0;

      }
    else
      {

        if( negative_muon_is_same_handedness )
          {

            same_handedness_factor = 1.0;

          }
        else
          {

            same_handedness_factor = -1.0;

          }

      }


    MINtoLK_segment.add_term( &MINtoLK_inv );
    MINtoLK_segment.add_term( &MINtoLK_const );
    MINtoLK_segment.add_term( &MINtoLK_log );
    MINtoLK_segment.add_term( &MINtoLK_lin );
    MINtoLK_segment.add_term( &MINtoLK_linlog );
    MINtoLK_segment.add_term( &MINtoLK_sq );

    segments.push_back( &MINtoLK_segment );


    LKtoHK_segment.add_term( &LKtoHK_const );
    LKtoHK_segment.add_term( &LKtoHK_lin );
    LKtoHK_segment.add_term( &LKtoHK_sq );

    segments.push_back( &LKtoHK_segment );


    HKtoMAX_segment.add_term( &HKtoMAX_inv );
    HKtoMAX_segment.add_term( &HKtoMAX_const );
    HKtoMAX_segment.add_term( &HKtoMAX_log );
    HKtoMAX_segment.add_term( &HKtoMAX_lin );
    HKtoMAX_segment.add_term( &HKtoMAX_linlog );
    HKtoMAX_segment.add_term( &HKtoMAX_sq );

    segments.push_back( &HKtoMAX_segment );

  }

  Z_handed_muon::~Z_handed_muon()
  {

    // does nothing.

  }


  void
  Z_handed_muon::calculateCoefficients()
  {

    /* it doesn't matter if the neutralino mass is negative, since only the
     * absolute square appears, & the couplings are only affected by it as an
     * overall phase.  I (BOL) assume no CP violation, so the neutralino mixing
     * matrices are strictly real. just as a quirk of the way I did my
     * workings, axial_coupling_factor needs to be +1 when the Z couples as a
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

        axial_coupling_factor = -1.0;

      }
    else
      {

        axial_coupling_factor = 1.0;

      }

    mQ = firstMass;
    mQsq = ( mQ * mQ );
    mQcu = ( mQ * mQsq );
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


    // MIN to LK:
    MINtoLK_segment.set_segment_range( minimumEnergy,
                                       Elk );
    // inv:

    // symmetric bit positive for both vector & axial coupling:
    // 0.0
    // tick!

    // symmetric bit positive for vector & negative for axial coupling:
    // 0.0
    // tick!

    if( 0.0 != same_handedness_factor )
      {

        /* antisymmetric bit positive for both vector & axial coupling (this is
         * the sum of 2 terms in the Mathematica notebook: VvPlusJjAaAllSq &
         * VvSqPlusAaSq):
         */
        current_term_coefficient
        = ( same_handedness_factor
            * ( 4.0 * ElMinsq * mCsq * mQcu * mVsq ) );
        // tick! tick!

        normalization = ( current_term_coefficient * ( lnElk - lnEmin ) );
        // current_term_coefficient should not be zero at this point.

        // antisymmetric bit positive for vector & negative for axial coupling:
        // 0.0
        // tick!

      }
    else
      {

        current_term_coefficient = 0.0;
        normalization = 0.0;

      }

    MINtoLK_inv.set_coefficient( current_term_coefficient );


    // const:

    // symmetric bit positive for both vector & axial coupling:
    current_term_coefficient
    = ( ElMin * mQsq * mQCsqDiff
        * ( EightmCEVElMin * ( lnElMin - lnmQC ) * mVsq
            - 4.0 * ElMin * mC * ( EightmCEVElMin - 4.0 * ElMinsq * mC ) ) );
    // tick!

    // symmetric bit positive for vector & negative for axial coupling:
    current_term_coefficient
    += ( axial_coupling_factor
         * ( EightmCXElMinmQsq * ElMin * ( lnmQC - lnElMin )
             * mQCsqDiff * mVsq ) );
    // tick!

    if( 0.0 != same_handedness_factor )
      {

        /* antisymmetric bit positive for both vector & axial coupling (this is
         * the sum of 2 terms in the Mathematica notebook: VvPlusJjAaAllSq &
         * VvSqPlusAaSq):
         */
        /*current_term_coefficient
        += ( same_handedness_factor
             * ( ( ElMin * mQsq * mVsq * mQCsqDiff
                   * ( EightmCEVElMin * ( lnElMin - lnmQC )
                       - 16.0 * ElMinsq * mC ) )

                 + ElMin * mQsq * mVsq
                   * ( 32.0 * ElMinsq * mC * mQsq
                       + ( EightmCEVElMin * mQsq
                           + 4.0 * ElMin * mCsq * mQCsqSum )
                         * 2.0 * ( lnmQC - lnElMin )
                       + 2.0 * mCsq * ( EightmCEVElMin + 4.0 * ElMin * mCsq )
                       - 4.0 * mC * mQCsqSum * mVsq ) ) );*/
        current_term_coefficient
        += ( same_handedness_factor
             * ( ( ElMin * mQsq * mVsq
                   * ( mQCsqDiff * ( EightmCEVElMin * ( lnElMin - lnmQC )
                                     - 16.0 * ElMinsq * mC ) )
                       + ( 32.0 * ElMinsq * mC * mQsq
                           + ( EightmCEVElMin * mQsq
                               + 4.0 * ElMin * mCsq * mQCsqSum )
                             * 2.0 * ( lnmQC - lnElMin )
                           + 2.0 * mCsq * ( EightmCEVElMin
                                            + 8.0 * ElMin * mCsq )
                           - 4.0 * mC * mQCsqSum * mVsq ) ) ) );
        // tick! tick!

        // antisymmetric bit positive for vector & negative for axial coupling:
        current_term_coefficient
        += ( same_handedness_factor * axial_coupling_factor
             * ( EightmCXElMinmQsq * ElMin * mVsq
                 * ( 2.0 * mCsq + ( lnmQC - lnElMin ) * mQCsqSum ) ) );
        // tick!

      }

    MINtoLK_const.set_coefficient( current_term_coefficient );
    if( 0.0 != current_term_coefficient )
      {

        normalization
        += ( current_term_coefficient * ( Elk - minimumEnergy ) );

      }


    // log:

    // symmetric bit positive for both vector & axial coupling:
    current_term_coefficient
    = ( -EightmCEVElMin * ElMin * mQsq * mQCsqDiff * mVsq );
    // tick!

    // symmetric bit positive for vector & negative for axial coupling:
    current_term_coefficient
    += ( axial_coupling_factor
         * ( EightmCXElMinmQsq * ElMin * mQCsqDiff * mVsq ) );
    // tick!

    if( 0.0 != same_handedness_factor )
      {

        /* antisymmetric bit positive for both vector & axial coupling (this is
         * the sum of 2 terms in the Mathematica notebook: VvPlusJjAaAllSq &
         * VvSqPlusAaSq):
         */
        /*current_term_coefficient
        += ( same_handedness_factor
             * ( -EightmCEVElMin * ElMin * mQsq * mQCsqDiff * mVsq
                 + 2.0 * ElMin * mQsq * mVsq
                   * ( EightmCEVElMin * mQsq
                       + 4.0 * ElMin * mCsq * mQCsqSum ) ) );*/
        current_term_coefficient
        += ( same_handedness_factor
             * ( 2.0 * ElMin * mQsq * mCsq * mVsq
                   * ( EightmCEVElMin + 4.0 * ElMin * mQCsqSum ) ) );
        // tick! tick!

        // antisymmetric bit positive for vector & negative for axial coupling:
        current_term_coefficient
        += ( same_handedness_factor * axial_coupling_factor
             * ( EightmCXElMinmQsq * ElMin * mQCsqSum * mVsq ) );
        // tick!

      }

    MINtoLK_log.set_coefficient( current_term_coefficient );
    if( 0.0 != current_term_coefficient )
      {

        normalization
        += ( current_term_coefficient
             * ( ( Elk * ( lnElk + 1.0 ) )
                 - ( minimumEnergy *  ( lnEmin + 1.0 ) ) ) );

      }


    // lin:

    // symmetric bit positive for both vector & axial coupling:
    current_term_coefficient
    = ( 4.0 * EightmCEVElMin * ElMin * mQcu * mQCsqDiff );
    // tick!

    // symmetric bit positive for vector & negative for axial coupling:
    // 0.0
    // tick!

    if( 0.0 != same_handedness_factor )
      {

        /* antisymmetric bit positive for both vector & axial coupling (this is
         * the sum of 2 terms in the Mathematica notebook: VvPlusJjAaAllSq &
         * VvSqPlusAaSq):
         */
        /*current_term_coefficient
        += ( same_handedness_factor
             * ( 16.0 * ElMinsq * mQcu * mQCsqDiff * mVsq
                 - 2.0 * mQcu * mVsq
                   * ( EightmCEVElMin * mC
                       + 8.0 * ElMin * ( mCsq * mC
                                         + 2.0 * ElMin
                                           * ( ( lnElMin - lnmQC ) * mCsq
                                               + mQsq ) )
                       - 2.0 * mCsq * mVsq ) ) );*/
        current_term_coefficient
        += ( same_handedness_factor
             * ( 2.0 * mQcu * mVsq
                 * ( 8.0 * ElMinsq * mQCsqDiff
                     - EightmCEVElMin * mC
                     + 8.0 * ElMin
                       * ( mC * mCsq
                           + 2.0 * ElMin * ( ( lnElMin - lnmQC ) * mCsq
                                             + mQsq ) )
                           - 2.0 * mCsq * mVsq ) ) );
        // tick! tick!

        // antisymmetric bit positive for vector & negative for axial coupling:
        current_term_coefficient
        += ( same_handedness_factor * axial_coupling_factor
             * ( -2.0 * EightmCXElMinmQsq * mC * mQ * mVsq ) );
        // tick!

      }

    MINtoLK_lin.set_coefficient( current_term_coefficient );
    if( 0.0 != current_term_coefficient )
      {

        normalization
        += ( current_term_coefficient
             * ( 0.5 * ( Elk * Elk - minimumEnergy * minimumEnergy ) ) );

      }


    // linlog:

    // symmetric bit positive for both vector & axial coupling:
    // 0.0
    // tick!

    // symmetric bit positive for vector & negative for axial coupling:
    // 0.0
    // tick!

    if( 0.0 != same_handedness_factor )
      {

        /* antisymmetric bit positive for both vector & axial coupling (this is
         * the sum of 2 terms in the Mathematica notebook: VvPlusJjAaAllSq &
         * VvSqPlusAaSq):
         */
        /* 0.0
         * + ( same_handedness_factor
               * ( 32.0 * ElMinsq * mCsq * mQcu * mVsq ) )*/
        current_term_coefficient
        = ( same_handedness_factor * ( 32.0 * ElMinsq * mCsq * mQcu * mVsq ) );
        // tick! tick!
        normalization
        += ( current_term_coefficient
             * ( ( Elk * Elk * ( lnElk + 0.5 ) )
                 - ( minimumEnergy * minimumEnergy * ( lnEmin + 0.5 ) ) ) );
        // current_term_coefficient should not be zero at this point.


        // antisymmetric bit positive for vector & negative for axial coupling:
        // 0.0
        // tick!

      }
    else
      {

        current_term_coefficient = 0.0;

      }

    MINtoLK_linlog.set_coefficient( current_term_coefficient );


    // sq:

    // symmetric bit positive for both vector & axial coupling:
    current_term_coefficient
    = ( -16.0 * ElMinsq * mQsq * mQsq * mQCsqDiff );
    // tick!

    // symmetric bit positive for vector & negative for axial coupling:
    // 0.0
    // tick!

    /*if( 0.0 != same_handedness_factor )
      {

        /* antisymmetric bit positive for both vector & axial coupling (this is
         * the sum of 2 terms in the Mathematica notebook: VvPlusJjAaAllSq &
         * VvSqPlusAaSq):
         *//*
        // 0.0
        // tick! tick!

        // antisymmetric bit positive for vector & negative for axial coupling:
        // 0.0
        // tick!

      }*/

    MINtoLK_sq.set_coefficient( current_term_coefficient );
    normalization
    += ( current_term_coefficient
         * ( ( Elk * Elk * Elk
               - minimumEnergy * minimumEnergy * minimumEnergy ) / 3.0 ) );
    // current_term_coefficient should not be zero at this point.


    // LK to HK:
    LKtoHK_segment.set_segment_range( Elk,
                                      Ehk );

    if( energy_limited_cosines )
      // if the cosine limits are limited by kinematics...
      {

        // const:

        // symmetric bit positive for both vector & axial coupling:
        current_term_coefficient
        = ( mQsq * mQCsqDiff
            * ( EightmCEVElMin * mVsq
                * ( ElMin * ( lnElMin - lnVsqOverFourElMinsq ) + mC )
                - 4.0 * EightmCEVElMin * ElMinsq * mC
                + mCsq * ( 16.0 * ElMinsq * ElMinsq - mVsq * mVsq ) ) );
        // tick!

        // symmetric bit positive for vector & negative for axial coupling:
        current_term_coefficient
        += ( axial_coupling_factor
             * ( EightmCXElMinmQsq * ElMin * mQCsqDiff * mVsq
                 * ( lnVsqOverFourElMinsq - lnElMin ) ) );
        // tick!

        if( 0.0 != same_handedness_factor )
          {

            /* antisymmetric bit positive for both vector & axial coupling
             * (this is the sum of 2 terms in the Mathematica notebook:
             * VvPlusJjAaAllSq & VvSqPlusAaSq):
             */
            /*current_term_coefficient
            += ( same_handedness_factor
                 * ( ElMin * mQsq * mQCsqDiff * mVsq
                     * ( EightmCEVElMin * ( lnElMin - lnVsqOverFourElMinsq )
                         + 4.0 * mC * ( mVsq - 4.0 * ElMinsq ) )
                     + 2.0 * ElMin * mQsq * mVsq
                       * ( EightmCEVElMin * mQsq
                           * ( lnVsqOverFourElMinsq - lnElMin )
                           + 8.0 * ElMinsq * mC * ( mCsq + 3.0 * mQsq )
                           - 4.0 * ElMin * mCsq * mQCsqSum
                             * ( lnElMin - lnVsqOverFourElMinsq )
                           - 2.0 * mC * ( mCsq + 3.0 * mQsq ) * mVsq ) ) );*/
            current_term_coefficient
            += ( same_handedness_factor
                 * ( ElMin * mQsq * mQCsqSum * mVsq
                     * ( ( EightmCEVElMin + 8.0 * ElMin * mCsq )
                         * ( lnVsqOverFourElMinsq - lnElMin )
                         + 8.0 * mC * ( 4.0 * ElMinsq - mVsq ) ) ) );
            // tick! tick!

            // antisymmetric bit positive for vector & negative for axial
            // coupling:
            current_term_coefficient
            += ( same_handedness_factor * axial_coupling_factor
                 * ( EightmCXElMinmQsq * ElMin * mQCsqSum * mVsq
                     * ( lnVsqOverFourElMinsq - lnElMin ) ) );
            // tick!

          }

        LKtoHK_const.set_coefficient( current_term_coefficient );


        // lin:

        // symmetric bit positive for both vector & axial coupling:
        // 0.0
        // tick!

        // symmetric bit positive for vector & negative for axial coupling:
        // 0.0
        // tick!

        if( 0.0 != same_handedness_factor )
          {

            /* antisymmetric bit positive for both vector & axial coupling
             * (this is the sum of 2 terms in the Mathematica notebook:
             * VvPlusJjAaAllSq & VvSqPlusAaSq):
             */
            /* 0.0
               + ( same_handedness_factor
                 * ( 2.0 * mC * mQcu
                     * ( 4.0 * ElMinsq
                         * ( EightmCEVElMin
                             + 8.0 * ElMin * mC * ( mC - ElMin ) )
                         - mVsq * ( EightmCEVElMin
                                    + 8.0 * ElMin * mC
                                      * ( 2.0 * ElMin
                                          * ( lnElMin - lnVsqOverFourElMinsq )
                                          + mC ) )
                       - 2.0 * mC * mVsq * mVsq ) ) )
             */
            current_term_coefficient
            = ( same_handedness_factor
                * ( 2.0 * mC * mQcu
                    * ( 4.0 * ElMinsq
                        * ( EightmCEVElMin
                            + 8.0 * ElMin * mC * ( mC - ElMin ) )
                        - mVsq * ( EightmCEVElMin
                                   + 8.0 * ElMin * mC
                                     * ( 2.0 * ElMin
                                         * ( lnElMin - lnVsqOverFourElMinsq )
                                             + mC ) )
                        + 2.0 * mC * mVsq * mVsq ) ) );
            // tick! tick!

            // antisymmetric bit positive for vector & negative for axial coupling:
            current_term_coefficient
            += ( same_handedness_factor * axial_coupling_factor
                 * ( 2.0 * EightmCXElMinmQsq * mC * mQ
                     * (4.0 * ElMinsq - mVsq ) ) );
            // tick!

          }
        else
          {

            current_term_coefficient = 0.0;

          }

        LKtoHK_lin.set_coefficient( current_term_coefficient );


        // sq:

        // symmetric bit positive for both vector & axial coupling:
        // 0.0
        // tick!

        // symmetric bit positive for vector & negative for axial coupling:
        // 0.0
        // tick!

        /*if( 0.0 != same_handedness_factor )
          {

            /* antisymmetric bit positive for both vector & axial coupling
             * (this is the sum of 2 terms in the Mathematica notebook:
             * VvPlusJjAaAllSq & VvSqPlusAaSq):
             *//*
            // 0.0
            // tick! tick!

            // antisymmetric bit positive for vector & negative for axial coupling:
            // 0.0
            // tick!

          }*/

        LKtoHK_sq.set_coefficient( 0.0 );

      }
    else
      // otherwise, if the cosine limits are not limited by kinematics...
      {

        // const:

        // symmetric bit positive for both vector & axial coupling:
        current_term_coefficient
        = ( -2.0 * EightmCEVElMin * ElMin * lnmQC * mQsq * mQCsqDiff * mVsq );
        // tick!

        // symmetric bit positive for vector & negative for axial coupling:
        current_term_coefficient
        += ( axial_coupling_factor
            * ( 2.0 * EightmCXElMinmQsq * ElMin * lnmQC * mQCsqDiff * mVsq ) );
        // tick!

        if( 0.0 != same_handedness_factor )
          {

            /* antisymmetric bit positive for both vector & axial coupling
             * (this is the sum of 2 terms in the Mathematica notebook:
             * VvPlusJjAaAllSq & VvSqPlusAaSq):
             */
            /*current_term_coefficient
            += ( same_handedness_factor
                 * ( -2.0 * EightmCEVElMin* ElMin * lnmQC * mQsq * mQCsqDiff
                     * mVsq
                     + 2.0 * ElMin * mQsq * mVsq
                       * ( EightmCEVElMin * mCsq
                           + 8.0 * ElMin * ( 1.0 + lnmQC ) * mCsq * mCsq
                           + ( EightmCEVElMin * ( 2.0 * lnmQC - 1.0 )
                               + 8.0 * ElMin * mQsq
                                 * ( lnmQC - 1.0 ) * mCsq ) ) ) );*/
            current_term_coefficient
            += ( same_handedness_factor
                 * ( 2.0 * ElMin * mQsq * mVsq
                     * ( EightmCEVElMin + 8.0 * ElMin * mCsq )
                     * ( lnmQC * mQCsqSum - mQCsqDiff ) ) );
            // tick! tick!

            // antisymmetric bit positive for vector & negative for axial
            // coupling:
            current_term_coefficient
            += ( same_handedness_factor * axial_coupling_factor
                 * ( 2.0 * EightmCXElMinmQsq * ElMin * mVsq
                     * ( lnmQC * mQCsqSum - mQCsqDiff ) ) );
            // tick!

          }

        LKtoHK_const.set_coefficient( current_term_coefficient );


        // lin:

        // symmetric bit positive for both vector & axial coupling:
        current_term_coefficient
        = ( 4.0 * EightmCEVElMin * ElMin * mQ * mQCsqDiff * mQCsqDiff );
        // tick!

        // symmetric bit positive for vector & negative for axial coupling:
        // 0.0
        // tick!

        if( 0.0 != same_handedness_factor )
          {

            /* antisymmetric bit positive for both vector & axial coupling
             * (this is the sum of 2 terms in the Mathematica notebook:
             * VvPlusJjAaAllSq & VvSqPlusAaSq):
             */
            /*current_term_coefficient
            += ( same_handedness_factor
                 * ( 16.0 * ElMinsq * mQ * mQCsqDiff * mQCsqDiff * mVsq
                     + 32.0 * ElMinsq * mQcu * mVsq
                       * ( 2.0 * mCsq * lnmQC - mQCsqDiff ) ) );*/
            current_term_coefficient
            += ( same_handedness_factor
                 * ( 16.0 * ElMinsq * mQ * mVsq
                     * ( 4.0 * mQsq * mCsq * lnmQC
                         - mQCsqDiff * mQCsqSum ) ) );
            // tick! tick!

            // antisymmetric bit positive for vector & negative for axial
            // coupling:
            // 0.0
            // tick!

          }

        LKtoHK_lin.set_coefficient( current_term_coefficient );


        // sq:

        // symmetric bit positive for both vector & axial coupling:
        current_term_coefficient
        = ( -16.0 * ElMinsq * mQCsqDiff * mQCsqDiff * mQCsqSum );
        // tick!

        // symmetric bit positive for vector & negative for axial coupling:
        // 0.0
        // tick!

        /*if( 0.0 != same_handedness_factor )
          {

            /* antisymmetric bit positive for both vector & axial coupling
             * (this is the sum of 2 terms in the Mathematica notebook:
             * VvPlusJjAaAllSq & VvSqPlusAaSq):
             *//*
            // 0.0
            // tick! tick!

            // antisymmetric bit positive for vector & negative for axial coupling:
            // 0.0
            // tick!

          }*/

        LKtoHK_sq.set_coefficient( current_term_coefficient );

      }
    current_term_coefficient = LKtoHK_sq.get_coefficient();
    if( 0.0 != current_term_coefficient )
      {

        normalization
        += ( current_term_coefficient
            * ( ( Ehk * Ehk * Ehk - Elk * Elk * Elk ) / 3.0 ) );

      }
    current_term_coefficient = LKtoHK_lin.get_coefficient();
    if( 0.0 != current_term_coefficient )
      {

        normalization
        += ( current_term_coefficient * ( 0.5 * ( Ehk * Ehk - Elk * Elk ) ) );

      }
    current_term_coefficient = LKtoHK_const.get_coefficient();
    if( 0.0 != current_term_coefficient )
      {

        normalization
        += ( current_term_coefficient * ( Ehk - Elk ) );

      }


    // HK to MAX:
    HKtoMAX_segment.set_segment_range( Ehk,
                                       maximumEnergy );

    // inv:

    // symmetric bit positive for both vector & axial coupling:
    // 0.0
    // tick!

    // symmetric bit positive for vector & negative for axial coupling:
    // 0.0
    // tick!

    if( 0.0 != same_handedness_factor )
      {

        /* antisymmetric bit positive for both vector & axial coupling (this is
         * the sum of 2 terms in the Mathematica notebook: VvPlusJjAaAllSq &
         * VvSqPlusAaSq):
         */
        /* 0.0
           + ( same_handedness_factor
               * ( -4.0 * ElMinsq * mCsq * mQcu * mVsq * mVsq ) );*/
        current_term_coefficient
        = ( same_handedness_factor
            * ( -4.0 * ElMinsq * mCsq * mQcu * mVsq * mVsq ) );
        // tick! tick!
        normalization
        += ( current_term_coefficient * ( lnEmax - lnEhk ) );
        // current_term_coefficient should not be zero at this point.

        // antisymmetric bit positive for vector & negative for axial coupling:
        // 0.0
        // tick!

      }
    else
      {

        current_term_coefficient = 0.0;

      }

    HKtoMAX_inv.set_coefficient( current_term_coefficient );


    // const:

    // symmetric bit positive for both vector & axial coupling:
    current_term_coefficient
    = ( mQsq * mQCsqDiff * mVsq
        * ( EightmCEVElMin * ( mC - ElMin * ( lnmQC + lnVsqOverFourElMinsq ) )
            - mCsq * mVsq ) );
    // tick!

    // symmetric bit positive for vector & negative for axial coupling:
    current_term_coefficient
    += ( axial_coupling_factor
        * ( EightmCXElMinmQsq * ElMin * mQCsqDiff * mVsq
            * ( lnmQC + lnVsqOverFourElMinsq ) ) );
    // tick!

    if( 0.0 != same_handedness_factor )
      {

        /* antisymmetric bit positive for both vector & axial coupling (this is
         * the sum of 2 terms in the Mathematica notebook: VvPlusJjAaAllSq &
         * VvSqPlusAaSq):
         */
        /*current_term_coefficient
        += ( same_handedness_factor
            * ( ElMin * mQsq * mQCsqDiff * mVsq
                * ( 4.0 * mC * mVsq
                    - EightmCEVElMin * ( lnmQC + lnVsqOverFourElMinsq ) )
                + 2.0 * ElMin * mQsq * mVsq
                  * ( 4.0 * ElMin * mCsq
                      * ( ( lnmQC + lnVsqOverFourElMinsq ) * mCsq
                          + ( lnmQC - 2.0 + lnVsqOverFourElMinsq ) * mQsq )
                      + 8.0 * ElMinsq * mC * mQCsqSum
                      + mQsq * ( EightmCEVElMin
                                 * ( lnmQC - 1.0 + lnVsqOverFourElMinsq )
                                 - 4.0 * mC * mVsq ) ) ) );*/
        current_term_coefficient
        += ( same_handedness_factor
            * ( ElMin * mQsq * mVsq
                * ( ( lnmQC + lnVsqOverFourElMinsq ) * mQCsqSum
                    * ( EightmCEVElMin + 8.0 * ElMin * mCsq )
                    + 4.0 * mC
                      * ( 4.0 * ElMin * ( ElMin * mQCsqSum - mC * mQsq )
                          - mQCsqSum * mVsq ) ) ) );
        // tick! tick!

        // antisymmetric bit positive for vector & negative for axial coupling:
        current_term_coefficient
        += ( same_handedness_factor * axial_coupling_factor
            * ( EightmCXElMinmQsq * ElMin * mVsq
                * ( ( lnmQC + lnVsqOverFourElMinsq ) * mQCsqSum
                    - 2.0 * mQsq ) ) );
        // tick!

      }

    HKtoMAX_const.set_coefficient( current_term_coefficient );
    if( 0.0 != current_term_coefficient )
      {

        normalization
        += ( current_term_coefficient * ( maximumEnergy - Ehk ) );

      }


    // log:

    // symmetric bit positive for both vector & axial coupling:
    current_term_coefficient
    = ( EightmCEVElMin * ElMin * mQsq * mQCsqDiff * mVsq );
    // tick!

    // symmetric bit positive for vector & negative for axial coupling:
    current_term_coefficient
    += ( axial_coupling_factor
         * ( -EightmCXElMinmQsq * ElMin * mQCsqDiff * mVsq ) );
    // tick!

    if( 0.0 != same_handedness_factor )
      {

        /* antisymmetric bit positive for both vector & axial coupling (this is
         * the sum of 2 terms in the Mathematica notebook: VvPlusJjAaAllSq &
         * VvSqPlusAaSq):
         */
        /*current_term_coefficient
        += ( same_handedness_factor
             * ( EightmCEVElMin * ElMin * mQsq * mQCsqDiff * mVsq
                 - 2.0 * ElMin * mQsq * mVsq
                   * ( EightmCEVElMin * mQsq
                       + 4.0 * ElMin * mCsq * mQCsqSum ) ) );*/
        current_term_coefficient
        += ( same_handedness_factor
             * ( -ElMin * mQsq * mQCsqSum * mVsq
                 * ( EightmCEVElMin + 8.0 * ElMin * mCsq ) ) );
        // tick! tick!

        // antisymmetric bit positive for vector & negative for axial coupling:
        current_term_coefficient
        += ( same_handedness_factor * axial_coupling_factor
             * ( -EightmCXElMinmQsq * ElMin * mQCsqSum * mVsq ) );
        // tick!

      }

    HKtoMAX_log.set_coefficient( current_term_coefficient );
    if( 0.0 != current_term_coefficient )
      {

        normalization
        += ( current_term_coefficient
             * ( ( maximumEnergy *  ( lnEmax + 1.0 ) )
                 - ( Ehk * ( lnEhk + 1.0 ) ) ) );

      }


    // lin:

    // symmetric bit positive for both vector & axial coupling:
    current_term_coefficient
    = ( -4.0 * EightmCEVElMin * ElMin * mCsq * mQ * mQCsqDiff );
    // tick!

    // symmetric bit positive for vector & negative for axial coupling:
    // 0.0
    // tick!

    if( 0.0 != same_handedness_factor )
      {

        /* antisymmetric bit positive for both vector & axial coupling (this is
         * the sum of 2 terms in the Mathematica notebook: VvPlusJjAaAllSq &
         * VvSqPlusAaSq):
         */
        /*current_term_coefficient
        += ( same_handedness_factor
            * ( -16.0 * ElMinsq * mCsq * mQ * mQCsqDiff * mVsq
                + 8.0 * ElMinsq * mC * mQcu
                  * ( EightmCEVElMin + 8.0 * ElMin * mC * ( mC - ElMin )
                      + 4.0 * mC * mVsq
                        * ( 1.0 + lnmQC + lnVsqOverFourElMinsq ) ) ) );*/
        current_term_coefficient
        += ( same_handedness_factor
            * ( 8.0 * ElMinsq * mC * mQ
                * ( 2.0 * mC * mVsq
                    * ( 2.0 * mQsq * ( lnmQC + lnVsqOverFourElMinsq )
                        + mQCsqSum )
                    + 8.0 * mQsq
                      * ( EightmCEVElMin
                          + 8.0 * ElMin * mC * ( mC - ElMin ) ) ) ) );
        // tick! tick!

        // antisymmetric bit positive for vector & negative for axial coupling:
        current_term_coefficient
        += ( same_handedness_factor * axial_coupling_factor
             * ( 8.0 * EightmCXElMinmQsq * ElMinsq * mC * mQ ) );
        // tick!

      }

    HKtoMAX_lin.set_coefficient( current_term_coefficient );
    if( 0.0 != current_term_coefficient )
      {

        normalization
        += ( current_term_coefficient
             * ( 0.5 * ( maximumEnergy * maximumEnergy - Ehk * Ehk ) ) );

      }


    // linlog:

    // symmetric bit positive for both vector & axial coupling:
    // 0.0
    // tick!

    // symmetric bit positive for vector & negative for axial coupling:
    // 0.0
    // tick!

    if( 0.0 != same_handedness_factor )
      {

        /* antisymmetric bit positive for both vector & axial coupling (this is
         * the sum of 2 terms in the Mathematica notebook: VvPlusJjAaAllSq &
         * VvSqPlusAaSq):
         */
        /* 0.0
           + ( same_handedness_factor
               * ( -32.0 * ElMinsq * mCsq * mQcu * mVsq ) );*/
        current_term_coefficient
        = ( same_handedness_factor
            * ( -32.0 * ElMinsq * mCsq * mQcu * mVsq ) );
        // tick! tick!
        normalization
        += ( current_term_coefficient
             * ( ( maximumEnergy * maximumEnergy * ( lnEmax + 0.5 ) )
                 - ( Ehk * Ehk * ( lnEhk + 0.5 ) ) ) );
        // current_term_coefficient should not be zero at this point.

        // antisymmetric bit positive for vector & negative for axial coupling:
        // 0.0
        // tick!

      }
    else
      {

        current_term_coefficient = 0.0;

      }

    HKtoMAX_linlog.set_coefficient( current_term_coefficient );


    // sq:

    // symmetric bit positive for both vector & axial coupling:
    current_term_coefficient = ( 16.0 * ElMinsq * mCsq * mCsq * mQCsqDiff );
    // tick!

    // symmetric bit positive for vector & negative for axial coupling:
    // 0.0
    // tick!

    /*if( 0.0 != same_handedness_factor )
      {

        /* antisymmetric bit positive for both vector & axial coupling (this is
         * the sum of 2 terms in the Mathematica notebook: VvPlusJjAaAllSq &
         * VvSqPlusAaSq):
         *//*
        // 0.0
    // tick! tick!

        // antisymmetric bit positive for vector & negative for axial coupling:
        // 0.0
    // tick!

      }*/

    HKtoMAX_sq.set_coefficient( current_term_coefficient );
    normalization
    += ( current_term_coefficient
         * ( ( maximumEnergy * maximumEnergy * maximumEnergy
               - Ehk * Ehk * Ehk ) / 3.0 ) );
    // current_term_coefficient should not be zero at this point.

    // now we normalize the distribution:
    if( 0.0 != normalization )
      {

        normalizeCoefficients();

      }

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



  W_minus_handed_muon::W_minus_handed_muon(
                                    readier_for_new_point* const given_readier,
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
    MINtoLK_inv( -1,
                 0,
                 CppSLHA::CppSLHA_global::really_wrong_value ),
    MINtoLK_const( 0,
                   0,
                   CppSLHA::CppSLHA_global::really_wrong_value ),
    MINtoLK_log( 0,
                 1,
                 CppSLHA::CppSLHA_global::really_wrong_value ),
    MINtoLK_lin( 1,
                 0,
                 CppSLHA::CppSLHA_global::really_wrong_value ),
    MINtoLK_linlog( 1,
                    1,
                    CppSLHA::CppSLHA_global::really_wrong_value ),
    MINtoLK_sq( 2,
                0,
                CppSLHA::CppSLHA_global::really_wrong_value ),
    LKtoHK_const( 0,
                  0,
                  CppSLHA::CppSLHA_global::really_wrong_value ),
    LKtoHK_lin( 1,
                0,
                CppSLHA::CppSLHA_global::really_wrong_value ),
    LKtoHK_sq( 2,
               0,
               CppSLHA::CppSLHA_global::really_wrong_value ),
    HKtoMAX_inv( -1,
                 0,
                 CppSLHA::CppSLHA_global::really_wrong_value ),
    HKtoMAX_const( 0,
                   0,
                   CppSLHA::CppSLHA_global::really_wrong_value ),
    HKtoMAX_log( 0,
                 1,
                 CppSLHA::CppSLHA_global::really_wrong_value ),
    HKtoMAX_lin( 1,
                 0,
                 CppSLHA::CppSLHA_global::really_wrong_value ),
    HKtoMAX_linlog( 1,
                    1,
                    CppSLHA::CppSLHA_global::really_wrong_value ),
    HKtoMAX_sq( 2,
                0,
                CppSLHA::CppSLHA_global::really_wrong_value )
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

    MINtoLK_segment.add_term( &MINtoLK_inv );
    MINtoLK_segment.add_term( &MINtoLK_const );
    MINtoLK_segment.add_term( &MINtoLK_log );
    MINtoLK_segment.add_term( &MINtoLK_lin );
    MINtoLK_segment.add_term( &MINtoLK_linlog );
    MINtoLK_segment.add_term( &MINtoLK_sq );

    segments.push_back( &MINtoLK_segment );


    LKtoHK_segment.add_term( &LKtoHK_const );
    LKtoHK_segment.add_term( &LKtoHK_lin );
    LKtoHK_segment.add_term( &LKtoHK_sq );

    segments.push_back( &LKtoHK_segment );


    HKtoMAX_segment.add_term( &HKtoMAX_inv );
    HKtoMAX_segment.add_term( &HKtoMAX_const );
    HKtoMAX_segment.add_term( &HKtoMAX_log );
    HKtoMAX_segment.add_term( &HKtoMAX_lin );
    HKtoMAX_segment.add_term( &HKtoMAX_linlog );
    HKtoMAX_segment.add_term( &HKtoMAX_sq );

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

    normalization = 0.0;
    lnEmin = log( minimumEnergy );
    lnElk = log( Elk );
    lnEhk = log( Ehk );
    lnEmax = log( maximumEnergy );

    // MIN to LK:
    MINtoLK_segment.set_segment_range( minimumEnergy,
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

    MINtoLK_inv.set_coefficient( current_term_coefficient );
    if( 0.0 != current_term_coefficient )
      {

        normalization = ( current_term_coefficient * ( lnElk - lnEmin ) );

      }
    else
      {

        normalization = 0.0;

      }


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

    MINtoLK_const.set_coefficient( current_term_coefficient );
    if( 0.0 != current_term_coefficient )
      {

        normalization
        += ( current_term_coefficient * ( Elk - minimumEnergy ) );

      }


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

    MINtoLK_log.set_coefficient( current_term_coefficient );
    if( 0.0 != current_term_coefficient )
      {

        normalization
        += ( current_term_coefficient
             * ( ( Elk * ( lnElk + 1.0 ) )
                 - ( minimumEnergy *  ( lnEmin + 1.0 ) ) ) );

      }


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

    MINtoLK_lin.set_coefficient( current_term_coefficient );
    if( 0.0 != current_term_coefficient )
      {

        normalization
        += ( current_term_coefficient
             * ( 0.5 * ( Elk * Elk - minimumEnergy * minimumEnergy ) ) );

      }


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

    MINtoLK_linlog.set_coefficient( current_term_coefficient );
    if( 0.0 != current_term_coefficient )
      {

        normalization
        += ( current_term_coefficient
             * ( ( Elk * Elk * ( lnElk + 0.5 ) )
                 - ( minimumEnergy * minimumEnergy * ( lnEmin + 0.5 ) ) ) );

      }


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

    MINtoLK_sq.set_coefficient( current_term_coefficient );
    if( 0.0 != current_term_coefficient )
      {

        normalization
        += ( current_term_coefficient
             * ( ( Elk * Elk * Elk
                   - minimumEnergy * minimumEnergy * minimumEnergy )
                   / 3.0 ) );

      }


    // LK to HK:
    LKtoHK_segment.set_segment_range( Elk,
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

        LKtoHK_const.set_coefficient( current_term_coefficient );


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

        LKtoHK_lin.set_coefficient( current_term_coefficient );


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

        LKtoHK_sq.set_coefficient( 0.0 );

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

        LKtoHK_const.set_coefficient( current_term_coefficient );


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

        LKtoHK_lin.set_coefficient( current_term_coefficient );


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

        LKtoHK_sq.set_coefficient( current_term_coefficient );

      }
    current_term_coefficient = LKtoHK_sq.get_coefficient();
    if( 0.0 != current_term_coefficient )
      {

        normalization
        += ( current_term_coefficient
            * ( ( Ehk * Ehk * Ehk - Elk * Elk * Elk ) / 3.0 ) );

      }
    current_term_coefficient = LKtoHK_lin.get_coefficient();
    if( 0.0 != current_term_coefficient )
      {

        normalization
        += ( current_term_coefficient * ( 0.5 * ( Ehk * Ehk - Elk * Elk ) ) );

      }
    current_term_coefficient = LKtoHK_const.get_coefficient();
    if( 0.0 != current_term_coefficient )
      {

        normalization
        += ( current_term_coefficient * ( Ehk - Elk ) );

      }


    // HK to MAX:
    HKtoMAX_segment.set_segment_range( Ehk,
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

    HKtoMAX_inv.set_coefficient( current_term_coefficient );
    if( 0.0 != current_term_coefficient )
      {

        normalization
        += ( current_term_coefficient * ( lnEmax - lnEhk ) );

      }


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

    HKtoMAX_const.set_coefficient( current_term_coefficient );
    if( 0.0 != current_term_coefficient )
      {

        normalization
        += ( current_term_coefficient * ( maximumEnergy - Ehk ) );

      }


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

    HKtoMAX_log.set_coefficient( current_term_coefficient );
    if( 0.0 != current_term_coefficient )
      {

        normalization
        += ( current_term_coefficient
             * ( ( maximumEnergy *  ( lnEmax + 1.0 ) )
                 - ( Ehk * ( lnEhk + 1.0 ) ) ) );

      }


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

    HKtoMAX_lin.set_coefficient( current_term_coefficient );
    if( 0.0 != current_term_coefficient )
      {

        normalization
        += ( current_term_coefficient
             * ( 0.5 * ( maximumEnergy * maximumEnergy - Ehk * Ehk ) ) );

      }


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

    HKtoMAX_linlog.set_coefficient( current_term_coefficient );
    if( 0.0 != current_term_coefficient )
      {

        normalization
        += ( current_term_coefficient
             * ( ( maximumEnergy * maximumEnergy * ( lnEmax + 0.5 ) )
                 - ( Ehk * Ehk * ( lnEhk + 0.5 ) ) ) );

      }


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

    HKtoMAX_sq.set_coefficient( current_term_coefficient );
    if( 0.0 != current_term_coefficient )
      {

        normalization
            += ( current_term_coefficient
                 * ( ( maximumEnergy * maximumEnergy * maximumEnergy
                       - Ehk * Ehk * Ehk ) / 3.0 ) );

      }


    // now we normalize the distribution:
    if( 0.0 != normalization )
      {

        normalizeCoefficients();

      }

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



  Z_direct_jet::Z_direct_jet( readier_for_new_point* const given_readier,
                              CppSLHA::CppSLHA0 const* const given_spectrum,
              CppSLHA::particle_property_set const* const given_first_particle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const given_second_particle,
              CppSLHA::particle_property_set const* const given_third_particle,
          CppSLHA::particle_property_set const* const given_fourth_particle ) :
    Z_handed_muon( given_readier,
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
                                    readier_for_new_point* const given_readier,
                                 CppSLHA::CppSLHA0 const* const given_spectrum,
              CppSLHA::particle_property_set const* const given_first_particle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
             CppSLHA::particle_property_set const* const given_second_particle,
              CppSLHA::particle_property_set const* const given_third_particle,
          CppSLHA::particle_property_set const* const given_fourth_particle ) :
    Z_handed_muon( given_readier,
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
                                    readier_for_new_point* const given_readier,
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
                                    readier_for_new_point* const given_readier,
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
                                    readier_for_new_point* const given_readier,
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
    MINtoMAX_const( 0,
                    0,
                    CppSLHA::CppSLHA_global::really_wrong_value )
  {

    MINtoMAX_segment.add_term( &MINtoMAX_const );

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


    // now we set up the coefficients of various terms, & use them to
    // determine the normalization:

    normalization = maximumEnergy;


    // now we set up the coefficients of various terms, & use them to determine
    // the normalization:

    // MIN to MAX segment:
    MINtoMAX_segment.set_segment_range( minimumEnergy,
                                        maximumEnergy );

    MINtoMAX_const.set_coefficient( ( 1.0 / normalization ) );

  }



  chargino_three_body_decay::chargino_three_body_decay(
                                    readier_for_new_point* const given_readier,
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
    MINtoMAX_const( 0,
                    0,
                    CppSLHA::CppSLHA_global::really_wrong_value )
  {

    MINtoMAX_segment.add_term( &MINtoMAX_const );

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


    // now we set up the coefficients of various terms, & use them to
    // determine the normalization:

    normalization = maximumEnergy;


    // now we set up the coefficients of various terms, & use them to determine
    // the normalization:

    // MIN to MAX segment:
    MINtoMAX_segment.set_segment_range( minimumEnergy,
                                        maximumEnergy );

    MINtoMAX_const.set_coefficient( ( 1.0 / normalization ) );

  }

}  // end of LHC_FASER namespace.

