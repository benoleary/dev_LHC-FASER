/*
 * LHC-FASER_global_stuff.cpp
 *
 *  Created on: 12 Jan 2011
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

#include "LHC-FASER_global_stuff.hpp"


namespace LHC_FASER
{

  // LHC-FASER_global_stuff non-inline functions:

  double
  LHC_FASER_global::unit_linear_interpolation( double const x_fraction,
                                               double const left_value,
                                               double const right_value )
  /* this does a linear interpolation from the values at the given 2
   * ends of a unit line segment for the given point in the segment.
   */
  {

    return ( left_value + ( right_value - left_value ) * x_fraction );

  }

  double
  LHC_FASER_global::square_bilinear_interpolation( double const x_fraction,
                                                   double const y_fraction,
                                                 double const lower_left_value,
                                                double const lower_right_value,
                                                double const upper_right_value,
                                                double const upper_left_value )
  {

    if( 0.0 == x_fraction )
      {

        return unit_linear_interpolation( y_fraction,
                                          lower_left_value,
                                          upper_left_value );

      }
    else if( 1.0 == x_fraction )
      {

        return unit_linear_interpolation( y_fraction,
                                          lower_right_value,
                                          upper_right_value );

      }
    else if( 0.0 == y_fraction )
      {

        return unit_linear_interpolation( x_fraction,
                                          lower_left_value,
                                          lower_right_value );

      }
    else if( 1.0 == y_fraction )
      {

        return unit_linear_interpolation( x_fraction,
                                          upper_left_value,
                                          upper_right_value );

      }
    else
      {

        return unit_linear_interpolation( x_fraction,
                                         unit_linear_interpolation( y_fraction,
                                                              lower_left_value,
                                                            upper_left_value ),
                                         unit_linear_interpolation( y_fraction,
                                                             lower_right_value,
                                                         upper_right_value ) );

      }

  }

  // this is where we decide what value marks when we throw away cascades
  // because the acceptance times branching ratio dropped too low:
  double const LHC_FASER_global::negligible_BR = 0.001;
  //double const LHC_FASER_global::negligible_BR = 0.000001;
  /* 1 millionth seems conservative enough. even at design luminosity, no
   * MSSM point could generate more than a few million electroweakinos in a
   * year.
   */

  /* this is where we decide what value marks when we throw away channels
   * because the cross-section times acceptance times branching ratio dropped
   * too low:
   */
  double const LHC_FASER_global::negligible_sigma = 0.001;
  /* internally, picobarns are used as the unit of cross-section. if we would
   * not see at least 1 event from this channel in the 1st 1/fb of data (so
   * 1000 events in the 1st 1/ab of data, 3 years of LHC running at design
   * luminosity), the channel's cross-section is deemed negligible.
   */



  readier_for_new_point::readier_for_new_point()
  {

    // does nothing.

  }

  readier_for_new_point::~readier_for_new_point()
  {

    // the readier_for_new_point has to let all its obsevers know that its
    // destructor is being called:
    for( std::list< bool* >::iterator
         observed_iterator = observer_bools_for_readier_existence.begin();
         observer_bools_for_readier_existence.end() != observed_iterator;
         ++observed_iterator )
      // go through the list of observed bools...
      {

        *(*observed_iterator) = false;
        // this should let any remaining observers know that they should *not*
        // attempt to de-register themselves from this readier_for_new_point.

      }

  }


  void
  readier_for_new_point::remove_me( bool* observer_bool_for_readier_existence,
                                    bool* observer_bool_for_readying )
  /* this removes the observer from the readier's list of observers which is
   * used when the readier's destructor is called to let the observers know
   * to stop asking to be reset.
   */
  {

    bool found_requested = false;

    for( std::list< bool* >::iterator
         bool_iterator = observer_bools_for_readier_existence.begin();
         ( ( observer_bools_for_readier_existence.end() != bool_iterator )
           &&
           !found_requested );
         ++bool_iterator )
      // look through the list for observer knowledge of readier existence...
      {

        if( observer_bool_for_readier_existence == *bool_iterator )
          // if we find the requested pointer...
          {

            bool_iterator
            = observer_bools_for_readier_existence.erase( bool_iterator );
            // remove the requested pointer from the list.
            found_requested = true;
            // stop looking (there should only be 1 instance of this pointer).

          }

      }

    found_requested = false;

    for( std::list< bool* >::iterator
         bool_iterator = observer_bools_for_readying.begin();
         ( ( observer_bools_for_readying.end() != bool_iterator )
           &&
           !found_requested );
         ++bool_iterator )
      // look through the list for readying bools...
      {

        if( observer_bool_for_readying == *bool_iterator )
          // if we find the requested pointer...
          {

            bool_iterator
            = observer_bools_for_readying.erase( bool_iterator );
            // remove the requested pointer from the list.
            found_requested = true;
            // stop looking (there should only be 1 instance of this pointer).

          }

      }

  }


  void
  readier_for_new_point::ready_for_new_point()
  // this goes through the list of pointers to bools created by observers
  // calling reset_me_on_next_update() & sets all the bools to true.
  {

    // debugging:
    /*std::cout
    << std::endl
    << "debugging: readier_for_new_point::readier_for_new_point is readying "
    << observer_bools_for_readying.size() << " bool pairs.";
    std::cout << std::endl;*/

    for( std::list< bool* >::iterator
         bool_iterator = observer_bools_for_readying.begin();
         observer_bools_for_readying.end() != bool_iterator;
         ++bool_iterator )
      // go through the list of observed bools...
      {

        // set the observer's bool to true so that it knows that an update has
        // happened:
        *(*bool_iterator) = true;

      }

    observer_bools_for_readying.clear();
    // the list needs to be reset afterwards. pointers to bools are added in
    // again by observers calling reset_me_on_next_update().

  }



  readied_for_new_point::readied_for_new_point(
                                 readier_for_new_point* const given_readier ) :
    needs_to_prepare_flag( true ),
    readier_exists_flag( true ),
    readier( given_readier )
  /* the constructor gives a pointer to a bool pair to given_readier which
   * has values that given_readier changes. this object keeps a pointer to
   * given_readier so that it can de-register when its destructor is called.
   */
  {

    given_readier->include_me( &readier_exists_flag );

  }

  readied_for_new_point::~readied_for_new_point()
  // the destructor tells readier, if it still exists, to stop modifying its
  // bools.
  {

    if( readier_exists_flag )
      // if the readier still exists...
      {

        readier->remove_me( &readier_exists_flag,
                            &needs_to_prepare_flag );

      }

  }

}  // end of LHC_FASER namespace.

