/*
 * LHC-FASER_cross-section_stuff.cpp
 *
 *  Created on: 15 Nov 2010
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

#include "LHC-FASER_cross-section_stuff.hpp"

namespace LHC_FASER
{

  square_grid::square_grid( std::string const* const grid_file_location,
                            std::string const* const given_grid_name,
                            square_grid const* const scaling_grid ) :
    grid_name( *given_grid_name )
  /* this constructor reads in a grid file, assumed to be in the format
   * x_coordinate y_coordinate value
   * in ascending order, y_coordinate varying first
   * (e.g.
   * 200.0 200.0 123.4
   * 200.0 210.0 567.8
   * ...
   * 210.0 200.0 876.5
   * 210.0 210.0 432.1
   * ...).
   * if provided with a non-NULL pointer to another square_grid, it stores
   * the its values scaled by an interpolated value from this other
   * square_grid.
   */
  {

    read_in( grid_file_location,
             scaling_grid );

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: square_grid() finished, value at ( 200.0, 200.0 ) = "
    << value_at( 200.0,
                 200.0 );
    std::cout << std::endl;**/

  }

  square_grid::square_grid( std::string base_grid_file_location,
                            std::string given_grid_name,
                            std::string scaling_grid_file_location ) :
    grid_name( given_grid_name )
  // this constructor uses the other constructor to make another square_grid
  // with the scaling factors, which is then used to construct this instance.
  {

    square_grid* scaling_grid = NULL;

    if( !(scaling_grid_file_location.empty()) )
      {

        scaling_grid = new square_grid( &scaling_grid_file_location,
                                        &given_grid_name,
                                        NULL );

      }

    read_in( &base_grid_file_location,
             scaling_grid );

    if( NULL != scaling_grid )
      {

        delete scaling_grid;

      }

  }

  square_grid::~square_grid()
  {

    // delete all the vectors held as pointers by the vector of vectors:

    for( std::vector< std::vector< double >* >::iterator
           deletion_iterator = values.begin();
         values.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;
        *deletion_iterator = NULL;  // just to be safe.

      }

  }

  void
  square_grid::read_in( std::string const* const grid_file_location,
                        square_grid const* const scaling_grid )
  // this actually does most of the job of the constructors.
  {

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: square_grid::read_in( "
    << *grid_file_location << ", " << scaling_grid << " ) called";
    std::cout << std::endl;**/

    if( NULL != grid_file_location )
      {

        // open the file:
        int const max_length_of_grid_file = 100000;
        CppSLHA::hash_commented_file_reader
        input_file_reader( grid_file_location,
                           max_length_of_grid_file,
                           false );
        /* a hundred thousand lines should be long enough for any grid file. if
         * read_in() is reading in more lines than that, it's probably because
         * of some file corruption.
         * also, if it's not enough, change it here.
         */

        std::string input_line_as_string;
        // prepare a string to parse the line.
        std::stringstream input_line_as_stream;
        // prepare a stringstream to parse the line.

        // the masses from the last line read & the currently-read line:
        double last_x_coordinate = 0.0;
        double last_y_coordinate = 0.0;
        double current_x_coordinate = 0.0;
        double current_y_coordinate = 0.0;

        double current_value;
        // the value currently being read in.

        bool grid_size_still_unknown = true;
        // we start out not knowing the step size for the grid (but we do
        // assume that it is a square grid).

        std::vector< double >* current_vector = NULL;
        // this holds all the values for a set of points with the same
        // x_coordinate.

        while( input_file_reader.read_line( &input_line_as_string ) )
          // the evaluation of the conditional reads in the next line.
          {

            // read in the line for interpretting:
            input_line_as_stream.clear();
            input_line_as_stream.str( input_line_as_string );

            // update what the masses for the last line were:
            last_y_coordinate = current_y_coordinate;
            last_x_coordinate = current_x_coordinate;

            input_line_as_stream >> current_x_coordinate
                                 >> current_y_coordinate
                                 >> current_value;

            if( current_x_coordinate > last_x_coordinate )
              // if the x_coordinate has changed, we need a new vector for the
              // values for varying y_coordinates.
              {

                // 1st record the last set of points with the same
                // x_coordinate:
                if( NULL != current_vector )
                  // if we have a vector to add...
                  {

                    values.push_back( current_vector );

                    // debugging:
                    /**std::cout
                    << std::endl
                    << "pushed back a vector of size "
                    << current_vector->size();**/

                  }

                // prepare a new vector for the new x_coordinate:
                current_vector = new std::vector< double >();

              }

            if( NULL != scaling_grid )
              // if we have to scale this value...
              {

                // debugging:
                /**std::cout
                << std::endl
                << "debugging: scaling point at ( "
                << current_x_coordinate << ", " << current_y_coordinate
                << " ) with value " << current_value << " by factor "
                << scaling_grid->value_at( current_x_coordinate,
                                           current_y_coordinate );
                std::cout << std::endl;**/

                current_value *= scaling_grid->value_at( current_x_coordinate,
                                                        current_y_coordinate );

              }

            current_vector->push_back( current_value );

            if( grid_size_still_unknown
                &&
                ( 0.0 < last_x_coordinate )
                &&
                ( 0.0 < last_y_coordinate ) )
              // if we still need to work out the grid step size & we have read
              // in enough to do so...
              {

                // first we take note of the 1st point, which should have the
                // lowest co-ordinate values:
                lowest_x_coordinate = last_x_coordinate;
                lowest_y_coordinate = last_y_coordinate;

                // then we work out the grid step size:
                if( current_x_coordinate > last_x_coordinate )
                  {

                    grid_step_size = ( current_x_coordinate
                                       - last_x_coordinate );
                    grid_size_still_unknown = false;

                  }
                else if( current_y_coordinate > last_y_coordinate )
                  {

                    grid_step_size = ( current_y_coordinate
                                       - last_y_coordinate );
                    grid_size_still_unknown = false;

                  }


                // debugging:
                /**std::cout
                << std::endl
                << "grid_step_size = " << grid_step_size;**/

              }

          }  // end of while loop going over the lines of the file.

        // we still have to add the last vector of values for the last x
        // co-ordinate:
        values.push_back( current_vector );

        // debugging:
        /**std::cout
        << std::endl
        << "pushed back last vector of size "
        << current_vector->size();
        std::cout << std::endl;**/

        // the last co-ordinates should be the maximal values:
        highest_x_coordinate = current_x_coordinate;
        highest_y_coordinate = current_y_coordinate;

      }
    else
      // otherwise, a NULL pointer was given instead of a pointer to the
      // grid file name string.
      {

        std::cout
        << std::endl
        << "LHC-FASER::error! asked to open a file with a NULL pointer!";
        std::cout << std::endl;

      }

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: square_grid::read_in() finished: "
    << std::endl << "lowest_x_coordinate = " << lowest_x_coordinate
    << std::endl << "highest_x_coordinate = " << highest_x_coordinate
    << std::endl << "lowest_y_coordinate = " << lowest_y_coordinate
    << std::endl << "highest_y_coordinate = " << highest_y_coordinate
    << std::endl << "value at ( " << lowest_x_coordinate << ", "
    << lowest_y_coordinate << " ) is " << value_at( lowest_x_coordinate,
                                                    lowest_y_coordinate )
    << std::endl << "value at ( "
    << ( 0.5 * ( lowest_x_coordinate + highest_x_coordinate ) + 0.2 ) << ", "
    << ( 0.5 * ( lowest_y_coordinate + highest_y_coordinate ) + 0.3 )
    << " ) is "
    << value_at( ( 0.5 * ( lowest_x_coordinate + highest_x_coordinate )
                   + 0.2 ),
              ( 0.5 * ( lowest_x_coordinate + highest_x_coordinate ) + 0.3 ) );
    std::cout << std::endl;**/

  }


  double
  square_grid::value_at( double const x_coordinate,
                         double const y_coordinate )
  const
  /* this finds the grid square which the given point is in, & then uses
   * LHC_FASER_global::square_bilinear_interpolation to get an interpolated
   * value. if the requested co-ordinates are outside the grid, or on the
   * edges of the grid corresponding to the maximum x or y values,
   * CppSLHA::CppSLHA_global::really_wrong_value is returned.
   * N.B.: cross-section grids use
   * x_coordinate == SQUARK mass, y_coordinate == GLUINO mass!
   */
  {

    if( ( 0.0 < grid_step_size )
        &&
        ( x_coordinate >= lowest_x_coordinate )
        &&
        ( y_coordinate >= lowest_y_coordinate ) )
      {

        double
        x_steps = ( ( x_coordinate - lowest_x_coordinate ) / grid_step_size );
        int lower_left_x = (int)x_steps;

        if( values.size() > ( lower_left_x + 1 ) )
          // if the x co-ordinate is less than its maximal grid value...
          {

            double
            y_steps
            = ( ( y_coordinate - lowest_y_coordinate ) / grid_step_size );
            int lower_left_y = (int)y_steps;

            if( values.at( lower_left_x )->size() > ( lower_left_y + 1 ) )
              // if the y co-ordinate is less than its maximal grid value...
              {

                return LHC_FASER_global::square_bilinear_interpolation(
                                                    ( x_steps - lower_left_x ),
                                                    ( y_steps - lower_left_y ),
                                 values.at( lower_left_x )->at( lower_left_y ),
                         values.at( ( lower_left_x + 1 ) )->at( lower_left_y ),
                 values.at( ( lower_left_x + 1 ) )->at( ( lower_left_y + 1 ) ),
                       values.at( lower_left_x )->at( ( lower_left_y + 1 ) ) );

              }
            else if( ( values.at( lower_left_x )->size()
                       == ( lower_left_y + 1 ) )
                     &&
                     ( (double)lower_left_y == y_steps ) )
              // otherwise, if it's on the maximal y edge...
              {

                return LHC_FASER_global::unit_linear_interpolation(
                                                    ( x_steps - lower_left_x ),
                                 values.at( lower_left_x )->at( lower_left_y ),
                       values.at( ( lower_left_x + 1 ) )->at( lower_left_y ) );

              }
            else
              // otherwise, it's off the grid:
              {

                return CppSLHA::CppSLHA_global::really_wrong_value;

              }

          }
        else if( ( values.size() == ( lower_left_x + 1 ) )
                 &&
                 ( (double)lower_left_x == x_steps ) )
          // otherwise, if it's on the maximal x edge...
          {

            double
            y_steps
            = ( ( y_coordinate - lowest_y_coordinate ) / grid_step_size );
            int lower_left_y = (int)y_steps;

            if( values.at( lower_left_x )->size() > ( lower_left_y + 1 ) )
              // if the y co-ordinate is less than its maximal grid value...
              {

                return LHC_FASER_global::unit_linear_interpolation(
                                                    ( y_steps - lower_left_y ),
                                 values.at( lower_left_x )->at( lower_left_y ),
                       values.at( lower_left_x )->at( ( lower_left_y + 1 ) ) );

              }
            else if( ( values.at( lower_left_x )->size()
                       == ( lower_left_y + 1 ) )
                     &&
                     ( (double)lower_left_y == y_steps ) )
              // otherwise, if it's on the maximal x & y corner...
              {

                return values.at( lower_left_x )->at( lower_left_y );

              }
            else
              // otherwise, it's off the grid:
              {

                return CppSLHA::CppSLHA_global::really_wrong_value;

              }

          }
        else
          // otherwise, it's off the grid:
          {

            return CppSLHA::CppSLHA_global::really_wrong_value;

          }

      }
    else
      // otherwise, it's off the grid (or the grid step size was not positive,
      // which is still a problem...):
      {

        return CppSLHA::CppSLHA_global::really_wrong_value;

      }

  }



  cross_section_table::cross_section_table(
                                   square_grid const* const given_lookup_table,
                         signed_particle_shortcut_pair const* const given_pair,
                                              double const given_flavor_factor,
                                   input_handler const* const given_shortcuts ) :
    readied_for_new_point( given_shortcuts->get_readier() ),
    shortcut( given_shortcuts ),
    scolored_pair( given_pair ),
    lookup_table( given_lookup_table ),
    flavor_factor( given_flavor_factor )
  {

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: cross_section_table( "
    << given_lookup_table << ", { "
    << given_pair->get_first_pointer()->get_name() << ": "
    << given_pair->is_first_particle() << "; "
    << given_pair->get_second_pointer()->get_name() << ": "
    << given_pair->is_second_particle() << " }, "
    << given_flavor_factor << ", " << given_shortcuts << " ) called";
    std::cout << std::endl;**/

    // we take the 1st squark found to ask for its mass in looking up the
    // cross-section. if no squark was provided, we keep the pointer as NULL.
    if( shortcut->is_in( given_pair->get_first_pointer()->get_PDG_code(),
                         shortcut->get_squarks() ) )
      {

        directly_produced_squark = given_pair->get_first_pointer();

      }
    else if( shortcut->is_in( given_pair->get_second_pointer()->get_PDG_code(),
                              shortcut->get_squarks() ) )
      {

        directly_produced_squark = given_pair->get_second_pointer();

      }
    else
      {

        directly_produced_squark = NULL;

      }

  }

  cross_section_table::~cross_section_table()
  {

    // does nothing.

  }


  double
  cross_section_table::get_value()
  /* this checks to see if it needs to recalculate the value for this point, &
   * if so, it does, & if the point lies outside the grid, it takes the value
   * of the nearest point on the grid instead.
   */
  {

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: cross_section_table::get_value() called."
    << " needs_to_prepare_for_this_point() = "
    << needs_to_prepare_for_this_point() << ", "
    << "flavor_factor = " << flavor_factor;
    std::cout << std::endl;**/

    if( needs_to_prepare_for_this_point() )
      {

        if( 0.0 < flavor_factor )
          {

            double squark_mass_to_use;

            if( NULL != directly_produced_squark )
              // if we have a directly-produced squark...
              {

                squark_mass_to_use
                = directly_produced_squark->get_absolute_mass();

              }
            else
              // otherwise, we use the average squark mass...
              {

                squark_mass_to_use = shortcut->get_average_squarks4_mass();

              }

            if( squark_mass_to_use
                < lookup_table->get_lowest_x() )
              {

                std::cout
                << std::endl
                << "LHC-FASER::warning! cross-section requested for a squark"
                << " mass ( " << squark_mass_to_use << " ) lower than the"
                << " lowest squark mass of the lookup table (";
                squark_mass_to_use = lookup_table->get_lowest_x();
                std::cout
                << squark_mass_to_use << "); using this lowest mass instead.";
                std::cout << std::endl;

              }
            else if( squark_mass_to_use
                     > lookup_table->get_highest_x() )
              {

                std::cout
                << std::endl
                << "LHC-FASER::warning! cross-section requested for a squark"
                << " mass ( " << squark_mass_to_use << " ) lower than the"
                << " highest squark mass of the lookup table ( ";
                squark_mass_to_use = lookup_table->get_highest_x();
                std::cout
                << squark_mass_to_use << " ); using this highest mass"
                << " instead.";
                std::cout << std::endl;

              }

            if( shortcut->get_gluino()->get_absolute_mass()
                < lookup_table->get_lowest_y() )
              {

                std::cout
                << std::endl
                << "LHC-FASER::warning! cross-section requested for a gluino"
                << " mass ( " << shortcut->get_gluino()->get_absolute_mass()
                << " ) lower than the lowest gluino mass of the lookup table"
                << " ( ";
                stored_value = lookup_table->value_at( squark_mass_to_use,
                                                lookup_table->get_lowest_y() );
                std::cout
                << lookup_table->get_lowest_y()
                << " ); using this lowest mass instead.";
                std::cout << std::endl;

              }
            else if( shortcut->get_gluino()->get_absolute_mass()
                     > lookup_table->get_highest_y() )
              {

                std::cout
                << std::endl
                << "LHC-FASER::warning! cross-section requested for a gluino"
                << " mass ( " << shortcut->get_gluino()->get_absolute_mass()
                << " ) lower than the highest gluino mass of the lookup table"
                << " ( ";
                stored_value = lookup_table->value_at( squark_mass_to_use,
                                               lookup_table->get_highest_y() );
                std::cout
                << lookup_table->get_highest_y()
                << " ); using this highest mass instead.";
                std::cout << std::endl;

              }
            else
              {
                stored_value = lookup_table->value_at( squark_mass_to_use,
                                 shortcut->get_gluino()->get_absolute_mass() );

              }

            stored_value *= flavor_factor;

          }
        else
          {

            stored_value = 0.0;

          }

        finish_preparing_for_this_point();

      }

    return stored_value;

  }



  cross_section_table_set::cross_section_table_set(
                                 std::string const* const given_grid_directory,
                                                    int const given_LHC_energy,
                                 input_handler const* const given_shortcuts ) :
    grid_directory( *given_grid_directory ),
    LHC_energy( given_LHC_energy ),
    shortcut( given_shortcuts )
  {

    // just an initialization list.

  }

  cross_section_table_set::~cross_section_table_set()
  {

    for( std::vector< square_grid* >::iterator
         deletion_iterator = grids.begin();
         grids.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

    for( std::vector< cross_section_table* >::iterator
         deletion_iterator = tables.begin();
         tables.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

  }


  double
  cross_section_table_set::prepare_grid_name( std::string* const grid_name,
                        signed_particle_shortcut_pair const* const given_pair )
  // this returns the flavor factor for the requested pair while putting the
  // square_grid's associated string in grid_name.
  {

    /* here we make the approximations about how the cross-sections are in
     * squark flavors.
     *
     * gluino+gluino is assumed to depend on the average squark mass.
     *
     * gluino+squark is assumed to be non-zero only for sdown & sup,
     * & that the ratio of sup to sdown is 2:1 (twice as many ups in the
     * other proton for the gluon from its proton to fuse with as downs).
     *
     * squark+antisquark is assumed to be non-zero only for both squarks
     * being the same flavor. all the _L+_L/_R+_R (_one+_one/_two+_two) entries
     * are assumed to be equal for the different flavors. for squark+antisquark
     * pairs of different chiralities, only sdown & sup have non-zero entries &
     * sdown_L+antisdown_R:sdown_L+antisup_R:sup_L+antisdown_R:sup_L+antisup_L
     * is assumed to be 1:1:2:2 & likewise for _L <-> _R (twice as many ups in
     * the other proton for the antiquark from its proton to fuse with as
     * downs).
     *
     * squark+squark is assumed to be non-zero only for both sdown, both
     * sup, or sdown+sup, & only for both being the same chirality (LO sQCD t-
     * & u-channels are identically zero if the squarks are opposite
     * chirality).
     * the ratios for sdown+sdown : sdown+sup : sup+sup are assumed to be
     * 1 * 8/3 * 1/2 : 4 * 1 * 1 : 4 * 8/3 * 1/2
     * (ratios of getting the various valence quarks from both
     * protons * matrix-element-squared in the heavy gluino mass limit
     * (t-channel + u-channel + interference terms with 1/N_c suppression, * 2
     * for t*u + u*t) * phase space halving for identical final-state
     * particles) = 4:12:16 = 1:3:4
     */

    if( ( CppSLHA::PDG_code::gluino
          == given_pair->get_first_pointer()->get_PDG_code() )
        &&
        ( CppSLHA::PDG_code::gluino
          == given_pair->get_second_pointer()->get_PDG_code() ) )
      {

        grid_name->assign( "gluino+gluino" );
        return 1.0;

      }
    else if( shortcut->is_in( given_pair->get_first_pointer()->get_PDG_code(),
                              shortcut->get_squarks() )
             &&
             ( CppSLHA::PDG_code::gluino
               == given_pair->get_second_pointer()->get_PDG_code() ) )
      {

        grid_name->assign( "sdown_L+gluino" );

        if( CppSLHA::PDG_code::sdown_L
            == given_pair->get_first_pointer()->get_PDG_code() )
          {

            return 1.0;

          }
        else if( CppSLHA::PDG_code::sdown_R
                 == given_pair->get_first_pointer()->get_PDG_code() )
          {

            return 1.0;

          }
        else if( CppSLHA::PDG_code::sup_L
                 == given_pair->get_first_pointer()->get_PDG_code() )
          {

            return 2.0;

          }
        else if( CppSLHA::PDG_code::sup_R
                 ==  given_pair->get_first_pointer()->get_PDG_code() )
          {

            return 2.0;

          }
        else
          {

            return 0.0;

          }

      }
    else if( shortcut->is_in( given_pair->get_second_pointer()->get_PDG_code(),
                              shortcut->get_squarks() )
             &&
             ( CppSLHA::PDG_code::gluino
               == given_pair->get_first_pointer()->get_PDG_code() ) )
      {

        grid_name->assign( "sdown_L+gluino" );

        if( CppSLHA::PDG_code::sdown_L
            == given_pair->get_second_pointer()->get_PDG_code() )
          {

            return 1.0;

          }
        else if( CppSLHA::PDG_code::sdown_R
                 == given_pair->get_second_pointer()->get_PDG_code() )
          {

            return 1.0;

          }
        else if( CppSLHA::PDG_code::sup_L
                 == given_pair->get_second_pointer()->get_PDG_code() )
          {

            return 2.0;

          }
        else if( CppSLHA::PDG_code::sup_R
                 ==  given_pair->get_second_pointer()->get_PDG_code() )
          {

            return 2.0;

          }
        else
          {

            return 0.0;

          }

      }
    else if( shortcut->is_in( given_pair->get_first_pointer()->get_PDG_code(),
                              shortcut->get_squarks() )
             &&
             shortcut->is_in( given_pair->get_second_pointer()->get_PDG_code(),
                              shortcut->get_squarks() ) )
      {

        if( ( given_pair->first_is_not_antiparticle()
              &&
              !given_pair->second_is_not_antiparticle() )
            ||
            ( !given_pair->first_is_not_antiparticle()
              &&
              given_pair->second_is_not_antiparticle() ) )
          {

            if( ( CppSLHA::PDG_code::sdown_L
                == given_pair->get_first_pointer()->get_PDG_code() )
                &&
                ( CppSLHA::PDG_code::sdown_L
                == given_pair->get_second_pointer()->get_PDG_code() ) )
              {

                grid_name->assign( "sdown_L+antisdown_L" );
                return 1.0;

              }
            else if( ( CppSLHA::PDG_code::sdown_R
                       == given_pair->get_first_pointer()->get_PDG_code() )
                      &&
                      ( CppSLHA::PDG_code::sdown_R
                        == given_pair->get_second_pointer()->get_PDG_code() ) )
              {

                grid_name->assign( "sdown_L+antisdown_L" );
                return 1.0;

              }
            else if( ( CppSLHA::PDG_code::sdown_L
                       == given_pair->get_first_pointer()->get_PDG_code() )
                      &&
                      ( CppSLHA::PDG_code::sdown_R
                        == given_pair->get_second_pointer()->get_PDG_code() ) )
              {

                grid_name->assign( "sdown_L+antisdown_R" );
                return 1.0;

              }
            else if( ( CppSLHA::PDG_code::sdown_R
                       == given_pair->get_first_pointer()->get_PDG_code() )
                      &&
                      ( CppSLHA::PDG_code::sdown_L
                        == given_pair->get_second_pointer()->get_PDG_code() ) )
              {

                grid_name->assign( "sdown_L+antisdown_R" );
                return 1.0;

              }
            else if( ( CppSLHA::PDG_code::sdown_L
                       == given_pair->get_first_pointer()->get_PDG_code() )
                      &&
                      ( CppSLHA::PDG_code::sup_R
                        == given_pair->get_second_pointer()->get_PDG_code() ) )
              {

                grid_name->assign( "sdown_L+antisdown_R" );
                return 1.0;

              }
            else if( ( CppSLHA::PDG_code::sdown_R
                       == given_pair->get_first_pointer()->get_PDG_code() )
                      &&
                      ( CppSLHA::PDG_code::sup_L
                        == given_pair->get_second_pointer()->get_PDG_code() ) )
              {

                grid_name->assign( "sdown_L+antisdown_R" );
                return 1.0;

              }
            else if( ( CppSLHA::PDG_code::sup_L
                       == given_pair->get_first_pointer()->get_PDG_code() )
                      &&
                      ( CppSLHA::PDG_code::sdown_R
                        == given_pair->get_second_pointer()->get_PDG_code() ) )
              {

                grid_name->assign( "sdown_L+antisdown_R" );
                return 2.0;

              }
            else if( ( CppSLHA::PDG_code::sup_R
                       == given_pair->get_first_pointer()->get_PDG_code() )
                      &&
                      ( CppSLHA::PDG_code::sdown_L
                        == given_pair->get_second_pointer()->get_PDG_code() ) )
              {

                grid_name->assign( "sdown_L+antisdown_R" );
                return 2.0;

              }
            else if( ( CppSLHA::PDG_code::sup_L
                       == given_pair->get_first_pointer()->get_PDG_code() )
                      &&
                      ( CppSLHA::PDG_code::sup_L
                        == given_pair->get_second_pointer()->get_PDG_code() ) )
              {

                grid_name->assign( "sdown_L+antisdown_L" );
                return 1.0;

              }
            else if( ( CppSLHA::PDG_code::sup_R
                       == given_pair->get_first_pointer()->get_PDG_code() )
                      &&
                      ( CppSLHA::PDG_code::sup_R
                        == given_pair->get_second_pointer()->get_PDG_code() ) )
              {

                grid_name->assign( "sdown_L+antisdown_L" );
                return 1.0;

              }
            else if( ( CppSLHA::PDG_code::sup_L
                       == given_pair->get_first_pointer()->get_PDG_code() )
                      &&
                      ( CppSLHA::PDG_code::sup_R
                        == given_pair->get_second_pointer()->get_PDG_code() ) )
              {

                grid_name->assign( "sdown_L+antisdown_R" );
                return 2.0;

              }
            else if( ( CppSLHA::PDG_code::sup_R
                       == given_pair->get_first_pointer()->get_PDG_code() )
                      &&
                      ( CppSLHA::PDG_code::sup_L
                        == given_pair->get_second_pointer()->get_PDG_code() ) )
              {

                grid_name->assign( "sdown_L+antisdown_R" );
                return 2.0;

              }
            else if( ( CppSLHA::PDG_code::sstrange_L
                       == given_pair->get_first_pointer()->get_PDG_code() )
                      &&
                      ( CppSLHA::PDG_code::sstrange_L
                        == given_pair->get_second_pointer()->get_PDG_code() ) )
              {

                grid_name->assign( "sdown_L+antisdown_L" );
                return 1.0;

              }
            else if( ( CppSLHA::PDG_code::sstrange_R
                       == given_pair->get_first_pointer()->get_PDG_code() )
                      &&
                      ( CppSLHA::PDG_code::sstrange_R
                        == given_pair->get_second_pointer()->get_PDG_code() ) )
              {

                grid_name->assign( "sdown_L+antisdown_L" );
                return 1.0;

              }
            else if( ( CppSLHA::PDG_code::scharm_L
                       == given_pair->get_first_pointer()->get_PDG_code() )
                      &&
                      ( CppSLHA::PDG_code::scharm_L
                        == given_pair->get_second_pointer()->get_PDG_code() ) )
              {

                grid_name->assign( "sdown_L+antisdown_L" );
                return 1.0;

              }
            else if( ( CppSLHA::PDG_code::scharm_R
                       == given_pair->get_first_pointer()->get_PDG_code() )
                      &&
                      ( CppSLHA::PDG_code::scharm_R
                        == given_pair->get_second_pointer()->get_PDG_code() ) )
              {

                grid_name->assign( "sdown_L+antisdown_L" );
                return 1.0;

              }
            else if( ( CppSLHA::PDG_code::sbottom_one
                       == given_pair->get_first_pointer()->get_PDG_code() )
                      &&
                      ( CppSLHA::PDG_code::sbottom_one
                        == given_pair->get_second_pointer()->get_PDG_code() ) )
              {

                grid_name->assign( "sdown_L+antisdown_L" );
                return 1.0;

              }
            else if( ( CppSLHA::PDG_code::sbottom_two
                       == given_pair->get_first_pointer()->get_PDG_code() )
                      &&
                      ( CppSLHA::PDG_code::sbottom_two
                        == given_pair->get_second_pointer()->get_PDG_code() ) )
              {

                grid_name->assign( "sdown_L+antisdown_L" );
                return 1.0;

              }
            else if( ( CppSLHA::PDG_code::stop_one
                       == given_pair->get_first_pointer()->get_PDG_code() )
                      &&
                      ( CppSLHA::PDG_code::stop_one
                        == given_pair->get_second_pointer()->get_PDG_code() ) )
              {

                grid_name->assign( "sdown_L+antisdown_L" );
                return 1.0;

              }
            else if( ( CppSLHA::PDG_code::stop_two
                       == given_pair->get_first_pointer()->get_PDG_code() )
                      &&
                      ( CppSLHA::PDG_code::stop_two
                        == given_pair->get_second_pointer()->get_PDG_code() ) )
              {

                grid_name->assign( "sdown_L+antisdown_L" );
                return 1.0;

              }
            else
              {

                return 0.0;

              }

          }
        else if( given_pair->first_is_not_antiparticle()
                 &&
                 given_pair->second_is_not_antiparticle() )
          {

            if( ( CppSLHA::PDG_code::sdown_L
                  == given_pair->get_first_pointer()->get_PDG_code() )
                &&
                ( CppSLHA::PDG_code::sdown_L
                  == given_pair->get_second_pointer()->get_PDG_code() ) )
              {

                grid_name->assign( "sdown_L+sdown_L" );
                return 1.0;

              }
            else if( ( CppSLHA::PDG_code::sdown_R
                        == given_pair->get_first_pointer()->get_PDG_code() )
                     &&
                     ( CppSLHA::PDG_code::sdown_R
                       == given_pair->get_second_pointer()->get_PDG_code() ) )
              {

                grid_name->assign( "sdown_L+sdown_L" );
                return 1.0;

              }
            else if( ( CppSLHA::PDG_code::sdown_L
                       == given_pair->get_first_pointer()->get_PDG_code() )
                     &&
                     ( CppSLHA::PDG_code::sup_L
                       == given_pair->get_second_pointer()->get_PDG_code() ) )
              {

                grid_name->assign( "sdown_L+sdown_L" );
                return 3.0;

              }
            else if( ( CppSLHA::PDG_code::sdown_R
                       == given_pair->get_first_pointer()->get_PDG_code() )
                     &&
                     ( CppSLHA::PDG_code::sup_R
                       == given_pair->get_second_pointer()->get_PDG_code() ) )
              {

                grid_name->assign( "sdown_L+sdown_L" );
                return 3.0;

              }
            else if( ( CppSLHA::PDG_code::sup_L
                       == given_pair->get_first_pointer()->get_PDG_code() )
                     &&
                     ( CppSLHA::PDG_code::sdown_L
                       == given_pair->get_second_pointer()->get_PDG_code() ) )
              {

                grid_name->assign( "sdown_L+sdown_L" );
                return 3.0;

              }
            else if( ( CppSLHA::PDG_code::sup_R
                       == given_pair->get_first_pointer()->get_PDG_code() )
                     &&
                     ( CppSLHA::PDG_code::sdown_R
                       == given_pair->get_second_pointer()->get_PDG_code() ) )
              {

                grid_name->assign( "sdown_L+sdown_L" );
                return 3.0;

              }
            else if( ( CppSLHA::PDG_code::sup_L
                       == given_pair->get_first_pointer()->get_PDG_code() )
                     &&
                     ( CppSLHA::PDG_code::sup_L
                       == given_pair->get_second_pointer()->get_PDG_code() ) )
              {

                grid_name->assign( "sdown_L+sdown_L" );
                return 4.0;

              }
            else if( ( CppSLHA::PDG_code::sup_R
                       == given_pair->get_first_pointer()->get_PDG_code() )
                     &&
                     ( CppSLHA::PDG_code::sup_R
                       == given_pair->get_second_pointer()->get_PDG_code() ) )
              {

                grid_name->assign( "sdown_L+sdown_L" );
                return 4.0;

              }
            else
              {

                return 0.0;

              }

          }
        else
          {

            return 0.0;

          }

      }
    else
      {

        return 0.0;

      }

  }


  square_grid const*
  cross_section_table_set::get_grid( std::string const* const grid_name )
  // this either finds the appropriate square_grid in grids or makes a new
  // instance, stores it, & returns its pointer.
  {

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: cross_section_table_set::get_grid( "
    << *grid_name << " ) called";
    std::cout << std::endl;**/

    square_grid* return_pointer = NULL;
    // this starts as NULL so that we know if it wasn't found among the
    // existing instances.

    for( std::vector< square_grid* >::iterator
         grid_iterator = grids.begin();
         grids.end() > grid_iterator;
         ++grid_iterator )
      {

        if( 0 == (*grid_iterator)->get_name()->compare( *grid_name ) )
          // if a table with this name already exists...
          {

            return *grid_iterator;
            // note the table.
            grid_iterator = grids.end();
            // stop looking.

          }

      }

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: return_pointer = " << return_pointer;
    std::cout << std::endl;**/

    if( NULL == return_pointer )
      // if we didn't find an existing match, we make a new instance, store
      // it & return a pointer to it:
      {

        std::string LO_grid_file_name( grid_directory );
        LO_grid_file_name.append( "/" );
        LO_grid_file_name.append( *grid_name );
        LO_grid_file_name.append( "_LO.dat" );

        // debugging:
        /**std::cout
        << std::endl
        << "debugging: LO_grid_file_name = " << LO_grid_file_name;
        std::cout << std::endl;**/

        // we check for whether we should use NLO K-factors:
        square_grid* NLO_grid = NULL;
        if( shortcut->using_NLO() )
          {

            std::string NLO_grid_file_name( grid_directory );
            NLO_grid_file_name.append( "/" );
            NLO_grid_file_name.append( *grid_name );
            NLO_grid_file_name.append( "_K.dat" );

            // debugging:
            /**std::cout
            << std::endl
            << "debugging: NLO_grid_file_name = " << NLO_grid_file_name;
            std::cout << std::endl;**/

            NLO_grid = new square_grid( &NLO_grid_file_name,
                                        grid_name,
                                        NULL );
            // the grid's name isn't important because it is only temporary.

            // debugging:
            /**std::cout
            << std::endl
            << "debugging: NLO_grid = " << NLO_grid
            << ", value_at( 200.0, 200.0 ) = "
            << NLO_grid->value_at( 200.0,
                                   200.0 );
            std::cout << std::endl;**/

          }

        // we make a new table:
        return_pointer = new square_grid( &LO_grid_file_name,
                                          grid_name,
                                          NLO_grid );

        if( NULL != NLO_grid )
          // if we needed a K-factor grid...
          {

            // we can let it go now:
            delete NLO_grid;

          }

        // debugging:
        /**std::cout
        << std::endl
        << "debugging: return_pointer = " << return_pointer;
        std::cout << std::endl;**/

        grids.push_back( return_pointer );

      }

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: returning = " << return_pointer;
    std::cout << std::endl;**/

    return return_pointer;

  }


  cross_section_table*
  cross_section_table_set::get_table(
                 signed_particle_shortcut_pair const* const requested_channel )
  // this returns the cross_section_table for the requested pair.
  {

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: cross_section_table_set::get_table( { "
    << *(requested_channel->get_first_pointer()->get_name_or_antiname(
                              requested_channel->first_is_not_antiparticle() ))
    << "; "
    << *(requested_channel->get_second_pointer()->get_name_or_antiname(
                             requested_channel->second_is_not_antiparticle() ))
    << " } ) called";
    std::cout << std::endl;**/

    cross_section_table* return_pointer = NULL;
    // this starts as NULL so that we know if it wasn't found among the
    // existing instances.

    for( std::vector< cross_section_table* >::iterator
         table_iterator = tables.begin();
         tables.end() > table_iterator;
         ++table_iterator )
      {

        if( requested_channel == (*table_iterator)->get_pair() )
          // if we find the requested pair...
          {

            // debugging:
            /**std::cout
            << std::endl
            << "found existing table.";
            std::cout << std::endl;**/

            return_pointer = *table_iterator;
            // note the table to return.
            table_iterator = tables.end();
            // stop looking.

          }

      }


    // debugging:
    /**std::cout
    << std::endl
    << "debugging: return_pointer = " << return_pointer;
    std::cout << std::endl;**/

    if( NULL == return_pointer )
      // if we didn't find an existing match, we make a new instance, store
      // it & return a pointer to it:
      {

        // debugging:
        /**std::cout
        << std::endl
        << "did not find existing table.";
        std::cout << std::endl;**/

        // 1st we need to find out which square_grid to use & what the flavor
        // factor is:
        std::string grid_name;
        double flavor_factor = prepare_grid_name( &grid_name,
                                                  requested_channel );

        // debugging:
        /**std::cout
        << std::endl
        << "debugging: grid_path = " << grid_name
        << std::endl
        << "debugging: flavor_factor = " << flavor_factor;
        std::cout << std::endl;**/
        if( flavor_factor > 0.0 )
          // if it is a non-zero channel as far as this code is concerned...
          {

            return_pointer = new cross_section_table( get_grid( &grid_name ),
                                                      requested_channel,
                                                      flavor_factor,
                                                      shortcut );

          }
        else
          {

            return_pointer = new cross_section_table( NULL,
                                                      requested_channel,
                                                      0.0,
                                                      shortcut );
            /* a cross_section_table initialized with a NULL square_grid
             * pointer *should* be OK if its flavor factor is 0.0, so that it
             * never calls any of the square_grid's member functions.
             */

          }

        tables.push_back( return_pointer );

      }

    return return_pointer;

  }



  cross_section_handler::cross_section_handler(
                                 input_handler const* const given_shortcuts ) :
    shortcut( given_shortcuts )
  {

    // just an initialization list

  }

  cross_section_handler::~cross_section_handler()
  {

    for( std::vector< cross_section_table_set* >::iterator
         deletion_iterator = table_sets.begin();
         table_sets.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

  }


  cross_section_table_set*
  cross_section_handler::get_table_set( int const LHC_energy_in_TeV )
  {

    cross_section_table_set* return_pointer = NULL;
    // this starts as NULL so that we know if it wasn't found among the
    // existing instances.

    for( std::vector< cross_section_table_set* >::iterator
         set_iterator = table_sets.begin();
         table_sets.end() > set_iterator;
         ++set_iterator )
      {

        if( LHC_energy_in_TeV == (*set_iterator)->get_LHC_energy() )
          {

            return_pointer = *set_iterator;
            set_iterator = table_sets.end();
            // stop looking.

          }

      }

    if( NULL == return_pointer )
      // if we didn't find an existing match, we make a new instance, store
      // it & return a pointer to it:
      {

        std::string
        energy_name( *(shortcut->inspect_path_to_cross_section_grids()) );
        std::stringstream energy_stream( "" );
        energy_stream << "/" << LHC_energy_in_TeV << "TeV";
        energy_name.append( energy_stream.str() );

        return_pointer = new cross_section_table_set( &energy_name,
                                                      LHC_energy_in_TeV,
                                                      shortcut );
        table_sets.push_back( return_pointer );

      }

    return return_pointer;

  }

}  // end of LHC_FASER namespace.
