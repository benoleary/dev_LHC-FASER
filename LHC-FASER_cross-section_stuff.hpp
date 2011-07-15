/*
 * LHC-FASER_cross-section_stuff.hpp
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

#ifndef LHC_FASER_CROSS_SECTION_STUFF_HPP_
#define LHC_FASER_CROSS_SECTION_STUFF_HPP_

#include "LHC-FASER_global_stuff.hpp"
#include "LHC-FASER_input_handling_stuff.hpp"

namespace LHC_FASER
{

  /* this class reads in a file in the assumed format, stores it, & gives out
   * interpolated values. it was written with cross-sections in mind, with data
   * files in the format
   * squark_mass gluino_mass cross-section newline
   * or K-factors, in a similar format.
   */
  class square_grid
  {

  public:

    square_grid( std::string const* const grid_file_location,
                 std::string const* const given_grid_name,
                 square_grid const* const scaling_grid )
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
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    square_grid( std::string base_grid_file_location,
                 std::string given_grid_name,
                 std::string scaling_grid_file_location )
    // this constructor uses the other constructor to make another square_grid
    // with the scaling factors, which is then used to construct this instance.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~square_grid()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    std::string const*
    get_name()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    value_at( double const x_coordinate,
              double const y_coordinate )
    const
    /* this finds the grid square which the given point is in, & then uses
     * LHC_FASER_global::square_bilinear_interpolation to get an interpolated
     * value.
     * N.B.: cross-section grids use
     * x_coordinate == SQUARK mass, y_coordinate == GLUINO mass!
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_lowest_x()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_highest_x()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_lowest_y()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_highest_y()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    std::string grid_name;
    double grid_step_size;
    double lowest_x_coordinate;
    double highest_x_coordinate;
    double lowest_y_coordinate;
    double highest_y_coordinate;
    std::vector< std::vector< double >* > values;

    void
    read_in( std::string const* const grid_file_location,
             square_grid const* const scaling_grid )
    // this actually does most of the job of the constructors.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };


  /* this class holds a square_grid instance as a cross-section lookup table
   * for a particular colored sparticle combination, & returns the
   * interpolated value by looking up the masses of its sparticles.
   * this class is small enough to be kept entirely in this header file.
   * N.B.: currently, this class assumes that only mass-degenerate
   * squark-squark or squark-antisquark pairs are produced. this should be fine
   * at the accuracy we work at unless there are large splittings among the
   * sdown & sup mass eigenstates.
   */
  class cross_section_table : public readied_for_new_point
  {

  public:

    cross_section_table( square_grid const* const given_lookup_table,
                         signed_particle_shortcut_pair const* const given_pair,
                         double const given_flavor_factor,
                         input_handler const* const given_shortcuts )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~cross_section_table()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    signed_particle_shortcut_pair const*
    get_pair()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    double
    get_value()
    /* this checks to see if it needs to recalculate the value for this point,
     * & if so, it does, & if the point lies outside the grid, it takes the
     * value of the nearest point on the grid instead.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    signed_particle_shortcut_pair const* scolored_pair;
    CppSLHA::particle_property_set const* directly_produced_squark;
    input_handler const* const shortcut;

    double stored_value;

    square_grid const* const lookup_table;
    double const flavor_factor;
    /* this is used when the flavor of squark for this table uses the grid for
     * a different flavor of squark, in the approximation that it has the same
     * cross-section scaled by some factor (e.g. the gluino+sup cross-section
     * is twice the gluino+sdown cross-section for the same squark mass).
     */

  };


  /* this class holds several square_grid instances as cross-section lookup
   * tables for various colored sparticle combinations, & returns their
   * interpolated values for given signed_particle_shortcut_pair pointers.
   */
  class cross_section_table_set
  {

  public:

    cross_section_table_set( std::string const* const given_grid_directory,
                             int const given_LHC_energy,
                             input_handler const* const given_shortcuts )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~cross_section_table_set()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    int
    get_LHC_energy()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    cross_section_table*
    get_table( signed_particle_shortcut_pair const* const requested_channel )
    // this returns the cross_section_table for the requested pair.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    input_handler const* const shortcut;

    std::string const grid_directory;
    int const LHC_energy;

    std::vector< square_grid* > grids;
    // this holds square_grid instances with unique names.
    std::vector< cross_section_table* > tables;
    // this holds cross_section_table instances with appropriate
    // signed_particle_shortcut_pair pointers.

    double
    prepare_grid_name( std::string* const grid_name,
                       signed_particle_shortcut_pair const* const given_pair )
    // this returns the flavor factor for the requested pair while putting the
    // square_grid's associated string in grid_name.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    square_grid const*
    get_grid( std::string const* const grid_name )
    // this either finds the appropriate square_grid in grids or makes a new
    // instance, stores it, & returns its pointer.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };


  class cross_section_handler
  {

  public:

    cross_section_handler( input_handler const* const given_shortcuts )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~cross_section_handler()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    cross_section_table_set*
    get_table_set( int const LHC_energy_in_TeV )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    cross_section_table*
    get_table( int const LHC_energy_in_TeV,
               signed_particle_shortcut_pair const* const requested_channel )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    input_handler const* const shortcut;
    std::vector< cross_section_table_set* > table_sets;
    // this is for holding all the cross-sections which are looked up by the
    // signal handlers.

  };




  // inline functions:


  inline std::string const*
  square_grid::get_name()
  const
  {

    return &grid_name;

  }

  inline double
  square_grid::get_lowest_x()
  const
  {

    return lowest_x_coordinate;

  }

  inline double
  square_grid::get_highest_x()
  const
  {

    return highest_x_coordinate;

  }

  inline double
  square_grid::get_lowest_y()
  const
  {

    return lowest_y_coordinate;

  }

  inline double
  square_grid::get_highest_y()
  const
  {

    return highest_y_coordinate;

  }



  inline signed_particle_shortcut_pair const*
  cross_section_table::get_pair()
  const
  {

    return scolored_pair;

  }



  inline int
  cross_section_table_set::get_LHC_energy()
  const
  {

    return LHC_energy;

  }



  inline cross_section_table*
  cross_section_handler::get_table( int const LHC_energy_in_TeV,
                 signed_particle_shortcut_pair const* const requested_channel )
  {

    return get_table_set( LHC_energy_in_TeV )->get_table( requested_channel );

  }

}  // end of LHC_FASER namespace.

#endif /* LHC_FASER_CROSS_SECTION_STUFF_HPP_ */
