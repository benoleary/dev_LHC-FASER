/*
 * CppSLHA_EW_scale_spectrum.hpp
 *
 *  Created on: 12 Sep 2009
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2009 Ben O'Leary
 *
 *      This file is part of CppSLHA.
 *
 *      CppSLHA is free software: you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation, either version 3 of the License, or
 *      (at your option) any later version.
 *
 *      CppSLHA is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *
 *      You should have received a copy of the GNU General Public License
 *      along with CppSLHA.  If not, see <http://www.gnu.org/licenses/>.
 *
 *      The GNU General Public License should be in GNU_public_license.txt
 *      the files of CppSLHA are:
 *      CppSLHA.hpp
 *      CppSLHA.cpp
 *      CppSLHA_BLOCK.hpp
 *      CppSLHA_BLOCK.cpp
 *      CppSLHA_EW_scale_spectrum.hpp
 *      CppSLHA_EW_scale_spectrum.cpp
 *      CppSLHA_file_reading_stuff.hpp
 *      CppSLHA_file_reading_stuff.cpp
 *      CppSLHA_global_stuff.hpp
 *      CppSLHA_global_stuff.cpp
 *      CppSLHA_PDG_codes_and_data.hpp
 *      CppSLHA_PDG_codes_and_data.cpp
 *      CppSLHA_particle_data.hpp
 *      CppSLHA_particle_data.cpp
 *      CppSLHA_waiting_on_subprocess_executor.hpp
 *      CppSLHA_waiting_on_subprocess_executor.cpp
 *      and ./README.CppSLHA.txt which describes the package.
 *
 */

/* this is the header file for a set of C++ classes that are basically just
 * the set of masses of a particle (+ sparticle) spectrum, stored as a
 * variable-sized vector array of particle_property_set instances
 * (see particle_data.hpp).
 * there is an internal ordering, but data are sought by means of the PDG
 * particle codes for particles.
 * (there is a possibility that I will implement lookup-by-name, but that's
 * probably too much work to be worth it.)
 *
 * the PDG particle codes are defined in the CppSLHA_particle_codes.hpp file.
 */

/*********************************************************************
 *                                                                   *
 * IMPORTANT!                                                        *
 *                                                                   *
 * welcome, person who has decided to read this file.                *
 *                                                                   *
 * I hope that you will find these classes easy to use. I have       *
 * tried to make the variable names self-explanatory, & hope that    *
 * the comments are adequate to explain the logic of the code.       *
 *                                                                   *
 * all dimensionful quantities in the code are assumed to be in GeV. *
 *                                                                   *
 * following some stuff which I read in a book, I have tried to      *
 * stick to object-oriented principles, & thus tried to make my      *
 * classes closed to modification, but open to extension. hence, if  *
 * you are intending to modify this code to accommodate a spectrum   *
 * for a different extension of the Standard Model, please have a    *
 * look at how the EW_scale_MSSM_spectrum class extends the          *
 * EW_scale_SM_spectrum class, & how the                             *
 * EW_scale_NMSSM_spectrum class extends the                         *
 * EW_scale_MSSM_spectrum class. these are found after the base      *
 * EW_scale_spectrum class, so from line 1110 of this file.          *
 * PLEASE NOTE that I have just used the PDG particle codes as they  *
 * are given in the SLHA2 paper:                                     *
 * B. C. Allanach et al., Comp.Phys.Commun.180:8-25,2009,            *
 * http://arxiv.org/abs/0801.0045 online.                            *
 * there the authors have given the 2 extra scalar boson states      *
 * numbers 45 & 46 (neutral scalar & pseudoscalar), & the 5th        *
 * neutralino number 1000045.                                        *
 * please also note that there are separate classes for the          *
 * flavor violation sparticle numbering scheme, called               *
 * EW_scale_FVMSSM_spectrum & EW_scale_FVNMSSM_spectrum which have   *
 * the MSSM & NMSSM particle content respectively, but relabelled    *
 * sfermions.                                                        *
 *                                                                   *
 * also, a note on antiparticles: the EW_scale_spectrum class only   *
 * holds data for the particles with positive PDG particle codes.    *
 * the antiparticles, which have the negative versions of the PDG    *
 * particle codes for the particles, are assumed to have exactly the *
 * same mass and decay width, and to decay into exactly the same     *
 * decays but with particles replaced by their antiparticles and     *
 * antiparticles replaced with their particles, as their particles.  *
 * accessing a property set with the antiparticle PDG particle code  *
 * will result in receiving the property set of the particle with    *
 * positive PDG particle code. it is up to the user to deal with any *
 * antiparticle-specific stuff. e.g., maybe going down a sbottom     *
 * cascade decay chain where it decays to a negative chargino        *
 * (-1000024) and then looking for the chargino decays to a          *
 * sneutrino and the lepton with the correct charge. particles which *
 * are also their own antiparticles have                             *
 * counts_as_self_conjugate_flag as true, and all the functions      *
 * should correctly relabel negative PDG particle codes for these    *
 * particles back to the positive versions.                          *
 *                                                                   *
 *********************************************************************/

/*********************************************************************
 *                                                                   *
 * important!                                                        *
 *                                                                   *
 * the internal ordering of the particles IS NOT THE ORDER IN WHICH  *
 * THE PARTICLES APPEAR in CppSLHA_PDG_codes_and_data.fpp for the    *
 * particle codes. the ordering is defined by the constructor and    *
 * doesn't actually matter, unless one has modified the code to look *
 * in the std::vector arrays directly, which is BAD PRACTICE!        *
 * all these classes were built around the idea of looking up by     *
 * PDG particle code.                                                *
 *                                                                   *
 *********************************************************************/


#ifndef CPPSLHA_EW_SCALE_SPECTRUM_HPP_
#define CPPSLHA_EW_SCALE_SPECTRUM_HPP_
// this should prevent this class being included multiple times.


#include "CppSLHA_global_stuff.hpp"
#include "CppSLHA_PDG_codes_and_data.hpp"
#include "CppSLHA_particle_data.hpp"
#include "CppSLHA_waiting_on_subprocess_executor.hpp"


namespace CppSLHA
{

  // this class holds a string (the label to draw) & a double (the height of
  // the label).
  class drawing_label
  {

  public:

    static bool
    compare_pointer_heights( drawing_label const* const first_label,
                             drawing_label const* const second_label )
    // this returns true if the pointers are ordered so that their heights are
    // in ascending order or equal, false otherwise.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    drawing_label( std::string const* const given_label,
                   double const given_height )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~drawing_label()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    std::string const*
    get_label()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    set_label( std::string const* const given_label )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_label_height()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    set_label_height( double const given_height )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    std::string label_string;
    double label_height;

  };  // end of drawing_label class.

  // this class holds a string (the color to draw) & a pair of doubles (the
  // height of the line & the height of the label to draw to).
  class drawing_line
  {

  public:

    drawing_line( std::string const* const given_color,
                  double const given_line_height,
                  drawing_label const* const given_label )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~drawing_line()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    std::string const*
    get_color()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_line_height()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_label_height()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    std::string color_string;
    double line_height;
    drawing_label const* line_label;

  };  // end of drawing_line class.

  // this class holds a set of particles which are to be drawn together in a
  // column, with their color & a few helper numbers.
  class drawing_column
  {

  public:

    drawing_column()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~drawing_column()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    static bool
    compare_particles_list_pointers(
                                 std::pair< particle_property_set const* const,
                                          std::string* > const* first_particle,
                                 std::pair< particle_property_set const* const,
                                        std::string* > const* second_particle )
    // this returns true if the pointers are ordered so that their absolute
    // heights are in ascending order or equal, false otherwise.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    bool
    is_empty()
    const
    // this returns true if there are no lines to draw, false otherwise.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    add_particle( particle_property_set const* const given_particle,
                  std::string const* const given_color )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    std::list< drawing_line* > const*
    get_lines()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    std::list< drawing_label* > const*
    get_labels()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    std::list< std::pair< particle_property_set const* const,
                          std::string* >* >*
    get_particles()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    double
    get_highest_height()
    const
    // this goes through all the particles stored in the column & finds the
    // highest height for drawing a particle's line.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    prepare_lines_and_labels( double const vertical_scale )
    /* this fills the lists of lines & labels from the list of particles,
     * ensuring that the labels are in the right order & not overlapping.
     * it needs to know the overall vertical scale though.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    std::list< std::pair< particle_property_set const* const,
                          std::string* >* > particles;
    std::list< drawing_line* > lines;
    std::list< drawing_label* > labels;

    double const label_merging_distance;
    /* labels for particles when drawn get put side-by-side instead of atop
     * each other if they are within label_merging_distance of the vertical
     * scale from each other.
     */
    double const label_clashing_distance;
    // labels for particles are kept at least label_clashing_distance of the
    // vertical scale away from each other.
    double const label_floating_distance;
    /* labels for particles are moved to label_floating_distance of the
     * vertical scale away from each other if they are within
     * label_clashing_distance of the vertical scale away from each other.
     */

  };  // end of drawing_column class.

  // this class organizes the columns & draws the spectrum.
  class spectrum_drawer
  {

  public:

    spectrum_drawer()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~spectrum_drawer()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    reset()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    double
    get_highest_height()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    std::vector< drawing_column* > const*
    get_column_set()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    add_particle( particle_property_set const* const particle_to_be_added,
                  int const column_to_draw_in,
                  std::string const* const color_to_draw_with )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    double
    find_highest_height()
    // this goes through all the particles stored in all the columns & finds
    // the highest height for drawing a particle's line.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    draw_spectrum( std::string const* const output_eps_filename,
                   std::string const* const path_to_gnuplot,
                   std::string const* const path_to_latex,
                   std::string const* const path_to_dvips,
                   std::string const* const path_to_ps2eps )
    /* this actually draws the spectrum using gnuplot, first setting up the
     * lines & labels. the "path_to_X" strings should include the executable
     * name, e.g. /usr/bin/gnuplot for path_to_gnuplot.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    std::vector< drawing_column* > column_set;
    double highest_height;

    double const vertical_scale_factor;
    // the vertical scale is vertical_scale_factor * the largest particle mass
    // to be drawn.
    double const line_gap;
    // how much of a gap to leave between mass line column pairs.
    double const line_length;
    // how long each mass line is.
    double const line_offset;
    // how much every 2nd column is drawn past its partner.
    double const label_size;
    // how far left to shunt a left-side label.

  };  // end of spectrum_drawer class.


  /* instances of this class are basically vector arrays of instances of the
   * particle_property_set class. this class handles unexpected input by
   * adding extra particles for PDG codes that it hadn't known about so far
   * (calling them "unknown_particle_1", "unknown_particle_2",
   * "unknown_particle_3"...).
   */
  class EW_scale_spectrum
  {

  public:

    EW_scale_spectrum()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    EW_scale_spectrum( EW_scale_spectrum const* const copy_source )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~EW_scale_spectrum()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    become_copy_of( EW_scale_spectrum const* const copy_source )
    /* this copies in all the data from the given EW_scale_spectrum, *keeping*
     * all *this* EW_scale_spectrum's pointers, updating *their* data, so that
     * for instance any shortcut pointers to this EW_scale_spectrum's particle
     * property sets will continue to work.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    int
    get_number_of_unknown_particles()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    std::vector< particle_property_set* > const*
    get_particle_property_sets()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set*
    get_particle_property_set( int sought_PDG_code )
    /* this returns the pointer to the particle_property_set with the
     * sought PDG code.  it returns NULL if the PDG code was not found. if
     * provided with a negative PDG code, it returns a pointer to the particle
     * property set with the positive PDG code, if it is found.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    inspect_particle_property_set( int sought_PDG_code )
    const
    /* this returns the pointer to the particle_property_set with the
     * sought PDG code.  it returns NULL if the PDG code was not found. if
     * provided with a negative PDG code, it returns a pointer to the particle
     * property set with the positive PDG code, if it is found.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_PDG_coded_mass( int const input_PDG_code )
    const
    /* this checks to see if there is a particle property set with the input
     * PDG code, or the positive PDG code if a negative PDG code is provided,
     * then if there is, returns its mass, otherwise
     * CppSLHA_global::really_wrong_value().
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    std::string const*
    get_PDG_coded_name( int const input_PDG_code )
    const
    /* this checks to see if there is a particle property set with the input
     * PDG code, then if there is, returns its name, otherwise
     * "error!".  this function accounts for negative PDG codes by
     * returning the antiparticle's name.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    std::string const*
    get_PDG_coded_LaTeX_name( int const input_PDG_code )
    const
    /* this checks to see if there is a particle property set with the input
     * PDG code, then if there is, returns its LaTeX name, otherwise
     * "error!".  this function accounts for negative PDG codes by
     * returning the antiparticle's LaTeX name.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    record_PDG_coded_mass( int const input_PDG_code,
                           double const input_mass )
    /* this records a mass for a particle specified by its PDG code. if the
     * antiparticle's PDG code is given, the mass is instead recorded for the
     * particle.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    particle_property_set*
    add_particle_property_set( int const particle_PDG_code,
                               std::string const* const particle_name,
                               bool const counts_as_stable_invisible,
                               bool const counts_as_jet,
                               bool const counts_as_light_lepton,
                               bool const counts_as_self_conjugate )
    /* this creates a new instance of a particle_property_set, adds a
     * pointer to it to the list of pointers to instances of the
     * particle_property_set class, & returns the pointer. if it receives a
     * negative PDG code, it records the positive PDG code & swaps the names it
     * received, e.g. if it received PDG code -24 & particle name "W_minus", it
     * would record the particle with code +24, name "anti-W_minus" &
     * antiparticle name "W_minus".
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set*
    add_particle_property_set( int const particle_PDG_code,
                               std::string const* const particle_name,
                               std::string const* const antiparticle_name,
                               bool const counts_as_stable_invisible,
                               bool const counts_as_jet,
                               bool const counts_as_light_lepton,
                               bool const counts_as_self_conjugate )
    /* this creates a new instance of a particle_property_set, adds a pointer
     * to it to the list of pointers to instances of the particle_property_set
     * class, & returns the pointer. if it receives a negative PDG code, it
     * records the positive PDG code & swaps the names it received, e.g. if it
     * received PDG code -24 & particle name "W_minus" & antiparticle name
     * "W_plus", it would record the particle with code +24, particle name
     * "W_plus" & antiparticle name "W_minus".
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set*
    add_particle_property_set( int const particle_PDG_code,
                               std::string const* const particle_name,
                               std::string const* const antiparticle_name,
                               std::string const* const LaTeX_name,
                               std::string const* const LaTeX_antiname,
                               int const column_to_draw_in,
                               std::string const* const color_to_draw_with,
                               bool const counts_as_stable_invisible,
                               bool const counts_as_jet,
                               bool const counts_as_light_lepton,
                               bool const counts_as_self_conjugate )
    /* this creates a new instance of a particle_property_set, adds a pointer
     * to it to the list of pointers to instances of the particle_property_set
     * class, & returns the pointer. if it receives a negative PDG code, it
     * records the positive PDG code & swaps the names it received, e.g. if it
     * received PDG code -24 & particle name "W_minus" & antiparticle name
     * "W_plus", it would record the particle with code +24, particle name
     * "W_plus" & antiparticle name "W_minus". it also swaps the LaTeX name
     * & antiname.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set*
    add_particle_property_set( int const particle_PDG_code,
                               std::string const* const particle_name,
                               std::string const* const antiparticle_name,
                               std::string const* const LaTeX_name,
                               std::string const* const LaTeX_antiname,
                               bool const counts_as_stable_invisible,
                               bool const counts_as_jet,
                               bool const counts_as_light_lepton,
                               bool const counts_as_self_conjugate )
    /* this creates a new instance of a particle_property_set, adds a pointer
     * to it to the list of pointers to instances of the particle_property_set
     * class, & returns the pointer. if it receives a negative PDG code, it
     * records the positive PDG code & swaps the names it received, e.g. if it
     * received PDG code -24 & particle name "W_minus" & antiparticle name
     * "W_plus", it would record the particle with code +24, particle name
     * "W_plus" & antiparticle name "W_minus". it also swaps the LaTeX name
     * & antiname. this version doesn't assign a drawing column or color
     * because it is used by the copy constructor, which finds out the drawing
     * properties by copying the spectrum drawer.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set*
    add_unknown_particle_property_set( int const particle_PDG_code )
    /* this creates a new instance of a particle_property_set, with a
     * name indicating that its PDG code was unknown to this point, not
     * counting as stable-&-invisible, as a jet, as a light lepton or
     * self-conjugate.  it then adds a pointer to it to the list of pointers
     * to instances of the particle_property_set class, & returns the
     * pointer.  if given a negative PDG code, this function swaps its sign.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    spectrum_drawer const*
    get_spectrum_drawer()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    draw_spectrum( std::string const output_eps_filename,
                   std::string const path_to_gnuplot,
                   std::string const path_to_latex,
                   std::string const path_to_dvips,
                   std::string const path_to_ps2eps )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    record_PDG_coded_decay_width( int const input_PDG_code,
                                  double const input_decay_width )
    /* this records a decay width for a particle specified by its PDG code.
     * if the antiparticle's PDG code is given, the decay width is instead
     * recorded for the particle.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    record_PDG_coded_direct_decay( int const decayer_PDG_code,
                                   int const product_PDG_code_one,
                                   int const product_PDG_code_two,
                                   double const input_branching_ratio )
    /* this adds a 2-body decay to a particle property set with the specified
     * PDG code, or if no particle property set with that PDG code exists, it
     * creates a new property set & adds the decay to that. if an
     * antiparticle's PDG code is given, the decay for the particle (with
     * decay product particles replaced by antiparticles & decay product
     * antiparticles replaced by particles) is recorded instead.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    record_PDG_coded_direct_decay( int const decayer_PDG_code,
                                   int const product_PDG_code_one,
                                   int const product_PDG_code_two,
                                   int const product_PDG_code_three,
                                   double const input_branching_ratio )
    /* this adds a 3-body decay to a particle property set with the specified
     * PDG code, or if no particle property set with that PDG code exists, it
     * creates a new property set & adds the decay to that. if an
     * antiparticle's PDG code is given, the decay for the particle (with
     * decay product particles replaced by antiparticles & decay product
     * antiparticles replaced by particles) is recorded instead.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    std::string
    get_decay_as_string( particle_decay const* const decay_to_have_as_string )
    const
    // this returns a string that is the decay with names as well as PDG codes.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    int
    find_all_cascade_decays( double const minimum_branching_ratio_to_keep )
    /* this goes through each possible set of final stable particles from each
     * particle's decay. if the branching ratio for any particular cascade
     * decay goes below miminum_branching_ratio_to_keep, the cascade is
     * deleted.
     * N.B. it is up to the user to ensure an appropriate minimum branching
     * ratio is applied! if there are several identical final-state cascades
     * that are just below the threshold but add together to be over the
     * threshold, they will not be kept! so the user needs to give a bit of
     * leeway, such as picking an order of magnitude smaller a threshold.
     * I'd put this into the particle_property_set class, but each
     * particle needs to know about all the other particles.
     * actually, currently it's not particularly useful, since experimental
     * cuts will affect different cascades differently.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    std::vector< particle_property_set* > particle_property_sets;
    // this is a pointer to a vector array of pointers to
    // particle_property_set instances.

    int number_of_unknown_particles;

    spectrum_drawer spectrum_drawing;


    void
    calculate_all_cascades( particle_property_set* given_particle,
                            double const minimum_branching_ratio_to_keep )
    /* this goes through each possible set of final stable particles for the
     * given particle's decay. if the branching ratio for any particular
     * cascade decay goes below miminum_branching_ratio_to_keep, the cascade
     * is deleted. it returns the number of cascade decays the particle has.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    bool
    decay_unstable_particle_in_cascade( particle_property_set* given_particle,
                                 double const minimum_branching_ratio_to_keep )
    /* this looks for the 1st unstable particle in the cascade decay set of
     * given_particle, then adds in every possible cascade decay based on that
     * of given_particle with the unstable particle's direct decay products
     * replacing the particle, then removes the original decay. it returns
     * true if there was an unstable particle which was replaced, false if
     * there was no unstable particle to replace.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    add_unstable_particle_substituted_with_decay_products(
                                    particle_decay_set_handler* given_cascades,
                   std::list< particle_decay* >::iterator cascade_with_decayer,
                                                   int const decaying_particle,
                          std::list< particle_decay* > const* decays_to_insert,
                                 double const minimum_branching_ratio_to_keep )
    /* this adds to given_cascades a set of cascades which are all its own
     * cascades which have decaying_particle, but with decaying_particle
     * replaced by each decay in decays_to_insert, except any decays which end
     * up with BRs smaller than minimum_branching_ratio_to_keep. it also checks
     * that decaying_particle is not among its own decay products.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    substitute_decaying_particle_with_products(
                                             particle_decay* cascade_duplicate,
                                  std::list< int >::iterator decaying_particle,
                                                particle_decay* decay_products,
                                 double const minimum_branching_ratio_to_keep )
    // this substitutes the decayer's PDG code with the given decay set.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };

  class EW_scale_SM_spectrum : public EW_scale_spectrum
  {

  public:

    // constructor (creates an instance of the EW_scale_spectrum class
    // & then adds in the SM particles):
    EW_scale_SM_spectrum()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~EW_scale_SM_spectrum()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    particle_property_set const*
    get_down()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_up()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_strange()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_charm()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_bottom()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_top()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_electron()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_electron_neutrino()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_muon()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_muon_neutrino()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_tau_lepton()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_tau_neutrino()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_gluon()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_photon()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_Z()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_W_plus()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_EWSB_scalar()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    particle_property_set* down_pointer;
    particle_property_set* up_pointer;
    particle_property_set* strange_pointer;
    particle_property_set* charm_pointer;
    particle_property_set* bottom_pointer;
    particle_property_set* top_pointer;
    particle_property_set* electron_pointer;
    particle_property_set* electron_neutrino_pointer;
    particle_property_set* muon_pointer;
    particle_property_set* muon_neutrino_pointer;
    particle_property_set* tau_lepton_pointer;
    particle_property_set* tau_neutrino_pointer;
    particle_property_set* gluon_pointer;
    particle_property_set* photon_pointer;
    particle_property_set* Z_pointer;
    particle_property_set* W_plus_pointer;
    particle_property_set* EWSB_scalar_pointer;

  };

  class EW_scale_MSSM_spectrum : public EW_scale_SM_spectrum
  {

  public:

    // constructor (this constructs a EW_scale_SM_spectrum & then
    // adds the extra MSSM particles):
    EW_scale_MSSM_spectrum()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~EW_scale_MSSM_spectrum()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    particle_property_set*
    get_lighter_particle( particle_property_set* const compared_particle,
                          int const comparison_PDG_code )
    // this compares the masses of 2 particles & returns a pointer to the
    // lighter particle's property set.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    particle_property_set*
    update_stable_invisible()
    /* this checks the sneutrinos & lightest neutralino for which is lightest,
     * & returns that, while also setting its is_stable_invisible to true, if
     * this returned sparticle has an empty decay list. this *can* be used to
     * identify if the LSP is cosmologically allowed, by running
     * ->counts_as_stable_invisible() on the pointer which this method
     * returns,
     * e.g. if update_stable_invisible()->counts_as_stable_invisible() is
     * false, then the lightest cosmologically-acceptable (in the sense of
     * being neutral to the color & electromagnetism forces) dark matter
     * candidate is unstable, & thus cannot be dark matter (because it decays,
     * into something, which, assuming R-parity is conserved in all the
     * recorded decays, is presumably an electromagnetically- and/or
     * color-charged sparticle as well as SM particles).
     * actually, there is a crazy fringe case where that won't work: if for
     * example the lightest sneutrino is lighter than the neutralinos, but the
     * LSP is actually stop_1, and the lightest sneutrino mass is less than
     * the charm mass heavier the stop_1 mass, and whichever spectrum
     * generator that worked out the branching ratios didn't catch the 5-body
     * decay into a neutrino, stop_1, positron, another neutrino, and an up
     * quark through a virtual top, or the 3-body decay through the
     * sneutrino-electron-chargino lepton-flavor-violating vertex followed by
     * the V_td-suppressed chargino-stop-down vertex. honestly though, I don't
     * think it's the job of this code to catch such psychotic cases.
     * this function could be more comprehensive. I might sort it out later.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    particle_property_set const*
    get_light_neutral_EWSB_scalar()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_heavy_neutral_EWSB_scalar()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_neutral_EWSB_pseudoscalar()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_charged_EWSB_scalar()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_gluino()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_sdown_L()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_sup_L()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_sstrange_L()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_scharm_L()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_sbottom_one()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_stop_one()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_sdown_R()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_sup_R()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_sstrange_R()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_scharm_R()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_sbottom_two()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_stop_two()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_selectron_L()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_electron_sneutrino_L()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_smuon_L()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_muon_sneutrino_L()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_stau_one()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_tau_sneutrino_L()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_selectron_R()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_electron_sneutrino_R()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_smuon_R()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_muon_sneutrino_R()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_stau_two()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_tau_sneutrino_R()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_neutralino_one()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_neutralino_two()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_neutralino_three()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_neutralino_four()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_chargino_one()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_chargino_two()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    particle_property_set* light_neutral_EWSB_scalar_pointer;
    particle_property_set* heavy_neutral_EWSB_scalar_pointer;
    particle_property_set* neutral_EWSB_pseudoscalar_pointer;
    particle_property_set* charged_EWSB_scalar_pointer;
    particle_property_set* gluino_pointer;
    particle_property_set* sdown_L_pointer;
    particle_property_set* sup_L_pointer;
    particle_property_set* sstrange_L_pointer;
    particle_property_set* scharm_L_pointer;
    particle_property_set* sbottom_one_pointer;
    particle_property_set* stop_one_pointer;
    particle_property_set* sdown_R_pointer;
    particle_property_set* sup_R_pointer;
    particle_property_set* sstrange_R_pointer;
    particle_property_set* scharm_R_pointer;
    particle_property_set* sbottom_two_pointer;
    particle_property_set* stop_two_pointer;
    particle_property_set* selectron_L_pointer;
    particle_property_set* electron_sneutrino_L_pointer;
    particle_property_set* smuon_L_pointer;
    particle_property_set* muon_sneutrino_L_pointer;
    particle_property_set* stau_one_pointer;
    particle_property_set* tau_sneutrino_L_pointer;
    particle_property_set* selectron_R_pointer;
    particle_property_set* electron_sneutrino_R_pointer;
    particle_property_set* smuon_R_pointer;
    particle_property_set* muon_sneutrino_R_pointer;
    particle_property_set* stau_two_pointer;
    particle_property_set* tau_sneutrino_R_pointer;
    particle_property_set* neutralino_one_pointer;
    particle_property_set* neutralino_two_pointer;
    particle_property_set* neutralino_three_pointer;
    particle_property_set* neutralino_four_pointer;
    particle_property_set* chargino_one_pointer;
    particle_property_set* chargino_two_pointer;

  };

  class EW_scale_FVMSSM_spectrum : public EW_scale_SM_spectrum
  {

  public:

    // constructor (this constructs a EW_scale_MSSM_spectrum & then
    // adds the extra NMSSM particles):
    EW_scale_FVMSSM_spectrum()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~EW_scale_FVMSSM_spectrum()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    particle_property_set const*
    get_light_neutral_EWSB_scalar()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_heavy_neutral_EWSB_scalar()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_neutral_EWSB_pseudoscalar()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_charged_EWSB_scalar()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_gluino()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_sdown_one()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_sup_one()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_sdown_two()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_sup_two()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_sdown_three()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_sup_three()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_sdown_four()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_sup_four()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_sdown_five()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_sup_five()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_sdown_six()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_sup_six()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_selectron_one()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_sneutrino_one()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_selectron_two()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_sneutrino_two()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_selectron_three()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_sneutrino_three()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_selectron_four()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_sneutrino_four()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_selectron_five()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_sneutrino_five()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_selectron_six()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_sneutrino_six()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_neutralino_one()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_neutralino_two()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_neutralino_three()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_neutralino_four()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_chargino_one()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_chargino_two()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    particle_property_set* light_neutral_EWSB_scalar_pointer;
    particle_property_set* heavy_neutral_EWSB_scalar_pointer;
    particle_property_set* neutral_EWSB_pseudoscalar_pointer;
    particle_property_set* charged_EWSB_scalar_pointer;
    particle_property_set* gluino_pointer;
    particle_property_set* sdown_one_pointer;
    particle_property_set* sdown_two_pointer;
    particle_property_set* sdown_three_pointer;
    particle_property_set* sdown_four_pointer;
    particle_property_set* sdown_five_pointer;
    particle_property_set* sdown_six_pointer;
    particle_property_set* sup_one_pointer;
    particle_property_set* sup_two_pointer;
    particle_property_set* sup_three_pointer;
    particle_property_set* sup_four_pointer;
    particle_property_set* sup_five_pointer;
    particle_property_set* sup_six_pointer;
    particle_property_set* selectron_one_pointer;
    particle_property_set* selectron_two_pointer;
    particle_property_set* selectron_three_pointer;
    particle_property_set* selectron_four_pointer;
    particle_property_set* selectron_five_pointer;
    particle_property_set* selectron_six_pointer;
    particle_property_set* sneutrino_one_pointer;
    particle_property_set* sneutrino_two_pointer;
    particle_property_set* sneutrino_three_pointer;
    particle_property_set* sneutrino_four_pointer;
    particle_property_set* sneutrino_five_pointer;
    particle_property_set* sneutrino_six_pointer;
    particle_property_set* neutralino_one_pointer;
    particle_property_set* neutralino_two_pointer;
    particle_property_set* neutralino_three_pointer;
    particle_property_set* neutralino_four_pointer;
    particle_property_set* chargino_one_pointer;
    particle_property_set* chargino_two_pointer;

  };

  class EW_scale_NMSSM_spectrum : public EW_scale_MSSM_spectrum
  {

  public:

    // constructor (this constructs a EW_scale_MSSM_spectrum & then
    // adds the extra NMSSM particles):
    EW_scale_NMSSM_spectrum()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~EW_scale_NMSSM_spectrum()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    particle_property_set const*
    get_singlet_scalar()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_singlet_pseudoscalar()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_neutralino_five()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    particle_property_set* singlet_scalar_pointer;
    particle_property_set* singlet_pseudoscalar_pointer;
    particle_property_set* neutralino_five_pointer;

  };

  class EW_scale_FVNMSSM_spectrum : public EW_scale_FVMSSM_spectrum
  {

  public:

    // constructor (this constructs a EW_scale_FVMSSM_spectrum & then
    // adds the extra FVNMSSM particles):
    EW_scale_FVNMSSM_spectrum()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~EW_scale_FVNMSSM_spectrum()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    particle_property_set const*
    get_singlet_scalar()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_singlet_pseudoscalar()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    particle_property_set const*
    get_neutralino_five()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    particle_property_set* singlet_scalar_pointer;
    particle_property_set* singlet_pseudoscalar_pointer;
    particle_property_set* neutralino_five_pointer;

  };





  // inline functions:

  inline std::string const*
  drawing_label::get_label()
  const
  {

    return &label_string;

  }

  inline void
  drawing_label::set_label( std::string const* const given_label )
  {

    label_string.assign( *given_label );

  }

  inline double
  drawing_label::get_label_height()
  const
  {

    return label_height;

  }

  inline void
  drawing_label::set_label_height( double const given_height )
  {

    label_height = given_height;

  }



  inline std::string const*
  drawing_line::get_color()
  const
  {

    return &color_string;

  }

  inline double
  drawing_line::get_line_height()
  const
  {

    return line_height;

  }

  inline double
  drawing_line::get_label_height()
  const
  {

    return line_label->get_label_height();

  }



  inline bool
  drawing_column::is_empty()
  const
  // this returns true if there are no lines to draw, false otherwise.
  {

    if( lines.empty()
        &&
        particles.empty() )
      {

        return true;

      }
    else
      {

        return false;

      }

  }


  inline void
  drawing_column::add_particle(
                             particle_property_set const* const given_particle,
                                std::string const* const given_color )
  {

    std::string* color_string = new std::string( *given_color );

    particles.push_back( new std::pair< particle_property_set const* const,
                                        std::string* >( given_particle,
                                                        color_string ) );

  }

  inline std::list< drawing_line* > const*
  drawing_column::get_lines()
  const
  {

    return &lines;

  }

  inline std::list< drawing_label* > const*
  drawing_column::get_labels()
  const
  {

    return &labels;

  }

  inline std::list< std::pair< particle_property_set const* const,
                               std::string* >* >*
  drawing_column::get_particles()
  {

    return &particles;

  }


  inline double
  drawing_column::get_highest_height()
  const
  // this goes through all the particles stored in the column & finds the
  // highest height for drawing a particle's line.
  {

    double column_height;
    double highest_height = 0.0;

    for( std::list< std::pair< particle_property_set const* const,
                               std::string* >* >::const_iterator
         particle_iterator = particles.begin();
         particles.end() != particle_iterator;
         ++particle_iterator )
      {

        column_height = (*particle_iterator)->first->get_mass();

        if( column_height > highest_height )
          {

            highest_height = column_height;

          }

      }

    return highest_height;

  }



  inline void
  spectrum_drawer::reset()
  {

    for( std::vector< drawing_column* >::iterator
         deletion_iterator = column_set.begin();
         column_set.end() != deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

    column_set.clear();

    highest_height = 0.0;

  }


  inline double
  spectrum_drawer::get_highest_height()
  const
  {

    return highest_height;

  }

  inline std::vector< drawing_column* > const*
  spectrum_drawer::get_column_set()
  const
  {

    return &column_set;

  }

  inline void
  spectrum_drawer::add_particle(
                       particle_property_set const* const particle_to_be_added,
                                 int const column_to_draw_in,
                                 std::string const* const color_to_draw_with )
  {

    // debugging:
    /*std::cout
    << std::endl
    << "debugging: spectrum_drawer::add_particle( "
    << particle_to_be_added->getName() << ", " << column_to_draw_in
    << ", " << color_to_draw_with << " ) has been called. prior to this,"
    << " drawing column 1 has ";
    if( 1 < column_set.size() )
      {

        std::cout << column_set.at( 1 )->get_particles()->size();

      }
    else
      {

        std::cout << "0";

      }
    std::cout << " particles.";
    std::cout << std::endl;
    */

    if( 0 <= column_to_draw_in )
      {

        while( column_set.size() <= (unsigned int)column_to_draw_in )
          {

            column_set.push_back( new drawing_column() );

          }
        column_set.at( column_to_draw_in )->add_particle( particle_to_be_added,
                                                          color_to_draw_with );

      }

    if( particle_to_be_added->has_mass_been_recorded() )
      {

        if( particle_to_be_added->get_mass() > highest_height )
          {

            highest_height = particle_to_be_added->get_mass();

          }

      }

  }


  inline double
  spectrum_drawer::find_highest_height()
  // this goes through all the particles stored in all the columns & finds
  // the highest height for drawing a particle's line.
  {

    double column_height;

    for( std::vector< drawing_column* >::const_iterator
         column_iterator = column_set.begin();
         column_set.end() != column_iterator;
         ++column_iterator )
      {

        column_height = (*column_iterator)->get_highest_height();

        if( column_height > highest_height )
          {

            highest_height = column_height;

          }

      }

    return highest_height;

  }



  inline int
  EW_scale_spectrum::get_number_of_unknown_particles()
  const
  {

    return number_of_unknown_particles;

  }

  inline std::vector< particle_property_set* > const*
  EW_scale_spectrum::get_particle_property_sets()
  const
  {

    return &particle_property_sets;

  }


  inline spectrum_drawer const*
  EW_scale_spectrum::get_spectrum_drawer()
  const
  {

    return &spectrum_drawing;

  }

  inline void
  EW_scale_spectrum::draw_spectrum( std::string const output_eps_filename,
                                    std::string const path_to_gnuplot,
                                    std::string const path_to_latex,
                                    std::string const path_to_dvips,
                                    std::string const path_to_ps2eps )
  {

    spectrum_drawing.draw_spectrum( &output_eps_filename,
                                    &path_to_gnuplot,
                                    &path_to_latex,
                                    &path_to_dvips,
                                    &path_to_ps2eps );

  }


  inline int
  EW_scale_spectrum::find_all_cascade_decays(
                                 double const minimum_branching_ratio_to_keep )
  /* this goes through each possible set of final stable particles from each
   * particle's decay. if the branching ratio for any particular cascade
   * decay goes below miminum_branching_ratio_to_keep, the cascade is
   * deleted.
   * N.B. it is up to the user to ensure an appropriate minimum branching
   * ratio is applied! if there are several identical final-state cascades
   * that are just below the threshold but add together to be over the
   * threshold, they will not be kept! so the user needs to give a bit of
   * leeway, such as picking an order of magnitude smaller a threshold.
   * I'd put this into the particle_property_set class, but each
   * particle needs to know about all the other particles.
   * actually, currently it's not particularly useful, since experimental
   * cuts will affect different cascades differently.
   */
  {

    if( 0.00001 > minimum_branching_ratio_to_keep )
    // if this function was called with an insane argument...
      {

        // alert the user, out of politeness.
        std::cout
        << std::endl
        << "CppSLHA::warning! find_all_cascade_decays was called with a"
        << " minimum_branching_ratio_to_keep of less than 0.00001! this"
        << " could result in a VERY long wait, since there could be up to"
        << " 10^5 or 10^6 decay chains for the gluino for example, and each"
        << " of these 10^6 chains has to be compared to the other 10^6 - 1"
        << " chains leading to something like 10^12 comparisons. in BOL's"
        << " humble opinion, this is a real waste of time...";
        std::cout
        << std::endl
        << "(you might want to press ctrl-c to kill this program, and"
        << " re-write the relevant bit of code.)";
        std::cout << std::endl;

      }

    int total_number_of_cascade_decays = 0;
    // this keeps track of how many cascade decays we are recording.

    for( std::vector< particle_property_set* >::iterator
         cascader_iterator = particle_property_sets.begin();
         particle_property_sets.end() > cascader_iterator;
         ++cascader_iterator )
      // for each particle in the spectrum...
      {

        calculate_all_cascades( *cascader_iterator,
                                minimum_branching_ratio_to_keep );

        total_number_of_cascade_decays
        += (*cascader_iterator)->get_cascade_decays()->size();

      }

    return total_number_of_cascade_decays;

  }

  inline void
  EW_scale_spectrum::calculate_all_cascades(
                                         particle_property_set* given_particle,
                                 double const minimum_branching_ratio_to_keep )
  /* this goes through each possible set of final stable particles for the
   * given particle's decay. if the branching ratio for any particular
   * cascade decay goes below miminum_branching_ratio_to_keep, the cascade
   * is deleted. it returns the number of cascade decays the particle has.
   */
  {

    // if no cascade decays for the given particle have been yet worked out,
    // start with the direct decays & carry on from there.
    if( given_particle->get_cascade_decays()->empty() )
      {

        given_particle->get_cascade_handler()->become_copy_of(
                                  given_particle->get_direct_decay_handler() );

      }

    while( decay_unstable_particle_in_cascade( given_particle,
                                            minimum_branching_ratio_to_keep ) )
      // until all the unstable particles have been replaced by their decay
      // products, keep replacing them.
      {

        // all the action for the while loop actually happens in the
        // condition expression.

      }

  }



  inline particle_property_set const*
  EW_scale_SM_spectrum::get_down()
  const
  {

    return down_pointer;

  }

  inline particle_property_set const*
  EW_scale_SM_spectrum::get_up()
  const
  {

    return up_pointer;

  }

  inline particle_property_set const*
  EW_scale_SM_spectrum::get_strange()
  const
  {

    return strange_pointer;

  }

  inline particle_property_set const*
  EW_scale_SM_spectrum::get_charm()
  const
  {

    return charm_pointer;

  }

  inline particle_property_set const*
  EW_scale_SM_spectrum::get_bottom()
  const
  {

    return bottom_pointer;

  }

  inline particle_property_set const*
  EW_scale_SM_spectrum::get_top()
  const
  {

    return top_pointer;

  }

  inline particle_property_set const*
  EW_scale_SM_spectrum::get_electron()
  const
  {

    return electron_pointer;

  }

  inline particle_property_set const*
  EW_scale_SM_spectrum::get_electron_neutrino()
  const
  {

    return electron_neutrino_pointer;

  }

  inline particle_property_set const*
  EW_scale_SM_spectrum::get_muon()
  const
  {

    return muon_pointer;

  }

  inline particle_property_set const*
  EW_scale_SM_spectrum::get_muon_neutrino()
  const
  {

    return muon_neutrino_pointer;

  }

  inline particle_property_set const*
  EW_scale_SM_spectrum::get_tau_lepton()
  const
  {

    return tau_lepton_pointer;

  }

  inline particle_property_set const*
  EW_scale_SM_spectrum::get_tau_neutrino()
  const
  {

    return tau_neutrino_pointer;

  }

  inline particle_property_set const*
  EW_scale_SM_spectrum::get_gluon()
  const
  {

    return gluon_pointer;

  }

  inline particle_property_set const*
  EW_scale_SM_spectrum::get_photon()
  const
  {

    return photon_pointer;

  }

  inline particle_property_set const*
  EW_scale_SM_spectrum::get_Z()
  const
  {

    return Z_pointer;

  }

  inline particle_property_set const*
  EW_scale_SM_spectrum::get_W_plus()
  const
  {

    return W_plus_pointer;

  }

  inline particle_property_set const*
  EW_scale_SM_spectrum::get_EWSB_scalar()
  const
  {

    return EWSB_scalar_pointer;

  }



  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_light_neutral_EWSB_scalar()
  const
  {

    return light_neutral_EWSB_scalar_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_heavy_neutral_EWSB_scalar()
  const
  {

    return heavy_neutral_EWSB_scalar_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_neutral_EWSB_pseudoscalar()
  const
  {

    return neutral_EWSB_pseudoscalar_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_charged_EWSB_scalar()
  const
  {

    return charged_EWSB_scalar_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_gluino()
  const
  {

    return gluino_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_sdown_L()
  const
  {

    return sdown_L_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_sup_L()
  const
  {

    return sup_L_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_sstrange_L()
  const
  {

    return sstrange_L_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_scharm_L()
  const
  {

    return scharm_L_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_sbottom_one()
  const
  {

    return sbottom_one_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_stop_one()
  const
  {

    return stop_one_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_sdown_R()
  const
  {

    return sdown_R_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_sup_R()
  const
  {

    return sup_R_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_sstrange_R()
  const
  {

    return sstrange_R_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_scharm_R()
  const
  {

    return scharm_R_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_sbottom_two()
  const
  {

    return sbottom_two_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_stop_two()
  const
  {

    return stop_two_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_selectron_L()
  const
  {

    return selectron_L_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_electron_sneutrino_L()
  const
  {

    return electron_sneutrino_L_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_smuon_L()
  const
  {

    return smuon_L_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_muon_sneutrino_L()
  const
  {

    return muon_sneutrino_L_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_stau_one()
  const
  {

    return stau_one_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_tau_sneutrino_L()
  const
  {

    return tau_sneutrino_L_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_selectron_R()
  const
  {

    return selectron_R_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_electron_sneutrino_R()
  const
  {

    return electron_sneutrino_R_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_smuon_R()
  const
  {

    return smuon_R_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_muon_sneutrino_R()
  const
  {

    return muon_sneutrino_R_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_stau_two()
  const
  {

    return stau_two_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_tau_sneutrino_R()
  const
  {

    return tau_sneutrino_R_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_neutralino_one()
  const
  {

    return neutralino_one_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_neutralino_two()
  const
  {

    return neutralino_two_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_neutralino_three()
  const
  {

    return neutralino_three_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_neutralino_four()
  const
  {

    return neutralino_four_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_chargino_one()
  const
  {

    return chargino_one_pointer;

  }

  inline particle_property_set const*
  EW_scale_MSSM_spectrum::get_chargino_two()
  const
  {

    return chargino_two_pointer;

  }



  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_light_neutral_EWSB_scalar()
  const
  {

    return light_neutral_EWSB_scalar_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_heavy_neutral_EWSB_scalar()
  const
  {

    return heavy_neutral_EWSB_scalar_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_neutral_EWSB_pseudoscalar()
  const
  {

    return neutral_EWSB_pseudoscalar_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_charged_EWSB_scalar()
  const
  {

    return charged_EWSB_scalar_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_gluino()
  const
  {

    return gluino_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_sdown_one()
  const
  {

    return sdown_one_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_sup_one()
  const
  {

    return sup_one_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_sdown_two()
  const
  {

    return sdown_two_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_sup_two()
  const
  {

    return sup_two_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_sdown_three()
  const
  {

    return sdown_three_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_sup_three()
  const
  {

    return sup_three_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_sdown_four()
  const
  {

    return sdown_four_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_sup_four()
  const
  {

    return sup_four_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_sdown_five()
  const
  {

    return sdown_five_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_sup_five()
  const
  {

    return sup_five_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_sdown_six()
  const
  {

    return sdown_six_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_sup_six()
  const
  {

    return sup_six_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_selectron_one()
  const
  {

    return selectron_one_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_sneutrino_one()
  const
  {

    return sneutrino_one_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_selectron_two()
  const
  {

    return selectron_two_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_sneutrino_two()
  const
  {

    return sneutrino_two_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_selectron_three()
  const
  {

    return selectron_three_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_sneutrino_three()
  const
  {

    return sneutrino_three_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_selectron_four()
  const
  {

    return selectron_four_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_sneutrino_four()
  const
  {

    return sneutrino_four_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_selectron_five()
  const
  {

    return selectron_five_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_sneutrino_five()
  const
  {

    return sneutrino_five_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_selectron_six()
  const
  {

    return selectron_six_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_sneutrino_six()
  const
  {

    return sneutrino_six_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_neutralino_one()
  const
  {

    return neutralino_one_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_neutralino_two()
  const
  {

    return neutralino_two_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_neutralino_three()
  const
  {

    return neutralino_three_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_neutralino_four()
  const
  {

    return neutralino_four_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_chargino_one()
  const
  {

    return chargino_one_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVMSSM_spectrum::get_chargino_two()
  const
  {

    return chargino_two_pointer;

  }



  inline particle_property_set const*
  EW_scale_NMSSM_spectrum::get_singlet_scalar()
  const
  {

    return singlet_scalar_pointer;

  }

  inline particle_property_set const*
  EW_scale_NMSSM_spectrum::get_singlet_pseudoscalar()
  const
  {

    return singlet_pseudoscalar_pointer;

  }

  inline particle_property_set const*
  EW_scale_NMSSM_spectrum::get_neutralino_five()
  const
  {

    return neutralino_five_pointer;

  }



  inline particle_property_set const*
  EW_scale_FVNMSSM_spectrum::get_singlet_scalar()
  const
  {

    return singlet_scalar_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVNMSSM_spectrum::get_singlet_pseudoscalar()
  const
  {

    return singlet_pseudoscalar_pointer;

  }

  inline particle_property_set const*
  EW_scale_FVNMSSM_spectrum::get_neutralino_five()
  const
  {

    return neutralino_five_pointer;

  }

} // end of CppSLHA namespace.

#endif /* CPPSLHA_EW_SCALE_SPECTRUM_HPP_ */
