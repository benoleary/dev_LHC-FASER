/*
 * CppSLHA.hpp
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

#ifndef CPPSLHA_HPP_
#define CPPSLHA_HPP_
// this should prevent this class being included multiple times.

/*
 * Ben O'Leary, originally 8th July 2009. this code is based in part on what
 * I learned from reading some private code by Martin Wiebusch
 * (mwiebusch@physik.rwth-aachen.de), along with examples from the
 * cplusplus.com website.
 */

/*
 * this is the header file for a C++ class that reads a supplied SLHA file &
 * stores the masses & decays of the sparticles in a
 * EW_scale_MSSM_spectrum instance & the rest in SLHA_BLOCK instances.
 * it has been organized so that it should be just this header file in the
 * directory where it is used, with the rest in the ./CppSLHA_stuff/
 * directory.
 */

/*********************************************************************
 *                                                                   *
 * IMPORTANT!                                                        *
 *                                                                   *
 * welcome, person who has decided to read this file.                *
 * please see README.CppSLHA.txt for information about the structure *
 * & intent of this code.                                            *
 *                                                                   *
 *********************************************************************/


// includes:
#include "CppSLHA_global_stuff.hpp"
#include "CppSLHA_EW_scale_spectrum.hpp"
#include "CppSLHA_BLOCK.hpp"
#include "CppSLHA_file_reading_stuff.hpp"


namespace CppSLHA
{

  // this class reads in a SLHA format file & stores the data internally,
  // in the manner explained above.
  class CppSLHA0
  {

  public:


    /* the constructor requires a file to open, & sets all entries to NaN or
     * not even assigned, in the case of std::vectors of entries & also
     * allocates the memory for the particle spectrum.
     */
    CppSLHA0( std::string const given_SLHA_file_name )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    /* this constructor requires a file to open, & this one takes a
     * pre-existing EW_scale_spectrum, which is how one can implement the NMSSM
     * or similar.
     */
    CppSLHA0( std::string const given_SLHA_file_name,
              EW_scale_spectrum* const given_particle_spectrum )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    /* there is no copy constructor because it would be a pain to write.
     * I cannot envision a situation where one could not simply open the
     * SLHA file again with a new instance if it were so important, though
     * I suppose it would be slow to write the SLHA file & then read it
     * in again. however, if someone really needs to do that, they can write
     * their own copy constructor...
     */

    virtual
    ~CppSLHA0()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    set_verbosity( bool const given_verbosity )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    particle_decay_set_handler const*
    get_DECAY( int const decayer_code )
    const
    // this returns a pointer to an instance of the particle_decay_set_handler
    // class (see particle_data.hpp).
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_branching_ratio( int const decayer_code,
                         int const product_one_code,
                         int const product_two_code )
    const
    /* this returns the branching ratio for a specified 2-body decay of a
     * particle with PDG MC code decayer_code into decay products with PDG
     * particle codes product_one_code & product_two_code.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_branching_ratio( int const decayer_code,
                         int const product_one_code,
                         int const product_two_code,
                         int const product_three_code )
    const
    /* this returns the branching ratio for a specified 3-body decay of a
     * particle with PDG particle code decayer_code into decay products with
     * PDG particle codes product_one_code, product_two_code &
     * product_three_code.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    EW_scale_spectrum*
    get_particle_spectrum()
    // this returns the particle spectrum as a pointer to an instance of the
    // EW_scale_spectrum class (see EW_scale_spectrum.hpp).
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    EW_scale_spectrum const*
    inspect_particle_spectrum()
    const
    // this returns the particle spectrum as a const pointer to an instance of
    // the EW_scale_spectrum class (see EW_scale_spectrum.hpp).
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    std::string const*
    get_SLHA_file_name()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_MASS( int const requested_element )
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    get_absolute_MASS( int const requested_element )
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    set_MASS( int const requested_element,
              double const set_value )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    draw_spectrum( std::string const output_eps_filename )
    /* this draws the spectrum using gnuplot, latex, dvips, & ps2eps. it
     * assumes that all the required executables are found in the specified
     * paths below.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    draw_spectrum( std::string const output_eps_filename,
                   std::string const path_to_gnuplot,
                   std::string const path_to_latex,
                   std::string const path_to_dvips,
                   std::string const path_to_ps2eps )
    /* this draws the spectrum using gnuplot, latex, dvips, & ps2eps. the
     * "path_to_X" strings should include the executable name, e.g.
     * /usr/bin/gnuplot for path_to_gnuplot.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    reset()
    // this resets all the entries apart from the SLHA file name.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    SLHA_BLOCK*
    get_BLOCK( std::string requested_BLOCK )
    // the input string is not const because it is converted to all uppercase
    // before using it as a search key.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    SLHA_BLOCK const*
    inspect_BLOCK( std::string requested_BLOCK )
    const
    // the input string is not const because it is converted to all uppercase
    // before using it as a search key.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    std::string
    get_MASS_as_string()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    std::string
    get_DECAY_as_string()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    bool
    read_file()
    // this reads in the SLHA file and returns true if it was successful.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    bool
    read_file( std::string const SLHA_file_name )
    // this reads in the specified SLHA file and returns true if it was
    // successful.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    bool
    write_file( std::string const output_SLHA_file_name )
    // this writes out all the BLOCKs to a file with the given name
    // (DECAY counts as a BLOCK name in this case).
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    bool
    write_file( std::string const output_SLHA_file_name,
                std::list< std::string > const* const list_of_BLOCKS )
    // this writes out the BLOCKs from the supplied list to a file with
    // the given name (DECAY counts as a BLOCK name in this case).
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    // STUFF FOR TESTING PURPOSES:

    int
    write_particle_spectrum( double const minimum_branching_ratio_to_keep,
                             std::string const output_file_name )
    // this writes all the stored spectrum data to a file with the provided
    // name. it returns the number of particle property sets written.
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    // END OF TESTING METHOD!


  protected:

    std::vector< SLHA_BLOCK* > BLOCKs;

    EW_scale_spectrum* particle_spectrum;
    EW_scale_spectrum* copy_of_given_spectrum;
    /* see CppSLHA_EW_scale_spectrum.hpp for details, but this class should
     * take care of any normal usage without worrying about the details of how
     * that class stores the particle masses.  actually, normal use that
     * requires decays probably does need to know the details of how the
     * decays are stored.  this class does provide basic access to the
     * decays, in the sense that it can return the branching ratio for a
     * specified decay (so 0.0 is returned for decays that do not happen),
     * but it's relatively restrictive, whereas the
     * particle_decay_set_handler class & the classes it relies on
     * allow for much more flexibility.
     * the reason for copy_of_given_spectrum is to allow easy resetting of
     * particle_spectrum. in the case of a non-MSSM particle spectrum, the
     * constructor is passed a pre-constructed EW_scale_spectrum object,
     * which is copied into copy_of_given_spectrum which never gets altered,
     * & copy_of_given_spectrum is used to initialize & to reset
     * particle_spectrum.
     */

    // masses: contained in particle_spectrum.
    // decays: contained in particle_spectrum.

    std::string SLHA_file_name;

    bool did_allocate_spectrum_memory;
    // this is to keep track of whether the destructor needs to deallocate
    // memory or not.

    bool is_verbose;

  }; // end of CppSLHA0 class.


  /* this derived class acknowledges all the BLOCKs described in the SLHA (1)
   * paper:
   * MODSEL, SMINPUTS, MINPAR, EXTPAR, NMIX, UMIX, VMIX, STOPMIX, SBOTMIX,
   * STAUMIX, ALPHA, HMIX, GAUGE, MSOFT, AU, AD, AE, YU, YD, YE & SPINFO.
   */
  class CppSLHA1 : public CppSLHA0
  {

  public:

    // the constructor requires a file to open.
    CppSLHA1( std::string const given_SLHA_file_name )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    /* this constructor requires a file to open, & this one takes a
     * pre-existing EW_scale_spectrum, which is how one can implement
     * the NMSSM or similar.
     */
    CppSLHA1( std::string const given_SLHA_file_name,
              EW_scale_spectrum* const given_particle_spectrum )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~CppSLHA1()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    initialize()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  }; // end of CppSLHA1 class


  /* this derived class acknowledges all the BLOCKs described in the SLHA (1)
   * paper & the SLHA2 paper.
   */
  class CppSLHA2 : public CppSLHA1
  {

  public:

    // the constructor requires a file to open.
    CppSLHA2( std::string const given_SLHA_file_name )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    /* this constructor requires a file to open, & this one takes a
     * pre-existing EW_scale_spectrum, which is how one can implement
     * the NMSSM or similar.
     */
    CppSLHA2( std::string const given_SLHA_file_name,
              EW_scale_spectrum* const given_particle_spectrum )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~CppSLHA2()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    void
    initialize()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  }; // end of CppSLHA2 class


    /* this derived class acknowledges all the BLOCKs described in the SLHA2
     * paper & also the SPHENOLOWENERGY BLOCK produced by SPheno with
     * low-energy observables, such as the anomalous magnetic moment of the
     * muon.
     */
    class CppSLHA2_SPheno : public CppSLHA2
    {

    public:

      // the constructor requires a file to open.
      CppSLHA2_SPheno( std::string const given_SLHA_file_name )
      /* code after the classes in this .hpp file, or in the .cpp file. */;

      /* this constructor requires a file to open, & this one takes a
       * pre-existing EW_scale_spectrum, which is how one can implement
       * the NMSSM or similar.
       */
      CppSLHA2_SPheno( std::string const given_SLHA_file_name,
                       EW_scale_spectrum* const given_particle_spectrum )
      /* code after the classes in this .hpp file, or in the .cpp file. */;

      ~CppSLHA2_SPheno()
      /* code after the classes in this .hpp file, or in the .cpp file. */;


      SPheno_cross_sections_BLOCK*
      get_SPheno_cross_sections_BLOCK()
      /* code after the classes in this .hpp file, or in the .cpp file. */;


      void
      initialize()
      /* code after the classes in this .hpp file, or in the .cpp file. */;


    protected:

      SPheno_cross_sections_BLOCK* shortcut_to_cross_sections;

    }; // end of CppSLHA2_SPheno class





    // inline functions:

    inline void
    CppSLHA0::set_verbosity( bool const given_verbosity )
      {

        is_verbose = given_verbosity;

      }


    inline particle_decay_set_handler const*
    CppSLHA0::get_DECAY( int const decayer_code )
    const
    // this returns a pointer to an instance of the particle_decay_set_handler
    // class (see particle_data.hpp).
    {

      return particle_spectrum->get_particle_property_set( decayer_code
                                                 )->get_direct_decay_handler();

    }

    inline double
    CppSLHA0::get_branching_ratio( int const decayer_code,
                                   int const product_one_code,
                                   int const product_two_code )
    const
    /* this returns the branching ratio for a specified 2-body decay of a
     * particle with PDG MC code decayer_code into decay products with PDG
     * particle codes product_one_code & product_two_code.
     */
    {

      return particle_spectrum->get_particle_property_set( decayer_code
                                                   )->get_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                                              product_one_code,
                                                            product_two_code );

    }

    inline double
    CppSLHA0::get_branching_ratio( int const decayer_code,
                                   int const product_one_code,
                                   int const product_two_code,
                                   int const product_three_code )
    const
    /* this returns the branching ratio for a specified 3-body decay of a
     * particle with PDG particle code decayer_code into decay products with
     * PDG particle codes product_one_code, product_two_code &
     * product_three_code.
     */
    {

      return particle_spectrum->get_particle_property_set( decayer_code
                                                   )->get_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                                              product_one_code,
                                                              product_two_code,
                                                          product_three_code );

    }

    inline EW_scale_spectrum*
    CppSLHA0::get_particle_spectrum()
    // this returns the particle spectrum as a pointer to an instance of the
    // EW_scale_spectrum class (see EW_scale_spectrum.hpp).
    {

      return particle_spectrum;

    }

    inline EW_scale_spectrum const*
    CppSLHA0::inspect_particle_spectrum()
    const
    // this returns the particle spectrum as a const pointer to an instance of
    // the EW_scale_spectrum class (see EW_scale_spectrum.hpp).
    {

      return particle_spectrum;

    }

    inline std::string const*
    CppSLHA0::get_SLHA_file_name()
    const
    {

      return &SLHA_file_name;

    }

    inline double
    CppSLHA0::get_MASS( int const requested_element )
    const
    {

      return particle_spectrum->get_PDG_coded_mass( requested_element );

    }

    inline double
    CppSLHA0::get_absolute_MASS( int const requested_element )
    const
    {

      double return_mass = get_MASS( requested_element );

      if( 0.0 > return_mass )
        {

          return -return_mass;

        }
      else
        {

          return return_mass;

        }

    }


    inline void
    CppSLHA0::set_MASS( int const requested_element,
                        double const set_value )
    {

      particle_spectrum->record_PDG_coded_mass( requested_element,
                                                set_value );

    }


    inline void
    CppSLHA0::draw_spectrum( std::string const output_eps_filename )
    /* this draws the spectrum using gnuplot, latex, dvips, & ps2eps. it
     * assumes that all the required executables are found in the specified
     * paths below.
     */
    {

      particle_spectrum->draw_spectrum( output_eps_filename,
                                        "/usr/bin/gnuplot",
                                        "/usr/bin/latex",
                                        "/usr/bin/dvips",
                                        "/usr/bin/ps2eps" );

    }

    inline void
    CppSLHA0::draw_spectrum( std::string const output_eps_filename,
                             std::string const path_to_gnuplot,
                             std::string const path_to_latex,
                             std::string const path_to_dvips,
                             std::string const path_to_ps2eps )
    /* this draws the spectrum using gnuplot, latex, dvips, & ps2eps. the
     * "path_to_X" strings should include the executable name, e.g.
     * /usr/bin/gnuplot for path_to_gnuplot.
     */
    {

      particle_spectrum->draw_spectrum( output_eps_filename,
                                        path_to_gnuplot,
                                        path_to_latex,
                                        path_to_dvips,
                                        path_to_ps2eps );

    }



    inline SPheno_cross_sections_BLOCK*
    CppSLHA2_SPheno::get_SPheno_cross_sections_BLOCK()
    {

      return shortcut_to_cross_sections;

    }

} // end of CppSLHA namespace.

#endif /* CPPSLHA_HPP_ */
