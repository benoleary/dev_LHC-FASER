/*
 * CppSLHA.cpp
 *
 *  Created on: 12 Jan 2010
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2010 Ben O'Leary
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

#include "CppSLHA.hpp"

namespace CppSLHA
{

  /* the constructor requires a file to open, & sets all entries to NaN or not
   * even assigned, in the case of std::vectors of entries & also allocates the
   * memory for the particle spectrum.
   */
  CppSLHA0::CppSLHA0( std::string const given_SLHA_file_name ) :
    SLHA_file_name( given_SLHA_file_name ),
    did_allocate_spectrum_memory( true ),
    // note that this instance allocated memory, so needs to deallocate it
    // with its destructor.
    is_verbose( false )
    // by default, the CppSLHA does not print out all its warnings.
  {

    particle_spectrum = new EW_scale_MSSM_spectrum();
    copy_of_given_spectrum = new EW_scale_spectrum( particle_spectrum );
    // in this case, copy_of_given_spectrum is just set up with an MSSM
    // spectrum.

    BLOCKs.push_back( new MASS_BLOCK( particle_spectrum ) );

  }

  /* this constructor requires a file to open, & this one takes a
   * pre-existing EW_scale_spectrum, which is how one can implement the NMSSM
   * or similar.
   */
  CppSLHA0::CppSLHA0( std::string const given_SLHA_file_name,
                      EW_scale_spectrum* const given_particle_spectrum ) :
    SLHA_file_name( given_SLHA_file_name ),
    did_allocate_spectrum_memory( false )
    // note that this instance did not allocate memory, so
    // should not deallocate it with its destructor.
  {

    if( NULL == given_particle_spectrum )
    /* this constructor prepares an instance that is meant to read in the SLHA
     * file data into a EW_scale_spectrum instance which is held
     * somewhere else in the code. thus it does not want to be involved in the
     * allocation & deallocation of the memory for the
     * EW_scale_spectrum instance, & so it checks that its memory has
     * already been allocated, forcing the programmer to think about memory
     * allocation & deallocation.
     */
      {

        std::cout
        << std::endl
        << "CppSLHA::error! this CppSLHA0 constructor, with given file name \""
        << SLHA_file_name
        << "\" and given pointer to a EW_scale_spectrum instance"
        << " found that the pointer was a NULL pointer.";
        std::cout
        << std::endl
        << "this constructor really doesn't want to do the memory"
        << " management when given such a pointer."
        << " instead, it is taking the easy way out & calling"
        << " exit( EXIT_FAILURE )."
        << " sorry about the inevitable memory leaks...";
        std::cout << std::endl;

        exit( EXIT_FAILURE );

      }
    else
      {

        particle_spectrum = given_particle_spectrum;
        copy_of_given_spectrum
        = new EW_scale_spectrum( given_particle_spectrum );

        BLOCKs.push_back( new MASS_BLOCK( particle_spectrum ) );

      }

  }

  CppSLHA0::~CppSLHA0()
  {

    for( std::vector< SLHA_BLOCK* >::iterator
         deletion_iterator = BLOCKs.begin();
         BLOCKs.end() > deletion_iterator;
         ++deletion_iterator )
      {

        // debugging:
        /*std::cout
        << std::endl
        << "debugging: ~CppSLHA0 deleting "
        << *((*deletion_iterator)->getName());
        std::cout << std::endl;*/
        delete *deletion_iterator;

      }

    // debugging:
    /*std::cout
    << std::endl
    << "debugging: ~CppSLHA0 successfully deleted all its BLOCKs";
    std::cout << std::endl;*/

    // deallocate the memory:
    delete copy_of_given_spectrum;
    if( did_allocate_spectrum_memory )
      {

        // debugging:
        /*std::cout
        << std::endl
        << "debugging: ~CppSLHA0 is about to delete its particle_spectrum.";
        std::cout << std::endl;*/

        delete particle_spectrum;

      }

  }


  void
  CppSLHA0::reset()
  // this resets all the entries apart from the SLHA file name.
  {

    particle_spectrum->become_copy_of( copy_of_given_spectrum );
    // reset particle_spectrum using copy_of_given_spectrum.

    for( std::vector< SLHA_BLOCK* >::iterator
           BLOCK_iterator = BLOCKs.begin();
         BLOCKs.end() > BLOCK_iterator;
         ++BLOCK_iterator )
      {

        (*BLOCK_iterator)->reset();

      }

  }


  SLHA_BLOCK*
  CppSLHA0::get_BLOCK( std::string requested_BLOCK )
  {

    // debugging:
    /*std::cout << std::endl << "debugging: CppSLHA::get_BLOCK( "
                           << requested_BLOCK << " ) called.";
    std::cout << std::endl;*/

    SLHA_BLOCK* return_BLOCK = NULL;

    CppSLHA_global::transform_to_uppercase( &requested_BLOCK );

    // debugging:
    /*std::cout << std::endl << "after transformation: " << requested_BLOCK;
    std::cout << std::endl;*/

    for( std::vector< SLHA_BLOCK* >::iterator
           BLOCK_iterator = BLOCKs.begin();
         BLOCKs.end() > BLOCK_iterator;
         ++BLOCK_iterator )
      {

        // debugging:
        /*std::cout << std::endl << "comparing against: "
                               << (*BLOCK_iterator)->getName();
        std::cout << std::endl;*/

        if( 0 == (*BLOCK_iterator)->get_name()->compare( requested_BLOCK ) )
          {

            return_BLOCK = (*BLOCK_iterator);

          }

      }

    // debugging:
    /*std::cout << std::endl << "returning " << return_BLOCK;
    std::cout << std::endl;*/

    return return_BLOCK;

  }


  SLHA_BLOCK const*
  CppSLHA0::inspect_BLOCK( std::string requested_BLOCK )
  const
  {

    SLHA_BLOCK const* return_BLOCK = NULL;

    CppSLHA_global::transform_to_uppercase( &requested_BLOCK );

    for( std::vector< SLHA_BLOCK* >::const_iterator
           BLOCK_iterator = BLOCKs.begin();
         BLOCKs.end() > BLOCK_iterator;
         ++BLOCK_iterator )
      {

        if( 0 == (*BLOCK_iterator)->get_name()->compare( requested_BLOCK ) )
          {

            return_BLOCK = (*BLOCK_iterator);

          }

      }

    return return_BLOCK;

  }


  std::string
  CppSLHA0::get_MASS_as_string()
  const
  {

    std::stringstream return_string_builder;

    int entries_written = 0;

    return_string_builder << "BLOCK MASS   #" << std::endl;

    for( std::vector< particle_property_set* >::const_iterator
         particle_set_iterator
         = particle_spectrum->get_particle_property_sets()->begin();
         particle_spectrum->get_particle_property_sets()->end()
         > particle_set_iterator;
         ++particle_set_iterator )
    // for each stored particle property set...
      {

        if( (*particle_set_iterator)->has_mass_been_recorded() )
          {

            ++entries_written;

            return_string_builder
            << " "
         << CppSLHA_global::SLHA_int( (*particle_set_iterator)->get_PDG_code(),
                                      9 )
            << "   "
         << CppSLHA_global::SLHA_double( (*particle_set_iterator)->get_mass() )
            << "   # " << *((*particle_set_iterator)->get_name()) << std::endl;

          }

      } // end of for loop over entries.

    if( 1 > entries_written )
      {

        return_string_builder
            << "# nothing recorded from original file for some reason."
            << std::endl;

      }

    return return_string_builder.str();

  }

  std::string
  CppSLHA0::get_DECAY_as_string()
  const
  {

    std::stringstream return_string_builder;

    return_string_builder
    << "# the DECAY section follows:" << std::endl
    << "# PDG = PDG particle ID" << std::endl
    << "# Width = decay width in GeV" << std::endl
    << "# BR = branching ratio" << std::endl
    << "# NDA = number of decay products" << std::endl
    << "# ID1 = PDG particle ID of 1st decay product" << std::endl
    << "# ID2 = PDG particle ID of 2nd decay product" << std::endl
    << "# ID3 = PDG particle ID of 3rd decay product" << std::endl;

    for( std::vector< particle_property_set* >::const_iterator
         particle_set_iterator
         = particle_spectrum->get_particle_property_sets()->begin();
         particle_spectrum->get_particle_property_sets()->end()
         > particle_set_iterator;
         ++particle_set_iterator )
    // for each stored particle property set...
      {

        if( (*particle_set_iterator)->has_decay_width_been_recorded() )
          {

            return_string_builder
            << "#" << std::endl
            << "#         PDG            Width" << std::endl << "DECAY "
         << CppSLHA_global::SLHA_int( (*particle_set_iterator)->get_PDG_code(),
                                      9 ) << "   "
  << CppSLHA_global::SLHA_double( (*particle_set_iterator)->get_decay_width() )
            << "   # " << *((*particle_set_iterator)->get_name()) << " decays"
            << std::endl
            << "#   BR               NDA         ID1       ID2       ID3"
            << std::endl;
            //  xxx+1.23456789E+123   12   123456789 123456789 123456789  #

            if( 0 < (*particle_set_iterator)->get_direct_decays()->size() )
              {

                for( std::list< particle_decay* >::const_iterator
                     direct_decay_iterator
                     = (*particle_set_iterator)->get_direct_decays()->begin();
                     (*particle_set_iterator)->get_direct_decays()->end()
                     != direct_decay_iterator;
                     ++direct_decay_iterator )
                  // for each decay channel of this particle...
                  {

                    return_string_builder
                    << "   "
                    << CppSLHA_global::SLHA_double(
                              (*direct_decay_iterator)->get_branching_ratio() )
                    // write the BR
                    << "   "
                    << CppSLHA_global::SLHA_int(
                         (*direct_decay_iterator)->get_decay_product_PDG_codes(
                                                                     )->size(),
                                                 2 )
                    // & the number of decay products.
                    << "   ";

                    for( std::list< int >::const_iterator
                         decay_product_iterator
                       = (*direct_decay_iterator)->get_decay_product_PDG_codes(
                                                                    )->begin();
                         (*direct_decay_iterator)->get_decay_product_PDG_codes(
                                                                       )->end()
                         != decay_product_iterator;
                         ++decay_product_iterator )
                      // for each decay product...
                      {

                        return_string_builder
                        << CppSLHA_global::SLHA_int( *decay_product_iterator,
                                                     9 )
                        << " ";
                        // print the PDG code.

                      }

                    return_string_builder
                    << "  # BR for " << *((*particle_set_iterator)->get_name())
                    << " ->";
                    // comment the human-readable version.

                    for( std::list< int >::const_iterator
                         decay_product_iterator
                       = (*direct_decay_iterator)->get_decay_product_PDG_codes(
                                                                    )->begin();
                         (*direct_decay_iterator)->get_decay_product_PDG_codes(
                                                                       )->end()
                           != decay_product_iterator;
                         ++decay_product_iterator )
                      // for each decay product...
                      {

                        return_string_builder
                        << " "
                        << *(particle_spectrum->get_PDG_coded_name(
                                                    *decay_product_iterator ));

                      }

                    return_string_builder << std::endl;

                  }

              }
            else
              {

                return_string_builder
                << "# no decays recorded for this particle."
                << std::endl;

              }

          }

      }

    return_string_builder
    << "#" << std::endl << "# end of the DECAY section"
    << std::endl << "#" << std::endl;

    return return_string_builder.str();

  }


  bool
  CppSLHA0::read_file()
  // this reads in the SLHA file and returns true if it was successful.
  {

    // debugging:
    /*std::cout << std::endl
    << "debugging: about to read in a file, but right now the BLOCKs are:";
    std::cout << std::endl;
    for( std::vector< SLHA_BLOCK* >::iterator
           BLOCK_iterator = BLOCKs.begin();
         BLOCKs.end() > BLOCK_iterator;
         ++BLOCK_iterator )
      {
        std::cout << (*BLOCK_iterator)->getName();
        std::cout << std::endl;
      }
    std::cout << std::endl << "that's all of them.";
    std::cout << std::endl;
    */

    reset();

    bool successful_reading = false;

    std::cout
    << std::endl
    << "CppSLHA::read_file() (using \""
    << SLHA_file_name << "\") has been called" << std::endl;

    // debugging:
    /*std::cout
    << std::endl << "debugging: beforehand though, MASS BLOCK now:"
    << std::endl << get_MASS_as_string();
    std::cout << std::endl;
    */

    int const max_length_of_SLHA_file = 10000;
    hash_commented_file_reader input_file_reader( &SLHA_file_name,
                                                  max_length_of_SLHA_file,
                                                  false );
    /* ten thousand lines should be long enough for any SLHA file. if read_file()
     * is reading in more lines than that, it's probably because of some file
     * corruption.
     * also, if it's not enough, change it here.
     */

    std::string input_line_as_string;
    // prepare a string to parse the line.
    std::stringstream input_line_as_stream;
    // prepare a stringstream to parse the line.

    // a set of integers & doubles for holding the data from a line or
    // lines before recording the data, plus a dummy string for ditching
    // stuff.
    int formatting_integer_one;
    int formatting_integer_two;
    int formatting_integer_three;
    int currently_decaying_particle;
    int current_number_of_decay_products;
    double formatting_double;
    std::string formatting_string;

    int current_BLOCK_type = 0;
    /* this is to record what kind of BLOCK is being recorded, since there
     * are 2 special cases, MASS & DECAY.
     * 0: no BLOCK being read
     * 1: MASS being read
     * 2: DECAY being read
     * 3: "normal" BLOCK being read.
     */
    SLHA_BLOCK* current_BLOCK = NULL;
    // this is for sending input lines to the correct BLOCK.

    while( input_file_reader.read_line( &input_line_as_string ) )
      // the evaluation of the conditional reads in the next line.
      {

        // debugging:
        /*std::cout << std::endl << "debugging: read line:";
        std::cout << std::endl << input_line_as_string;
        std::cout << std::endl;
        */

        CppSLHA_global::transform_to_uppercase( &input_line_as_string );
        // convert all the lowercase letters in the line to uppercase.

        input_line_as_stream.clear();
        input_line_as_stream.str( input_line_as_string );

        if( 0 == input_line_as_string.compare( 0,
                                               5,
                                               "BLOCK" ) )
          // if we have found the declaration line of a BLOCK...
          {

            // debugging:
            /*std::cout
            << std::endl
            << "debugging: found BLOCK!";
            std::cout << std::endl;
            */

            input_line_as_stream.ignore( 5 );
            // skip the BLOCK of this stringstream (leading whitespace has already
            // been removed).

            input_line_as_stream >> formatting_string;

            // debugging:
            /*std::cout
            << "after finding \"BLOCK\", formatting_string = "
            << formatting_string;
            std::cout << std::endl;
            */

            if( 0 == formatting_string.compare( "MASS" ) )
              // if the BLOCK is MASS...
              {

                current_BLOCK_type = 1;
                // note that we are reading the MASS BLOCK.

                // debugging:
                /*std::cout
                << std::endl << "debugging: found MASS BLOCK!";
                std::cout << std::endl;
                */

              }
            else
              // otherwise, we are reading a normal BLOCK.
              {

                current_BLOCK_type = 3;
                current_BLOCK = get_BLOCK( formatting_string );
                if( NULL == current_BLOCK )
                  {

                    if( is_verbose )
                      {

                        std::cout
                        << std::endl
                        << "CppSLHA::warning! ignoring unknown BLOCK \""
                        << formatting_string << "\".";

                      }

                    current_BLOCK_type = 0;

                  }
                else
                  {

                    // debugging:
                    /*std::cout
                    << std::endl << "debugging: found "
                    << current_BLOCK->getName() << " BLOCK!" << std::endl
                    << "with pointer " << current_BLOCK;
                    if( !(current_BLOCK->has_recorded()) )
                      {

                        std::cout
                        << " which has not been recorded";

                      }
                    std::cout << std::endl;
                    std::cout
                    << "has_recorded? " << current_BLOCK->has_recorded()
                    << std::endl
                    << "get_extra_BLOCK() = "
                    << current_BLOCK->get_extra_BLOCK();
                    std::cout << std::endl;
                    */

                    int sanity_checker = 0;
                    // sanity_checker is for aborting possible infinite loops.
                    while( ( current_BLOCK->has_recorded() )
                           &&
                           ( NULL != current_BLOCK->get_extra_BLOCK() )
                           &&
                           ( sanity_checker < max_length_of_SLHA_file ) )
                      // if this BLOCK has already recorded & it has an extra
                      // BLOCK already, go to that BLOCK instead.
                      {

                        // debugging:
                        /*std::cout
                        << std::endl << "moving from pointer "
                        << current_BLOCK;
                        */

                        ++sanity_checker;

                        current_BLOCK = current_BLOCK->get_extra_BLOCK();

                        // debugging:
                        /*std::cout
                        << std::endl << " to pointer " << current_BLOCK;
                        std::cout << std::endl;
                        */

                      }
                    if( ( current_BLOCK->has_recorded() )
                        &&
                        ( NULL == current_BLOCK->get_extra_BLOCK() ) )
                      /* if this BLOCK has already recorded but it does not
                       * have an extra BLOCK already, make the extra BLOCK & go
                       * to it.
                       */
                      {

                        current_BLOCK = current_BLOCK->make_extra_BLOCK();

                      }

                    // debugging:
                    /*std::cout
                    << std::endl << "debugging: about to look for Q= for "
                    << current_BLOCK->getName() << " BLOCK!"
                    << std::endl << "with pointer " << current_BLOCK;
                    std::cout << std::endl;
                    */

                    input_line_as_stream >> formatting_string;

                    // if the BLOCK's name is followed by "Q=" or "Q" then
                    // whitespace then "=", then whitespace, note the Q value:
                    if( 0 == formatting_string.compare( "Q=" ) )
                      {

                        input_line_as_stream >> formatting_double;

                        current_BLOCK->set_Q( formatting_double );

                      }
                    else if( 0 == formatting_string.compare( "Q" ) )
                      {

                        input_line_as_stream >> formatting_string;

                        if( 0 == formatting_string.compare( "=" ) )
                          {

                            input_line_as_stream >> formatting_double;

                            current_BLOCK->set_Q( formatting_double );

                          }

                      }

                    current_BLOCK->note_recorded();

                  }

              }

          }
        else if( 0 == input_line_as_string.compare( 0,
                                                    5,
                                                    "DECAY" ) )
          // otherwise, maybe we are about to read in a decay...
          {

            current_BLOCK_type = 2;

            input_line_as_stream.ignore( 5 );
            // skip the DECAY of this stringstream (leading whitespace has
            // already been removed).

            input_line_as_stream >> currently_decaying_particle
                                 >> formatting_double;

            particle_spectrum->record_PDG_coded_decay_width(
                                                   currently_decaying_particle,
                                                         formatting_double );

          }
        else if( 1 == current_BLOCK_type )
          // otherwise, we could be reading the MASS BLOCK...
          {

            input_line_as_stream >> formatting_integer_one
                                 >> formatting_double;

            particle_spectrum->record_PDG_coded_mass( formatting_integer_one,
                                                      formatting_double );
            // record the mass according to its PDG code.

          }
        else if( 3 == current_BLOCK_type )
          // otherwise, we could be reading a normal BLOCK...
          {

            // debugging:
            /*std::cout
            << std::endl << "debugging: reading the following line for "
            << current_BLOCK->getName();
            std::cout << std::endl << input_line_as_stream.str();
            std::cout << std::endl;
            */

            current_BLOCK->read_input( &input_line_as_stream );

          }
        else if( 2 == current_BLOCK_type )
          // or we could be reading a DECAY.
          {

            input_line_as_stream >> formatting_double
                                 >> current_number_of_decay_products;
            // read in as far as knowing how many decay products there are.

            if( 2 == current_number_of_decay_products )
              // if it is a 2-body decay...
              {

                input_line_as_stream >> formatting_integer_one
                                     >> formatting_integer_two;
                // read in the decay products.

                particle_spectrum->record_PDG_coded_direct_decay(
                                                   currently_decaying_particle,
                                                        formatting_integer_one,
                                                        formatting_integer_two,
                                                           formatting_double );
                // record a 2-body decay.

              }
            else if( 3 == current_number_of_decay_products )
              // if it is a 3-body decay...
              {

                input_line_as_stream >> formatting_integer_one
                                     >> formatting_integer_two
                                     >> formatting_integer_three;
                // read in the decay products.

                particle_spectrum->record_PDG_coded_direct_decay(
                                                   currently_decaying_particle,
                                                        formatting_integer_one,
                                                        formatting_integer_two,
                                                      formatting_integer_three,
                                                           formatting_double );
                // record a 3-body decay.

              }
            else
              /* otherwise it was a 1-body-decay (nonsensical in this context) or
               * a 4-or-more-body decay, & coping with this has not been
               * implemented.
               */
              {

                if( is_verbose )
                  {

                    std::cout
                    << std::endl
                    << "CppSLHA::Warning! read_file() has encountered a "
                    << current_number_of_decay_products
                    << "-body decay in the SLHA file \"" << SLHA_file_name
                    << "\" and has not recorded it." << std::endl;

                  }

              }

          } // end of recording a DECAY.

      }  // end of while loop going through all the lines of the SLHA file.

    successful_reading = true;

    // debugging:
    /*std::cout
    << std::endl << "debugging: MASS BLOCK now:" << std::endl
    << get_MASS_as_string();
    std::cout << std::endl;
    */

    return input_file_reader.was_successful();

  }

  bool
  CppSLHA0::read_file( std::string const given_SLHA_file_name )
  // this reads in the specified SLHA file and returns true if it was
  // successful.
  {

    if( is_verbose )
      {

        std::cout
        << std::endl
        << "CppSLHA::read_file( \"" << given_SLHA_file_name << "\" ) called";
        std::cout
        << " - warning, all previously held data in this instance has been"
        << " discarded.";
        std::cout << std::endl;

      }

    SLHA_file_name.assign( given_SLHA_file_name );

    return read_file();

  }


  bool
  CppSLHA0::write_file( std::string const output_SLHA_file_name )
  // this writes out all the BLOCKs to a file with the given name (DECAY
  // counts as a BLOCK name in this case).
  {

    bool return_bool = false;

    std::ofstream* output_file = new std::ofstream();

    output_file->open( output_SLHA_file_name.c_str() );

    if( output_file->is_open() )
      {

        return_bool = true;

        *output_file
        << "# this was generated automatically by a program using BOL's"
        << " CppSLHA #" << std::endl;

        for( std::vector< SLHA_BLOCK* >::iterator
               BLOCK_iterator = BLOCKs.begin();
             BLOCKs.end() > BLOCK_iterator;
             ++BLOCK_iterator )
          {
            *output_file << "#-------------------------------------------#"
                << std::endl << *((*BLOCK_iterator)->get_as_string());
          }
        //*output_file << "#-------------------------------------------#"
        //    << std::endl << get_MASS_as_string();
        *output_file << "#-------------------------------------------#"
            << std::endl << get_DECAY_as_string();
        *output_file << "#-------------------------------------------#"
            << std::endl;

        *output_file << std::endl;
        // put an extra newline at the end, just for kicks.

        // now close the file.
        output_file->close();

      }

    // clean up.
    delete output_file;
    output_file = NULL;  // just to be safe.

    return return_bool;

  }

  bool
  CppSLHA0::write_file( std::string const output_SLHA_file_name,
                        std::list< std::string > const* const list_of_BLOCKS )
  // this writes out the BLOCKs from the supplied list to a file with the
  // given name (DECAY counts as a BLOCK name in this case).
  {

    bool return_bool = false;

    // use the list of strings to make a vector of pointers to BLOCKs:
    std::vector< SLHA_BLOCK* > BLOCKs_to_write;
    for( std::list< std::string >::const_iterator
           BLOCK_iterator = list_of_BLOCKS->begin();
         list_of_BLOCKS->end() != BLOCK_iterator;
         ++BLOCK_iterator )
      {

        BLOCKs_to_write.push_back( get_BLOCK(*BLOCK_iterator) );

      }

    std::ofstream* output_file = new std::ofstream();

    output_file->open( output_SLHA_file_name.c_str() );

    if( output_file->is_open() )
      {

        return_bool = true;

        for( std::vector< SLHA_BLOCK* >::iterator
               BLOCK_iterator = BLOCKs_to_write.begin();
             BLOCKs_to_write.end() > BLOCK_iterator;
             ++BLOCK_iterator )
        // go through the supplied BLOCKs:
          {

            *output_file << *((*BLOCK_iterator)->get_as_string());
            // write each BLOCK.

          }

        *output_file << std::endl;
        // put an extra newline at the end, just for kicks.

      }

    // now close the file:
    output_file->close();

    // clean up.
    delete output_file;
    output_file = NULL;  // just to be safe.

    return return_bool;

  }


  // STUFF FOR TESTING PURPOSES!

  int
  CppSLHA0::write_particle_spectrum(
                                  double const minimum_branching_ratio_to_keep,
                                     std::string const output_file_name )
  // this writes all the stored spectrum data to a file with the provided
  // name.  it returns the number of particle property sets written.
  {

    std::cout << std::endl
        << "CppSLHA::write_particle_spectrum() has been called" << std::endl;

    std::cout
    << std::endl << "found "
    << particle_spectrum->find_all_cascade_decays(
                                              minimum_branching_ratio_to_keep )
    << " cascade decays" << std::endl;

    std::ofstream* output_file = new std::ofstream();

    output_file->open( output_file_name.c_str() );

    *output_file << "# this is a file with the data from a CppSLHA instance #"
                 << std::endl;

    int property_sets_written = 0;

    for( std::vector< particle_property_set* >::const_iterator
           particle_set_iterator
             = particle_spectrum->get_particle_property_sets()->begin();
         particle_spectrum->get_particle_property_sets()->end()
           > particle_set_iterator;
         ++particle_set_iterator )
    // for each stored particle property set...
      {

        ++property_sets_written;
        // increment the number of property sets written.

        // write out the properties:

        *output_file
        << std::endl
        << " stored particle " << property_sets_written
        << ": PDG code " << (*particle_set_iterator)->get_PDG_code()
        << ", name " << *((*particle_set_iterator)->get_name())
        << ", mass " << (*particle_set_iterator)->get_mass()
        << ", decay width " << (*particle_set_iterator)->get_decay_width();

        if( (*particle_set_iterator)->counts_as_stable_invisible() )
          {

            *output_file << ", counts as stable and invisible";

          }
        else
          {

            *output_file << ", does not count as stable and invisible";

          }
        if( (*particle_set_iterator)->counts_as_jet() )
          {

            *output_file << ", counts as a jet";

          }
        else
          {

            *output_file << ", does not count as a jet";

          }
        if( (*particle_set_iterator)->counts_as_light_lepton() )
          {

            *output_file << ", counts as a light lepton";

          }
        else
          {

            *output_file << ", does not count as a light lepton";

          }
        if( (*particle_set_iterator)->counts_as_self_conjugate() )
          {

            *output_file << ", counts as self-conjugate";

          }
        else
          {

            *output_file << ", does not count as self-conjugate";

          }

        *output_file
        << std::endl
        << " direct decays ("
        << (*particle_set_iterator)->get_direct_decays()->size() << "):"
        << std::endl;

        for( std::list< particle_decay* >::const_iterator
               direct_decay_iterator
                 = (*particle_set_iterator)->get_direct_decays()->begin();
             (*particle_set_iterator)->get_direct_decays()->end()
               != direct_decay_iterator;
             ++direct_decay_iterator )
          {

            *output_file
            << "("
            << (*direct_decay_iterator)->get_decay_product_PDG_codes()->size()
            << "-body)";

            for( std::list< int >::const_iterator
                   decay_product_iterator
                     = (*direct_decay_iterator)->get_decay_product_PDG_codes(
                                                   )->begin();
                 (*direct_decay_iterator)->get_decay_product_PDG_codes(
                                             )->end()
                   != decay_product_iterator;
                 ++decay_product_iterator )
              {

                *output_file
                << " " << *decay_product_iterator << " ("
                << *(particle_spectrum->get_PDG_coded_name(
                                                     *decay_product_iterator ))
                << ")";

              }

            *output_file
            << " BR:" << (*direct_decay_iterator)->get_branching_ratio()
            << std::endl;

          }

        *output_file
        << std::endl
        << " cascade decays ("
        << (*particle_set_iterator)->get_cascade_decays()->size() << "):"
        << std::endl;

        for( std::list< particle_decay* >::const_iterator
               cascade_decay_iterator
                 = (*particle_set_iterator)->get_cascade_decays()->begin();
             (*particle_set_iterator)->get_cascade_decays()->end()
               != cascade_decay_iterator;
             ++cascade_decay_iterator )
          {

            *output_file
            << "("
            << (*cascade_decay_iterator)->get_decay_product_PDG_codes(
                                            )->size()
            << "-body)";

            for( std::list< int >::const_iterator
                   decay_product_iterator
                   = (*cascade_decay_iterator)->get_decay_product_PDG_codes(
                                                  )->begin();
                 (*cascade_decay_iterator)->get_decay_product_PDG_codes(
                                              )->end()
                   != decay_product_iterator;
                 ++decay_product_iterator )
              {

                *output_file
                << " " << *decay_product_iterator << " ("
                << *(particle_spectrum->get_PDG_coded_name(
                                                     *decay_product_iterator ))
                << ")";

              }

            *output_file
            << " BR:" << (*cascade_decay_iterator)->get_branching_ratio()
            << std::endl;

          }

      }

    std::cout
    << std::endl
    << "finished spectrum data, now writing the extra stuff.";
    std::cout << std::endl;

    // debugging:
    /*std::cout << std::endl << "debugging:";
     std::cout << std::endl;
     */

    *output_file << std::endl << std::endl;
    for( std::vector< SLHA_BLOCK* >::const_iterator
           BLOCK_iterator = BLOCKs.begin();
         BLOCKs.end() > BLOCK_iterator;
         ++BLOCK_iterator )
      {
        // debugging:
        /*std::cout << "should be \"" << (*BLOCK_iterator)->getName() << "\"";
         std::cout << std::endl;
         std::cout
         << "#-------------------------------------------#"
         << std::endl << (*BLOCK_iterator)->get_as_string();
         */
        *output_file
        << "#-------------------------------------------#"
        << std::endl << *((*BLOCK_iterator)->get_as_string());
      }
    *output_file
    << "#-------------------------------------------#"
    << std::endl << get_MASS_as_string();
    *output_file
    << "#-------------------------------------------#"
    << std::endl << get_DECAY_as_string();
    *output_file
    << "#-------------------------------------------#"
    << std::endl;

    // now close the file.
    output_file->close();

    // clean up.
    delete output_file;
    output_file = NULL;  // just to be safe.

    return property_sets_written;

  }
  // END OF TESTING FUNCTION!


  /* the derived classes build off each other. I have tried to keep it all
   * modular. CppSLHA1 builds off CppSLHA directly. CppSLHA2 uses the CppSLHA1
   * constructor to add the BLOCKs from SLHA1 & then adds in the extra SLHA2
   * BLOCKs. User-programmed extension should be easy, following these
   * examples & the others below.
   */

  CppSLHA1::CppSLHA1( std::string const given_SLHA_file_name ) :
    CppSLHA0( given_SLHA_file_name )
  // the constructor requires a file to open.
  {

    initialize();

  }

  CppSLHA1::CppSLHA1( std::string const given_SLHA_file_name,
                      EW_scale_spectrum* const given_particle_spectrum ) :
    CppSLHA0( given_SLHA_file_name,
              given_particle_spectrum )
  /* this constructor requires a file to open, & this one takes a pre-existing
   * EW_scale_spectrum, which is how one can implement the NMSSM or
   * similar.
   */
  {

    initialize();

  }

  CppSLHA1::~CppSLHA1()
  {

    // does nothing.

  }


  void
  CppSLHA1::initialize()
  {

    // debugging:
    /*std::cout << std::endl << "debugging: CppSLHA1::initialize() called";
     std::cout << std::endl;
     */

    SLHA_BLOCK*
    BLOCK_adder
    = new sparse_one_index_BLOCK( "MODSEL",
                                  "model selection" );
    BLOCK_adder->write_entry_as_integer( true );
    // SPheno at least requires these entries to be printed as integers.
    BLOCKs.push_back( BLOCK_adder );
    BLOCK_adder->set_comment( 1,
                              "choice of SUSY breaking model" );
    BLOCK_adder->set_comment( 3,
                              "choice of particle content" );
    BLOCK_adder->set_comment( 11,
                              "number of BLOCK copies with log-spaced Q" );
    BLOCK_adder->set_comment( 12,
                              "maximum Q" );
    BLOCK_adder->set_comment( 21,
                   "PDG code of particle with mass to be taken as Q" );

    BLOCK_adder = new one_index_BLOCK( "SMINPUTS",
                                       "standard model inputs" );
    BLOCKs.push_back( BLOCK_adder );
    BLOCK_adder->set_comment( 1,
                   "1/alpha_EM at the Z mass in the MSbar scheme" );
    BLOCK_adder->set_comment( 2,
                              "the Fermi constant G_F in GeV^(-2)" );
    BLOCK_adder->set_comment( 3,
                              "1/alpha_s at the Z mass in the MSbar scheme" );
    BLOCK_adder->set_comment( 4,
                              "the pole mass of the Z boson" );
    BLOCK_adder->set_comment( 5,
                   "the bottom quark mass at itself in the MSbar scheme" );
    BLOCK_adder->set_comment( 6,
                              "the pole mass of the top quark" );
    BLOCK_adder->set_comment( 7,
                              "the pole mass of the tau lepton" );

    BLOCK_adder = new sparse_one_index_BLOCK( "MINPAR",
                        "minimal set of model input parameters" );
    BLOCKs.push_back( BLOCK_adder );
    BLOCK_adder->set_comment( 1,
      "in the mSUGRA-inspired CMSSM, the common scalar soft mass" );
    BLOCK_adder->set_comment( 2,
      "in the mSUGRA-inspired CMSSM, the common gaugino soft mass" );
    BLOCK_adder->set_comment( 3,
      "in the mSUGRA-inspired CMSSM, the tangent of beta at the Z mass" );
    BLOCK_adder->set_comment( 4,
      "in the mSUGRA-inspired CMSSM, the sign of the mu parameter" );
    BLOCK_adder->set_comment( 5,
      "in the mSUGRA-inspired CMSSM, the common scalar soft trilinear mass" );

    BLOCK_adder = new sparse_one_index_BLOCK( "EXTPAR",
                                              "extra model input parameters" );
    BLOCKs.push_back( BLOCK_adder );
    BLOCK_adder->set_comment( 0,
                              "the scale for these parameters" );
    BLOCK_adder->set_comment( 1,
                              "M_1" );
    BLOCK_adder->set_comment( 2,
                              "M_2" );
    BLOCK_adder->set_comment( 3,
                              "M_3" );
    BLOCK_adder->set_comment( 21,
                              "M^2_(H,d)" );
    BLOCK_adder->set_comment( 22,
                              "M^2_(H,u)" );
    BLOCK_adder->set_comment( 31,
                              "M_(L,11)" );
    BLOCK_adder->set_comment( 32,
                              "M_(L,22)" );
    BLOCK_adder->set_comment( 33,
                              "M_(L,33)" );
    BLOCK_adder->set_comment( 34,
                              "M_(E,11)" );
    BLOCK_adder->set_comment( 35,
                              "M_(E,22)" );
    BLOCK_adder->set_comment( 36,
                              "M_(E,33)" );
    BLOCK_adder->set_comment( 41,
                              "M_(Q,11)" );
    BLOCK_adder->set_comment( 42,
                              "M_(Q,22)" );
    BLOCK_adder->set_comment( 43,
                              "M_(Q,33)" );
    BLOCK_adder->set_comment( 44,
                              "M_(U,11)" );
    BLOCK_adder->set_comment( 45,
                              "M_(U,22)" );
    BLOCK_adder->set_comment( 46,
                              "M_(U,33)" );
    BLOCK_adder->set_comment( 47,
                              "M_(D,11)" );
    BLOCK_adder->set_comment( 48,
                              "M_(D,22)" );
    BLOCK_adder->set_comment( 49,
                              "M_(D,33)" );
    BLOCK_adder->set_comment( 51,
                              "N_1" );
    BLOCK_adder->set_comment( 52,
                              "N_2" );
    BLOCK_adder->set_comment( 53,
                              "N_3" );

    BLOCK_adder = new two_index_BLOCK( "NMIX",
                                       "neutralino mixing matrix N" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "UMIX",
                                       "chargino mixing matrix U" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "VMIX",
                                       "chargino mixing matrix V" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "STOPMIX",
                                       "stop mixing matrix R_stop" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "SBOTMIX",
                                       "sbottom mixing matrix R_sbottom" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "STAUMIX",
                                       "stau mixing matrix R_stau" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new zero_index_BLOCK( "ALPHA",
                                        "effective Higgs mixing angle" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new sparse_one_index_BLOCK( "HMIX",
                                              "Higgs mixing parameters" );
    BLOCKs.push_back( BLOCK_adder );
    BLOCK_adder->set_comment( 1,
                              "mu" );
    BLOCK_adder->set_comment( 2,
                              "tangent of beta at scale Q" );
    BLOCK_adder->set_comment( 3,
                              "v(Q)" );
    BLOCK_adder->set_comment( 4,
                              "m^2_A(Q)" );

    BLOCK_adder = new one_index_BLOCK( "GAUGE",
                                       "gauge couplings" );
    BLOCKs.push_back( BLOCK_adder );
    BLOCK_adder->set_comment( 1,
                              "g'(Q) in the DRbar scheme" );
    BLOCK_adder->set_comment( 2,
                              "g(Q) in the DRbar scheme" );
    BLOCK_adder->set_comment( 3,
                              "g_3(Q) in the DRbar scheme" );

    BLOCK_adder = new sparse_one_index_BLOCK( "MSOFT",
                         "soft mass parameters at scale Q" );
    BLOCKs.push_back( BLOCK_adder );
    BLOCK_adder->set_comment( 1,
                              "M_1" );
    BLOCK_adder->set_comment( 2,
                              "M_2" );
    BLOCK_adder->set_comment( 3,
                              "M_3" );
    BLOCK_adder->set_comment( 21,
                              "M^2_(H,d)" );
    BLOCK_adder->set_comment( 22,
                              "M^2_(H,u)" );
    BLOCK_adder->set_comment( 31,
                              "M_(L,11)" );
    BLOCK_adder->set_comment( 32,
                              "M_(L,22)" );
    BLOCK_adder->set_comment( 33,
                              "M_(L,33)" );
    BLOCK_adder->set_comment( 34,
                              "M_(E,11)" );
    BLOCK_adder->set_comment( 35,
                              "M_(E,22)" );
    BLOCK_adder->set_comment( 36,
                              "M_(E,33)" );
    BLOCK_adder->set_comment( 41,
                              "M_(Q,11)" );
    BLOCK_adder->set_comment( 42,
                              "M_(Q,22)" );
    BLOCK_adder->set_comment( 43,
                              "M_(Q,33)" );
    BLOCK_adder->set_comment( 44,
                              "M_(U,11)" );
    BLOCK_adder->set_comment( 45,
                              "M_(U,22)" );
    BLOCK_adder->set_comment( 46,
                              "M_(U,33)" );
    BLOCK_adder->set_comment( 47,
                              "M_(D,11)" );
    BLOCK_adder->set_comment( 48,
                              "M_(D,22)" );
    BLOCK_adder->set_comment( 49,
                              "M_(D,33)" );

    BLOCK_adder = new two_index_BLOCK( "AU",
                        "soft trilinear parameters for the up-type squarks" );
    BLOCKs.push_back( BLOCK_adder );
    BLOCK_adder->set_comment( 1,
                              1,
                              "A_u at scale Q in the DRbar scheme" );
    BLOCK_adder->set_comment( 2,
                              2,
                              "A_c at scale Q in the DRbar scheme" );
    BLOCK_adder->set_comment( 3,
                              3,
                              "A_t at scale Q in the DRbar scheme" );

    BLOCK_adder = new two_index_BLOCK( "AD",
                      "soft trilinear parameters for the down-type squarks" );
    BLOCKs.push_back( BLOCK_adder );
    BLOCK_adder->set_comment( 1,
                              1,
                              "A_d at scale Q in the DRbar scheme" );
    BLOCK_adder->set_comment( 2,
                              2,
                              "A_s at scale Q in the DRbar scheme" );
    BLOCK_adder->set_comment( 3,
                              3,
                              "A_b at scale Q in the DRbar scheme" );

    BLOCK_adder = new two_index_BLOCK( "AE",
                      "soft trilinear parameters for the charged sleptons" );
    BLOCKs.push_back( BLOCK_adder );
    BLOCK_adder->set_comment( 1,
                              1,
                              "A_e at scale Q in the DRbar scheme" );
    BLOCK_adder->set_comment( 2,
                              2,
                              "A_mu at scale Q in the DRbar scheme" );
    BLOCK_adder->set_comment( 3,
                              3,
                              "A_tau at scale Q in the DRbar scheme" );

    BLOCK_adder = new two_index_BLOCK( "YU",
                        "Yukawa couplings for the up-type quarks" );
    BLOCKs.push_back( BLOCK_adder );
    BLOCK_adder->set_comment( 1,
                              1,
                              "Y_u at scale Q in the DRbar scheme" );
    BLOCK_adder->set_comment( 2,
                              2,
                              "Y_c at scale Q in the DRbar scheme" );
    BLOCK_adder->set_comment( 3,
                              3,
                              "Y_t at scale Q in the DRbar scheme" );

    BLOCK_adder = new two_index_BLOCK( "YD",
                        "Yukawa couplings for the down-type quarks" );
    BLOCKs.push_back( BLOCK_adder );
    BLOCK_adder->set_comment( 1,
                              1,
                              "Y_d at scale Q in the DRbar scheme" );
    BLOCK_adder->set_comment( 2,
                              2,
                              "Y_s at scale Q in the DRbar scheme" );
    BLOCK_adder->set_comment( 3,
                              3,
                              "Y_b at scale Q in the DRbar scheme" );

    BLOCK_adder = new two_index_BLOCK( "YE",
                        "Yukawa couplings for the charged leptons" );
    BLOCKs.push_back( BLOCK_adder );
    BLOCK_adder->set_comment( 1,
                              1,
                              "Y_e at scale Q in the DRbar scheme" );
    BLOCK_adder->set_comment( 2,
                              2,
                              "Y_mu at scale Q in the DRbar scheme" );
    BLOCK_adder->set_comment( 3,
                              3,
                              "Y_tau at scale Q in the DRbar scheme" );

    BLOCK_adder
    = new multiple_same_one_index_string_BLOCK( "SPINFO",
                                                "spectrum generator info" );
    BLOCKs.push_back( BLOCK_adder );
    BLOCK_adder->set_comment( 1,
                              1,
                              "spectrum calculator name" );
    BLOCK_adder->set_comment( 2,
                              1,
                              "spectrum calculator version" );
    BLOCK_adder->set_comment( 3,
                              1,
                              "warning message" );
    BLOCK_adder->set_comment( 4,
                              1,
                              "error message" );

  }



  CppSLHA2::CppSLHA2( std::string const given_SLHA_file_name ) :
    CppSLHA1( given_SLHA_file_name )
  // the constructor requires a file to open.
  {

    initialize();

  }

  CppSLHA2::CppSLHA2( std::string const given_SLHA_file_name,
                      EW_scale_spectrum* const given_particle_spectrum ) :
    CppSLHA1( given_SLHA_file_name,
              given_particle_spectrum )
  /* this constructor requires a file to open, & this one takes a pre-existing
   * EW_scale_spectrum, which is how one can implement the NMSSM or similar.
   */
  {

    initialize();

  }

  CppSLHA2::~CppSLHA2()
  {

    // does nothing.

  }


  void
  CppSLHA2::initialize()
  {

    // the BLOCKs from SLHA1 have been added by the CppSLHA1 constructor part
    // of the CppSLHA2 constructor.

    // comments for extra parts of the SLHA1 BLOCKs:

    SLHA_BLOCK* BLOCK_adder = get_BLOCK( "MODSEL" );
    BLOCK_adder->write_entry_as_integer( true );
    // SPheno at least requires these entries to be printed as integers.
    if( NULL != BLOCK_adder )
      {

        BLOCK_adder->set_comment( 4,
                                  "choice of R-parity conservation" );
        BLOCK_adder->set_comment( 5,
                                  "choice of CP conservation" );
        BLOCK_adder->set_comment( 6,
                                  "choice of flavor conservation" );

      }
    else
      {

        if( is_verbose )
          {

            std::cout
            << std::endl
            << "CppSLHA::warning! CppSLHA2::initialize() could not add"
            << " comments to MODSEL because the BLOCK does not (yet) exist!";
            std::cout << std::endl;

          }

      }
    BLOCK_adder = get_BLOCK( "SMINPUTS" );
    if( NULL != BLOCK_adder )
      {

        BLOCK_adder->set_comment( 8,
                     "the pole mass of the 3rd-lightest (heaviest) neutrino" );
        BLOCK_adder->set_comment( 11,
                                  "the pole mass of the electron" );
        BLOCK_adder->set_comment( 12,
                                  "the pole mass of the lightest neutrino" );
        BLOCK_adder->set_comment( 13,
                                  "the pole mass of the muon" );
        BLOCK_adder->set_comment( 14,
                       "the pole mass of the 2nd-lightest neutrino'" );
        BLOCK_adder->set_comment( 21,
                       "the down quark mass at 2 GeV in the MSbar scheme" );
        BLOCK_adder->set_comment( 22,
                       "the up quark mass at 2 GeV in the MSbar scheme" );
        BLOCK_adder->set_comment( 23,
                       "the strange quark mass at 2 GeV in the MSbar scheme" );
        BLOCK_adder->set_comment( 24,
                       "the charm quark mass at itself in the MSbar scheme" );

      }
    else
      {

        if( is_verbose )
          {

            std::cout
            << std::endl
            << "CppSLHA::warning! CppSLHA2::initialize() could not add"
            << " comments to SMINPUTS because the BLOCK does not (yet) exist!";
            std::cout << std::endl;

          }

      }


    // extra non-model-specific BLOCKs:

    BLOCK_adder = new two_index_BLOCK( "QEXTPAR",
            "the input scales for the EXTPAR parameters (index for index)" );
    BLOCKs.push_back( BLOCK_adder );


    // Charge-Parity violation:
    // imaginary parts of the SLHA1 BLOCKs & extra scalar mixing BLOCKs:

    BLOCK_adder = get_BLOCK( "MINPAR" );
    if( NULL != BLOCK_adder )
      {

        BLOCK_adder->set_comment( 4,
 "in the mSUGRA-inspired CMSSM, the cosine of the phase of the mu parameter" );

      }
    else
      {

        if( is_verbose )
          {

            std::cout
            << std::endl
            << "CppSLHA::warning! CppSLHA2::initialize() could not add"
            << " comments to MINPAR because the BLOCK does not (yet) exist!";
            std::cout << std::endl;

          }

      }

    BLOCK_adder = new sparse_one_index_BLOCK( "IMMINPAR",
                                              "imaginary parts of MINPAR" );
    BLOCKs.push_back( BLOCK_adder );
    BLOCK_adder->set_comment( 4,
  "in the mSUGRA-inspired CMSSM, the sine of the phase of the mu parameter" );

    BLOCK_adder = new sparse_one_index_BLOCK( "IMEXTPAR",
                                              "imaginary parts of EXTPAR" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMNMIX",
                                       "imaginary parts of NMIX" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMUMIX",
                                       "imaginary parts of UMIX" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMVMIX",
                                       "imaginary parts of VMIX" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMSTOPMIX",
                                       "imaginary parts of STOPMIX" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMSBOTMIX",
                                       "imaginary parts of SBOTMIX" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMSTAUMIX",
                                       "imaginary parts of STAUMIX" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new zero_index_BLOCK( "IMALPHA",
                                        "imaginary part of ALPHA" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new sparse_one_index_BLOCK( "IMHMIX",
                                              "imaginary parts of HMIX" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new one_index_BLOCK( "IMGAUGE",
                                       "imaginary parts of GAUGE" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new sparse_one_index_BLOCK( "IMMSOFT",
                                              "imaginary parts of MSOFT" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMAU",
                                       "imaginary parts of AU" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMAD",
                                       "imaginary parts of AD" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMAE",
                                       "imaginary parts of AE" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMYU",
                                       "imaginary parts of YU" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMYD",
                                       "imaginary parts of YD" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMYE",
                                       "imaginary parts of YE" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "CVHMIX",
                                       "neutral scalar mixing matrix" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMCVHMIX",
                                       "imaginary parts of CVHMIX" );
    BLOCKs.push_back( BLOCK_adder );


    // flavor:
    // flavor-mixing matrices

    BLOCK_adder = new one_index_BLOCK( "VCKMIN",
  "input-scale CKM quark mixing matrix in the Wolfenstein parameterization" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new one_index_BLOCK( "IMVCKMIN",
                                       "imaginary parts of VCKMIN" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new one_index_BLOCK( "VCKM",
    "CKM quark mixing matrix in the Wolfenstein parameterization at scale Q" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new one_index_BLOCK( "IMVCKM",
                                       "imaginary parts of VCKM" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new one_index_BLOCK( "UPMNSIN",
    "input-scale UPMNS neutrino mixing matrix in the PDG parameterization" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new one_index_BLOCK( "IMUPMNSIN",
                                       "imaginary parts of UPMNSIN" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new one_index_BLOCK( "UPMNS",
    "UPMNS neutrino mixing matrix in the PDG parameterization at scale Q" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new one_index_BLOCK( "IMUPMNS",
                                       "imaginary parts of UPMNS" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "MSQ2IN",
                        "input-scale m^2_(Q,ij) in the super-CKM basis" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMMSQ2IN",
                                       "imaginary parts of MSQ2IN" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "MSQ2",
                        "m^2_(Q,ij)(Q) in the super-CKM basis" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMMSQ2",
                                       "imaginary parts of MSQ2" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "MSU2IN",
                        "input-scale m^2_(U,ij) in the super-CKM basis" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMMSU2IN",
                                       "imaginary parts of MSU2IN" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "MSU2",
                        "m^2_(U,ij)(Q) in the super-CKM basis" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMMSU2",
                                       "imaginary parts of MSU2" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "MSD2IN",
                        "input-scale m^2_(D,ij) in the super-CKM basis" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMMSD2IN",
                                       "imaginary parts of MSD2IN" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "MSD2",
                        "m^2_(D,ij)(Q) in the super-CKM basis" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMMSD2",
                                       "imaginary parts of MSD2" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "MSL2IN",
                        "input-scale m^2_(L,ij) in the super-PMNS basis" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMMSL2IN",
                                       "imaginary parts of MSL2IN" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "MSL2",
                        "m^2_(L,ij)(Q) in the super-PMNS basis" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMMSL2",
                                       "imaginary parts of MSL2" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "MSE2IN",
                        "input-scale m^2_(E,ij) in the super-PMNS basis" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMMSE2IN",
                                       "imaginary parts of MSE2IN" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "MSE2",
                        "m^2_(E,ij)(Q) in the super-PMNS basis" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMMSE2",
                                       "imaginary parts of MSE2" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "MSN2IN",
                        "input-scale m^2_(N,ij) in the super-PMNS basis" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMMSN2IN",
                                       "imaginary parts of MSN2IN" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "MSN2",
                        "m^2_(N,ij)(Q) in the super-PMNS basis" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMMSN2",
                                       "imaginary parts of MSN2" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "TUIN",
                        "input-scale T_(U,ij) in the super-CKM basis" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMTUIN",
                                       "imaginary parts of TUIN" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "TU",
                        "T_(U,ij)(Q) in the super-CKM basis" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMTU",
                                       "imaginary parts of TU" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "TDIN",
                        "input-scale T_(D,ij) in the super-CKM basis" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMTDIN",
                                       "imaginary parts of TDIN" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "TD",
                        "T_(D,ij)(Q) in the super-CKM basis" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMTD",
                                       "imaginary parts of TD" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "TEIN",
                        "input-scale T_(E,ij) in the super-PMNS basis" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMTEIN",
                                       "imaginary parts of TEIN" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "TE",
                        "T_(E,ij)(Q) in the super-PMNS basis" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMTE",
                                       "imaginary parts of TE" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "TNIN",
                        "input-scale T_(N,ij) in the super-PMNS basis" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMTNIN",
                                       "imaginary parts of TNIN" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "TN",
                        "T_(N,ij)(Q) in the super-PMNS basis" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMTN",
                                       "imaginary parts of TN" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "YN",
                                       "Yukawa couplings for the neutrinos" );
    BLOCKs.push_back( BLOCK_adder );
    BLOCK_adder->set_comment( 1,
                              1,
                              "Y_nu_1 at scale Q in the DRbar scheme" );
    BLOCK_adder->set_comment( 2,
                              2,
                              "Y_nu_2 at scale Q in the DRbar scheme" );
    BLOCK_adder->set_comment( 3,
                              3,
                              "Y_nu_3 at scale Q in the DRbar scheme" );

    BLOCK_adder = new two_index_BLOCK( "USQMIX",
                                       "up-type squark mixing matrix R_u" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMUSQMIX",
                                       "imaginary parts of USQMIX" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "DSQMIX",
                                       "down-type squark mixing matrix R_d" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMDSQMIX",
                                       "imaginary parts of DSQMIX" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "SELMIX",
                                       "charged slepton mixing matrix R_e" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMSELMIX",
                                       "imaginary parts of SELMIX" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "SNUMIX",
                                       "sneutrino mixing matrix R_nu" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMSNUMIX",
                                       "imaginary parts of SNUMIX" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "SNSMIX",
                                       "scalar sneutrino mixing matrix" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMSNSMIX",
                                       "imaginary parts of SNSMIX" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "SNAMIX",
                        "pseudoscalar sneutrino mixing matrix" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMSNAMIX",
                                       "imaginary parts of SNAMIX" );
    BLOCKs.push_back( BLOCK_adder );


    // R-parity violation:
    // trilinear & bilinear couplings & soft SUSY-breaking terms:

    BLOCK_adder = new three_index_BLOCK( "RVLAMLLEIN",
                      "input-scale R-parity violating couplings lambda_ijk" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new three_index_BLOCK( "IMRVLAMLLEIN",
                                         "imaginary parts of RVLAMLLEIN" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new three_index_BLOCK( "RVLAMLLE",
                      "R-parity violating couplings lambda_ijk at scale Q" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new three_index_BLOCK( "IMRVLAMLLE",
                                         "imaginary parts of RVLAMLLE" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new three_index_BLOCK( "RVLAMLQDIN",
                      "input-scale R-parity violating couplings lambda'_ijk" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new three_index_BLOCK( "IMRVLAMLQDIN",
                                         "imaginary parts of RVLAMLQDIN" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new three_index_BLOCK( "RVLAMLQD",
                      "R-parity violating couplings lambda'_ijk at scale Q" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new three_index_BLOCK( "IMRVLAMLQD",
                                         "imaginary parts of RVLAMLQD" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new three_index_BLOCK( "RVLAMUDDIN",
                  "input-scale R-parity violating couplings lambda''_ijk" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new three_index_BLOCK( "IMRVLAMUDDIN",
                                         "imaginary parts of RVLAMUDDIN" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new three_index_BLOCK( "RVLAMUDD",
                      "R-parity violating couplings lambda''_ijk at scale Q" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new three_index_BLOCK( "IMRVLAMUDD",
                                         "imaginary parts of RVLAMLQD" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new three_index_BLOCK( "RVTLLEIN",
      "input-scale R-parity violating trilinear soft terms T_ijk" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new three_index_BLOCK( "IMRVTLLEIN",
                                         "imaginary parts of RVTLLEIN" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new three_index_BLOCK( "RVTLLE",
      "R-parity violating trilinear soft terms T_ijk at scale Q" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new three_index_BLOCK( "IMRVTLLE",
                                         "imaginary parts of RVTLLE" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new three_index_BLOCK( "RVTLQDIN",
      "input-scale R-parity violating  trilinear soft terms T'_ijk" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new three_index_BLOCK( "IMRVTLQDIN",
                                         "imaginary parts of RVTLQDIN" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new three_index_BLOCK( "RVTLQD",
      "R-parity violating trilinear soft terms T'_ijk at scale Q" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new three_index_BLOCK( "IMRVTLQD",
                                         "imaginary parts of RVTLQD" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new three_index_BLOCK( "RVTUDDIN",
      "input-scale R-parity violating trilinear soft terms T''_ijk" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new three_index_BLOCK( "IMRVTUDDIN",
                                         "imaginary parts of RVTUDDIN" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new three_index_BLOCK( "RVTUDD",
      "R-parity violating trilinear soft terms T''_ijk at scale Q" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new three_index_BLOCK( "IMRVTUDD",
                                         "imaginary parts of RVTLQD" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new one_index_BLOCK( "RVKAPPAIN",
                        "input-scale R-parity violating couplings kappa_i" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new one_index_BLOCK( "IMRVKAPPAIN",
                                       "imaginary parts of IMRVKAPPAIN" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new one_index_BLOCK( "RVKAPPA",
                        "R-parity violating couplings kappa_i at scale Q" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new one_index_BLOCK( "IMRVKAPPA",
                                       "imaginary parts of IMRVKAPPA" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new one_index_BLOCK( "RVDIN",
                        "input-scale R-parity violating couplings D_i" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new one_index_BLOCK( "IMRVDIN",
                                       "imaginary parts of IMRVDIN" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new one_index_BLOCK( "RVD",
                        "R-parity violating couplings D_i at scale Q" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new one_index_BLOCK( "IMRVD",
                                       "imaginary parts of IMRVD" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new one_index_BLOCK( "RVM2LH1IN",
      "input-scale R-parity violating soft squared masses m^2_LiH1" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new one_index_BLOCK( "IMRVM2LH1IN",
                                       "imaginary parts of IMRVM2LH1IN" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new one_index_BLOCK( "RVM2LH1",
      "R-parity violating soft squared masses m^2_LiH1 at scale Q" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new one_index_BLOCK( "IMRVM2LH1",
                                       "imaginary parts of IMRVM2LH1" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new one_index_BLOCK( "RVSNVEVIN",
                        "input-scale left-handed sneutrino vevs" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new one_index_BLOCK( "IMRVSNVEVIN",
                                       "imaginary parts of IMRVSNVEVIN" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new one_index_BLOCK( "RVSNVEV",
                        "left-handed sneutrino vevs at scale Q" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new one_index_BLOCK( "IMRVSNVEV",
                                       "imaginary parts of IMRVSNVEV" );
    BLOCKs.push_back( BLOCK_adder );

    /* the silly buggers specified a new NMIX even though NMIX is not
     * specified to be 4 * 4, so the old NMIX could do this job perfectly well
     * without another BLOCK. oh well, it's the standard now...
     */
    BLOCK_adder = new two_index_BLOCK( "RVNMIX",
                                       "neutrino-neutralino mixing matrix N" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMRVNMIX",
                                       "imaginary parts of RVNMIX" );
    BLOCKs.push_back( BLOCK_adder );

    /* the silly buggers specified a new UMIX even though UMIX is not
     * specified to be 2 * 2, so the old UMIX could do this job perfectly well
     * without another BLOCK. oh well, it's the standard now...
     */
    BLOCK_adder = new two_index_BLOCK( "RVUMIX",
                        "charged lepton-chargino mixing matrix U" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMRVUMIX",
                                       "imaginary parts of RVUMIX" );
    BLOCKs.push_back( BLOCK_adder );

    /* the silly buggers specified a new VMIX even though VMIX is not
     * specified to be 2 * 2, so the old VMIX could do this job perfectly well
     * without another BLOCK. oh well, it's the standard now...
     */
    BLOCK_adder = new two_index_BLOCK( "RVVMIX",
                        "charged lepton-chargino mixing matrix V" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMRVVMIX",
                                       "imaginary parts of RVVMIX" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "RVHMIX",
                                       "neutral scalar mixing matrix" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMRVHMIX",
                                       "imaginary parts of RVHMIX" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "RVAMIX",
                                       "neutral pseudoscalar mixing matrix" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMRVAMIX",
                                       "imaginary parts of RVAMIX" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "RVLMIX",
                                       "charged scalar mixing matrix" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMRVLMIX",
                                       "imaginary parts of RVLMIX" );
    BLOCKs.push_back( BLOCK_adder );


    // NMSSM (N.B.: remember to use a CppSLHA::EW_scale_NMSSM_spectrum!)
    // output parameters & mixing:

    BLOCK_adder = get_BLOCK( "EXTPAR" );
    if( NULL != BLOCK_adder )
      {

        BLOCK_adder->set_comment( 61,
                       "input-scale trilinear H_u-H_d-S coupling lambda" );
        BLOCK_adder->set_comment( 62,
                                  "input-scale cubic coupling kappa" );
        BLOCK_adder->set_comment( 63,
                       "input-scale soft trilinear H_u-H_d-S term A_lambda" );
        BLOCK_adder->set_comment( 64,
                                  "input-scale soft cubic term A_kappa" );
        BLOCK_adder->set_comment( 65,
                                  "input-scale singlet vev times lambda" );
        BLOCK_adder->set_comment( 66,
                                  "input-scale linear coupling xi_F" );
        BLOCK_adder->set_comment( 67,
                                  "input-scale soft linear term xi_S" );
        BLOCK_adder->set_comment( 68,
                                  "input-scale quadratic coupling mu'" );
        BLOCK_adder->set_comment( 69,
                                  "input-scale soft quadratic term m'_S^2" );
        BLOCK_adder->set_comment( 70,
                                  "input-scale soft singlet mass^2 m_S^2" );

      }
    else
      {

        if( is_verbose )
          {

            std::cout
            << std::endl
            << "CppSLHA::warning! CppSLHA2::initialize() could not add"
            << " comments to MINPAR because the BLOCK does not (yet) exist!";
            std::cout << std::endl;

          }

      }

    BLOCK_adder = new one_index_BLOCK( "NMSSMRUN",
                                       "NMSSM parameters at scale Q" );
    BLOCKs.push_back( BLOCK_adder );
    BLOCK_adder->set_comment( 1,
                   "the trilinear H_u-H_d-S coupling lambda at scale Q" );
    BLOCK_adder->set_comment( 2,
                              "the cubic coupling kappa at scale Q" );
    BLOCK_adder->set_comment( 3,
                   "the soft trilinear H_u-H_d-S term A_lambda at scale Q" );
    BLOCK_adder->set_comment( 4,
                              "the soft cubic term A_kappa at scale Q" );
    BLOCK_adder->set_comment( 5,
                              "the singlet vev times lambda at scale Q" );
    BLOCK_adder->set_comment( 6,
                              "the linear coupling xi_F at scale Q" );
    BLOCK_adder->set_comment( 7,
                              "the soft linear term xi_S at scale Q" );
    BLOCK_adder->set_comment( 8,
                              "the quadratic coupling mu' at scale Q" );
    BLOCK_adder->set_comment( 9,
                              "the soft quadratic term m'_S^2 at scale Q" );
    BLOCK_adder->set_comment( 10,
                              "the soft singlet mass^2 m_S^2 at scale Q" );

    BLOCK_adder = new two_index_BLOCK( "NMHMIX",
                                       "neutral scalar mixing matrix" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMNMHMIX",
                                       "imaginary parts of NMHMIX" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "NMAMIX",
                                       "neutral pseudoscalar mixing matrix" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMNMAMIX",
                                       "imaginary parts of NMAMIX" );
    BLOCKs.push_back( BLOCK_adder );

    /* the silly buggers specified a new NMIX even though NMIX is not
     * specified to be 4 * 4, so the old NMIX could do this job perfectly well
     * without another BLOCK. oh well, it's the standard now...
     */
    BLOCK_adder = new two_index_BLOCK( "NMNMIX",
                                       "neutrino-neutralino mixing matrix N" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder = new two_index_BLOCK( "IMNMNMIX",
                                       "imaginary parts of NMNMIX" );
    BLOCKs.push_back( BLOCK_adder );

  }



  CppSLHA2_SPheno::CppSLHA2_SPheno( std::string const given_SLHA_file_name ) :
    CppSLHA2( given_SLHA_file_name )
  // the constructor requires a file to open.
  {

    // the SLHA1 BLOCKs have already been added by the CppSLHA1 constructor.
    initialize();

  }

  CppSLHA2_SPheno::CppSLHA2_SPheno( std::string const given_SLHA_file_name,
                           EW_scale_spectrum* const given_particle_spectrum ) :
    CppSLHA2( given_SLHA_file_name,
              given_particle_spectrum )
  /* this constructor requires a file to open, & this one takes a pre-existing
   * EW_scale_spectrum, which is how one can implement the NMSSM or similar.
   */
  {

    // the SLHA1 BLOCKs have already been added by the CppSLHA1 constructor.
    initialize();

  }

  CppSLHA2_SPheno::~CppSLHA2_SPheno()
  {

    // does nothing.

  }


  void
  CppSLHA2_SPheno::initialize()
  {

    SLHA_BLOCK*
    BLOCK_adder
    = new sparse_one_index_BLOCK( "SPHENOINFO",
                                  "SPheno-specific information" );
    BLOCKs.push_back( BLOCK_adder );
    BLOCK_adder->set_comment( 1,
                              "error level" );
    BLOCK_adder->set_comment( 2,
                              "SPA conventions" );
    BLOCK_adder->set_comment( 11,
                              "branching ratios calculation switch" );
    BLOCK_adder->set_comment( 12,
                   "write only branching ratios larger than this value" );
    BLOCK_adder->set_comment( 21,
                              "cross-sections calculation switch" );

    BLOCK_adder = new SPheno_input_BLOCK( "SPHENOINPUT",
                                          "SPheno-specific input" );
    BLOCKs.push_back( BLOCK_adder );

    BLOCK_adder
    = new SPheno_cross_sections_BLOCK( "SPHENOCROSSSECTIONS",
                        "cross-sections as calculated by SPheno",
                        particle_spectrum );
    //BLOCK_adder = (SLHA_BLOCK*)shortcut_to_cross_sections;
    BLOCKs.push_back( BLOCK_adder );
    shortcut_to_cross_sections = (SPheno_cross_sections_BLOCK*)BLOCK_adder;

    BLOCK_adder = new sparse_one_index_BLOCK( "SPHENOLOWENERGY",
                        "low energy observables as calculated by SPheno" );
    BLOCKs.push_back( BLOCK_adder );
    BLOCK_adder->set_comment( 1,
                              "BR(b -> s gamma)" );
    BLOCK_adder->set_comment( 2,
                              "BR(b -> s mu+ mu-)" );
    BLOCK_adder->set_comment( 3,
                              "BR(b -> s nu nu)" );
    BLOCK_adder->set_comment( 4,
                              "BR(Bs -> mu+ mu-)" );
    BLOCK_adder->set_comment( 5,
                              "BR(B_u -> tau nu)" );
    BLOCK_adder->set_comment( 6,
                              "|Delta(M_Bd)| [ps^-1]" );
    BLOCK_adder->set_comment( 7,
                              "|Delta(M_Bs)| [ps^-1]" );
    BLOCK_adder->set_comment( 10,
                              "Delta(g-2)_electron" );
    BLOCK_adder->set_comment( 11,
                              "Delta(g-2)_muon" );
    BLOCK_adder->set_comment( 12,
                              "Delta(g-2)_tau" );
    BLOCK_adder->set_comment( 13,
                              "electric dipole moment of the electron" );
    BLOCK_adder->set_comment( 14,
                              "electric dipole moment of the muon" );
    BLOCK_adder->set_comment( 15,
                              "electric dipole moment of the tau lepton" );
    BLOCK_adder->set_comment( 16,
                              "Br(mu -> e gamma)" );
    BLOCK_adder->set_comment( 17,
                              "Br(tau -> e gamma)" );
    BLOCK_adder->set_comment( 18,
                              "Br(tau -> mu gamma)" );
    BLOCK_adder->set_comment( 19,
                              "Br(mu -> 3 e)" );
    BLOCK_adder->set_comment( 20,
                              "Br(tau -> 3 e)" );
    BLOCK_adder->set_comment( 21,
                              "Br(tau -> 3 mu)" );
    BLOCK_adder->set_comment( 30,
                              "Delta(rho_parameter)" );
    BLOCK_adder->set_comment( 40,
                              "BR(Z -> e mu)" );
    BLOCK_adder->set_comment( 41,
                              "BR(Z -> e tau)" );
    BLOCK_adder->set_comment( 42,
                              "BR(Z -> mu tau)" );

  }

}  // end of CppSLHA namespace.

