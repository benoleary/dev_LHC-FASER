/*
 * CppSLHA_EW_scale_spectrum.cpp
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

#include "CppSLHA_EW_scale_spectrum.hpp"

namespace CppSLHA
{

  bool
  drawing_label::compare_pointer_heights(
                                        drawing_label const* const first_label,
                                      drawing_label const* const second_label )
  // this returns true if the pointers are ordered so that their heights are in
  // ascending order or equal, false otherwise
  {

    if( first_label->get_label_height() <= second_label->get_label_height() )
      {

        return true;

      }
    else
      {

        return false;

      }

  }


  drawing_label::drawing_label( std::string const* const given_label,
                                double const given_height ) :
    label_string( *given_label ),
    label_height( given_height )
  {

    // just an initialization list.

  }

  drawing_label::~drawing_label()
  {

    // does nothing.

  }



  drawing_line::drawing_line( std::string const* const given_color,
                              double const given_line_height,
                              drawing_label const* const given_label ) :
    color_string( *given_color ),
    line_height( given_line_height ),
    line_label( given_label )
  {

    // just an initialization list.

  }

  drawing_line::~drawing_line()
  {

    // does nothing.

  }



  bool
  drawing_column::compare_particles_list_pointers(
                                 std::pair< particle_property_set const* const,
                                          std::string* > const* first_particle,
                                 std::pair< particle_property_set const* const,
                                        std::string* > const* second_particle )
  // this returns true if the pointers are ordered so that their absolute
  // heights are in ascending order or equal, false otherwise.
  {

    if( first_particle->first->get_absolute_mass()
        <= second_particle->first->get_absolute_mass() )
      {

        return true;

      }
    else
      {

        return false;

      }

  }

  drawing_column::drawing_column() :
    label_merging_distance( 0.005 ),
    label_clashing_distance( 0.06 ),
    label_floating_distance( 0.062 )
  {

    // just an initialization list. here is where the formatting constants
    // should be changed, along with the constructor for spectrum_drawer.

  }

  drawing_column::~drawing_column()
  {

    for( std::list< std::pair< particle_property_set const* const,
                               std::string* >* >::iterator
         deletion_iterator = particles.begin();
         particles.end() != deletion_iterator;
         ++deletion_iterator )
      {

        // 1st we delete the string:
        delete (*deletion_iterator)->second;

        // now we can delete the pair:
        delete *deletion_iterator;

      }

    for( std::list< drawing_label* >::iterator
         deletion_iterator = labels.begin();
         labels.end() != deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

    for( std::list< drawing_line* >::iterator deletion_iterator = lines.begin();
         lines.end() != deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

  }


  void
  drawing_column::prepare_lines_and_labels( double const vertical_scale )
  /* this fills the lists of lines & labels from the list of particles,
   * ensuring that the labels are in the right order & not overlapping.
   * it needs to know the overall vertical scale though.
   */
  {

    // debugging:
    /*std::cout
    << std::endl << "debugging: column has " << particles.size()
    << " particles:";
    std::cout << std::endl;*/

    // copy the particles into a list, rejecting any without recorded masses:
    std::list< std::pair< particle_property_set const* const,
                          std::string* >* > massed_particles;
    for( std::list< std::pair< particle_property_set const* const,
                               std::string* >* >::const_iterator
         particle_iterator = particles.begin();
         particles.end() != particle_iterator;
         ++particle_iterator )
      {

        // debugging:
        /*std::cout
        << std::endl << *((*particle_iterator)->first->get_name())
        << ", " << (*particle_iterator)->second << ", recorded mass: "
        << (*particle_iterator)->first->has_mass_been_recorded();
        std::cout << std::endl;*/

        if( (*particle_iterator)->first->has_mass_been_recorded() )
          {

            massed_particles.push_back( *particle_iterator );

          }

      }


    // clear any existing lines & labels:

    for( std::list< drawing_label* >::iterator
         deletion_iterator = labels.begin();
         labels.end() != deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }
    labels.clear();

    for( std::list< drawing_line* >::iterator
         deletion_iterator = lines.begin();
         lines.end() != deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }
    lines.clear();

    if( 1 == massed_particles.size() )
      // if there is only 1 particle to draw...
      {

        // create a label for the only element:
        drawing_label*
        particle_label
        = new drawing_label(
                        (*(massed_particles.begin()))->first->get_LaTeX_name(),
                            (*(massed_particles.begin()))->first->get_mass() );
        labels.push_back( particle_label );

        // make a drawing_line for it:
        lines.push_back(
                       new drawing_line( (*(massed_particles.begin()))->second,
                              (*(massed_particles.begin()))->first->get_mass(),
                                         particle_label ) );

        // remove it from the list:
        massed_particles.erase( massed_particles.begin() );

      }
    else if( 1 < massed_particles.size() )
      // if there are more than 1 particles to draw...
      {

        massed_particles.sort( &compare_particles_list_pointers );
        // sort the particles so that they are in ascending order by
        // absolute mass.

        // debugging:
        /*std::cout
        << std::endl << "sorted " << massed_particles.size()
        << " particles:";
        std::cout << std::endl;
        for( std::list< std::pair< particle_property_set const* const,
                                   std::string const >* >::const_iterator
               particle_iterator = massed_particles.begin();
             massed_particles.end() != particle_iterator;
             ++particle_iterator )
          {

            // debugging:
            std::cout
            << std::endl << *((*particle_iterator)->first->get_name())
            << ", " << (*particle_iterator)->second;
            std::cout << std::endl;

          }*/

        // prepare a pair of iterators for ease of popping out the closest
        // pairs of particles in mass:
        std::list< std::pair< particle_property_set const* const,
                              std::string* >* >::iterator
        first_particle_of_closest_pair = massed_particles.begin();
        std::list< std::pair< particle_property_set const* const,
                              std::string* >* >::iterator
        second_particle_of_closest_pair = massed_particles.begin();
        ++second_particle_of_closest_pair;
        std::list< std::pair< particle_property_set const* const,
                              std::string* >* >::iterator
        first_particle_iterator;
        std::list< std::pair< particle_property_set const* const,
                              std::string* >* >::iterator
        second_particle_iterator;
        double smallest_mass_difference;
        double current_mass_difference;

        bool finished_merging_labels_flag = false;

        while( !(massed_particles.empty()) )
          // while there are still particles to pop out...
          {

            // debugging:
            /*std::cout
            << std::endl << "!(massed_particles.empty())"
            << ", massed_particles.size() = " << massed_particles.size()
            << ", massed_particles.begin() is "
            << *((*(massed_particles.begin()))->first->get_name());*/

            if( !finished_merging_labels_flag
                &&
                ( massed_particles.size() > 1 ) )
              // if we are still looking for pairs of particles to merge
              // their labels...
              {

                first_particle_iterator = massed_particles.begin();
                second_particle_iterator = massed_particles.begin();
                ++second_particle_iterator;
                // start at the lowest-mass pair.
                smallest_mass_difference
                = ( (*second_particle_iterator)->first->get_absolute_mass()
                    - (*first_particle_iterator)->first->get_absolute_mass() );
                // this should be positive because the particles have already
                // been sorted by absolute mass.
                first_particle_of_closest_pair = first_particle_iterator;
                second_particle_of_closest_pair = second_particle_iterator;

                while( massed_particles.end() != second_particle_iterator )
                  // until we reach the last pairing...
                  {

                    // debugging:
                    /*std::cout
                    << std::endl
                    << "first_particle_iterator is "
                    << *((*first_particle_iterator)->first->get_name())
                    << ", second_particle_iterator is "
                    << *((*second_particle_iterator)->first->get_name());
                    std::cout << std::endl;*/

                    current_mass_difference
                    = ( (*second_particle_iterator)->first->get_absolute_mass()
                    - (*first_particle_iterator)->first->get_absolute_mass() );

                    // debugging:
                    /*std::cout
                    << std::endl
                    << "current_mass_difference = "
                    << current_mass_difference
                    << ", smallest_mass_difference = "
                    << smallest_mass_difference;
                    std::cout << std::endl;*/

                    if( current_mass_difference < smallest_mass_difference )
                      // if we find a pair with a smaller mass difference
                      // than the previous smallest difference, note them:
                      {

                        first_particle_of_closest_pair
                        = first_particle_iterator;
                        second_particle_of_closest_pair
                        = second_particle_iterator;
                        smallest_mass_difference = current_mass_difference;

                      }

                    // debugging:
                    /*std::cout
                    << std::endl
                    << "first_particle_of_closest_pair is "
                    << *((*first_particle_of_closest_pair)->first->get_name())
                    << ", second_particle_of_closest_pair is "
                   << *((*second_particle_of_closest_pair)->first->get_name());
                    std::cout << std::endl;*/

                    // move on to the next pair of particles:
                    ++first_particle_iterator;
                    ++second_particle_iterator;
                    // debugging:
                    /*std::cout
                    << std::endl
                    << "after iteration, first_particle_iterator is "
                    << *((*first_particle_iterator)->first->get_name())
                    << ", ( massed_particles.end()"
                    << " != second_particle_iterator ) = "
                    << ( massed_particles.end() != second_particle_iterator );
                    std::cout << std::endl;
                    if( massed_particles.end() != second_particle_iterator )
                      {

                        std::cout
                        << std::endl
                        << "after iteration, second_particle_iterator is "
                        << *((*second_particle_iterator)->first->get_name());
                        std::cout << std::endl;

                      }*/

                  }  // end of for loop over particle pairs.

                // debugging:
                /*std::cout
                << std::endl << "after"
                << " while( massed_particles.end()"
                << " != second_particle_iterator ),"
                << " first_particle_of_closest_pair is "
                << *((*first_particle_of_closest_pair)->first->get_name())
                << ", second_particle_of_closest_pair is "
                << *((*second_particle_of_closest_pair)->first->get_name())
                << std::endl
                << "( label_merging_distance * vertical_scale ) = "
                << ( label_merging_distance * vertical_scale )
                << ", smallest_mass_difference = " << smallest_mass_difference;
                std::cout << std::endl;*/


                // now we check whether the particles closest in mass need
                // to merge their labels:

                if( ( label_merging_distance * vertical_scale )
                    > smallest_mass_difference )
                  {

                    std::string
                    merged_label_string(
               *((*first_particle_of_closest_pair)->first->get_LaTeX_name()) );
                    merged_label_string.append( "," );
                    merged_label_string.append(
              *((*second_particle_of_closest_pair)->first->get_LaTeX_name()) );

                    drawing_label*
                    merged_label
                    = new drawing_label( &merged_label_string,
                                         ( 0.5
             * ( (*second_particle_of_closest_pair)->first->get_absolute_mass()
         + (*first_particle_of_closest_pair)->first->get_absolute_mass() ) ) );

                    labels.push_back( merged_label );

                    // debugging:
                    /*std::cout
                    << std::endl << "pushed back merged label.";
                    std::cout << std::endl;*/

                    // now we make drawing_line objects for them:
                    lines.push_back(
                   new drawing_line( (*first_particle_of_closest_pair)->second,
                 (*first_particle_of_closest_pair)->first->get_absolute_mass(),
                                     merged_label ) );
                    lines.push_back(
                  new drawing_line( (*second_particle_of_closest_pair)->second,
                (*second_particle_of_closest_pair)->first->get_absolute_mass(),
                                    merged_label ) );

                    // debugging:
                    /*std::cout
                    << std::endl << "about to erase pair.";
                    std::cout << std::endl;*/

                    // now we remove the pair of particles from the list:
                    massed_particles.erase( first_particle_of_closest_pair );
                    massed_particles.erase( second_particle_of_closest_pair );

                  }  // end of if the labels are close enough to merge.
                else
                  // otherwise...
                  {

                    // note that the closest particles in mass didn't need
                    // their labels merged, so no more pairs do:
                    finished_merging_labels_flag = true;

                    // make separate labels for each:
                    drawing_label*
                    first_label
                    = new drawing_label(
                    (*first_particle_of_closest_pair)->first->get_LaTeX_name(),
               (*first_particle_of_closest_pair)->first->get_absolute_mass() );
                    labels.push_back( first_label );
                    drawing_label*
                    second_label
                    = new drawing_label(
                   (*second_particle_of_closest_pair)->first->get_LaTeX_name(),
              (*second_particle_of_closest_pair)->first->get_absolute_mass() );
                    labels.push_back( second_label );

                    lines.push_back(
                   new drawing_line( (*first_particle_of_closest_pair)->second,
                 (*first_particle_of_closest_pair)->first->get_absolute_mass(),
                                     first_label ) );
                    lines.push_back(
                  new drawing_line( (*second_particle_of_closest_pair)->second,
                (*second_particle_of_closest_pair)->first->get_absolute_mass(),
                                    second_label ) );

                    // now we remove the pair of particles from the list:
                    massed_particles.erase( first_particle_of_closest_pair );
                    massed_particles.erase( second_particle_of_closest_pair );

                  }  // end of if-else for merging labels.

              }  // end of if we are still looking for pairs of particles.
            else
              {

                while( !massed_particles.empty() )
                  // until the list is emptied...
                  {

                    // create a label for the 1st element:
                    drawing_label*
                    particle_label
                    = new drawing_label(
                        (*(massed_particles.begin()))->first->get_LaTeX_name(),
                   (*(massed_particles.begin()))->first->get_absolute_mass() );
                    labels.push_back( particle_label );

                    // make a drawing_line for it:
                    lines.push_back(
                      new drawing_line( (*(massed_particles.begin()))->second,
                     (*(massed_particles.begin()))->first->get_absolute_mass(),
                                        particle_label ) );

                    // remove it from the list:
                    massed_particles.erase( massed_particles.begin() );

                  }  // end of for loop over remaining particles.

              }  // end of if-else looking for pairs of particles.

          }
        // end of while loop making all the drawing_labels & drawing_lines.

        // now order the labels from lowest to highest:
        labels.sort( &(drawing_label::compare_pointer_heights) );

        std::list< drawing_label* >::iterator
        first_label_iterator = labels.begin();
        std::list< drawing_label* >::iterator
        second_label_iterator = labels.begin();
        ++second_label_iterator;
        bool finished_floating_labels_flag = false;
        bool moved_labels_flag = false;
        double label_gap;
        double label_average;
        double
        floating_distance = ( 0.5 * label_floating_distance * vertical_scale );

        while( !finished_floating_labels_flag )
          {

            std::list< drawing_label* >::iterator
            first_label_iterator = labels.begin();
            std::list< drawing_label* >::iterator
            second_label_iterator = labels.begin();
            ++second_label_iterator;

            moved_labels_flag = false;

            while( !moved_labels_flag
                   &&
                   ( labels.end() != second_label_iterator ) )
              // until we move a pair of labels or we reach the end of the
              // list...
              {

                label_gap
                  = ( (*second_label_iterator)->get_label_height()
                      - (*first_label_iterator)->get_label_height() );
                // the label heights are all positive by this point.

                if( ( label_clashing_distance * vertical_scale ) > label_gap )
                  // if the labels are too close together (or in the wrong
                  // order because of an earlier move)...
                  {

                    label_average
                    = ( 0.5
                        * ( (*second_label_iterator)->get_label_height()
                            + (*first_label_iterator)->get_label_height() ) );

                    // move the labels away from each other:
                    (*first_label_iterator)->set_label_height(
                                       ( label_average - floating_distance ) );
                    (*second_label_iterator)->set_label_height(
                                       ( label_average + floating_distance ) );

                    moved_labels_flag = true;
                    // note that a pair of labels had to move apart, so we
                    // should look through the list again.

                  }
                else
                  // otherwise...
                  {

                    // look at the next pair of labels:
                    ++first_label_iterator;
                    ++second_label_iterator;

                  }

              }
            // end of while loop looking for a pair of labels to float apart.

            if( !moved_labels_flag )
              // if no labels were moved this time, it means that all the
              // labels in the list are far enough apart from each other.
              {

                finished_floating_labels_flag = true;

              }

          }  // end of while loop floating the labels apart.

      }  // end of if-else checking the number of particles.

  }



  spectrum_drawer::spectrum_drawer() :
    highest_height( 0.0 ),
    vertical_scale_factor( 1.2 ),
    line_gap( 24.0 ),
    line_length( 8.0 ),
    line_offset( 2.0 ),
    label_size( 3.0 )
  {

    // just an initialization list. here is where the formatting constants
    // should be changed, along with the constructor for drawing_column.

  }

  spectrum_drawer::~spectrum_drawer()
  {

    for( std::vector< drawing_column* >::iterator
         deletion_iterator = column_set.begin();
         column_set.end() != deletion_iterator;
         ++deletion_iterator )
      {

        // debugging:
        /*std::cout
        << std::endl
        << "debugging: ~spectrum_drawer deleting drawing_column pointer "
        << *deletion_iterator;
        std::cout << std::endl;*/
        delete *deletion_iterator;

      }

    // debugging:
    /*std::cout
    << std::endl
    << "debugging: ~spectrum_drawer successfully deleted its"
    << " drawing_column pointers.";
    std::cout << std::endl;*/

  }


  void
  spectrum_drawer::draw_spectrum( std::string const* const output_eps_filename,
                                  std::string const* const path_to_gnuplot,
                                  std::string const* const path_to_latex,
                                  std::string const* const path_to_dvips,
                                  std::string const* const path_to_ps2eps )
  /* this actually draws the spectrum using gnuplot, first setting up the
   * lines & labels. the "path_to_X" strings should include the executable
   * name, e.g. /usr/bin/gnuplot for path_to_gnuplot.
   */
  {

    std::string data_filename( "CppSLHA_gnuplot_data_file.data" );
    std::string command_filename( "CppSLHA_gnuplot_input_file.input" );
    std::string gnuplot_LaTeX_name_stem( "CppSLHA_gnuplot_LaTeX" );
    std::string gnuplot_label;

    double vertical_scale = ( vertical_scale_factor * find_highest_height() );
    double line_left_edge;
    double line_right_edge = 0.0;
    double label_x_value;
    int gnuplot_line_index = 0;

    // prepare a data file for gnuplot:
    std::ofstream* gnuplot_data_file = new std::ofstream();
    gnuplot_data_file->open( data_filename.c_str() );

    // prepare a file to tell gnuplot how to plot the data:
    std::ofstream* gnuplot_command_file = new std::ofstream();
    gnuplot_command_file->open( command_filename.c_str() );

    *gnuplot_command_file
    << "set term epslatex color solid"
    << std::endl << "set output \"" << gnuplot_LaTeX_name_stem << ".eps\""
    << std::endl << "set format x \"\""
    << std::endl;


    // go through all the columns, setting them up & printing their data
    // to the gnuplot files per column:
    int column_counter = 0;
    for( std::vector< drawing_column* >::iterator
         column_iterator = column_set.begin();
         column_set.end() > column_iterator;
         ++column_iterator )
      {

        // debugging:
        /*std::cout
        << std::endl << "debugging: telling column " << column_counter
        << " to prepare_lines_and_labels( " << vertical_scale << " )";
        std::cout << std::endl;*/

        (*column_iterator)->prepare_lines_and_labels( vertical_scale );

        if( 0 == ( column_counter % 2 ) )
          // if the column gets its labels to its right...
          {

            line_left_edge
            = ( ( column_counter / 2 ) * ( line_length + line_gap
                                           + line_offset )
                + line_offset );

            if( column_set.at( 0 )->is_empty() )
              // if we're not drawing column 0...
              {

                line_left_edge -= line_length;

              }

            line_right_edge = ( line_left_edge + line_length );
            label_x_value = ( line_right_edge + line_offset );

          }
        else
          {

            line_left_edge
            = ( ( ( column_counter + 1 ) / 2 ) * ( line_length + line_gap
                                                   + line_offset ) );

            if( column_set.at( 0 )->is_empty() )
              {

                line_left_edge -= line_length;

              }

            line_right_edge = ( line_left_edge + line_length );
            label_x_value = ( line_left_edge - line_offset );

          }

        for( std::list< drawing_line* >::const_iterator
             line_iterator = (*column_iterator)->get_lines()->begin();
             (*column_iterator)->get_lines()->end() != line_iterator;
             ++line_iterator )
          {

            if( 0 == ( column_counter % 2 ) )
              // if the column gets its labels to its right...
              {

                *gnuplot_data_file
                << line_left_edge << " "
                << (*line_iterator)->get_line_height()
                << std::endl
                << line_right_edge << " "
                << (*line_iterator)->get_line_height()
                << std::endl
                << label_x_value << " "
                << (*line_iterator)->get_label_height()
                << std::endl;

              }
            else
              {

                *gnuplot_data_file
                << label_x_value << " "
                << (*line_iterator)->get_label_height()
                << std::endl
                << line_left_edge << " "
                << (*line_iterator)->get_line_height()
                << std::endl
                << line_right_edge << " "
                << (*line_iterator)->get_line_height()
                << std::endl;

              }

            // gnuplot treats double blank lines as meaning that the data
            // above has a separate "index".
            *gnuplot_data_file << std::endl << std::endl;

            *gnuplot_command_file
            << "set style line " << ( gnuplot_line_index + 1 ) << " lt rgb \""
            << *((*line_iterator)->get_color()) << "\" lw 3"
            << std::endl;

            // note that we have finished with a line.
            ++gnuplot_line_index;

          }  // end of loop over lines in the current column.

        // now deal with labels for this column:

        for( std::list< drawing_label* >::const_iterator
             label_iterator = (*column_iterator)->get_labels()->begin();
             (*column_iterator)->get_labels()->end() != label_iterator;
             ++label_iterator )
          {

            if( 0 == ( column_counter % 2 ) )
              // if the column gets its labels to its right...
              {

                gnuplot_label.assign( "\\tiny " );
                gnuplot_label.append( *((*label_iterator)->get_label()) );
                gnuplot_label.append( "\\normalsize" );

                *gnuplot_command_file
                << "set label '" << gnuplot_label << "' at "
                << label_x_value << ", "
                << (*label_iterator)->get_label_height() << std::endl;


              }
            else
              {

                gnuplot_label.assign( "\\tiny" );
                gnuplot_label.append( " \\settowidth{\\labelmover}{" );
                gnuplot_label.append( *((*label_iterator)->get_label()) );
                gnuplot_label.append( "} \\hspace{-\\labelmover} " );
                gnuplot_label.append( *((*label_iterator)->get_label()) );
                gnuplot_label.append( " \\normalsize" );

                *gnuplot_command_file
                << "set label '" << gnuplot_label << "' at "
                << ( label_x_value - label_size ) << ", "
                << (*label_iterator)->get_label_height() << std::endl;

              }

          }  // end of loop over labels for the current column.

        // move on to the next column:
        ++column_counter;

      }  // end of loop over columns.

    gnuplot_data_file->close();
    // close the gnuplot data file.


    // finally construct the plot command:
    *gnuplot_command_file
    << "plot [0:" << ( line_right_edge + line_gap )
    << "] [0:" << vertical_scale << "] '" << data_filename
    << "' index 0 notitle with lines ls 1";
    for( int line_counter = 1;
         gnuplot_line_index > line_counter;
         ++line_counter )
      {

        *gnuplot_command_file
        << ", '" << data_filename
        << "' index " << line_counter
        << " notitle with lines ls " << ( line_counter + 1 );

      }
    *gnuplot_command_file << std::endl;

    gnuplot_command_file->close();
    // close the gnuplot command file.

    std::vector< std::string const* > executable_arguments( 1,
                                                           &command_filename );

    waiting_on_subprocess_executor gnuplot_executor( *path_to_gnuplot,
                                                     &executable_arguments,
                                                     10000 );
    gnuplot_executor.fork_and_execv_and_wait();
    // "gnuplot CppSLHA_gnuplot_input_file.input"

    gnuplot_executor.set_executable_name( *path_to_latex );
    executable_arguments.clear();
    std::string gnuplot_LaTeX_filename( gnuplot_LaTeX_name_stem );
    gnuplot_LaTeX_filename.append( ".tex" );
    std::string proper_LaTeX_filename( gnuplot_LaTeX_name_stem );
    proper_LaTeX_filename.append( ".proper.tex" );
    std::ofstream* proper_LaTeX_file = new std::ofstream();
    proper_LaTeX_file->open( proper_LaTeX_filename.c_str() );
    *proper_LaTeX_file
    << "\\documentclass{article}" << std::endl
    << "\\usepackage{graphics}" << std::endl
    << "\\newlength\\labelmover" << std::endl
    << "\\begin{document}" << std::endl
    << "\\pagestyle{empty}" << std::endl
    << "\\begin{center}" << std::endl
    << "\\input{./" << gnuplot_LaTeX_filename << "}" << std::endl
    << "\\end{center}" << std::endl
    << "\\end{document}" << std::endl;
    proper_LaTeX_file->close();
    delete proper_LaTeX_file;
    executable_arguments.push_back( &proper_LaTeX_filename );
    gnuplot_executor.set_arguments( &executable_arguments );
    gnuplot_executor.fork_and_execv_and_wait();
    // "latex CppSLHA_gnuplot_LaTeX.proper.tex"

    gnuplot_executor.set_executable_name( *path_to_dvips );
    executable_arguments.clear();
    //std::string dvips_argument_one( "-f" );
    //executable_arguments.push_back( &dvips_argument_one );
    std::string dvips_argument_two( gnuplot_LaTeX_name_stem );
    dvips_argument_two.append( ".proper.dvi" );
    executable_arguments.push_back( &dvips_argument_two );
    //std::string dvips_argument_three( ">" );
    //executable_arguments.push_back( &dvips_argument_three );
    //std::string dvips_argument_four( gnuplot_LaTeX_name_stem );
    //dvips_argument_four.append( ".proper.ps" );
    //executable_arguments.push_back( &dvips_argument_four );
    gnuplot_executor.set_arguments( &executable_arguments );
    gnuplot_executor.fork_and_execv_and_wait();
//"dvips -f CppSLHA_gnuplot_LaTeX.proper.dvi > CppSLHA_gnuplot_LaTeX.proper.ps"
    // actually, now just "dvips CppSLHA_gnuplot_LaTeX.proper.dvi"

#ifndef CppSLHA_NO_CLEANUP
    gnuplot_executor.set_executable_name( "/bin/rm" );
    executable_arguments.clear();
    std::string rm_argument( gnuplot_LaTeX_name_stem );
    rm_argument.append( ".proper.eps" );
    executable_arguments.push_back( &rm_argument );
    gnuplot_executor.set_arguments( &executable_arguments );
    gnuplot_executor.fork_and_execv_and_wait();
    // "rm CppSLHA_gnuplot_LaTeX.proper.eps"
#endif

    gnuplot_executor.set_executable_name( *path_to_ps2eps );
    executable_arguments.clear();
    std::string ps2eps_argument( gnuplot_LaTeX_name_stem );
    ps2eps_argument.append( ".proper.ps" );
    executable_arguments.push_back( &ps2eps_argument );
    gnuplot_executor.set_arguments( &executable_arguments );
    gnuplot_executor.fork_and_execv_and_wait();
    // "ps2eps CppSLHA_gnuplot_LaTeX.proper.ps"

    /* the following is to allow for the use of both ps2eps & ps2epsi.
     * as long as the user sets the right path to the executable, either
     *  CppSLHA_gnuplot_LaTeX.proper.eps or CppSLHA_gnuplot_LaTeX.proper.epsi
     *  will be produced. trying to cp both files to output_eps_filename
     *  should then correctly produce the required final file.
     */
    gnuplot_executor.set_executable_name( "/bin/cp" );
    executable_arguments.clear();
    std::string cp_argument_one( gnuplot_LaTeX_name_stem );
    cp_argument_one.append( ".proper.epsi" );
    executable_arguments.push_back( &cp_argument_one );
    std::string cp_argument_two( *output_eps_filename );
    executable_arguments.push_back( &cp_argument_two );
    gnuplot_executor.set_arguments( &executable_arguments );
    gnuplot_executor.fork_and_execv_and_wait();
    // "cp CppSLHA_gnuplot_LaTeX.proper.epsi [whatever output_eps_filename is]"

    // the executable name stays the same.
    executable_arguments.clear();
    cp_argument_one.assign( gnuplot_LaTeX_name_stem );
    cp_argument_one.append( ".proper.eps" );
    executable_arguments.push_back( &cp_argument_one );
    // cp argument 2 stays the same.
    executable_arguments.push_back( &cp_argument_two );
    gnuplot_executor.set_arguments( &executable_arguments );
    gnuplot_executor.fork_and_execv_and_wait();
    // "cp CppSLHA_gnuplot_LaTeX.proper.eps [whatever output_eps_filename is]"

#ifndef CppSLHA_NO_CLEANUP
    gnuplot_executor.set_executable_name( "/bin/rm" );
    executable_arguments.clear();
    rm_argument.assign( gnuplot_LaTeX_name_stem );
    rm_argument.append( ".eps" );
    executable_arguments.push_back( &rm_argument );
    gnuplot_executor.set_arguments( &executable_arguments );
    gnuplot_executor.fork_and_execv_and_wait();
    // "rm CppSLHA_gnuplot_LaTeX.eps"

    // the executable name stays the same.
    executable_arguments.clear();
    rm_argument.assign( gnuplot_LaTeX_name_stem );
    rm_argument.append( ".tex" );
    executable_arguments.push_back( &rm_argument );
    gnuplot_executor.set_arguments( &executable_arguments );
    gnuplot_executor.fork_and_execv_and_wait();
    // "rm CppSLHA_gnuplot_LaTeX.tex"

    // the executable name stays the same.
    executable_arguments.clear();
    rm_argument.assign( gnuplot_LaTeX_name_stem );
    rm_argument.append( ".proper.aux" );
    executable_arguments.push_back( &rm_argument );
    gnuplot_executor.set_arguments( &executable_arguments );
    gnuplot_executor.fork_and_execv_and_wait();
    // "rm CppSLHA_gnuplot_LaTeX.proper.aux"

    // the executable name stays the same.
    executable_arguments.clear();
    rm_argument.assign( gnuplot_LaTeX_name_stem );
    rm_argument.append( ".proper.dvi" );
    executable_arguments.push_back( &rm_argument );
    gnuplot_executor.set_arguments( &executable_arguments );
    gnuplot_executor.fork_and_execv_and_wait();
    // "rm CppSLHA_gnuplot_LaTeX.proper.dvi"

    // the executable name stays the same.
    executable_arguments.clear();
    rm_argument.assign( gnuplot_LaTeX_name_stem );
    rm_argument.append( ".proper.epsi" );
    executable_arguments.push_back( &rm_argument );
    gnuplot_executor.set_arguments( &executable_arguments );
    gnuplot_executor.fork_and_execv_and_wait();
    // "rm CppSLHA_gnuplot_LaTeX.proper.epsi"

    // the executable name stays the same.
    executable_arguments.clear();
    rm_argument.assign( gnuplot_LaTeX_name_stem );
    rm_argument.append( ".proper.eps" );
    executable_arguments.push_back( &rm_argument );
    gnuplot_executor.set_arguments( &executable_arguments );
    gnuplot_executor.fork_and_execv_and_wait();
    // "rm CppSLHA_gnuplot_LaTeX.proper.eps"

    // the executable name stays the same.
    executable_arguments.clear();
    rm_argument.assign( gnuplot_LaTeX_name_stem );
    rm_argument.append( ".proper.log" );
    executable_arguments.push_back( &rm_argument );
    gnuplot_executor.set_arguments( &executable_arguments );
    gnuplot_executor.fork_and_execv_and_wait();
    // "rm CppSLHA_gnuplot_LaTeX.proper.log"

    // the executable name stays the same.
    executable_arguments.clear();
    rm_argument.assign( gnuplot_LaTeX_name_stem );
    rm_argument.append( ".proper.ps" );
    executable_arguments.push_back( &rm_argument );
    gnuplot_executor.set_arguments( &executable_arguments );
    gnuplot_executor.fork_and_execv_and_wait();
    // "rm CppSLHA_gnuplot_LaTeX.proper.ps"

    // the executable name stays the same.
    executable_arguments.clear();
    rm_argument.assign( gnuplot_LaTeX_name_stem );
    rm_argument.append( ".proper.tex" );
    executable_arguments.push_back( &rm_argument );
    gnuplot_executor.set_arguments( &executable_arguments );
    gnuplot_executor.fork_and_execv_and_wait();
    // "rm CppSLHA_gnuplot_LaTeX.proper.tex"
#endif

  }



  EW_scale_spectrum::EW_scale_spectrum() :
    number_of_unknown_particles( 0 )
  {

    // just an initialization list.

  }

  EW_scale_spectrum::EW_scale_spectrum(
                                 EW_scale_spectrum const* const copy_source ) :
    number_of_unknown_particles( copy_source->number_of_unknown_particles )
  {

    // debugging:
    /*std::cout
    << std::endl
    << "EW_scale_spectrum( " << copy_source << " ) called.";
    std::cout << std::endl;*/

    for( std::vector< particle_property_set* >::const_iterator
         copying_iterator = copy_source->get_particle_property_sets()->begin();
         copy_source->get_particle_property_sets()->end() > copying_iterator;
         ++copying_iterator )
    // for each particle in the spectrum to be copied...
      {

        // debugging:
        /*std::cout
        << "about to copy " << (*copying_iterator)->get_name();
        std::cout << std::endl;*/

        particle_property_set*
        copied_property_set
        = add_particle_property_set( (*copying_iterator)->get_PDG_code(),
                                     (*copying_iterator)->get_name(),
                                     (*copying_iterator)->get_antiname(),
                                     (*copying_iterator)->get_LaTeX_name(),
                                     (*copying_iterator)->get_LaTeX_antiname(),
                             (*copying_iterator)->counts_as_stable_invisible(),
                                     (*copying_iterator)->counts_as_jet(),
                                 (*copying_iterator)->counts_as_light_lepton(),
                             (*copying_iterator)->counts_as_self_conjugate() );

        if( (*copying_iterator)->has_mass_been_recorded() )
          {

            copied_property_set->set_mass( (*copying_iterator)->get_mass() );

          }
        if( (*copying_iterator)->has_decay_width_been_recorded() )
          {

            copied_property_set->set_decay_width(
                                      (*copying_iterator)->get_decay_width() );
            copied_property_set->get_direct_decay_handler()->become_copy_of(
                             (*copying_iterator)->get_direct_decay_handler() );
            copied_property_set->get_cascade_handler()->become_copy_of(
                                  (*copying_iterator)->get_cascade_handler() );

          }

      }

    int column_counter = 0;
    for( std::vector< drawing_column* >::const_iterator
         column_copier
         = copy_source->get_spectrum_drawer()->get_column_set()->begin();
         copy_source->get_spectrum_drawer()->get_column_set()->end()
         > column_copier;
         ++column_copier )
      {

        for( std::list< std::pair< particle_property_set const* const,
                                   std::string* >* >::const_iterator
             copy_iterator = (*column_copier)->get_particles()->begin();
             (*column_copier)->get_particles()->end() != copy_iterator;
             ++copy_iterator )
          {

            spectrum_drawing.add_particle( get_particle_property_set(
                                     (*copy_iterator)->first->get_PDG_code() ),
                                           column_counter,
                                           (*copy_iterator)->second );

          }

        ++column_counter;

      }

  }

  EW_scale_spectrum::~EW_scale_spectrum()
  {

    for( std::vector< particle_property_set* >::iterator
         deletion_iterator = particle_property_sets.begin();
         particle_property_sets.end() > deletion_iterator;
         ++deletion_iterator )
      // go through each pointer in the std::vector.
      {

        // debugging:
        /*std::cout
        << std::endl
        << "debugging: deleting " << *((*deletion_iterator)->get_name());
        std::cout << std::endl;*/
        // delete the allocated memory pointed at by each pointer.
        delete *deletion_iterator;

      }

    // debugging:
    /*std::cout
    << std::endl
    << "debugging: ~EW_scale_spectrum successfully deleted all its"
    << " particle_property_sets.";
    std::cout << std::endl;*/

  }


  void
  EW_scale_spectrum::become_copy_of(
                                   EW_scale_spectrum const* const copy_source )
  /* this copies in all the data from the given EW_scale_spectrum, *keeping*
   * all *this* EW_scale_spectrum's pointers, updating *their* data, so that
   * for instance any shortcut pointers to this EW_scale_spectrum's particle
   * property sets will continue to work.
   */
  {

    // debugging:
    /*std::cout
    << std::endl
    << "EW_scale_spectrum::become_copy_of( " << copy_source << " ) called.";
    std::cout << std::endl;*/

    particle_property_set* copying_property_set;

    for( std::vector< particle_property_set* >::const_iterator
         copying_iterator = copy_source->get_particle_property_sets()->begin();
         copy_source->get_particle_property_sets()->end() > copying_iterator;
         ++copying_iterator )
      // for each particle in the spectrum to be copied...
      {

        // debugging:
        /*std::cout
        << "about to copy " << (*copying_iterator)->get_name();
        std::cout << std::endl;*/

        // find the pointer in *this* EW_scale_spectrum with the same PDG code
        // as *copying_iterator has:
        copying_property_set
        = get_particle_property_set( (*copying_iterator)->get_PDG_code() );

        if( NULL != copying_property_set )
          // if *this* EW_scale_spectrum has a pointer for a
          // particle_property_set with the same PDG code...
          {

            // copy the property set data from the source:
            copying_property_set->become_copy_of( *copying_iterator );

          }
        else
          // otherwise, add it in:
          {

            copying_property_set
            = add_particle_property_set( (*copying_iterator)->get_PDG_code(),
                                     (*copying_iterator)->get_name(),
                                     (*copying_iterator)->get_antiname(),
                                     (*copying_iterator)->get_LaTeX_name(),
                                     (*copying_iterator)->get_LaTeX_antiname(),
                             (*copying_iterator)->counts_as_stable_invisible(),
                                     (*copying_iterator)->counts_as_jet(),
                                 (*copying_iterator)->counts_as_light_lepton(),
                             (*copying_iterator)->counts_as_self_conjugate() );

          }

      }

    /* now all the property sets from copy_source have been copied, we check
     * for extra property sets in *this* EW_scale_spectrum which were *not* in
     * copy_source, & remove them.
     */
    if( 0 != ( particle_property_sets.size()
               - copy_source->get_particle_property_sets()->size() ) )
      {

        std::vector< particle_property_set* >::iterator
        checking_iterator = particle_property_sets.begin();
        while( particle_property_sets.end() > checking_iterator )
          // for each particle in the spectrum...
          {

            if( NULL
                == copy_source->inspect_particle_property_set(
                                       (*checking_iterator)->get_PDG_code() ) )
              // if copy_source does not have this particle...
              {

                // delete it & remove it from the set:
                delete *checking_iterator;
                checking_iterator
                = particle_property_sets.erase( checking_iterator );
                // now checking_iterator points to the element after that just
                // removed.

              }
            else
              // otherwise...
              {

                // check the next particle:
                ++checking_iterator;

              }

          }

      }

    // now the spectrum_drawer is copied:
    spectrum_drawing.reset();
    int column_counter = 0;
    for( std::vector< drawing_column* >::const_iterator
         column_copier
         = copy_source->get_spectrum_drawer()->get_column_set()->begin();
         copy_source->get_spectrum_drawer()->get_column_set()->end()
         > column_copier;
         ++column_copier )
      {

        for( std::list< std::pair< particle_property_set const* const,
                                   std::string* >* >::const_iterator
             copy_iterator = (*column_copier)->get_particles()->begin();
             (*column_copier)->get_particles()->end() != copy_iterator;
             ++copy_iterator )
          {

            spectrum_drawing.add_particle( get_particle_property_set(
                                     (*copy_iterator)->first->get_PDG_code() ),
                                           column_counter,
                                           (*copy_iterator)->second );

          }

        ++column_counter;

      }

  }


  particle_property_set*
  EW_scale_spectrum::get_particle_property_set( int sought_PDG_code )
  /* this returns the pointer to the particle_property_set with the
   * sought PDG code.  it returns NULL if the PDG code was not found. if
   * provided with a negative PDG code, it returns a pointer to the particle
   * property set with the positive PDG code, if it is found.
   */
  {

    // debugging:
    /*std::cout
    << std::endl
    << "get_particle_property_set( " << sought_PDG_code << " ) called.";
    std::cout << std::endl;*/

    particle_property_set* pointer_to_return = NULL;

    if( 0 > sought_PDG_code )
    // if told to look for an antiparticle, look for its particle instead.
      {

        sought_PDG_code = -sought_PDG_code;

      }

    for( std::vector< particle_property_set* >::iterator
         search_iterator = particle_property_sets.begin();
         particle_property_sets.end() > search_iterator;
         ++search_iterator )
    // go through each particle property set...
      {

        // debugging:
        /*std::cout
        << std::endl
        << "comparing against " << (*search_iterator)->get_PDG_code();
        std::cout << std::endl;*/

        if( sought_PDG_code == (*search_iterator)->get_PDG_code() )
        // if we find the match...
          {

            pointer_to_return = *search_iterator;
            // note the pointer to return.

            search_iterator = particle_property_sets.end();
            // stop searching.

          }

      }

    return pointer_to_return;

  }

  particle_property_set const*
  EW_scale_spectrum::inspect_particle_property_set( int sought_PDG_code )
  const
  /* this returns the pointer to the particle_property_set with the
   * sought PDG code.  it returns NULL if the PDG code was not found. if
   * provided with a negative PDG code, it returns a pointer to the particle
   * property set with the positive PDG code, if it is found.
   */
  {

    particle_property_set const* pointer_to_return = NULL;

    if( 0 > sought_PDG_code )
      // if told to look for an antiparticle, look for its particle instead.
      {

        sought_PDG_code = -sought_PDG_code;

      }

    for( std::vector< particle_property_set* >::const_iterator
         search_iterator = particle_property_sets.begin();
         particle_property_sets.end() > search_iterator;
         ++search_iterator )
      // go through each particle property set...
      {

        if( sought_PDG_code == (*search_iterator)->get_PDG_code() )
        // if we find the match...
          {

            pointer_to_return = *search_iterator;
            // note the pointer to return.

            search_iterator = particle_property_sets.end();
            // stop searching.

          }

      }

    return pointer_to_return;

  }


  double
  EW_scale_spectrum::get_PDG_coded_mass( int const input_PDG_code )
  const
  /* this checks to see if there is a particle property set with the input PDG
   * code, or the positive PDG code if a negative PDG code is provided, then
   * if there is, returns its mass, otherwise
   * CppSLHA_global::really_wrong_value.
   */
  {

    double mass_to_return = CppSLHA_global::really_wrong_value;
    // assume that the PDG code will not be found.

    particle_property_set const*
    sought_particle_property_set
    = inspect_particle_property_set( input_PDG_code );
    // find the particle property set.

    if( NULL != sought_particle_property_set )
    // if there is a particle with this PDG code...
      {

        mass_to_return = sought_particle_property_set->get_mass();
        // get the mass.

      }

    return mass_to_return;

  }

  std::string const*
  EW_scale_spectrum::get_PDG_coded_name( int const input_PDG_code )
  const
  /* this checks to see if there is a particle property set with the input
   * PDG code, then if there is, returns its name, otherwise
   * "error!".  this function accounts for negative PDG codes by
   * returning the antiparticle's name.
   */
  {

    std::string const* name_to_return = NULL;
    // assume that the PDG code will not be found.

    particle_property_set const*
    sought_particle_property_set
    = inspect_particle_property_set( input_PDG_code );
    // find the particle property set. (negative codes are already handled by
    // get_particle_property_set( input_PDG_code ) properly.)

    if( NULL != sought_particle_property_set )
    // if there is a particle with this PDG code...
      {

        if( 0 > input_PDG_code )
          // if an antiparticle's name is being sought...
          {

            name_to_return = sought_particle_property_set->get_antiname();
            // return the antiparticle's name.

          }
        else
          {

            name_to_return = sought_particle_property_set->get_name();
            // get the particle's name.

          }

      }

    return name_to_return;

  }

  std::string const*
  EW_scale_spectrum::get_PDG_coded_LaTeX_name( int const input_PDG_code )
  const
  /* this checks to see if there is a particle property set with the input
   * PDG code, then if there is, returns its LaTeX name, otherwise
   * "error!".  this function accounts for negative PDG codes by
   * returning the antiparticle's LaTeX name.
   */
  {

    std::string const* name_to_return = NULL;
    // assume that the PDG code will not be found.

    particle_property_set const*
    sought_particle_property_set
      = inspect_particle_property_set( input_PDG_code );
    // find the particle property set. (negative codes are already handled by
    // get_particle_property_set( input_PDG_code ) properly.)

    if( NULL != sought_particle_property_set )
    // if there is a particle with this PDG code...
      {

        if( 0 > input_PDG_code )
          // if an antiparticle's name is being sought...
          {

            name_to_return
            = sought_particle_property_set->get_LaTeX_antiname();
            // return the antiparticle's LaTeX name.

          }
        else
          {

            name_to_return = sought_particle_property_set->get_LaTeX_name();
            // get the particle's LaTeX name.

          }

      }

    return name_to_return;

  }


  void
  EW_scale_spectrum::record_PDG_coded_mass( int const input_PDG_code,
                                            double const input_mass )
  /* this records a mass for a particle specified by its PDG code.  if the
   * antiparticle's PDG code is given, the mass is instead recorded for the
   * particle.
   */
  {

    // debugging:
    /*std::cout
    << std::endl << "debugging:"
    << " EW_scale_spectrum::record_PDG_coded_mass( "
    << input_PDG_code << ", " << input_mass << " ) called.";
    std::cout << std::endl;
    */

    particle_property_set*
    sought_property_set = get_particle_property_set( input_PDG_code );
    // find the particle property set.

    if( NULL == sought_property_set )
    // if there is no particle with this PDG code, create a new particle
    // property set for it:
      {

        sought_property_set
        = add_unknown_particle_property_set( input_PDG_code );

      }

    sought_property_set->set_mass( input_mass );
    // set the mass now that the property set definitely exists.

  }


  particle_property_set*
  EW_scale_spectrum::add_particle_property_set( int const particle_PDG_code,
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
  {

    // sort out the antiparticle_name:
    std::string antiparticle_name( "anti-" );
    if( counts_as_self_conjugate )
      {

        antiparticle_name.assign( *particle_name );

      }
    else
      {

        antiparticle_name.append( *particle_name );

      }


    // use the PDG code as the LaTeX name:
    std::ostringstream LaTeX_name_stream;
    if( 0 > particle_PDG_code )
      {

        LaTeX_name_stream << -particle_PDG_code;

      }
    else
      {

        LaTeX_name_stream << particle_PDG_code;

      }
    std::string LaTeX_name( LaTeX_name_stream.str() );
    std::string LaTeX_antiname( "-" );
    if( counts_as_self_conjugate )
      {

        LaTeX_antiname.assign( LaTeX_name );

      }
    else
      {

        LaTeX_antiname.append( LaTeX_name );

      }

    particle_property_set*
    pointer_to_particle_to_be_added;

    if( 0 > particle_PDG_code )
      {

        pointer_to_particle_to_be_added
        = new particle_property_set( -particle_PDG_code,
                                     &antiparticle_name,
                                     particle_name,
                                     &LaTeX_name,
                                     &LaTeX_antiname,
                                     counts_as_stable_invisible,
                                     counts_as_jet,
                                     counts_as_light_lepton,
                                     counts_as_self_conjugate );

      }
    else
      {

        pointer_to_particle_to_be_added
        = new particle_property_set( particle_PDG_code,
                                     particle_name,
                                     &antiparticle_name,
                                     &LaTeX_name,
                                     &LaTeX_antiname,
                                     counts_as_stable_invisible,
                                     counts_as_jet,
                                     counts_as_light_lepton,
                                     counts_as_self_conjugate );

      }

    particle_property_sets.push_back( pointer_to_particle_to_be_added );

    std::string default_unknown_color( "black" );
    spectrum_drawing.add_particle( pointer_to_particle_to_be_added,
                                   0,
                                   &default_unknown_color );

    return pointer_to_particle_to_be_added;

  }

  particle_property_set*
  EW_scale_spectrum::add_particle_property_set( int const particle_PDG_code,
                                        std::string const* const particle_name,
                                    std::string const* const antiparticle_name,
                                         bool const counts_as_stable_invisible,
                                                bool const counts_as_jet,
                                             bool const counts_as_light_lepton,
                                          bool const counts_as_self_conjugate )
  /* this creates a new instance of a particle_property_set, adds a
   * pointer to it to the list of pointers to instances of the
   * particle_property_set class, & returns the pointer. if it receives a
   * negative PDG code, it records the positive PDG code & swaps the names it
   * received, e.g. if it received PDG code -24 & particle name "W_minus" &
   * antiparticle name "W_plus", it would record the particle with code +24,
   * name "W_plus" & antiparticle name "W_minus".
   */
  {

    // use the PDG code as the LaTeX name:
    std::ostringstream LaTeX_name_stream;
    if( 0 > particle_PDG_code )
      {

        LaTeX_name_stream << -particle_PDG_code;

      }
    else
      {

        LaTeX_name_stream << particle_PDG_code;

      }
    std::string LaTeX_name( LaTeX_name_stream.str() );
    std::string LaTeX_antiname( "-" );
    if( counts_as_self_conjugate )
      {

        LaTeX_antiname.assign( LaTeX_name );

      }
    else
      {

        LaTeX_antiname.append( LaTeX_name );

      }


    particle_property_set*
    pointer_to_particle_to_be_added;

    if( 0 > particle_PDG_code )
      {

        pointer_to_particle_to_be_added
        = new particle_property_set( -particle_PDG_code,
                                     antiparticle_name,
                                     particle_name,
                                     &LaTeX_name,
                                     &LaTeX_antiname,
                                     counts_as_stable_invisible,
                                     counts_as_jet,
                                     counts_as_light_lepton,
                                     counts_as_self_conjugate );

      }
    else
      {

        pointer_to_particle_to_be_added
        = new particle_property_set( particle_PDG_code,
                                     particle_name,
                                     antiparticle_name,
                                     &LaTeX_name,
                                     &LaTeX_antiname,
                                     counts_as_stable_invisible,
                                     counts_as_jet,
                                     counts_as_light_lepton,
                                     counts_as_self_conjugate );

      }

    particle_property_sets.push_back( pointer_to_particle_to_be_added );

    std::string default_unknown_color( "black" );
    spectrum_drawing.add_particle( pointer_to_particle_to_be_added,
                                   0,
                                   &default_unknown_color );

    return pointer_to_particle_to_be_added;

  }

  particle_property_set*
  EW_scale_spectrum::add_particle_property_set( int const particle_PDG_code,
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
  /* this creates a new instance of a particle_property_set, adds a
   * pointer to it to the list of pointers to instances of the
   * particle_property_set class, & returns the pointer. if it receives a
   * negative PDG code, it records the positive PDG code & swaps the names it
   * received, e.g. if it received PDG code -24 & particle name "W_minus" &
   * antiparticle name "W_plus", it would record the particle with code +24,
   * name "W_plus" & antiparticle name "W_minus". it also swaps the LaTeX name
   * & antiname.
   */
  {

    // debugging:
    /*std::cout
    << std::endl
    << "add_particle_property_set( " << particle_PDG_code << ", "
                                     << particle_name << ", "
                                     << antiparticle_name << ", "
                                     << LaTeX_name << ", "
                                     << LaTeX_antiname << ", "
                                     << column_to_draw_in << ", "
                                     << color_to_draw_with << ", "
                                     << counts_as_stable_invisible << ", "
                                     << counts_as_jet << ", "
                                     << counts_as_light_lepton << ", "
                                     << counts_as_self_conjugate
                                     << " ) called.";
    std::cout << std::endl;*/

    particle_property_set*
    pointer_to_particle_to_be_added;

    if( 0 > particle_PDG_code )
      {

        pointer_to_particle_to_be_added
        = new particle_property_set( -particle_PDG_code,
                                     antiparticle_name,
                                     particle_name,
                                     LaTeX_antiname,
                                     LaTeX_name,
                                     counts_as_stable_invisible,
                                     counts_as_jet,
                                     counts_as_light_lepton,
                                     counts_as_self_conjugate );

      }
    else
      {

        pointer_to_particle_to_be_added
        = new particle_property_set( particle_PDG_code,
                                     particle_name,
                                     antiparticle_name,
                                     LaTeX_name,
                                     LaTeX_antiname,
                                     counts_as_stable_invisible,
                                     counts_as_jet,
                                     counts_as_light_lepton,
                                     counts_as_self_conjugate );

      }

    particle_property_sets.push_back( pointer_to_particle_to_be_added );

    spectrum_drawing.add_particle( pointer_to_particle_to_be_added,
                                   column_to_draw_in,
                                   color_to_draw_with );

    return pointer_to_particle_to_be_added;

  }

  particle_property_set*
  EW_scale_spectrum::add_particle_property_set( int const particle_PDG_code,
                                        std::string const* const particle_name,
                                    std::string const* const antiparticle_name,
                                           std::string const* const LaTeX_name,
                                       std::string const* const LaTeX_antiname,
                                         bool const counts_as_stable_invisible,
                                                bool const counts_as_jet,
                                             bool const counts_as_light_lepton,
                                          bool const counts_as_self_conjugate )
  /* this creates a new instance of a particle_property_set, adds a
   * pointer to it to the list of pointers to instances of the
   * particle_property_set class, & returns the pointer. if it receives a
   * negative PDG code, it records the positive PDG code & swaps the names it
   * received, e.g. if it received PDG code -24 & particle name "W_minus" &
   * antiparticle name "W_plus", it would record the particle with code +24,
   * name "W_plus" & antiparticle name "W_minus". it also swaps the LaTeX name
   * & antiname. this version doesn't assign a drawing column or color because
   * it is used by the copy constructor, which finds out the drawing properties
   * by copying the spectrum drawer.
   */
  {

    particle_property_set*
    pointer_to_particle_to_be_added;

    if( 0 > particle_PDG_code )
      {

        pointer_to_particle_to_be_added
        = new particle_property_set( -particle_PDG_code,
                                     antiparticle_name,
                                     particle_name,
                                     LaTeX_antiname,
                                     LaTeX_name,
                                     counts_as_stable_invisible,
                                     counts_as_jet,
                                     counts_as_light_lepton,
                                     counts_as_self_conjugate );

      }
    else
      {

        pointer_to_particle_to_be_added
        = new particle_property_set( particle_PDG_code,
                                     particle_name,
                                     antiparticle_name,
                                     LaTeX_name,
                                     LaTeX_antiname,
                                     counts_as_stable_invisible,
                                     counts_as_jet,
                                     counts_as_light_lepton,
                                     counts_as_self_conjugate );

      }

    particle_property_sets.push_back( pointer_to_particle_to_be_added );

    return pointer_to_particle_to_be_added;

  }


  particle_property_set*
  EW_scale_spectrum::add_unknown_particle_property_set(
                                                    int constructing_PDG_code )
  /* this creates a new instance of a particle_property_set, with a
   * name indicating that its PDG code was unknown to this point, not counting
   * as stable-&-invisible, as a jet, as a light lepton or self-conjugate. it
   * then adds a pointer to it to the list of pointers to instances of the
   * particle_property_set class, & returns the pointer. if given a
   * negative PDG code, this function swaps its sign.
   */
  {

    ++number_of_unknown_particles;
    // 1st, note that there is now another unknown particle in the spectrum.

    // sort out the name:
    std::ostringstream unknown_name_stream;
    unknown_name_stream << number_of_unknown_particles;

    std::string unknown_name_string( "unknown_particle_" );
    unknown_name_string.append( unknown_name_stream.str() );
    std::string unknown_antiname_string( "unknown_antiparticle_" );
    unknown_antiname_string.append( unknown_name_stream.str() );

    if( 0 > constructing_PDG_code )
      // if told to add an unknown antiparticle, add the particle instead.
      {

        constructing_PDG_code = -constructing_PDG_code;

      }

    // actually create the new particle property set & return a pointer to it:
    return add_particle_property_set( constructing_PDG_code,
                                      &unknown_name_string,
                                      &unknown_antiname_string,
                                      false,
              /* unknown particles don't get to count as stable-&-invisible. */
                                      false,
                            /* unknown particles don't get to count as jets. */
                                      false,
                   /* unknown particles don't get to count as light leptons. */
                                      false
                /* unknown particles don't get to count as self-conjugate. */);

  }


  void
  EW_scale_spectrum::record_PDG_coded_decay_width( int const input_PDG_code,
                                               double const input_decay_width )
  /* this records a decay width for a particle specified by its PDG code. if
   * the antiparticle's PDG code is given, the decay width is instead recorded
   * for the particle.
   */
  {

    particle_property_set*
    sought_property_set
    = get_particle_property_set( input_PDG_code );
    // find the particle property set.

    if( NULL == sought_property_set )
      // if there is no particle with this PDG code, create a new particle
      // property set for it:
      {

        sought_property_set
        = add_unknown_particle_property_set( input_PDG_code );

      }

    sought_property_set->set_decay_width( input_decay_width );
    // set the decay width now that the property set definitely exists.

  }

  void
  EW_scale_spectrum::record_PDG_coded_direct_decay( int decayer_PDG_code,
                                                    int product_PDG_code_one,
                                                    int product_PDG_code_two,
                                           double const input_branching_ratio )
  /* this adds a 2-body decay to a particle property set with the specified
   * PDG code, or if no particle property set with that PDG code exists, it
   * creates a new property set & adds the decay to that. if an antiparticle's
   * PDG code is given, the decay for the particle (with decay product
   * particles replaced by antiparticles & decay product antiparticles
   * replaced by particles) is recorded instead.
   */
  {

    if( 0 > decayer_PDG_code )
    // if an antiparticle was provided...
      {

        // flip the codes of all particles involved:
        decayer_PDG_code = -decayer_PDG_code;
        product_PDG_code_one = -product_PDG_code_one;
        product_PDG_code_two = -product_PDG_code_two;

      }

    particle_property_set*
    product_one_set = get_particle_property_set( product_PDG_code_one );
    particle_property_set*
    product_two_set = get_particle_property_set( product_PDG_code_two );
    particle_property_set*
    decayer_property_set = get_particle_property_set( decayer_PDG_code );

    // check to see if the spectrum knows the decay products & if they are
    // self-conjugate, in which case ensure that they have positive sign again:
    if( NULL != product_one_set )
      {

        if( ( 0 > product_PDG_code_one )
            &&
            product_one_set->counts_as_self_conjugate() )
          {

            product_PDG_code_one = -product_PDG_code_one;

          }

      }
    if( NULL != product_two_set )
      {

        if( ( 0 > product_PDG_code_two )
            &&
            product_two_set->counts_as_self_conjugate() )
          {

            product_PDG_code_two = -product_PDG_code_two;

          }

      }

    if( NULL == decayer_property_set )
    // if there is no particle with this PDG code, create a new particle
    // property set for it:
      {

        decayer_property_set
            = add_unknown_particle_property_set( decayer_PDG_code );

      }

    decayer_property_set->add_direct_decay( product_PDG_code_one,
                                            product_PDG_code_two,
                                            input_branching_ratio );
    // now that the decaying particle's property set definitely exists, add
    // the decay to the particle's properties.

  }

  void
  EW_scale_spectrum::record_PDG_coded_direct_decay( int decayer_PDG_code,
                                                    int product_PDG_code_one,
                                                    int product_PDG_code_two,
                                                    int product_PDG_code_three,
                                           double const input_branching_ratio )
  /* this adds a 3-body decay to a particle property set with the specified
   * PDG code, or if no particle property set with that PDG code exists, it
   * creates a new property set & adds the decay to that. if an antiparticle's
   * PDG code is given, the decay for the particle (with decay product
   * particles replaced by antiparticles & decay product antiparticles
   * replaced by particles) is recorded instead.
   */
  {

    if( 0 > decayer_PDG_code )
    // if an antiparticle was provided...
      {

        // flip the codes of all particles involved:
        decayer_PDG_code = -decayer_PDG_code;
        product_PDG_code_one = -product_PDG_code_one;
        product_PDG_code_two = -product_PDG_code_two;
        product_PDG_code_three = -product_PDG_code_three;

      }

    particle_property_set*
    product_one_set = get_particle_property_set( product_PDG_code_one );
    particle_property_set*
    product_two_set = get_particle_property_set( product_PDG_code_two );
    particle_property_set*
    product_three_set = get_particle_property_set( product_PDG_code_three );
    particle_property_set*
    decayer_property_set = get_particle_property_set( decayer_PDG_code );

    // check to see if the spectrum knows the decay products & if they are
    // self-conjugate, in which case ensure that they have positive sign again:
    if( NULL != product_one_set )
      {

        if( ( 0 > product_PDG_code_one )
            &&
            product_one_set->counts_as_self_conjugate() )
          {

            product_PDG_code_one = -product_PDG_code_one;

          }

      }
    if( NULL != product_two_set )
      {

        if( ( 0 > product_PDG_code_two )
            &&
            product_two_set->counts_as_self_conjugate() )
          {

            product_PDG_code_two = -product_PDG_code_two;

          }

      }
    if( NULL != product_three_set )
      {

        if( ( 0 > product_PDG_code_three )
            &&
            product_three_set->counts_as_self_conjugate() )
          {

            product_PDG_code_three = -product_PDG_code_three;

          }

      }

    if( NULL == decayer_property_set )
    // if there is no particle with this PDG code, create a new particle
    // property set for it:
      {

        decayer_property_set
            = add_unknown_particle_property_set( decayer_PDG_code );

      }

    decayer_property_set->add_direct_decay( product_PDG_code_one,
                                            product_PDG_code_two,
                                            product_PDG_code_three,
                                            input_branching_ratio );
    // now that the decaying particle's property set definitely exists, add
    // the decay to the particle's properties.

  }


  std::string
  EW_scale_spectrum::get_decay_as_string(
                          particle_decay const* const decay_to_have_as_string )
  const
  // this returns a string that is the decay with names as well as PDG codes.
  {

    std::ostringstream decay_as_stringstream;

    decay_as_stringstream << "-> { ";

    int comma_count_to_insert
    = ( decay_to_have_as_string->inspect_decay_product_PDG_codes()->size() - 1 );

    for( std::list< int >::const_iterator
           decay_product_iterator
             = decay_to_have_as_string->inspect_decay_product_PDG_codes()->begin();
         decay_to_have_as_string->inspect_decay_product_PDG_codes()->end()
           != decay_product_iterator;
         ++decay_product_iterator )
      {

        decay_as_stringstream
        << get_PDG_coded_name( *decay_product_iterator )
        << " (" << *decay_product_iterator << ")";

        if( 0 < comma_count_to_insert )
          {

            decay_as_stringstream << ", ";

          }

        --comma_count_to_insert;

      }

    decay_as_stringstream << " }";

    return decay_as_stringstream.str();

  }


  bool
  EW_scale_spectrum::decay_unstable_particle_in_cascade(
                                         particle_property_set* given_particle,
                                 double const minimum_branching_ratio_to_keep )
  /* this looks for the 1st unstable particle in the cascade decay set of
   * given_particle, then adds in every possible cascade decay based on that
   * of given_particle with the unstable particle's direct decay products
   * replacing the particle, then removes the original decay. it returns
   * true if there was an unstable particle which was replaced, false if
   * there was no unstable particle to replace.
   */
  {

    // debugging:
    /*std::cout
    << std::endl
    << "decay_unstable_particle_in_cascade( " << given_particle
    << " [" << given_particle->get_name() << "], "
    << minimum_branching_ratio_to_keep << " ) called.";
    std::cout << std::endl;*/

    particle_property_set* current_decay_product;

    bool found_unstable_particle = false;

    for( std::list< particle_decay* >::iterator
         cascade_iterator = given_particle->get_cascade_decays()->begin();
         ( ( given_particle->get_cascade_decays()->end() != cascade_iterator )
           &&
           !found_unstable_particle );
         ++cascade_iterator )
      // look through the set of cascade decays for the 1st decay still with
      // an unstable particle.
      {

        // debugging:
        /*std::cout
        << std::endl
        << "looking at *cascade_iterator = " << *cascade_iterator
        << ", products:";
        std::cout << std::endl;*/

        for( std::list< int >::const_iterator
             decay_iterator
             = (*cascade_iterator)->get_decay_product_PDG_codes()->begin();
             ( ( (*cascade_iterator)->get_decay_product_PDG_codes()->end()
                 != decay_iterator )
               &&
               !found_unstable_particle );
             ++decay_iterator )
          // look through the list of particles in the cascade decay for the
          // 1st unstable particle.
          {

            // debugging:
            /*std::cout << *decay_iterator;
            std::cout << std::endl;*/

            current_decay_product
            = get_particle_property_set( *decay_iterator );

            if( !(current_decay_product->get_direct_decays()->empty()) )
              /* if this particle has decays, so is not stable... (there seems
               * to be no good check to prevent bad input SLHA files which
               * have A -> B + X, B -> A + Y causing infinite loops. the best
               * I can think of is to check to see if at any point, the
               * particle lists itself as a decay product. this check is done
               * in add_unstable_particle_subsituted_with_decay_products
               * below.)
               */
              {

                found_unstable_particle = true;
                /* note that we have found an unstable particle so we should
                 * stop looking for more in this particular set of decay
                 * products.
                 */

                // debugging:
                /*std::cout
                << std::endl
                << "found unstable particle " << *decay_iterator << " ("
                << current_decay_product->get_name() << ").";
                std::cout << std::endl;*/

                std::list< particle_decay* > const* decays_to_insert;

                if( !(current_decay_product->get_cascade_decays()->empty()) )
                  // if this decay product already has had its
                  // cascade decays worked out...
                  {

                    decays_to_insert
                    = current_decay_product->get_cascade_decays();
                    // use the cascade decays as replacements for the unstable
                    // particle.

                  }
                else
                  // otherwise just use the direct decays.
                  {

                    decays_to_insert
                    = current_decay_product->get_direct_decays();

                  }

                // add the set of cascades which are copies of this cascade,
                // substituting current_decay_product with its decay products:
                add_unstable_particle_substituted_with_decay_products(
                                         given_particle->get_cascade_handler(),
                                                              cascade_iterator,
                                                               *decay_iterator,
                                                              decays_to_insert,
                                             minimum_branching_ratio_to_keep );

                given_particle->get_cascade_decays()->erase(
                                                            cascade_iterator );

              }  // end of if current_decay_product can decay.

          }  // end of for loop over decay products in this cascade decay.

      }  // end of for loop over cascade decays of given_particle.

    return found_unstable_particle;

  }

  void
  EW_scale_spectrum::add_unstable_particle_substituted_with_decay_products(
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
  {

    // debugging:
    /*std::cout
    << std::endl
    << "add_unstable_particle_substituted_with_decay_products( "
    << given_cascades << ", "
    << *cascade_with_decayer << ", "
    << decaying_particle << " [" << get_PDG_coded_name( decaying_particle )
    << "], " << decays_to_insert << ", "
    << minimum_branching_ratio_to_keep << " ) called.";
    std::cout << std::endl;*/

    double new_branching_ratio;
    particle_decay* cascade_duplicate;
    bool already_replaced;

    /*
    // I had planned to replace all instances of decaying_particle in a single
    // sweep, but then realized that it's probably not worth it.
    bool cascade_has_decayer = false;
    for( std::list< particle_decay* >::iterator
           replacement_iterator = given_cascades->end();
         first_cascade_with_decayer != replacement_iterator;
         --replacement_iterator )
      // go backwards from the last cascade to the 1st with the unstable
      // particle which is being replaced by its decay products.
      {

        for( std::list< int >::const_iterator
               product_iterator
             = (*replacement_iterator)->get_decay_product_PDG_codes()->begin();
             (*replacement_iterator)->get_decay_product_PDG_codes()->end()
                                                           != product_iterator;
             ++product_iterator )
          // look through the list of particles in the cascade decay for the
          // unstable particle.
          {

            if( decaying_particle == *product_iterator )
              {

                // stuff would have gone here if I kept with this idea.

              }

          }

      }
      */

    for( std::list< particle_decay* >::const_iterator
         substitution_iterator = decays_to_insert->begin();
         decays_to_insert->end() != substitution_iterator;
         ++substitution_iterator )
      // for each way to replace the decaying particle...
      {

        new_branching_ratio
        = ( (*cascade_with_decayer)->get_branching_ratio()
            * (*substitution_iterator)->get_branching_ratio() );
        // calculate the BR for the new cascade.

        // debugging:
        /*std::cout
        << std::endl
        << "new_branching_ratio = " << new_branching_ratio;
        std::cout << std::endl;*/

        if( minimum_branching_ratio_to_keep <= new_branching_ratio )
          // if the branching ratio for this cascade decay has not gotten too
          // small...
          {

            // prepare the new cascade:
            cascade_duplicate
            = given_cascades->add_decay( *cascade_with_decayer );
            cascade_duplicate->set_branching_ratio( new_branching_ratio );

            // replace the PDG code for the unstable particle in the decay
            // with the set of PDG codes for its decay products:
            std::list< int >::iterator
            replacing_iterator
            = cascade_duplicate->get_decay_product_PDG_codes()->begin();
            already_replaced = false;
            while( !already_replaced
                   &&
                   ( cascade_duplicate->get_decay_product_PDG_codes()->end()
                     != replacing_iterator ) )
              {

                if( *replacing_iterator == decaying_particle )
                // when we find the 1st occurrance of the unstable particle's
                // PDG code in the new cascade decay...
                  {

                    // substitute the decayer's PDG code with the decay set:
                    substitute_decaying_particle_with_products(
                                                             cascade_duplicate,
                                                            replacing_iterator,
                                                        *substitution_iterator,
                                             minimum_branching_ratio_to_keep );

                    already_replaced = true;
                    /* note that this particular cascade had an unstable
                     * particle substituted by one of its sets of decay
                     * products.
                     */

                  }
                // end of if looking for the decaying particle in the cascade.

                // move on to the next decay product:
                ++replacing_iterator;

              }
            // end of while loop looking for the decaying particle in the
            // cascade.

          }
          // end of if the BR is large enough to bother with the replacement.

      }  // end of for loop over ways of replacing the decaying particle.

  }

  void
  EW_scale_spectrum::substitute_decaying_particle_with_products(
                                             particle_decay* cascade_duplicate,
                                  std::list< int >::iterator decaying_particle,
                                                particle_decay* decay_products,
                                 double const minimum_branching_ratio_to_keep )
  // this substitutes the decayer's PDG code with the given decay set.
  {

    int const decaying_particle_code = *decaying_particle;
    // this extra int is necessary because decaying_particle no longer points
    // at the decaying particle's code after the erase below.

    // debugging:
    /*std::cout
    << std::endl
    << "substitute_decaying_particle_with_products( "
    << cascade_duplicate << ", "
    << *decaying_particle << " [" << get_PDG_coded_name( *decaying_particle )
    << "], " << decay_products << ", "
    << minimum_branching_ratio_to_keep << " ) called.";
    std::cout << std::endl;*/

    // erase the unstable particle's PDG code:
    cascade_duplicate->get_decay_product_PDG_codes(
                                                 )->erase( decaying_particle );

    // add the decay products of the decay which is replacing the unstable
    // particle:
    for( std::list< int >::iterator
         replacement_iterator
          = decay_products->get_decay_product_PDG_codes()->begin();
         decay_products->get_decay_product_PDG_codes()->end()
         != replacement_iterator;
         ++replacement_iterator )
      {

        // debugging:
        /*std::cout
        << std::endl
        << "*decaying_particle = " << *decaying_particle
        << ", *replacement_iterator = " << *replacement_iterator;
        std::cout << std::endl;*/

        if( ( decaying_particle_code == *replacement_iterator )
            ||
            ( decaying_particle_code == -(*replacement_iterator) ) )
          // if the particle is in its own cascade decay...
          {

            std::cout
            << std::endl
            << "CppSLHA::error! somehow the SLHA file listed decays which"
            << " allow a particle to decay (" << *decaying_particle
            <<  ", " << get_PDG_coded_name( *decaying_particle )
            << ") back to itself (possibly via some intermediate particle or"
            << " particles). this has been discovered by"
            << " find_all_cascade_decays( " << minimum_branching_ratio_to_keep
            << " ) and is serious enough to warrant bringing the whole program"
            << " to a halt, in BOL's opinion. calling exit( EXIT_FAILURE )!";
            std::cout << std::endl;
            exit( EXIT_FAILURE );

          }
        else if( ( 0 > *decaying_particle )
            // if we are replacing an antiparticle...
                 &&
                 !(get_particle_property_set( *replacement_iterator
                                               )->counts_as_self_conjugate()) )
          // & this decay product is not self-conjugate...
          {

            cascade_duplicate->get_decay_product_PDG_codes(
                                      )->push_back( -(*replacement_iterator) );
            // add the appropriate decay product for the antiparticle.

          }
        else
          // otherwise add the decay product as it is in the
          // replacing decay:
          {

            cascade_duplicate->get_decay_product_PDG_codes(
                                         )->push_back( *replacement_iterator );

          }

      }  // end of for loop over the replacement PDG codes.

  }



  // constructor (creates an instance of the EW_scale_spectrum class
  // & then adds in the SM particles):
  EW_scale_SM_spectrum::EW_scale_SM_spectrum() :
    EW_scale_spectrum()
  {

    // debugging:
    /*std::cout
    << std::endl
    << "EW_scale_SM_spectrum() called.";
    std::cout << std::endl;*/

    std::string name_enterer;
    std::string antiname_enterer;
    std::string LaTeX_enterer;
    std::string antiLaTeX_enterer;
    std::string color_enterer;

    // add down quark:
    name_enterer.assign( "down_quark" );
    antiname_enterer.assign( "down_antiquark" );
    LaTeX_enterer.assign( "$d$" );
    antiLaTeX_enterer.assign( "${\\bar{d}}$" );
    color_enterer.assign( "black" );
    // not drawn, so doesn't matter.
    down_pointer =
    add_particle_property_set( PDG_code::down,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               -1, /* don't draw it on the mass plot. */
                               &color_enterer,
                               false, /* not counted as stable-&-invisible. */
                               true, /* is counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */ );
    down_pointer->set_mass( PDG_data::down_mass );
    // this constructor doesn't account for quarks decaying.

    // add up quark:
    name_enterer.assign( "up_quark" );
    antiname_enterer.assign( "up_antiquark" );
    LaTeX_enterer.assign( "$u$" );
    antiLaTeX_enterer.assign( "${\\bar{u}}$" );
    up_pointer =
    add_particle_property_set( PDG_code::up,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               -1, /* don't draw it on the mass plot. */
                               &color_enterer, /* not drawn, so no matter. */
                               false, /* not counted as stable-&-invisible. */
                               true, /* is counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */ );
    up_pointer->set_mass( PDG_data::up_mass );
    // this constructor doesn't account for quarks decaying.

    // add strange quark:
    name_enterer.assign( "strange_quark" );
    antiname_enterer.assign( "strange_antiquark" );
    LaTeX_enterer.assign( "$s$" );
    antiLaTeX_enterer.assign( "${\\bar{s}}$" );
    strange_pointer =
    add_particle_property_set( PDG_code::strange,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               -1, /* don't draw it on the mass plot. */
                               &color_enterer, /* not drawn, so no matter. */
                               false, /* not counted as stable-&-invisible. */
                               true, /* is counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */ );
    strange_pointer->set_mass( PDG_data::strange_mass );
    // this constructor doesn't account for quarks decaying.

    // add charm quark:
    name_enterer.assign( "charm_quark" );
    antiname_enterer.assign( "charm_antiquark" );
    LaTeX_enterer.assign( "$c$" );
    antiLaTeX_enterer.assign( "${\\bar{c}}$" );
    charm_pointer =
    add_particle_property_set( PDG_code::charm,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               -1, /* don't draw it on the mass plot. */
                               &color_enterer, /* not drawn, so no matter. */
                               false, /* not counted as stable-&-invisible. */
                               true, /* is counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */ );
    charm_pointer->set_mass( PDG_data::charm_mass );
    // this constructor doesn't account for quarks decaying.

    // add bottom quark:
    name_enterer.assign( "bottom_quark" );
    antiname_enterer.assign( "bottom_antiquark" );
    LaTeX_enterer.assign( "$b$" );
    antiLaTeX_enterer.assign( "${\\bar{b}}$" );
    bottom_pointer =
    add_particle_property_set( PDG_code::bottom,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               -1, /* don't draw it on the mass plot. */
                               &color_enterer, /* not drawn, so no matter. */
                               false, /* not counted as stable-&-invisible. */
                               true, /* is counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */ );
    bottom_pointer->set_mass( PDG_data::bottom_mass );
    // this constructor doesn't account for quarks decaying.

    // add top quark:
    name_enterer.assign( "top_quark" );
    antiname_enterer.assign( "top_antiquark" );
    LaTeX_enterer.assign( "$t$" );
    antiLaTeX_enterer.assign( "${\\bar{t}}$" );
    top_pointer =
    add_particle_property_set( PDG_code::top,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               -1, /* don't draw it on the mass plot. */
                               &color_enterer, /* not drawn, so no matter. */
                               false, /* not counted as stable-&-invisible. */
                               true, /* is counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */ );
    top_pointer->set_mass( PDG_data::top_mass );
    // this constructor doesn't account for quarks decaying.

    // add electron:
    name_enterer.assign( "electron" );
    antiname_enterer.assign( "positron" );
    LaTeX_enterer.assign( "$e$" );
    antiLaTeX_enterer.assign( "${\\bar{e}}$" );
    electron_pointer =
    add_particle_property_set( PDG_code::electron,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               -1, /* don't draw it on the mass plot. */
                               &color_enterer, /* not drawn, so no matter. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               true, /* is counted as a light lepton. */
                               false /* not counted as self-conjugate. */ );
    electron_pointer->set_mass( PDG_data::electron_mass );
    get_particle_property_set( PDG_code::electron )->set_decay_width( 0.0 );
    // electrons do not decay, as far as I know. ;]

    // add electron neutrino (currently, this code treats neutrinos as
    // Dirac fermions, so are have self_conjugate_flag = false):
    name_enterer.assign( "electron_neutrino" );
    antiname_enterer.assign( "electron_antineutrino" );
    LaTeX_enterer.assign( "${\\nu}_{e}$" );
    antiLaTeX_enterer.assign( "${\\bar{{\\nu}}}_{e}$" );
    electron_neutrino_pointer =
    add_particle_property_set( PDG_code::electron_neutrino,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               -1, /* don't draw it on the mass plot. */
                               &color_enterer, /* not drawn, so no matter. */
                               true, /* is counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */ );
    electron_neutrino_pointer->set_mass( PDG_data::electron_neutrino_mass );
    electron_neutrino_pointer->set_decay_width( 0.0 );
    // this code also doesn't care about the mismatch between flavor & mass
    // eigenstates for neutrinos, nor the decays of the heavier ones.

    // add muon:
    name_enterer.assign( "muon" );
    antiname_enterer.assign( "antimuon" );
    LaTeX_enterer.assign( "${\\mu}$" );
    antiLaTeX_enterer.assign( "${\\bar{{\\mu}}}$" );
    muon_pointer =
    add_particle_property_set( PDG_code::muon,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               -1, /* don't draw it on the mass plot. */
                               &color_enterer, /* not drawn, so no matter. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               true, /* is counted as a light lepton. */
                               false /* not counted as self-conjugate. */ );
    muon_pointer->set_mass( PDG_data::muon_mass );
    // this constructor doesn't account for muons decaying.

    // add muon neutrino (currently, this code treats neutrinos as Dirac
    // fermions, so are have self_conjugate_flag = false):
    name_enterer.assign( "muon_neutrino" );
    antiname_enterer.assign( "muon_antineutrino" );
    LaTeX_enterer.assign( "${\\nu}_{{\\mu}}$" );
    antiLaTeX_enterer.assign( "${\\bar{{\\nu}}}_{{\\mu}}$" );
    muon_neutrino_pointer =
    add_particle_property_set( PDG_code::muon_neutrino,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               -1, /* don't draw it on the mass plot. */
                               &color_enterer, /* not drawn, so no matter. */
                               true, /* is counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */ );
    muon_neutrino_pointer->set_mass( PDG_data::muon_neutrino_mass );
    muon_neutrino_pointer->set_decay_width( 0.0 );
    // this code also doesn't care about the mismatch between flavor & mass
    // eigenstates for neutrinos, nor the decays of the heavier ones.

    // add tau lepton:
    name_enterer.assign( "tau_lepton" );
    antiname_enterer.assign( "tau_antilepton" );
    LaTeX_enterer.assign( "${\\tau}$" );
    antiLaTeX_enterer.assign( "${\\bar{{\\tau}}}$" );
    tau_lepton_pointer =
    add_particle_property_set( PDG_code::tau_lepton,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               -1, /* don't draw it on the mass plot. */
                               &color_enterer, /* not drawn, so no matter. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */ );
    tau_lepton_pointer->set_mass( PDG_data::tau_lepton_mass );
    // this code doesn't account for tau leptons decaying.

    // add tau neutrino (currently, this code treats neutrinos as Dirac
    // fermions, so are have self_conjugate_flag = false):
    name_enterer.assign( "tau_neutrino" );
    antiname_enterer.assign( "tau_antineutrino" );
    LaTeX_enterer.assign( "${\\nu}_{{\\tau}}$" );
    antiLaTeX_enterer.assign( "${\\bar{{\\nu}}}_{{\\tau}}$" );
    tau_neutrino_pointer =
    add_particle_property_set( PDG_code::tau_neutrino,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               -1, /* don't draw it on the mass plot. */
                               &color_enterer, /* not drawn, so no matter. */
                               true, /* is counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */ );
    tau_neutrino_pointer->set_mass( PDG_data::tau_neutrino_mass );
    tau_neutrino_pointer->set_decay_width( 0.0 );
    // this code also doesn't care about the mismatch between flavor & mass
    // eigenstates for neutrinos, nor the decays of the heavier ones.

    // add gluon:
    name_enterer.assign( "gluon" );
    antiname_enterer.assign( "gluon" );
    LaTeX_enterer.assign( "$g$" );
    antiLaTeX_enterer.assign( "$g$" );
    gluon_pointer =
    add_particle_property_set( PDG_code::gluon,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               -1, /* don't draw it on the mass plot. */
                               &color_enterer, /* not drawn, so no matter. */
                               false, /* not counted as stable-&-invisible. */
                               true, /* is counted as a jet. */
                               false, /* not counted as a light lepton. */
                               true /* is counted as self-conjugate. */ );
    gluon_pointer->set_mass( PDG_data::gluon_mass );
    gluon_pointer->set_decay_width( 0.0 );
    // gluons do not decay, though it's not really a well-defined concept in
    // this case.

    // add photon:
    name_enterer.assign( "photon" );
    antiname_enterer.assign( "photon" );
    LaTeX_enterer.assign( "${\\gamma}$" );
    antiLaTeX_enterer.assign( "${\\gamma}$" );
    photon_pointer =
    add_particle_property_set( PDG_code::photon,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               -1, /* don't draw it on the mass plot. */
                               &color_enterer, /* not drawn, so no matter. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               true /* is counted as self-conjugate. */ );
    photon_pointer->set_mass( PDG_data::photon_mass );
    photon_pointer->set_decay_width( 0.0 );
    // photons do not decay, as far as I know. ;]

    // add Z boson (all its decays were taken from the PDG on 2009-11-08):
    name_enterer.assign( "Z_boson" );
    antiname_enterer.assign( "Z_boson" );
    LaTeX_enterer.assign( "$Z$" );
    antiLaTeX_enterer.assign( "$Z$" );
    Z_pointer =
    add_particle_property_set( PDG_code::Z,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               -1, /* don't draw it on the mass plot. */
                               &color_enterer, /* not drawn, so no matter. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               true /* is counted as self-conjugate. */ );
    Z_pointer->set_mass( PDG_data::Z_mass );
    Z_pointer->set_decay_width( PDG_data::Z_decay_width );
    record_PDG_coded_direct_decay( PDG_code::Z,
                                   PDG_code::electron,
                                   -PDG_code::electron,
                                   PDG_data::Z_to_electron_antielectron_BR );
    record_PDG_coded_direct_decay( PDG_code::Z,
                                   PDG_code::muon,
                                   -PDG_code::muon,
                                   PDG_data::Z_to_muon_antimuon_BR );
    record_PDG_coded_direct_decay( PDG_code::Z,
                                   PDG_code::tau_lepton,
                                   -PDG_code::tau_lepton,
                                   PDG_data::Z_to_tau_lepton_antilepton_BR );
    record_PDG_coded_direct_decay( PDG_code::Z,
                                   PDG_code::electron_neutrino,
                                   -PDG_code::electron_neutrino,
                          PDG_data::Z_to_electron_neutrino_antineutrino_BR );
    record_PDG_coded_direct_decay( PDG_code::Z,
                                   PDG_code::muon_neutrino,
                                   -PDG_code::muon_neutrino,
                              PDG_data::Z_to_muon_neutrino_antineutrino_BR );
    record_PDG_coded_direct_decay( PDG_code::Z,
                                   PDG_code::tau_neutrino,
                                   -PDG_code::tau_neutrino,
                               PDG_data::Z_to_tau_neutrino_antineutrino_BR );
    record_PDG_coded_direct_decay( PDG_code::Z,
                                   PDG_code::charm,
                                   -PDG_code::charm,
                                   PDG_data::Z_to_charm_anticharm_BR );
    record_PDG_coded_direct_decay( PDG_code::Z,
                                   PDG_code::bottom,
                                   -PDG_code::bottom,
                                   PDG_data::Z_to_bottom_antibottom_BR );
    record_PDG_coded_direct_decay( PDG_code::Z,
                                   PDG_code::down,
                                   -PDG_code::down,
                                   PDG_data::Z_to_down_antidown_BR );
    record_PDG_coded_direct_decay( PDG_code::Z,
                                   PDG_code::up,
                                   -PDG_code::up,
                                   PDG_data::Z_to_up_antiup_BR );
    record_PDG_coded_direct_decay( PDG_code::Z,
                                   PDG_code::strange,
                                   -PDG_code::strange,
                                   PDG_data::Z_to_strange_antistrange_BR );

    // add W boson (all its decays were taken from the PDG on 2009-11-08):
    name_enterer.assign( "W_plus_boson" );
    antiname_enterer.assign( "W_minus_boson" );
    LaTeX_enterer.assign( "$W^{+}$" );
    LaTeX_enterer.assign( "$W^{-}$" );
    W_plus_pointer =
    add_particle_property_set( PDG_code::W_plus,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               -1, /* don't draw it on the mass plot. */
                               &color_enterer, /* not drawn, so no matter. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */ );
    W_plus_pointer->set_mass( PDG_data::W_plus_mass );
    W_plus_pointer->set_decay_width( PDG_data::W_plus_decay_width );
    record_PDG_coded_direct_decay( PDG_code::W_plus,
                                   PDG_code::electron_neutrino,
                                   -PDG_code::electron,
                                PDG_data::W_plus_to_neutrino_antielectron_BR );
    record_PDG_coded_direct_decay( PDG_code::W_plus,
                                   PDG_code::muon_neutrino,
                                   -PDG_code::muon,
                                   PDG_data::W_plus_to_neutrino_antimuon_BR );
    record_PDG_coded_direct_decay( PDG_code::W_plus,
                                   PDG_code::tau_neutrino,
                                   -PDG_code::tau_lepton,
                              PDG_data::W_plus_to_neutrino_tau_antilepton_BR );
    record_PDG_coded_direct_decay( PDG_code::W_plus,
                                   PDG_code::charm,
                                   -PDG_code::down,
                                   PDG_data::W_plus_to_charm_antidown_BR );
    record_PDG_coded_direct_decay( PDG_code::W_plus,
                                   PDG_code::charm,
                                   -PDG_code::strange,
                                  PDG_data::W_plus_to_charm_antistrange_BR );
    record_PDG_coded_direct_decay( PDG_code::W_plus,
                                   PDG_code::up,
                                   -PDG_code::down,
                                   PDG_data::W_plus_to_up_antidown_BR );
    record_PDG_coded_direct_decay( PDG_code::W_plus,
                                   PDG_code::up,
                                   -PDG_code::strange,
                                   PDG_data::W_plus_to_up_antistrange_BR );

    // add electroweak-symmetry-breaking scalar boson.
    name_enterer.assign( "Higgs_boson" );
    antiname_enterer.assign( "Higgs_boson" );
    LaTeX_enterer.assign( "$h^{0}$" );
    antiLaTeX_enterer.assign( "$h^{0}$" );
    color_enterer.assign( "gold" );
    EWSB_scalar_pointer =
    add_particle_property_set( PDG_code::EWSB_scalar,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               1, /* draw it in column 1 on the mass plot. */
                               &color_enterer, /* to be drawn in gold. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               true /* is counted as self-conjugate. */);
    // no known mass as yet.

  }

  EW_scale_SM_spectrum::~EW_scale_SM_spectrum()
  {

    // does nothing.

  }



  // constructor (this constructs a EW_scale_SM_spectrum & then adds the
  // extra particles):
  EW_scale_MSSM_spectrum::EW_scale_MSSM_spectrum() :
    EW_scale_SM_spectrum()
  {

    // debugging:
    /*std::cout
    << std::endl
    << "EW_scale_MSSM_spectrum() called.";
    std::cout << std::endl;*/

    std::string name_enterer;
    std::string antiname_enterer;
    std::string LaTeX_enterer;
    std::string antiLaTeX_enterer;
    std::string color_enterer;

    // rename SM EWSB scalar:
    name_enterer.assign( "light_neutral_scalar_Higgs_boson" );
    antiname_enterer.assign( "light_neutral_scalar_Higgs_boson" );
    get_particle_property_set(
      PDG_code::light_neutral_EWSB_scalar )->set_name( name_enterer );
    light_neutral_EWSB_scalar_pointer = EWSB_scalar_pointer;

    // add heavy neutral electroweak-symmetry-breaking scalar boson:
    name_enterer.assign( "heavy_neutral_scalar_Higgs_boson" );
    antiname_enterer.assign( "heavy_neutral_scalar_Higgs_boson" );
    LaTeX_enterer.assign( "$H^{0}$" );
    antiLaTeX_enterer.assign( "$H^{0}$" );
    color_enterer.assign( "gold" );
    heavy_neutral_EWSB_scalar_pointer =
    add_particle_property_set( PDG_code::heavy_neutral_EWSB_scalar,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               1, /* draw it in column 1 on the mass plot. */
                               &color_enterer, /* to be drawn in brown. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               true /* is counted as self-conjugate. */);
    // no known mass as yet.

    // add neutral electroweak-symmetry-breaking pseudoscalar boson:
    name_enterer.assign( "neutral_pseudoscalar_Higgs_boson" );
    antiname_enterer.assign( "heavy_neutral_scalar_Higgs_boson" );
    LaTeX_enterer.assign( "$A^{0}$" );
    antiLaTeX_enterer.assign( "$A^{0}$" );
    neutral_EWSB_pseudoscalar_pointer =
    add_particle_property_set( PDG_code::neutral_EWSB_pseudoscalar,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               1, /* draw it in column 1 on the mass plot. */
                               &color_enterer, /* to be drawn in brown. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               true /* is counted as self-conjugate. */);
    // no known mass as yet.

    // add charged electroweak-symmetry-breaking scalar boson:
    name_enterer.assign( "positively_charged_scalar_Higgs_boson" );
    antiname_enterer.assign( "negatively_charged_scalar_Higgs_boson" );
    LaTeX_enterer.assign( "$H^{+}$" );
    antiLaTeX_enterer.assign( "$H^{-}$" );
    color_enterer.assign( "orangered" );
    charged_EWSB_scalar_pointer =
    add_particle_property_set( PDG_code::charged_EWSB_scalar,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               2, /* draw it in column 2 on the mass plot. */
                               &color_enterer, /* to be drawn in orangered. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add down squark associated with left-handed quarks:
    name_enterer.assign( "down_squark_L" );
    antiname_enterer.assign( "down_antisquark_L" );
    LaTeX_enterer.assign( "${\\tilde{d}}_{L}$" );
    antiLaTeX_enterer.assign( "${\\tilde{d}}_{L}^{{\\ast}}$" );
    color_enterer.assign( "blue" );
    sdown_L_pointer =
    add_particle_property_set( PDG_code::sdown_L,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               7, /* draw it in column 7 on the mass plot. */
                               &color_enterer, /* to be drawn in blue. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add down squark associated with right-handed quarks:
    name_enterer.assign( "down_squark_R" );
    antiname_enterer.assign( "down_antisquark_R" );
    LaTeX_enterer.assign( "${\\tilde{d}}_{R}$" );
    antiLaTeX_enterer.assign( "${\\tilde{d}}_{R}^{{\\ast}}$" );
    sdown_R_pointer =
    add_particle_property_set( PDG_code::sdown_R,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               7, /* draw it in column 7 on the mass plot. */
                               &color_enterer, /* to be drawn in blue. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add up squark associated with left-handed quarks:
    name_enterer.assign( "up_squark_L" );
    antiname_enterer.assign( "up_antisquark_L" );
    LaTeX_enterer.assign( "${\\tilde{u}}_{L}$" );
    antiLaTeX_enterer.assign( "${\\tilde{u}}_{L}^{{\\ast}}$" );
    sup_L_pointer =
    add_particle_property_set( PDG_code::sup_L,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               7, /* draw it in column 7 on the mass plot. */
                               &color_enterer, /* to be drawn in blue. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add up squark associated with right-handed quarks:
    name_enterer.assign( "up_squark_R" );
    antiname_enterer.assign( "up_antisquark_R" );
    LaTeX_enterer.assign( "${\\tilde{u}}_{R}$" );
    antiLaTeX_enterer.assign( "${\\tilde{u}}_{R}^{{\\ast}}$" );
    sup_R_pointer =
    add_particle_property_set( PDG_code::sup_R,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               7, /* draw it in column 7 on the mass plot. */
                               &color_enterer, /* to be drawn in blue. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add strange squark associated with left-handed quarks.
    name_enterer.assign( "strange_squark_L" );
    antiname_enterer.assign( "strange_antisquark_L" );
    LaTeX_enterer.assign( "${\\tilde{s}}_{L}$" );
    antiLaTeX_enterer.assign( "${\\tilde{s}}_{L}^{{\\ast}}$" );
    sstrange_L_pointer =
    add_particle_property_set( PDG_code::sstrange_L,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               7, /* draw it in column 7 on the mass plot. */
                               &color_enterer, /* to be drawn in blue. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add strange squark associated with right-handed quarks:
    name_enterer.assign( "strange_squark_R" );
    antiname_enterer.assign( "strange_antisquark_R" );
    LaTeX_enterer.assign( "${\\tilde{s}}_{R}$" );
    antiLaTeX_enterer.assign( "${\\tilde{s}}_{R}^{{\\ast}}$" );
    sstrange_R_pointer =
    add_particle_property_set( PDG_code::sstrange_R,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               7, /* draw it in column 7 on the mass plot. */
                               &color_enterer, /* to be drawn in blue. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add charm squark associated with left-handed quarks:
    name_enterer.assign( "charm_squark_L" );
    antiname_enterer.assign( "charm_antisquark_L" );
    LaTeX_enterer.assign( "${\\tilde{c}}_{L}$" );
    antiLaTeX_enterer.assign( "${\\tilde{c}}_{L}^{{\\ast}}$" );
    scharm_L_pointer =
    add_particle_property_set( PDG_code::scharm_L,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               7, /* draw it in column 7 on the mass plot. */
                               &color_enterer, /* to be drawn in blue. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add charm squark associated with right-handed quarks:
    name_enterer.assign( "charm_squark_R" );
    antiname_enterer.assign( "charm_antisquark_R" );
    LaTeX_enterer.assign( "${\\tilde{c}}_{R}$" );
    antiLaTeX_enterer.assign( "${\\tilde{c}}_{R}^{{\\ast}}$" );
    scharm_R_pointer =
    add_particle_property_set( PDG_code::scharm_R,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               7, /* draw it in column 7 on the mass plot. */
                               &color_enterer, /* to be drawn in blue. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add lightest bottom squark:
    name_enterer.assign( "bottom_squark_one" );
    antiname_enterer.assign( "bottom_antisquark_one" );
    LaTeX_enterer.assign( "${\\tilde{b}}_{1}$" );
    antiLaTeX_enterer.assign( "${\\tilde{b}}_{1}^{{\\ast}}$" );
    color_enterer.assign( "cyan" );
    sbottom_one_pointer =
    add_particle_property_set( PDG_code::sbottom_one,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               8, /* draw it in column 8 on the mass plot. */
                               &color_enterer, /* to be drawn in cyan. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add heaviest bottom squark:
    name_enterer.assign( "bottom_squark_two" );
    antiname_enterer.assign( "bottom_antisquark_two" );
    LaTeX_enterer.assign( "${\\tilde{b}}_{2}$" );
    antiLaTeX_enterer.assign( "${\\tilde{b}}_{2}^{{\\ast}}$" );
    sbottom_two_pointer =
    add_particle_property_set( PDG_code::sbottom_two,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               8, /* draw it in column 8 on the mass plot. */
                               &color_enterer, /* to be drawn in cyan. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add lightest top squark:
    name_enterer.assign( "top_squark_one" );
    antiname_enterer.assign( "top_antisquark_one" );
    LaTeX_enterer.assign( "${\\tilde{t}}_{1}$" );
    antiLaTeX_enterer.assign( "${\\tilde{t}}_{1}^{{\\ast}}$" );
    stop_one_pointer =
    add_particle_property_set( PDG_code::stop_one,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               8, /* draw it in column 8 on the mass plot. */
                               &color_enterer, /* to be drawn in cyan. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add heaviest top squark:
    name_enterer.assign( "top_squark_two" );
    antiname_enterer.assign( "top_antisquark_two" );
    LaTeX_enterer.assign( "${\\tilde{t}}_{2}$" );
    antiLaTeX_enterer.assign( "${\\tilde{t}}_{2}^{{\\ast}}$" );
    stop_two_pointer =
    add_particle_property_set( PDG_code::stop_two,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               8, /* draw it in column 8 on the mass plot. */
                               &color_enterer, /* to be drawn in cyan. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add selectron associated with left-handed electrons:
    name_enterer.assign( "selectron_L" );
    antiname_enterer.assign( "antiselectron_L" );
    LaTeX_enterer.assign( "${\\tilde{e}}_{L}$" );
    antiLaTeX_enterer.assign( "${\\tilde{e}}_{L}^{{\\ast}}$" );
    color_enterer.assign( "green" );
    selectron_L_pointer =
    add_particle_property_set( PDG_code::selectron_L,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               5, /* draw it in column 5 on the mass plot. */
                               &color_enterer, /* to be drawn in green. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.


    // add selectron associated with right-handed electrons:
    name_enterer.assign( "selectron_R" );
    antiname_enterer.assign( "antiselectron_R" );
    LaTeX_enterer.assign( "${\\tilde{e}}_{R}$" );
    antiLaTeX_enterer.assign( "${\\tilde{e}}_{R}^{{\\ast}}$" );
    selectron_R_pointer =
    add_particle_property_set( PDG_code::selectron_R,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               5, /* draw it in column 5 on the mass plot. */
                               &color_enterer, /* to be drawn in green. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    /* add electron sneutrino associated with left-handed electron neutrinos
     * (this constructor assumes that the LSP is not a sneutrino;
     * currently, this code treats neutrinos as Dirac fermions, so sneutrinos
     * also have self_conjugate_flag = false):
     */
    name_enterer.assign( "electron_sneutrino_L" );
    antiname_enterer.assign( "electron_antisneutrino_L" );
    LaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{eL}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{eL}^{{\\ast}}$" );
    electron_sneutrino_L_pointer =
    add_particle_property_set( PDG_code::electron_sneutrino_L,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               5, /* draw it in column 5 on the mass plot. */
                               &color_enterer, /* to be drawn in green. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    /* add electron sneutrino associated with right-handed electron neutrinos
     * (this constructor assumes that the LSP is not a sneutrino;
     * currently, this code treats neutrinos as Dirac fermions, so sneutrinos
     * also have self_conjugate_flag = false):
     */
    name_enterer.assign( "electron_sneutrino_R" );
    antiname_enterer.assign( "electron_antisneutrino_R" );
    LaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{eR}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{eR}^{{\\ast}}$" );
    electron_sneutrino_R_pointer =
    add_particle_property_set( PDG_code::electron_sneutrino_R,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               5, /* draw it in column 5 on the mass plot. */
                               &color_enterer, /* to be drawn in green. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add smuon associated with left-handed muons:
    name_enterer.assign( "smuon_L" );
    antiname_enterer.assign( "antismuon_L" );
    LaTeX_enterer.assign( "${\\tilde{{\\mu}}}_{L}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\mu}}}_{L}^{{\\ast}}$" );
    smuon_L_pointer =
    add_particle_property_set( PDG_code::smuon_L,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               5, /* draw it in column 5 on the mass plot. */
                               &color_enterer, /* to be drawn in green. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add smuon associated with right-handed muons:
    name_enterer.assign( "smuon_R" );
    antiname_enterer.assign( "antismuon_R" );
    LaTeX_enterer.assign( "${\\tilde{{\\mu}}}_{R}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\mu}}}_{R}^{{\\ast}}$" );
    smuon_R_pointer =
    add_particle_property_set( PDG_code::smuon_R,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               5, /* draw it in column 5 on the mass plot. */
                               &color_enterer, /* to be drawn in green. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    /* add muon sneutrino associated with left-handed muon neutrinos
     * (this constructor assumes that the LSP is not a sneutrino;
     * currently, this code treats neutrinos as Dirac fermions, so sneutrinos
     * also have self_conjugate_flag = false):
     */
    name_enterer.assign( "muon_sneutrino_L" );
    antiname_enterer.assign( "electron_antisneutrino_L" );
    LaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{{\\mu}L}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{{\\mu}L}^{{\\ast}}$" );
    muon_sneutrino_L_pointer =
    add_particle_property_set( PDG_code::muon_sneutrino_L,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               5, /* draw it in column 5 on the mass plot. */
                               &color_enterer, /* to be drawn in green. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    /* add muon sneutrino associated with right-handed muon neutrinos
     * (this constructor assumes that the LSP is not a sneutrino;
     * currently, this code treats neutrinos as Dirac fermions, so sneutrinos
     * also have self_conjugate_flag = false):
     */
    name_enterer.assign( "muon_sneutrino_R" );
    antiname_enterer.assign( "muon_antisneutrino_R" );
    LaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{{\\mu}R}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{{\\mu}R}^{{\\ast}}$" );
    muon_sneutrino_R_pointer =
    add_particle_property_set( PDG_code::muon_sneutrino_R,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               5, /* draw it in column 5 on the mass plot. */
                               &color_enterer, /* to be drawn in green. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add lightest tau slepton:
    name_enterer.assign( "stau_one" );
    antiname_enterer.assign( "antistau_one" );
    LaTeX_enterer.assign( "${\\tilde{{\\tau}}}_{1}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\tau}}}_{1}^{{\\ast}}$" );
    color_enterer.assign( "turquoise" );
    stau_one_pointer =
    add_particle_property_set( PDG_code::stau_one,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               6, /* draw it in column 6 on the mass plot. */
                               &color_enterer, /* to be drawn in turquoise. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add heaviest tau slepton:
    name_enterer.assign( "stau_two" );
    antiname_enterer.assign( "antistau_two" );
    LaTeX_enterer.assign( "${\\tilde{{\\tau}}}_{2}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\tau}}}_{2}^{{\\ast}}$" );
    stau_two_pointer =
    add_particle_property_set( PDG_code::stau_two,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               6, /* draw it in column 6 on the mass plot. */
                               &color_enterer, /* to be drawn in turquoise. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    /* add tau sneutrino associated with left-handed tau sneutrinos
     * (this constructor assumes that the LSP is not a sneutrino);
     * currently, this code treats neutrinos as Dirac fermions, so sneutrinos
     * also have self_conjugate_flag = false):
     */
    name_enterer.assign( "tau_sneutrino_L" );
    antiname_enterer.assign( "tau_antisneutrino_L" );
    LaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{{\\tau}L}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{{\\tau}L}^{{\\ast}}$" );
    tau_sneutrino_L_pointer =
    add_particle_property_set( PDG_code::tau_sneutrino_L,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               6, /* draw it in column 6 on the mass plot. */
                               &color_enterer, /* to be drawn in turquoise. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    /* add tau sneutrino associated with right-handed tau sneutrinos
     * (this constructor assumes that the LSP is not a sneutrino;
     * currently, this code treats neutrinos as Dirac fermions, so sneutrinos
     * also have self_conjugate_flag = false);
     */
    name_enterer.assign( "tau_sneutrino_R" );
    antiname_enterer.assign( "tau_antisneutrino_R" );
    LaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{{\\tau}R}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{{\\tau}R}^{{\\ast}}$" );
    tau_sneutrino_R_pointer =
    add_particle_property_set( PDG_code::tau_sneutrino_R,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               6, /* draw it in column 6 on the mass plot. */
                               &color_enterer, /* to be drawn in turquoise. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add gluino:
    name_enterer.assign( "gluino" );
    antiname_enterer.assign( "gluino" );
    LaTeX_enterer.assign( "${\\tilde{g}}$" );
    antiLaTeX_enterer.assign( "${\\tilde{g}}$" );
    color_enterer.assign( "purple" );
    gluino_pointer =
    add_particle_property_set( PDG_code::gluino,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               4, /* draw it in column 4 on the mass plot. */
                               &color_enterer, /* to be drawn in purple. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               true /* is counted as self-conjugate. */);
    // no known mass as yet.

    // add lightest neutralino (this constructor assumes that the LSP is the
    // lightest neutralino):
    name_enterer.assign( "neutralino_1" );
    antiname_enterer.assign( "neutralino_1" );
    LaTeX_enterer.assign( "${\\tilde{{\\chi}}}_{1}^{0}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\chi}}}_{1}^{0}$" );
    color_enterer.assign( "orange" );
    neutralino_one_pointer =
    add_particle_property_set( PDG_code::neutralino_one,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               3, /* draw it in column 3 on the mass plot. */
                               &color_enterer, /* to be drawn in orange. */
                               true, /* is counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               true /* is counted as self-conjugate. */);
    // no known mass as yet.

    // add 2nd-lightest neutralino:
    name_enterer.assign( "neutralino_2" );
    antiname_enterer.assign( "neutralino_2" );
    LaTeX_enterer.assign( "${\\tilde{{\\chi}}}_{2}^{0}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\chi}}}_{2}^{0}$" );
    neutralino_two_pointer =
    add_particle_property_set( PDG_code::neutralino_two,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               3, /* draw it in column 3 on the mass plot. */
                               &color_enterer, /* to be drawn in orange. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               true /* is counted as self-conjugate. */);
    // no known mass as yet.

    // add 3rd-lightest neutralino:
    name_enterer.assign( "neutralino_3" );
    antiname_enterer.assign( "neutralino_3" );
    LaTeX_enterer.assign( "${\\tilde{{\\chi}}}_{3}^{0}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\chi}}}_{3}^{0}$" );
    neutralino_three_pointer =
    add_particle_property_set( PDG_code::neutralino_three,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               3, /* draw it in column 3 on the mass plot. */
                               &color_enterer, /* to be drawn in orange. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               true /* is counted as self-conjugate. */);
    // no known mass as yet.

    // add 4th-lightest (i.e heaviest in the MSSM) neutralino:
    name_enterer.assign( "neutralino_4" );
    antiname_enterer.assign( "neutralino_4" );
    LaTeX_enterer.assign( "${\\tilde{{\\chi}}}_{4}^{0}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\chi}}}_{4}^{0}$" );
    neutralino_four_pointer =
    add_particle_property_set( PDG_code::neutralino_four,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               3, /* draw it in column 3 on the mass plot. */
                               &color_enterer, /* to be drawn in orange. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               true /* is counted as self-conjugate. */);
    // no known mass as yet.

    // add lightest chargino:
    name_enterer.assign( "positive_chargino_1" );
    antiname_enterer.assign( "negative_chargino_1" );
    LaTeX_enterer.assign( "${\\tilde{{\\chi}}}_{1}^{+}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\chi}}}_{1}^{-}$" );
    color_enterer.assign( "red" );
    chargino_one_pointer =
    add_particle_property_set( PDG_code::chargino_one,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               4, /* draw it in column 4 on the mass plot. */
                               &color_enterer, /* to be drawn in red. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add 2nd-lightest (i.e heaviest in the MSSM) chargino:
    name_enterer.assign( "positive_chargino_2" );
    antiname_enterer.assign( "negative_chargino_2" );
    LaTeX_enterer.assign( "${\\tilde{{\\chi}}}_{2}^{+}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\chi}}}_{2}^{-}$" );
    chargino_two_pointer =
    add_particle_property_set( PDG_code::chargino_two,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               4, /* draw it in column 4 on the mass plot. */
                               &color_enterer, /* to be drawn in red. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

  }

  EW_scale_MSSM_spectrum::~EW_scale_MSSM_spectrum()
  {

    // does nothing.

  }


  particle_property_set*
  EW_scale_MSSM_spectrum::get_lighter_particle(
                                particle_property_set* const compared_particle,
                                                int const comparison_PDG_code )
  // this compares the masses of 2 particles & returns a pointer to the
  // lighter particle's property set.
  {

    particle_property_set* assumed_lighter = compared_particle;

    double assumed_lighter_mass = assumed_lighter->get_absolute_mass();

    particle_property_set*
    assumed_heavier = get_particle_property_set( comparison_PDG_code );

    double assumed_heavier_mass = assumed_heavier->get_absolute_mass();

    if( assumed_lighter_mass > assumed_heavier_mass )
      // if the guess as to which was lighter was wrong...
      {

        return assumed_heavier;
        // correct the pointer to the lighter particle.

      }
    else
      {

        return assumed_lighter;

      }

  }

  particle_property_set*
  EW_scale_MSSM_spectrum::update_stable_invisible()
  /* this checks the sneutrinos & lightest neutralino for which is lightest, &
   * returns that, while also setting its is_stable_invisible to true, if this
   * returned sparticle has an empty decay list. this *can* be used to
   * identify if the LSP is cosmologically allowed, by running
   * ->counts_as_stable_invisible() on the pointer which this function
   * returns, e.g. if update_stable_invisible()->counts_as_stable_invisible()
   * is false, then the lightest cosmologically-acceptable (in the sense of
   * being neutral to the color & electromagnetism forces) dark matter
   * candidate is unstable, & thus cannot be dark matter (because it decays,
   * into something, which, assuming R-parity is conserved in all the recorded
   * decays, is presumably an electromagnetically- and/or color-charged
   * sparticle as well as SM particles).
   * actually, there is a crazy fringe case where that won't work: if for
   * example the lightest sneutrino is lighter than the neutralinos, but the
   * LSP is actually stop_1, and the lightest sneutrino mass is less than the
   * charm mass heavier the stop_1 mass, and whichever spectrum generator that
   * worked out the branching ratios didn't catch the 5-body decay into a
   * neutrino, stop_1, positron, another neutrino, and an up quark through a
   * virtual top, or the 3-body decay through the sneutrino-electron-chargino
   * lepton-flavor-violating vertex followed by the V_td-suppressed
   * chargino-stop-down vertex. honestly though, I don't think it's the job of
   * this code to catch such psychotic cases.
   * this function could be more comprehensive. I might sort it out later.
   */
  {

    // start by assuming that the lightest neutralino is lighter than all the
    // sneutrinos:
    particle_property_set* assumed_LSP
    = get_particle_property_set( PDG_code::neutralino_one );

    // compare against each sneutrino, noting the lighter particle each time:
    assumed_LSP = get_lighter_particle( assumed_LSP,
                                        PDG_code::electron_sneutrino_L );
    assumed_LSP = get_lighter_particle( assumed_LSP,
                                        PDG_code::muon_sneutrino_L );
    assumed_LSP = get_lighter_particle( assumed_LSP,
                                        PDG_code::tau_sneutrino_L );
    assumed_LSP = get_lighter_particle( assumed_LSP,
                                        PDG_code::electron_sneutrino_R );
    assumed_LSP = get_lighter_particle( assumed_LSP,
                                        PDG_code::muon_sneutrino_R );
    assumed_LSP = get_lighter_particle( assumed_LSP,
                                        PDG_code::tau_sneutrino_R );

    if( PDG_code::neutralino_one != assumed_LSP->get_PDG_code() )
    // if, after that, it turns out that the lightest neutralino is not the
    // LSP...
      {

        get_particle_property_set(
          PDG_code::neutralino_one )->make_count_as_stable_invisible( false );
        // flag the lightest neutralino as no longer escaping the detector.

      }

    if( !( assumed_LSP->get_direct_decays()->empty() ) )
      // if this supposed LSP has a list of decays...
      {

        assumed_LSP->make_count_as_stable_invisible( false );
        // flag it as no longer escaping the detector.

      }
    else
      {

        assumed_LSP->make_count_as_stable_invisible( true );
        // flag it as escaping the detector.

      }

    return assumed_LSP;

  }



  // constructor (this constructs a EW_scale_SM_spectrum & then adds the
  // extra particles):
  EW_scale_FVMSSM_spectrum::EW_scale_FVMSSM_spectrum() :
    EW_scale_SM_spectrum()
  {

    std::string name_enterer;
    std::string antiname_enterer;
    std::string LaTeX_enterer;
    std::string antiLaTeX_enterer;
    std::string color_enterer;

    // rename SM EWSB scalar:
    name_enterer.assign( "light_neutral_scalar_Higgs_boson" );
    antiname_enterer.assign( "light_neutral_scalar_Higgs_boson" );
    get_particle_property_set(
               PDG_code::light_neutral_EWSB_scalar )->set_name( name_enterer );

    // add heavy neutral electroweak-symmetry-breaking scalar boson:
    name_enterer.assign( "heavy_neutral_scalar_Higgs_boson" );
    antiname_enterer.assign( "heavy_neutral_scalar_Higgs_boson" );
    LaTeX_enterer.assign( "$H^{0}$" );
    antiLaTeX_enterer.assign( "$H^{0}$" );
    color_enterer.assign( "gold" );
    heavy_neutral_EWSB_scalar_pointer =
    add_particle_property_set( PDG_code::heavy_neutral_EWSB_scalar,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               1, /* draw it in column 1 on the mass plot. */
                               &color_enterer, /* to be drawn in brown. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               true /* is counted as self-conjugate. */);
    // no known mass as yet.

    // add neutral electroweak-symmetry-breaking pseudoscalar boson:
    name_enterer.assign( "neutral_pseudoscalar_Higgs_boson" );
    antiname_enterer.assign( "heavy_neutral_scalar_Higgs_boson" );
    LaTeX_enterer.assign( "$A^{0}$" );
    antiLaTeX_enterer.assign( "$A^{0}$" );
    neutral_EWSB_pseudoscalar_pointer =
    add_particle_property_set( PDG_code::neutral_EWSB_pseudoscalar,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               1, /* draw it in column 1 on the mass plot. */
                               &color_enterer, /* to be drawn in brown. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               true /* is counted as self-conjugate. */);
    // no known mass as yet.

    // add charged electroweak-symmetry-breaking scalar boson:
    name_enterer.assign( "positively_charged_scalar_Higgs_boson" );
    antiname_enterer.assign( "negatively_charged_scalar_Higgs_boson" );
    LaTeX_enterer.assign( "$H^{+}$" );
    antiLaTeX_enterer.assign( "$H^{-}$" );
    color_enterer.assign( "orangered" );
    charged_EWSB_scalar_pointer =
    add_particle_property_set( PDG_code::charged_EWSB_scalar,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               2, /* draw it in column 2 on the mass plot. */
                               &color_enterer, /* to be drawn in orangered. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add lightest down squark:
    name_enterer.assign( "down_squark_one" );
    antiname_enterer.assign( "down_antisquark_one" );
    LaTeX_enterer.assign( "${\\tilde{d}}_{1}$" );
    antiLaTeX_enterer.assign( "${\\tilde{d}}_{1}^{{\\ast}}$" );
    color_enterer.assign( "blue" );
    sdown_one_pointer =
    add_particle_property_set( PDG_code::sdown_one,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               7, /* draw it in column 7 on the mass plot. */
                               &color_enterer, /* to be drawn in blue. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add 2nd-lightest down squark:
    name_enterer.assign( "down_squark_two" );
    antiname_enterer.assign( "down_antisquark_two" );
    LaTeX_enterer.assign( "${\\tilde{d}}_{2}$" );
    antiLaTeX_enterer.assign( "${\\tilde{d}}_{2}^{{\\ast}}$" );
    sdown_two_pointer =
    add_particle_property_set( PDG_code::sdown_two,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               7, /* draw it in column 7 on the mass plot. */
                               &color_enterer, /* to be drawn in blue. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add 3rd-lightest down squark:
    name_enterer.assign( "down_squark_three" );
    antiname_enterer.assign( "down_antisquark_three" );
    LaTeX_enterer.assign( "${\\tilde{d}}_{3}$" );
    antiLaTeX_enterer.assign( "${\\tilde{d}}_{3}^{{\\ast}}$" );
    sdown_three_pointer =
    add_particle_property_set( PDG_code::sdown_three,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               7, /* draw it in column 7 on the mass plot. */
                               &color_enterer, /* to be drawn in blue. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add 4th-lightest down squark:
    name_enterer.assign( "down_squark_four" );
    antiname_enterer.assign( "down_antisquark_four" );
    LaTeX_enterer.assign( "${\\tilde{d}}_{4}$" );
    antiLaTeX_enterer.assign( "${\\tilde{d}}_{4}^{{\\ast}}$" );
    sdown_four_pointer =
    add_particle_property_set( PDG_code::sdown_four,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               7, /* draw it in column 7 on the mass plot. */
                               &color_enterer, /* to be drawn in blue. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add 5th-lightest down squark:
    name_enterer.assign( "down_squark_five" );
    antiname_enterer.assign( "down_antisquark_five" );
    LaTeX_enterer.assign( "${\\tilde{d}}_{5}$" );
    antiLaTeX_enterer.assign( "${\\tilde{d}}_{5}^{{\\ast}}$" );
    sdown_five_pointer =
    add_particle_property_set( PDG_code::sdown_five,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               7, /* draw it in column 7 on the mass plot. */
                               &color_enterer, /* to be drawn in blue. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add 6th-lightest down squark:
    name_enterer.assign( "down_squark_six" );
    antiname_enterer.assign( "down_antisquark_six" );
    LaTeX_enterer.assign( "${\\tilde{d}}_{6}$" );
    antiLaTeX_enterer.assign( "${\\tilde{d}}_{6}^{{\\ast}}$" );
    sdown_six_pointer =
    add_particle_property_set( PDG_code::sdown_six,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               7, /* draw it in column 7 on the mass plot. */
                               &color_enterer, /* to be drawn in blue. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add lightest up squark:
    name_enterer.assign( "up_squark_one" );
    antiname_enterer.assign( "up_antisquark_one" );
    LaTeX_enterer.assign( "${\\tilde{u}}_{1}$" );
    antiLaTeX_enterer.assign( "${\\tilde{u}}_{1}^{{\\ast}}$" );
    color_enterer.assign( "cyan" );
    sup_one_pointer =
    add_particle_property_set( PDG_code::sup_one,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               8, /* draw it in column 8 on the mass plot. */
                               &color_enterer, /* to be drawn in cyan. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add 2nd-lightest up squark:
    name_enterer.assign( "up_squark_two" );
    antiname_enterer.assign( "up_antisquark_two" );
    LaTeX_enterer.assign( "${\\tilde{u}}_{2}$" );
    antiLaTeX_enterer.assign( "${\\tilde{u}}_{2}^{{\\ast}}$" );
    sup_two_pointer =
    add_particle_property_set( PDG_code::sup_two,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               8, /* draw it in column 8 on the mass plot. */
                               &color_enterer, /* to be drawn in cyan. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add 3rd-lightest up squark:
    name_enterer.assign( "up_squark_three" );
    antiname_enterer.assign( "up_antisquark_three" );
    LaTeX_enterer.assign( "${\\tilde{u}}_{3}$" );
    antiLaTeX_enterer.assign( "${\\tilde{u}}_{3}^{{\\ast}}$" );
    sup_three_pointer =
    add_particle_property_set( PDG_code::sup_three,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               8, /* draw it in column 8 on the mass plot. */
                               &color_enterer, /* to be drawn in cyan. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add 4th-lightest up squark:
    name_enterer.assign( "up_squark_four" );
    antiname_enterer.assign( "up_antisquark_four" );
    LaTeX_enterer.assign( "${\\tilde{u}}_{4}$" );
    antiLaTeX_enterer.assign( "${\\tilde{u}}_{4}^{{\\ast}}$" );
    sup_four_pointer =
    add_particle_property_set( PDG_code::sup_four,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               8, /* draw it in column 8 on the mass plot. */
                               &color_enterer, /* to be drawn in cyan. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add 5th-lightest up squark:
    name_enterer.assign( "up_squark_five" );
    antiname_enterer.assign( "up_antisquark_five" );
    LaTeX_enterer.assign( "${\\tilde{u}}_{5}$" );
    antiLaTeX_enterer.assign( "${\\tilde{u}}_{5}^{{\\ast}}$" );
    sup_five_pointer =
    add_particle_property_set( PDG_code::sup_five,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               8, /* draw it in column 8 on the mass plot. */
                               &color_enterer, /* to be drawn in cyan. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add 6th-lightest up squark:
    name_enterer.assign( "up_squark_six" );
    antiname_enterer.assign( "up_antisquark_six" );
    LaTeX_enterer.assign( "${\\tilde{u}}_{6}$" );
    antiLaTeX_enterer.assign( "${\\tilde{u}}_{6}^{{\\ast}}$" );
    sup_six_pointer =
    add_particle_property_set( PDG_code::sup_six,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               8, /* draw it in column 8 on the mass plot. */
                               &color_enterer, /* to be drawn in cyan. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add lightest charged slepton:
    name_enterer.assign( "selectron_one" );
    antiname_enterer.assign( "antiselectron_one" );
    LaTeX_enterer.assign( "${\\tilde{e}}_{1}$" );
    antiLaTeX_enterer.assign( "${\\tilde{e}}_{1}^{{\\ast}}$" );
    color_enterer.assign( "green" );
    selectron_one_pointer =
    add_particle_property_set( PDG_code::selectron_one,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               5, /* draw it in column 5 on the mass plot. */
                               &color_enterer, /* to be drawn in green. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add 2nd-lightest charged slepton:
    name_enterer.assign( "selectron_two" );
    antiname_enterer.assign( "antiselectron_two" );
    LaTeX_enterer.assign( "${\\tilde{e}}_{2}$" );
    antiLaTeX_enterer.assign( "${\\tilde{e}}_{2}^{{\\ast}}$" );
    selectron_two_pointer =
    add_particle_property_set( PDG_code::selectron_two,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               5, /* draw it in column 5 on the mass plot. */
                               &color_enterer, /* to be drawn in green. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add 3rd-lightest charged slepton:
    name_enterer.assign( "selectron_three" );
    antiname_enterer.assign( "antiselectron_three" );
    LaTeX_enterer.assign( "${\\tilde{e}}_{3}$" );
    antiLaTeX_enterer.assign( "${\\tilde{e}}_{3}^{{\\ast}}$" );
    selectron_three_pointer =
    add_particle_property_set( PDG_code::selectron_three,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               5, /* draw it in column 5 on the mass plot. */
                               &color_enterer, /* to be drawn in green. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add 4th-lightest charged slepton:
    name_enterer.assign( "selectron_four" );
    antiname_enterer.assign( "antiselectron_four" );
    LaTeX_enterer.assign( "${\\tilde{e}}_{4}$" );
    antiLaTeX_enterer.assign( "${\\tilde{e}}_{4}^{{\\ast}}$" );
    selectron_four_pointer =
    add_particle_property_set( PDG_code::selectron_four,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               5, /* draw it in column 5 on the mass plot. */
                               &color_enterer, /* to be drawn in green. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add 5th-lightest charged slepton:
    name_enterer.assign( "selectron_five" );
    antiname_enterer.assign( "antiselectron_five" );
    LaTeX_enterer.assign( "${\\tilde{e}}_{5}$" );
    antiLaTeX_enterer.assign( "${\\tilde{e}}_{5}^{{\\ast}}$" );
    selectron_five_pointer =
    add_particle_property_set( PDG_code::selectron_five,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               5, /* draw it in column 5 on the mass plot. */
                               &color_enterer, /* to be drawn in green. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add 6th-lightest charged slepton:
    name_enterer.assign( "selectron_six" );
    antiname_enterer.assign( "antiselectron_six" );
    LaTeX_enterer.assign( "${\\tilde{e}}_{6}$" );
    antiLaTeX_enterer.assign( "${\\tilde{e}}_{6}^{{\\ast}}$" );
    selectron_six_pointer =
    add_particle_property_set( PDG_code::selectron_six,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               5, /* draw it in column 5 on the mass plot. */
                               &color_enterer, /* to be drawn in green. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    /* add lightest sneutrino (this constructor assumes that the LSP is not a
     * sneutrino); currently, this code treats neutrinos as Dirac fermions, so
     * sneutrinos also have self_conjugate_flag = false):
     */
    name_enterer.assign( "sneutrino_one" );
    antiname_enterer.assign( "antisneutrino_one" );
    LaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{1}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{1}^{{\\ast}}$" );
    color_enterer.assign( "turquoise" );
    sneutrino_one_pointer =
    add_particle_property_set( PDG_code::sneutrino_one,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               6, /* draw it in column 6 on the mass plot. */
                               &color_enterer, /* to be drawn in turquoise. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    /* add 2nd-lightest sneutrino (this constructor assumes that the LSP is
     * not a sneutrino); currently, this code treats neutrinos as Dirac
     * fermions, so sneutrinos also have self_conjugate_flag = false):
     */
    name_enterer.assign( "sneutrino_two" );
    antiname_enterer.assign( "antisneutrino_two" );
    LaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{2}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{2}^{{\\ast}}$" );
    sneutrino_two_pointer =
    add_particle_property_set( PDG_code::sneutrino_two,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               6, /* draw it in column 6 on the mass plot. */
                               &color_enterer, /* to be drawn in turquoise. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    /* add 3rd-lightest sneutrino (this constructor assumes that the LSP is
     * not a sneutrino); currently, this code treats neutrinos as Dirac
     * fermions, so sneutrinos also have self_conjugate_flag = false):
     */
    name_enterer.assign( "sneutrino_three" );
    antiname_enterer.assign( "antisneutrino_three" );
    LaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{3}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{3}^{{\\ast}}$" );
    sneutrino_three_pointer =
    add_particle_property_set( PDG_code::sneutrino_three,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               6, /* draw it in column 6 on the mass plot. */
                               &color_enterer, /* to be drawn in turquoise. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    /* add 4th-lightest sneutrino (this constructor assumes that the LSP is
     * not a sneutrino); currently, this code treats neutrinos as Dirac
     * fermions, so sneutrinos also have self_conjugate_flag = false):
     */
    name_enterer.assign( "sneutrino_four" );
    antiname_enterer.assign( "antisneutrino_four" );
    LaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{4}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{4}^{{\\ast}}$" );
    sneutrino_four_pointer =
    add_particle_property_set( PDG_code::sneutrino_four,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               6, /* draw it in column 6 on the mass plot. */
                               &color_enterer, /* to be drawn in turquoise. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    /* add 5th-lightest sneutrino (this constructor assumes that the LSP is
     * not a sneutrino); currently, this code treats neutrinos as Dirac
     * fermions, so sneutrinos also have self_conjugate_flag = false):
     */
    name_enterer.assign( "sneutrino_five" );
    antiname_enterer.assign( "antisneutrino_five" );
    LaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{5}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{5}^{{\\ast}}$" );
    sneutrino_five_pointer =
    add_particle_property_set( PDG_code::sneutrino_five,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               6, /* draw it in column 6 on the mass plot. */
                               &color_enterer, /* to be drawn in turquoise. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    /* add 6th-lightest sneutrino (this constructor assumes that the LSP is
     * not a sneutrino); currently, this code treats neutrinos as Dirac
     * fermions, so sneutrinos also have self_conjugate_flag = false):
     */
    name_enterer.assign( "sneutrino_six" );
    antiname_enterer.assign( "antisneutrino_six" );
    LaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{6}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{6}^{{\\ast}}$" );
    sneutrino_six_pointer =
    add_particle_property_set( PDG_code::sneutrino_six,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               6, /* draw it in column 6 on the mass plot. */
                               &color_enterer, /* to be drawn in turquoise. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    /* add lightest pseudoscalar sneutrino, for the case of CP-violating mixing
     * between the components, in which case the sneutrinos above are the
     * scalar components.
     */
    name_enterer.assign( "sneutrino_pseudoscalar_one" );
    antiname_enterer.assign( "sneutrino_pseudoscalar_one" );
    LaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{P1}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{P1}$" );
    color_enterer.assign( "darkgreen" );
    add_particle_property_set( PDG_code::sneutrino_pseudoscalar_one,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               6, /* draw it in column 6 on the mass plot. */
                               &color_enterer, /* to be drawn in darkgreen. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               true /* is counted as self-conjugate. */);

    /* add 2nd-lightest pseudoscalar sneutrino, for the case of CP-violating
     * mixing between the components, in which case the sneutrinos above are
     * the scalar components.
     */
    name_enterer.assign( "sneutrino_pseudoscalar_two" );
    antiname_enterer.assign( "sneutrino_pseudoscalar_two" );
    LaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{P2}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{P2}$" );
    add_particle_property_set( PDG_code::sneutrino_pseudoscalar_two,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               6, /* draw it in column 6 on the mass plot. */
                               &color_enterer, /* to be drawn in darkgreen. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               true /* is counted as self-conjugate. */);

    /* add 3rd-lightest pseudoscalar sneutrino, for the case of CP-violating
     * mixing between the components, in which case the sneutrinos above are
     * the scalar components.
     */
    name_enterer.assign( "sneutrino_pseudoscalar_three" );
    antiname_enterer.assign( "sneutrino_pseudoscalar_three" );
    LaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{P3}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{P3}$" );
    add_particle_property_set( PDG_code::sneutrino_pseudoscalar_three,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               6, /* draw it in column 6 on the mass plot. */
                               &color_enterer, /* to be drawn in darkgreen. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               true /* is counted as self-conjugate. */);

    /* add 4th-lightest pseudoscalar sneutrino, for the case of CP-violating
     * mixing between the components, in which case the sneutrinos above are
     * the scalar components.
     */
    name_enterer.assign( "sneutrino_pseudoscalar_four" );
    antiname_enterer.assign( "sneutrino_pseudoscalar_four" );
    LaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{P4}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{P4}$" );
    add_particle_property_set( PDG_code::sneutrino_pseudoscalar_four,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               6, /* draw it in column 6 on the mass plot. */
                               &color_enterer, /* to be drawn in darkgreen. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               true /* is counted as self-conjugate. */);

    /* add 5th-lightest pseudoscalar sneutrino, for the case of CP-violating
     * mixing between the components, in which case the sneutrinos above are
     * the scalar components.
     */
    name_enterer.assign( "sneutrino_pseudoscalar_five" );
    antiname_enterer.assign( "sneutrino_pseudoscalar_five" );
    LaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{P5}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{P5}$" );
    add_particle_property_set( PDG_code::sneutrino_pseudoscalar_five,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               6, /* draw it in column 6 on the mass plot. */
                               &color_enterer, /* to be drawn in darkgreen. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               true /* is counted as self-conjugate. */);

    /* add 6th-lightest pseudoscalar sneutrino, for the case of CP-violating
     * mixing between the components, in which case the sneutrinos above are
     * the scalar components.
     */
    name_enterer.assign( "sneutrino_pseudoscalar_six" );
    antiname_enterer.assign( "sneutrino_pseudoscalar_six" );
    LaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{P6}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\nu}}}_{P6}$" );
    add_particle_property_set( PDG_code::sneutrino_pseudoscalar_six,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               6, /* draw it in column 6 on the mass plot. */
                               &color_enterer, /* to be drawn in darkgreen. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               true /* is counted as self-conjugate. */);

    // add gluino:
    name_enterer.assign( "gluino" );
    antiname_enterer.assign( "gluino" );
    LaTeX_enterer.assign( "${\\tilde{g}}$" );
    antiLaTeX_enterer.assign( "${\\tilde{g}}$" );
    color_enterer.assign( "purple" );
    gluino_pointer =
    add_particle_property_set( PDG_code::gluino,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               4, /* draw it in column 4 on the mass plot. */
                               &color_enterer, /* to be drawn in purple. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               true /* is counted as self-conjugate. */);
    // no known mass as yet.

    // add lightest neutralino (this constructor assumes that the LSP is the
    // lightest neutralino):
    name_enterer.assign( "neutralino_1" );
    antiname_enterer.assign( "neutralino_1" );
    LaTeX_enterer.assign( "${\\tilde{{\\chi}}}_{1}^{0}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\chi}}}_{1}^{0}$" );
    color_enterer.assign( "orange" );
    neutralino_one_pointer =
    add_particle_property_set( PDG_code::neutralino_one,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               3, /* draw it in column 3 on the mass plot. */
                               &color_enterer, /* to be drawn in orange. */
                               true, /* is counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               true /* is counted as self-conjugate. */);
    // no known mass as yet.

    // add 2nd-lightest neutralino:
    name_enterer.assign( "neutralino_2" );
    antiname_enterer.assign( "neutralino_2" );
    LaTeX_enterer.assign( "${\\tilde{{\\chi}}}_{2}^{0}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\chi}}}_{2}^{0}$" );
    neutralino_two_pointer =
    add_particle_property_set( PDG_code::neutralino_two,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               3, /* draw it in column 3 on the mass plot. */
                               &color_enterer, /* to be drawn in orange. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               true /* is counted as self-conjugate. */);
    // no known mass as yet.

    // add 3rd-lightest neutralino:
    name_enterer.assign( "neutralino_3" );
    antiname_enterer.assign( "neutralino_3" );
    LaTeX_enterer.assign( "${\\tilde{{\\chi}}}_{3}^{0}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\chi}}}_{3}^{0}$" );
    neutralino_three_pointer =
    add_particle_property_set( PDG_code::neutralino_three,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               3, /* draw it in column 3 on the mass plot. */
                               &color_enterer, /* to be drawn in orange. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               true /* is counted as self-conjugate. */);
    // no known mass as yet.

    // add 4th-lightest (i.e heaviest in the MSSM) neutralino:
    name_enterer.assign( "neutralino_4" );
    antiname_enterer.assign( "neutralino_4" );
    LaTeX_enterer.assign( "${\\tilde{{\\chi}}}_{4}^{0}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\chi}}}_{4}^{0}$" );
    neutralino_four_pointer =
    add_particle_property_set( PDG_code::neutralino_four,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               3, /* draw it in column 3 on the mass plot. */
                               &color_enterer, /* to be drawn in orange. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               true /* is counted as self-conjugate. */);
    // no known mass as yet.

    // add lightest chargino:
    name_enterer.assign( "positive_chargino_1" );
    antiname_enterer.assign( "negative_chargino_1" );
    LaTeX_enterer.assign( "${\\tilde{{\\chi}}}_{1}^{+}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\chi}}}_{1}^{-}$" );
    color_enterer.assign( "red" );
    chargino_one_pointer =
    add_particle_property_set( PDG_code::chargino_one,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               4, /* draw it in column 4 on the mass plot. */
                               &color_enterer, /* to be drawn in red. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

    // add 2nd-lightest (i.e heaviest in the MSSM) chargino:
    name_enterer.assign( "positive_chargino_2" );
    antiname_enterer.assign( "negative_chargino_2" );
    LaTeX_enterer.assign( "${\\tilde{{\\chi}}}_{2}^{+}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\chi}}}_{2}^{-}$" );
    chargino_two_pointer =
    add_particle_property_set( PDG_code::chargino_two,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               4, /* draw it in column 4 on the mass plot. */
                               &color_enterer, /* to be drawn in red. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               false /* not counted as self-conjugate. */);
    // no known mass as yet.

  }

  EW_scale_FVMSSM_spectrum::~EW_scale_FVMSSM_spectrum()
  {

    // does nothing.

  }



  // constructor (this constructs a EW_scale_MSSM_spectrum & then adds the
  // extra NMSSM particles):
  EW_scale_NMSSM_spectrum::EW_scale_NMSSM_spectrum() :
    EW_scale_MSSM_spectrum()
  {

    /* when adding new particles, either remember their particle code as a
     * number, or add them to the static ints in the PDG_code class in
     * CppSLHA_particle_codes.hpp, or, more conveniently but maybe a bit more
     * error-prone, given the global nature of #define, you can do something
     * like:
     */
    //#define PDG_singlet_scalar_CODE 45
    //#define PDG_singlet_pseudoscalar_CODE 46
    //#define PDG_neutralino_five_CODE 1000045
    /* however, the NMSSM particle codes are already in the PDG_code class, so
     * I am not redefining them here.
     */

    std::string name_enterer;
    std::string antiname_enterer;
    std::string LaTeX_enterer;
    std::string antiLaTeX_enterer;
    std::string color_enterer;

    /* 1st we add an extra EWSB scalar. note that if the scalar sector has a
     * state that is mainly singlet, it is unlikely to be intermediate in mass
     * to the states that are mainly not singlet, but all that matters is that
     * now there need to be 3 scalar bosons, & 2 have already had property
     * sets allocated to them, so another needs to be allocated. the actual
     * proportions of flavor eigenstates can be implemented by adding mixing
     * matrices to property sets, but that is not yet implemented.
     */
    name_enterer.assign( "heavy_scalar_Higgs_boson" );
    antiname_enterer.assign( "heavy_scalar_Higgs_boson" );
    LaTeX_enterer.assign( "$h_{3}^{0}$" );
    antiLaTeX_enterer.assign( "$h_{3}^{0}$" );
    color_enterer.assign( "gold" );
    singlet_scalar_pointer =
    add_particle_property_set( PDG_code::singlet_scalar,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               1, /* draw it in column 1 on the mass plot. */
                               &color_enterer, /* to be drawn in gold. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               true /* is counted as self-conjugate. */);
    // no known mass as yet.

    // rename the MSSM EWSB scalars:
    LaTeX_enterer.assign( "$h_{1}^{0}$" );
    antiLaTeX_enterer.assign( "$h_{1}^{0}$" );
    get_particle_property_set(
        PDG_code::light_neutral_EWSB_scalar )->set_LaTeX_name( LaTeX_enterer );
    get_particle_property_set(
        PDG_code::light_neutral_EWSB_scalar )->set_LaTeX_antiname(
                                                           antiLaTeX_enterer );
    LaTeX_enterer.assign( "$h_{2}^{0}$" );
    antiLaTeX_enterer.assign( "$h_{2}^{0}$" );
    get_particle_property_set(
        PDG_code::heavy_neutral_EWSB_scalar )->set_LaTeX_name( LaTeX_enterer );
    get_particle_property_set(
        PDG_code::heavy_neutral_EWSB_scalar )->set_LaTeX_antiname(
                                                           antiLaTeX_enterer );

    // add an extra EWSB pseudoscalar (see above for comment on the
    // irrelevance of flavor eigenstates):
    name_enterer.assign( "heavy_pseudoscalar_Higgs_boson" );
    antiname_enterer.assign( "heavy_pseudoscalar_Higgs_boson" );
    LaTeX_enterer.assign( "$a_{2}^{0}$" );
    antiLaTeX_enterer.assign( "$a_{2}^{0}$" );
    singlet_pseudoscalar_pointer =
    add_particle_property_set( PDG_code::singlet_pseudoscalar,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               1, /* draw it in column 1 on the mass plot. */
                               &color_enterer, /* to be drawn in gold. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               true /* is counted as self-conjugate. */);
    // no known mass as yet.

    // rename the MSSM EWSB pseudoscalar:
    name_enterer.assign( "light_pseudoscalar_Higgs_boson" );
    antiname_enterer.assign( "light_pseudoscalar_Higgs_boson" );
    LaTeX_enterer.assign( "$a_{1}^{0}$" );
    antiLaTeX_enterer.assign( "$a_{1}^{0}$" );
    get_particle_property_set(
      PDG_code::neutral_EWSB_pseudoscalar )->set_name( name_enterer );
    get_particle_property_set(
      PDG_code::neutral_EWSB_pseudoscalar )->set_antiname( antiname_enterer );
    get_particle_property_set(
      PDG_code::neutral_EWSB_pseudoscalar )->set_LaTeX_name( LaTeX_enterer );
    get_particle_property_set(
      PDG_code::neutral_EWSB_pseudoscalar )->set_LaTeX_antiname(
                                                           antiLaTeX_enterer );

    /* add a fifth neutralino (since with the addition of the singlino the
     * neutralino sector has 5 neutralinos) (see above for comment on the
     * irrelevance of flavor eigenstates):
     */
    name_enterer.assign( "neutralino_5" );
    antiname_enterer.assign( "neutralino_5" );
    LaTeX_enterer.assign( "${\\tilde{{\\chi}}}_{5}^{0}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\chi}}}_{5}^{0}$" );
    color_enterer.assign( "orange" );
    neutralino_five_pointer =
    add_particle_property_set( PDG_code::neutralino_five,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               3, /* draw it in column 3 on the mass plot. */
                               &color_enterer, /* to be drawn in orange. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               true /* is counted as self-conjugate. */);
    // no known mass as yet.

  }

  EW_scale_NMSSM_spectrum::~EW_scale_NMSSM_spectrum()
  {

    // does nothing.

  }



  // constructor (this constructs a EW_scale_FVMSSM_spectrum & then adds the
  // extra FVNMSSM particles):
  EW_scale_FVNMSSM_spectrum::EW_scale_FVNMSSM_spectrum() :
    EW_scale_FVMSSM_spectrum()
  {

    std::string name_enterer;
    std::string antiname_enterer;
    std::string LaTeX_enterer;
    std::string antiLaTeX_enterer;
    std::string color_enterer;

    /* 1st we add an extra EWSB scalar. note that if the scalar sector has a
     * state that is mainly singlet, it is unlikely to be intermediate in mass
     * to the states that are mainly not singlet, but all that matters is that
     * now there need to be 3 scalar bosons, & 2 have already had property
     * sets allocated to them, so another needs to be allocated. the actual
     * proportions of flavor eigenstates can be implemented by adding mixing
     * matrices to property sets, but that is not yet implemented.
     */
    name_enterer.assign( "heavy_scalar_Higgs_boson" );
    antiname_enterer.assign( "heavy_scalar_Higgs_boson" );
    LaTeX_enterer.assign( "$h_{3}^{0}$" );
    antiLaTeX_enterer.assign( "$h_{3}^{0}$" );
    color_enterer.assign( "gold" );
    singlet_scalar_pointer =
    add_particle_property_set( PDG_code::singlet_scalar,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               1, /* draw it in column 1 on the mass plot. */
                               &color_enterer, /* to be drawn in gold. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               true /* is counted as self-conjugate. */);
    // no known mass as yet.

    // rename the MSSM EWSB scalars:
    LaTeX_enterer.assign( "$h_{1}^{0}$" );
    antiLaTeX_enterer.assign( "$h_{1}^{0}$" );
    get_particle_property_set(
        PDG_code::light_neutral_EWSB_scalar )->set_LaTeX_name( LaTeX_enterer );
    get_particle_property_set(
        PDG_code::light_neutral_EWSB_scalar )->set_LaTeX_antiname(
                                                           antiLaTeX_enterer );
    LaTeX_enterer.assign( "$h_{2}^{0}$" );
    antiLaTeX_enterer.assign( "$h_{2}^{0}$" );
    get_particle_property_set(
        PDG_code::heavy_neutral_EWSB_scalar )->set_LaTeX_name( LaTeX_enterer );
    get_particle_property_set(
        PDG_code::heavy_neutral_EWSB_scalar )->set_LaTeX_antiname(
                                                           antiLaTeX_enterer );

    // add an extra EWSB pseudoscalar (see above for comment on the
    // irrelevance of flavor eigenstates):
    name_enterer.assign( "heavy_pseudoscalar_Higgs_boson" );
    antiname_enterer.assign( "heavy_pseudoscalar_Higgs_boson" );
    LaTeX_enterer.assign( "$a_{2}^{0}$" );
    antiLaTeX_enterer.assign( "$a_{2}^{0}$" );
    singlet_pseudoscalar_pointer =
    add_particle_property_set( PDG_code::singlet_pseudoscalar,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               1, /* draw it in column 1 on the mass plot. */
                               &color_enterer, /* to be drawn in gold. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               true /* is counted as self-conjugate. */);
    // no known mass as yet.

    // rename the MSSM EWSB pseudoscalar:
    name_enterer.assign( "light_pseudoscalar_Higgs_boson" );
    antiname_enterer.assign( "light_pseudoscalar_Higgs_boson" );
    LaTeX_enterer.assign( "$a_{1}^{0}$" );
    antiLaTeX_enterer.assign( "$a_{1}^{0}$" );
    get_particle_property_set(
      PDG_code::neutral_EWSB_pseudoscalar )->set_name( name_enterer );
    get_particle_property_set(
      PDG_code::neutral_EWSB_pseudoscalar )->set_antiname( antiname_enterer );
    get_particle_property_set(
      PDG_code::neutral_EWSB_pseudoscalar )->set_LaTeX_name( LaTeX_enterer );
    get_particle_property_set(
      PDG_code::neutral_EWSB_pseudoscalar )->set_LaTeX_antiname(
                                                           antiLaTeX_enterer );

    /* add a fifth neutralino (since with the addition of the singlino the
     * neutralino sector has 5 neutralinos) (see above for comment on the
     * irrelevance of flavor eigenstates):
     */
    name_enterer.assign( "neutralino_5" );
    antiname_enterer.assign( "neutralino_5" );
    LaTeX_enterer.assign( "${\\tilde{{\\chi}}}_{5}^{0}$" );
    antiLaTeX_enterer.assign( "${\\tilde{{\\chi}}}_{5}^{0}$" );
    color_enterer.assign( "orange" );
    neutralino_five_pointer =
    add_particle_property_set( PDG_code::neutralino_five,
                               &name_enterer,
                               &antiname_enterer,
                               &LaTeX_enterer,
                               &antiLaTeX_enterer,
                               3, /* draw it in column 3 on the mass plot. */
                               &color_enterer, /* to be drawn in orange. */
                               false, /* not counted as stable-&-invisible. */
                               false, /* not counted as a jet. */
                               false, /* not counted as a light lepton. */
                               true /* is counted as self-conjugate. */);
    // no known mass as yet.

  }

  EW_scale_FVNMSSM_spectrum::~EW_scale_FVNMSSM_spectrum()
  {

    // does nothing.

  }

} // end of CppSLHA namespace.
