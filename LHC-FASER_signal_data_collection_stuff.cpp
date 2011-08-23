/*
 * LHC-FASER_signal_data_collection_stuff.cpp
 *
 *  Created on: 21 Mar 2011
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

#include "LHC-FASER_signal_data_collection_stuff.hpp"

namespace LHC_FASER
{

  sQCD_to_EWino::sQCD_to_EWino(
                signed_particle_shortcut_pair const* const given_scolored_pair,
                                input_handler const* const given_shortcuts,
                               std::list< int > const* const given_jets_to_tag,
                         std::list< int > const* const given_jets_to_reject ) :
    getsReadiedForNewPoint( given_shortcuts->get_readier() ),
    shortcut( given_shortcuts ),
    scolored_pair( given_scolored_pair ),
    jets_to_tag( given_jets_to_tag ),
    jets_to_reject( given_jets_to_reject ),
    first_scolored_cascade( given_shortcuts->get_colored_cascade_set(
                                  given_scolored_pair->get_first_pointer() ) ),
    second_scolored_cascade( given_shortcuts->get_colored_cascade_set(
                                  given_scolored_pair->get_second_pointer() ) )
  {

    if( scolored_pair->get_first_pointer()->get_PDG_code()
        != CppSLHA::PDG_code::gluino )
      {

        first_scolored_to_gluino_BR
        = shortcut->get_exclusive_BR_calculator(
                                            scolored_pair->get_first_pointer(),
                                                 shortcut->get_gluino(),
                                                 true,
                                                 given_jets_to_reject );

      }
    else
      {

        first_scolored_to_gluino_BR = NULL;

      }

    if( scolored_pair->get_second_pointer()->get_PDG_code()
        != CppSLHA::PDG_code::gluino )
      {

        second_scolored_to_gluino_BR
        = shortcut->get_exclusive_BR_calculator(
                                           scolored_pair->get_second_pointer(),
                                                 shortcut->get_gluino(),
                                                 true,
                                                 given_jets_to_reject );

      }
    else
      {

        second_scolored_to_gluino_BR = NULL;

      }

    std::vector< exclusive_BR_calculator* >* scolored_to_EWino_vector;
    std::vector< exclusive_BR_calculator* >* scolored_to_antiEWino_vector;

    for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
         scolored_iterator = shortcut->get_colored_sparticles()->begin();
         shortcut->get_colored_sparticles()->end() > scolored_iterator;
         ++scolored_iterator )
      {

        // gluinos can decay into squarks.
        gluino_to_scolored_BRs.push_back(
                                         shortcut->get_exclusive_BR_calculator(
                                                        shortcut->get_gluino(),
                                                            *scolored_iterator,
                                                                          true,
                                                      given_jets_to_reject ) );
        /* even if get_exclusive_BR_calculator returns NULL, we push it into
         * the vector because we need a 1-to-1 correspondence with
         * shortcut->get_colored_sparticles(). also because of needing this
         * correspondence, we don't explicitly deal with the
         * gluino->antisquark decays. they are assumed to be the same as the
         * gluino->squark decays, because we assume no CP violation in this
         * code.
         */

        scolored_to_EWino_vector
        = new std::vector< exclusive_BR_calculator* >();
        scolored_to_antiEWino_vector
        = new std::vector< exclusive_BR_calculator* >();

        if( (*scolored_iterator)->get_PDG_code()
            == scolored_pair->get_first_pointer()->get_PDG_code() )
          // if the scolored particle we're looking at also happens to be in
          // scoloredPair...
          {

            first_scolored_to_EWino_BRs = scolored_to_EWino_vector;
            first_scolored_to_antiEWino_BRs = scolored_to_antiEWino_vector;
            // note that the BRs of the scoloredPair are to be found at the
            // same place.

          }

        if( (*scolored_iterator)->get_PDG_code()
            == scolored_pair->get_second_pointer()->get_PDG_code() )
          // if the scolored particle we're looking at also happens to be in
          // scoloredPair...
          {

            second_scolored_to_EWino_BRs = scolored_to_EWino_vector;
            second_scolored_to_antiEWino_BRs = scolored_to_antiEWino_vector;
            // note that the BRs of the scoloredPair are to be found at the
            // same place.

          }

        if( CppSLHA::PDG_code::gluino == (*scolored_iterator)->get_PDG_code() )
          // if the scolored particle we're looking at also happens to be the
          // gluino...
          {

            gluino_to_EWino_BRs = scolored_to_EWino_vector;
            // note that the BRs of the gluino are to be found at the same
            // place.

          }

        for( std::vector< CppSLHA::particle_property_set const*
                                                              >::const_iterator
             EWino_iterator = shortcut->get_electroweakinos()->begin();
             shortcut->get_electroweakinos()->end() > EWino_iterator;
             ++EWino_iterator )
          {

            scolored_to_EWino_vector->push_back(
                     shortcut->get_exclusive_BR_calculator( *scolored_iterator,
                                                            *EWino_iterator,
                                                            true,
                                                      given_jets_to_reject ) );

            if( (*EWino_iterator)->counts_as_self_conjugate()
                ||
                (*scolored_iterator)->counts_as_self_conjugate() )
              // if the antielectroweakino or *scolored_iterator is itself
              // (i.e. *scolored_iterator is the gluino)...
              {

                // don't overcount:
                scolored_to_antiEWino_vector->push_back( NULL );

              }
            else
              {

                scolored_to_antiEWino_vector->push_back(
                     shortcut->get_exclusive_BR_calculator( *scolored_iterator,
                                                            *EWino_iterator,
                                                            false,
                                                      given_jets_to_reject ) );

              }

          }  // end of loop over electroweakinos.

        scolored_to_EWino_BRs.push_back( scolored_to_EWino_vector );
        scolored_to_antiEWino_BRs.push_back( scolored_to_antiEWino_vector );

      }  // end of loop over colored sparticles.

    if( ( CppSLHA::PDG_code::gluino
          == given_scolored_pair->get_first_pointer()->get_PDG_code() )
        &&
        ( CppSLHA::PDG_code::gluino
          == given_scolored_pair->get_second_pointer()->get_PDG_code() ) )
      {

        channel_type = gg;

      }
    else if( ( CppSLHA::PDG_code::gluino
               != given_scolored_pair->get_first_pointer()->get_PDG_code() )
             &&
             ( CppSLHA::PDG_code::gluino
               == given_scolored_pair->get_second_pointer()->get_PDG_code() ) )
      {

        channel_type = sg;

      }
    else if( ( CppSLHA::PDG_code::gluino
               == given_scolored_pair->get_first_pointer()->get_PDG_code() )
             &&
             ( CppSLHA::PDG_code::gluino
               != given_scolored_pair->get_second_pointer()->get_PDG_code() ) )
      {

        channel_type = sg;

      }
    else
      {

        channel_type = ss;

      }

  }

  sQCD_to_EWino::~sQCD_to_EWino()
  {

    for( std::vector< scoloreds_to_EWinos* >::iterator
         deletion_iterator = sQCD_end.begin();
         sQCD_end.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

    for( std::vector< std::vector< exclusive_BR_calculator* >* >::iterator
         deletion_iterator = scolored_to_EWino_BRs.begin();
         scolored_to_EWino_BRs.end()  > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

    for( std::vector< std::vector< exclusive_BR_calculator* >* >::iterator
         deletion_iterator = scolored_to_antiEWino_BRs.begin();
         scolored_to_antiEWino_BRs.end()  > deletion_iterator;
         ++deletion_iterator )
     {

       delete *deletion_iterator;

     }

  }


  void
  sQCD_to_EWino::add_EWinos_to_sQCD_end(
                                    colored_cascade const* const first_cascade,
                             bool const first_end_scolored_is_not_antiparticle,
                                    bool const first_EWino_is_not_antiparticle,
                                   colored_cascade const* const second_cascade,
                                   double const second_by_scolored_to_EWino_BR,
                               double const second_by_scolored_to_antiEWino_BR,
                               double const second_by_antiscolored_to_EWino_BR,
                           double const second_by_antiscolored_to_antiEWino_BR,
                                          double const first_sQCD_BR )
  // this adds the appropriate combinations of squarks &/or antisquarks with
  // electroweakinos &/or antielectroweakinos to sQCD_end.
  {

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: sQCD_to_EWino::add_EWinos_to_sQCD_end() called."
    << std::endl
    << "1st cascade:"
    << *(first_cascade->get_decayer()->get_name_or_antiname(
                                  scoloredPair->first_is_not_antiparticle() ))
    << "->"
    << *(first_cascade->get_product()->get_name_or_antiname(
                                      first_end_scolored_is_not_antiparticle ))
    << "->"
    << *(first_cascade->get_EWino()->get_name_or_antiname(
                                             first_EWino_is_not_antiparticle ))
    << std::endl
    << "2nd cascade:"
    << *(second_cascade->get_decayer()->get_name_or_antiname(
                                 scoloredPair->second_is_not_antiparticle() ))
    << "->(anti)"
    << *(second_cascade->get_product()->getName())
    << "->(anti)"
    << *(second_cascade->get_EWino()->getName())
    << std::endl
    << "BRs: 1st: " << first_sQCD_BR
    << ", 11: " << second_by_scolored_to_EWino_BR
    << ", 10: " << second_by_scolored_to_antiEWino_BR
    << ", 01: " << second_by_antiscolored_to_EWino_BR
    << ", 00: " << second_by_antiscolored_to_antiEWino_BR;
    std::cout << std::endl;**/

    if( lhcFaserGlobal::negligibleBr
        < ( first_sQCD_BR * second_by_scolored_to_EWino_BR ) )
      // if the squark to electroweakino combination is not negligible...
      {

        // debugging:
        /**std::cout
        << std::endl
        << "adding "
        << first_end_scolored_is_not_antiparticle
        << first_EWino_is_not_antiparticle
        << ",11";
        std::cout << std::endl;**/

        sQCD_end.push_back( new scoloreds_to_EWinos( scolored_pair,
                                                     first_cascade,
                                        first_end_scolored_is_not_antiparticle,
                                               first_EWino_is_not_antiparticle,
                                                     second_cascade,
                                                     true,
                                                     true,
                                                     ( first_sQCD_BR
                                        * second_by_scolored_to_EWino_BR ) ) );

      }

    if( lhcFaserGlobal::negligibleBr
        < ( first_sQCD_BR * second_by_scolored_to_antiEWino_BR ) )
      // if the squark to antielectroweakino combination is not negligible...
      {

        // debugging:
        /**std::cout
        << std::endl
        << "adding "
        << first_end_scolored_is_not_antiparticle
        << first_EWino_is_not_antiparticle
        << ",10";
        std::cout << std::endl;**/

        sQCD_end.push_back( new scoloreds_to_EWinos( scolored_pair,
                                                     first_cascade,
                                        first_end_scolored_is_not_antiparticle,
                                               first_EWino_is_not_antiparticle,
                                                     second_cascade,
                                                     true,
                                                     false,
                                                     ( first_sQCD_BR
                                    * second_by_scolored_to_antiEWino_BR ) ) );

      }


    if( lhcFaserGlobal::negligibleBr
        < ( first_sQCD_BR * second_by_antiscolored_to_EWino_BR ) )
      // if the antisquark to electroweakino combination is not negligible...
      {

        // debugging:
        /**std::cout
        << std::endl
        << "adding "
        << first_end_scolored_is_not_antiparticle
        << first_EWino_is_not_antiparticle
        << ",01";
        std::cout << std::endl;**/

        sQCD_end.push_back( new scoloreds_to_EWinos( scolored_pair,
                                                     first_cascade,
                                        first_end_scolored_is_not_antiparticle,
                                               first_EWino_is_not_antiparticle,
                                                     second_cascade,
                                                     false,
                                                     true,
                                                     ( first_sQCD_BR
                                    * second_by_antiscolored_to_EWino_BR ) ) );

      }


    if( lhcFaserGlobal::negligibleBr
        < ( first_sQCD_BR * second_by_antiscolored_to_antiEWino_BR ) )
      // if the antisquark to antielectroweakino combination is not
      // negligible...
      {

        // debugging:
        /**std::cout
        << std::endl
        << "adding "
        << first_end_scolored_is_not_antiparticle
        << first_EWino_is_not_antiparticle
        << ",00";
        std::cout << std::endl;**/

        sQCD_end.push_back( new scoloreds_to_EWinos( scolored_pair,
                                                     first_cascade,
                                        first_end_scolored_is_not_antiparticle,
                                               first_EWino_is_not_antiparticle,
                                                     second_cascade,
                                                     false,
                                                     false,
                                                     ( first_sQCD_BR
                                * second_by_antiscolored_to_antiEWino_BR ) ) );

      }

  }


  bool
  sQCD_to_EWino::cascade_is_not_negligible(
                                          colored_cascade const* const cascade,
                                       bool const scolored_is_not_antiparticle,
                                            double* const sc_ew_BR,
                                            double* const sc_aew_BR,
                                            double* const asc_ew_BR,
                                            double* const asc_aew_BR,
                              exclusive_BR_calculator* const scolored_to_EWino,
                          exclusive_BR_calculator* const scolored_to_antiEWino,
                                exclusive_BR_calculator* const gluino_to_EWino,
                             exclusive_BR_calculator* const scolored_to_gluino,
                              exclusive_BR_calculator* const gluino_to_product,
                               exclusive_BR_calculator* const product_to_EWino,
                          exclusive_BR_calculator* const product_to_antiEWino )
  /* this checks to see if the colored_cascade is not negligible given the
   * provided branching ratio calculators. it returns true if at least 1 of the
   * final possibilities (whether the scolored before the electroweakino is a
   * particle or antiparticle, & whether the electroweakino is a particle or
   * antiparticle) is not negligible, as well as assigning the relevant values
   * to the given doubles.
   */
  {

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: sQCD_to_EWino::cascade_is_not_negligible() called.";
    std::cout << std::endl;**/

    // set the BRs to 0.0 by default:
    *sc_ew_BR = 0.0;
    *sc_aew_BR = 0.0;
    *asc_ew_BR = 0.0;
    *asc_aew_BR = 0.0;

    bool return_value = false;
    double BR_to_record;

    if( cascade->is_sx() )
      // if it's a direct squark decay to the electroweakino...
      {

        // debugging:
        /**std::cout
        << std::endl
        << "cascade is sx.";
        std::cout << std::endl;**/

        if( NULL != scolored_to_EWino )
          {

            BR_to_record = scolored_to_EWino->get_BR();

            if( lhcFaserGlobal::negligibleBr < BR_to_record )
              {

                if( scolored_is_not_antiparticle )
                  {

                    *sc_ew_BR = BR_to_record;

                  }
                else
                  {

                    if( cascade->get_EWino()->counts_as_self_conjugate() )
                      {

                        *asc_ew_BR = BR_to_record;

                      }
                    else
                      {

                        *asc_aew_BR = BR_to_record;

                      }

                  }

                return_value = true;

              }

          }

        if( NULL != scolored_to_antiEWino )
          // neutralinos should never pass this conditional...
          {

            BR_to_record = scolored_to_antiEWino->get_BR();

            if( lhcFaserGlobal::negligibleBr < BR_to_record )
              {

                if( scolored_is_not_antiparticle )
                  {

                    *sc_aew_BR = BR_to_record;

                  }
                else
                  {

                    *asc_ew_BR = BR_to_record;

                  }

                return_value = true;

              }

          }

      }
    else if( cascade->is_gx() )
      // if it's a direct gluino decay to the electroweakino...
      {

        // debugging:
        /**std::cout
        << std::endl
        << "cascade is gx.";
        std::cout << std::endl;**/

        if( NULL != scolored_to_EWino )
          {

            BR_to_record = scolored_to_EWino->get_BR();

            if( lhcFaserGlobal::negligibleBr < BR_to_record )
              {

                *sc_ew_BR = BR_to_record;
                if( !(cascade->get_EWino()->counts_as_self_conjugate()) )
                  {

                    *sc_aew_BR = BR_to_record;

                  }
                /* we take all gluinos as particles rather than allowing for
                 * antigluinos. also, we assume that the charge conjugate decay
                 * is the same.
                 */

                return_value = true;

                // this function should set gluinos to have their
                // "antiscolored_to_(anti)EWino" value set to 0.0 always.

              }

          }

      }
    else if( cascade->is_sgx() )
      // otherwise if it's a squark 2-stage decay to the electroweakino...
      {

        // debugging:
        /**std::cout
        << std::endl
        << "cascade is sgx.";
        std::cout << std::endl;**/

        if( ( NULL != scolored_to_gluino )
            &&
            ( NULL != gluino_to_EWino ) )
          {

            BR_to_record
            = ( scolored_to_gluino->get_BR() * gluino_to_EWino->get_BR() );

            if( lhcFaserGlobal::negligibleBr
                < BR_to_record )
              {

                *sc_ew_BR = BR_to_record;
                if( !(cascade->get_EWino()->counts_as_self_conjugate()) )
                  {

                    *sc_aew_BR = BR_to_record;

                  }
                /* we take all gluinos as particles rather than allowing for
                 * antigluinos. also, we assume that the charge conjugate decay
                 * is the same.
                 */

                return_value = true;

              }

          }

      }
    else if( cascade->is_gsx() )
      // otherwise if it's a gluino 2-stage decay to the electroweakino...
      {

        // debugging:
        /**std::cout
        << std::endl
        << "cascade is gsx; gluino_to_product = " << gluino_to_product
        << ", product_to_EWino = " << product_to_EWino;
        std::cout << std::endl;**/

        if( NULL != gluino_to_product )
          {

            if( NULL != product_to_EWino )
              {

                BR_to_record
                = ( gluino_to_product->get_BR() * product_to_EWino->get_BR() );

                // debugging:
                /**std::cout
                << std::endl
                << "gluino_to_product->get_BR() = "
                << gluino_to_product->get_BR()
                << ", product_to_EWino->get_BR() = "
                << product_to_EWino->get_BR();
                std::cout << std::endl;**/

                if( lhcFaserGlobal::negligibleBr < BR_to_record )
                  {

                    *sc_ew_BR = BR_to_record;
                    if( cascade->get_EWino()->counts_as_self_conjugate() )
                      {

                        *asc_ew_BR = BR_to_record;

                      }
                    else
                      {

                        *asc_aew_BR = BR_to_record;

                      }
                    /* we take all gluinos as particles rather than allowing
                     * for antigluinos. also, we assume that the charge
                     * conjugate decay is the same.
                     */

                    return_value = true;

                    // this function should set gluinos to have their
                    // "antiscolored_to_(anti)EWino" value set to 0.0 always.

                  }

              }

            if( NULL != product_to_antiEWino )
              // neutralinos should never pass this conditional...
              {

                BR_to_record = ( gluino_to_product->get_BR()
                                 * product_to_antiEWino->get_BR() );

                if( lhcFaserGlobal::negligibleBr < BR_to_record )
                  {

                    *sc_aew_BR = BR_to_record;
                    *asc_ew_BR = BR_to_record;
                    /* we take all gluinos as particles rather than allowing
                     * for antigluinos. also, we assume that the charge
                     * conjugate decay is the same.
                     */

                    return_value = true;

                  }

              }

          }

      }
    else
      // otherwise it's a 3-stage decay to the electroweakino (though
      // not necessarily valid)...
      {

        // debugging:
        /**std::cout
        << std::endl
        << "cascade is sgsx.";
        std::cout << std::endl;**/

        if( ( NULL != scolored_to_gluino )
            &&
            ( NULL != gluino_to_product ) )
          {

            if( NULL != product_to_EWino )
              {

                BR_to_record = ( scolored_to_gluino->get_BR()
                                 * gluino_to_product->get_BR()
                                 * product_to_EWino->get_BR() );

                if( lhcFaserGlobal::negligibleBr < BR_to_record )
                  {

                    *sc_ew_BR = BR_to_record;
                    if( !(cascade->get_EWino()->counts_as_self_conjugate()) )
                      {

                        *asc_aew_BR = BR_to_record;

                      }
                    /* we take all gluinos as particles rather than allowing
                     * for antigluinos. also, we assume that the charge
                     * conjugate decay is the same.
                     */

                    return_value = true;

                  }

              }

            if( NULL != product_to_antiEWino )
              // neutralinos should never pass this conditional...
              {

                BR_to_record = ( scolored_to_gluino->get_BR()
                                 * gluino_to_product->get_BR()
                                 * product_to_antiEWino->get_BR() );

                if( lhcFaserGlobal::negligibleBr < BR_to_record )
                  {

                    *sc_aew_BR = BR_to_record;
                    *asc_ew_BR = BR_to_record;
                    /* we take all gluinos as particles rather than allowing
                     * for antigluinos. also, we assume that the charge
                     * conjugate decay is the same.
                     */

                    return_value = true;

                  }

              }

          }

      }

    // debugging:
    /**std::cout
    << std::endl
    << "sQCD_to_EWino::cascade_is_not_negligible() returning "
    << return_value;
    std::cout << std::endl;**/

    return return_value;

  }


  void
  sQCD_to_EWino::prepare_channels()
  /* this checks to see if any of the kinematically-allowed cascade
   * combinations following from the channel survive being multiplied by their
   * branching ratios. those that do survive get put into sQCD_end for later
   * processing by the cascade handler.
   */
  {

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: sQCD_to_EWino::prepare_channels() called.";
    std::cout << std::endl;**/

    // 1st, reset sQCD_end:
    for( std::vector< scoloreds_to_EWinos* >::iterator
         deletion_iterator = sQCD_end.begin();
         sQCD_end.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }
    sQCD_end.clear();

    // debugging:
    /**std::cout
    << std::endl
    << "sQCD_end.clear() called successfully.";
    std::cout << std::endl;**/

    double first_by_scolored_to_EWino_BR;
    double first_by_scolored_to_antiEWino_BR;
    double first_by_antiscolored_to_EWino_BR;
    double first_by_antiscolored_to_antiEWino_BR;
    double second_by_scolored_to_EWino_BR;
    double second_by_scolored_to_antiEWino_BR;
    double second_by_antiscolored_to_EWino_BR;
    double second_by_antiscolored_to_antiEWino_BR;
    double identical_cascade_BR_product;


    for( std::vector< colored_cascade* >::const_iterator
         first_sQCD_iterator
         = first_scolored_cascade->get_colored_cascades()->begin();
         first_scolored_cascade->get_colored_cascades()->end()
         > first_sQCD_iterator;
         ++first_sQCD_iterator )
      {

        // debugging:
        /**std::cout
        << std::endl
        << "1st cascade is "
        << *((*first_sQCD_iterator)->get_decayer()->getName()) << " => "
        << *((*first_sQCD_iterator)->get_product()->getName()) << " => "
        << *((*first_sQCD_iterator)->get_EWino()->getName());
        std::cout << std::endl;**/

        // debugging:
        /**std::cout
        << std::endl
        << "about to call cascade_is_not_negligible() for 1st cascade.";
        std::cout << std::endl;**/

        if( cascade_is_not_negligible( *first_sQCD_iterator,
                                    scolored_pair->first_is_not_antiparticle(),
                                       &first_by_scolored_to_EWino_BR,
                                       &first_by_scolored_to_antiEWino_BR,
                                       &first_by_antiscolored_to_EWino_BR,
                                       &first_by_antiscolored_to_antiEWino_BR,
                                       first_scolored_to_EWino_BRs->at(
                                (*first_sQCD_iterator)->get_EWino_position() ),
                                       first_scolored_to_antiEWino_BRs->at(
                                (*first_sQCD_iterator)->get_EWino_position() ),
                                       gluino_to_EWino_BRs->at(
                                (*first_sQCD_iterator)->get_EWino_position() ),
                                       first_scolored_to_gluino_BR,
                                       gluino_to_scolored_BRs.at(
                              (*first_sQCD_iterator)->get_product_position() ),
                                       scolored_to_EWino_BRs.at(
                          (*first_sQCD_iterator)->get_product_position() )->at(
                                (*first_sQCD_iterator)->get_EWino_position() ),
                                       scolored_to_antiEWino_BRs.at(
                          (*first_sQCD_iterator)->get_product_position() )->at(
                             (*first_sQCD_iterator)->get_EWino_position() ) ) )
          {

            // debugging:
            /**std::cout
            << std::endl
            << "cascade_is_not_negligible() for 1st cascade true.";
            std::cout << std::endl;**/

            if( ( scolored_pair->first_is_not_antiparticle()
                  == scolored_pair->second_is_not_antiparticle() )
                &&
                ( scolored_pair->get_first_pointer()
                  == scolored_pair->get_second_pointer() ) )
              // if the 2 cascades start identically...
              {

                std::vector< colored_cascade* >::const_iterator
                second_iterator_start = first_sQCD_iterator;
                ++second_iterator_start;

                // debugging:
                /**std::cout
                << std::endl
                << "2nd cascade is "
                << *((*first_sQCD_iterator)->get_decayer()->getName())
                << " => "
                << *((*first_sQCD_iterator)->get_product()->getName())
                << " => "
                << *((*first_sQCD_iterator)->get_EWino()->getName());
                std::cout << std::endl;**/


                // we start with the identical cascades case (we know that
                // the 2nd cascade is also not neglible, as the 1st is not):

                // debugging:
                /**std::cout
                << std::endl
                << "cascade_is_not_negligible() for 2nd cascade true.";
                std::cout << std::endl;**/

                identical_cascade_BR_product
                = ( first_by_scolored_to_EWino_BR
                    * first_by_scolored_to_EWino_BR );
                if( lhcFaserGlobal::negligibleBr
                    < identical_cascade_BR_product )
                  {

                    sQCD_end.push_back( new scoloreds_to_EWinos( scolored_pair,
                                                          *first_sQCD_iterator,
                                                                 true,
                                                                 true,
                                                          *first_sQCD_iterator,
                                                                 true,
                                                                 true,
                                              identical_cascade_BR_product ) );

                  }
                //debugging:
                /**else if( ( CppSLHA::PDG_code::gluino
                           == scoloredPair->get_first_pointer(
                                                            )->get_PDG_code() )
                         &&
                         ( CppSLHA::PDG_code::gluino
                           == scoloredPair->get_second_pointer(
                                                          )->get_PDG_code() ) )
                  {

                    std::cout
                    << std::endl
                    << "++++ rejected: cascade BR = "
                    << identical_cascade_BR_product;
                    std::cout << std:: endl;

                  }**/

                // debugging:
                /**if( (*first_sQCD_iterator)->is_gx() )
                  {

                    std::cout
                    << std::endl
                    << "debugging: gx has first_by_antiscolored_to_EWino_BR = "
                    << first_by_antiscolored_to_EWino_BR
                    << ", first_by_antiscolored_to_antiEWino_BR = "
                    << first_by_antiscolored_to_antiEWino_BR;

                  }**/

                identical_cascade_BR_product
                = ( 2.0 * first_by_scolored_to_EWino_BR
                        * first_by_antiscolored_to_EWino_BR );
                // first_by_antiscolored_to_EWino_BR is 0.0 for gx.
                if( lhcFaserGlobal::negligibleBr
                    < identical_cascade_BR_product )
                  {

                    sQCD_end.push_back( new scoloreds_to_EWinos( scolored_pair,
                                                          *first_sQCD_iterator,
                                                                 true,
                                                                 true,
                                                          *first_sQCD_iterator,
                                                                 false,
                                                                 true,
                                              identical_cascade_BR_product ) );

                  }
                //debugging:
                /**else if( ( CppSLHA::PDG_code::gluino
                           == scoloredPair->get_first_pointer(
                                                            )->get_PDG_code() )
                         &&
                         ( CppSLHA::PDG_code::gluino
                           == scoloredPair->get_second_pointer(
                                                          )->get_PDG_code() ) )
                  {

                    std::cout
                    << std::endl
                    << "++-+ rejected: cascade BR = "
                    << identical_cascade_BR_product;
                    std::cout << std:: endl;

                  }**/

                identical_cascade_BR_product
                = ( first_by_antiscolored_to_EWino_BR
                    * first_by_antiscolored_to_EWino_BR );
                // first_by_antiscolored_to_EWino_BR is 0.0 for gx.
                if( lhcFaserGlobal::negligibleBr
                    < identical_cascade_BR_product )
                  {

                    sQCD_end.push_back( new scoloreds_to_EWinos( scolored_pair,
                                                          *first_sQCD_iterator,
                                                                 false,
                                                                 true,
                                                          *first_sQCD_iterator,
                                                                 false,
                                                                 true,
                                              identical_cascade_BR_product ) );

                  }
                //debugging:
                /**else if( ( CppSLHA::PDG_code::gluino
                           == scoloredPair->get_first_pointer(
                                                            )->get_PDG_code() )
                         &&
                         ( CppSLHA::PDG_code::gluino
                           == scoloredPair->get_second_pointer(
                                                          )->get_PDG_code() ) )
                  {

                    std::cout
                    << std::endl
                    << "-+-+ rejected: cascade BR = "
                    << identical_cascade_BR_product;
                    std::cout << std:: endl;

                  }**/

                if( !((*first_sQCD_iterator)->get_EWino(
                                               )->counts_as_self_conjugate()) )
                  // if the cascades end with charginos, we have to allow for
                  // both charges...
                  {

                    identical_cascade_BR_product
                    = ( 2.0 * first_by_scolored_to_EWino_BR
                            * first_by_scolored_to_antiEWino_BR );
                    if( lhcFaserGlobal::negligibleBr
                        < identical_cascade_BR_product )
                      {

                        sQCD_end.push_back( new scoloreds_to_EWinos(
                                                                 scolored_pair,
                                                          *first_sQCD_iterator,
                                                                     true,
                                                                     true,
                                                          *first_sQCD_iterator,
                                                                     true,
                                                                     false,
                                              identical_cascade_BR_product ) );

                      }
                    //debugging:
                    /**else if( ( CppSLHA::PDG_code::gluino
                               == scoloredPair->get_first_pointer(
                                                            )->get_PDG_code() )
                           &&
                           ( CppSLHA::PDG_code::gluino
                             == scoloredPair->get_second_pointer(
                                                          )->get_PDG_code() ) )
                      {

                        std::cout
                        << std::endl
                        << "+++- rejected: cascade BR = "
                        << identical_cascade_BR_product;
                        std::cout << std:: endl;

                      }**/

                    identical_cascade_BR_product
                    = ( first_by_scolored_to_antiEWino_BR
                        * first_by_scolored_to_antiEWino_BR );
                    if( lhcFaserGlobal::negligibleBr
                        < identical_cascade_BR_product )
                      {

                        sQCD_end.push_back( new scoloreds_to_EWinos(
                                                                 scolored_pair,
                                                          *first_sQCD_iterator,
                                                                     true,
                                                                     false,
                                                          *first_sQCD_iterator,
                                                                     true,
                                                                     false,
                                              identical_cascade_BR_product ) );

                      }
                    //debugging:
                    /**else if( ( CppSLHA::PDG_code::gluino
                               == scoloredPair->get_first_pointer(
                                                            )->get_PDG_code() )
                           &&
                           ( CppSLHA::PDG_code::gluino
                             == scoloredPair->get_second_pointer(
                                                          )->get_PDG_code() ) )
                      {

                        std::cout
                        << std::endl
                        << "+-+- rejected: cascade BR = "
                        << identical_cascade_BR_product;
                        std::cout << std:: endl;

                      }**/

                    identical_cascade_BR_product
                    = ( 2.0 * first_by_scolored_to_antiEWino_BR
                            * first_by_antiscolored_to_EWino_BR );
                    if( lhcFaserGlobal::negligibleBr
                        < identical_cascade_BR_product )
                      {

                        sQCD_end.push_back( new scoloreds_to_EWinos(
                                                                 scolored_pair,
                                                          *first_sQCD_iterator,
                                                                     true,
                                                                     false,
                                                          *first_sQCD_iterator,
                                                                     false,
                                                                     true,
                                              identical_cascade_BR_product ) );

                      }
                    //debugging:
                    /**else if( ( CppSLHA::PDG_code::gluino
                               == scoloredPair->get_first_pointer(
                                                            )->get_PDG_code() )
                             &&
                             ( CppSLHA::PDG_code::gluino
                               == scoloredPair->get_second_pointer(
                                                          )->get_PDG_code() ) )
                      {

                        std::cout
                        << std::endl
                        << "+--+ rejected: cascade BR = "
                        << identical_cascade_BR_product;
                        std::cout << std:: endl;

                      }**/

                    identical_cascade_BR_product
                    = ( 2.0 * first_by_scolored_to_EWino_BR
                            * first_by_antiscolored_to_antiEWino_BR );
                    if( lhcFaserGlobal::negligibleBr
                        < identical_cascade_BR_product )
                      {

                        sQCD_end.push_back( new scoloreds_to_EWinos(
                                                                 scolored_pair,
                                                          *first_sQCD_iterator,
                                                                     true,
                                                                     true,
                                                          *first_sQCD_iterator,
                                                                     false,
                                                                     false,
                                              identical_cascade_BR_product ) );

                      }
                    //debugging:
                    /**else if( ( CppSLHA::PDG_code::gluino
                               == scoloredPair->get_first_pointer(
                                                            )->get_PDG_code() )
                             &&
                             ( CppSLHA::PDG_code::gluino
                               == scoloredPair->get_second_pointer(
                                                          )->get_PDG_code() ) )
                      {

                        std::cout
                        << std::endl
                        << "++-- rejected: cascade BR = "
                        << identical_cascade_BR_product;
                        std::cout << std:: endl;

                      }**/

                    identical_cascade_BR_product
                    = ( 2.0 * first_by_scolored_to_antiEWino_BR
                            * first_by_antiscolored_to_antiEWino_BR );
                    if( lhcFaserGlobal::negligibleBr
                        < identical_cascade_BR_product )
                      {

                        sQCD_end.push_back( new scoloreds_to_EWinos(
                                                                 scolored_pair,
                                                          *first_sQCD_iterator,
                                                                     true,
                                                                     false,
                                                          *first_sQCD_iterator,
                                                                     false,
                                                                     false,
                                              identical_cascade_BR_product ) );

                      }
                    //debugging:
                    /**else if( ( CppSLHA::PDG_code::gluino
                               == scoloredPair->get_first_pointer(
                                                            )->get_PDG_code() )
                             &&
                             ( CppSLHA::PDG_code::gluino
                               == scoloredPair->get_second_pointer(
                                                          )->get_PDG_code() ) )
                      {

                        std::cout
                        << std::endl
                        << "+--- rejected: cascade BR = "
                        << identical_cascade_BR_product;
                        std::cout << std:: endl;

                      }**/

                    identical_cascade_BR_product
                    = ( 2.0 * first_by_antiscolored_to_EWino_BR
                            * first_by_antiscolored_to_antiEWino_BR );
                    if( lhcFaserGlobal::negligibleBr
                        < identical_cascade_BR_product )
                      {

                        sQCD_end.push_back( new scoloreds_to_EWinos(
                                                                 scolored_pair,
                                                          *first_sQCD_iterator,
                                                                     false,
                                                                     true,
                                                          *first_sQCD_iterator,
                                                                     false,
                                                                     false,
                                              identical_cascade_BR_product ) );

                      }
                    //debugging:
                    /**else if( ( CppSLHA::PDG_code::gluino
                               == scoloredPair->get_first_pointer(
                                                            )->get_PDG_code() )
                             &&
                             ( CppSLHA::PDG_code::gluino
                               == scoloredPair->get_second_pointer(
                                                          )->get_PDG_code() ) )
                      {

                        std::cout
                        << std::endl
                        << "-+-- rejected: cascade BR = "
                        << identical_cascade_BR_product;
                        std::cout << std:: endl;

                      }**/

                    identical_cascade_BR_product
                    = ( first_by_antiscolored_to_antiEWino_BR
                        * first_by_antiscolored_to_antiEWino_BR );
                    if( lhcFaserGlobal::negligibleBr
                        < identical_cascade_BR_product )
                      {

                        sQCD_end.push_back( new scoloreds_to_EWinos(
                                                                 scolored_pair,
                                                          *first_sQCD_iterator,
                                                                     false,
                                                                     false,
                                                          *first_sQCD_iterator,
                                                                     false,
                                                                     false,
                                              identical_cascade_BR_product ) );

                      }
                    //debugging:
                    /**else if( ( CppSLHA::PDG_code::gluino
                               == scoloredPair->get_first_pointer(
                                                            )->get_PDG_code() )
                             &&
                             ( CppSLHA::PDG_code::gluino
                               == scoloredPair->get_second_pointer(
                                                          )->get_PDG_code() ) )
                      {

                        std::cout
                        << std::endl
                        << "---- rejected: cascade BR = "
                        << identical_cascade_BR_product;
                        std::cout << std:: endl;

                      }**/

                  }


                // now we deal with the non-identical cases, doubling their
                // values to account for not adding the swapped versions later:

                first_by_scolored_to_EWino_BR
                += first_by_scolored_to_EWino_BR;
                first_by_antiscolored_to_EWino_BR
                += first_by_antiscolored_to_EWino_BR;
                first_by_scolored_to_antiEWino_BR
                += first_by_scolored_to_antiEWino_BR;
                first_by_antiscolored_to_antiEWino_BR
                += first_by_antiscolored_to_antiEWino_BR;

                for( std::vector< colored_cascade* >::const_iterator
                     second_sQCD_iterator = second_iterator_start;
                     second_scolored_cascade->get_colored_cascades()->end()
                     > second_sQCD_iterator;
                     ++second_sQCD_iterator )
                  // we avoid overcounting any previous identical decays...
                  {

                    // debugging:
                    /**std::cout
                    << std::endl
                    << "2nd cascade is "
                    << *((*second_sQCD_iterator)->get_decayer()->getName())
                    << " => "
                    << *((*second_sQCD_iterator)->get_product()->getName())
                    << " => "
                    << *((*second_sQCD_iterator)->get_EWino()->getName());
                    std::cout << std::endl;**/

                    // debugging:
                    /**std::cout
                    << std::endl
                    << "about to call cascade_is_not_negligible() for 2nd"
                    << "cascade.";
                    std::cout << std::endl;**/

                    if( cascade_is_not_negligible( *second_sQCD_iterator,
                                   scolored_pair->second_is_not_antiparticle(),
                                               &second_by_scolored_to_EWino_BR,
                                           &second_by_scolored_to_antiEWino_BR,
                                           &second_by_antiscolored_to_EWino_BR,
                                       &second_by_antiscolored_to_antiEWino_BR,
                                              second_scolored_to_EWino_BRs->at(
                               (*second_sQCD_iterator)->get_EWino_position() ),
                                          second_scolored_to_antiEWino_BRs->at(
                               (*second_sQCD_iterator)->get_EWino_position() ),
                                                   gluino_to_EWino_BRs->at(
                               (*second_sQCD_iterator)->get_EWino_position() ),
                                                  second_scolored_to_gluino_BR,
                                                   gluino_to_scolored_BRs.at(
                             (*second_sQCD_iterator)->get_product_position() ),
                                                   scolored_to_EWino_BRs.at(
                         (*second_sQCD_iterator)->get_product_position() )->at(
                               (*second_sQCD_iterator)->get_EWino_position() ),
                                                  scolored_to_antiEWino_BRs.at(
                         (*second_sQCD_iterator)->get_product_position() )->at(
                            (*second_sQCD_iterator)->get_EWino_position() ) ) )
                      {

                        // debugging:
                        /**std::cout
                        << std::endl
                        << "cascade_is_not_negligible() for 2nd cascade true.";
                        std::cout << std::endl;**/

                        add_EWinos_to_sQCD_end( *first_sQCD_iterator,
                                                true,
                                                true,
                                                *second_sQCD_iterator,
                                                second_by_scolored_to_EWino_BR,
                                            second_by_scolored_to_antiEWino_BR,
                                            second_by_antiscolored_to_EWino_BR,
                                        second_by_antiscolored_to_antiEWino_BR,
                                               first_by_scolored_to_EWino_BR );

                        add_EWinos_to_sQCD_end( *first_sQCD_iterator,
                                                false,
                                                true,
                                                *second_sQCD_iterator,
                                                second_by_scolored_to_EWino_BR,
                                            second_by_scolored_to_antiEWino_BR,
                                            second_by_antiscolored_to_EWino_BR,
                                        second_by_antiscolored_to_antiEWino_BR,
                                           first_by_antiscolored_to_EWino_BR );

                        if( !((*first_sQCD_iterator)->get_EWino(
                                               )->counts_as_self_conjugate()) )
                          // if the 1st cascade ends with a chargino, we have
                          // to allow for both charges...
                          {

                            add_EWinos_to_sQCD_end( *first_sQCD_iterator,
                                                    true,
                                                    false,
                                                    *second_sQCD_iterator,
                                                second_by_scolored_to_EWino_BR,
                                            second_by_scolored_to_antiEWino_BR,
                                            second_by_antiscolored_to_EWino_BR,
                                        second_by_antiscolored_to_antiEWino_BR,
                                           first_by_scolored_to_antiEWino_BR );

                            add_EWinos_to_sQCD_end( *first_sQCD_iterator,
                                                    false,
                                                    false,
                                                    *second_sQCD_iterator,
                                                second_by_scolored_to_EWino_BR,
                                            second_by_scolored_to_antiEWino_BR,
                                            second_by_antiscolored_to_EWino_BR,
                                        second_by_antiscolored_to_antiEWino_BR,
                                       first_by_antiscolored_to_antiEWino_BR );

                          }

                      }  // end of if 2nd sQCD cascade is not negligible.
                    //debugging:
                    /**else if( ( CppSLHA::PDG_code::gluino
                               == scoloredPair->get_first_pointer(
                                                            )->get_PDG_code() )
                             &&
                             ( CppSLHA::PDG_code::gluino
                               == scoloredPair->get_second_pointer(
                                                          )->get_PDG_code() ) )
                      {

                        std::cout
                        << std::endl
                        << "rejected: cascade BR = { "
                        << first_by_scolored_to_EWino_BR << " / "
                        << first_by_scolored_to_antiEWino_BR << " / "
                        << first_by_antiscolored_to_EWino_BR << " / "
                        << first_by_antiscolored_to_antiEWino_BR << " } * { "
                        << second_by_scolored_to_EWino_BR << " / "
                        << second_by_scolored_to_antiEWino_BR << " / "
                        << second_by_antiscolored_to_EWino_BR << " / "
                        << second_by_antiscolored_to_antiEWino_BR << " }";
                        std::cout << std:: endl;

                      }**/

                  }  // end of loop over 2nd sQCD cascade.

              } // end of if the cascades are identical
            else
              {

                for( std::vector< colored_cascade* >::const_iterator
                     second_sQCD_iterator
                     = second_scolored_cascade->get_colored_cascades(
                                                                    )->begin();
                     second_scolored_cascade->get_colored_cascades()->end()
                     > second_sQCD_iterator;
                     ++second_sQCD_iterator )
                  {

                    // debugging:
                    /**std::cout
                    << std::endl
                    << "2nd cascade is "
                    << *((*second_sQCD_iterator)->get_decayer()->getName())
                    << " => "
                    << *((*second_sQCD_iterator)->get_product()->getName())
                    << " => "
                    << *((*second_sQCD_iterator)->get_EWino()->getName());
                    std::cout << std::endl;**/

                    // debugging:
                    /**std::cout
                    << std::endl
                    << "about to call cascade_is_not_negligible() for 2nd"
                    << " cascade.";
                    std::cout << std::endl;**/

                    if( cascade_is_not_negligible( *second_sQCD_iterator,
                                   scolored_pair->second_is_not_antiparticle(),
                                               &second_by_scolored_to_EWino_BR,
                                           &second_by_scolored_to_antiEWino_BR,
                                           &second_by_antiscolored_to_EWino_BR,
                                       &second_by_antiscolored_to_antiEWino_BR,
                                              second_scolored_to_EWino_BRs->at(
                               (*second_sQCD_iterator)->get_EWino_position() ),
                                          second_scolored_to_antiEWino_BRs->at(
                               (*second_sQCD_iterator)->get_EWino_position() ),
                                                   gluino_to_EWino_BRs->at(
                               (*second_sQCD_iterator)->get_EWino_position() ),
                                                  second_scolored_to_gluino_BR,
                                                   gluino_to_scolored_BRs.at(
                             (*second_sQCD_iterator)->get_product_position() ),
                                                   scolored_to_EWino_BRs.at(
                         (*second_sQCD_iterator)->get_product_position() )->at(
                               (*second_sQCD_iterator)->get_EWino_position() ),
                                                  scolored_to_antiEWino_BRs.at(
                         (*second_sQCD_iterator)->get_product_position() )->at(
                            (*second_sQCD_iterator)->get_EWino_position() ) ) )
                      {

                        // debugging:
                        /**std::cout
                        << std::endl
                        << "cascade_is_not_negligible() for 2nd cascade true.";
                        std::cout << std::endl;**/

                        add_EWinos_to_sQCD_end( *first_sQCD_iterator,
                                                true,
                                                true,
                                                *second_sQCD_iterator,
                                                second_by_scolored_to_EWino_BR,
                                            second_by_scolored_to_antiEWino_BR,
                                            second_by_antiscolored_to_EWino_BR,
                                        second_by_antiscolored_to_antiEWino_BR,
                                               first_by_scolored_to_EWino_BR );

                        add_EWinos_to_sQCD_end( *first_sQCD_iterator,
                                                false,
                                                true,
                                                *second_sQCD_iterator,
                                                second_by_scolored_to_EWino_BR,
                                            second_by_scolored_to_antiEWino_BR,
                                            second_by_antiscolored_to_EWino_BR,
                                        second_by_antiscolored_to_antiEWino_BR,
                                           first_by_antiscolored_to_EWino_BR );

                        if( !((*first_sQCD_iterator)->get_EWino(
                                               )->counts_as_self_conjugate()) )
                          // if the 1st cascade ends with a chargino, we have
                          // to allow for both charges...
                          {

                            add_EWinos_to_sQCD_end( *first_sQCD_iterator,
                                                    true,
                                                    false,
                                                    *second_sQCD_iterator,
                                                second_by_scolored_to_EWino_BR,
                                            second_by_scolored_to_antiEWino_BR,
                                            second_by_antiscolored_to_EWino_BR,
                                        second_by_antiscolored_to_antiEWino_BR,
                                           first_by_scolored_to_antiEWino_BR );

                            add_EWinos_to_sQCD_end( *first_sQCD_iterator,
                                                    false,
                                                    false,
                                                    *second_sQCD_iterator,
                                                second_by_scolored_to_EWino_BR,
                                            second_by_scolored_to_antiEWino_BR,
                                            second_by_antiscolored_to_EWino_BR,
                                        second_by_antiscolored_to_antiEWino_BR,
                                       first_by_antiscolored_to_antiEWino_BR );

                          }

                      }  // end of if 2nd sQCD cascade is not negligible.
                    //debugging:
                    /**else if( ( CppSLHA::PDG_code::gluino
                               == scoloredPair->get_first_pointer(
                                                            )->get_PDG_code() )
                             &&
                             ( CppSLHA::PDG_code::gluino
                               == scoloredPair->get_second_pointer(
                                                          )->get_PDG_code() ) )
                      {

                        std::cout
                        << std::endl
                        << "rejected: cascade BR = { "
                        << first_by_scolored_to_EWino_BR << " / "
                        << first_by_scolored_to_antiEWino_BR << " / "
                        << first_by_antiscolored_to_EWino_BR << " / "
                        << first_by_antiscolored_to_antiEWino_BR << " } * { "
                        << second_by_scolored_to_EWino_BR << " / "
                        << second_by_scolored_to_antiEWino_BR << " / "
                        << second_by_antiscolored_to_EWino_BR << " / "
                        << second_by_antiscolored_to_antiEWino_BR << " }";
                        std::cout << std:: endl;

                      }**/

                  }  // end of loop over 2nd sQCD cascade.

              }  // end of else the cascades did not start identically.

          }  // end of if 1st sQCD cascade is not negligible.
        //debugging:
        /**else if( ( CppSLHA::PDG_code::gluino
                   == scoloredPair->get_first_pointer()->get_PDG_code() )
                 &&
                 ( CppSLHA::PDG_code::gluino
                   == scoloredPair->get_second_pointer()->get_PDG_code() ) )
          {

            std::cout
            << std::endl
            << "rejected: cascade BR = "
            << first_by_scolored_to_EWino_BR << " / "
            << first_by_scolored_to_antiEWino_BR << " / "
            << first_by_antiscolored_to_EWino_BR << " / "
            << first_by_antiscolored_to_antiEWino_BR;
            std::cout << std:: endl;

          }**/

      }  // end of loop over 1st sQCD cascade.

  }



  sQCD_to_EWino_handler::sQCD_to_EWino_handler(
                                 input_handler const* const given_shortcuts ) :
     shortcut( given_shortcuts )
   {

     // just an initialization list.

   }

  sQCD_to_EWino_handler::~sQCD_to_EWino_handler()
  {

    for( std::vector< sQCD_to_EWino* >::iterator
         deletion_iterator = sQCD_to_EWinos.begin();
         sQCD_to_EWinos.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

  }


  sQCD_to_EWino*
  sQCD_to_EWino_handler::get_sQCD_to_EWino(
                signed_particle_shortcut_pair const* const given_scolored_pair,
                               std::list< int > const* const given_jets_to_tag,
                           std::list< int > const* const given_jets_to_reject )
  {

    sQCD_to_EWino* return_pointer = NULL;

    for( std::vector< sQCD_to_EWino* >::const_iterator
         checker_iterator = sQCD_to_EWinos.begin();
         sQCD_to_EWinos.end() > checker_iterator;
         ++checker_iterator )
      {

        if( (*checker_iterator)->is_requested( given_scolored_pair,
                                               given_jets_to_tag,
                                               given_jets_to_reject ) )
          // if we find a perfect match...
          {

            return_pointer = *checker_iterator;
            // return this exclusive_BR_calculator.
            checker_iterator = sQCD_to_EWinos.end();
            // stop looking.

          }

      }

    if( NULL == return_pointer )
      // if we didn't find a pre-existing match...
      {

        return_pointer = new sQCD_to_EWino( given_scolored_pair,
                                            shortcut,
                                            given_jets_to_tag,
                                            given_jets_to_reject );
        sQCD_to_EWinos.push_back( return_pointer );
        // add this calculator to the collection.

      }

    return return_pointer;

  }



  signal_shortcuts::signal_shortcuts(
                              input_handler const* const given_input_shortcuts,
                      crossSectionHandler* const given_cross_section_pointer,
                            kinematics_handler* const given_kinematics_pointer,
                               cascade_handler* const given_cascade_pointer ) :
    input_shortcuts( given_input_shortcuts ),
    sQCD_to_EWinos( given_input_shortcuts ),
    cross_section_pointer( given_cross_section_pointer ),
    kinematics_pointer( given_kinematics_pointer ),
    cascade_pointer( given_cascade_pointer )
  {

    // just an initialization list.

  }

  signal_shortcuts::~signal_shortcuts()
  {

    // does nothing.

  }


  double
  signal_shortcuts::get_uncertainty()
  /* this looks at the various mass differences in the hierarchy & then
   * guesses how uncertain the signal value is (as a multiplicative factor
   * > 1.0).
   */
  {

    // actually, right now it returns a blanket 1.3 factor...

    return 1.3;

  }


  productionChannelPointerSet::productionChannelPointerSet(
                                  signalDefinitionSet* const signalDefinitions,
                    signed_particle_shortcut_pair const* const scoloredPair ) :
    crossSection( signalDefinitions->getCrossSections()->getTable(
        scoloredPair ) ),
    firstCascadeSet( signalDefinitions->getShortcuts()->get_cascade_handler(
                                                            )->get_cascade_set(
                                         scoloredPair->get_first_pointer() ) ),
    secondCascadeSet( signalDefinitions->getShortcuts()->get_cascade_handler(
                                                            )->get_cascade_set(
                                         scoloredPair->get_second_pointer() ) )
  {
    // just an initialization list.
  }

  productionChannelPointerSet::~productionChannelPointerSet()
  {
    // does nothing.
  }



  signalDefinitionSet::signalDefinitionSet(
                                           signal_shortcuts* const shortcut ) :
  acceptanceCutSet(),
  shortcut( shortcut ),
  crossSections( NULL ),
  leptonKinematics( NULL ),
  jetPlusMetAcceptance( NULL )
  {
    // just an initialization list.
  }

  signalDefinitionSet::signalDefinitionSet(
                                     signalDefinitionSet* const copySource  ) :
  acceptanceCutSet( (acceptanceCutSet)copySource ),
  shortcut( copySource->getShortcuts() ),
  crossSections( copySource->getCrossSections() ),
  leptonKinematics( copySource->getLeptonKinematics() ),
  jetPlusMetAcceptance( copySource->getJetPlusMetAcceptance() )
  {
    // just an initialization list.
  }

  signalDefinitionSet::~signalDefinitionSet()
  {
    // does nothing.
  }

}  // end of LHC_FASER namespace
