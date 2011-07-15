/*
 * LHC-FASER_sparticle_decay_stuff.cpp
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

#include "LHC-FASER_sparticle_decay_stuff.hpp"

namespace LHC_FASER
{

  signed_particle_shortcut_pair::signed_particle_shortcut_pair(
             CppSLHA::particle_property_set const* const given_first_pointer,
                               bool const given_first_is_not_antiparticle_flag,
            CppSLHA::particle_property_set const* const given_second_pointer,
                             bool const given_second_is_not_antiparticle_flag ) :
    first_particle_pointer( given_first_pointer ),
    first_is_not_antiparticle_flag( given_first_is_not_antiparticle_flag ),
    second_particle_pointer( given_second_pointer ),
    second_is_not_antiparticle_flag( given_second_is_not_antiparticle_flag )
  {

    // just an initialization list.

  }

  signed_particle_shortcut_pair::signed_particle_shortcut_pair(
                     signed_particle_shortcut_pair const* const copy_source ) :
    first_particle_pointer( copy_source->get_first_pointer() ),
    first_is_not_antiparticle_flag( copy_source->first_is_not_antiparticle() ),
    second_particle_pointer( copy_source->get_second_pointer() ),
    second_is_not_antiparticle_flag(
                                    copy_source->second_is_not_antiparticle() )
  {

    // just an initialization list.

  }

  signed_particle_shortcut_pair::~signed_particle_shortcut_pair()
  {

    // does nothing.

  }



  direct_decay_checker::direct_decay_checker(
                     CppSLHA::particle_property_set const* const given_decayer,
              std::vector< CppSLHA::particle_property_set const* > const* const
                                              given_products,
                                 readier_for_new_point* const given_readier ) :
    readied_for_new_point( given_readier ),
    decayer( given_decayer ),
    products( given_products ),
    decayer_does_decay_to_at_least_one_product( false )
  {

    // just an initialization list.

  }

  direct_decay_checker::~direct_decay_checker()
  {

    // does nothing.

  }


  bool
  direct_decay_checker::check_for_decay()
  // this checks the particles to see if the decay exists for this point.
  {

    bool return_bool = false;

    // check each possible decay mode:
    for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
         product_iterator = products->begin();
         products->end() > product_iterator;
         ++product_iterator )
      {

        if( decayer->get_absolute_mass()
            > (*product_iterator)->get_absolute_mass() )
          // if the decayer can decay into this product...
          {

            // note that the decayer decays to at least 1 of the products in
            // products with a 2-body decay:
            return_bool = true;
            product_iterator = products->end();
            // stop searching.

          }

      }

    return return_bool;

  }



  decay_checker_handler::decay_checker_handler(
                                 readier_for_new_point* const given_readier ) :
    readier( given_readier )
  {

    // just an initialization list.

  }

  decay_checker_handler::~decay_checker_handler()
  {

    for( std::vector< direct_decay_checker* >::iterator
         deletion_iterator = decay_checkers.begin();
         decay_checkers.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

  }


  direct_decay_checker*
  decay_checker_handler::get_decay_checker(
                     CppSLHA::particle_property_set const* const given_decayer,
              std::vector< CppSLHA::particle_property_set const* > const* const
                                            given_products )
  {

    direct_decay_checker* return_pointer = NULL;

    for( std::vector< direct_decay_checker* >::const_iterator
         checker_iterator = decay_checkers.begin();
         decay_checkers.end() > checker_iterator;
         ++checker_iterator )
      {

        if( ( given_decayer == (*checker_iterator)->get_decayer() )
            &&
            ( given_products == (*checker_iterator)->get_products() ) )
          // if we find a perfect match...
          {

            return_pointer = *checker_iterator;
            // return this exclusive_BR_calculator.
            checker_iterator = decay_checkers.end();
            // stop looking.

          }

      }

    if( NULL == return_pointer )
      // if we didn't find a pre-existing match...
      {

        return_pointer = add( given_decayer,
                              given_products );
        // add this calculator to the collection.

      }

    return return_pointer;

  }



  exclusive_BR_calculator::exclusive_BR_calculator(
                     CppSLHA::particle_property_set const* const given_decayer,
                     CppSLHA::particle_property_set const* const given_product,
                                   bool const product_is_not_antiparticle_flag,
                                  std::list< int > const* const given_excluded,
                                 readier_for_new_point* const given_readier ) :
    readied_for_new_point( given_readier ),
    decayer( given_decayer ),
    product( given_product ),
    exclusion_list( given_excluded ),
    BR_value( 0.0 )
  {

    if( product_is_not_antiparticle_flag )
      {

        product_code = new std::list< int >( 1,
                                             given_product->get_PDG_code() );

      }
    else
      {

        product_code = new std::list< int >( 1,
                                            -(given_product->get_PDG_code()) );

      }

    // debugging:
    /**std::cout
    << std::endl << "debugging: product code is a list of:" << std::endl;
    for( std::list< int >::const_iterator
         list_iterator = product_code->begin();
         product_code->end() != list_iterator;
         ++list_iterator )
      {
        std::cout << " " << *list_iterator;
      }
    std::cout << std::endl;**/

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: exclusive_BR_calculator::exclusive_BR_calculator("
    << *(given_decayer->get_name()) << ", "
    << *(given_product->get_name_or_antiname(
                                            product_is_not_antiparticle_flag ))
    << ", excluded:{";
    for( std::list< int >::const_iterator
         list_iterator = given_excluded->begin();
         given_excluded->end() != list_iterator;
         ++list_iterator )
      {
        std::cout << " " << *list_iterator;
      }
    std::cout << " } ) called.";**/

  }

  exclusive_BR_calculator::~exclusive_BR_calculator()
  {

    delete product_code;

  }



  exclusive_BR_handler::exclusive_BR_handler(
      std::vector< CppSLHA::particle_property_set const* > const* const sdowns,
        std::vector< CppSLHA::particle_property_set const* > const* const sups,
                                 readier_for_new_point* const given_readier ) :
    readier( given_readier )
  {

    for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
         sdown_iterator = sdowns->begin();
         sdowns->end() > sdown_iterator;
         ++sdown_iterator )
      {

        always_neglected_decays.push_back( new std::pair< int,
                                                          int >(
                                             (*sdown_iterator)->get_PDG_code(),
                                           CppSLHA::PDG_code::chargino_one ) );
        always_neglected_decays.push_back( new std::pair< int,
                                                          int >(
                                             (*sdown_iterator)->get_PDG_code(),
                                           CppSLHA::PDG_code::chargino_two ) );

      }

    for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
         sup_iterator = sups->begin();
         sups->end() > sup_iterator;
         ++sup_iterator )
      {

        always_neglected_decays.push_back( new std::pair< int,
                                                          int >(
                                               (*sup_iterator)->get_PDG_code(),
                                          -CppSLHA::PDG_code::chargino_one ) );
        always_neglected_decays.push_back( new std::pair< int,
                                                          int >(
                                               (*sup_iterator)->get_PDG_code(),
                                          -CppSLHA::PDG_code::chargino_two ) );

      }

  }

  exclusive_BR_handler::~exclusive_BR_handler()
  {

    for( std::vector< exclusive_BR_calculator* >::iterator
         deletion_iterator = exclusive_BRs.begin();
         exclusive_BRs.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

    for( std::vector< std::pair< int,
                                 int >* >::iterator
         deletion_iterator = always_neglected_decays.begin();
         always_neglected_decays.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

  }


  exclusive_BR_calculator*
  exclusive_BR_handler::get_BR_calculator(
                     CppSLHA::particle_property_set const* const given_decayer,
                     CppSLHA::particle_property_set const* const given_product,
                                   bool const product_is_not_antiparticle_flag,
                                 std::list< int > const* const given_excluded )
  {

    exclusive_BR_calculator* return_pointer = NULL;

    int given_product_code;
    if( product_is_not_antiparticle_flag )
      {

        given_product_code = given_product->get_PDG_code();

      }
    else
      {

        given_product_code = -(given_product->get_PDG_code());

      }

    bool always_neglected = false;

    for( std::vector< std::pair< int,
                                 int >* >::iterator
         neglection_iterator = always_neglected_decays.begin();
         always_neglected_decays.end() > neglection_iterator;
         ++neglection_iterator )
      // check whether this decay is always neglected because it cannot happen
      // often enough no matter what the scenario:
      {

        if( given_decayer->get_PDG_code() == (*neglection_iterator)->first )
          {

            if( given_product_code == (*neglection_iterator)->second )
              {

                always_neglected = true;
                neglection_iterator = always_neglected_decays.end();
                // stop looking.

              }

          }

      }

    if( !always_neglected )
      {

        for( std::vector< exclusive_BR_calculator* >::iterator
             BR_iterator = exclusive_BRs.begin();
             exclusive_BRs.end() > BR_iterator;
             ++BR_iterator )
          {

            if( ( given_decayer == (*BR_iterator)->get_decayer() )
                &&
                ( given_product == (*BR_iterator)->get_product() )
                &&
                ( given_excluded == (*BR_iterator)->get_excluded() ) )
              // if we find a perfect match...
              {

                return_pointer = *BR_iterator;
                // return this exclusive_BR_calculator.
                BR_iterator = exclusive_BRs.end();
                // stop looking.

              }

          }

        if( NULL == return_pointer )
          // if we didn't find a pre-existing match...
          {

            return_pointer = add( given_decayer,
                                  given_product,
                                  product_is_not_antiparticle_flag,
                                  given_excluded );
            // add this calculator to the collection.

          }

      }
    // end of if( !always_neglected ) - if always_neglected == true, then
    // NULL is returned. we should always check whether we got a NULL!

    return return_pointer;

  }



  colored_cascade::colored_cascade(
                     CppSLHA::particle_property_set const* const given_decayer,
                     CppSLHA::particle_property_set const* const given_product,
                                    int const given_product_position,
                       CppSLHA::particle_property_set const* const given_EWino,
                                    int const given_EWino_position ) :
    decayer( given_decayer ),
    product( given_product ),
    product_position( given_product_position ),
    EWino( given_EWino ),
    EWino_position( given_EWino_position )
  {

    if( ( CppSLHA::PDG_code::gluino != given_decayer->get_PDG_code() )
        &&
        ( given_decayer->get_PDG_code() == given_product->get_PDG_code() ) )
      {

        cascade_code = sx;

      }
    else if( ( CppSLHA::PDG_code::gluino == given_decayer->get_PDG_code() )
             &&
             ( CppSLHA::PDG_code::gluino == given_product->get_PDG_code() ) )
      {

        cascade_code = gx;

      }
    else if( ( CppSLHA::PDG_code::gluino != given_decayer->get_PDG_code() )
             &&
             ( CppSLHA::PDG_code::gluino == given_product->get_PDG_code() ) )
      {

        cascade_code = sgx;

      }
    else if( ( CppSLHA::PDG_code::gluino == given_decayer->get_PDG_code() )
             &&
             ( CppSLHA::PDG_code::gluino != given_product->get_PDG_code() ) )
      {

        cascade_code = gsx;

      }
    else
      {

        cascade_code = sgsx;

      }

  }

  colored_cascade::~colored_cascade()
  {

    // does nothing.

  }



  colored_cascade_set::colored_cascade_set(
                     CppSLHA::particle_property_set const* const given_decayer,
                      CppSLHA::particle_property_set const* const given_gluino,
              std::vector< CppSLHA::particle_property_set const* > const* const
                                            given_scoloreds,
              std::vector< CppSLHA::particle_property_set const* > const* const
                                            given_EWinos,
                                 readier_for_new_point* const given_readier ) :
    readied_for_new_point( given_readier ),
    decayer( given_decayer ),
    gluino_pointer( given_gluino ),
    scoloreds( given_scoloreds ),
    EWinos( given_EWinos )
  {

    // just an initialization list.

  }

  colored_cascade_set::~colored_cascade_set()
  {

    for( std::vector< colored_cascade* >::iterator
         deletion_iterator = colored_cascades.begin();
         colored_cascades.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

  }


  void
  colored_cascade_set::prepare_open_channels()
  /* this fills colored_cascades with "good" combinations, ignoring "bad"
   * combinations, of possible decays of a colored sparticle to any
   * electroweakino.
   * in the following, h is a squark heavier than the gluino (or equal mass),
   * g is the gluino, & l is a squark lighter than the gluino (definitely
   * lower mass).
   * the "bad" combinations are:
   * - g/l, h (l & g cannot decay into h, duh)
   * - l, g (if there is an l, g decays into it)
   * - h/l, h/l, if the decay mode squark is *not* the same as the the decay
   *   product squark, to avoid overcounting.
   * the "good" combinations are:
   * - h/l, h/l, but only if the decay mode squark is the same as the decay
   *   product squark, to avoid overcounting.
   * - h, g (then we assume that h -> g, & that g decays via a 3-body decay).
   * - g, g (then we assume that g decays via a 3-body decay).
   * - g, l (then we assume that g -> l).
   * - h, l (then we assume that h -> g & then g -> l).
   */
  {

    // 1st, empty colored_cascades:
    for( std::vector< colored_cascade* >::iterator
         deletion_iterator = colored_cascades.begin();
         colored_cascades.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }
    colored_cascades.clear();


    // next, check through all the potential colored decay products:
    for( int scolored_counter = 0;
         scoloreds->size() > scolored_counter;
         ++scolored_counter )
      {

        if( ( scoloreds->at( scolored_counter )->get_PDG_code()
              == decayer->get_PDG_code() )
            // if it should be a direct decay to an electroweakino...
            ||
            ( ( ( CppSLHA::PDG_code::gluino == decayer->get_PDG_code() )
                ||
                ( CppSLHA::PDG_code::gluino
                  == scoloreds->at( scolored_counter )->get_PDG_code() ) )
              &&
              ( scoloreds->at( scolored_counter )->get_absolute_mass()
                < decayer->get_absolute_mass() ) )
            // or if it's a gluino decay to a squark or a squark decay to a
            // gluino, then directly to an electroweakino...
            ||
            ( ( scoloreds->at( scolored_counter )->get_absolute_mass()
                < gluino_pointer->get_absolute_mass() )
              &&
              ( gluino_pointer->get_absolute_mass()
                < decayer->get_absolute_mass() ) ) )
            // or if it's the decay of a squark to a gluino then to a lighter
            // squark which then decays to an electroweakino
          {

            // go through all the possible electroweakinos which can come from
            // scoloreds->at( scolored_counter ):
            for( int EWino_counter = 0;
                 EWinos->size() > EWino_counter;
                 ++EWino_counter )
              {

                if( EWinos->at( EWino_counter )->get_absolute_mass()
                    < scoloreds->at( scolored_counter )->get_absolute_mass() )
                  // if the decay to this electroweakino is kinematically
                  // allowed...
                  {

                    colored_cascades.push_back( new colored_cascade( decayer,
                                             scoloreds->at( scolored_counter ),
                                                              scolored_counter,
                                                   EWinos->at( EWino_counter ),
                                                             EWino_counter ) );

                  }

              }  // end of loop over electroweakinos.

          }  // end of if the scolored particle mass hierarchy is correct.

      }  // end of loop over colored sparticles.

  }



  colored_cascade_handler::colored_cascade_handler(
                    CppSLHA::particle_property_set const* const given_gluino,
            std::vector< CppSLHA::particle_property_set const* > const* const
                                                    given_scoloreds,
            std::vector< CppSLHA::particle_property_set const* > const* const
                                                    given_EWinos,
                                 readier_for_new_point* const given_readier ) :
    gluino_pointer( given_gluino ),
    scoloreds( given_scoloreds ),
    EWinos( given_EWinos ),
    readier( given_readier )
  {

    // just an initialization list;

  }

  colored_cascade_handler::~colored_cascade_handler()
  {

    for( std::vector< colored_cascade_set* >::iterator
         deletion_iterator = colored_cascade_sets.begin();
         colored_cascade_sets.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

  }


  colored_cascade_set*
  colored_cascade_handler::get_colored_cascade_set(
                    CppSLHA::particle_property_set const* const given_decayer )
  {

    colored_cascade_set* return_pointer = NULL;

    // look to see if we already have a colored_cascade_set for this decayer:
    for( std::vector< colored_cascade_set* >::iterator
         search_iterator = colored_cascade_sets.begin();
         colored_cascade_sets.end() > search_iterator;
         ++search_iterator )
      {

        if( (*search_iterator)->get_decayer() == given_decayer )
          {

            return_pointer = *search_iterator;
            // note the pointer to return.
            search_iterator = colored_cascade_sets.end();
            // stop looking.

          }

      }

    if( NULL == return_pointer )
      // if we didn't have a colored_cascade_set for the decayer, make 1:
      {

        return_pointer = new colored_cascade_set( given_decayer,
                                                  gluino_pointer,
                                                  scoloreds,
                                                  EWinos,
                                                  readier );

        colored_cascade_sets.push_back( return_pointer );

      }

    return return_pointer;

  }



  scoloreds_to_EWinos::scoloreds_to_EWinos(
                signed_particle_shortcut_pair const* const given_scolored_pair,
                              colored_cascade const* const given_first_cascade,
                       bool const given_first_end_scolored_is_not_antiparticle,
                              bool const given_first_EWino_is_not_antiparticle,
                             colored_cascade const* const given_second_cascade,
                      bool const given_second_end_scolored_is_not_antiparticle,
                             bool const given_second_EWino_is_not_antiparticle,
                                            double given_BR_product ) :
    scolored_pair( given_scolored_pair ),
    first_cascade( given_first_cascade ),
    first_end_scolored_is_not_antiparticle_flag(
                                given_first_end_scolored_is_not_antiparticle ),
    first_EWino_is_not_antiparticle_flag(
                                       given_first_EWino_is_not_antiparticle ),
    second_cascade( given_second_cascade ),
    second_end_scolored_is_not_antiparticle_flag(
                               given_second_end_scolored_is_not_antiparticle ),
    second_EWino_is_not_antiparticle_flag(
                                      given_second_EWino_is_not_antiparticle ),
    BR_product( given_BR_product )
  {

    // just an initialization list.

  }

  scoloreds_to_EWinos::~scoloreds_to_EWinos()
  {

    // does nothing.

  }


}  // end of LHC_FASER namespace
