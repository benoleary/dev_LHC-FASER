/*
 * LHC-FASER_signal_data_collection_stuff.hpp
 *
 *  Created on: 17 Dec 2010
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

#ifndef LHC_FASER_SIGNAL_DATA_COLLECTION_STUFF_HPP_
#define LHC_FASER_SIGNAL_DATA_COLLECTION_STUFF_HPP_

#include "LHC-FASER_input_handling_stuff.hpp"
#include "LHC-FASER_cross-section_stuff.hpp"
#include "LHC-FASER_kinematics_stuff.hpp"
#include "LHC-FASER_full_cascade_stuff.hpp"

namespace LHC_FASER
{

  /* this holds each production channel for a combination of colored sparticles
   * as a set of pointers to the signed_particle_shortcut_pair, the gluino
   * decay product, & the electroweakino end products of each colored sparticle
   * decay, along with the product of the branching ratios for the decay set.
   */
  class sQCD_to_EWino : public readied_for_new_point
  {

  public:

    sQCD_to_EWino(
                signed_particle_shortcut_pair const* const given_scolored_pair,
                   input_handler const* const given_shortcuts,
                   std::list< int > const* const given_jets_to_tag,
                   std::list< int > const* const given_jets_to_reject )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~sQCD_to_EWino()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    bool
    is_requested(
                signed_particle_shortcut_pair const* const given_scolored_pair,
                  std::list< int > const* const given_jets_to_tag,
                  std::list< int > const* const given_jets_to_reject )
    const
    /* this checks to see if the pointers are the same: the assumption is that
     * the lists of jets to tag or reject are kept unique, with pointers to
     * them being kept by all the objects which use them.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    std::vector< scoloreds_to_EWinos* > const*
    get_sQCD_end()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    signed_particle_shortcut_pair const*
    get_scolored_pair()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:

    enum basic_channel
    { gg,
      sg,
      ss }
    channel_type;

    input_handler const* const shortcut;
    signed_particle_shortcut_pair const* const scolored_pair;

    std::vector< scoloreds_to_EWinos* > sQCD_end;
    int number_of_sQCD_steps_in_first_channel;
    int number_of_sQCD_steps_in_second_channel;
    bool channel_is_valid_flag;
    std::list< int > const* const jets_to_tag;
    std::list< int > const* const jets_to_reject;

    colored_cascade_set* first_scolored_cascade;
    colored_cascade_set* second_scolored_cascade;
    std::vector< exclusive_BR_calculator* > gluino_to_scolored_BRs;
    std::vector< std::vector< exclusive_BR_calculator* >* >
    scolored_to_EWino_BRs;
    // scolored_to_EWino_BRs.at( i )->at( j ) is the BR for scolored i to decay
    // to EWino j.
    std::vector< std::vector< exclusive_BR_calculator* >* >
    scolored_to_antiEWino_BRs;
    // scolored_to_antiEWino_BRs.at( i )->at( j ) is the BR for scolored i to
    // decay to antiEWino j.
    exclusive_BR_calculator* first_scolored_to_gluino_BR;
    exclusive_BR_calculator* second_scolored_to_gluino_BR;
    std::vector< exclusive_BR_calculator* >* gluino_to_EWino_BRs;
    std::vector< exclusive_BR_calculator* >* first_scolored_to_EWino_BRs;
    std::vector< exclusive_BR_calculator* >* first_scolored_to_antiEWino_BRs;
    std::vector< exclusive_BR_calculator* >* second_scolored_to_EWino_BRs;
    std::vector< exclusive_BR_calculator* >* second_scolored_to_antiEWino_BRs;


    bool
    cascade_is_not_negligible( colored_cascade const* const cascade,
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
     * provided branching ratio calculators. it returns true if at least 1 of
     * the final possibilities (whether the squark before the electroweakino is
     * a particle or antiparticle, & whether the electroweakino is a particle
     * or antiparticle) is not negligible, as well as assigning the relevant
     * values to the given doubles.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    add_EWinos_to_sQCD_end( colored_cascade const* const first_cascade,
                            bool const first_end_scolored_is_not_antiparticle,
                            bool const first_EWino_is_not_antiparticle,
                            colored_cascade const* const second_cascade,
                            double const second_end_scolored_to_EWino_BR,
                            double const second_end_scolored_to_antiEWino_BR,
                            double const second_end_antiscolored_to_EWino_BR,
                          double const second_end_antiscolored_to_antiEWino_BR,
                            double const first_sQCD_BR )
    // this adds the appropriate combinations of squarks &/or antisquarks with
    // electroweakinos &/or antielectroweakinos to sQCD_end.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    prepare_channels()
    /* this checks to see if any of the kinematically-allowed cascade
     * combinations following from the channel survive being multiplied by
     * their branching ratios. those that do survive get put into sQCD_end for
     * later processing by the cascade handler.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };  // end of sQCD_to_EWino class.

  // this class holds a set of direct_decay_checkers & passes out pointers
  // to them.
  class sQCD_to_EWino_handler
  {

  public:

    sQCD_to_EWino_handler( input_handler const* const given_shortcuts )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~sQCD_to_EWino_handler()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    sQCD_to_EWino*
    get_sQCD_to_EWino(
                signed_particle_shortcut_pair const* const given_scolored_pair,
                       std::list< int > const* const given_jets_to_tag,
                       std::list< int > const* const given_jets_to_reject )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    std::vector< sQCD_to_EWino* > sQCD_to_EWinos;
    input_handler const* const shortcut;

  };  // end of decay_checker_handler class.


  /* this class holds a pointer to an input_handler instance along with
   * various pointers to cross-section, kinematics, & cascade handlers, for
   * ease of reference. it also keeps an instance of the sQCD_to_EWino_handler
   * class.
   */
  class signal_shortcuts
  {

  public:

    signal_shortcuts( input_handler const* const given_input_shortcuts,
                      cross_section_handler* const given_cross_section_pointer,
                      kinematics_handler* const given_kinematics_pointer,
                      cascade_handler* const given_cascade_pointer )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~signal_shortcuts()
    /* code after the classes in this .hpp file, or in the .cpp file. */;



    input_handler const*
    get_input_shortcuts()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    cross_section_handler*
    get_cross_section_handler()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    kinematics_handler*
    get_kinematics_handler()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    cascade_handler*
    get_cascade_handler()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    sQCD_to_EWino*
    get_sQCD_to_EWinos(
                signed_particle_shortcut_pair const* const given_scolored_pair,
                        std::list< int > const* const given_jets_to_tag,
                        std::list< int > const* const given_jets_to_reject )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    double
    get_uncertainty()
    /* this looks at the various mass differences in the hierarchy & then
     * guesses how uncertain the signal value is (as a multiplicative factor
     * > 1.0).
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    input_handler const* const input_shortcuts;
    sQCD_to_EWino_handler sQCD_to_EWinos;
    cross_section_handler* const cross_section_pointer;
    kinematics_handler* const kinematics_pointer;
    cascade_handler* const cascade_pointer;

  };  // end of signal_shortcuts class.


  /* this holds together the various doubles & pointers that store information
   * that is particular to a signal. it doesn't hold all the information
   * required by every single signal type, & it does hold a little bit more
   * than is necessary for the least elaborate signal, but it should be a
   * convenient package of information. it builds on the acceptanceCutSet
   * class.
   */
  class signalDefinitionSet : public acceptanceCutSet
  {

  public:
    signalDefinitionSet( signal_shortcuts* const shortcut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    signalDefinitionSet( signalDefinitionSet* const copySource )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~signalDefinitionSet()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    signal_shortcuts*
    getShortcuts()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setBeamEnergy( double const inputValue )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    cross_section_table_set*
    getCrossSections()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    lepton_acceptance_table*
    getLeptonKinematics()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    jet_acceptance_table*
    getJetPlusMetAcceptance()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setJetPlusMetAcceptance( jet_acceptance_table* const inputTable )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    signal_shortcuts* const shortcut;

    // each signal needs to look up specific tables, based on the beam energy:
    cross_section_table_set* crossSections;
    lepton_acceptance_table* leptonKinematics;
    jet_acceptance_table* jetPlusMetAcceptance;
    // by default each signal has only a single jet+MET table to look up; if a
    // signal needs to look up more tables, it should store the extras itself.
  };  // end of signal_definitions class.


  /* this holds each production channel for a combination of colored sparticles
   * as a set of pointers to the cross_section_table (which holds a pointer to
   * the signed_particle_shortcut_pair) & the relevant cascade_sets for the
   * sparticles.
   */
  class productionChannelPointerSet
  {
  public:
    productionChannelPointerSet( signalDefinitionSet* const signalDefinitions,
                      signed_particle_shortcut_pair const* const scoloredPair )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~productionChannelPointerSet()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    signed_particle_shortcut_pair const*
    getScoloredPair()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getCrossSection()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    fullCascadeSet*
    getFirstCascadeSet()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    fullCascadeSet*
    getSecondCascadeSet()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    //signal_shortcuts* const shortcut;
    cross_section_table* crossSection;
    fullCascadeSet* firstCascadeSet;
    fullCascadeSet* secondCascadeSet;
  };  // end of productionChannelPointerSet class.





  // inline functions:


  inline bool
  sQCD_to_EWino::is_requested(
                signed_particle_shortcut_pair const* const given_scolored_pair,
                               std::list< int > const* const given_jets_to_tag,
                           std::list< int > const* const given_jets_to_reject )
  const
  /* this checks to see if the pointers are the same: the assumption is that
   * the lists of jets to tag or reject are kept unique, with pointers to them
   * being kept by all the objects which use them.
   */
  {

    if( ( scolored_pair == given_scolored_pair )
        &&
        ( jets_to_tag == given_jets_to_tag )
        &&
        ( jets_to_reject == given_jets_to_reject ) )
      {

        return true;

      }
    else
      {

        return false;

      }

  }

  inline std::vector< scoloreds_to_EWinos* > const*
  sQCD_to_EWino::get_sQCD_end()
  {

    if( needs_to_prepare_for_this_point() )
      {

        prepare_channels();
        finish_preparing_for_this_point();

      }

    return &sQCD_end;

  }

  inline signed_particle_shortcut_pair const*
  sQCD_to_EWino::get_scolored_pair()
  const
  {

    return scolored_pair;

  }



  inline input_handler const*
  signal_shortcuts::get_input_shortcuts()
  const
  {

    return input_shortcuts;

  }


  inline cross_section_handler*
  signal_shortcuts::get_cross_section_handler()
  const
  {

    return cross_section_pointer;

  }

  inline kinematics_handler*
  signal_shortcuts::get_kinematics_handler()
  const
  {

    return kinematics_pointer;

  }

  inline cascade_handler*
  signal_shortcuts::get_cascade_handler()
  const
  {

    return cascade_pointer;

  }


  inline sQCD_to_EWino*
  signal_shortcuts::get_sQCD_to_EWinos(
              signed_particle_shortcut_pair const* const given_scolored_pair,
                      std::list< int > const* const given_jets_to_tag,
                      std::list< int > const* const given_jets_to_reject )
  {

    return sQCD_to_EWinos.get_sQCD_to_EWino( given_scolored_pair,
                                             given_jets_to_tag,
                                             given_jets_to_reject );

  }



  inline signal_shortcuts*
  signalDefinitionSet::getShortcuts()
  {
    return shortcut;
  }

  inline void
  signalDefinitionSet::setBeamEnergy( double const inputValue )
  {
    beamEnergy = inputValue;
    crossSections
    = shortcut->get_cross_section_handler()->get_table_set( inputValue );
    leptonKinematics
    = shortcut->get_kinematics_handler()->get_lepton_acceptance_table(
                                                                  inputValue );
  }

  inline cross_section_table_set*
  signalDefinitionSet::getCrossSections()
  {
    return crossSections;
  }

  inline lepton_acceptance_table*
  signalDefinitionSet::getLeptonKinematics()
  {
    return leptonKinematics;
  }

  inline jet_acceptance_table*
  signalDefinitionSet::getJetPlusMetAcceptance()
  // this returns the 1st element of the vector. it is intended for when the
  // vector has only 1 entry.
  {
    return jetPlusMetAcceptance;
  }

  inline void
  signalDefinitionSet::setJetPlusMetAcceptance(
                                      jet_acceptance_table* const inputTable )
  // this sets the 1st element of the vector. it is intended for when the
  // vector has only 1 entry.
  {
    jetPlusMetAcceptance = inputTable;
  }



  inline signed_particle_shortcut_pair const*
  productionChannelPointerSet::getScoloredPair()
  const
  {

    return crossSection->get_pair();

  }

  inline double
  productionChannelPointerSet::getCrossSection()
  {
    return crossSection->get_value();
  }

  inline fullCascadeSet*
  productionChannelPointerSet::getFirstCascadeSet()
  {
    return firstCascadeSet;
  }

  inline fullCascadeSet*
  productionChannelPointerSet::getSecondCascadeSet()
  {
    return secondCascadeSet;
  }
}  // end of LHC_FASER namespace

#endif /* LHC_FASER_SIGNAL_DATA_COLLECTION_STUFF_HPP_ */
