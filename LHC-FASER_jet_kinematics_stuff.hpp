/*
 * LHC-FASER_kinematics_stuff.hpp
 *
 *  Created on: 15 Dec 2010
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

#ifndef LHC_FASER_JET_KINEMATICS_STUFF_HPP_
#define LHC_FASER_JET_KINEMATICS_STUFF_HPP_

#include "LHC-FASER_input_handling_stuff.hpp"
#include "LHC-FASER_full_cascade_stuff.hpp"


/* this file contains all the code relevant to obtaining the acceptance for
 * jet transverse momentum plus missing transverse momentum cuts.
 *
 * if the format of the grid files changes, a lot of the code here has to
 * change too. the following are format-dependent:
 * - all of the leptonAcceptanceParameterSet class.
 * - all of the acceptanceGrid class.
 * - all of the leptonAcceptanceGrid class, in particular
 *   -- leptonAcceptanceGrid::getValue()
 */

namespace LHC_FASER
{
  // this class holds an acceptanceGrid with a string identifying the type of
  // jet+MET signal which the acceptances are for.
  class jetAcceptanceTable
  {
  public:
    enum usedTypes
    {
      gx = 0,
      sx = 1,
      gsx = 2,
      sgx = 3,
      sizeOfEnum = 4 /* should be useful */
    };
    // this could be expanded, provided grids that account for multiple squark
    // flavors are provided.

    jetAcceptanceTable( std::string const* const gridFilesLocation,
                        std::string const* const jetCutName,
                        int const acceptanceColumn,
                        input_handler const* const shortcut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~jetAcceptanceTable()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    std::string const*
    getCutName()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getAcceptance( fullCascade const* const firstCascade,
                   fullCascade const* const secondCascade )
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    std::string const jetCutName;
    int const acceptanceColumn;
    input_handler const* const shortcut;
    acceptanceGrid* gluinoGluinoGrid;
    acceptanceGrid* squarkGluinoGrid;
    acceptanceGrid* squarkAntisquarkGrid;
    acceptanceGrid* squarkSquarkGrid;
    acceptanceGrid* gridToUse;
    // this is decided by the types of cascade given.
    std::vector< std::vector< acceptanceGrid* >* > gridsMatrix;
    double squarkMassForGrid;
    double gluinoMassForGrid;
    double lighterNeutralinoMassForGrid;
    double heavierNeutralinoMassForGrid;
    usedTypes firstCascadeType;
    double firstCascadeSquarkMass;
    usedTypes secondCascadeType;
    double secondCascadeSquarkMass;

    int
    getIntForCascadeType( fullCascade const* const givenCascade,
                          double* const squarkMassFromCascade )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };



  /* this class holds pointers to a jet_acceptance_table & a
   * lepton_acceptance_table for looking up the jets plus missing transverse
   * momentum acceptance with an int for which column to use & the lepton
   * acceptances, & gives out pointers from the jet_acceptance_table & the
   * lepton_acceptance_table to the value-calculating instances they produce
   * (which are all released & re-created with every update for a new point).
   */
  class kinematics_table
  {

  public:

    kinematics_table( jet_acceptance_table* const given_jets_table,
                      int const given_jet_acceptance_column,
                      lepton_acceptance_table* const given_leptons_table )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~kinematics_table()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    int
    get_column()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    jet_acceptance_value*
    get_jet_acceptance( signed_particle_shortcut_pair const* const given_pair,
                        colored_cascade const* const given_first_sQCD_cascade,
                       colored_cascade const* const given_second_sQCD_cascade )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    lepton_acceptance_value*
    get_lepton_acceptance(
                         signed_particle_shortcut_pair const* const given_pair,
                         colored_cascade const* const given_first_sQCD_cascade,
                       colored_cascade const* const given_second_sQCD_cascade )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    jet_acceptance_table* const jets_table;
    int const jet_acceptance_column;
    lepton_acceptance_table* const leptons_table;

  };

  /* this class holds several kinematics_table instances as lookup tables for
   * various colored sparticle combinations, & returns pointers to them for
   * given kinematics_identifer pointers. it holds the lepton_acceptance_table
   * instances directly, giving pointers to them to its kinematics_table
   * instances.
   */
  class kinematics_table_set
  {

  public:

    kinematics_table_set( int const given_LHC_energy_in_TeV,
                          std::string const* const given_grid_directory,
                          std::string const* const given_jet_subdirectory,
                          lepton_acceptance_table* const given_lepton_table,
                          input_handler const* const given_shortcuts )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~kinematics_table_set()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    kinematics_table*
    get_table( int const given_acceptance_column )
    // this returns the kinematics_table for the requested column.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    bool
    is_requested( int const given_LHC_energy_in_TeV,
                  std::string const* const given_jet_subdirectory )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:

    input_handler const* const shortcut;

    int const LHC_energy_in_TeV;
    std::string const grid_directory;
    std::string jet_grid_subdirectory;

    lepton_acceptance_table* lepton_table;
    jet_acceptance_table* jet_table;
    std::vector< kinematics_table* > kinematics_tables;

  };


  // this holds a set of kinematics_table_set instances & returns pointers to
  // them based on the given LHC energy & signal name.
  class kinematics_handler
  {

  public:

    kinematics_handler( input_handler const* const given_shortcuts )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~kinematics_handler()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    kinematics_table*
    get_table( int const LHC_energy_in_TeV,
               std::string const* const given_jet_subdirectory,
               int const given_acceptance_column )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    lepton_acceptance_table*
    get_lepton_acceptance_table( int const LHC_energy_in_TeV )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    enum sQCD_cascade_type
    {

      sx /* direct squark to electroweakino */,
      sw /* direct squark to other squark plus W boson */,
      gx /* direct gluino to electroweakino */,
      sgx /* squark to gluino to electroweakino */,
      gsx /* gluino to squark to electroweakino */,
      gsw /* gluino to squark to other squark plus W boson */,
      sgsx /* squark to gluino to squark to electroweakino */,
      sgsw /* squark to gluino to squark to other squark plus W boson */

      /* we don't take into account squark -> other squark + W boson, other
       * squark -> gluino etc., because if the 1st decay is open, so is the
       * decay of the heavier squark directly into the gluino, & for the
       * purposes of this code, we assume that that sQCD decay trumps the EW
       * decay.
       */

    };


  protected:

    input_handler const* const shortcut;

    std::vector< kinematics_table_set* > table_sets;
    std::vector< lepton_acceptance_table* > lepton_tables;

    kinematics_table_set*
    get_table_set( int const LHC_energy_in_TeV,
                   std::string const* const given_jet_subdirectory )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };





  // inline functions:

  inline std::string const*
  jetAcceptanceTable::getCutName()
  const
  {
    return &jetCutName;
  }

  inline void
  acceptanceCutSet::becomeCopyOf( acceptanceCutSet* const copyPointer )
  // this copies the values from a given acceptanceCutSet.
  {
    beamEnergy = copyPointer->getBeamEnergy();
    primaryLeptonCut = copyPointer->getPrimaryLeptonCut();
    secondaryLeptonCut = copyPointer->getSecondaryLeptonCut();
    jetCut = copyPointer->getJetCut();
  }

  inline bool
  acceptanceCutSet::compareJetAndBothLeptonCuts(
                                          acceptanceCutSet* const firstPointer,
                                        acceptanceCutSet* const secondPointer )
  // this returns true if both acceptanceCutSets have the same jet cut & both
  // lepton cuts.
  {
    if( firstPointer->isSameAcceptanceCutSet( secondPointer ) )
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  inline bool
  acceptanceCutSet::compareJustBothLeptonCuts(
                                          acceptanceCutSet* const firstPointer,
                                        acceptanceCutSet* const secondPointer )
  {
    if( ( firstPointer->getPrimaryLeptonCut()
          == secondPointer->getPrimaryLeptonCut() )
        &&
        ( firstPointer->getSecondaryLeptonCut()
          == secondPointer->getSecondaryLeptonCut() ) )
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  inline bool
  acceptanceCutSet::compareJustPrimaryLeptonCut(
                                          acceptanceCutSet* const firstPointer,
                                        acceptanceCutSet* const secondPointer )
  {
    if( firstPointer->getPrimaryLeptonCut()
        == secondPointer->getPrimaryLeptonCut() )
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  inline bool
  acceptanceCutSet::justReturnTrue( acceptanceCutSet* const firstPointer,
                                    acceptanceCutSet* const secondPointer )
  {
    return true;
  }


  inline double
  acceptanceCutSet::getBeamEnergy()
  const
  {

    return beamEnergy;

  }


  inline double
  acceptanceCutSet::getPrimaryLeptonCut()
  const
  {

    return primaryLeptonCut;

  }

  inline void
  acceptanceCutSet::setPrimaryLeptonCut( double inputValue )
  {

    primaryLeptonCut = inputValue;

  }


  inline double
  acceptanceCutSet::getSecondaryLeptonCut()
  const
  {

    return secondaryLeptonCut;

  }

  inline void
  acceptanceCutSet::setSecondaryLeptonCut( double inputValue )
  {

    secondaryLeptonCut = inputValue;

  }


  inline double
  acceptanceCutSet::getLeptonCut()
  const
  {

    return getPrimaryLeptonCut();

  }

  inline void
  acceptanceCutSet::setLeptonCut( double inputValue )
  {

    setPrimaryLeptonCut( inputValue );
    setSecondaryLeptonCut( inputValue );

  }


  inline double
  acceptanceCutSet::getJetCut()
  const
  {

    return jetCut;

  }

  inline void
  acceptanceCutSet::setJetCut( double inputValue )
  {

    jetCut = inputValue;

  }


  inline bool
  acceptanceCutSet::isSameAcceptanceCutSet(
                              acceptanceCutSet const* const comparisonCuts )
  const
  {

    if( ( comparisonCuts->getBeamEnergy() == beamEnergy )
        &&
        ( comparisonCuts->getJetCut() == jetCut )
        &&
        ( comparisonCuts->getPrimaryLeptonCut() == primaryLeptonCut )
        &&
        ( comparisonCuts->getSecondaryLeptonCut()
          == secondaryLeptonCut ) )
    {

      return true;

    }
    else
    {

      return false;

    }

  }



  inline bool
  singleLeptonCut::isSameAcceptanceCutSet(
                              acceptanceCutSet const* const comparison_cuts )
  const
  {

    if( ( comparison_cuts->getBeamEnergy() == beamEnergy )
        &&
        ( comparison_cuts->getPrimaryLeptonCut() == primaryLeptonCut )  )
    {

      return true;

    }
    else
    {

      return false;

    }

  }



  inline double
  leptonAcceptanceParameterSet::getEffectiveSquarkMass()
  {
    // debugging:
    /**std::cout
    << std::endl
    << "debugging: leptonAcceptanceParameterSet::getEffectiveSquarkMass() called.";
    std::cout << std::endl;**/

    if( needs_to_prepare_for_this_point() )
    {
      resetValues();
    }
    return effectiveSquarkMass;
  }

  inline double
  leptonAcceptanceParameterSet::acceptanceAt( double const givenEnergy,
                                     double const givenCut )
  /* this checks to see if the acceptances need updating, then returns
   * calculateAcceptanceAt( givenEnergy,
   *                        givenCut ), which interpolates the values in
   * acceptance_bins to the requested value, or returns
   * pseudorapidityAcceptance if it's lower, scaled to the given value for
   * the transverse momentum cut.
   */
  {
    // debugging:
    /**std::cout
    << std::endl
    << "debugging: leptonAcceptanceParameterSet::acceptanceAt( " << givenEnergy
    << ", " << givenCut << " ) called.";
    std::cout << std::endl;**/

    if( needs_to_prepare_for_this_point() )
    {
      resetValues();
    }
    return calculateAcceptanceAt( givenEnergy,
                                  givenCut );
  }


  inline int
  acceptanceGrid::getNumberOfAcceptanceColumns()
  const
  {

    return acceptanceColumns;

  }

  inline double
  acceptanceGrid::getLowestSquarkMass()
  const
  {

    return lowestSquarkMass;

  }

  inline double
  acceptanceGrid::getHighestSquarkMass()
  const
  {

    return  highestSquarkMass;

  }

  inline double
  acceptanceGrid::getLowestGluinoMass()
  const
  {

    return  lowestGluinoMass;

  }

  inline double
  acceptanceGrid::getHighestGluinoMass()
  const
  {

    return  highestGluinoMass;

  }



  inline double
  jet_acceptance_value::get_acceptance()
  {

    if( not_already_calculated )
      {

        stored_acceptance = jet_table->getValue( scolored_pair,
                                                  first_sQCD_cascade,
                                                  second_sQCD_cascade,
                                                  acceptance_column,
                                                  false,
                                                  false );
        not_already_calculated = false;

      }

    return stored_acceptance;

  }


  inline bool
  jet_acceptance_value::is_requested( int const given_acceptance_column,
                         signed_particle_shortcut_pair const* const given_pair,
                         colored_cascade const* const given_first_sQCD_cascade,
                       colored_cascade const* const given_second_sQCD_cascade )
  const
  {

    if( ( given_acceptance_column == acceptance_column )
        &&
        ( given_pair == scolored_pair )
        &&
        ( given_first_sQCD_cascade == first_sQCD_cascade )
        &&
        ( given_second_sQCD_cascade == second_sQCD_cascade ) )
      {

        return true;

      }
    else
      {

        return false;

      }

  }



  inline bool
  lepton_acceptance_value::is_requested(
                         signed_particle_shortcut_pair const* const given_pair,
                         colored_cascade const* const given_first_sQCD_cascade,
                       colored_cascade const* const given_second_sQCD_cascade )
  const
  {

    if( ( given_pair == scolored_pair )
        &&
        ( given_first_sQCD_cascade == first_sQCD_cascade )
        &&
        ( given_second_sQCD_cascade == second_sQCD_cascade ) )
      {

        return true;

      }
    else
      {

        return false;

      }

  }



  inline std::string const*
  jet_acceptance_table::get_name()
  const
  {

    return &table_name;

  }



  inline int
  lepton_acceptance_table::get_energy()
  const
  {

    return LHC_energy_in_TeV;

  }



  inline jet_acceptance_value*
  kinematics_table::get_jet_acceptance(
                         signed_particle_shortcut_pair const* const given_pair,
                         colored_cascade const* const given_first_sQCD_cascade,
                       colored_cascade const* const given_second_sQCD_cascade )
  {

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: kinematics_table::get_jet_acceptance(...) called,";
    std::cout << std::endl;**/

    return jets_table->get_acceptance( jet_acceptance_column,
                                       given_pair,
                                       given_first_sQCD_cascade,
                                       given_second_sQCD_cascade );

  }

  inline lepton_acceptance_value*
  kinematics_table::get_lepton_acceptance(
                       signed_particle_shortcut_pair const* const given_pair,
                       colored_cascade const* const given_first_sQCD_cascade,
                     colored_cascade const* const given_second_sQCD_cascade )
  {

    return leptons_table->get_acceptance( given_pair,
                                          given_first_sQCD_cascade,
                                          given_second_sQCD_cascade );

  }


  inline int
  kinematics_table::get_column()
  const
  {

    return jet_acceptance_column;

  }



  inline bool
  kinematics_table_set::is_requested( int const given_LHC_energy_in_TeV,
                              std::string const* const given_jet_subdirectory )
  {

    if( ( given_LHC_energy_in_TeV == LHC_energy_in_TeV )
        &&
        ( 0 == given_jet_subdirectory->compare( jet_grid_subdirectory ) ) )
      {

        return true;

      }
    else
      {

        return false;

      }

  }



  inline kinematics_table*
  kinematics_handler::get_table( int const LHC_energy_in_TeV,
                               std::string const* const given_jet_subdirectory,
                                 int const given_acceptance_column )
  {

    return get_table_set( LHC_energy_in_TeV,
                given_jet_subdirectory )->get_table( given_acceptance_column );

  }

}  // end of LHC_FASER namespace.

#endif /* LHC_FASER_JET_KINEMATICS_STUFF_HPP_ */
