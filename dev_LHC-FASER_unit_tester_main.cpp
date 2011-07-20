/*
 * dev_LHC-FASER_unit_tester_main.cpp
 *
 *  Created on: 21 Jan 2011
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

/* this program is to test various bits of LHC-FASER in isolation (or as much
 * isolation as practical).
 *
 * I should really test each bit as it is written, BUT the whole structure of
 * LHC-FASER is subject to changes, even major changes, leaving good chances
 * of testing units that do not find their way into even the 1st prototype
 * compilable version.
 *
 *
 * the current ordering for testing is:
 *
 * (0: CppSLHA - done.)
 * 1: lepton distributions - sorta redone, but for just a subset of cases
 *                           (framework works though, being the point!).
 *                         - @Carsten: please check that the new W & Z
 *                           distributions work.
 * 2: readied_for_new_point + readier_for_new_point - not yet redone, but
 *                                                    should not be a problem,
 *                                                    since I don't *think*
 *                                                    that I changed anything
 *                                                    here since 0.1 (which
 *                                                    works)...
 * 3: cross_section_table - again, not yet redone, but I don't recall changing
 *                          anything here since the working 0.1 version.
 * 4: kinematics_table - not yet redone, & it'll be quite different in 0.2,
 *                       unfortunately.
 * 5: lepton acceptance for given cascade - not yet redone
 * 6: an entire signal - not yet redone
 *
 * after that, I think that we'll be ready to test the whole thing.
 */


// includes:

#include "CppSLHA/CppSLHA.hpp"
#include "LHC-FASER_derived_lepton_distributions.hpp"
#include "LHC-FASER_input_handling_stuff.hpp"
#include "LHC-FASER_cross-section_stuff.hpp"
#include "LHC-FASER_kinematics_stuff.hpp"
//#include "LHC-FASER_derived_electroweak_cascade_stuff.hpp"
//#include "LHC-FASER_signal_data_collection_stuff.hpp"
//#include "LHC-FASER_signal_calculator_stuff.hpp"
//#include "LHC-FASER.hpp"
#include <sys/time.h>


// this holds the set of distributions & handles adding them.
class distribution_set
{

public:

  distribution_set( LHC_FASER::readier_for_new_point* const readier ) :
    readier( readier ),
    hard_muon( new LHC_FASER::hard_muon_from_tau() ),
    soft_muon( new LHC_FASER::soft_muon_from_tau() ),
    hard_pion( new LHC_FASER::hard_pion_from_tau() ),
    soft_pion( new LHC_FASER::soft_pion_from_tau() )
  {

    // just an initialization list.

  }

  virtual
  ~distribution_set()
  {

    for( std::vector< LHC_FASER::leptonEnergyDistribution* >::iterator
           deletion_iterator = distributions_vector.begin();
         distributions_vector.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

    for( std::vector< std::string* >::iterator
           deletion_iterator = names_vector.begin();
         names_vector.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

  }


  void
  add_distributions(
         LHC_FASER::leptonEnergyDistribution* const given_distribution,
                     std::string given_name )
  {

    std::string* test_name = new std::string( given_name );
    test_name->append( "_ll" );
    distributions_vector.push_back( given_distribution );
    names_vector.push_back( test_name );

    LHC_FASER::leptonEnergyDistribution* tau_distribution
    = new LHC_FASER::visible_tau_decay_product( readier,
                                                given_distribution,
                                                hard_muon );
    test_name = new std::string( given_name );
    test_name->append( "_tLm" );
    distributions_vector.push_back( tau_distribution );
    names_vector.push_back( test_name );
    tau_distribution
    = new LHC_FASER::visible_tau_decay_product( readier,
                                                given_distribution,
                                                soft_muon );
    test_name = new std::string( given_name );
    test_name->append( "_tRm" );
    distributions_vector.push_back( tau_distribution );
    names_vector.push_back( test_name );
    tau_distribution
    = new LHC_FASER::visible_tau_decay_product( readier,
                                                given_distribution,
                                                soft_pion );
    test_name = new std::string( given_name );
    test_name->append( "_tLp" );
    distributions_vector.push_back( tau_distribution );
    names_vector.push_back( test_name );
    tau_distribution
    = new LHC_FASER::visible_tau_decay_product( readier,
                                                given_distribution,
                                                hard_pion );
    test_name = new std::string( given_name );
    test_name->append( "_tRp" );
    distributions_vector.push_back( tau_distribution );
    names_vector.push_back( test_name );

  }

  std::vector< LHC_FASER::leptonEnergyDistribution* > const*
  get_distributions_vector()
  const
  {

    return &distributions_vector;

  }

  std::vector< std::string* > const*
  get_names_vector()
  const
  {

    return &names_vector;

  }


protected:
  LHC_FASER::readier_for_new_point* readier;
  LHC_FASER::tau_decay_coefficient const* const hard_muon;
  LHC_FASER::tau_decay_coefficient const* const soft_muon;
  LHC_FASER::tau_decay_coefficient const* const hard_pion;
  LHC_FASER::tau_decay_coefficient const* const soft_pion;

  std::vector< LHC_FASER::leptonEnergyDistribution* >
  distributions_vector;
  std::vector< std::string* > names_vector;

};

class squarkMassReturner : public LHC_FASER::effectiveSquarkMassHolder
{
public:
  squarkMassReturner( CppSLHA::particle_property_set const* const squark ) :
    squark( squark )
  {
    // just an initialization list.
  }

  virtual
  ~squarkMassReturner()
  {
    // does nothing.
  }

  double
  getEffectiveSquarkMass()
  {
    return squark->get_absolute_mass();
  }

protected:
  CppSLHA::particle_property_set const* const squark;
};


// the main program:

int main( int argument_count,
          char* argument_strings[] )
{

  int recalculation_amount = 10;

  if( 2 != argument_count )
    // if the input was definitely wrong...
    {

      std::cout
      << std::endl
      << "error! this program requires an integer as its argument, as the"
      << " number of times to re-calculate all the points. defaulting to 10.";
      std::cout << std::endl;  // let the user know the format.

    }
  else  // otherwise, do what is supposed to be done.
    {

      std::string argument_interpretter( argument_strings[ 1 ] );
      std::stringstream argument_stream( argument_interpretter );
      argument_stream >> recalculation_amount;

    }

  CppSLHA::CppSLHA2 SLHA_data( "SPS1a_spectrum.out" );
  SLHA_data.read_file();
  LHC_FASER::readier_for_new_point test_readier;
  LHC_FASER::input_handler test_input_handler( &SLHA_data,
                                         SLHA_data.get_particle_spectrum(),
                                         "./grids",
                                         &test_readier );

  std::string bluh( "bluh" );

  /* testing the lepton energy distributions in the squark rest frame: (slight changes of input) *//**/

  distribution_set::distribution_set distributions( &test_readier );

  CppSLHA::EW_scale_spectrum* const
  spectrum = SLHA_data.get_particle_spectrum();
  CppSLHA::particle_property_set const* const
  sdown_L_pointer
  = spectrum->inspect_particle_property_set( CppSLHA::PDG_code::sdown_L );
  squarkMassReturner sdownLMass( sdown_L_pointer );
  CppSLHA::particle_property_set const* const
  neutralino_four_pointer
  = spectrum->inspect_particle_property_set(
                                          CppSLHA::PDG_code::neutralino_four );
  CppSLHA::particle_property_set* const
  neutralino_three_pointer
  = spectrum->get_particle_property_set( CppSLHA::PDG_code::neutralino_three );
  CppSLHA::particle_property_set const* const
  neutralino_two_pointer
  = spectrum->inspect_particle_property_set(
                                           CppSLHA::PDG_code::neutralino_two );
  CppSLHA::particle_property_set const* const
  chargino_two_pointer
  = spectrum->inspect_particle_property_set( CppSLHA::PDG_code::chargino_two );
  CppSLHA::particle_property_set const* const
  chargino_one_pointer
  = spectrum->inspect_particle_property_set( CppSLHA::PDG_code::chargino_one );
  CppSLHA::particle_property_set const* const
  smuon_L_pointer
  = spectrum->inspect_particle_property_set( CppSLHA::PDG_code::smuon_L );
  CppSLHA::particle_property_set const* const
  smuon_R_pointer
  = spectrum->inspect_particle_property_set( CppSLHA::PDG_code::smuon_R );
  CppSLHA::particle_property_set const* const
  muon_sneutrino_L_pointer
  = spectrum->inspect_particle_property_set(
                                         CppSLHA::PDG_code::muon_sneutrino_L );
  CppSLHA::particle_property_set const* const
  muon_sneutrino_R_pointer
  = spectrum->inspect_particle_property_set(
                                         CppSLHA::PDG_code::muon_sneutrino_R );
  CppSLHA::particle_property_set const* const
  light_neutral_EWSB_scalar_pointer
  = spectrum->inspect_particle_property_set(
                                CppSLHA::PDG_code::light_neutral_EWSB_scalar );
  CppSLHA::particle_property_set* const
  charged_EWSB_scalar_pointer
  = spectrum->get_particle_property_set(
                                      CppSLHA::PDG_code::charged_EWSB_scalar );
  CppSLHA::particle_property_set const* const
  Z_pointer
  = spectrum->inspect_particle_property_set( CppSLHA::PDG_code::Z );
  CppSLHA::particle_property_set const* const
  W_pointer
  = spectrum->inspect_particle_property_set( CppSLHA::PDG_code::W_plus );
  CppSLHA::particle_property_set const* const
  neutralino_one_pointer
  = spectrum->inspect_particle_property_set(
                                           CppSLHA::PDG_code::neutralino_one );

  // here we set n3 to have enough mass to decay to a h/Z + n1:
  neutralino_three_pointer->set_mass( ( 20.0
                       + light_neutral_EWSB_scalar_pointer->get_absolute_mass()
                             + neutralino_one_pointer->get_absolute_mass() ) );

  // here we set H+ to be light enough for x2 to decay into it plus n1:
  charged_EWSB_scalar_pointer->set_mass( (
                                      chargino_two_pointer->get_absolute_mass()
                                  - neutralino_one_pointer->get_absolute_mass()
                                           - 20.0 ) );

  LHC_FASER::leptonEnergyDistribution*
  currentDistribution
  = new LHC_FASER::flat_near_muon_plus_antimuon( &test_readier,
                                                 &SLHA_data,
                                                 sdown_L_pointer,
                                                 &sdownLMass,
                                                 neutralino_two_pointer,
                                                 smuon_R_pointer );
  distributions.add_distributions( currentDistribution,
                                   "flat_near" );
  currentDistribution
  = new LHC_FASER::same_chirality_near_muon( &test_readier,
                                             &SLHA_data,
                                             sdown_L_pointer,
                                             &sdownLMass,
                                             neutralino_two_pointer,
                                             smuon_R_pointer );
  distributions.add_distributions( currentDistribution,
                                   "same_near" );
  currentDistribution
  = new LHC_FASER::opposite_chirality_near_muon( &test_readier,
                                                 &SLHA_data,
                                                 sdown_L_pointer,
                                                 &sdownLMass,
                                                 neutralino_two_pointer,
                                                 smuon_R_pointer );
  distributions.add_distributions( currentDistribution,
                                   "opp_near" );
  currentDistribution
  = new LHC_FASER::flat_far_muon_plus_antimuon( &test_readier,
                                                &SLHA_data,
                                                sdown_L_pointer,
                                                &sdownLMass,
                                                neutralino_two_pointer,
                                                smuon_R_pointer,
                                                neutralino_one_pointer );
  distributions.add_distributions( currentDistribution,
                                   "flat_far" );
  currentDistribution
  = new LHC_FASER::same_chirality_far_muon( &test_readier,
                                            &SLHA_data,
                                            sdown_L_pointer,
                                            &sdownLMass,
                                            neutralino_two_pointer,
                                            smuon_R_pointer,
                                            neutralino_one_pointer );
  distributions.add_distributions( currentDistribution,
                                   "same_far" );
  currentDistribution
  = new LHC_FASER::opposite_chirality_far_muon( &test_readier,
                                                &SLHA_data,
                                                sdown_L_pointer,
                                                &sdownLMass,
                                                neutralino_two_pointer,
                                                smuon_R_pointer,
                                                neutralino_one_pointer );
  distributions.add_distributions( currentDistribution,
                                   "opp_far" );
  currentDistribution
  = new LHC_FASER::Higgs_muon_plus_antimuon( &test_readier,
                                             &SLHA_data,
                                             sdown_L_pointer,
                                             &sdownLMass,
                                             neutralino_four_pointer,
                                             light_neutral_EWSB_scalar_pointer,
                                             neutralino_one_pointer );
  distributions.add_distributions( currentDistribution,
                                   "hn4" );
  currentDistribution
  = new LHC_FASER::Higgs_muon_plus_antimuon( &test_readier,
                                             &SLHA_data,
                                             sdown_L_pointer,
                                             &sdownLMass,
                                             neutralino_three_pointer,
                                             light_neutral_EWSB_scalar_pointer,
                                             neutralino_one_pointer );
  distributions.add_distributions( currentDistribution,
                                   "hn3" );
  currentDistribution
  = new LHC_FASER::Higgs_muon_plus_antimuon( &test_readier,
                                             &SLHA_data,
                                             sdown_L_pointer,
                                             &sdownLMass,
                                             chargino_two_pointer,
                                             charged_EWSB_scalar_pointer,
                                             neutralino_one_pointer );
  distributions.add_distributions( currentDistribution,
                                   "H+" );
  /*currentDistribution
  = new LHC_FASER::Z_muon_plus_antimuon( &test_readier,
                                         &SLHA_data,
                                         sdown_L_pointer,
                                         &sdownLMass,
                                         neutralino_three_pointer,
                                         Z_pointer,
                                         neutralino_one_pointer );
  distributions.add_distributions( currentDistribution,
                                   "Zn4" );
  currentDistribution
  = new LHC_FASER::Z_muon_plus_antimuon( &test_readier,
                                         &SLHA_data,
                                         sdown_L_pointer,
                                         &sdownLMass,
                                         neutralino_four_pointer,
                                         Z_pointer,
                                         neutralino_one_pointer );
  distributions.add_distributions( currentDistribution,
                                   "Zn3" );
  currentDistribution
  = new LHC_FASER::W_minus_plus_up_L_type_muon( &test_readier,
                                                &SLHA_data,
                                                sdown_L_pointer,
                                                &sdownLMass,
                                                chargino_one_pointer,
                                                W_pointer,
                                                neutralino_one_pointer );
  distributions.add_distributions( currentDistribution,
                                   "LWx1" );
  currentDistribution
  = new LHC_FASER::W_minus_plus_up_R_type_muon( &test_readier,
                                                &SLHA_data,
                                                sdown_L_pointer,
                                                &sdownLMass,
                                                chargino_one_pointer,
                                                W_pointer,
                                                neutralino_one_pointer );
  distributions.add_distributions( currentDistribution,
                                   "RWx1" );
  currentDistribution
  = new LHC_FASER::W_minus_plus_up_L_type_muon( &test_readier,
                                                &SLHA_data,
                                                sdown_L_pointer,
                                                &sdownLMass,
                                                chargino_two_pointer,
                                                W_pointer,
                                                neutralino_one_pointer );
  distributions.add_distributions( currentDistribution,
                                   "LWx2" );
  currentDistribution
  = new LHC_FASER::W_minus_plus_up_R_type_muon( &test_readier,
                                                &SLHA_data,
                                                sdown_L_pointer,
                                                &sdownLMass,
                                                chargino_two_pointer,
                                                W_pointer,
                                                neutralino_one_pointer );
  distributions.add_distributions( currentDistribution,
                                   "RWx2" );
  currentDistribution
  = new LHC_FASER::neutralino_three_body_decay( &test_readier,
                                                &SLHA_data,
                                                sdown_L_pointer,
                                                &sdownLMass,
                                                neutralino_two_pointer,
                                                neutralino_one_pointer,
                                                smuon_L_pointer,
                                                smuon_R_pointer );
  distributions.add_distributions( currentDistribution,
                                   "nn3b" );
  currentDistribution
  = new LHC_FASER::chargino_three_body_decay( &test_readier,
                                              &SLHA_data,
                                              sdown_L_pointer,
                                              &sdownLMass,
                                              chargino_one_pointer,
                                              neutralino_one_pointer,
                                              smuon_L_pointer,
                                              smuon_R_pointer,
                                              muon_sneutrino_L_pointer,
                                              muon_sneutrino_R_pointer );
  distributions.add_distributions( currentDistribution,
                                   "xx3b" );*/


  std::ofstream* gnuplot_file = new std::ofstream();
  gnuplot_file->open( "./lepton_distributions/test_plotter.input" );


  if( gnuplot_file->is_open() )
    // if the file was successfully opened...
    {

      *gnuplot_file << "set term postscript enhanced color solid" << std::endl;

      std::string file_name;

      for( int graph_maker = 0;
           distributions.get_distributions_vector()->size() > graph_maker;
           ++graph_maker )
        {

          // open a file for writing out the data points.
          std::ofstream* output_file = new std::ofstream();
          file_name.assign( "./lepton_distributions/" );
          file_name.append(
                      *(distributions.get_names_vector()->at( graph_maker )) );
          file_name.append( ".dat" );
          output_file->open( file_name.c_str() );

          if( output_file->is_open() )
            // if the file was successfully opened...
            {

              *gnuplot_file
              << "set output \""
              << *(distributions.get_names_vector()->at( graph_maker ))
              << ".eps\"" << std::endl
              << "set style line 1 lt rgb \"red\" lw 3" << std::endl
              << "plot '"
              << *(distributions.get_names_vector()->at( graph_maker ))
              << ".dat' index 0 notitle with lines ls 1" << std::endl;

              std::cout
              << std::endl
              << *(distributions.get_names_vector()->at( graph_maker ))
              << ": "
              << distributions.get_distributions_vector()->at( graph_maker
                                                          )->getMinimumEnergy()
              << " to "
              << distributions.get_distributions_vector()->at( graph_maker
                                                         )->getMaximumEnergy();
              std::cout << std::endl;


              int const number_of_bins = 100;

              double
              bin_size
              = ( ( distributions.get_distributions_vector()->at( graph_maker
                                                          )->getMaximumEnergy()
                    - distributions.get_distributions_vector()->at( graph_maker
                                                        )->getMinimumEnergy() )
                  / number_of_bins );

              double bin_position;

              double plot_value;

              double area = 0.0;

              for( int bin_counter = 0;
                   number_of_bins >= bin_counter;
                   ++bin_counter )
                {

                  bin_position
                  = ( (double)bin_counter * bin_size
                      + distributions.get_distributions_vector()->at(
                                           graph_maker )->getMinimumEnergy() );

                  plot_value
                  = distributions.get_distributions_vector()->at( graph_maker
                                                   )->valueAt( bin_position );

                  area += ( plot_value * bin_size );

                  *output_file
                  << bin_position << " " << plot_value << std::endl;

                }

              output_file->close();

              std::cout << "area = " << area;
              std::cout << std::endl;

            }
          else
            {

              std::cout
              << std::endl << "error!  output file could not be opened!";

              exit( EXIT_FAILURE );

            }

          delete output_file;

        }

      gnuplot_file->close();

    }
  else
    {

      std::cout << std::endl << "error! gnuplot file could not be opened!";

      exit( EXIT_FAILURE );

    }

  // clean up:

  delete gnuplot_file;
  /**/


  /* testing the means of checking if something has already been prepared for a
   * new point: (should work)*//**

  std::cout
  << std::endl
  << "is verbose? " << test_input_handler.is_verbose()
  << std::endl
  << "set to be not verbose.";
  test_input_handler.set_verbosity( false );
  std::cout
  << std::endl
  << "is verbose? " << test_input_handler.is_verbose()
  << std::endl
  << "set to be verbose.";
  test_input_handler.set_verbosity( true );
  std::cout
  << std::endl
  << "is verbose? " << test_input_handler.is_verbose()
  << std::endl
  << "using NLO? " << test_input_handler.using_NLO()
  << std::endl
  << "set to not use NLO.";
  test_input_handler.set_using_NLO( false );
  std::cout
  << std::endl
  << "using NLO? " << test_input_handler.using_NLO()
  << std::endl
  << "set to use NLO.";
  test_input_handler.set_using_NLO( true );
  std::cout
  << std::endl
  << "using NLO? " << test_input_handler.using_NLO()
  << std::endl
  << std::endl
  << "path_to_grids = "
  << *(test_input_handler.get_path_to_grids())
  << std::endl
  << "path_to_cross_section_grids = "
  << *(test_input_handler.get_path_to_cross_section_grids())
  << std::endl
  << "(inspect: path_to_cross_section_grids = "
  << *(test_input_handler.inspect_path_to_cross_section_grids()) << ")"
  << std::endl
  << "path_to_kinematics_grids = "
  << *(test_input_handler.get_path_to_kinematics_grids())
  << std::endl
  << "(inspect: path_to_kinematics_grids = "
  << *(test_input_handler.inspect_path_to_kinematics_grids()) << ")"
  << std::endl
  << std::endl
  << "SPS1a MASS:"
  << test_input_handler.get_CppSLHA()->get_MASS_as_string()
  << std::endl
  << std::endl
  << *(test_input_handler.get_spectrum()->get_PDG_coded_name(
                                                        CppSLHA::PDG_code::Z ))
  << " "
  << test_input_handler.get_spectrum()->get_PDG_coded_mass(
                                                         CppSLHA::PDG_code::Z )
  << " GeV"
  << std::endl
  << *(test_input_handler.get_W_plus()->get_name()) << " "
  << test_input_handler.get_W_plus()->get_mass() << " GeV"
  << std::endl
  << *(test_input_handler.get_top()->get_name()) << " "
  << test_input_handler.get_top()->get_mass() << " GeV"
  << std::endl
  << *(test_input_handler.get_gluino()->get_name()) << " "
  << test_input_handler.get_gluino()->get_mass() << " GeV"
  << std::endl
  << "average_squarks4_mass = "
  << test_input_handler.get_average_squarks4_mass()
  << std::endl
  << "scoloreds: " << std::endl;
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = test_input_handler.get_colored_sparticles()->begin();
       test_input_handler.get_colored_sparticles()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->get_name()) << std::endl;
    }
  std::cout
  << std::endl
  << "charginos: " << std::endl;
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = test_input_handler.get_charginos()->begin();
       test_input_handler.get_charginos()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->get_name()) << std::endl;
    }
  std::cout
  << std::endl
  << "neutralinos: " << std::endl;
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = test_input_handler.get_neutralinos()->begin();
       test_input_handler.get_neutralinos()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->get_name()) << std::endl;
    }
  std::cout
  << std::endl
  << "unstable neutralinos: " << std::endl;
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = test_input_handler.get_unstable_neutralinos()->begin();
       test_input_handler.get_unstable_neutralinos()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->get_name()) << std::endl;
    }
  std::cout
  << std::endl
  << "electroweakinos: " << std::endl;
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = test_input_handler.get_electroweakinos()->begin();
       test_input_handler.get_electroweakinos()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->get_name()) << std::endl;
    }
  std::cout
  << std::endl
  << "unstable electroweakinos: " << std::endl;
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator
       = test_input_handler.get_unstable_electroweakinos()->begin();
       test_input_handler.get_unstable_electroweakinos()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->get_name()) << std::endl;
    }
  std::cout
  << std::endl
  << "charged sleptons: " << std::endl;
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = test_input_handler.get_charged_sleptons()->begin();
       test_input_handler.get_charged_sleptons()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->get_name()) << std::endl;
    }
  std::cout
  << std::endl
  << "sneutrinos: " << std::endl;
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = test_input_handler.get_sneutrinos()->begin();
       test_input_handler.get_sneutrinos()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->get_name()) << std::endl;
    }
  std::cout
  << std::endl
  << "sdowns: " << std::endl;
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = test_input_handler.get_sdown_types()->begin();
       test_input_handler.get_sdown_types()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->get_name()) << std::endl;
    }
  std::cout
  << std::endl
  << "sups: " << std::endl;
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = test_input_handler.get_sup_types()->begin();
       test_input_handler.get_sup_types()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->get_name()) << std::endl;
    }
  std::cout
  << std::endl
  << "jets: " << std::endl;
  for( std::list< int >::const_iterator
       set_iterator = test_input_handler.get_jets()->begin();
       test_input_handler.get_jets()->end() != set_iterator;
       ++set_iterator )
    {
      std::cout
      << *(test_input_handler.get_spectrum()->get_PDG_coded_name(
                                                             (*set_iterator) ))
      << std::endl;
    }
  std::cout
  << std::endl
  << "jets4: " << std::endl;
  for( std::list< int >::const_iterator
       set_iterator = test_input_handler.get_jets4()->begin();
       test_input_handler.get_jets4()->end() != set_iterator;
       ++set_iterator )
    {
      std::cout
      << *(test_input_handler.get_spectrum()->get_PDG_coded_name(
                                                             (*set_iterator) ))
      << std::endl;
    }
  std::cout
  << std::endl
  << "jets5: " << std::endl;
  for( std::list< int >::const_iterator
       set_iterator = test_input_handler.get_jets5()->begin();
       test_input_handler.get_jets5()->end() != set_iterator;
       ++set_iterator )
    {
      std::cout
      << *(test_input_handler.get_spectrum()->get_PDG_coded_name(
                                                             (*set_iterator) ))
      << std::endl;
    }
  std::cout
  << std::endl
  << "not in jets4: " << std::endl;
  for( std::list< int >::const_iterator
       set_iterator = test_input_handler.get_not_in_jets4()->begin();
       test_input_handler.get_not_in_jets4()->end() != set_iterator;
       ++set_iterator )
    {
      std::cout
      << *(test_input_handler.get_spectrum()->get_PDG_coded_name(
                                                             (*set_iterator) ))
      << std::endl;
    }
  std::cout
  << std::endl
  << "not in jets5: " << std::endl;
  for( std::list< int >::const_iterator
       set_iterator = test_input_handler.get_not_in_jets5()->begin();
       test_input_handler.get_not_in_jets5()->end() != set_iterator;
       ++set_iterator )
    {
      std::cout
      << *(test_input_handler.get_spectrum()->get_PDG_coded_name(
                                                             (*set_iterator) ))
      << std::endl;
    }
  std::cout
  << std::endl
  << "squarks: " << std::endl;
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = test_input_handler.get_squarks()->begin();
       test_input_handler.get_squarks()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->get_name()) << std::endl;
    }
  std::cout
  << std::endl
  << "squarks4: " << std::endl;
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = test_input_handler.get_squarks4()->begin();
       test_input_handler.get_squarks4()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->get_name()) << std::endl;
    }
  std::cout
  << std::endl
  << "squarks5: " << std::endl;
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = test_input_handler.get_squarks5()->begin();
       test_input_handler.get_squarks5()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->get_name()) << std::endl;
    }
  std::cout
  << std::endl
  << "sbottom_1 is in squarks? "
  << test_input_handler.is_in(
                          test_input_handler.get_sbottom_one()->get_PDG_code(),
                               test_input_handler.get_squarks() )
  << std::endl
  << "sbottom_1 is in squarks4? "
  << test_input_handler.is_in(
                          test_input_handler.get_sbottom_one()->get_PDG_code(),
                               test_input_handler.get_squarks4() )
  << std::endl
  << std::endl
  << "colored cascade set: ";
  LHC_FASER::colored_cascade_set*
  gluino_cascade
  = test_input_handler.get_colored_cascade_set(
                                             test_input_handler.get_gluino() );
  LHC_FASER::colored_cascade_set*
  sbottom_cascade
  = test_input_handler.get_colored_cascade_set(
                                        test_input_handler.get_sbottom_one() );
  LHC_FASER::colored_cascade_set*
  gluino_cascade2
  = test_input_handler.get_colored_cascade_set(
                                             test_input_handler.get_gluino() );
  std::cout
  << *(gluino_cascade->get_decayer()->get_name())
  << ":" << std::endl
  << "pointer = " << gluino_cascade << ", pointer2 = " << gluino_cascade2
  << std::endl;
  for( std::vector< LHC_FASER::colored_cascade* >::const_iterator
       set_iterator = gluino_cascade->get_colored_cascades()->begin();
       gluino_cascade->get_colored_cascades()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout
      << *((*set_iterator)->get_decayer()->get_name())
      << " -> " << *((*set_iterator)->get_product()->get_name())
      << " ( " << (*set_iterator)->get_product_position() << " ) -> "
      << *((*set_iterator)->get_EWino()->get_name())
      << " ( " << (*set_iterator)->get_EWino_position() << " )";
      if( (*set_iterator)->is_sx() )
        {
          std::cout<< " => sx" << std::endl;
        }
      if( (*set_iterator)->is_gx() )
        {
          std::cout<< " => gx" << std::endl;
        }
      if( (*set_iterator)->is_sgx() )
        {
          std::cout<< " => sgx" << std::endl;
        }
      if( (*set_iterator)->is_gsx() )
        {
          std::cout<< " => gsx" << std::endl;
        }
      if( (*set_iterator)->is_sgsx() )
        {
          std::cout<< " => sgsx" << std::endl;
        }
    }
  std::cout
  << *(sbottom_cascade->get_decayer()->get_name())
  << ":" << std::endl
  << "pointer = " << sbottom_cascade
  << std::endl;
  for( std::vector< LHC_FASER::colored_cascade* >::const_iterator
       set_iterator = sbottom_cascade->get_colored_cascades()->begin();
       sbottom_cascade->get_colored_cascades()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout
      << *((*set_iterator)->get_decayer()->get_name())
      << " -> " << *((*set_iterator)->get_product()->get_name())
      << " ( " << (*set_iterator)->get_product_position() << " ) -> "
      << *((*set_iterator)->get_EWino()->get_name())
      << " ( " << (*set_iterator)->get_EWino_position() << " )";
      if( (*set_iterator)->is_sx() )
        {
          std::cout<< " => sx" << std::endl;
        }
      if( (*set_iterator)->is_gx() )
        {
          std::cout<< " => gx" << std::endl;
        }
      if( (*set_iterator)->is_sgx() )
        {
          std::cout<< " => sgx" << std::endl;
        }
      if( (*set_iterator)->is_gsx() )
        {
          std::cout<< " => gsx" << std::endl;
        }
      if( (*set_iterator)->is_sgsx() )
        {
          std::cout<< " => sgsx" << std::endl;
        }
    }
  std::cout
  << std::endl
  << std::endl
  << "sparticle production combinations: " << std::endl;
  for( std::vector< LHC_FASER::signed_particle_shortcut_pair* >::const_iterator
       set_iterator
       = test_input_handler.get_sparticle_production_combinations()->begin();
       test_input_handler.get_sparticle_production_combinations()->end()
       > set_iterator;
       ++set_iterator )
    {
      if( (*set_iterator)->is_first_particle() )
        {
          std::cout << *((*set_iterator)->get_first_pointer()->get_name());
        }
      else
        {
          std::cout << *((*set_iterator)->get_first_pointer()->get_antiname());
        }
      std::cout << " + ";
      if( (*set_iterator)->is_second_particle() )
        {
          std::cout << *((*set_iterator)->get_second_pointer()->get_name());
        }
      else
        {
          std::cout
          << *((*set_iterator)->get_second_pointer()->get_antiname());
        }
      std::cout << std::endl;
    }
  std::vector< CppSLHA::particle_property_set const* > sdown_LR;
  sdown_LR.push_back( test_input_handler.get_sdown_L() );
  sdown_LR.push_back( test_input_handler.get_sdown_R() );
  std::vector< CppSLHA::particle_property_set const* > sup_LR;
  sup_LR.push_back( test_input_handler.get_sup_L() );
  sup_LR.push_back( test_input_handler.get_sup_R() );
  LHC_FASER::direct_decay_checker*
  gluino_to_sdown_LR
  = test_input_handler.get_direct_decay_checker(
                                               test_input_handler.get_gluino(),
                                                 &sdown_LR );
  LHC_FASER::direct_decay_checker*
  gluino_to_sup_LR
  = test_input_handler.get_direct_decay_checker(
                                               test_input_handler.get_gluino(),
                                                 &sup_LR );
  LHC_FASER::direct_decay_checker*
  gluino_to_sdown_LR2
  = test_input_handler.get_direct_decay_checker(
                                               test_input_handler.get_gluino(),
                                                 &sdown_LR );
  std::cout
  << std::endl
  << "direct_decay_checker: gluino -> sdown_L, sdown_R"
  << std::endl << "pointer = "
  << gluino_to_sdown_LR
  << std::endl << *(gluino_to_sdown_LR->get_decayer()->get_name())
  << " -> ";
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = gluino_to_sdown_LR->get_products()->begin();
       gluino_to_sdown_LR->get_products()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->get_name()) << " ";
    }
  std::cout
  << std::endl
  << "has direct decay? " << gluino_to_sdown_LR->has_direct_decay();
  std::cout
  << std::endl
  << "2nd gluino -> sdown_L, sdown_R"
  << std::endl << "pointer = "
  << gluino_to_sdown_LR2
  << std::endl << *(gluino_to_sdown_LR2->get_decayer()->get_name())
  << " -> ";
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = gluino_to_sdown_LR2->get_products()->begin();
       gluino_to_sdown_LR2->get_products()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->get_name()) << " ";
    }
  std::cout
  << std::endl
  << "has direct decay? " << gluino_to_sdown_LR2->has_direct_decay();
  std::cout
  << std::endl
  << "direct_decay_checker: gluino -> sup_L, sup_R"
  << std::endl << "pointer = "
  << gluino_to_sup_LR
  << std::endl << *(gluino_to_sup_LR->get_decayer()->get_name())
  << " -> ";
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = gluino_to_sup_LR->get_products()->begin();
       gluino_to_sup_LR->get_products()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->get_name()) << " ";
    }
  std::cout
  << std::endl
  << "has direct decay? " << gluino_to_sup_LR->has_direct_decay();

  LHC_FASER::exclusive_BR_calculator*
  gluino_to_sdown_L_BR
  = test_input_handler.get_exclusive_BR_calculator(
                                               test_input_handler.get_gluino(),
                                              test_input_handler.get_sdown_L(),
                                                    true,
                                       test_input_handler.get_not_in_jets5() );
  LHC_FASER::exclusive_BR_calculator*
  gluino_to_sdown_R_BR
  = test_input_handler.get_exclusive_BR_calculator(
                                               test_input_handler.get_gluino(),
                                              test_input_handler.get_sdown_R(),
                                                    true,
                                       test_input_handler.get_not_in_jets5() );
  LHC_FASER::exclusive_BR_calculator*
  gluino_to_sdown_L_BR2
  = test_input_handler.get_exclusive_BR_calculator(
                                               test_input_handler.get_gluino(),
                                              test_input_handler.get_sdown_L(),
                                                    true,
                                       test_input_handler.get_not_in_jets5() );

  std::cout
  << std::endl
  << "exclusive BR checker: gluino -> sdown_L"
  << std::endl << "pointer = "
  << gluino_to_sdown_L_BR
  << std::endl << *(gluino_to_sdown_L_BR->get_decayer()->get_name())
  << " -> " << *(gluino_to_sdown_L_BR->get_product()->get_name())
  << " has BR " << gluino_to_sdown_L_BR->get_BR()
  << std::endl
  << "2nd gluino -> sdown_L"
  << std::endl << "pointer = "
  << gluino_to_sdown_L_BR2
  << std::endl << *(gluino_to_sdown_L_BR2->get_decayer()->get_name())
  << " -> " << *(gluino_to_sdown_L_BR2->get_product()->get_name())
  << " has BR " << gluino_to_sdown_L_BR2->get_BR()
  << std::endl
  << "gluino -> sdown_R"
  << std::endl << "pointer = "
  << gluino_to_sdown_R_BR
  << std::endl << *(gluino_to_sdown_R_BR->get_decayer()->get_name())
  << " -> " << *(gluino_to_sdown_R_BR->get_product()->get_name())
  << " has BR " << gluino_to_sdown_R_BR->get_BR()
  << std::endl
  << std::endl;


  SLHA_data.read_file( "SPS2_spectrum.out" );
  test_readier.ready_for_new_point();
  std::cout
  << std::endl
  << "SPS2 MASS:"
  << test_input_handler.get_CppSLHA()->get_MASS_as_string()
  << std::endl
  << std::endl
  << *(test_input_handler.get_gluino()->get_name()) << " "
  << test_input_handler.get_gluino()->get_mass() << " GeV"
  << std::endl
  << std::endl
  << "charged sleptons: " << std::endl;
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = test_input_handler.get_charged_sleptons()->begin();
       test_input_handler.get_charged_sleptons()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->get_name()) << std::endl;
    }
  std::cout
  << *(gluino_cascade->get_decayer()->get_name())
  << ":" << std::endl
  << "pointer = " << gluino_cascade << ", pointer2 = " << gluino_cascade2
  << std::endl;
  for( std::vector< LHC_FASER::colored_cascade* >::const_iterator
       set_iterator = gluino_cascade->get_colored_cascades()->begin();
       gluino_cascade->get_colored_cascades()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout
      << *((*set_iterator)->get_decayer()->get_name())
      << " -> " << *((*set_iterator)->get_product()->get_name())
      << " ( " << (*set_iterator)->get_product_position() << " ) -> "
      << *((*set_iterator)->get_EWino()->get_name())
      << " ( " << (*set_iterator)->get_EWino_position() << " )";
      if( (*set_iterator)->is_sx() )
        {
          std::cout<< " => sx" << std::endl;
        }
      if( (*set_iterator)->is_gx() )
        {
          std::cout<< " => gx" << std::endl;
        }
      if( (*set_iterator)->is_sgx() )
        {
          std::cout<< " => sgx" << std::endl;
        }
      if( (*set_iterator)->is_gsx() )
        {
          std::cout<< " => gsx" << std::endl;
        }
      if( (*set_iterator)->is_sgsx() )
        {
          std::cout<< " => sgsx" << std::endl;
        }
    }
  std::cout
  << *(sbottom_cascade->get_decayer()->get_name())
  << ":" << std::endl
  << "pointer = " << sbottom_cascade
  << std::endl;
  for( std::vector< LHC_FASER::colored_cascade* >::const_iterator
       set_iterator = sbottom_cascade->get_colored_cascades()->begin();
       sbottom_cascade->get_colored_cascades()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout
      << *((*set_iterator)->get_decayer()->get_name())
      << " -> " << *((*set_iterator)->get_product()->get_name())
      << " ( " << (*set_iterator)->get_product_position() << " ) -> "
      << *((*set_iterator)->get_EWino()->get_name())
      << " ( " << (*set_iterator)->get_EWino_position() << " )";
      if( (*set_iterator)->is_sx() )
        {
          std::cout<< " => sx" << std::endl;
        }
      if( (*set_iterator)->is_gx() )
        {
          std::cout<< " => gx" << std::endl;
        }
      if( (*set_iterator)->is_sgx() )
        {
          std::cout<< " => sgx" << std::endl;
        }
      if( (*set_iterator)->is_gsx() )
        {
          std::cout<< " => gsx" << std::endl;
        }
      if( (*set_iterator)->is_sgsx() )
        {
          std::cout<< " => sgsx" << std::endl;
        }
    }
  std::cout
  << std::endl
  << "direct_decay_checker: gluino -> sdown_L, sdown_R"
  << std::endl << "pointer = "
  << gluino_to_sdown_LR
  << std::endl << *(gluino_to_sdown_LR->get_decayer()->get_name())
  << " -> ";
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = gluino_to_sdown_LR->get_products()->begin();
       gluino_to_sdown_LR->get_products()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->get_name()) << " ";
    }
  std::cout
  << std::endl
  << "has direct decay? " << gluino_to_sdown_LR->has_direct_decay();
  std::cout
  << std::endl
  << "exclusive BR checker: gluino -> sdown_L"
  << std::endl << "pointer = "
  << gluino_to_sdown_L_BR
  << std::endl << *(gluino_to_sdown_L_BR->get_decayer()->get_name())
  << " -> " << *(gluino_to_sdown_L_BR->get_product()->get_name())
  << " has BR " << gluino_to_sdown_L_BR->get_BR()
  << std::endl;


  SLHA_data.read_file( "SPS1a_spectrum.out" );
  test_readier.ready_for_new_point();
  std::cout
  << std::endl
  << "SPS1a MASS:"
  << test_input_handler.get_CppSLHA()->get_MASS_as_string()
  << std::endl
  << std::endl
  << *(test_input_handler.get_gluino()->get_name()) << " "
  << test_input_handler.get_gluino()->get_mass() << " GeV"
  << std::endl
  << std::endl
  << "charged sleptons: " << std::endl;
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = test_input_handler.get_charged_sleptons()->begin();
       test_input_handler.get_charged_sleptons()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->get_name()) << std::endl;
    }
  std::cout
  << *(gluino_cascade->get_decayer()->get_name())
  << ":" << std::endl
  << "pointer = " << gluino_cascade << ", pointer2 = " << gluino_cascade2
  << std::endl;
  for( std::vector< LHC_FASER::colored_cascade* >::const_iterator
       set_iterator = gluino_cascade->get_colored_cascades()->begin();
       gluino_cascade->get_colored_cascades()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout
      << *((*set_iterator)->get_decayer()->get_name())
      << " -> " << *((*set_iterator)->get_product()->get_name())
      << " ( " << (*set_iterator)->get_product_position() << " ) -> "
      << *((*set_iterator)->get_EWino()->get_name())
      << " ( " << (*set_iterator)->get_EWino_position() << " )";
      if( (*set_iterator)->is_sx() )
        {
          std::cout<< " => sx" << std::endl;
        }
      if( (*set_iterator)->is_gx() )
        {
          std::cout<< " => gx" << std::endl;
        }
      if( (*set_iterator)->is_sgx() )
        {
          std::cout<< " => sgx" << std::endl;
        }
      if( (*set_iterator)->is_gsx() )
        {
          std::cout<< " => gsx" << std::endl;
        }
      if( (*set_iterator)->is_sgsx() )
        {
          std::cout<< " => sgsx" << std::endl;
        }
    }
  std::cout
  << *(sbottom_cascade->get_decayer()->get_name())
  << ":" << std::endl
  << "pointer = " << sbottom_cascade
  << std::endl;
  for( std::vector< LHC_FASER::colored_cascade* >::const_iterator
       set_iterator = sbottom_cascade->get_colored_cascades()->begin();
       sbottom_cascade->get_colored_cascades()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout
      << *((*set_iterator)->get_decayer()->get_name())
      << " -> " << *((*set_iterator)->get_product()->get_name())
      << " ( " << (*set_iterator)->get_product_position() << " ) -> "
      << *((*set_iterator)->get_EWino()->get_name())
      << " ( " << (*set_iterator)->get_EWino_position() << " )";
      if( (*set_iterator)->is_sx() )
        {
          std::cout<< " => sx" << std::endl;
        }
      if( (*set_iterator)->is_gx() )
        {
          std::cout<< " => gx" << std::endl;
        }
      if( (*set_iterator)->is_sgx() )
        {
          std::cout<< " => sgx" << std::endl;
        }
      if( (*set_iterator)->is_gsx() )
        {
          std::cout<< " => gsx" << std::endl;
        }
      if( (*set_iterator)->is_sgsx() )
        {
          std::cout<< " => sgsx" << std::endl;
        }
    }
  std::cout
  << std::endl
  << "direct_decay_checker: gluino -> sdown_L, sdown_R"
  << std::endl << "pointer = "
  << gluino_to_sdown_LR
  << std::endl << *(gluino_to_sdown_LR->get_decayer()->get_name())
  << " -> ";
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = gluino_to_sdown_LR->get_products()->begin();
       gluino_to_sdown_LR->get_products()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->get_name()) << " ";
    }
  std::cout
  << std::endl
  << "has direct decay? " << gluino_to_sdown_LR->has_direct_decay();
  std::cout
  << std::endl
  << "exclusive BR checker: gluino -> sdown_L"
  << std::endl << "pointer = "
  << gluino_to_sdown_L_BR
  << std::endl << *(gluino_to_sdown_L_BR->get_decayer()->get_name())
  << " -> " << *(gluino_to_sdown_L_BR->get_product()->get_name())
  << " has BR " << gluino_to_sdown_L_BR->get_BR()
  << std::endl;

  SLHA_data.read_file( "SPS2_spectrum.out" );
  test_readier.ready_for_new_point();
  std::cout
  << std::endl
  << "exclusive BR checker: gluino -> sdown_L"
  << std::endl << "pointer = "
  << gluino_to_sdown_L_BR
  << std::endl << *(gluino_to_sdown_L_BR->get_decayer()->get_name())
  << " -> " << *(gluino_to_sdown_L_BR->get_product()->get_name())
  << " has BR " << gluino_to_sdown_L_BR->get_BR()
  << std::endl;
  std::cout
  << std::endl
  << "exclusive BR checker: gluino -> sdown_L"
  << std::endl << "pointer = "
  << gluino_to_sdown_L_BR
  << std::endl << *(gluino_to_sdown_L_BR->get_decayer()->get_name())
  << " -> " << *(gluino_to_sdown_L_BR->get_product()->get_name())
  << " has BR " << gluino_to_sdown_L_BR->get_BR()
  << std::endl;
  SLHA_data.read_file( "SPS1a_spectrum.out" );
  test_readier.ready_for_new_point();
  std::cout
  << std::endl
  << "exclusive BR checker: gluino -> sdown_L"
  << std::endl << "pointer = "
  << gluino_to_sdown_L_BR
  << std::endl << *(gluino_to_sdown_L_BR->get_decayer()->get_name())
  << " -> " << *(gluino_to_sdown_L_BR->get_product()->get_name())
  << " has BR " << gluino_to_sdown_L_BR->get_BR()
  << std::endl;
  std::cout
  << std::endl
  << "exclusive BR checker: gluino -> sdown_L"
  << std::endl << "pointer = "
  << gluino_to_sdown_L_BR
  << std::endl << *(gluino_to_sdown_L_BR->get_decayer()->get_name())
  << " -> " << *(gluino_to_sdown_L_BR->get_product()->get_name())
  << " has BR " << gluino_to_sdown_L_BR->get_BR()
  << std::endl
  << std::endl
  << std::endl
  << "gluino_to_sdown_LR = " << gluino_to_sdown_LR
  << std::endl
  << "gluino pointer = " << test_input_handler.get_gluino()
  << std::endl
  << "sdown_L pointer = " << test_input_handler.get_sdown_L()
  << std::endl
  << "sdown_R pointer = " << test_input_handler.get_sdown_R()
  << std::endl;

  LHC_FASER::readier_for_new_point*
  test_readier2 = new LHC_FASER::readier_for_new_point();
  LHC_FASER::exclusive_BR_calculator*
  sdown_L_chargino_one_minus
  = new LHC_FASER::exclusive_BR_calculator( test_input_handler.get_sdown_L(),
                                         test_input_handler.get_chargino_one(),
                                            false,
                                         test_input_handler.get_not_in_jets4(),
                                            test_readier2 );
  LHC_FASER::exclusive_BR_calculator*
  sdown_L_chargino_one_plus
  = new LHC_FASER::exclusive_BR_calculator( test_input_handler.get_sdown_L(),
                                         test_input_handler.get_chargino_one(),
                                            true,
                                         test_input_handler.get_not_in_jets4(),
                                            test_readier2 );
  std::cout
  << std::endl << "sdown_L_chargino_one_minus BR = "
  << sdown_L_chargino_one_minus->get_BR()
  << std::endl << "sdown_L_chargino_one_plus BR = "
  << sdown_L_chargino_one_plus->get_BR();
  test_readier2->ready_for_new_point();
  delete sdown_L_chargino_one_plus;
  std::cout
  << std::endl << "sdown_L_chargino_one_minus BR = "
  << sdown_L_chargino_one_minus->get_BR();
  test_readier2->ready_for_new_point();
  delete test_readier2;
  std::cout
  << std::endl << "sdown_L_chargino_one_minus BR = "
  << sdown_L_chargino_one_minus->get_BR();
  delete sdown_L_chargino_one_minus;
  **/


  /* testing the cross-section grid stuff: (should work, may need slight changes of input)*//**

  std::cout
  << std::endl
  << "path_to_grids = "
  << *(test_input_handler.get_path_to_grids())
  << std::endl;

  std::string
  path_to_test_grid( *(test_input_handler.get_path_to_cross_section_grids()) );
  path_to_test_grid.append( "/7TeV/gluino+gluino_LO.dat" );
  LHC_FASER::square_grid test_grid( &path_to_test_grid,
                                    &bluh,
                                    NULL );
  LHC_FASER::cross_section_table gg_table( &test_grid,
           test_input_handler.get_sparticle_production_combinations()->at( 0 ),
                                           1.0,
                                           &test_input_handler );
  std::cout
  << std::endl << "x: " << test_grid.get_lowest_x() << " - "
  << test_grid.get_highest_x();
  std::cout
  << std::endl << "y: " << test_grid.get_lowest_y() << " - "
  << test_grid.get_highest_y();
  std::cout
  << std::endl << "value at { 200.0, 200.0 } = " << test_grid.value_at( 200.0,
                                                                       200.0 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at { 200.0, 205.0 } = " << test_grid.value_at( 200.0,
                                                                       205.0 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at { 205.0, 200.0 } = " << test_grid.value_at( 205.0,
                                                                       200.0 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at { 201.0, 201.0 } = " << test_grid.value_at( 201.0,
                                                                       201.0 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at { 210.0, 210.0 } = " << test_grid.value_at( 210.0,
                                                                       210.0 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at { 199.5, 201.0 } = " << test_grid.value_at( 199.5,
                                                                       201.0 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at { 199.5, 199.5 } = " << test_grid.value_at( 199.5,
                                                                       199.5 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at { 2000.0, 200.5 } = "
  << test_grid.value_at( 2000.0,
                         200.5 );
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout
  << std::endl << "value at { 2000.5, 200.5 } = "
  << test_grid.value_at( 2000.5,
                         200.5 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at { 200.5, 2000.0 } = "
  << test_grid.value_at( 200.5,
                         2000.0 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at { 200.5, 2000.5 } = "
  << test_grid.value_at( 200.5,
                         2000.5 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at { 1999.5, 1999.5 } = "
  << test_grid.value_at( 1999.5,
                         1999.5 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at { 2000.0, 2000.0 } = "
  << test_grid.value_at( 2000.0,
                         2000.0 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at { 2000.5, 2000.5 } = "
  << test_grid.value_at( 2000.5,
                         2000.5 );
  std::cout << std::endl;
  std::cout
  << std::endl << "ggtable value = "
  << gg_table.get_value();

  path_to_test_grid.assign(
                     *(test_input_handler.get_path_to_cross_section_grids()) );
  path_to_test_grid.append( "/7TeV/gluino+gluino_K.dat" );
  LHC_FASER::square_grid test_K( &path_to_test_grid,
                                 &bluh,
                                 NULL );
  std::cout
  << std::endl << "K value at { 200.5, 200.5 } = "
  << test_K.value_at( 200.5,
                      200.5 );
  std::cout << std::endl;


  LHC_FASER::cross_section_handler cross_section_handler_tester(
                                                         &test_input_handler );
  LHC_FASER::cross_section_table_set*
  test_cross_section_set = cross_section_handler_tester.get_table_set( 7 );
  std::cout
  << std::endl << "test_cross_section_set = " << test_cross_section_set;
  N.B.: for the next bit to work, cross_section_handler::get_table() has to be
  moved back to public.
  LHC_FASER::cross_section_table*
  sd_go_sigma
  = test_cross_section_set->get_table(
         test_input_handler.get_sparticle_production_combinations()->at( 1 ) );
  LHC_FASER::cross_section_table*
  su_go_sigma
  = test_cross_section_set->get_table(
        test_input_handler.get_sparticle_production_combinations()->at( 10 ) );
  LHC_FASER::cross_section_table_set*
  test_cross_section_set2 = cross_section_handler_tester.get_table_set( 7 );
  std::cout
  << std::endl << "test_cross_section_set2 = " << test_cross_section_set2;
  LHC_FASER::cross_section_table*
  sd_go_sigma2
  = test_cross_section_set2->get_table(
         test_input_handler.get_sparticle_production_combinations()->at( 1 ) );
  std::cout
  << std::endl << "no extra square_grid loading?";

  path_to_test_grid.assign(
                     *(test_input_handler.get_path_to_cross_section_grids()) );
  path_to_test_grid.append( "/7TeV/sdown_L+gluino_LO.dat" );
  LHC_FASER::square_grid test_sd_go_LO( &path_to_test_grid,
                                        &bluh,
                                        NULL );
  path_to_test_grid.assign(
                     *(test_input_handler.get_path_to_cross_section_grids()) );
  path_to_test_grid.append( "/7TeV/sdown_L+gluino_K.dat" );
  LHC_FASER::square_grid test_sd_go_K( &path_to_test_grid,
                                       &bluh,
                                       NULL );
  double
  sd_mass
  = test_input_handler.get_sparticle_production_combinations()->at( 1
                                   )->get_first_pointer()->get_absolute_mass();
  double
  go_mass
  = test_input_handler.get_sparticle_production_combinations()->at( 1
                                  )->get_second_pointer()->get_absolute_mass();
  std::cout << std::endl;
  std::cout
  << std::endl << "sd_mass = " << sd_mass
  << std::endl << "go_mass = " << go_mass
  << std::endl << "test_sd_go_LO value = " << test_sd_go_LO.value_at( sd_mass,
                                                                      go_mass )
  << std::endl << "test_sd_go_K value = " << test_sd_go_K.value_at( sd_mass,
                                                                    go_mass )
  << std::endl << "hence NLO value should be "
  << ( test_sd_go_LO.value_at( sd_mass,
                               go_mass )
       * test_sd_go_K.value_at( sd_mass,
                                go_mass ) )
  << std::endl << "sd_go_sigma = " << sd_go_sigma
  << std::endl << "sd_go_sigma value = "
  << sd_go_sigma->get_value()
  << std::endl << "sd_go_sigma2 = " << sd_go_sigma2
  << std::endl << "sd_go_sigma2 value = "
  << sd_go_sigma->get_value()
  << std::endl << "su_go_sigma = " << su_go_sigma
  << std::endl << "su_go_sigma value = "
  << su_go_sigma->get_value();

  **/


  /* testing the kinematics grid stuff: (ahhhh, tssss, might still work involes a lot of old classes that are not yet deleted)*//**

  std::string
  acceptance_grid_path( *(test_input_handler.get_path_to_kinematics_grids()) );
  std::cout
  << std::endl << "path to kinematics grids = " << acceptance_grid_path;
  std::cout << std::endl;
  acceptance_grid_path.append(
                        "/7TeV/jets/Atlas4jMET/squark+gluino_acceptance.dat" );
  LHC_FASER::acceptance_grid*
  jets_test_acceptance
  = new LHC_FASER::acceptance_grid( &acceptance_grid_path );
  std::cout
  << std::endl << "jets_test_acceptance->value_at( "
  << jets_test_acceptance->get_lowest_squark_mass() << ", "
  << jets_test_acceptance->get_lowest_gluino_mass() << ", 0.0, 0.0, 0 ) = "
  << jets_test_acceptance->value_at(
                                jets_test_acceptance->get_lowest_squark_mass(),
                                jets_test_acceptance->get_lowest_gluino_mass(),
                                0.0,
                                0.0,
                                0,
                                false )
  << std::endl
  << "jets_test_acceptance->value_at( 123.4, 567.8, 12.3, 45.6, 4 ) = "
  << jets_test_acceptance->value_at( 123.4,
                                     567.8,
                                     12.3,
                                     45.6,
                                     4,
                                     false )
  << std::endl
  << "jets_test_acceptance->value_at( 234.5, 678.9, 12.3, 45.6, 4 ) = "
  << jets_test_acceptance->value_at( 234.5,
                                     678.9,
                                     12.3,
                                     45.6,
                                     4,
                                     false )
  << std::endl
  << "jets_test_acceptance->value_at( 876.5, 432.1, 1.0, 1.0, 4 ) = "
  << jets_test_acceptance->value_at( 876.5,
                                     432.1,
                                     1.0,
                                     1.0,
                                     4,
                                     false )
  << std::endl
  << "jets_test_acceptance->value_at( 501.0, 501.0, 0.0, 0.0, 4 ) = "
  << jets_test_acceptance->value_at( 501.0,
                                     501.0,
                                     0.0,
                                     0.0,
                                     4,
                                     false );
  std::cout << std::endl;
  acceptance_grid_path.assign(
                        *(test_input_handler.get_path_to_kinematics_grids()) );
  std::cout
  << std::endl << "path to kinematics grids = " << acceptance_grid_path;
  std::cout << std::endl;
  acceptance_grid_path.append( "/7TeV/leptons/squark+gluino_acceptance.dat" );
  LHC_FASER::acceptance_grid*
  leptons_test_acceptance
  = new LHC_FASER::acceptance_grid( &acceptance_grid_path );
  std::cout
  << std::endl << "leptons_test_acceptance->value_at( "
  << leptons_test_acceptance->get_lowest_squark_mass() << ", "
  << leptons_test_acceptance->get_lowest_gluino_mass() << ", 0.0, 0.0, 0 ) = "
  << leptons_test_acceptance->value_at(
                             leptons_test_acceptance->get_lowest_squark_mass(),
                             leptons_test_acceptance->get_lowest_gluino_mass(),
                                        0.0,
                                        0.0,
                                        0,
                                        true )
  << std::endl
  << "leptons_test_acceptance->value_at( 123.4, 567.8, 12.3, 45.6, 4 ) = "
  << leptons_test_acceptance->value_at( 123.4,
                                        567.8,
                                        12.3,
                                        45.6,
                                        4,
                                        false )
  << std::endl
  << "leptons_test_acceptance->value_at( 234.5, 678.9, 12.3, 45.6, 4 ) = "
  << leptons_test_acceptance->value_at( 234.5,
                                        678.9,
                                        12.3,
                                        45.6,
                                        4,
                                        false )
  << std::endl
  << "leptons_test_acceptance->value_at( 876.5, 432.1, 1.0, 1.0, 4 ) = "
  << leptons_test_acceptance->value_at( 876.5,
                                        432.1,
                                        1.0,
                                        1.0,
                                        4,
                                        false )
  << std::endl
  << "leptons_test_acceptance->value_at( 501.0, 501.0, 0.0, 0.0, 4 ) = "
  << leptons_test_acceptance->value_at( 501.0,
                                        501.0,
                                        0.0,
                                        0.0,
                                        4,
                                        false );
  std::cout << std::endl;

  delete leptons_test_acceptance;

  LHC_FASER::colored_cascade_set*
  gluino_cascade3
  = test_input_handler.get_colored_cascade_set(
                                             test_input_handler.get_gluino() );
  LHC_FASER::colored_cascade_set*
  sdown_L_cascade
  = test_input_handler.get_colored_cascade_set(
                                            test_input_handler.get_sdown_L() );
  int const sdL_go = 1;
  acceptance_grid_path.assign(
                        *(test_input_handler.get_path_to_kinematics_grids()) );
  acceptance_grid_path.append( "/7TeV/jets/Atlas4jMET" );
  std::cout
  << std::endl << "path to jets kinematics grids = " << acceptance_grid_path;
  std::cout << std::endl;
  LHC_FASER::jet_acceptance_table test_jet_kinematics_table(
                                                         &acceptance_grid_path,
                                                             &bluh,
                                                         &test_input_handler );
  std::cout
  << std::endl << "test_jet_kinematics_table.get_name() = "
  << *(test_jet_kinematics_table.get_name())
  << std::endl << "test_jet_kinematics_table.get_acceptance( 4,"
  << std::endl
  << "[ ";
  if( test_input_handler.get_sparticle_production_combinations()->at( sdL_go
                                                       )->is_first_particle() )
    {
      std::cout
      << *(test_input_handler.get_sparticle_production_combinations(
                             )->at( sdL_go )->get_first_pointer()->get_name());
    }
  else
    {
      std::cout
      << *(test_input_handler.get_sparticle_production_combinations()->at(
                                sdL_go )->get_first_pointer()->get_antiname());
    }
  std::cout << " + ";
  if( test_input_handler.get_sparticle_production_combinations()->at( sdL_go
                                                      )->is_second_particle() )
    {
      std::cout
      << *(test_input_handler.get_sparticle_production_combinations()->at(
                                   sdL_go )->get_second_pointer()->get_name());
    }
  else
    {
      std::cout
      << *(test_input_handler.get_sparticle_production_combinations()->at(
                               sdL_go )->get_second_pointer()->get_antiname());
    }
  std::cout << " ]," << std::endl;
  std::cout
  << "[ "
  << *(gluino_cascade3->get_colored_cascades()->at( 1
                                                  )->get_decayer()->get_name())
  << " -> " << *(gluino_cascade3->get_colored_cascades()->at( 1
                                                  )->get_product()->get_name())
  << " ( " << gluino_cascade3->get_colored_cascades()->at( 1
                                          )->get_product_position() << " ) -> "
  << *(gluino_cascade3->get_colored_cascades()->at( 1
                                                    )->get_EWino()->get_name())
  << " ( " << gluino_cascade3->get_colored_cascades()->at( 1
                                               )->get_EWino_position() << " )";
  if( gluino_cascade3->get_colored_cascades()->at( 1 )->is_sx() )
    {
      std::cout<< " => sx";
    }
  if( gluino_cascade3->get_colored_cascades()->at( 1 )->is_gx() )
    {
      std::cout<< " => gx";
    }
  if( gluino_cascade3->get_colored_cascades()->at( 1 )->is_sgx() )
    {
      std::cout<< " => sgx";
    }
  if( gluino_cascade3->get_colored_cascades()->at( 1 )->is_gsx() )
    {
      std::cout<< " => gsx";
    }
  if( gluino_cascade3->get_colored_cascades()->at( 1 )->is_sgsx() )
    {
      std::cout<< " => sgsx";
    }
  std::cout
  << " ],"
  << std::endl;
  std::cout
  << "[ "
  << *(sdown_L_cascade->get_colored_cascades()->at( 1
                                                  )->get_decayer()->get_name())
  << " -> " << *(sdown_L_cascade->get_colored_cascades()->at( 1
                                                  )->get_product()->get_name())
  << " ( " << sdown_L_cascade->get_colored_cascades()->at( 1
                                          )->get_product_position() << " ) -> "
  << *(sdown_L_cascade->get_colored_cascades()->at( 1
                                                    )->get_EWino()->get_name())
  << " ( " << sdown_L_cascade->get_colored_cascades()->at( 1
                                               )->get_EWino_position() << " )";
  if( sdown_L_cascade->get_colored_cascades()->at( 1 )->is_sx() )
    {
      std::cout<< " => sx";
    }
  if( sdown_L_cascade->get_colored_cascades()->at( 1 )->is_gx() )
    {
      std::cout<< " => gx";
    }
  if( sdown_L_cascade->get_colored_cascades()->at( 1 )->is_sgx() )
    {
      std::cout<< " => sgx";
    }
  if( sdown_L_cascade->get_colored_cascades()->at( 1 )->is_gsx() )
    {
      std::cout<< " => gsx";
    }
  if( sdown_L_cascade->get_colored_cascades()->at( 1 )->is_sgsx() )
    {
      std::cout<< " => sgsx";
    }
  std::cout
  << " ] )->get_acceptance()"
  << std::endl << " = "
  << test_jet_kinematics_table.get_acceptance( 4,
      test_input_handler.get_sparticle_production_combinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
          sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance();
  std::cout << std::endl;

  acceptance_grid_path.assign(
                        *(test_input_handler.get_path_to_kinematics_grids()) );
  acceptance_grid_path.append( "/7TeV/leptons" );
  std::cout
  << std::endl
  << "path to leptons kinematics grids = " << acceptance_grid_path;
  std::cout << std::endl;
  LHC_FASER::lepton_acceptance_table test_lepton_kinematics_table(
                                                         &acceptance_grid_path,
                                                                   7,
                                                         &test_input_handler );
  std::cout
  << std::endl << "test_lepton_kinematics_table.get_energy() = "
  << test_lepton_kinematics_table.get_energy()
  << std::endl << "test_lepton_kinematics_table.get_acceptance( "
  << "[ ";
  if( test_input_handler.get_sparticle_production_combinations()->at( sdL_go
                                                       )->is_first_particle() )
    {
      std::cout
      << *(test_input_handler.get_sparticle_production_combinations()->at(
                                    sdL_go )->get_first_pointer()->get_name());
    }
  else
    {
      std::cout
      << *(test_input_handler.get_sparticle_production_combinations()->at(
                                sdL_go )->get_first_pointer()->get_antiname());
    }
  std::cout << " + ";
  if( test_input_handler.get_sparticle_production_combinations()->at( sdL_go
                                                      )->is_second_particle() )
    {
      std::cout
      << *(test_input_handler.get_sparticle_production_combinations()->at(
                                   sdL_go )->get_second_pointer()->get_name());
    }
  else
    {
      std::cout
      << *(test_input_handler.get_sparticle_production_combinations()->at(
                               sdL_go )->get_second_pointer()->get_antiname());
    }
  std::cout << " ]," << std::endl;
  std::cout
  << "[ "
  << *(gluino_cascade3->get_colored_cascades()->at( 1
                                                  )->get_decayer()->get_name())
  << " -> " << *(gluino_cascade3->get_colored_cascades()->at( 1
                                                  )->get_product()->get_name())
  << " ( " << gluino_cascade3->get_colored_cascades()->at( 1
                                          )->get_product_position() << " ) -> "
  << *(gluino_cascade3->get_colored_cascades()->at( 1
                                                    )->get_EWino()->get_name())
  << " ( " << gluino_cascade3->get_colored_cascades()->at( 1
                                               )->get_EWino_position() << " )";
  if( gluino_cascade3->get_colored_cascades()->at( 1 )->is_sx() )
    {
      std::cout<< " => sx";
    }
  if( gluino_cascade3->get_colored_cascades()->at( 1 )->is_gx() )
    {
      std::cout<< " => gx";
    }
  if( gluino_cascade3->get_colored_cascades()->at( 1 )->is_sgx() )
    {
      std::cout<< " => sgx";
    }
  if( gluino_cascade3->get_colored_cascades()->at( 1 )->is_gsx() )
    {
      std::cout<< " => gsx";
    }
  if( gluino_cascade3->get_colored_cascades()->at( 1 )->is_sgsx() )
    {
      std::cout<< " => sgsx";
    }
  std::cout
  << " ],"
  << std::endl;
  std::cout
  << "[ "
  << *(sdown_L_cascade->get_colored_cascades()->at( 1
                                                  )->get_decayer()->get_name())
  << " -> " << *(sdown_L_cascade->get_colored_cascades()->at( 1
                                                  )->get_product()->get_name())
  << " ( " << sdown_L_cascade->get_colored_cascades()->at( 1
                                          )->get_product_position() << " ) -> "
  << *(sdown_L_cascade->get_colored_cascades()->at( 1
                                                    )->get_EWino()->get_name())
  << " ( " << sdown_L_cascade->get_colored_cascades()->at( 1
                                               )->get_EWino_position() << " )";
  if( sdown_L_cascade->get_colored_cascades()->at( 1 )->is_sx() )
    {
      std::cout<< " => sx";
    }
  if( sdown_L_cascade->get_colored_cascades()->at( 1 )->is_gx() )
    {
      std::cout<< " => gx";
    }
  if( sdown_L_cascade->get_colored_cascades()->at( 1 )->is_sgx() )
    {
      std::cout<< " => sgx";
    }
  if( sdown_L_cascade->get_colored_cascades()->at( 1 )->is_gsx() )
    {
      std::cout<< " => gsx";
    }
  if( sdown_L_cascade->get_colored_cascades()->at( 1 )->is_sgsx() )
    {
      std::cout<< " => sgsx";
    }
  std::cout
  << " ] )->get_acceptance( 23.4, 20.0 )"
  << std::endl << " = "
  << test_lepton_kinematics_table.get_acceptance(
      test_input_handler.get_sparticle_production_combinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                        )->acceptance_at( 23.4,
                                                                        20.0 );
  std::cout << std::endl;
  std::cout
  << "...->get_acceptance( -1.0, 10.0 )"
  << std::endl << " = "
  << test_lepton_kinematics_table.get_acceptance(
      test_input_handler.get_sparticle_production_combinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                        )->acceptance_at( -1.0,
                                                                        10.0 );
  std::cout << std::endl;
  std::cout
  << "...->get_acceptance( 1.0, 10.0 )"
  << std::endl << " = "
  << test_lepton_kinematics_table.get_acceptance(
      test_input_handler.get_sparticle_production_combinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                        )->acceptance_at( 1.0,
                                                                        10.0 );
  std::cout << std::endl;
  std::cout
  << "...->get_acceptance( 10.0, 10.0 )"
  << std::endl << " = "
  << test_lepton_kinematics_table.get_acceptance(
      test_input_handler.get_sparticle_production_combinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                        )->acceptance_at( 10.0,
                                                                        10.0 );
  std::cout << std::endl;
  std::cout
  << "...->get_acceptance( 12.3, 10.0 )"
  << std::endl << " = "
  << test_lepton_kinematics_table.get_acceptance(
      test_input_handler.get_sparticle_production_combinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                        )->acceptance_at( 12.3,
                                                                        10.0 );
  std::cout << std::endl;
  std::cout
  << "...->get_acceptance( 20.0, 10.0 )"
  << std::endl << " = "
  << test_lepton_kinematics_table.get_acceptance(
      test_input_handler.get_sparticle_production_combinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                        )->acceptance_at( 20.0,
                                                                        10.0 );
  std::cout << std::endl;
  std::cout
  << "...->get_acceptance( 79.9, 10.0 )"
  << std::endl << " = "
  << test_lepton_kinematics_table.get_acceptance(
      test_input_handler.get_sparticle_production_combinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                        )->acceptance_at( 79.9,
                                                                        10.0 );
  std::cout << std::endl;
  std::cout
  << "...->get_acceptance( 80.0, 10.0 )"
  << std::endl << " = "
  << test_lepton_kinematics_table.get_acceptance(
      test_input_handler.get_sparticle_production_combinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                        )->acceptance_at( 80.0,
                                                                        10.0 );
  std::cout << std::endl;
  std::cout
  << "...->get_acceptance( 100.0, 10.0 )"
  << std::endl << " = "
  << test_lepton_kinematics_table.get_acceptance(
      test_input_handler.get_sparticle_production_combinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                       )->acceptance_at( 100.0,
                                                                        10.0 );
  std::cout << std::endl;
  std::cout
  << std::endl << "test_lepton_kinematics_table.get_acceptance( "
  << "[ ";
  if( test_input_handler.get_sparticle_production_combinations()->at( sdL_go
                                                       )->is_first_particle() )
    {
      std::cout
      << *(test_input_handler.get_sparticle_production_combinations()->at(
                                    sdL_go )->get_first_pointer()->get_name());
    }
  else
    {
      std::cout
      << *(test_input_handler.get_sparticle_production_combinations()->at(
                                sdL_go )->get_first_pointer()->get_antiname());
    }
  std::cout << " + ";
  if( test_input_handler.get_sparticle_production_combinations()->at( sdL_go
                                                      )->is_second_particle() )
    {
      std::cout
      << *(test_input_handler.get_sparticle_production_combinations()->at(
                                   sdL_go )->get_second_pointer()->get_name());
    }
  else
    {
      std::cout
      << *(test_input_handler.get_sparticle_production_combinations()->at(
                               sdL_go )->get_second_pointer()->get_antiname());
    }
  std::cout << " ]," << std::endl;
  std::cout
  << "[ "
  << *(gluino_cascade3->get_colored_cascades()->back(
                                                  )->get_decayer()->get_name())
  << " -> " << *(gluino_cascade3->get_colored_cascades()->back(
                                                  )->get_product()->get_name())
  << " ( " << gluino_cascade3->get_colored_cascades()->back(
                                          )->get_product_position() << " ) -> "
  << *(gluino_cascade3->get_colored_cascades()->back(
                                                    )->get_EWino()->get_name())
  << " ( " << gluino_cascade3->get_colored_cascades()->back(
                                               )->get_EWino_position() << " )";
  if( gluino_cascade3->get_colored_cascades()->back()->is_sx() )
    {
      std::cout<< " => sx";
    }
  if( gluino_cascade3->get_colored_cascades()->back()->is_gx() )
    {
      std::cout<< " => gx";
    }
  if( gluino_cascade3->get_colored_cascades()->back()->is_sgx() )
    {
      std::cout<< " => sgx";
    }
  if( gluino_cascade3->get_colored_cascades()->back()->is_gsx() )
    {
      std::cout<< " => gsx";
    }
  if( gluino_cascade3->get_colored_cascades()->back()->is_sgsx() )
    {
      std::cout<< " => sgsx";
    }
  std::cout
  << " ],"
  << std::endl;
  std::cout
  << "[ "
  << *(sdown_L_cascade->get_colored_cascades()->at( 1
                                                  )->get_decayer()->get_name())
  << " -> " << *(sdown_L_cascade->get_colored_cascades()->at( 1
                                                  )->get_product()->get_name())
  << " ( " << sdown_L_cascade->get_colored_cascades()->at( 1
                                          )->get_product_position() << " ) -> "
  << *(sdown_L_cascade->get_colored_cascades()->at( 1
                                                    )->get_EWino()->get_name())
  << " ( " << sdown_L_cascade->get_colored_cascades()->at( 1
                                               )->get_EWino_position() << " )";
  if( sdown_L_cascade->get_colored_cascades()->at( 1 )->is_sx() )
    {
      std::cout<< " => sx";
    }
  if( sdown_L_cascade->get_colored_cascades()->at( 1 )->is_gx() )
    {
      std::cout<< " => gx";
    }
  if( sdown_L_cascade->get_colored_cascades()->at( 1 )->is_sgx() )
    {
      std::cout<< " => sgx";
    }
  if( sdown_L_cascade->get_colored_cascades()->at( 1 )->is_gsx() )
    {
      std::cout<< " => gsx";
    }
  if( sdown_L_cascade->get_colored_cascades()->at( 1 )->is_sgsx() )
    {
      std::cout<< " => sgsx";
    }
  std::cout
  << " ] )->get_acceptance( 23.4, 20.0 )"
  << std::endl << " = "
  << test_lepton_kinematics_table.get_acceptance(
      test_input_handler.get_sparticle_production_combinations()->at( sdL_go ),
                               gluino_cascade3->get_colored_cascades()->back(),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                        )->acceptance_at( 23.4,
                                                                        20.0 );

  LHC_FASER::kinematics_table
  kinematics_table_tester( &test_jet_kinematics_table,
                           4,
                           &test_lepton_kinematics_table );
  std::cout
  << std::endl << "kinematics_table_tester.get_column() = "
  << kinematics_table_tester.get_column()
  << std::endl << "kinematics_table_tester.get_jet_acceptance( "
  << "[ ";
  if( test_input_handler.get_sparticle_production_combinations()->at( sdL_go
                                                       )->is_first_particle() )
    {
      std::cout
      << *(test_input_handler.get_sparticle_production_combinations()->at(
                                    sdL_go )->get_first_pointer()->get_name());
    }
  else
    {
      std::cout
      << *(test_input_handler.get_sparticle_production_combinations()->at(
                                sdL_go )->get_first_pointer()->get_antiname());
    }
  std::cout << " + ";
  if( test_input_handler.get_sparticle_production_combinations()->at( sdL_go
                                                      )->is_second_particle() )
    {
      std::cout
      << *(test_input_handler.get_sparticle_production_combinations()->at(
                                   sdL_go )->get_second_pointer()->get_name());
    }
  else
    {
      std::cout
      << *(test_input_handler.get_sparticle_production_combinations()->at(
                               sdL_go )->get_second_pointer()->get_antiname());
    }
  std::cout << " ]," << std::endl;
  std::cout
  << "[ "
  << *(gluino_cascade3->get_colored_cascades()->at( 1
                                                  )->get_decayer()->get_name())
  << " -> " << *(gluino_cascade3->get_colored_cascades()->at( 1
                                                  )->get_product()->get_name())
  << " ( " << gluino_cascade3->get_colored_cascades()->at( 1
                                          )->get_product_position() << " ) -> "
  << *(gluino_cascade3->get_colored_cascades()->at( 1
                                                    )->get_EWino()->get_name())
  << " ( " << gluino_cascade3->get_colored_cascades()->at( 1
                                               )->get_EWino_position() << " )";
  if( gluino_cascade3->get_colored_cascades()->at( 1 )->is_sx() )
    {
      std::cout<< " => sx";
    }
  if( gluino_cascade3->get_colored_cascades()->at( 1 )->is_gx() )
    {
      std::cout<< " => gx";
    }
  if( gluino_cascade3->get_colored_cascades()->at( 1 )->is_sgx() )
    {
      std::cout<< " => sgx";
    }
  if( gluino_cascade3->get_colored_cascades()->at( 1 )->is_gsx() )
    {
      std::cout<< " => gsx";
    }
  if( gluino_cascade3->get_colored_cascades()->at( 1 )->is_sgsx() )
    {
      std::cout<< " => sgsx";
    }
  std::cout
  << " ],"
  << std::endl;
  std::cout
  << "[ "
  << *(sdown_L_cascade->get_colored_cascades()->at( 1
                                                  )->get_decayer()->get_name())
  << " -> " << *(sdown_L_cascade->get_colored_cascades()->at( 1
                                                  )->get_product()->get_name())
  << " ( " << sdown_L_cascade->get_colored_cascades()->at( 1
                                          )->get_product_position() << " ) -> "
  << *(sdown_L_cascade->get_colored_cascades()->at( 1
                                                    )->get_EWino()->get_name())
  << " ( " << sdown_L_cascade->get_colored_cascades()->at( 1
                                               )->get_EWino_position() << " )";
  if( sdown_L_cascade->get_colored_cascades()->at( 1 )->is_sx() )
    {
      std::cout<< " => sx";
    }
  if( sdown_L_cascade->get_colored_cascades()->at( 1 )->is_gx() )
    {
      std::cout<< " => gx";
    }
  if( sdown_L_cascade->get_colored_cascades()->at( 1 )->is_sgx() )
    {
      std::cout<< " => sgx";
    }
  if( sdown_L_cascade->get_colored_cascades()->at( 1 )->is_gsx() )
    {
      std::cout<< " => gsx";
    }
  if( sdown_L_cascade->get_colored_cascades()->at( 1 )->is_sgsx() )
    {
      std::cout<< " => sgsx";
    }
  std::cout
  << " ] )->get_acceptance()"
  << std::endl << " = "
  << kinematics_table_tester.get_jet_acceptance(
      test_input_handler.get_sparticle_production_combinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
          sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance();
  std::cout << std::endl;
  std::cout
  << "kinematics_table_tester.get_lepton_acceptance(...)->get_acceptance("
  << " 12.34, 8.0 )"
  << std::endl << " = "
  << kinematics_table_tester.get_lepton_acceptance(
      test_input_handler.get_sparticle_production_combinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                       )->acceptance_at( 12.34,
                                                                         8.0 );

  std::string
  kinematics_at_energy_name(
                        *(test_input_handler.get_path_to_kinematics_grids()) );
  kinematics_at_energy_name.append( "/7TeV" );
  std::string jet_signal_name( "Atlas4jMET" );
  LHC_FASER::kinematics_table_set
  kinematics_table_set_tester( 7,
                               &kinematics_at_energy_name,
                               &jet_signal_name,
                               &test_lepton_kinematics_table,
                               &test_input_handler );
  std::cout << std::endl;
  std::cout
  << "kinematics_table_set_tester.get_table( 4"
  << " )->get_jet_acceptance(...)->get_acceptance() = "
  << kinematics_table_set_tester.get_table( 4 )->get_jet_acceptance(
      test_input_handler.get_sparticle_production_combinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
          sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance();
  std::cout << std::endl;
  std::cout
  << "kinematics_table_set_tester.get_table( 4"
  << " )->get_lepton_acceptance(...)->get_acceptance( 12.34, 8.0 ) = "
  << kinematics_table_set_tester.get_table( 4 )->get_lepton_acceptance(
      test_input_handler.get_sparticle_production_combinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                       )->acceptance_at( 12.34,
                                                                         8.0 );
  std::cout << std::endl;

  LHC_FASER::kinematics_handler kinematics_handler_tester(
                                                         &test_input_handler );
  LHC_FASER::kinematics_table*
  kinematics_handler_table_one = kinematics_handler_tester.get_table( 7,
                                                              &jet_signal_name,
                                                                      4 );
  LHC_FASER::kinematics_table*
  kinematics_handler_table_two = kinematics_handler_tester.get_table( 7,
                                                              &jet_signal_name,
                                                                      5 );
  LHC_FASER::kinematics_table*
  kinematics_handler_table_three = kinematics_handler_tester.get_table( 7,
                                                              &jet_signal_name,
                                                                        4 );
  std::cout << std::endl;
  std::cout
  << std::endl
  << "kinematics_handler_table_one = " << kinematics_handler_table_one
  << std::endl
  << "->get_jet_acceptance(...)->get_acceptance() = "
  << kinematics_handler_table_one->get_jet_acceptance(
      test_input_handler.get_sparticle_production_combinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
           sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance()
  << std::endl
  << "->get_lepton_acceptance(...)->get_acceptance( 12.34, 8.0 ) = "
  << kinematics_handler_table_one->get_lepton_acceptance(
      test_input_handler.get_sparticle_production_combinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                       )->acceptance_at( 12.34,
                                                                         8.0 )
  << std::endl
  << "kinematics_handler_table_two = " << kinematics_handler_table_two
  << std::endl
  << "->get_jet_acceptance(...)->get_acceptance() = "
  << kinematics_handler_table_two->get_jet_acceptance(
      test_input_handler.get_sparticle_production_combinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
           sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance()
  << std::endl
  << "->get_lepton_acceptance(...)->get_acceptance( 12.34, 8.0 ) = "
  << kinematics_handler_table_two->get_lepton_acceptance(
      test_input_handler.get_sparticle_production_combinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                       )->acceptance_at( 12.34,
                                                                         8.0 )
  << std::endl
  << "kinematics_handler_table_three = " << kinematics_handler_table_three
  << std::endl
  << "->get_jet_acceptance(...)->get_acceptance() = "
  << kinematics_handler_table_three->get_jet_acceptance(
      test_input_handler.get_sparticle_production_combinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
           sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance()
  << std::endl
  << "->get_lepton_acceptance(...)->get_acceptance( 12.34, 8.0 ) = "
  << kinematics_handler_table_three->get_lepton_acceptance(
      test_input_handler.get_sparticle_production_combinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                       )->acceptance_at( 12.34,
                                                                         8.0 )
  << std::endl;
  **/


  /* testing the cascade calculation stuff: (doesnt work, why? objects have been redone completly, open: base_electroweak_cascade_stuff look at electro weak cascade class. find analogous way to set it up (n2->smu->n1 should be possible). look also at derived_electroweak...)*//**

  LHC_FASER::kinematics_handler test_kinematics( &test_input_handler );
  std::string test_signal_name( "Atlas4jMET" );
  LHC_FASER::lepton_acceptance_value*
  test_lepton_acceptance_value
  = test_kinematics.get_table( 7,
                               &test_signal_name,
                               4 )->get_lepton_acceptance(
           test_input_handler.get_sparticle_production_combinations()->at( 1 ),
                               test_input_handler.get_colored_cascade_set(
            test_input_handler.get_gluino() )->get_colored_cascades()->at( 1 ),
                               test_input_handler.get_colored_cascade_set(
         test_input_handler.get_sdown_L() )->get_colored_cascades()->at( 1 ) );

  LHC_FASER::cascade_acceptance_value two_jets_no_leptons( 0.0 );
  LHC_FASER::cascade_acceptance_value one_jet_one_muon( 0.0 );
  LHC_FASER::cascade_acceptance_value one_jet_one_electron( 0.0 );
  LHC_FASER::cascade_acceptance_value no_jets_two_muons( 0.0 );
  LHC_FASER::cascade_acceptance_value no_jets_two_electrons( 0.0 );
  LHC_FASER::cascade_acceptance_value no_jets_one_muon_one_electron( 0.0 );
  LHC_FASER::cascade_acceptance_value OSSF_minus_OSDF_leptons( 0.0 );
  LHC_FASER::cascade_acceptance_value one_jet_no_leptons( 0.0 );
  LHC_FASER::cascade_acceptance_value no_jets_one_muon( 0.0 );
  LHC_FASER::cascade_acceptance_value no_jets_one_electron( 0.0 );
  LHC_FASER::cascade_acceptance_value no_jets_no_leptons( 0.0 );


  // the following use n1 kinematics but n4 cascades:

  LHC_FASER::channel_calculator*
  smuon_channel = new LHC_FASER::charge_summed_neutralino_to_semu( 10.0,
                                                                   10.0,
                                                  test_lepton_acceptance_value,
                                              test_input_handler.get_sdown_L(),
                                                                   true,
                                      test_input_handler.get_neutralino_four(),
                                                                   true,
                                              test_input_handler.get_smuon_R(),
                                                           &test_input_handler,
                                   OSSF_minus_OSDF_leptons.get_value_pointer(),
                                         no_jets_two_muons.get_value_pointer(),
                                          no_jets_one_muon.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
  OSSF_minus_OSDF_leptons.add_channel( smuon_channel );
  no_jets_two_muons.add_channel( smuon_channel );
  no_jets_one_muon.add_channel( smuon_channel );
  no_jets_no_leptons.add_channel( smuon_channel );

  LHC_FASER::channel_calculator*
  stau_channel = new LHC_FASER::charge_summed_neutralino_to_stau( 10.0,
                                                                  10.0,
                                                                  40.0,
                                                  test_lepton_acceptance_value,
                                              test_input_handler.get_sdown_L(),
                                                                  true,
                                      test_input_handler.get_neutralino_four(),
                                                                  true,
                                             test_input_handler.get_stau_one(),
                                                           &test_input_handler,
                                       two_jets_no_leptons.get_value_pointer(),
                                          one_jet_one_muon.get_value_pointer(),
                                      one_jet_one_electron.get_value_pointer(),
                                         no_jets_two_muons.get_value_pointer(),
                                     no_jets_two_electrons.get_value_pointer(),
                             no_jets_one_muon_one_electron.get_value_pointer(),
                                        one_jet_no_leptons.get_value_pointer(),
                                          no_jets_one_muon.get_value_pointer(),
                                      no_jets_one_electron.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
  two_jets_no_leptons.add_channel( stau_channel );
  one_jet_one_muon.add_channel( stau_channel );
  one_jet_one_electron.add_channel( stau_channel );
  no_jets_two_muons.add_channel( stau_channel );
  no_jets_two_electrons.add_channel( stau_channel );
  no_jets_one_muon_one_electron.add_channel( stau_channel );
  one_jet_no_leptons.add_channel( stau_channel );
  no_jets_one_muon.add_channel( stau_channel );
  no_jets_one_electron.add_channel( stau_channel );
  no_jets_no_leptons.add_channel( stau_channel );

  LHC_FASER::channel_calculator*
  Z_channel = new LHC_FASER::charge_summed_neutralino_to_Z( 10.0,
                                                            10.0,
                                                            40.0,
                                                  test_lepton_acceptance_value,
                                              test_input_handler.get_sdown_L(),
                                                            true,
                                      test_input_handler.get_neutralino_four(),
                                                            true,
                                                    test_input_handler.get_Z(),
                                                           &test_input_handler,
                                       two_jets_no_leptons.get_value_pointer(),
                                          one_jet_one_muon.get_value_pointer(),
                                      one_jet_one_electron.get_value_pointer(),
                                   OSSF_minus_OSDF_leptons.get_value_pointer(),
                                         no_jets_two_muons.get_value_pointer(),
                                     no_jets_two_electrons.get_value_pointer(),
                             no_jets_one_muon_one_electron.get_value_pointer(),
                                        one_jet_no_leptons.get_value_pointer(),
                                          no_jets_one_muon.get_value_pointer(),
                                      no_jets_one_electron.get_value_pointer(),
                                      no_jets_no_leptons.get_value_pointer() );
  two_jets_no_leptons.add_channel( Z_channel );
  one_jet_one_muon.add_channel( Z_channel );
  one_jet_one_electron.add_channel( Z_channel );
  OSSF_minus_OSDF_leptons.add_channel( Z_channel );
  no_jets_two_muons.add_channel( Z_channel );
  no_jets_two_electrons.add_channel( Z_channel );
  no_jets_one_muon_one_electron.add_channel( Z_channel );
  one_jet_no_leptons.add_channel( Z_channel );
  no_jets_one_muon.add_channel( Z_channel );
  no_jets_one_electron.add_channel( Z_channel );
  no_jets_no_leptons.add_channel( Z_channel );

  double
  neutralino_four_to_smuon_R
  = test_input_handler.get_exclusive_BR_calculator(
                                      test_input_handler.get_neutralino_four(),
                                              test_input_handler.get_smuon_R(),
                                              true,
                               test_input_handler.get_empty_list() )->get_BR();
  double
  smuon_R_to_neutralino_one
  = test_input_handler.get_exclusive_BR_calculator(
                                              test_input_handler.get_smuon_R(),
                                       test_input_handler.get_neutralino_one(),
                                                    true,
                               test_input_handler.get_empty_list() )->get_BR();
  double
  neutralino_four_to_stau_one
  = test_input_handler.get_exclusive_BR_calculator(
                                      test_input_handler.get_neutralino_four(),
                                             test_input_handler.get_stau_one(),
                                                    true,
                               test_input_handler.get_empty_list() )->get_BR();
  double
  stau_one_to_neutralino_one
  = test_input_handler.get_exclusive_BR_calculator(
                                             test_input_handler.get_stau_one(),
                                       test_input_handler.get_neutralino_one(),
                                                    true,
                               test_input_handler.get_empty_list() )->get_BR();
  double
  neutralino_four_to_Z
  = test_input_handler.get_neutralino_four()->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                    test_input_handler.get_Z()->get_PDG_code(),
                                           CppSLHA::PDG_code::neutralino_one );
  std::cout
  << std::endl
  << "BR n4 -> smu_R = " << neutralino_four_to_smuon_R
  << std::endl
  << "BR smu_R -> n1 = " << smuon_R_to_neutralino_one
  << std::endl
  << "so BR n4 -> smu_R -> n1 + charge conjugate should be "
  << ( 2.0 * neutralino_four_to_smuon_R * smuon_R_to_neutralino_one )
  << std::endl
  << std::endl
  << "BR n4 -> stau_1 = " << neutralino_four_to_stau_one
  << std::endl
  << "BR stau_1 -> n1 = " << stau_one_to_neutralino_one
  << std::endl
  << "so BR n4 -> stau_1 -> n1 + charge conjugate should be "
  << ( 2.0 * neutralino_four_to_stau_one * stau_one_to_neutralino_one )
  << std::endl
  << "& combined with square of BR tauon -> muon, we have "
  << ( 2.0 * neutralino_four_to_stau_one * stau_one_to_neutralino_one
       * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR
       * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR )
  << std::endl
  << "while combined with square of BR tauon -> pion, we have "
       << ( 2.0 * neutralino_four_to_stau_one * stau_one_to_neutralino_one
            * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR
            * CppSLHA::PDG_data::tau_lepton_to_neutrino_hadron_BR )
  << std::endl
  << std::endl
  << "BR n4 -> Z = " << neutralino_four_to_Z
  << std::endl
  << "BR Z -> muon pair = " << CppSLHA::PDG_data::Z_to_muon_antimuon_BR
  << std::endl
  << "hence BR n4 -> Z -> muon pair = "
  << ( neutralino_four_to_Z * CppSLHA::PDG_data::Z_to_muon_antimuon_BR )
  << std::endl
  << "BR Z -> tauon pair = "
  << CppSLHA::PDG_data::Z_to_tau_lepton_antilepton_BR
  << std::endl
  << "combined with square of tauon -> muon, we have"
  << " n4 -> Z -> tauon pair -> muon pair = "
  << ( neutralino_four_to_Z * CppSLHA::PDG_data::Z_to_tau_lepton_antilepton_BR
       * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR
       * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR )
  << std::endl
  << "hence the full \"BR\" into 2 muons is "
  << ( 2.0 * ( neutralino_four_to_smuon_R * smuon_R_to_neutralino_one
               + neutralino_four_to_stau_one * stau_one_to_neutralino_one
                 * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR
                 * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR )
       + neutralino_four_to_Z
         * ( CppSLHA::PDG_data::Z_to_muon_antimuon_BR
             + CppSLHA::PDG_data::Z_to_tau_lepton_antilepton_BR
               * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR
               * CppSLHA::PDG_data::tau_lepton_to_neutrinos_muon_BR ) )
  << std::endl
  << "smuon channel adder pointer = " << smuon_channel
  << std::endl
  << "stau channel adder pointer = " << stau_channel
  << std::endl
  << "Z channel adder pointer = " << Z_channel
  << std::endl
  << std::endl
  << "before asking for *any* values, "
  << std::endl << "two_jets_no_leptons = "
  << *(two_jets_no_leptons.get_value_pointer())
  << std::endl << "one_jet_one_muon = "
  << *(one_jet_one_muon.get_value_pointer())
  << std::endl << "one_jet_one_electron = "
  << *(one_jet_one_electron.get_value_pointer())
  << std::endl << "no_jets_two_muons = "
  << *(no_jets_two_muons.get_value_pointer())
  << std::endl << "no_jets_two_electrons = "
  << *(no_jets_two_electrons.get_value_pointer())
  << std::endl << "no_jets_one_muon_one_electron = "
  << *(no_jets_one_muon_one_electron.get_value_pointer())
  << std::endl << "OSSF_minus_OSDF_leptons = "
  << *(OSSF_minus_OSDF_leptons.get_value_pointer())
  << std::endl << "one_jet_no_leptons = "
  << *(one_jet_no_leptons.get_value_pointer())
  << std::endl << "no_jets_one_muon = "
  << *(no_jets_one_muon.get_value_pointer())
  << std::endl << "no_jets_one_electron = "
  << *(no_jets_one_electron.get_value_pointer())
  << std::endl << "no_jets_no_leptons = "
  << *(no_jets_no_leptons.get_value_pointer());
  std::cout << std::endl;
  std::cout
  << std::endl
  << "OSSF_minus_OSDF_leptons.get_value() = "
  << OSSF_minus_OSDF_leptons.get_value();
  std::cout << std::endl;
  std::cout
  << std::endl
  << "now, "
  << std::endl << "two_jets_no_leptons = "
  << *(two_jets_no_leptons.get_value_pointer())
  << std::endl << "one_jet_one_muon = "
  << *(one_jet_one_muon.get_value_pointer())
  << std::endl << "one_jet_one_electron = "
  << *(one_jet_one_electron.get_value_pointer())
  << std::endl << "no_jets_two_muons = "
  << *(no_jets_two_muons.get_value_pointer())
  << std::endl << "no_jets_two_electrons = "
  << *(no_jets_two_electrons.get_value_pointer())
  << std::endl << "no_jets_one_muon_one_electron = "
  << *(no_jets_one_muon_one_electron.get_value_pointer())
  << std::endl << "OSSF_minus_OSDF_leptons = "
  << *(OSSF_minus_OSDF_leptons.get_value_pointer())
  << std::endl << "one_jet_no_leptons = "
  << *(one_jet_no_leptons.get_value_pointer())
  << std::endl << "no_jets_one_muon = "
  << *(no_jets_one_muon.get_value_pointer())
  << std::endl << "no_jets_one_electron = "
  << *(no_jets_one_electron.get_value_pointer())
  << std::endl << "no_jets_no_leptons = "
  << *(no_jets_no_leptons.get_value_pointer())
  << std::endl;
  std::cout << std::endl;
  std::cout
  << std::endl
  << "no_jets_one_muon.get_value() = "
  << no_jets_one_muon.get_value();
  std::cout << std::endl;
  std::cout
  << std::endl
  << "now, "
  << std::endl << "two_jets_no_leptons = "
  << *(two_jets_no_leptons.get_value_pointer())
  << std::endl << "one_jet_one_muon = "
  << *(one_jet_one_muon.get_value_pointer())
  << std::endl << "one_jet_one_electron = "
  << *(one_jet_one_electron.get_value_pointer())
  << std::endl << "no_jets_two_muons = "
  << *(no_jets_two_muons.get_value_pointer())
  << std::endl << "no_jets_two_electrons = "
  << *(no_jets_two_electrons.get_value_pointer())
  << std::endl << "no_jets_one_muon_one_electron = "
  << *(no_jets_one_muon_one_electron.get_value_pointer())
  << std::endl << "OSSF_minus_OSDF_leptons = "
  << *(OSSF_minus_OSDF_leptons.get_value_pointer())
  << std::endl << "one_jet_no_leptons = "
  << *(one_jet_no_leptons.get_value_pointer())
  << std::endl << "no_jets_one_muon = "
  << *(no_jets_one_muon.get_value_pointer())
  << std::endl << "no_jets_one_electron = "
  << *(no_jets_one_electron.get_value_pointer())
  << std::endl << "no_jets_no_leptons = "
  << *(no_jets_no_leptons.get_value_pointer())
  << std::endl;
  std::cout << std::endl;


  std::cout
  << std::endl
  << "total = "
  << ( two_jets_no_leptons.get_value()
       + one_jet_one_muon.get_value()
       + one_jet_one_electron.get_value()
       + no_jets_two_muons.get_value()
       + no_jets_two_electrons.get_value()
       + no_jets_one_muon_one_electron.get_value()
       + one_jet_no_leptons.get_value()
       + no_jets_one_muon.get_value()
       + no_jets_one_electron.get_value()
       + no_jets_no_leptons.get_value() )
  << std::endl
  << "to be compared to cascade BR sum = "
  << ( 2.0 * ( neutralino_four_to_smuon_R * smuon_R_to_neutralino_one
               + neutralino_four_to_stau_one * stau_one_to_neutralino_one  )
       + neutralino_four_to_Z );
  std::cout << std::endl;

  delete smuon_channel;
  delete stau_channel;
  delete Z_channel;


  // now we call the values in a different order...

  LHC_FASER::cascade_acceptance_value two_jets_no_leptons2( 0.0 );
  LHC_FASER::cascade_acceptance_value one_jet_one_muon2( 0.0 );
  LHC_FASER::cascade_acceptance_value one_jet_one_electron2( 0.0 );
  LHC_FASER::cascade_acceptance_value no_jets_two_muons2( 0.0 );
  LHC_FASER::cascade_acceptance_value no_jets_two_electrons2( 0.0 );
  LHC_FASER::cascade_acceptance_value no_jets_one_muon_one_electron2( 0.0 );
  LHC_FASER::cascade_acceptance_value OSSF_minus_OSDF_leptons2( 0.0 );
  LHC_FASER::cascade_acceptance_value one_jet_no_leptons2( 0.0 );
  LHC_FASER::cascade_acceptance_value no_jets_one_muon2( 0.0 );
  LHC_FASER::cascade_acceptance_value no_jets_one_electron2( 0.0 );
  LHC_FASER::cascade_acceptance_value no_jets_no_leptons2( 0.0 );


  // the following use n1 kinematics but n4 cascades:

  smuon_channel = new LHC_FASER::charge_summed_neutralino_to_semu( 10.0,
                                                                   10.0,
                                                  test_lepton_acceptance_value,
                                              test_input_handler.get_sdown_L(),
                                                                   true,
                                      test_input_handler.get_neutralino_four(),
                                                                   true,
                                              test_input_handler.get_smuon_R(),
                                                           &test_input_handler,
                                  OSSF_minus_OSDF_leptons2.get_value_pointer(),
                                        no_jets_two_muons2.get_value_pointer(),
                                         no_jets_one_muon2.get_value_pointer(),
                                     no_jets_no_leptons2.get_value_pointer() );
  OSSF_minus_OSDF_leptons2.add_channel( smuon_channel );
  no_jets_two_muons2.add_channel( smuon_channel );
  no_jets_one_muon2.add_channel( smuon_channel );
  no_jets_no_leptons2.add_channel( smuon_channel );

  stau_channel = new LHC_FASER::charge_summed_neutralino_to_stau( 10.0,
                                                                  10.0,
                                                                  40.0,
                                                  test_lepton_acceptance_value,
                                              test_input_handler.get_sdown_L(),
                                                                    true,
                                      test_input_handler.get_neutralino_four(),
                                                                    true,
                                             test_input_handler.get_stau_one(),
                                                           &test_input_handler,
                                      two_jets_no_leptons2.get_value_pointer(),
                                         one_jet_one_muon2.get_value_pointer(),
                                     one_jet_one_electron2.get_value_pointer(),
                                        no_jets_two_muons2.get_value_pointer(),
                                    no_jets_two_electrons2.get_value_pointer(),
                            no_jets_one_muon_one_electron2.get_value_pointer(),
                                       one_jet_no_leptons2.get_value_pointer(),
                                         no_jets_one_muon2.get_value_pointer(),
                                     no_jets_one_electron2.get_value_pointer(),
                                     no_jets_no_leptons2.get_value_pointer() );
  two_jets_no_leptons2.add_channel( stau_channel );
  one_jet_one_muon2.add_channel( stau_channel );
  one_jet_one_electron2.add_channel( stau_channel );
  no_jets_two_muons2.add_channel( stau_channel );
  no_jets_two_electrons2.add_channel( stau_channel );
  no_jets_one_muon_one_electron2.add_channel( stau_channel );
  one_jet_no_leptons2.add_channel( stau_channel );
  no_jets_one_muon2.add_channel( stau_channel );
  no_jets_one_electron2.add_channel( stau_channel );
  no_jets_no_leptons2.add_channel( stau_channel );

  Z_channel = new LHC_FASER::charge_summed_neutralino_to_Z( 10.0,
                                                            10.0,
                                                            40.0,
                                                  test_lepton_acceptance_value,
                                              test_input_handler.get_sdown_L(),
                                                            true,
                                      test_input_handler.get_neutralino_four(),
                                                            true,
                                                    test_input_handler.get_Z(),
                                                           &test_input_handler,
                                      two_jets_no_leptons2.get_value_pointer(),
                                         one_jet_one_muon2.get_value_pointer(),
                                     one_jet_one_electron2.get_value_pointer(),
                                  OSSF_minus_OSDF_leptons2.get_value_pointer(),
                                        no_jets_two_muons2.get_value_pointer(),
                                    no_jets_two_electrons2.get_value_pointer(),
                            no_jets_one_muon_one_electron2.get_value_pointer(),
                                       one_jet_no_leptons2.get_value_pointer(),
                                         no_jets_one_muon2.get_value_pointer(),
                                     no_jets_one_electron2.get_value_pointer(),
                                     no_jets_no_leptons2.get_value_pointer() );
  two_jets_no_leptons2.add_channel( Z_channel );
  one_jet_one_muon2.add_channel( Z_channel );
  one_jet_one_electron2.add_channel( Z_channel );
  OSSF_minus_OSDF_leptons2.add_channel( Z_channel );
  no_jets_two_muons2.add_channel( Z_channel );
  no_jets_two_electrons2.add_channel( Z_channel );
  no_jets_one_muon_one_electron2.add_channel( Z_channel );
  one_jet_no_leptons2.add_channel( Z_channel );
  no_jets_one_muon2.add_channel( Z_channel );
  no_jets_one_electron2.add_channel( Z_channel );
  no_jets_no_leptons2.add_channel( Z_channel );

  std::cout
  << std::endl
  << "now the values are called in the reverse order..."
  << std::endl
  << "smuon channel adder pointer = " << smuon_channel
  << std::endl
  << "stau channel adder pointer = " << stau_channel
  << std::endl
  << "Z channel adder pointer = " << Z_channel
  << std::endl
  << std::endl
  << "before asking for *any* values, "
  << std::endl << "two_jets_no_leptons2 = "
  << *(two_jets_no_leptons2.get_value_pointer())
  << std::endl << "one_jet_one_muon2 = "
  << *(one_jet_one_muon2.get_value_pointer())
  << std::endl << "one_jet_one_electron2 = "
  << *(one_jet_one_electron2.get_value_pointer())
  << std::endl << "no_jets_two_muons2 = "
  << *(no_jets_two_muons2.get_value_pointer())
  << std::endl << "no_jets_two_electrons2 = "
  << *(no_jets_two_electrons2.get_value_pointer())
  << std::endl << "no_jets_one_muon_one_electron2 = "
  << *(no_jets_one_muon_one_electron2.get_value_pointer())
  << std::endl << "OSSF_minus_OSDF_leptons2 = "
  << *(OSSF_minus_OSDF_leptons2.get_value_pointer())
  << std::endl << "one_jet_no_leptons2 = "
  << *(one_jet_no_leptons2.get_value_pointer())
  << std::endl << "no_jets_one_muon2 = "
  << *(no_jets_one_muon2.get_value_pointer())
  << std::endl << "no_jets_one_electron2 = "
  << *(no_jets_one_electron2.get_value_pointer())
  << std::endl << "no_jets_no_leptons2 = "
  << *(no_jets_no_leptons2.get_value_pointer());
  std::cout << std::endl;
  std::cout
  << std::endl
  << "one_jet_one_muon2.get_value() = "
  << one_jet_one_muon2.get_value();
  std::cout << std::endl;
  std::cout
  << std::endl
  << "now, "
  << std::endl << "two_jets_no_leptons2 = "
  << *(two_jets_no_leptons2.get_value_pointer())
  << std::endl << "one_jet_one_muon2 = "
  << *(one_jet_one_muon2.get_value_pointer())
  << std::endl << "one_jet_one_electron2 = "
  << *(one_jet_one_electron2.get_value_pointer())
  << std::endl << "no_jets_two_muons2 = "
  << *(no_jets_two_muons2.get_value_pointer())
  << std::endl << "no_jets_two_electrons2 = "
  << *(no_jets_two_electrons2.get_value_pointer())
  << std::endl << "no_jets_one_muon_one_electron2 = "
  << *(no_jets_one_muon_one_electron2.get_value_pointer())
  << std::endl << "OSSF_minus_OSDF_leptons2 = "
  << *(OSSF_minus_OSDF_leptons2.get_value_pointer())
  << std::endl << "one_jet_no_leptons2 = "
  << *(one_jet_no_leptons2.get_value_pointer())
  << std::endl << "no_jets_one_muon2 = "
  << *(no_jets_one_muon2.get_value_pointer())
  << std::endl << "no_jets_one_electron2 = "
  << *(no_jets_one_electron2.get_value_pointer())
  << std::endl << "no_jets_no_leptons2 = "
  << *(no_jets_no_leptons2.get_value_pointer());
  std::cout << std::endl;
  std::cout
  << std::endl
  << "OSSF_minus_OSDF_leptons2.get_value() = "
  << OSSF_minus_OSDF_leptons2.get_value();
  std::cout << std::endl;
  std::cout
  << std::endl
  << "now, "
  << std::endl << "two_jets_no_leptons2 = "
  << *(two_jets_no_leptons2.get_value_pointer())
  << std::endl << "one_jet_one_muon2 = "
  << *(one_jet_one_muon2.get_value_pointer())
  << std::endl << "one_jet_one_electron2 = "
  << *(one_jet_one_electron2.get_value_pointer())
  << std::endl << "no_jets_two_muons2 = "
  << *(no_jets_two_muons2.get_value_pointer())
  << std::endl << "no_jets_two_electrons2 = "
  << *(no_jets_two_electrons2.get_value_pointer())
  << std::endl << "no_jets_one_muon_one_electron2 = "
  << *(no_jets_one_muon_one_electron2.get_value_pointer())
  << std::endl << "OSSF_minus_OSDF_leptons2 = "
  << *(OSSF_minus_OSDF_leptons2.get_value_pointer())
  << std::endl << "one_jet_no_leptons2 = "
  << *(one_jet_no_leptons2.get_value_pointer())
  << std::endl << "no_jets_one_muon2 = "
  << *(no_jets_one_muon2.get_value_pointer())
  << std::endl << "no_jets_one_electron2 = "
  << *(no_jets_one_electron2.get_value_pointer())
  << std::endl << "no_jets_no_leptons2 = "
  << *(no_jets_no_leptons2.get_value_pointer())
  << std::endl;
  std::cout << std::endl;


  std::cout
  << std::endl
  << "total = "
  << ( two_jets_no_leptons2.get_value()
       + one_jet_one_muon2.get_value()
       + one_jet_one_electron2.get_value()
       + no_jets_two_muons2.get_value()
       + no_jets_two_electrons2.get_value()
       + no_jets_one_muon_one_electron2.get_value()
       + one_jet_no_leptons2.get_value()
       + no_jets_one_muon2.get_value()
       + no_jets_one_electron2.get_value()
       + no_jets_no_leptons2.get_value() )
  << std::endl
  << "to be compared to cascade BR sum = "
  << ( 2.0 * ( neutralino_four_to_smuon_R * smuon_R_to_neutralino_one
               + neutralino_four_to_stau_one * stau_one_to_neutralino_one  )
       + neutralino_four_to_Z );
  std::cout << std::endl;

  delete smuon_channel;
  delete stau_channel;
  delete Z_channel;

  LHC_FASER::cascade_handler*
  cascade_handler_tester
   = new LHC_FASER::cascade_handler( &test_input_handler );

  LHC_FASER::charge_summed_cascade_calculator*
  charge_summed_test_10_10_40_a
  = cascade_handler_tester->get_charge_summed_calculator( 10.0,
                                                          10.0,
                                                          40.0,
                                                  test_lepton_acceptance_value,
                                              test_input_handler.get_sdown_L(),
                                                          true,
                                      test_input_handler.get_neutralino_four(),
                                                          true );
  LHC_FASER::charge_summed_cascade_calculator*
  charge_summed_test_20_10_40
  = cascade_handler_tester->get_charge_summed_calculator( 20.0,
                                                          10.0,
                                                          40.0,
                                                  test_lepton_acceptance_value,
                                              test_input_handler.get_sdown_L(),
                                                          true,
                                      test_input_handler.get_neutralino_four(),
                                                          true );
  LHC_FASER::charge_summed_cascade_calculator*
  charge_summed_test_10_10_40_b
  = cascade_handler_tester->get_charge_summed_calculator( 10.0,
                                                          10.0,
                                                          40.0,
                                                  test_lepton_acceptance_value,
                                              test_input_handler.get_sdown_L(),
                                                          true,
                                      test_input_handler.get_neutralino_four(),
                                                          true );

  std::cout
  << std::endl << "cascade handler testing:"
  << std::endl
  << "charge_summed_test_10_10_40_a = " << charge_summed_test_10_10_40_a
  << std::endl
  << "charge_summed_test_10_10_40_b = " << charge_summed_test_10_10_40_b
  << std::endl
  << "charge_summed_test_20_10_40 = " << charge_summed_test_20_10_40
  << std::endl;
  std::cout
  << std::endl
  << "charge_summed_test_10_10_40_a:"
  << std::endl << "OSSF-OSDF = "
  << charge_summed_test_10_10_40_a->get_OSSF_minus_OSDF_acceptance()
  << std::endl << "2j0l = "
  << charge_summed_test_10_10_40_a->get_acceptance( 2,
                                                    0 )
  << std::endl << "1j1l = "
  << charge_summed_test_10_10_40_a->get_acceptance( 1,
                                                    1 )
  << std::endl << "0j2l = "
  << charge_summed_test_10_10_40_a->get_acceptance( 0,
                                                    2 )
  << std::endl << "0j2mu = "
  << charge_summed_test_10_10_40_a->get_muon_acceptance( 0,
                                                         2 )
  << std::endl << "1j0l = "
  << charge_summed_test_10_10_40_a->get_acceptance( 1,
                                                    0 )
  << std::endl << "0j1l = "
  << charge_summed_test_10_10_40_a->get_acceptance( 0,
                                                    1 )
  << std::endl << "0j1mu = "
  << charge_summed_test_10_10_40_a->get_muon_acceptance( 0,
                                                         1 )
  << std::endl << "0j0l = "
  << charge_summed_test_10_10_40_a->get_acceptance( 0,
                                                    0 );
  std::cout
  << std::endl
  << "charge_summed_test_20_10_40:"
  << std::endl << "OSSF-OSDF = "
  << charge_summed_test_20_10_40->get_OSSF_minus_OSDF_acceptance()
  << std::endl << "2j0l = "
  << charge_summed_test_20_10_40->get_acceptance( 2,
                                                  0 )
  << std::endl << "1j1l = "
  << charge_summed_test_20_10_40->get_acceptance( 1,
                                                  1 )
  << std::endl << "0j2l = "
  << charge_summed_test_20_10_40->get_acceptance( 0,
                                                  2 )
  << std::endl << "0j2mu = "
  << charge_summed_test_20_10_40->get_muon_acceptance( 0,
                                                       2 )
  << std::endl << "1j0l = "
  << charge_summed_test_20_10_40->get_acceptance( 1,
                                                  0 )
  << std::endl << "0j1l = "
  << charge_summed_test_20_10_40->get_acceptance( 0,
                                                  1 )
  << std::endl << "0j1mu = "
  << charge_summed_test_20_10_40->get_muon_acceptance( 0,
                                                       1 )
  << std::endl << "0j0l = "
  << charge_summed_test_20_10_40->get_acceptance( 0,
                                                  0 );

  delete test_cascade_handler;

  **/

	

  /* testing the signal stuff: shouldnt work, because cascade calculation doesnt work (this requires all objects to work, which is currently not the case)*//**

  LHC_FASER::colored_cascade_set*
  gluino_BR_sum_tester
  = test_input_handler.get_colored_cascade_set(
                                             test_input_handler.get_gluino() );
  double gluino_BR_sum = 0.0;
  LHC_FASER::exclusive_BR_calculator* BR_sum_helper;
  double BR_sum_helper_double;

  for( std::vector< LHC_FASER::colored_cascade* >::const_iterator
       BR_iterator = gluino_BR_sum_tester->get_colored_cascades()->begin();
       gluino_BR_sum_tester->get_colored_cascades()->end() > BR_iterator;
       ++BR_iterator )
    {

      if( (*BR_iterator)->is_gx() )
        // if it's a direct gluino decay to the electroweakino...
        {

          BR_sum_helper
          = test_input_handler.get_exclusive_BR_calculator(
                                               test_input_handler.get_gluino(),
                                                   (*BR_iterator)->get_EWino(),
                                                            true,
                                         test_input_handler.get_empty_list() );

          if( NULL != BR_sum_helper )
            {

              gluino_BR_sum += BR_sum_helper->get_BR();

            }

          if( !((*BR_iterator)->get_EWino()->counts_as_self_conjugate()) )
            {

              BR_sum_helper
              = test_input_handler.get_exclusive_BR_calculator(
                                               test_input_handler.get_gluino(),
                                                   (*BR_iterator)->get_EWino(),
                                                                false,
                                         test_input_handler.get_empty_list() );

              if( NULL != BR_sum_helper )
                {

                  gluino_BR_sum += BR_sum_helper->get_BR();

                }

            }

        }
      else if( (*BR_iterator)->is_gsx() )
        // otherwise if it's a gluino 2-stage decay to the electroweakino...
        {

          BR_sum_helper
          = test_input_handler.get_exclusive_BR_calculator(
                                               test_input_handler.get_gluino(),
                                                 (*BR_iterator)->get_product(),
                                                            true,
                                         test_input_handler.get_empty_list() );
          if( NULL != BR_sum_helper )
            {

              BR_sum_helper_double = BR_sum_helper->get_BR();

              BR_sum_helper
              = test_input_handler.get_exclusive_BR_calculator(
                                                 (*BR_iterator)->get_product(),
                                                   (*BR_iterator)->get_EWino(),
                                                                true,
                                         test_input_handler.get_empty_list() );

              if( NULL != BR_sum_helper )
                {

                  gluino_BR_sum
                  += ( 2.0 * BR_sum_helper_double * BR_sum_helper->get_BR() );

                }

              if( !((*BR_iterator)->get_EWino()->counts_as_self_conjugate()) )
                {

                  BR_sum_helper
                  = test_input_handler.get_exclusive_BR_calculator(
                                                 (*BR_iterator)->get_product(),
                                                   (*BR_iterator)->get_EWino(),
                                                                    false,
                                         test_input_handler.get_empty_list() );

                  if( NULL != BR_sum_helper )
                    {

                      gluino_BR_sum += ( 2.0 * BR_sum_helper_double
                                             * BR_sum_helper->get_BR() );

                    }

                }

            }

        }

    }
  std::cout
  << std::endl
  << "gluino_BR_sum = " << gluino_BR_sum;
  std::cout << std::endl;



  LHC_FASER::cross_section_handler test_cross_section_handler(
                                                         &test_input_handler );

  LHC_FASER::kinematics_handler test_kinematics_handler( &test_input_handler );

  LHC_FASER::cascade_handler test_cascade_handler( &test_input_handler );
  LHC_FASER::signal_shortcuts test_signal_shortcuts( &test_input_handler,
                                                   &test_cross_section_handler,
                                                     &test_kinematics_handler,
                                                     &test_cascade_handler );

  LHC_FASER::signal_handler test_bad_signal( "testBadSignal",
                                             &test_signal_shortcuts );

  std::cout
  << std::endl
  << *(test_bad_signal.get_name()) << " has signal value "
  << test_bad_signal.get_value();
  std::cout << std::endl;

  LHC_FASER::signal_handler
  cross_section_breakdown_test( "sigmaBreakdownTest",
                                &test_signal_shortcuts );

  std::cout
  << std::endl
  << *(cross_section_breakdown_test.get_name()) << " has signal value "
  << cross_section_breakdown_test.get_value();
  std::cout << std::endl;

  LHC_FASER::signal_handler test_Atlas_signal( "Atlas4jMET0l7TeV",
                                               &test_signal_shortcuts );

  std::cout
  << std::endl
  << *(test_Atlas_signal.get_name()) << " constructed without crashing.";
  std::cout
  << std::endl;
  std::cout
  << std::endl
  << *(test_Atlas_signal.get_name()) << " has signal value "
  << test_Atlas_signal.get_value();
  std::cout << std::endl;

  SLHA_data.read_file( "SPS2_spectrum.out" );
  std::cout << std::endl;
  std::cout
  << std::endl << "without updating for new point, "
  << *(cross_section_breakdown_test.get_name()) << " has signal value "
  << cross_section_breakdown_test.get_value()
  << std::endl
  << *(test_Atlas_signal.get_name()) << " has signal value "
  << test_Atlas_signal.get_value();
  std::cout << std::endl;

  test_readier.ready_for_new_point();
  std::cout
  << std::endl << "after updating for new point, "
  << *(cross_section_breakdown_test.get_name()) << " has signal value "
  << cross_section_breakdown_test.get_value()
  << std::endl
  << *(test_Atlas_signal.get_name()) << " has signal value "
  << test_Atlas_signal.get_value();
  std::cout << std::endl;

  gluino_BR_sum = 0.0;
  for( std::vector< LHC_FASER::colored_cascade* >::const_iterator
       BR_iterator = gluino_BR_sum_tester->get_colored_cascades()->begin();
       gluino_BR_sum_tester->get_colored_cascades()->end() > BR_iterator;
       ++BR_iterator )
    {

      if( (*BR_iterator)->is_gx() )
        // if it's a direct gluino decay to the electroweakino...
        {

          BR_sum_helper
          = test_input_handler.get_exclusive_BR_calculator(
                                               test_input_handler.get_gluino(),
                                                   (*BR_iterator)->get_EWino(),
                                                            true,
                                         test_input_handler.get_empty_list() );

          if( NULL != BR_sum_helper )
            {

              gluino_BR_sum += BR_sum_helper->get_BR();

            }

          if( !((*BR_iterator)->get_EWino()->counts_as_self_conjugate()) )
            {

              BR_sum_helper
              = test_input_handler.get_exclusive_BR_calculator(
                                               test_input_handler.get_gluino(),
                                                   (*BR_iterator)->get_EWino(),
                                                                false,
                                         test_input_handler.get_empty_list() );

              if( NULL != BR_sum_helper )
                {

                  gluino_BR_sum += BR_sum_helper->get_BR();

                }

            }

        }
      else if( (*BR_iterator)->is_gsx() )
        // otherwise if it's a gluino 2-stage decay to the electroweakino...
        {

          BR_sum_helper
          = test_input_handler.get_exclusive_BR_calculator(
                                               test_input_handler.get_gluino(),
                                                 (*BR_iterator)->get_product(),
                                                            true,
                                         test_input_handler.get_empty_list() );
          if( NULL != BR_sum_helper )
            {

              BR_sum_helper_double = BR_sum_helper->get_BR();

              BR_sum_helper
              = test_input_handler.get_exclusive_BR_calculator(
                                                 (*BR_iterator)->get_product(),
                                                   (*BR_iterator)->get_EWino(),
                                                                true,
                                         test_input_handler.get_empty_list() );

              if( NULL != BR_sum_helper )
                {

                  gluino_BR_sum
                  += ( 2.0 * BR_sum_helper_double * BR_sum_helper->get_BR() );

                }

              if( !((*BR_iterator)->get_EWino()->counts_as_self_conjugate()) )
                {

                  BR_sum_helper
                  = test_input_handler.get_exclusive_BR_calculator(
                                                 (*BR_iterator)->get_product(),
                                                   (*BR_iterator)->get_EWino(),
                                                                    false,
                                         test_input_handler.get_empty_list() );

                  if( NULL != BR_sum_helper )
                    {

                      gluino_BR_sum += ( 2.0 * BR_sum_helper_double
                                             * BR_sum_helper->get_BR() );

                    }

                }

            }

        }

    }
  std::cout
  << std::endl
  << "gluino_BR_sum = " << gluino_BR_sum;
  std::cout << std::endl;

  std::cout
  << std::endl
  << "now trying alternating SPS1a & SPS2 4 times.";
  std::cout << std::endl;

  timeval starting_time;
  timeval ending_time;
  double elapsed_time;
  for( int reload_count = 1;
       4 >= reload_count;
       ++reload_count )
    {

      gettimeofday( &starting_time,
                    NULL );

      std::cout
      << std::endl << "count = " << reload_count;

      SLHA_data.read_file( "SPS1a_spectrum.out" );
      test_readier.ready_for_new_point();
      std::cout
      << std::endl << "SPS1a: "
      << *(test_Atlas_signal.get_name()) << " has signal value "
      << test_Atlas_signal.get_value();
      std::cout << std::endl;

      SLHA_data.read_file( "SPS2_spectrum.out" );
      gettimeofday( &ending_time,
                    NULL );
      elapsed_time
      = ( (double)( ending_time.tv_sec - starting_time.tv_sec )
          + 0.000001 * (double)( ending_time.tv_usec
                                 - starting_time.tv_usec ) );
      std::cout
      << std::endl
      << "took " << seconds_taken << " seconds.";
      std::cout << std::endl;

      gettimeofday( &starting_time,
                    NULL );

      test_readier.ready_for_new_point();
      std::cout
      << std::endl << "SPS2: "
      << *(test_Atlas_signal.get_name()) << " has signal value "
      << test_Atlas_signal.get_value();
      std::cout << std::endl;

      gettimeofday( &ending_time,
                    NULL );
      elapsed_time
      = ( (double)( ending_time.tv_sec - starting_time.tv_sec )
          + 0.000001 * (double)( ending_time.tv_usec
                                 - starting_time.tv_usec ) );
      std::cout
      << std::endl
      << "took " << elapsed_time << " seconds.";
      std::cout << std::endl;

    }


   **/


  /* finally, testing the LHC_FASER_UI itself: (not working)*//**

  LHC_FASER::LHC_FASER_UI SPS1a_UI( "SPS1a_spectrum.out" );
  LHC_FASER::LHC_FASER_UI SPS2_UI( "SPS2_spectrum.out" );

  SPS1a_UI.add_signal( "sigmaBreakdownTest" );
  SPS1a_UI.add_signal( "Atlas4jMET0l7TeV" );
  SPS2_UI.add_signal( "badSignalName" );
  SPS2_UI.add_signal( "Atlas4jMET0l7TeV" );

  LHC_FASER::signal_handler* sigma_SPS1a
  = SPS1a_UI.get_signal( "sigmaBreakdownTest" );
  LHC_FASER::signal_handler* Atlas_SPS1a
  = SPS1a_UI.get_signal( "Atlas4jMET0l7TeV" );
  LHC_FASER::signal_handler* Atlas_SPS2
  = SPS2_UI.get_signal( "Atlas4jMET0l7TeV" );
  LHC_FASER::signal_handler* bad_SPS2
  = SPS2_UI.get_signal( "badSignalName" );


  std::cout << std::endl;
  std::cout << std::endl << "all 4 test signals from test UIs have pointers.";
  std::cout << std::endl;

  SPS1a_UI.update_for_new_point();
  std::cout << std::endl;
  std::cout << std::endl << "SPS1a_UI.update_for_new_point() successful.";
  std::cout << std::endl;
  SPS2_UI.update_for_new_point();
  std::cout << std::endl;
  std::cout << std::endl << "SPS2_UI.update_for_new_point() successful.";
  std::cout << std::endl;
  std::cout
  << std::endl << "after updating for new point, "
  << *(sigma_SPS1a->get_name()) << " (sigma_SPS1a) has signal value "
  << sigma_SPS1a->get_value();
  std::cout << std::endl;
  std::cout
  << *(bad_SPS2->get_name()) << " (bad_SPS2) has signal value "
  << bad_SPS2->get_value();
  std::cout << std::endl;
  std::cout
  << *(Atlas_SPS1a->get_name()) << " (Atlas_SPS1a) has signal value "
  << Atlas_SPS1a->get_value();
  std::cout << std::endl;
  std::cout
  << *(Atlas_SPS2->get_name()) << " (Atlas_SPS2) has signal value "
  << Atlas_SPS2->get_value();
  std::cout << std::endl;

  **/


  /* testing adding new signals & how long it takes to load & calculate: *//**

  timeval start_time;
  timeval end_time;
  double seconds_taken;

  std::vector< LHC_FASER::signalHandler* > signals;
  std::string signal_name;

  //std::string grid_path( "/home/bol/BOL_work/projects/Eclipse_workspace/");
  //grid_path.append( "LHC-FASER_unit_tester/testing/grids/" );
  std::string grid_path( "./grids/");

  system( "cp LM1_spectrum.out test_spectrum.out" );

  gettimeofday( &start_time,
                NULL );
  std::cout
  << std::endl << "constructing test_UI";
  LHC_FASER::LHC_FASER_UI test_UI( "test_spectrum.out",
                                   grid_path,
                                   "fb" );
  gettimeofday( &end_time,
                NULL );
  seconds_taken
  = ( (double)( end_time.tv_sec - start_time.tv_sec )
      + 0.000001 * (double)( end_time.tv_usec - start_time.tv_usec ) );
  std::cout
  << " took " << seconds_taken << " seconds.";
  std::cout << std::endl;

  gettimeofday( &start_time,
                NULL );
  signal_name.assign( "Atlas4jMET0l7TeV" );
  std::cout
  << std::endl << "adding " << signal_name;
  test_UI.add_signal( signal_name );
  signals.push_back( test_UI.get_signal( signal_name ) );
  gettimeofday( &end_time,
                NULL );
  seconds_taken
  = ( (double)( end_time.tv_sec - start_time.tv_sec )
      + 0.000001 * (double)( end_time.tv_usec - start_time.tv_usec ) );
  std::cout
  << " took " << seconds_taken << " seconds.";
  std::cout << std::endl;

  gettimeofday( &start_time,
                NULL );
  signal_name.assign( "Atlas3jMET1l7TeV" );
  std::cout
  << std::endl << "adding " << signal_name;
  test_UI.add_signal( signal_name );
  signals.push_back( test_UI.get_signal( signal_name ) );
  gettimeofday( &end_time,
                NULL );
  seconds_taken
  = ( (double)( end_time.tv_sec - start_time.tv_sec )
      + 0.000001 * (double)( end_time.tv_usec - start_time.tv_usec ) );
  std::cout
  << " took " << seconds_taken << " seconds.";
  std::cout << std::endl;


  std::vector< std::string* > spectrum_files;
  spectrum_files.push_back( new std::string( "SPS1a_spectrum.out" ) );
  spectrum_files.push_back( new std::string( "SPS1b_spectrum.out" ) );
  spectrum_files.push_back( new std::string( "SPS2_spectrum.out" ) );
  spectrum_files.push_back( new std::string( "SPS3_spectrum.out" ) );
  //spectrum_files.push_back( new std::string( "SPS4_spectrum.out" ) );
  // SPS4 is borked.
  spectrum_files.push_back( new std::string( "SPS5_spectrum.out" ) );
  spectrum_files.push_back( new std::string( "SPS6_spectrum.out" ) );
  spectrum_files.push_back( new std::string( "SPS7_spectrum.out" ) );
  spectrum_files.push_back( new std::string( "SPS8_spectrum.out" ) );
  spectrum_files.push_back( new std::string( "SPS9_spectrum.out" ) );
  spectrum_files.push_back( new std::string( "LM1_spectrum.out" ) );
  spectrum_files.push_back( new std::string( "SU4_spectrum.out" ) );

  std::string system_call_string;
  for( int reload_count = 1;
       recalculation_amount >= reload_count;
       ++reload_count )
    {

      std::cout
      << std::endl << "count = " << reload_count;

      for( std::vector< std::string* >::const_iterator
           spectrum_iterator = spectrum_files.begin();
           spectrum_files.end() > spectrum_iterator;
           ++spectrum_iterator )
        {

          system_call_string.assign( "cp " );
          system_call_string.append( *(*spectrum_iterator) );
          system_call_string.append( " test_spectrum.out" );
          gettimeofday( &start_time,
                        NULL );
          std::cout
          << std::endl << "copying " << *(*spectrum_iterator);
          system( system_call_string.c_str() );
          gettimeofday( &end_time,
                        NULL );
          seconds_taken
          = ( (double)( end_time.tv_sec - start_time.tv_sec )
              + 0.000001 * (double)( end_time.tv_usec - start_time.tv_usec ) );
          std::cout
          << " took " << seconds_taken << " seconds.";
          std::cout << std::endl;

          gettimeofday( &start_time,
                        NULL );
          std::cout
          << std::endl << "reading it in & updating";
          test_UI.update_for_new_point();
          gettimeofday( &end_time,
                        NULL );
          seconds_taken
          = ( (double)( end_time.tv_sec - start_time.tv_sec )
              + 0.000001 * (double)( end_time.tv_usec - start_time.tv_usec ) );
          std::cout
          << " took " << seconds_taken << " seconds.";
          std::cout << std::endl;

          for( std::vector< LHC_FASER::signalHandler* >::iterator
               signal_iterator = signals.begin();
               signals.end() > signal_iterator;
               ++signal_iterator )
            {

              gettimeofday( &start_time,
                            NULL );
              std::cout
              << std::endl << "calculating the point's value for "
              << *((*signal_iterator)->getName()) << ": ";
              std::cout << (*signal_iterator)->getValue();
              gettimeofday( &end_time,
                            NULL );
              seconds_taken
              = ( (double)( end_time.tv_sec - start_time.tv_sec )
                  + 0.000001 * (double)( end_time.tv_usec
                                         - start_time.tv_usec ) );
              std::cout
              << "; took " << seconds_taken << " seconds.";
              std::cout << std::endl;

            }

        }

    }

  for( std::vector< std::string* >::const_iterator
       deletion_iterator = spectrum_files.begin();
       spectrum_files.end() > deletion_iterator;
       ++deletion_iterator )
    {

      delete *deletion_iterator;

    }

  **/


  std::cout << std::endl;
  std::cout << std::endl << "ended successfully, I hope";
  std::cout << std::endl;

  return EXIT_SUCCESS;

}
