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
 * 2: getsReadiedForNewPoint + readierForNewPoint - not yet redone, but
 *                                                    should not be a problem,
 *                                                    since I don't *think*
 *                                                    that I changed anything
 *                                                    here since 0.1 (which
 *                                                    works)...
 * 3: crossSectionTable - again, not yet redone, but I don't recall changing
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
//#include "LHC-FASER_derived_electroweak_cascade_stuff.hpp"
//#include "LHC-FASER_signal_data_collection_stuff.hpp"
//#include "LHC-FASER_signal_calculator_stuff.hpp"
//#include "LHC-FASER.hpp"
#include <sys/time.h>


// this holds the set of distributions & handles adding them.
class distribution_set
{

public:

  distribution_set( LHC_FASER::readierForNewPoint* const readier ) :
    readier( readier ),
    hard_muon( new LHC_FASER::hardMuonFromTau() ),
    soft_muon( new LHC_FASER::softMuonFromTau() ),
    hard_pion( new LHC_FASER::hardPionFromTau() ),
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
  LHC_FASER::readierForNewPoint* readier;
  LHC_FASER::tauDecayCoefficient const* const hard_muon;
  LHC_FASER::tauDecayCoefficient const* const soft_muon;
  LHC_FASER::tauDecayCoefficient const* const hard_pion;
  LHC_FASER::tauDecayCoefficient const* const soft_pion;

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

int main( int argumentCount,
          char* argumentStrings[] )
{

  int recalculation_amount = 10;

  if( 2 != argumentCount )
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

      std::string argument_interpretter( argumentStrings[ 1 ] );
      std::stringstream argument_stream( argument_interpretter );
      argument_stream >> recalculation_amount;

    }

  CppSLHA::CppSLHA2 SLHA_data( "SPS1a_spectrum.out" );
  SLHA_data.read_file();
  LHC_FASER::readierForNewPoint test_readier;
  LHC_FASER::inputHandler test_input_handler( &SLHA_data,
                                         SLHA_data.get_particle_spectrum(),
                                         "./grids",
                                         &test_readier );

  std::string bluh( "bluh" );

  /* testing the lepton energy distributions in the squark rest frame: (slight changes of input) *//**/

  distribution_set distributions( &test_readier );

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
                                         sdownLPointer,
                                         &sdownLMass,
                                         neutralinoThreePointer,
                                         zPointer,
                                         neutralinoOnePointer );
  distributions.add_distributions( currentDistribution,
                                   "Zn4" );
  currentDistribution
  = new LHC_FASER::Z_muon_plus_antimuon( &test_readier,
                                         &SLHA_data,
                                         sdownLPointer,
                                         &sdownLMass,
                                         neutralinoFourPointer,
                                         zPointer,
                                         neutralinoOnePointer );
  distributions.add_distributions( currentDistribution,
                                   "Zn3" );
  currentDistribution
  = new LHC_FASER::W_minus_plus_up_L_type_muon( &test_readier,
                                                &SLHA_data,
                                                sdownLPointer,
                                                &sdownLMass,
                                                charginoOnePointer,
                                                W_pointer,
                                                neutralinoOnePointer );
  distributions.add_distributions( currentDistribution,
                                   "LWx1" );
  currentDistribution
  = new LHC_FASER::W_minus_plus_up_R_type_muon( &test_readier,
                                                &SLHA_data,
                                                sdownLPointer,
                                                &sdownLMass,
                                                charginoOnePointer,
                                                W_pointer,
                                                neutralinoOnePointer );
  distributions.add_distributions( currentDistribution,
                                   "RWx1" );
  currentDistribution
  = new LHC_FASER::W_minus_plus_up_L_type_muon( &test_readier,
                                                &SLHA_data,
                                                sdownLPointer,
                                                &sdownLMass,
                                                charginoTwoPointer,
                                                W_pointer,
                                                neutralinoOnePointer );
  distributions.add_distributions( currentDistribution,
                                   "LWx2" );
  currentDistribution
  = new LHC_FASER::W_minus_plus_up_R_type_muon( &test_readier,
                                                &SLHA_data,
                                                sdownLPointer,
                                                &sdownLMass,
                                                charginoTwoPointer,
                                                W_pointer,
                                                neutralinoOnePointer );
  distributions.add_distributions( currentDistribution,
                                   "RWx2" );
  currentDistribution
  = new LHC_FASER::neutralino_three_body_decay( &test_readier,
                                                &SLHA_data,
                                                sdownLPointer,
                                                &sdownLMass,
                                                neutralinoTwoPointer,
                                                neutralinoOnePointer,
                                                smuonLPointer,
                                                smuonRPointer );
  distributions.add_distributions( currentDistribution,
                                   "nn3b" );
  currentDistribution
  = new LHC_FASER::chargino_three_body_decay( &test_readier,
                                              &SLHA_data,
                                              sdownLPointer,
                                              &sdownLMass,
                                              charginoOnePointer,
                                              neutralinoOnePointer,
                                              smuonLPointer,
                                              smuonRPointer,
                                              muonSneutrinoLPointer,
                                              muonSneutrinoRPointer );
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
  << "is verbose? " << test_input_handler.isVerbose()
  << std::endl
  << "set to be not verbose.";
  test_input_handler.setVerbosity( false );
  std::cout
  << std::endl
  << "is verbose? " << test_input_handler.isVerbose()
  << std::endl
  << "set to be verbose.";
  test_input_handler.setVerbosity( true );
  std::cout
  << std::endl
  << "is verbose? " << test_input_handler.isVerbose()
  << std::endl
  << "using NLO? " << test_input_handler.usingNlo()
  << std::endl
  << "set to not use NLO.";
  test_input_handler.setUsingNloFlag( false );
  std::cout
  << std::endl
  << "using NLO? " << test_input_handler.usingNlo()
  << std::endl
  << "set to use NLO.";
  test_input_handler.setUsingNloFlag( true );
  std::cout
  << std::endl
  << "using NLO? " << test_input_handler.usingNlo()
  << std::endl
  << std::endl
  << "pathToGrids = "
  << *(test_input_handler.getPathToGrids())
  << std::endl
  << "pathToCrossSectionGrids = "
  << *(test_input_handler.getPathToCrossSectionGrids())
  << std::endl
  << "(inspect: pathToCrossSectionGrids = "
  << *(test_input_handler.inspectPathToCrossSectionGrids()) << ")"
  << std::endl
  << "pathToKinematicsGrids = "
  << *(test_input_handler.getPathToKinematicsGrids())
  << std::endl
  << "(inspect: pathToKinematicsGrids = "
  << *(test_input_handler.inspectPathToKinematicsGrids()) << ")"
  << std::endl
  << std::endl
  << "SPS1a MASS:"
  << test_input_handler.getCppSlha()->get_MASS_as_string()
  << std::endl
  << std::endl
  << *(test_input_handler.getSpectrum()->get_PDG_coded_name(
                                                        CppSLHA::PDG_code::Z ))
  << " "
  << test_input_handler.getSpectrum()->get_PDG_coded_mass(
                                                         CppSLHA::PDG_code::Z )
  << " GeV"
  << std::endl
  << *(test_input_handler.getWPlus()->getName()) << " "
  << test_input_handler.getWPlus()->get_mass() << " GeV"
  << std::endl
  << *(test_input_handler.getTop()->getName()) << " "
  << test_input_handler.getTop()->get_mass() << " GeV"
  << std::endl
  << *(test_input_handler.getGluino()->getName()) << " "
  << test_input_handler.getGluino()->get_mass() << " GeV"
  << std::endl
  << "average_squarks4_mass = "
  << test_input_handler.getAverageSquarks4Mass()
  << std::endl
  << "scoloreds: " << std::endl;
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = test_input_handler.getColoredSparticles()->begin();
       test_input_handler.getColoredSparticles()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->getName()) << std::endl;
    }
  std::cout
  << std::endl
  << "charginos: " << std::endl;
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = test_input_handler.getCharginos()->begin();
       test_input_handler.getCharginos()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->getName()) << std::endl;
    }
  std::cout
  << std::endl
  << "neutralinos: " << std::endl;
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = test_input_handler.getNeutralinos()->begin();
       test_input_handler.getNeutralinos()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->getName()) << std::endl;
    }
  std::cout
  << std::endl
  << "unstable neutralinos: " << std::endl;
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = test_input_handler.getUnstableNeutralinos()->begin();
       test_input_handler.getUnstableNeutralinos()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->getName()) << std::endl;
    }
  std::cout
  << std::endl
  << "electroweakinos: " << std::endl;
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = test_input_handler.getElectroweakinos()->begin();
       test_input_handler.getElectroweakinos()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->getName()) << std::endl;
    }
  std::cout
  << std::endl
  << "unstable electroweakinos: " << std::endl;
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator
       = test_input_handler.getUnstableElectroweakinos()->begin();
       test_input_handler.getUnstableElectroweakinos()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->getName()) << std::endl;
    }
  std::cout
  << std::endl
  << "charged sleptons: " << std::endl;
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = test_input_handler.getChargedSleptons()->begin();
       test_input_handler.getChargedSleptons()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->getName()) << std::endl;
    }
  std::cout
  << std::endl
  << "sneutrinos: " << std::endl;
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = test_input_handler.getSneutrinos()->begin();
       test_input_handler.getSneutrinos()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->getName()) << std::endl;
    }
  std::cout
  << std::endl
  << "sdowns: " << std::endl;
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = test_input_handler.getSdownTypes()->begin();
       test_input_handler.getSdownTypes()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->getName()) << std::endl;
    }
  std::cout
  << std::endl
  << "sups: " << std::endl;
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = test_input_handler.getSupTypes()->begin();
       test_input_handler.getSupTypes()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->getName()) << std::endl;
    }
  std::cout
  << std::endl
  << "jets: " << std::endl;
  for( std::list< int >::const_iterator
       set_iterator = test_input_handler.getJets()->begin();
       test_input_handler.getJets()->end() != set_iterator;
       ++set_iterator )
    {
      std::cout
      << *(test_input_handler.getSpectrum()->get_PDG_coded_name(
                                                             (*set_iterator) ))
      << std::endl;
    }
  std::cout
  << std::endl
  << "jets4: " << std::endl;
  for( std::list< int >::const_iterator
       set_iterator = test_input_handler.getJets4()->begin();
       test_input_handler.getJets4()->end() != set_iterator;
       ++set_iterator )
    {
      std::cout
      << *(test_input_handler.getSpectrum()->get_PDG_coded_name(
                                                             (*set_iterator) ))
      << std::endl;
    }
  std::cout
  << std::endl
  << "jets5: " << std::endl;
  for( std::list< int >::const_iterator
       set_iterator = test_input_handler.getJets5()->begin();
       test_input_handler.getJets5()->end() != set_iterator;
       ++set_iterator )
    {
      std::cout
      << *(test_input_handler.getSpectrum()->get_PDG_coded_name(
                                                             (*set_iterator) ))
      << std::endl;
    }
  std::cout
  << std::endl
  << "not in jets4: " << std::endl;
  for( std::list< int >::const_iterator
       set_iterator = test_input_handler.getNotInJets4()->begin();
       test_input_handler.getNotInJets4()->end() != set_iterator;
       ++set_iterator )
    {
      std::cout
      << *(test_input_handler.getSpectrum()->get_PDG_coded_name(
                                                             (*set_iterator) ))
      << std::endl;
    }
  std::cout
  << std::endl
  << "not in jets5: " << std::endl;
  for( std::list< int >::const_iterator
       set_iterator = test_input_handler.getNotInJets5()->begin();
       test_input_handler.getNotInJets5()->end() != set_iterator;
       ++set_iterator )
    {
      std::cout
      << *(test_input_handler.getSpectrum()->get_PDG_coded_name(
                                                             (*set_iterator) ))
      << std::endl;
    }
  std::cout
  << std::endl
  << "squarks: " << std::endl;
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = test_input_handler.getSquarks()->begin();
       test_input_handler.getSquarks()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->getName()) << std::endl;
    }
  std::cout
  << std::endl
  << "squarks4: " << std::endl;
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = test_input_handler.getSquarks4()->begin();
       test_input_handler.getSquarks4()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->getName()) << std::endl;
    }
  std::cout
  << std::endl
  << "squarks5: " << std::endl;
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = test_input_handler.getSquarks5()->begin();
       test_input_handler.getSquarks5()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->getName()) << std::endl;
    }
  std::cout
  << std::endl
  << "sbottom_1 is in squarks? "
  << test_input_handler.isIn(
                          test_input_handler.getSbottomOne()->get_PDG_code(),
                               test_input_handler.getSquarks() )
  << std::endl
  << "sbottom_1 is in squarks4? "
  << test_input_handler.isIn(
                          test_input_handler.getSbottomOne()->get_PDG_code(),
                               test_input_handler.getSquarks4() )
  << std::endl
  << std::endl
  << "colored cascade set: ";
  LHC_FASER::colored_cascade_set*
  gluino_cascade
  = test_input_handler.get_colored_cascade_set(
                                             test_input_handler.getGluino() );
  LHC_FASER::colored_cascade_set*
  sbottom_cascade
  = test_input_handler.get_colored_cascade_set(
                                        test_input_handler.getSbottomOne() );
  LHC_FASER::colored_cascade_set*
  gluino_cascade2
  = test_input_handler.get_colored_cascade_set(
                                             test_input_handler.getGluino() );
  std::cout
  << *(gluino_cascade->getDecayer()->getName())
  << ":" << std::endl
  << "pointer = " << gluino_cascade << ", pointer2 = " << gluino_cascade2
  << std::endl;
  for( std::vector< LHC_FASER::colored_cascade* >::const_iterator
       set_iterator = gluino_cascade->get_colored_cascades()->begin();
       gluino_cascade->get_colored_cascades()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout
      << *((*set_iterator)->getDecayer()->getName())
      << " -> " << *((*set_iterator)->getProduct()->getName())
      << " ( " << (*set_iterator)->get_product_position() << " ) -> "
      << *((*set_iterator)->get_EWino()->getName())
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
  << *(sbottom_cascade->getDecayer()->getName())
  << ":" << std::endl
  << "pointer = " << sbottom_cascade
  << std::endl;
  for( std::vector< LHC_FASER::colored_cascade* >::const_iterator
       set_iterator = sbottom_cascade->get_colored_cascades()->begin();
       sbottom_cascade->get_colored_cascades()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout
      << *((*set_iterator)->getDecayer()->getName())
      << " -> " << *((*set_iterator)->getProduct()->getName())
      << " ( " << (*set_iterator)->get_product_position() << " ) -> "
      << *((*set_iterator)->get_EWino()->getName())
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
  for( std::vector< LHC_FASER::signedParticleShortcutPair* >::const_iterator
       set_iterator
       = test_input_handler.getScoloredProductionCombinations()->begin();
       test_input_handler.getScoloredProductionCombinations()->end()
       > set_iterator;
       ++set_iterator )
    {
      if( (*set_iterator)->is_first_particle() )
        {
          std::cout << *((*set_iterator)->getFirstParticle()->getName());
        }
      else
        {
          std::cout << *((*set_iterator)->getFirstParticle()->get_antiname());
        }
      std::cout << " + ";
      if( (*set_iterator)->is_second_particle() )
        {
          std::cout << *((*set_iterator)->getSecondParticle()->getName());
        }
      else
        {
          std::cout
          << *((*set_iterator)->getSecondParticle()->get_antiname());
        }
      std::cout << std::endl;
    }
  std::vector< CppSLHA::particle_property_set const* > sdown_LR;
  sdown_LR.push_back( test_input_handler.getSdownL() );
  sdown_LR.push_back( test_input_handler.getSdownR() );
  std::vector< CppSLHA::particle_property_set const* > sup_LR;
  sup_LR.push_back( test_input_handler.getSupL() );
  sup_LR.push_back( test_input_handler.getSupR() );
  LHC_FASER::directDecayChecker*
  gluino_to_sdown_LR
  = test_input_handler.getDirectDecayChecker(
                                               test_input_handler.getGluino(),
                                                 &sdown_LR );
  LHC_FASER::directDecayChecker*
  gluino_to_sup_LR
  = test_input_handler.getDirectDecayChecker(
                                               test_input_handler.getGluino(),
                                                 &sup_LR );
  LHC_FASER::directDecayChecker*
  gluino_to_sdown_LR2
  = test_input_handler.getDirectDecayChecker(
                                               test_input_handler.getGluino(),
                                                 &sdown_LR );
  std::cout
  << std::endl
  << "directDecayChecker: gluino -> sdown_L, sdown_R"
  << std::endl << "pointer = "
  << gluino_to_sdown_LR
  << std::endl << *(gluino_to_sdown_LR->getDecayer()->getName())
  << " -> ";
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = gluino_to_sdown_LR->getProducts()->begin();
       gluino_to_sdown_LR->getProducts()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->getName()) << " ";
    }
  std::cout
  << std::endl
  << "has direct decay? " << gluino_to_sdown_LR->hasDirectDecay();
  std::cout
  << std::endl
  << "2nd gluino -> sdown_L, sdown_R"
  << std::endl << "pointer = "
  << gluino_to_sdown_LR2
  << std::endl << *(gluino_to_sdown_LR2->getDecayer()->getName())
  << " -> ";
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = gluino_to_sdown_LR2->getProducts()->begin();
       gluino_to_sdown_LR2->getProducts()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->getName()) << " ";
    }
  std::cout
  << std::endl
  << "has direct decay? " << gluino_to_sdown_LR2->hasDirectDecay();
  std::cout
  << std::endl
  << "directDecayChecker: gluino -> sup_L, sup_R"
  << std::endl << "pointer = "
  << gluino_to_sup_LR
  << std::endl << *(gluino_to_sup_LR->getDecayer()->getName())
  << " -> ";
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = gluino_to_sup_LR->getProducts()->begin();
       gluino_to_sup_LR->getProducts()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->getName()) << " ";
    }
  std::cout
  << std::endl
  << "has direct decay? " << gluino_to_sup_LR->hasDirectDecay();

  LHC_FASER::exclusiveBrCalculator*
  gluino_to_sdown_L_BR
  = test_input_handler.getExclusiveBrCalculator(
                                               test_input_handler.getGluino(),
                                              test_input_handler.getSdownL(),
                                                    true,
                                       test_input_handler.getNotInJets5() );
  LHC_FASER::exclusiveBrCalculator*
  gluino_to_sdown_R_BR
  = test_input_handler.getExclusiveBrCalculator(
                                               test_input_handler.getGluino(),
                                              test_input_handler.getSdownR(),
                                                    true,
                                       test_input_handler.getNotInJets5() );
  LHC_FASER::exclusiveBrCalculator*
  gluino_to_sdown_L_BR2
  = test_input_handler.getExclusiveBrCalculator(
                                               test_input_handler.getGluino(),
                                              test_input_handler.getSdownL(),
                                                    true,
                                       test_input_handler.getNotInJets5() );

  std::cout
  << std::endl
  << "exclusive BR checker: gluino -> sdown_L"
  << std::endl << "pointer = "
  << gluino_to_sdown_L_BR
  << std::endl << *(gluino_to_sdown_L_BR->getDecayer()->getName())
  << " -> " << *(gluino_to_sdown_L_BR->getProduct()->getName())
  << " has BR " << gluino_to_sdown_L_BR->getBr()
  << std::endl
  << "2nd gluino -> sdown_L"
  << std::endl << "pointer = "
  << gluino_to_sdown_L_BR2
  << std::endl << *(gluino_to_sdown_L_BR2->getDecayer()->getName())
  << " -> " << *(gluino_to_sdown_L_BR2->getProduct()->getName())
  << " has BR " << gluino_to_sdown_L_BR2->getBr()
  << std::endl
  << "gluino -> sdown_R"
  << std::endl << "pointer = "
  << gluino_to_sdown_R_BR
  << std::endl << *(gluino_to_sdown_R_BR->getDecayer()->getName())
  << " -> " << *(gluino_to_sdown_R_BR->getProduct()->getName())
  << " has BR " << gluino_to_sdown_R_BR->getBr()
  << std::endl
  << std::endl;


  SLHA_data.read_file( "SPS2_spectrum.out" );
  test_readier.readyObserversForNewPoint();
  std::cout
  << std::endl
  << "SPS2 MASS:"
  << test_input_handler.getCppSlha()->get_MASS_as_string()
  << std::endl
  << std::endl
  << *(test_input_handler.getGluino()->getName()) << " "
  << test_input_handler.getGluino()->get_mass() << " GeV"
  << std::endl
  << std::endl
  << "charged sleptons: " << std::endl;
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = test_input_handler.getChargedSleptons()->begin();
       test_input_handler.getChargedSleptons()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->getName()) << std::endl;
    }
  std::cout
  << *(gluino_cascade->getDecayer()->getName())
  << ":" << std::endl
  << "pointer = " << gluino_cascade << ", pointer2 = " << gluino_cascade2
  << std::endl;
  for( std::vector< LHC_FASER::colored_cascade* >::const_iterator
       set_iterator = gluino_cascade->get_colored_cascades()->begin();
       gluino_cascade->get_colored_cascades()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout
      << *((*set_iterator)->getDecayer()->getName())
      << " -> " << *((*set_iterator)->getProduct()->getName())
      << " ( " << (*set_iterator)->get_product_position() << " ) -> "
      << *((*set_iterator)->get_EWino()->getName())
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
  << *(sbottom_cascade->getDecayer()->getName())
  << ":" << std::endl
  << "pointer = " << sbottom_cascade
  << std::endl;
  for( std::vector< LHC_FASER::colored_cascade* >::const_iterator
       set_iterator = sbottom_cascade->get_colored_cascades()->begin();
       sbottom_cascade->get_colored_cascades()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout
      << *((*set_iterator)->getDecayer()->getName())
      << " -> " << *((*set_iterator)->getProduct()->getName())
      << " ( " << (*set_iterator)->get_product_position() << " ) -> "
      << *((*set_iterator)->get_EWino()->getName())
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
  << "directDecayChecker: gluino -> sdown_L, sdown_R"
  << std::endl << "pointer = "
  << gluino_to_sdown_LR
  << std::endl << *(gluino_to_sdown_LR->getDecayer()->getName())
  << " -> ";
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = gluino_to_sdown_LR->getProducts()->begin();
       gluino_to_sdown_LR->getProducts()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->getName()) << " ";
    }
  std::cout
  << std::endl
  << "has direct decay? " << gluino_to_sdown_LR->hasDirectDecay();
  std::cout
  << std::endl
  << "exclusive BR checker: gluino -> sdown_L"
  << std::endl << "pointer = "
  << gluino_to_sdown_L_BR
  << std::endl << *(gluino_to_sdown_L_BR->getDecayer()->getName())
  << " -> " << *(gluino_to_sdown_L_BR->getProduct()->getName())
  << " has BR " << gluino_to_sdown_L_BR->getBr()
  << std::endl;


  SLHA_data.read_file( "SPS1a_spectrum.out" );
  test_readier.readyObserversForNewPoint();
  std::cout
  << std::endl
  << "SPS1a MASS:"
  << test_input_handler.getCppSlha()->get_MASS_as_string()
  << std::endl
  << std::endl
  << *(test_input_handler.getGluino()->getName()) << " "
  << test_input_handler.getGluino()->get_mass() << " GeV"
  << std::endl
  << std::endl
  << "charged sleptons: " << std::endl;
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = test_input_handler.getChargedSleptons()->begin();
       test_input_handler.getChargedSleptons()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->getName()) << std::endl;
    }
  std::cout
  << *(gluino_cascade->getDecayer()->getName())
  << ":" << std::endl
  << "pointer = " << gluino_cascade << ", pointer2 = " << gluino_cascade2
  << std::endl;
  for( std::vector< LHC_FASER::colored_cascade* >::const_iterator
       set_iterator = gluino_cascade->get_colored_cascades()->begin();
       gluino_cascade->get_colored_cascades()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout
      << *((*set_iterator)->getDecayer()->getName())
      << " -> " << *((*set_iterator)->getProduct()->getName())
      << " ( " << (*set_iterator)->get_product_position() << " ) -> "
      << *((*set_iterator)->get_EWino()->getName())
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
  << *(sbottom_cascade->getDecayer()->getName())
  << ":" << std::endl
  << "pointer = " << sbottom_cascade
  << std::endl;
  for( std::vector< LHC_FASER::colored_cascade* >::const_iterator
       set_iterator = sbottom_cascade->get_colored_cascades()->begin();
       sbottom_cascade->get_colored_cascades()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout
      << *((*set_iterator)->getDecayer()->getName())
      << " -> " << *((*set_iterator)->getProduct()->getName())
      << " ( " << (*set_iterator)->get_product_position() << " ) -> "
      << *((*set_iterator)->get_EWino()->getName())
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
  << "directDecayChecker: gluino -> sdown_L, sdown_R"
  << std::endl << "pointer = "
  << gluino_to_sdown_LR
  << std::endl << *(gluino_to_sdown_LR->getDecayer()->getName())
  << " -> ";
  for( std::vector< CppSLHA::particle_property_set const* >::const_iterator
       set_iterator = gluino_to_sdown_LR->getProducts()->begin();
       gluino_to_sdown_LR->getProducts()->end() > set_iterator;
       ++set_iterator )
    {
      std::cout << *((*set_iterator)->getName()) << " ";
    }
  std::cout
  << std::endl
  << "has direct decay? " << gluino_to_sdown_LR->hasDirectDecay();
  std::cout
  << std::endl
  << "exclusive BR checker: gluino -> sdown_L"
  << std::endl << "pointer = "
  << gluino_to_sdown_L_BR
  << std::endl << *(gluino_to_sdown_L_BR->getDecayer()->getName())
  << " -> " << *(gluino_to_sdown_L_BR->getProduct()->getName())
  << " has BR " << gluino_to_sdown_L_BR->getBr()
  << std::endl;

  SLHA_data.read_file( "SPS2_spectrum.out" );
  test_readier.readyObserversForNewPoint();
  std::cout
  << std::endl
  << "exclusive BR checker: gluino -> sdown_L"
  << std::endl << "pointer = "
  << gluino_to_sdown_L_BR
  << std::endl << *(gluino_to_sdown_L_BR->getDecayer()->getName())
  << " -> " << *(gluino_to_sdown_L_BR->getProduct()->getName())
  << " has BR " << gluino_to_sdown_L_BR->getBr()
  << std::endl;
  std::cout
  << std::endl
  << "exclusive BR checker: gluino -> sdown_L"
  << std::endl << "pointer = "
  << gluino_to_sdown_L_BR
  << std::endl << *(gluino_to_sdown_L_BR->getDecayer()->getName())
  << " -> " << *(gluino_to_sdown_L_BR->getProduct()->getName())
  << " has BR " << gluino_to_sdown_L_BR->getBr()
  << std::endl;
  SLHA_data.read_file( "SPS1a_spectrum.out" );
  test_readier.readyObserversForNewPoint();
  std::cout
  << std::endl
  << "exclusive BR checker: gluino -> sdown_L"
  << std::endl << "pointer = "
  << gluino_to_sdown_L_BR
  << std::endl << *(gluino_to_sdown_L_BR->getDecayer()->getName())
  << " -> " << *(gluino_to_sdown_L_BR->getProduct()->getName())
  << " has BR " << gluino_to_sdown_L_BR->getBr()
  << std::endl;
  std::cout
  << std::endl
  << "exclusive BR checker: gluino -> sdown_L"
  << std::endl << "pointer = "
  << gluino_to_sdown_L_BR
  << std::endl << *(gluino_to_sdown_L_BR->getDecayer()->getName())
  << " -> " << *(gluino_to_sdown_L_BR->getProduct()->getName())
  << " has BR " << gluino_to_sdown_L_BR->getBr()
  << std::endl
  << std::endl
  << std::endl
  << "gluino_to_sdown_LR = " << gluino_to_sdown_LR
  << std::endl
  << "gluino pointer = " << test_input_handler.getGluino()
  << std::endl
  << "sdown_L pointer = " << test_input_handler.getSdownL()
  << std::endl
  << "sdown_R pointer = " << test_input_handler.getSdownR()
  << std::endl;

  LHC_FASER::readierForNewPoint*
  test_readier2 = new LHC_FASER::readierForNewPoint();
  LHC_FASER::exclusiveBrCalculator*
  sdown_L_chargino_one_minus
  = new LHC_FASER::exclusiveBrCalculator( test_input_handler.getSdownL(),
                                         test_input_handler.getCharginoOne(),
                                            false,
                                         test_input_handler.getNotInJets4(),
                                            test_readier2 );
  LHC_FASER::exclusiveBrCalculator*
  sdown_L_chargino_one_plus
  = new LHC_FASER::exclusiveBrCalculator( test_input_handler.getSdownL(),
                                         test_input_handler.getCharginoOne(),
                                            true,
                                         test_input_handler.getNotInJets4(),
                                            test_readier2 );
  std::cout
  << std::endl << "sdown_L_chargino_one_minus BR = "
  << sdown_L_chargino_one_minus->getBr()
  << std::endl << "sdown_L_chargino_one_plus BR = "
  << sdown_L_chargino_one_plus->getBr();
  test_readier2->readyObserversForNewPoint();
  delete sdown_L_chargino_one_plus;
  std::cout
  << std::endl << "sdown_L_chargino_one_minus BR = "
  << sdown_L_chargino_one_minus->getBr();
  test_readier2->readyObserversForNewPoint();
  delete test_readier2;
  std::cout
  << std::endl << "sdown_L_chargino_one_minus BR = "
  << sdown_L_chargino_one_minus->getBr();
  delete sdown_L_chargino_one_minus;
  **/


  /* testing the cross-section grid stuff: (should work, may need slight changes of input)*//**

  std::cout
  << std::endl
  << "pathToGrids = "
  << *(test_input_handler.getPathToGrids())
  << std::endl;

  std::string
  path_to_test_grid( *(test_input_handler.getPathToCrossSectionGrids()) );
  path_to_test_grid.append( "/7TeV/gluino+gluino_LO.dat" );
  LHC_FASER::squareGrid test_grid( &path_to_test_grid,
                                    &bluh,
                                    NULL );
  LHC_FASER::crossSectionTable gg_table( &test_grid,
           test_input_handler.getScoloredProductionCombinations()->at( 0 ),
                                           1.0,
                                           &test_input_handler );
  std::cout
  << std::endl << "x: " << test_grid.getLowestX() << " - "
  << test_grid.getHighestX();
  std::cout
  << std::endl << "y: " << test_grid.getLowestY() << " - "
  << test_grid.getHighestY();
  std::cout
  << std::endl << "value at { 200.0, 200.0 } = " << test_grid.valueAt( 200.0,
                                                                       200.0 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at { 200.0, 205.0 } = " << test_grid.valueAt( 200.0,
                                                                       205.0 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at { 205.0, 200.0 } = " << test_grid.valueAt( 205.0,
                                                                       200.0 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at { 201.0, 201.0 } = " << test_grid.valueAt( 201.0,
                                                                       201.0 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at { 210.0, 210.0 } = " << test_grid.valueAt( 210.0,
                                                                       210.0 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at { 199.5, 201.0 } = " << test_grid.valueAt( 199.5,
                                                                       201.0 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at { 199.5, 199.5 } = " << test_grid.valueAt( 199.5,
                                                                       199.5 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at { 2000.0, 200.5 } = "
  << test_grid.valueAt( 2000.0,
                         200.5 );
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout
  << std::endl << "value at { 2000.5, 200.5 } = "
  << test_grid.valueAt( 2000.5,
                         200.5 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at { 200.5, 2000.0 } = "
  << test_grid.valueAt( 200.5,
                         2000.0 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at { 200.5, 2000.5 } = "
  << test_grid.valueAt( 200.5,
                         2000.5 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at { 1999.5, 1999.5 } = "
  << test_grid.valueAt( 1999.5,
                         1999.5 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at { 2000.0, 2000.0 } = "
  << test_grid.valueAt( 2000.0,
                         2000.0 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at { 2000.5, 2000.5 } = "
  << test_grid.valueAt( 2000.5,
                         2000.5 );
  std::cout << std::endl;
  std::cout
  << std::endl << "ggtable value = "
  << gg_table.getValue();

  path_to_test_grid.assign(
                     *(test_input_handler.getPathToCrossSectionGrids()) );
  path_to_test_grid.append( "/7TeV/gluino+gluino_K.dat" );
  LHC_FASER::squareGrid test_K( &path_to_test_grid,
                                 &bluh,
                                 NULL );
  std::cout
  << std::endl << "K value at { 200.5, 200.5 } = "
  << test_K.valueAt( 200.5,
                      200.5 );
  std::cout << std::endl;


  LHC_FASER::crossSectionHandler cross_section_handler_tester(
                                                         &test_input_handler );
  LHC_FASER::crossSectionTableSet*
  test_cross_section_set = cross_section_handler_tester.getTableSet( 7 );
  std::cout
  << std::endl << "test_cross_section_set = " << test_cross_section_set;
  N.B.: for the next bit to work, crossSectionHandler::getTable() has to be
  moved back to public.
  LHC_FASER::crossSectionTable*
  sd_go_sigma
  = test_cross_section_set->getTable(
         test_input_handler.getScoloredProductionCombinations()->at( 1 ) );
  LHC_FASER::crossSectionTable*
  su_go_sigma
  = test_cross_section_set->getTable(
        test_input_handler.getScoloredProductionCombinations()->at( 10 ) );
  LHC_FASER::crossSectionTableSet*
  test_cross_section_set2 = cross_section_handler_tester.getTableSet( 7 );
  std::cout
  << std::endl << "test_cross_section_set2 = " << test_cross_section_set2;
  LHC_FASER::crossSectionTable*
  sd_go_sigma2
  = test_cross_section_set2->getTable(
         test_input_handler.getScoloredProductionCombinations()->at( 1 ) );
  std::cout
  << std::endl << "no extra squareGrid loading?";

  path_to_test_grid.assign(
                     *(test_input_handler.getPathToCrossSectionGrids()) );
  path_to_test_grid.append( "/7TeV/sdown_L+gluino_LO.dat" );
  LHC_FASER::squareGrid test_sd_go_LO( &path_to_test_grid,
                                        &bluh,
                                        NULL );
  path_to_test_grid.assign(
                     *(test_input_handler.getPathToCrossSectionGrids()) );
  path_to_test_grid.append( "/7TeV/sdown_L+gluino_K.dat" );
  LHC_FASER::squareGrid test_sd_go_K( &path_to_test_grid,
                                       &bluh,
                                       NULL );
  double
  sd_mass
  = test_input_handler.getScoloredProductionCombinations()->at( 1
                                   )->getFirstParticle()->get_absolute_mass();
  double
  go_mass
  = test_input_handler.getScoloredProductionCombinations()->at( 1
                                  )->getSecondParticle()->get_absolute_mass();
  std::cout << std::endl;
  std::cout
  << std::endl << "sd_mass = " << sd_mass
  << std::endl << "go_mass = " << go_mass
  << std::endl << "test_sd_go_LO value = " << test_sd_go_LO.valueAt( sd_mass,
                                                                      go_mass )
  << std::endl << "test_sd_go_K value = " << test_sd_go_K.valueAt( sd_mass,
                                                                    go_mass )
  << std::endl << "hence NLO value should be "
  << ( test_sd_go_LO.valueAt( sd_mass,
                               go_mass )
       * test_sd_go_K.valueAt( sd_mass,
                                go_mass ) )
  << std::endl << "sd_go_sigma = " << sd_go_sigma
  << std::endl << "sd_go_sigma value = "
  << sd_go_sigma->getValue()
  << std::endl << "sd_go_sigma2 = " << sd_go_sigma2
  << std::endl << "sd_go_sigma2 value = "
  << sd_go_sigma->getValue()
  << std::endl << "su_go_sigma = " << su_go_sigma
  << std::endl << "su_go_sigma value = "
  << su_go_sigma->getValue();

  **/


  /* testing the kinematics grid stuff: (ahhhh, tssss, might still work involes a lot of old classes that are not yet deleted)*//**

  std::string
  acceptance_grid_path( *(test_input_handler.getPathToKinematicsGrids()) );
  std::cout
  << std::endl << "path to kinematics grids = " << acceptance_grid_path;
  std::cout << std::endl;
  acceptance_grid_path.append(
                        "/7TeV/jets/Atlas4jMET/squark+gluino_acceptance.dat" );
  LHC_FASER::acceptanceGrid*
  jets_test_acceptance
  = new LHC_FASER::acceptanceGrid( &acceptance_grid_path );
  std::cout
  << std::endl << "jets_test_acceptance->valueAt( "
  << jets_test_acceptance->getLowestSquarkMass() << ", "
  << jets_test_acceptance->getLowestGluinoMass() << ", 0.0, 0.0, 0 ) = "
  << jets_test_acceptance->valueAt(
                                jets_test_acceptance->getLowestSquarkMass(),
                                jets_test_acceptance->getLowestGluinoMass(),
                                0.0,
                                0.0,
                                0,
                                false )
  << std::endl
  << "jets_test_acceptance->valueAt( 123.4, 567.8, 12.3, 45.6, 4 ) = "
  << jets_test_acceptance->valueAt( 123.4,
                                     567.8,
                                     12.3,
                                     45.6,
                                     4,
                                     false )
  << std::endl
  << "jets_test_acceptance->valueAt( 234.5, 678.9, 12.3, 45.6, 4 ) = "
  << jets_test_acceptance->valueAt( 234.5,
                                     678.9,
                                     12.3,
                                     45.6,
                                     4,
                                     false )
  << std::endl
  << "jets_test_acceptance->valueAt( 876.5, 432.1, 1.0, 1.0, 4 ) = "
  << jets_test_acceptance->valueAt( 876.5,
                                     432.1,
                                     1.0,
                                     1.0,
                                     4,
                                     false )
  << std::endl
  << "jets_test_acceptance->valueAt( 501.0, 501.0, 0.0, 0.0, 4 ) = "
  << jets_test_acceptance->valueAt( 501.0,
                                     501.0,
                                     0.0,
                                     0.0,
                                     4,
                                     false );
  std::cout << std::endl;
  acceptance_grid_path.assign(
                        *(test_input_handler.getPathToKinematicsGrids()) );
  std::cout
  << std::endl << "path to kinematics grids = " << acceptance_grid_path;
  std::cout << std::endl;
  acceptance_grid_path.append( "/7TeV/leptons/squark+gluino_acceptance.dat" );
  LHC_FASER::acceptanceGrid*
  leptons_test_acceptance
  = new LHC_FASER::acceptanceGrid( &acceptance_grid_path );
  std::cout
  << std::endl << "leptons_test_acceptance->valueAt( "
  << leptons_test_acceptance->getLowestSquarkMass() << ", "
  << leptons_test_acceptance->getLowestGluinoMass() << ", 0.0, 0.0, 0 ) = "
  << leptons_test_acceptance->valueAt(
                             leptons_test_acceptance->getLowestSquarkMass(),
                             leptons_test_acceptance->getLowestGluinoMass(),
                                        0.0,
                                        0.0,
                                        0,
                                        true )
  << std::endl
  << "leptons_test_acceptance->valueAt( 123.4, 567.8, 12.3, 45.6, 4 ) = "
  << leptons_test_acceptance->valueAt( 123.4,
                                        567.8,
                                        12.3,
                                        45.6,
                                        4,
                                        false )
  << std::endl
  << "leptons_test_acceptance->valueAt( 234.5, 678.9, 12.3, 45.6, 4 ) = "
  << leptons_test_acceptance->valueAt( 234.5,
                                        678.9,
                                        12.3,
                                        45.6,
                                        4,
                                        false )
  << std::endl
  << "leptons_test_acceptance->valueAt( 876.5, 432.1, 1.0, 1.0, 4 ) = "
  << leptons_test_acceptance->valueAt( 876.5,
                                        432.1,
                                        1.0,
                                        1.0,
                                        4,
                                        false )
  << std::endl
  << "leptons_test_acceptance->valueAt( 501.0, 501.0, 0.0, 0.0, 4 ) = "
  << leptons_test_acceptance->valueAt( 501.0,
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
                                             test_input_handler.getGluino() );
  LHC_FASER::colored_cascade_set*
  sdown_L_cascade
  = test_input_handler.get_colored_cascade_set(
                                            test_input_handler.getSdownL() );
  int const sdL_go = 1;
  acceptance_grid_path.assign(
                        *(test_input_handler.getPathToKinematicsGrids()) );
  acceptance_grid_path.append( "/7TeV/jets/Atlas4jMET" );
  std::cout
  << std::endl << "path to jets kinematics grids = " << acceptance_grid_path;
  std::cout << std::endl;
  LHC_FASER::jet_acceptance_table test_jet_kinematics_table(
                                                         &acceptance_grid_path,
                                                             &bluh,
                                                         &test_input_handler );
  std::cout
  << std::endl << "test_jet_kinematics_table.getName() = "
  << *(test_jet_kinematics_table.getName())
  << std::endl << "test_jet_kinematics_table.get_acceptance( 4,"
  << std::endl
  << "[ ";
  if( test_input_handler.getScoloredProductionCombinations()->at( sdL_go
                                                       )->is_first_particle() )
    {
      std::cout
      << *(test_input_handler.getScoloredProductionCombinations(
                             )->at( sdL_go )->getFirstParticle()->getName());
    }
  else
    {
      std::cout
      << *(test_input_handler.getScoloredProductionCombinations()->at(
                                sdL_go )->getFirstParticle()->get_antiname());
    }
  std::cout << " + ";
  if( test_input_handler.getScoloredProductionCombinations()->at( sdL_go
                                                      )->is_second_particle() )
    {
      std::cout
      << *(test_input_handler.getScoloredProductionCombinations()->at(
                                   sdL_go )->getSecondParticle()->getName());
    }
  else
    {
      std::cout
      << *(test_input_handler.getScoloredProductionCombinations()->at(
                               sdL_go )->getSecondParticle()->get_antiname());
    }
  std::cout << " ]," << std::endl;
  std::cout
  << "[ "
  << *(gluino_cascade3->get_colored_cascades()->at( 1
                                                  )->getDecayer()->getName())
  << " -> " << *(gluino_cascade3->get_colored_cascades()->at( 1
                                                  )->getProduct()->getName())
  << " ( " << gluino_cascade3->get_colored_cascades()->at( 1
                                          )->get_product_position() << " ) -> "
  << *(gluino_cascade3->get_colored_cascades()->at( 1
                                                    )->get_EWino()->getName())
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
                                                  )->getDecayer()->getName())
  << " -> " << *(sdown_L_cascade->get_colored_cascades()->at( 1
                                                  )->getProduct()->getName())
  << " ( " << sdown_L_cascade->get_colored_cascades()->at( 1
                                          )->get_product_position() << " ) -> "
  << *(sdown_L_cascade->get_colored_cascades()->at( 1
                                                    )->get_EWino()->getName())
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
      test_input_handler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
          sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance();
  std::cout << std::endl;

  acceptance_grid_path.assign(
                        *(test_input_handler.getPathToKinematicsGrids()) );
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
  if( test_input_handler.getScoloredProductionCombinations()->at( sdL_go
                                                       )->is_first_particle() )
    {
      std::cout
      << *(test_input_handler.getScoloredProductionCombinations()->at(
                                    sdL_go )->getFirstParticle()->getName());
    }
  else
    {
      std::cout
      << *(test_input_handler.getScoloredProductionCombinations()->at(
                                sdL_go )->getFirstParticle()->get_antiname());
    }
  std::cout << " + ";
  if( test_input_handler.getScoloredProductionCombinations()->at( sdL_go
                                                      )->is_second_particle() )
    {
      std::cout
      << *(test_input_handler.getScoloredProductionCombinations()->at(
                                   sdL_go )->getSecondParticle()->getName());
    }
  else
    {
      std::cout
      << *(test_input_handler.getScoloredProductionCombinations()->at(
                               sdL_go )->getSecondParticle()->get_antiname());
    }
  std::cout << " ]," << std::endl;
  std::cout
  << "[ "
  << *(gluino_cascade3->get_colored_cascades()->at( 1
                                                  )->getDecayer()->getName())
  << " -> " << *(gluino_cascade3->get_colored_cascades()->at( 1
                                                  )->getProduct()->getName())
  << " ( " << gluino_cascade3->get_colored_cascades()->at( 1
                                          )->get_product_position() << " ) -> "
  << *(gluino_cascade3->get_colored_cascades()->at( 1
                                                    )->get_EWino()->getName())
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
                                                  )->getDecayer()->getName())
  << " -> " << *(sdown_L_cascade->get_colored_cascades()->at( 1
                                                  )->getProduct()->getName())
  << " ( " << sdown_L_cascade->get_colored_cascades()->at( 1
                                          )->get_product_position() << " ) -> "
  << *(sdown_L_cascade->get_colored_cascades()->at( 1
                                                    )->get_EWino()->getName())
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
      test_input_handler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                        )->acceptance_at( 23.4,
                                                                        20.0 );
  std::cout << std::endl;
  std::cout
  << "...->get_acceptance( -1.0, 10.0 )"
  << std::endl << " = "
  << test_lepton_kinematics_table.get_acceptance(
      test_input_handler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                        )->acceptance_at( -1.0,
                                                                        10.0 );
  std::cout << std::endl;
  std::cout
  << "...->get_acceptance( 1.0, 10.0 )"
  << std::endl << " = "
  << test_lepton_kinematics_table.get_acceptance(
      test_input_handler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                        )->acceptance_at( 1.0,
                                                                        10.0 );
  std::cout << std::endl;
  std::cout
  << "...->get_acceptance( 10.0, 10.0 )"
  << std::endl << " = "
  << test_lepton_kinematics_table.get_acceptance(
      test_input_handler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                        )->acceptance_at( 10.0,
                                                                        10.0 );
  std::cout << std::endl;
  std::cout
  << "...->get_acceptance( 12.3, 10.0 )"
  << std::endl << " = "
  << test_lepton_kinematics_table.get_acceptance(
      test_input_handler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                        )->acceptance_at( 12.3,
                                                                        10.0 );
  std::cout << std::endl;
  std::cout
  << "...->get_acceptance( 20.0, 10.0 )"
  << std::endl << " = "
  << test_lepton_kinematics_table.get_acceptance(
      test_input_handler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                        )->acceptance_at( 20.0,
                                                                        10.0 );
  std::cout << std::endl;
  std::cout
  << "...->get_acceptance( 79.9, 10.0 )"
  << std::endl << " = "
  << test_lepton_kinematics_table.get_acceptance(
      test_input_handler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                        )->acceptance_at( 79.9,
                                                                        10.0 );
  std::cout << std::endl;
  std::cout
  << "...->get_acceptance( 80.0, 10.0 )"
  << std::endl << " = "
  << test_lepton_kinematics_table.get_acceptance(
      test_input_handler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                        )->acceptance_at( 80.0,
                                                                        10.0 );
  std::cout << std::endl;
  std::cout
  << "...->get_acceptance( 100.0, 10.0 )"
  << std::endl << " = "
  << test_lepton_kinematics_table.get_acceptance(
      test_input_handler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                       )->acceptance_at( 100.0,
                                                                        10.0 );
  std::cout << std::endl;
  std::cout
  << std::endl << "test_lepton_kinematics_table.get_acceptance( "
  << "[ ";
  if( test_input_handler.getScoloredProductionCombinations()->at( sdL_go
                                                       )->is_first_particle() )
    {
      std::cout
      << *(test_input_handler.getScoloredProductionCombinations()->at(
                                    sdL_go )->getFirstParticle()->getName());
    }
  else
    {
      std::cout
      << *(test_input_handler.getScoloredProductionCombinations()->at(
                                sdL_go )->getFirstParticle()->get_antiname());
    }
  std::cout << " + ";
  if( test_input_handler.getScoloredProductionCombinations()->at( sdL_go
                                                      )->is_second_particle() )
    {
      std::cout
      << *(test_input_handler.getScoloredProductionCombinations()->at(
                                   sdL_go )->getSecondParticle()->getName());
    }
  else
    {
      std::cout
      << *(test_input_handler.getScoloredProductionCombinations()->at(
                               sdL_go )->getSecondParticle()->get_antiname());
    }
  std::cout << " ]," << std::endl;
  std::cout
  << "[ "
  << *(gluino_cascade3->get_colored_cascades()->back(
                                                  )->getDecayer()->getName())
  << " -> " << *(gluino_cascade3->get_colored_cascades()->back(
                                                  )->getProduct()->getName())
  << " ( " << gluino_cascade3->get_colored_cascades()->back(
                                          )->get_product_position() << " ) -> "
  << *(gluino_cascade3->get_colored_cascades()->back(
                                                    )->get_EWino()->getName())
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
                                                  )->getDecayer()->getName())
  << " -> " << *(sdown_L_cascade->get_colored_cascades()->at( 1
                                                  )->getProduct()->getName())
  << " ( " << sdown_L_cascade->get_colored_cascades()->at( 1
                                          )->get_product_position() << " ) -> "
  << *(sdown_L_cascade->get_colored_cascades()->at( 1
                                                    )->get_EWino()->getName())
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
      test_input_handler.getScoloredProductionCombinations()->at( sdL_go ),
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
  if( test_input_handler.getScoloredProductionCombinations()->at( sdL_go
                                                       )->is_first_particle() )
    {
      std::cout
      << *(test_input_handler.getScoloredProductionCombinations()->at(
                                    sdL_go )->getFirstParticle()->getName());
    }
  else
    {
      std::cout
      << *(test_input_handler.getScoloredProductionCombinations()->at(
                                sdL_go )->getFirstParticle()->get_antiname());
    }
  std::cout << " + ";
  if( test_input_handler.getScoloredProductionCombinations()->at( sdL_go
                                                      )->is_second_particle() )
    {
      std::cout
      << *(test_input_handler.getScoloredProductionCombinations()->at(
                                   sdL_go )->getSecondParticle()->getName());
    }
  else
    {
      std::cout
      << *(test_input_handler.getScoloredProductionCombinations()->at(
                               sdL_go )->getSecondParticle()->get_antiname());
    }
  std::cout << " ]," << std::endl;
  std::cout
  << "[ "
  << *(gluino_cascade3->get_colored_cascades()->at( 1
                                                  )->getDecayer()->getName())
  << " -> " << *(gluino_cascade3->get_colored_cascades()->at( 1
                                                  )->getProduct()->getName())
  << " ( " << gluino_cascade3->get_colored_cascades()->at( 1
                                          )->get_product_position() << " ) -> "
  << *(gluino_cascade3->get_colored_cascades()->at( 1
                                                    )->get_EWino()->getName())
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
                                                  )->getDecayer()->getName())
  << " -> " << *(sdown_L_cascade->get_colored_cascades()->at( 1
                                                  )->getProduct()->getName())
  << " ( " << sdown_L_cascade->get_colored_cascades()->at( 1
                                          )->get_product_position() << " ) -> "
  << *(sdown_L_cascade->get_colored_cascades()->at( 1
                                                    )->get_EWino()->getName())
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
      test_input_handler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
          sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance();
  std::cout << std::endl;
  std::cout
  << "kinematics_table_tester.get_lepton_acceptance(...)->get_acceptance("
  << " 12.34, 8.0 )"
  << std::endl << " = "
  << kinematics_table_tester.get_lepton_acceptance(
      test_input_handler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                       )->acceptance_at( 12.34,
                                                                         8.0 );

  std::string
  kinematics_at_energy_name(
                        *(test_input_handler.getPathToKinematicsGrids()) );
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
  << "kinematics_table_set_tester.getTable( 4"
  << " )->get_jet_acceptance(...)->get_acceptance() = "
  << kinematics_table_set_tester.getTable( 4 )->get_jet_acceptance(
      test_input_handler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
          sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance();
  std::cout << std::endl;
  std::cout
  << "kinematics_table_set_tester.getTable( 4"
  << " )->get_lepton_acceptance(...)->get_acceptance( 12.34, 8.0 ) = "
  << kinematics_table_set_tester.getTable( 4 )->get_lepton_acceptance(
      test_input_handler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                       )->acceptance_at( 12.34,
                                                                         8.0 );
  std::cout << std::endl;

  LHC_FASER::kinematics_handler kinematics_handler_tester(
                                                         &test_input_handler );
  LHC_FASER::kinematics_table*
  kinematics_handler_table_one = kinematics_handler_tester.getTable( 7,
                                                              &jet_signal_name,
                                                                      4 );
  LHC_FASER::kinematics_table*
  kinematics_handler_table_two = kinematics_handler_tester.getTable( 7,
                                                              &jet_signal_name,
                                                                      5 );
  LHC_FASER::kinematics_table*
  kinematics_handler_table_three = kinematics_handler_tester.getTable( 7,
                                                              &jet_signal_name,
                                                                        4 );
  std::cout << std::endl;
  std::cout
  << std::endl
  << "kinematics_handler_table_one = " << kinematics_handler_table_one
  << std::endl
  << "->get_jet_acceptance(...)->get_acceptance() = "
  << kinematics_handler_table_one->get_jet_acceptance(
      test_input_handler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
           sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance()
  << std::endl
  << "->get_lepton_acceptance(...)->get_acceptance( 12.34, 8.0 ) = "
  << kinematics_handler_table_one->get_lepton_acceptance(
      test_input_handler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                       )->acceptance_at( 12.34,
                                                                         8.0 )
  << std::endl
  << "kinematics_handler_table_two = " << kinematics_handler_table_two
  << std::endl
  << "->get_jet_acceptance(...)->get_acceptance() = "
  << kinematics_handler_table_two->get_jet_acceptance(
      test_input_handler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
           sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance()
  << std::endl
  << "->get_lepton_acceptance(...)->get_acceptance( 12.34, 8.0 ) = "
  << kinematics_handler_table_two->get_lepton_acceptance(
      test_input_handler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                       )->acceptance_at( 12.34,
                                                                         8.0 )
  << std::endl
  << "kinematics_handler_table_three = " << kinematics_handler_table_three
  << std::endl
  << "->get_jet_acceptance(...)->get_acceptance() = "
  << kinematics_handler_table_three->get_jet_acceptance(
      test_input_handler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
           sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance()
  << std::endl
  << "->get_lepton_acceptance(...)->get_acceptance( 12.34, 8.0 ) = "
  << kinematics_handler_table_three->get_lepton_acceptance(
      test_input_handler.getScoloredProductionCombinations()->at( sdL_go ),
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
  = test_kinematics.getTable( 7,
                               &test_signal_name,
                               4 )->get_lepton_acceptance(
           test_input_handler.getScoloredProductionCombinations()->at( 1 ),
                               test_input_handler.get_colored_cascade_set(
            test_input_handler.getGluino() )->get_colored_cascades()->at( 1 ),
                               test_input_handler.get_colored_cascade_set(
         test_input_handler.getSdownL() )->get_colored_cascades()->at( 1 ) );

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
                                              test_input_handler.getSdownL(),
                                                                   true,
                                      test_input_handler.getNeutralinoFour(),
                                                                   true,
                                              test_input_handler.getSmuonR(),
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
                                              test_input_handler.getSdownL(),
                                                                  true,
                                      test_input_handler.getNeutralinoFour(),
                                                                  true,
                                             test_input_handler.getStauOne(),
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
                                              test_input_handler.getSdownL(),
                                                            true,
                                      test_input_handler.getNeutralinoFour(),
                                                            true,
                                                    test_input_handler.getZ(),
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
  = test_input_handler.getExclusiveBrCalculator(
                                      test_input_handler.getNeutralinoFour(),
                                              test_input_handler.getSmuonR(),
                                              true,
                               test_input_handler.getEmptyList() )->getBr();
  double
  smuon_R_to_neutralino_one
  = test_input_handler.getExclusiveBrCalculator(
                                              test_input_handler.getSmuonR(),
                                       test_input_handler.getNeutralinoOne(),
                                                    true,
                               test_input_handler.getEmptyList() )->getBr();
  double
  neutralino_four_to_stau_one
  = test_input_handler.getExclusiveBrCalculator(
                                      test_input_handler.getNeutralinoFour(),
                                             test_input_handler.getStauOne(),
                                                    true,
                               test_input_handler.getEmptyList() )->getBr();
  double
  stau_one_to_neutralino_one
  = test_input_handler.getExclusiveBrCalculator(
                                             test_input_handler.getStauOne(),
                                       test_input_handler.getNeutralinoOne(),
                                                    true,
                               test_input_handler.getEmptyList() )->getBr();
  double
  neutralino_four_to_Z
  = test_input_handler.getNeutralinoFour()->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                    test_input_handler.getZ()->get_PDG_code(),
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
  << "OSSF_minus_OSDF_leptons.getValue() = "
  << OSSF_minus_OSDF_leptons.getValue();
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
  << "no_jets_one_muon.getValue() = "
  << no_jets_one_muon.getValue();
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
  << ( two_jets_no_leptons.getValue()
       + one_jet_one_muon.getValue()
       + one_jet_one_electron.getValue()
       + no_jets_two_muons.getValue()
       + no_jets_two_electrons.getValue()
       + no_jets_one_muon_one_electron.getValue()
       + one_jet_no_leptons.getValue()
       + no_jets_one_muon.getValue()
       + no_jets_one_electron.getValue()
       + no_jets_no_leptons.getValue() )
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
                                              test_input_handler.getSdownL(),
                                                                   true,
                                      test_input_handler.getNeutralinoFour(),
                                                                   true,
                                              test_input_handler.getSmuonR(),
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
                                              test_input_handler.getSdownL(),
                                                                    true,
                                      test_input_handler.getNeutralinoFour(),
                                                                    true,
                                             test_input_handler.getStauOne(),
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
                                              test_input_handler.getSdownL(),
                                                            true,
                                      test_input_handler.getNeutralinoFour(),
                                                            true,
                                                    test_input_handler.getZ(),
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
  << "one_jet_one_muon2.getValue() = "
  << one_jet_one_muon2.getValue();
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
  << "OSSF_minus_OSDF_leptons2.getValue() = "
  << OSSF_minus_OSDF_leptons2.getValue();
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
  << ( two_jets_no_leptons2.getValue()
       + one_jet_one_muon2.getValue()
       + one_jet_one_electron2.getValue()
       + no_jets_two_muons2.getValue()
       + no_jets_two_electrons2.getValue()
       + no_jets_one_muon_one_electron2.getValue()
       + one_jet_no_leptons2.getValue()
       + no_jets_one_muon2.getValue()
       + no_jets_one_electron2.getValue()
       + no_jets_no_leptons2.getValue() )
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
                                              test_input_handler.getSdownL(),
                                                          true,
                                      test_input_handler.getNeutralinoFour(),
                                                          true );
  LHC_FASER::charge_summed_cascade_calculator*
  charge_summed_test_20_10_40
  = cascade_handler_tester->get_charge_summed_calculator( 20.0,
                                                          10.0,
                                                          40.0,
                                                  test_lepton_acceptance_value,
                                              test_input_handler.getSdownL(),
                                                          true,
                                      test_input_handler.getNeutralinoFour(),
                                                          true );
  LHC_FASER::charge_summed_cascade_calculator*
  charge_summed_test_10_10_40_b
  = cascade_handler_tester->get_charge_summed_calculator( 10.0,
                                                          10.0,
                                                          40.0,
                                                  test_lepton_acceptance_value,
                                              test_input_handler.getSdownL(),
                                                          true,
                                      test_input_handler.getNeutralinoFour(),
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
                                             test_input_handler.getGluino() );
  double gluino_BR_sum = 0.0;
  LHC_FASER::exclusiveBrCalculator* BR_sum_helper;
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
          = test_input_handler.getExclusiveBrCalculator(
                                               test_input_handler.getGluino(),
                                                   (*BR_iterator)->get_EWino(),
                                                            true,
                                         test_input_handler.getEmptyList() );

          if( NULL != BR_sum_helper )
            {

              gluino_BR_sum += BR_sum_helper->getBr();

            }

          if( !((*BR_iterator)->get_EWino()->counts_as_self_conjugate()) )
            {

              BR_sum_helper
              = test_input_handler.getExclusiveBrCalculator(
                                               test_input_handler.getGluino(),
                                                   (*BR_iterator)->get_EWino(),
                                                                false,
                                         test_input_handler.getEmptyList() );

              if( NULL != BR_sum_helper )
                {

                  gluino_BR_sum += BR_sum_helper->getBr();

                }

            }

        }
      else if( (*BR_iterator)->is_gsx() )
        // otherwise if it's a gluino 2-stage decay to the electroweakino...
        {

          BR_sum_helper
          = test_input_handler.getExclusiveBrCalculator(
                                               test_input_handler.getGluino(),
                                                 (*BR_iterator)->getProduct(),
                                                            true,
                                         test_input_handler.getEmptyList() );
          if( NULL != BR_sum_helper )
            {

              BR_sum_helper_double = BR_sum_helper->getBr();

              BR_sum_helper
              = test_input_handler.getExclusiveBrCalculator(
                                                 (*BR_iterator)->getProduct(),
                                                   (*BR_iterator)->get_EWino(),
                                                                true,
                                         test_input_handler.getEmptyList() );

              if( NULL != BR_sum_helper )
                {

                  gluino_BR_sum
                  += ( 2.0 * BR_sum_helper_double * BR_sum_helper->getBr() );

                }

              if( !((*BR_iterator)->get_EWino()->counts_as_self_conjugate()) )
                {

                  BR_sum_helper
                  = test_input_handler.getExclusiveBrCalculator(
                                                 (*BR_iterator)->getProduct(),
                                                   (*BR_iterator)->get_EWino(),
                                                                    false,
                                         test_input_handler.getEmptyList() );

                  if( NULL != BR_sum_helper )
                    {

                      gluino_BR_sum += ( 2.0 * BR_sum_helper_double
                                             * BR_sum_helper->getBr() );

                    }

                }

            }

        }

    }
  std::cout
  << std::endl
  << "gluino_BR_sum = " << gluino_BR_sum;
  std::cout << std::endl;



  LHC_FASER::crossSectionHandler test_cross_section_handler(
                                                         &test_input_handler );

  LHC_FASER::kinematics_handler test_kinematics_handler( &test_input_handler );

  LHC_FASER::cascade_handler test_cascade_handler( &test_input_handler );
  LHC_FASER::signalShortcuts test_signal_shortcuts( &test_input_handler,
                                                   &test_cross_section_handler,
                                                     &test_kinematics_handler,
                                                     &test_cascade_handler );

  LHC_FASER::signal_handler test_bad_signal( "testBadSignal",
                                             &test_signal_shortcuts );

  std::cout
  << std::endl
  << *(test_bad_signal.getName()) << " has signal value "
  << test_bad_signal.getValue();
  std::cout << std::endl;

  LHC_FASER::signal_handler
  cross_section_breakdown_test( "sigmaBreakdownTest",
                                &test_signal_shortcuts );

  std::cout
  << std::endl
  << *(cross_section_breakdown_test.getName()) << " has signal value "
  << cross_section_breakdown_test.getValue();
  std::cout << std::endl;

  LHC_FASER::signal_handler test_Atlas_signal( "Atlas4jMET0l7TeV",
                                               &test_signal_shortcuts );

  std::cout
  << std::endl
  << *(test_Atlas_signal.getName()) << " constructed without crashing.";
  std::cout
  << std::endl;
  std::cout
  << std::endl
  << *(test_Atlas_signal.getName()) << " has signal value "
  << test_Atlas_signal.getValue();
  std::cout << std::endl;

  SLHA_data.read_file( "SPS2_spectrum.out" );
  std::cout << std::endl;
  std::cout
  << std::endl << "without updating for new point, "
  << *(cross_section_breakdown_test.getName()) << " has signal value "
  << cross_section_breakdown_test.getValue()
  << std::endl
  << *(test_Atlas_signal.getName()) << " has signal value "
  << test_Atlas_signal.getValue();
  std::cout << std::endl;

  test_readier.readyObserversForNewPoint();
  std::cout
  << std::endl << "after updating for new point, "
  << *(cross_section_breakdown_test.getName()) << " has signal value "
  << cross_section_breakdown_test.getValue()
  << std::endl
  << *(test_Atlas_signal.getName()) << " has signal value "
  << test_Atlas_signal.getValue();
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
          = test_input_handler.getExclusiveBrCalculator(
                                               test_input_handler.getGluino(),
                                                   (*BR_iterator)->get_EWino(),
                                                            true,
                                         test_input_handler.getEmptyList() );

          if( NULL != BR_sum_helper )
            {

              gluino_BR_sum += BR_sum_helper->getBr();

            }

          if( !((*BR_iterator)->get_EWino()->counts_as_self_conjugate()) )
            {

              BR_sum_helper
              = test_input_handler.getExclusiveBrCalculator(
                                               test_input_handler.getGluino(),
                                                   (*BR_iterator)->get_EWino(),
                                                                false,
                                         test_input_handler.getEmptyList() );

              if( NULL != BR_sum_helper )
                {

                  gluino_BR_sum += BR_sum_helper->getBr();

                }

            }

        }
      else if( (*BR_iterator)->is_gsx() )
        // otherwise if it's a gluino 2-stage decay to the electroweakino...
        {

          BR_sum_helper
          = test_input_handler.getExclusiveBrCalculator(
                                               test_input_handler.getGluino(),
                                                 (*BR_iterator)->getProduct(),
                                                            true,
                                         test_input_handler.getEmptyList() );
          if( NULL != BR_sum_helper )
            {

              BR_sum_helper_double = BR_sum_helper->getBr();

              BR_sum_helper
              = test_input_handler.getExclusiveBrCalculator(
                                                 (*BR_iterator)->getProduct(),
                                                   (*BR_iterator)->get_EWino(),
                                                                true,
                                         test_input_handler.getEmptyList() );

              if( NULL != BR_sum_helper )
                {

                  gluino_BR_sum
                  += ( 2.0 * BR_sum_helper_double * BR_sum_helper->getBr() );

                }

              if( !((*BR_iterator)->get_EWino()->counts_as_self_conjugate()) )
                {

                  BR_sum_helper
                  = test_input_handler.getExclusiveBrCalculator(
                                                 (*BR_iterator)->getProduct(),
                                                   (*BR_iterator)->get_EWino(),
                                                                    false,
                                         test_input_handler.getEmptyList() );

                  if( NULL != BR_sum_helper )
                    {

                      gluino_BR_sum += ( 2.0 * BR_sum_helper_double
                                             * BR_sum_helper->getBr() );

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
      test_readier.readyObserversForNewPoint();
      std::cout
      << std::endl << "SPS1a: "
      << *(test_Atlas_signal.getName()) << " has signal value "
      << test_Atlas_signal.getValue();
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

      test_readier.readyObserversForNewPoint();
      std::cout
      << std::endl << "SPS2: "
      << *(test_Atlas_signal.getName()) << " has signal value "
      << test_Atlas_signal.getValue();
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


  /* finally, testing the lhcFaser itself: (not working)*//**

  LHC_FASER::lhcFaser SPS1a_UI( "SPS1a_spectrum.out" );
  LHC_FASER::lhcFaser SPS2_UI( "SPS2_spectrum.out" );

  SPS1a_UI.addSignal( "sigmaBreakdownTest" );
  SPS1a_UI.addSignal( "Atlas4jMET0l7TeV" );
  SPS2_UI.addSignal( "badSignalName" );
  SPS2_UI.addSignal( "Atlas4jMET0l7TeV" );

  LHC_FASER::signal_handler* sigma_SPS1a
  = SPS1a_UI.getSignal( "sigmaBreakdownTest" );
  LHC_FASER::signal_handler* Atlas_SPS1a
  = SPS1a_UI.getSignal( "Atlas4jMET0l7TeV" );
  LHC_FASER::signal_handler* Atlas_SPS2
  = SPS2_UI.getSignal( "Atlas4jMET0l7TeV" );
  LHC_FASER::signal_handler* bad_SPS2
  = SPS2_UI.getSignal( "badSignalName" );


  std::cout << std::endl;
  std::cout << std::endl << "all 4 test signals from test UIs have pointers.";
  std::cout << std::endl;

  SPS1a_UI.updateForNewSlha();
  std::cout << std::endl;
  std::cout << std::endl << "SPS1a_UI.updateForNewSlha() successful.";
  std::cout << std::endl;
  SPS2_UI.updateForNewSlha();
  std::cout << std::endl;
  std::cout << std::endl << "SPS2_UI.updateForNewSlha() successful.";
  std::cout << std::endl;
  std::cout
  << std::endl << "after updating for new point, "
  << *(sigma_SPS1a->getName()) << " (sigma_SPS1a) has signal value "
  << sigma_SPS1a->getValue();
  std::cout << std::endl;
  std::cout
  << *(bad_SPS2->getName()) << " (bad_SPS2) has signal value "
  << bad_SPS2->getValue();
  std::cout << std::endl;
  std::cout
  << *(Atlas_SPS1a->getName()) << " (Atlas_SPS1a) has signal value "
  << Atlas_SPS1a->getValue();
  std::cout << std::endl;
  std::cout
  << *(Atlas_SPS2->getName()) << " (Atlas_SPS2) has signal value "
  << Atlas_SPS2->getValue();
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
  LHC_FASER::lhcFaser test_UI( "test_spectrum.out",
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
  test_UI.addSignal( signal_name );
  signals.push_back( test_UI.getSignal( signal_name ) );
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
  test_UI.addSignal( signal_name );
  signals.push_back( test_UI.getSignal( signal_name ) );
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
          test_UI.updateForNewSlha();
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
