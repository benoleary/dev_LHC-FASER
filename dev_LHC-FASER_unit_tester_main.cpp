/*
 * dev_LHC-FASER_unit_tester_main.cpp
 *
 *  Created on: 21 Jan 2011
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
 *      LHC-FASER_base_electroweak_cascade_stuff.hpp
 *      LHC-FASER_base_electroweak_cascade_stuff.cpp
 *      LHC-FASER_base_kinematics_stuff.hpp
 *      LHC-FASER_base_kinematics_stuff.cpp
 *      LHC-FASER_base_lepton_distribution_stuff.hpp
 *      LHC-FASER_base_lepton_distribution_stuff.cpp
 *      LHC-FASER_charged_electroweak_cascade_stuff.hpp
 *      LHC-FASER_charged_electroweak_cascade_stuff.cpp
 *      LHC-FASER_cross-section_stuff.hpp
 *      LHC-FASER_cross-section_stuff.cpp
 *      LHC-FASER_derived_lepton_distributions.hpp
 *      LHC-FASER_derived_lepton_distributions.cpp
 *      LHC-FASER_electroweak_cascade_collection_stuff.hpp
 *      LHC-FASER_electroweak_cascade_collection_stuff.cpp
 *      LHC-FASER_full_cascade_stuff.hpp
 *      LHC-FASER_full_cascade_stuff.cpp
 *      LHC-FASER_global_stuff.hpp
 *      LHC-FASER_global_stuff.cpp
 *      LHC-FASER_input_handling_stuff.hpp
 *      LHC-FASER_input_handling_stuff.cpp
 *      LHC-FASER_jet_kinematics_stuff.hpp
 *      LHC-FASER_jet_kinematics_stuff.cpp
 *      LHC-FASER_lepton_kinematics_stuff.hpp
 *      LHC-FASER_lepton_kinematics_stuff.cpp
 *      LHC-FASER_neutral_electroweak_cascade_stuff.hpp
 *      LHC-FASER_neutral_electroweak_cascade_stuff.cpp
 *      LHC-FASER_signal_calculator_stuff.hpp
 *      LHC-FASER_signal_calculator_stuff.cpp
 *      LHC-FASER_signal_data_collection_stuff.hpp
 *      LHC-FASER_signal_data_collection_stuff.cpp
 *      LHC-FASER_sparticle_decay_stuff.hpp
 *      LHC-FASER_sparticle_decay_stuff.cpp
 *      LHC-FASER_template_classes.hpp
 *      and README.LHC-FASER.txt which describes the package.
 *
 *      LHC-FASER also requires CppSLHA. It should be found in a subdirectory
 *      included with this package.
 *
 *      LHC-FASER also requires grids of lookup values. These should also be
 *      found in a subdirectory included with this package.
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
 * done so far:
 *      LHC-FASER_template_classes.hpp
 *      LHC-FASER_global_stuff.hpp
 *      LHC-FASER_global_stuff.cpp
 *      LHC-FASER_base_lepton_distribution_stuff.hpp
 *      LHC-FASER_base_lepton_distribution_stuff.cpp
 *      LHC-FASER_derived_lepton_distributions.hpp
 *      LHC-FASER_derived_lepton_distributions.cpp
 *      LHC-FASER_sparticle_decay_stuff.hpp
 *      LHC-FASER_sparticle_decay_stuff.cpp
 *      LHC-FASER_input_handling_stuff.hpp
 *      LHC-FASER_input_handling_stuff.cpp
 *      LHC-FASER_cross-section_stuff.hpp
 *      LHC-FASER_cross-section_stuff.cpp
 *      LHC-FASER_base_kinematics_stuff.hpp
 *      LHC-FASER_base_kinematics_stuff.cpp
 *
 *
 * still to test:
 *      LHC-FASER_lepton_kinematics_stuff.hpp
 *      LHC-FASER_lepton_kinematics_stuff.cpp
 *      LHC-FASER_base_electroweak_cascade_stuff.hpp
 *      LHC-FASER_base_electroweak_cascade_stuff.cpp
 *      LHC-FASER_neutral_electroweak_cascade_stuff.hpp
 *      LHC-FASER_neutral_electroweak_cascade_stuff.cpp
 *      LHC-FASER_charged_electroweak_cascade_stuff.hpp
 *      LHC-FASER_charged_electroweak_cascade_stuff.cpp
 *      LHC-FASER_electroweak_cascade_collection_stuff.hpp
 *      LHC-FASER_electroweak_cascade_collection_stuff.cpp
 *      LHC-FASER_full_cascade_stuff.hpp
 *      LHC-FASER_full_cascade_stuff.cpp
 *      LHC-FASER_jet_kinematics_stuff.hpp
 *      LHC-FASER_jet_kinematics_stuff.cpp
 *      LHC-FASER_signal_data_collection_stuff.hpp
 *      LHC-FASER_signal_data_collection_stuff.cpp
 *      LHC-FASER_signal_calculator_stuff.hpp
 *      LHC-FASER_signal_calculator_stuff.cpp
 *      LHC-FASER.hpp
 *      LHC-FASER.cpp
 */


// includes:

#include <sys/time.h>
//#include "CppSLHA/CppSLHA.hpp"
//#include "LHC-FASER_template_classes.hpp"
//#include "LHC-FASER_global_stuff.hpp"
//#include "LHC-FASER_base_lepton_distribution_stuff.hpp"
#include "LHC-FASER_derived_lepton_distributions.hpp"
//#include "LHC-FASER_sparticle_decay_stuff.hpp"
//#include "LHC-FASER_input_handling_stuff.hpp"
#include "LHC-FASER_cross-section_stuff.hpp"
//#include "LHC-FASER_base_kinematics_stuff.hpp"
#include "LHC-FASER_lepton_kinematics_stuff.hpp"

// future includes:
//#include "LHC-FASER_base_electroweak_cascade_stuff.hpp"
//#include "LHC-FASER_neutral_electroweak_cascade_stuff.hpp"
//#include "LHC-FASER_charged_electroweak_cascade_stuff.hpp"
//#include "LHC-FASER_electroweak_cascade_collection_stuff.hpp"
//#include "LHC-FASER_full_cascade_stuff.hpp"
//#include "LHC-FASER_jet_kinematics_stuff.hpp"
//#include "LHC-FASER_signal_data_collection_stuff.hpp"
//#include "LHC-FASER_signal_calculator_stuff.hpp"
//#include "LHC-FASER.hpp"


// this holds the set of distributions & handles adding them.
class distributionSet
{
public:
  distributionSet( LHC_FASER::readierForNewPoint* const readier ) :
    readier( readier ),
    hardMuon( new LHC_FASER::hardMuonFromTau() ),
    softMuon( new LHC_FASER::softMuonFromTau() ),
    hardPion( new LHC_FASER::hardPionFromTau() ),
    softPion( new LHC_FASER::softPionFromTau() )
  {
    // just an initialization list.
  }
  ~distributionSet()
  {
    for( std::vector< LHC_FASER::leptonEnergyDistribution* >::iterator
         deletionIterator( distributionsVector.begin() );
         distributionsVector.end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
    for( std::vector< std::string* >::iterator
         deletionIterator( namesVector.begin() );
         namesVector.end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
  }


  void
  addDistributions(
                  LHC_FASER::leptonEnergyDistribution* const inputDistribution,
                    std::string distributionName )
  {
    std::string* testName = new std::string( distributionName );
    testName->append( "_ll" );
    distributionsVector.push_back( inputDistribution );
    namesVector.push_back( testName );
    LHC_FASER::leptonEnergyDistribution*
    tauDistribution = new LHC_FASER::visibleTauDecayProduct( readier,
                                                             inputDistribution,
                                                             hardMuon );
    testName = new std::string( distributionName );
    testName->append( "_tLm" );
    distributionsVector.push_back( tauDistribution );
    namesVector.push_back( testName );
    tauDistribution = new LHC_FASER::visibleTauDecayProduct( readier,
                                                             inputDistribution,
                                                             softMuon );
    testName = new std::string( distributionName );
    testName->append( "_tRm" );
    distributionsVector.push_back( tauDistribution );
    namesVector.push_back( testName );
    tauDistribution = new LHC_FASER::visibleTauDecayProduct( readier,
                                                             inputDistribution,
                                                             softPion );
    testName = new std::string( distributionName );
    testName->append( "_tLp" );
    distributionsVector.push_back( tauDistribution );
    namesVector.push_back( testName );
    tauDistribution = new LHC_FASER::visibleTauDecayProduct( readier,
                                                             inputDistribution,
                                                             hardPion );
    testName = new std::string( distributionName );
    testName->append( "_tRp" );
    distributionsVector.push_back( tauDistribution );
    namesVector.push_back( testName );
  }

  std::vector< LHC_FASER::leptonEnergyDistribution* > const*
  getDistributionsVector()
  const
  {
    return &distributionsVector;
  }

  std::vector< std::string* > const*
  get_namesVector()
  const
  {
    return &namesVector;
  }

protected:
  LHC_FASER::readierForNewPoint* readier;
  LHC_FASER::tauDecayCoefficient const* const hardMuon;
  LHC_FASER::tauDecayCoefficient const* const softMuon;
  LHC_FASER::tauDecayCoefficient const* const hardPion;
  LHC_FASER::tauDecayCoefficient const* const softPion;
  std::vector< LHC_FASER::leptonEnergyDistribution* >
  distributionsVector;
  std::vector< std::string* > namesVector;
};

class squarkMassReturner : public LHC_FASER::effectiveSquarkMassHolder
{
public:
  squarkMassReturner( LHC_FASER::particlePointer const squark ) :
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
  LHC_FASER::particlePointer const squark;
};


// the main program:

int main( int argumentCount,
          char* argumentStrings[] )
{
  int recalculationAmount( 10 );

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
      std::string argumentInterpretter( argumentStrings[ 1 ] );
      std::stringstream argumentStream( argumentInterpretter );
      argumentStream >> recalculationAmount;
    }

  CppSLHA::CppSLHA2 slhaData( "SPS1a_spectrum.out" );
  slhaData.read_file();
  LHC_FASER::readierForNewPoint testReadier;
  LHC_FASER::inputHandler testInputHandler( &slhaData,
                                            slhaData.get_particle_spectrum(),
                                            "./grids",
                                            &testReadier );
  testInputHandler.setVerbosity( true );
  std::string bluh( "bluh" );


  /* testing the lepton energy distributions in the squark rest frame:
   * (works) *//**

  distributionSet distributions( &testReadier );

  CppSLHA::EW_scale_spectrum* const
  spectrum( slhaData.get_particle_spectrum() );
  LHC_FASER::particlePointer const
  sdownLPointer( spectrum->inspect_particle_property_set(
                                                CppSLHA::PDG_code::sdown_L ) );
  squarkMassReturner sdownLMass( sdownLPointer );
  CppSLHA::particle_property_set* const
  scharmLPointer(
          spectrum->get_particle_property_set( CppSLHA::PDG_code::scharm_L ) );
  LHC_FASER::particlePointer const
  neutralinoFourPointer( spectrum->inspect_particle_property_set(
                                        CppSLHA::PDG_code::neutralino_four ) );
  CppSLHA::particle_property_set* const
  neutralinoThreePointer( spectrum->get_particle_property_set(
                                       CppSLHA::PDG_code::neutralino_three ) );
  LHC_FASER::particlePointer const
  neutralinoTwoPointer( spectrum->inspect_particle_property_set(
                                         CppSLHA::PDG_code::neutralino_two ) );
  LHC_FASER::particlePointer const
  charginoTwoPointer( spectrum->inspect_particle_property_set(
                                           CppSLHA::PDG_code::chargino_two ) );
  LHC_FASER::particlePointer const
  charginoOnePointer( spectrum->inspect_particle_property_set(
                                           CppSLHA::PDG_code::chargino_one ) );
  LHC_FASER::particlePointer const
  smuonLPointer( spectrum->inspect_particle_property_set(
                                                CppSLHA::PDG_code::smuon_L ) );
  LHC_FASER::particlePointer const
  smuonRPointer( spectrum->inspect_particle_property_set(
                                                CppSLHA::PDG_code::smuon_R ) );
  LHC_FASER::particlePointer const
  muonSneutrinoLPointer( spectrum->inspect_particle_property_set(
                                       CppSLHA::PDG_code::muon_sneutrino_L ) );
  LHC_FASER::particlePointer const
  muonSneutrinoRPointer( spectrum->inspect_particle_property_set(
                                       CppSLHA::PDG_code::muon_sneutrino_R ) );
  LHC_FASER::particlePointer const
  lightNeutralEwsbScalarPointer( spectrum->inspect_particle_property_set(
                              CppSLHA::PDG_code::light_neutral_EWSB_scalar ) );
  CppSLHA::particle_property_set* const
  chargedEwsbScalarPointer( spectrum->get_particle_property_set(
                                    CppSLHA::PDG_code::charged_EWSB_scalar ) );
  LHC_FASER::particlePointer const
  zPointer( spectrum->inspect_particle_property_set( CppSLHA::PDG_code::Z ) );
  LHC_FASER::particlePointer const
  wPointer( spectrum->inspect_particle_property_set(
                                                 CppSLHA::PDG_code::W_plus ) );
  LHC_FASER::particlePointer const
  neutralinoOnePointer( spectrum->inspect_particle_property_set(
                                         CppSLHA::PDG_code::neutralino_one ) );

  // here we set n3 to have enough mass to decay to a h/Z + n1:
  neutralinoThreePointer->set_mass( ( 20.0
                           + lightNeutralEwsbScalarPointer->get_absolute_mass()
                               + neutralinoOnePointer->get_absolute_mass() ) );

  // here we set H+ to be light enough for x2 to decay into it plus n1:
  chargedEwsbScalarPointer->set_mass( (
                                      charginoTwoPointer->get_absolute_mass()
                                  - neutralinoOnePointer->get_absolute_mass()
                                           - 20.0 ) );
  // here we set scL to be light enough for sdL to decay to it + W-:
  scharmLPointer->set_mass( sdownLPointer->get_absolute_mass()
                            - wPointer->get_absolute_mass()
                            - 20.0 );

  LHC_FASER::leptonEnergyDistribution*
  currentDistribution( new LHC_FASER::flatNearMuonPlusAntimuon( &testReadier,
                                                                &slhaData,
                                                                sdownLPointer,
                                                                &sdownLMass,
                                                          neutralinoTwoPointer,
                                                             smuonRPointer ) );
  distributions.addDistributions( currentDistribution,
                                  "flat_near" );
  currentDistribution = new LHC_FASER::sameChiralityNearMuon( &testReadier,
                                                              &slhaData,
                                                              sdownLPointer,
                                                              &sdownLMass,
                                                          neutralinoTwoPointer,
                                                              smuonRPointer );
  distributions.addDistributions( currentDistribution,
                                  "same_near" );
  currentDistribution = new LHC_FASER::oppositeChiralityNearMuon( &testReadier,
                                                                  &slhaData,
                                                                 sdownLPointer,
                                                                  &sdownLMass,
                                                          neutralinoTwoPointer,
                                                               smuonRPointer );
  distributions.addDistributions( currentDistribution,
                                  "opp_near" );
  currentDistribution = new LHC_FASER::flatFarMuonPlusAntimuon( &testReadier,
                                                                &slhaData,
                                                                sdownLPointer,
                                                                &sdownLMass,
                                                          neutralinoTwoPointer,
                                                                smuonRPointer,
                                                        neutralinoOnePointer );
  distributions.addDistributions( currentDistribution,
                                  "flat_far" );
  currentDistribution = new LHC_FASER::sameChiralityFarMuon( &testReadier,
                                                             &slhaData,
                                                             sdownLPointer,
                                                             &sdownLMass,
                                                          neutralinoTwoPointer,
                                                             smuonRPointer,
                                                        neutralinoOnePointer );
  distributions.addDistributions( currentDistribution,
                                  "same_far" );
  currentDistribution = new LHC_FASER::oppositeChiralityFarMuon( &testReadier,
                                                                 &slhaData,
                                                                 sdownLPointer,
                                                                 &sdownLMass,
                                                          neutralinoTwoPointer,
                                                                 smuonRPointer,
                                                        neutralinoOnePointer );
  distributions.addDistributions( currentDistribution,
                                  "opp_far" );
  currentDistribution = new LHC_FASER::HiggsMuonPlusAntimuon( &testReadier,
                                                              &slhaData,
                                                              sdownLPointer,
                                                              &sdownLMass,
                                                         neutralinoFourPointer,
                                                 lightNeutralEwsbScalarPointer,
                                                        neutralinoOnePointer );
  distributions.addDistributions( currentDistribution,
                                  "hn4" );
  currentDistribution = new LHC_FASER::HiggsMuonPlusAntimuon( &testReadier,
                                                              &slhaData,
                                                              sdownLPointer,
                                                              &sdownLMass,
                                                        neutralinoThreePointer,
                                                 lightNeutralEwsbScalarPointer,
                                                        neutralinoOnePointer );
  distributions.addDistributions( currentDistribution,
                                  "hn3" );
  currentDistribution = new LHC_FASER::HiggsMuonPlusAntimuon( &testReadier,
                                                              &slhaData,
                                                              sdownLPointer,
                                                              &sdownLMass,
                                                            charginoTwoPointer,
                                                      chargedEwsbScalarPointer,
                                                        neutralinoOnePointer );
  distributions.addDistributions( currentDistribution,
                                  "H+" );
  currentDistribution = new LHC_FASER::zHandedMuon( &testReadier,
                                                    &slhaData,
                                                    sdownLPointer,
                                                    &sdownLMass,
                                                    neutralinoThreePointer,
                                                    zPointer,
                                                    neutralinoOnePointer,
                                                    true,
                                                    false );
  distributions.addDistributions( currentDistribution,
                                  "Zn3L" );
  currentDistribution = new LHC_FASER::zHandedMuon( &testReadier,
                                                    &slhaData,
                                                    sdownLPointer,
                                                    &sdownLMass,
                                                    neutralinoThreePointer,
                                                    zPointer,
                                                    neutralinoOnePointer,
                                                    false,
                                                    false );
  distributions.addDistributions( currentDistribution,
                                  "Zn3R" );
  currentDistribution = new LHC_FASER::zHandedMuon( &testReadier,
                                                    &slhaData,
                                                    sdownLPointer,
                                                    &sdownLMass,
                                                    neutralinoThreePointer,
                                                    zPointer,
                                                    neutralinoOnePointer,
                                                    false,
                                                    true );
  distributions.addDistributions( currentDistribution,
                                  "Zn3summed" );
  currentDistribution = new LHC_FASER::zHandedMuon( &testReadier,
                                                    &slhaData,
                                                    sdownLPointer,
                                                    &sdownLMass,
                                                    neutralinoFourPointer,
                                                    zPointer,
                                                    neutralinoOnePointer,
                                                    true,
                                                    false );
  distributions.addDistributions( currentDistribution,
                                  "Zn4L" );
  currentDistribution = new LHC_FASER::zHandedMuon( &testReadier,
                                                    &slhaData,
                                                    sdownLPointer,
                                                    &sdownLMass,
                                                    neutralinoFourPointer,
                                                    zPointer,
                                                    neutralinoOnePointer,
                                                    false,
                                                    false );
  distributions.addDistributions( currentDistribution,
                                  "Zn4R" );
  currentDistribution = new LHC_FASER::zHandedMuon( &testReadier,
                                                    &slhaData,
                                                    sdownLPointer,
                                                    &sdownLMass,
                                                    neutralinoFourPointer,
                                                    zPointer,
                                                    neutralinoOnePointer,
                                                    false,
                                                    true );
  distributions.addDistributions( currentDistribution,
                                  "Zn4summed" );
  currentDistribution = new LHC_FASER::vectorFromSquarkToMuon( &testReadier,
                                                               &slhaData,
                                                               sdownLPointer,
                                                               &sdownLMass,
                                                               scharmLPointer,
                                                               wPointer );
  distributions.addDistributions( currentDistribution,
                                  "sdscW" );
  currentDistribution = new LHC_FASER::wMinusHandedMuon( &testReadier,
                                                         &slhaData,
                                                         sdownLPointer,
                                                         &sdownLMass,
                                                         charginoOnePointer,
                                                         wPointer,
                                                         neutralinoOnePointer,
                                                         true );
  distributions.addDistributions( currentDistribution,
                                  "Wx1L" );
  currentDistribution = new LHC_FASER::wMinusHandedMuon( &testReadier,
                                                         &slhaData,
                                                         sdownLPointer,
                                                         &sdownLMass,
                                                         charginoOnePointer,
                                                         wPointer,
                                                         neutralinoOnePointer,
                                                         false );
  distributions.addDistributions( currentDistribution,
                                  "Wx1R" );
  currentDistribution = new LHC_FASER::wMinusHandedMuon( &testReadier,
                                                         &slhaData,
                                                         sdownLPointer,
                                                         &sdownLMass,
                                                         charginoTwoPointer,
                                                         wPointer,
                                                         neutralinoOnePointer,
                                                         true );
  distributions.addDistributions( currentDistribution,
                                  "Wx2L" );
  currentDistribution = new LHC_FASER::wMinusHandedMuon( &testReadier,
                                                         &slhaData,
                                                         sdownLPointer,
                                                         &sdownLMass,
                                                         charginoTwoPointer,
                                                         wPointer,
                                                         neutralinoOnePointer,
                                                         false );
  distributions.addDistributions( currentDistribution,
                                  "Wx2R" );
  currentDistribution = new LHC_FASER::neutralinoThreeBodyDecay( &testReadier,
                                                                 &slhaData,
                                                                 sdownLPointer,
                                                                 &sdownLMass,
                                                          neutralinoTwoPointer,
                                                          neutralinoOnePointer,
                                                                 smuonLPointer,
                                                               smuonRPointer );
  distributions.addDistributions( currentDistribution,
                                  "nn3b" );
  currentDistribution = new LHC_FASER::charginoThreeBodyDecay( &testReadier,
                                                               &slhaData,
                                                               sdownLPointer,
                                                               &sdownLMass,
                                                            charginoOnePointer,
                                                          neutralinoOnePointer,
                                                               smuonLPointer,
                                                               smuonRPointer,
                                                         muonSneutrinoLPointer,
                                                       muonSneutrinoRPointer );
  distributions.addDistributions( currentDistribution,
                                  "xx3b" );

  std::ofstream* gnuplotCommandFile( new std::ofstream() );
  gnuplotCommandFile->open( "./lepton_distributions/test_plotter.input" );
  if( gnuplotCommandFile->is_open() )
    // if the file was successfully opened...
  {
    *gnuplotCommandFile
    << "set term postscript enhanced color solid" << std::endl;
    std::string gnuplotDataFileName;
    for( int graphCount( 0 );
         distributions.getDistributionsVector()->size()
         > (unsigned int)graphCount;
         ++graphCount )
    {
      // open a file for writing out the data points.
      std::ofstream* gnuplotDataFile( new std::ofstream() );
      gnuplotDataFileName.assign( "./lepton_distributions/" );
      gnuplotDataFileName.append(
                         *(distributions.get_namesVector()->at( graphCount )) );
      gnuplotDataFileName.append( ".dat" );
      gnuplotDataFile->open( gnuplotDataFileName.c_str() );
      int const numberOfBins( 100 );
      double binSize;
      double binPosition;
      double plotValue;
      double areaSum;
      if( gnuplotDataFile->is_open() )
        // if the file was successfully opened...
      {
        *gnuplotCommandFile
        << "set output \""
        << *(distributions.get_namesVector()->at( graphCount ))
        << ".eps\"" << std::endl
        << "set style line 1 lt rgb \"red\" lw 3" << std::endl
        << "plot '"
        << *(distributions.get_namesVector()->at( graphCount ))
        << ".dat' index 0 notitle with lines ls 1" << std::endl;

        std::cout
        << std::endl
        << *(distributions.get_namesVector()->at( graphCount ))
        << ": "
        << distributions.getDistributionsVector()->at( graphCount
                                                          )->getMinimumEnergy()
        << " to "
        << distributions.getDistributionsVector()->at( graphCount
                                                         )->getMaximumEnergy();
        std::cout << std::endl;

        binSize = ( ( distributions.getDistributionsVector()->at( graphCount
                                                          )->getMaximumEnergy()
                      - distributions.getDistributionsVector()->at( graphCount
                                                        )->getMinimumEnergy() )
                    / numberOfBins );
        areaSum = 0.0;
        for( int binCounter( 0 );
             numberOfBins >= binCounter;
             ++binCounter )
        {
          binPosition
          = ( (double)binCounter * binSize
              + distributions.getDistributionsVector()->at(
                                            graphCount )->getMinimumEnergy() );
          plotValue
          = distributions.getDistributionsVector()->at( graphCount
                                                     )->valueAt( binPosition );
          areaSum += ( plotValue * binSize );
          *gnuplotDataFile << binPosition << " " << plotValue << std::endl;
        }
        gnuplotDataFile->close();
        std::cout << "area = " << areaSum;
        std::cout << std::endl;
      }
      else
      {
        std::cout
        << std::endl << "error!  output file could not be opened!";
        exit( EXIT_FAILURE );
      }
      delete gnuplotDataFile;
    }
    gnuplotCommandFile->close();
  }
  else
  {
    std::cout << std::endl << "error! gnuplot file could not be opened!";
    exit( EXIT_FAILURE );
  }
  // clean up:
  delete gnuplotCommandFile;
  **/


  /* testing the means of checking if something has already been prepared for a
   * new point: (works)*//**

  std::cout
  << std::endl
  << "is verbose? " << testInputHandler.isVerbose()
  << std::endl
  << "set to be not verbose.";
  testInputHandler.setVerbosity( false );
  std::cout
  << std::endl
  << "is verbose? " << testInputHandler.isVerbose()
  << std::endl
  << "set to be verbose.";
  testInputHandler.setVerbosity( true );
  std::cout
  << std::endl
  << "is verbose? " << testInputHandler.isVerbose()
  << std::endl
  << "using NLO? " << testInputHandler.usingNlo()
  << std::endl
  << "set to not use NLO.";
  testInputHandler.setUsingNloFlag( false );
  std::cout
  << std::endl
  << "using NLO? " << testInputHandler.usingNlo()
  << std::endl
  << "set to use NLO.";
  testInputHandler.setUsingNloFlag( true );
  std::cout
  << std::endl
  << "using NLO? " << testInputHandler.usingNlo()
  << std::endl
  << std::endl
  << "pathToGrids = "
  << *(testInputHandler.getPathToGrids())
  << std::endl
  << "pathToCrossSectionGrids = "
  << *(testInputHandler.getPathToCrossSectionGrids())
  << std::endl
  << "(inspect: pathToCrossSectionGrids = "
  << *(testInputHandler.inspectPathToCrossSectionGrids()) << ")"
  << std::endl
  << "pathToKinematicsGrids = "
  << *(testInputHandler.getPathToKinematicsGrids())
  << std::endl
  << "(inspect: pathToKinematicsGrids = "
  << *(testInputHandler.inspectPathToKinematicsGrids()) << ")"
  << std::endl
  << std::endl
  << "g_EM = "
  << testInputHandler.getElectromagneticFineStructureConstant()
  << std::endl
  << std::endl
  << "SPS1a MASS:"
  << testInputHandler.getCppSlha()->get_MASS_as_string()
  << std::endl
  << std::endl
  << *(testInputHandler.getSpectrum()->get_PDG_coded_name(
                                                        CppSLHA::PDG_code::Z ))
  << " "
  << testInputHandler.getSpectrum()->get_PDG_coded_mass( CppSLHA::PDG_code::Z )
  << " GeV"
  << std::endl
  << *(testInputHandler.getWPlus()->get_name()) << " "
  << testInputHandler.getWPlus()->get_mass() << " GeV"
  << std::endl
  << *(testInputHandler.getTop()->get_name()) << " "
  << testInputHandler.getTop()->get_mass() << " GeV"
  << std::endl
  << *(testInputHandler.getGluino()->get_name()) << " "
  << testInputHandler.getGluino()->get_mass() << " GeV"
  << std::endl
  << "averageSquarks4Mass = "
  << testInputHandler.getAverageSquarks4Mass()
  << std::endl
  << "scoloreds: " << std::endl;
  for( std::vector< LHC_FASER::particlePointer >::const_iterator
       setIterator( testInputHandler.getColoredSparticles()->begin() );
       testInputHandler.getColoredSparticles()->end() > setIterator;
       ++setIterator )
    {
      std::cout << *((*setIterator)->get_name()) << std::endl;
    }
  std::cout
  << std::endl
  << "charginos: " << std::endl;
  for( std::vector< LHC_FASER::particlePointer >::const_iterator
       setIterator( testInputHandler.getCharginos()->begin() );
       testInputHandler.getCharginos()->end() > setIterator;
       ++setIterator )
    {
      std::cout << *((*setIterator)->get_name()) << std::endl;
    }
  std::cout
  << std::endl
  << "neutralinos: " << std::endl;
  for( std::vector< LHC_FASER::particlePointer >::const_iterator
       setIterator( testInputHandler.getNeutralinos()->begin() );
       testInputHandler.getNeutralinos()->end() > setIterator;
       ++setIterator )
    {
      std::cout << *((*setIterator)->get_name()) << std::endl;
    }
  std::cout
  << std::endl
  << "unstable neutralinos: " << std::endl;
  for( std::vector< LHC_FASER::particlePointer >::const_iterator
       setIterator( testInputHandler.getUnstableNeutralinos()->begin() );
       testInputHandler.getUnstableNeutralinos()->end() > setIterator;
       ++setIterator )
    {
      std::cout << *((*setIterator)->get_name()) << std::endl;
    }
  std::cout
  << std::endl
  << "electroweakinos: " << std::endl;
  for( std::vector< LHC_FASER::particlePointer >::const_iterator
       setIterator( testInputHandler.getElectroweakinos()->begin() );
       testInputHandler.getElectroweakinos()->end() > setIterator;
       ++setIterator )
    {
      std::cout << *((*setIterator)->get_name()) << std::endl;
    }
  std::cout
  << std::endl
  << "unstable electroweakinos: " << std::endl;
  for( std::vector< LHC_FASER::particlePointer >::const_iterator
       setIterator( testInputHandler.getUnstableElectroweakinos()->begin() );
       testInputHandler.getUnstableElectroweakinos()->end() > setIterator;
       ++setIterator )
    {
      std::cout << *((*setIterator)->get_name()) << std::endl;
    }
  std::cout
  << std::endl
  << "charged sleptons: " << std::endl;
  for( std::vector< LHC_FASER::particlePointer >::const_iterator
       setIterator( testInputHandler.getChargedSleptons()->begin() );
       testInputHandler.getChargedSleptons()->end() > setIterator;
       ++setIterator )
    {
      std::cout << *((*setIterator)->get_name()) << std::endl;
    }
  std::cout
  << std::endl
  << "sneutrinos: " << std::endl;
  for( std::vector< LHC_FASER::particlePointer >::const_iterator
       setIterator( testInputHandler.getSneutrinos()->begin() );
       testInputHandler.getSneutrinos()->end() > setIterator;
       ++setIterator )
    {
      std::cout << *((*setIterator)->get_name()) << std::endl;
    }
  std::cout
  << std::endl
  << "sdowns: " << std::endl;
  for( std::vector< LHC_FASER::particlePointer >::const_iterator
       setIterator( testInputHandler.getSdownTypes()->begin() );
       testInputHandler.getSdownTypes()->end() > setIterator;
       ++setIterator )
    {
      std::cout << *((*setIterator)->get_name()) << std::endl;
    }
  std::cout
  << std::endl
  << "sups: " << std::endl;
  for( std::vector< LHC_FASER::particlePointer >::const_iterator
       setIterator( testInputHandler.getSupTypes()->begin() );
       testInputHandler.getSupTypes()->end() > setIterator;
       ++setIterator )
    {
      std::cout << *((*setIterator)->get_name()) << std::endl;
    }
  std::cout
  << std::endl
  << "jets: " << std::endl;
  for( std::list< int >::const_iterator
       setIterator( testInputHandler.getJets()->begin() );
       testInputHandler.getJets()->end() != setIterator;
       ++setIterator )
    {
      std::cout
      << *(testInputHandler.getSpectrum()->get_PDG_coded_name(
                                                              (*setIterator) ))
      << std::endl;
    }
  std::cout
  << std::endl
  << "jets4: " << std::endl;
  for( std::list< int >::const_iterator
       setIterator( testInputHandler.getJets4()->begin() );
       testInputHandler.getJets4()->end() != setIterator;
       ++setIterator )
    {
      std::cout
      << *(testInputHandler.getSpectrum()->get_PDG_coded_name(
                                                              (*setIterator) ))
      << std::endl;
    }
  std::cout
  << std::endl
  << "jets5: " << std::endl;
  for( std::list< int >::const_iterator
       setIterator( testInputHandler.getJets5()->begin() );
       testInputHandler.getJets5()->end() != setIterator;
       ++setIterator )
    {
      std::cout
      << *(testInputHandler.getSpectrum()->get_PDG_coded_name(
                                                              (*setIterator) ))
      << std::endl;
    }
  std::cout
  << std::endl
  << "not in jets4: " << std::endl;
  for( std::list< int >::const_iterator
       setIterator( testInputHandler.getNotInJets4()->begin() );
       testInputHandler.getNotInJets4()->end() != setIterator;
       ++setIterator )
    {
      std::cout
      << *(testInputHandler.getSpectrum()->get_PDG_coded_name(
                                                              (*setIterator) ))
      << std::endl;
    }
  std::cout
  << std::endl
  << "not in jets5: " << std::endl;
  for( std::list< int >::const_iterator
       setIterator( testInputHandler.getNotInJets5()->begin() );
       testInputHandler.getNotInJets5()->end() != setIterator;
       ++setIterator )
    {
      std::cout
      << *(testInputHandler.getSpectrum()->get_PDG_coded_name(
                                                              (*setIterator) ))
      << std::endl;
    }
  std::cout
  << std::endl
  << "squarks: " << std::endl;
  for( std::vector< LHC_FASER::particlePointer >::const_iterator
       setIterator( testInputHandler.getSquarks()->begin() );
       testInputHandler.getSquarks()->end() > setIterator;
       ++setIterator )
    {
      std::cout << *((*setIterator)->get_name()) << std::endl;
    }
  std::cout
  << std::endl
  << "squarks4: " << std::endl;
  for( std::vector< LHC_FASER::particlePointer >::const_iterator
       setIterator( testInputHandler.getSquarks4()->begin() );
       testInputHandler.getSquarks4()->end() > setIterator;
       ++setIterator )
    {
      std::cout << *((*setIterator)->get_name()) << std::endl;
    }
  std::cout
  << std::endl
  << "squarks5: " << std::endl;
  for( std::vector< LHC_FASER::particlePointer >::const_iterator
       setIterator( testInputHandler.getSquarks5()->begin() );
       testInputHandler.getSquarks5()->end() > setIterator;
       ++setIterator )
    {
      std::cout << *((*setIterator)->get_name()) << std::endl;
    }
  std::cout
  << std::endl
  << "sbottom_1 is in squarks? "
  << testInputHandler.isIn( testInputHandler.getSbottomOne()->get_PDG_code(),
                            testInputHandler.getSquarks() )
  << std::endl
  << "sbottom_1 is in squarks4? "
  << testInputHandler.isIn( testInputHandler.getSbottomOne()->get_PDG_code(),
                            testInputHandler.getSquarks4() )
  << std::endl;
  std::cout
  << std::endl
  << std::endl
  << "sparticle production combinations: " << std::endl;
  for( std::vector< LHC_FASER::signedParticleShortcutPair* >::const_iterator
       setIterator( testInputHandler.getScoloredProductionCombinations(
                                                                  )->begin() );
       testInputHandler.getScoloredProductionCombinations()->end()
       > setIterator;
       ++setIterator )
    {
      if( (*setIterator)->firstIsNotAntiparticle() )
        {
          std::cout << *((*setIterator)->getFirstParticle()->get_name());
        }
      else
        {
          std::cout << *((*setIterator)->getFirstParticle()->get_antiname());
        }
      std::cout << " + ";
      if( (*setIterator)->secondIsNotAntiparticle() )
        {
          std::cout << *((*setIterator)->getSecondParticle()->get_name());
        }
      else
        {
          std::cout
          << *((*setIterator)->getSecondParticle()->get_antiname());
        }
      std::cout << std::endl;
    }

  LHC_FASER::effectiveSquarkMassHolder*
  sdownLEffectiveMass( testInputHandler.getOnShellEffectiveSquarkMass(
                                              testInputHandler.getSdownL() ) );
  LHC_FASER::effectiveSquarkMassHolder*
  sbottomOneEffectiveMass( testInputHandler.getOnShellEffectiveSquarkMass(
                                          testInputHandler.getSbottomOne() ) );
  LHC_FASER::effectiveSquarkMassHolder*
  sdownLEffectiveMass2( testInputHandler.getOnShellEffectiveSquarkMass(
                                              testInputHandler.getSdownL() ) );
  std::cout
  << std::endl
  << "sdown_L effective mass holder = "
  << sdownLEffectiveMass
  << ", mass = " << sdownLEffectiveMass->getEffectiveSquarkMass()
  << std::endl
  << "2nd sdown_L effective mass holder = "
  << sdownLEffectiveMass2
  << std::endl
  << "sbottom_1 effective mass holder = "
  << sbottomOneEffectiveMass
  << ", mass = " << sbottomOneEffectiveMass->getEffectiveSquarkMass()
  << std::endl;


  std::cout
  << std::endl
  << "sdown_L to neutralino_1 fermion handedness = "
  << testInputHandler.quarkOrLeptonLeftHandedness( CppSLHA::PDG_code::sdown_L,
                                            CppSLHA::PDG_code::neutralino_one )
  << std::endl
  << "sdown_R to neutralino_1 fermion handedness = "
  << testInputHandler.quarkOrLeptonLeftHandedness( CppSLHA::PDG_code::sdown_R,
                                            CppSLHA::PDG_code::neutralino_one )
  << std::endl
  << "sbottom_1 to neutralino_1 fermion handedness = "
  << testInputHandler.quarkOrLeptonLeftHandedness(
                                                CppSLHA::PDG_code::sbottom_one,
                                            CppSLHA::PDG_code::neutralino_one )
  << std::endl
  << "sbottom_2 to neutralino_1 fermion handedness = "
  << testInputHandler.quarkOrLeptonLeftHandedness(
                                                CppSLHA::PDG_code::sbottom_two,
                                            CppSLHA::PDG_code::neutralino_one )
  << std::endl;


  LHC_FASER::exclusiveBrCalculator*
  gluinoToSdownLBr( testInputHandler.getExclusiveBrCalculator(
                                                  testInputHandler.getGluino(),
                                                  testInputHandler.getSdownL(),
                                                               true,
                                          testInputHandler.getNotInJets5() ) );
  LHC_FASER::exclusiveBrCalculator*
  gluinoToSdownRBr( testInputHandler.getExclusiveBrCalculator(
                                                  testInputHandler.getGluino(),
                                                  testInputHandler.getSdownR(),
                                                               true,
                                          testInputHandler.getNotInJets5() ) );
  LHC_FASER::exclusiveBrCalculator*
  sdownLToGluinoBr( testInputHandler.getExclusiveBrCalculator(
                                                  testInputHandler.getSdownL(),
                                                  testInputHandler.getGluino(),
                                                               true,
                                          testInputHandler.getNotInJets5() ) );
  LHC_FASER::exclusiveBrCalculator*
  sdownRToGluinoBr( testInputHandler.getExclusiveBrCalculator(
                                                  testInputHandler.getSdownR(),
                                                  testInputHandler.getGluino(),
                                                               true,
                                          testInputHandler.getNotInJets5() ) );
  LHC_FASER::exclusiveBrCalculator*
  gluinoToSdownLBr2( testInputHandler.getExclusiveBrCalculator(
                                                  testInputHandler.getGluino(),
                                                  testInputHandler.getSdownL(),
                                                                true,
                                          testInputHandler.getNotInJets5() ) );
  std::cout
  << std::endl
  << "exclusive BR checker: gluino -> sdown_L"
  << std::endl << "pointer = "
  << gluinoToSdownLBr
  << std::endl << *(gluinoToSdownLBr->getDecayer()->get_name())
  << " -> " << *(gluinoToSdownLBr->getProduct()->get_name())
  << " has BR " << gluinoToSdownLBr->getBr()
  << std::endl
  << "2nd gluino -> sdown_L"
  << std::endl << "pointer = "
  << gluinoToSdownLBr2
  << std::endl << *(gluinoToSdownLBr2->getDecayer()->get_name())
  << " -> " << *(gluinoToSdownLBr2->getProduct()->get_name())
  << " has BR " << gluinoToSdownLBr2->getBr()
  << std::endl
  << "exclusive BR checker: sdown_L -> gluino"
  << std::endl << "pointer = "
  << sdownLToGluinoBr
  << std::endl << *(sdownLToGluinoBr->getDecayer()->get_name())
  << " -> " << *(sdownLToGluinoBr->getProduct()->get_name())
  << " has BR " << sdownLToGluinoBr->getBr()
  << std::endl
  << "gluino -> sdown_R"
  << std::endl << "pointer = "
  << gluinoToSdownRBr
  << std::endl << *(gluinoToSdownRBr->getDecayer()->get_name())
  << " -> " << *(gluinoToSdownRBr->getProduct()->get_name())
  << " has BR " << gluinoToSdownRBr->getBr()
  << std::endl
  << "exclusive BR checker: sdown_R -> gluino"
  << std::endl << "pointer = "
  << sdownLToGluinoBr
  << std::endl << *(sdownRToGluinoBr->getDecayer()->get_name())
  << " -> " << *(sdownRToGluinoBr->getProduct()->get_name())
  << " has BR " << sdownRToGluinoBr->getBr()
  << std::endl
  << std::endl;


  slhaData.read_file( "SPS2_spectrum.out" );
  testReadier.readyObserversForNewPoint();
  std::cout
  << std::endl
  << "SPS2 MASS:"
  << testInputHandler.getCppSlha()->get_MASS_as_string()
  << std::endl
  << std::endl
  << *(testInputHandler.getGluino()->get_name()) << " "
  << testInputHandler.getGluino()->get_mass() << " GeV"
  << std::endl
  << std::endl
  << "charged sleptons: " << std::endl;
  for( std::vector< LHC_FASER::particlePointer >::const_iterator
       setIterator( testInputHandler.getChargedSleptons()->begin() );
       testInputHandler.getChargedSleptons()->end() > setIterator;
       ++setIterator )
    {
      std::cout << *((*setIterator)->get_name()) << std::endl;
    }
  std::cout
  << std::endl
  << "exclusive BR checker: gluino -> sdown_L"
  << std::endl << "pointer = "
  << gluinoToSdownLBr
  << std::endl << *(gluinoToSdownLBr->getDecayer()->get_name())
  << " -> " << *(gluinoToSdownLBr->getProduct()->get_name())
  << " has BR " << gluinoToSdownLBr->getBr()
  << std::endl
  << "2nd gluino -> sdown_L"
  << std::endl << "pointer = "
  << gluinoToSdownLBr2
  << std::endl << *(gluinoToSdownLBr2->getDecayer()->get_name())
  << " -> " << *(gluinoToSdownLBr2->getProduct()->get_name())
  << " has BR " << gluinoToSdownLBr2->getBr()
  << std::endl
  << "exclusive BR checker: sdown_L -> gluino"
  << std::endl << "pointer = "
  << sdownLToGluinoBr
  << std::endl << *(sdownLToGluinoBr->getDecayer()->get_name())
  << " -> " << *(sdownLToGluinoBr->getProduct()->get_name())
  << " has BR " << sdownLToGluinoBr->getBr()
  << std::endl
  << "gluino -> sdown_R"
  << std::endl << "pointer = "
  << gluinoToSdownRBr
  << std::endl << *(gluinoToSdownRBr->getDecayer()->get_name())
  << " -> " << *(gluinoToSdownRBr->getProduct()->get_name())
  << " has BR " << gluinoToSdownRBr->getBr()
  << std::endl
  << "exclusive BR checker: sdown_R -> gluino"
  << std::endl << "pointer = "
  << sdownLToGluinoBr
  << std::endl << *(sdownRToGluinoBr->getDecayer()->get_name())
  << " -> " << *(sdownRToGluinoBr->getProduct()->get_name())
  << " has BR " << sdownRToGluinoBr->getBr()
  << std::endl
  << std::endl;


  slhaData.read_file( "SPS1a_spectrum.out" );
  testReadier.readyObserversForNewPoint();
  std::cout
  << std::endl
  << "SPS1a MASS:"
  << testInputHandler.getCppSlha()->get_MASS_as_string()
  << std::endl
  << std::endl
  << *(testInputHandler.getGluino()->get_name()) << " "
  << testInputHandler.getGluino()->get_mass() << " GeV"
  << std::endl
  << std::endl
  << "charged sleptons: " << std::endl;
  for( std::vector< LHC_FASER::particlePointer >::const_iterator
       setIterator( testInputHandler.getChargedSleptons()->begin() );
       testInputHandler.getChargedSleptons()->end() > setIterator;
       ++setIterator )
    {
      std::cout << *((*setIterator)->get_name()) << std::endl;
    }
  std::cout
  << std::endl
  << "exclusive BR checker: gluino -> sdown_L"
  << std::endl << "pointer = "
  << gluinoToSdownLBr
  << std::endl << *(gluinoToSdownLBr->getDecayer()->get_name())
  << " -> " << *(gluinoToSdownLBr->getProduct()->get_name())
  << " has BR " << gluinoToSdownLBr->getBr()
  << std::endl
  << "2nd gluino -> sdown_L"
  << std::endl << "pointer = "
  << gluinoToSdownLBr2
  << std::endl << *(gluinoToSdownLBr2->getDecayer()->get_name())
  << " -> " << *(gluinoToSdownLBr2->getProduct()->get_name())
  << " has BR " << gluinoToSdownLBr2->getBr()
  << std::endl
  << "exclusive BR checker: sdown_L -> gluino"
  << std::endl << "pointer = "
  << sdownLToGluinoBr
  << std::endl << *(sdownLToGluinoBr->getDecayer()->get_name())
  << " -> " << *(sdownLToGluinoBr->getProduct()->get_name())
  << " has BR " << sdownLToGluinoBr->getBr()
  << std::endl
  << "gluino -> sdown_R"
  << std::endl << "pointer = "
  << gluinoToSdownRBr
  << std::endl << *(gluinoToSdownRBr->getDecayer()->get_name())
  << " -> " << *(gluinoToSdownRBr->getProduct()->get_name())
  << " has BR " << gluinoToSdownRBr->getBr()
  << std::endl
  << "exclusive BR checker: sdown_R -> gluino"
  << std::endl << "pointer = "
  << sdownLToGluinoBr
  << std::endl << *(sdownRToGluinoBr->getDecayer()->get_name())
  << " -> " << *(sdownRToGluinoBr->getProduct()->get_name())
  << " has BR " << sdownRToGluinoBr->getBr()
  << std::endl
  << std::endl;

  slhaData.read_file( "SPS2_spectrum.out" );
  testReadier.readyObserversForNewPoint();
  std::cout
  << std::endl
  << "exclusive BR checker: gluino -> sdown_L"
  << std::endl << "pointer = "
  << gluinoToSdownLBr
  << std::endl << *(gluinoToSdownLBr->getDecayer()->get_name())
  << " -> " << *(gluinoToSdownLBr->getProduct()->get_name())
  << " has BR " << gluinoToSdownLBr->getBr()
  << std::endl
  << "exclusive BR checker: sdown_L -> gluino"
  << std::endl << "pointer = "
  << sdownLToGluinoBr
  << std::endl << *(sdownLToGluinoBr->getDecayer()->get_name())
  << " -> " << *(sdownLToGluinoBr->getProduct()->get_name())
  << " has BR " << sdownLToGluinoBr->getBr()
  << std::endl
  << "gluino -> sdown_R"
  << std::endl << "pointer = "
  << gluinoToSdownRBr
  << std::endl << *(gluinoToSdownRBr->getDecayer()->get_name())
  << " -> " << *(gluinoToSdownRBr->getProduct()->get_name())
  << " has BR " << gluinoToSdownRBr->getBr()
  << std::endl
  << "exclusive BR checker: sdown_R -> gluino"
  << std::endl << "pointer = "
  << sdownLToGluinoBr
  << std::endl << *(sdownRToGluinoBr->getDecayer()->get_name())
  << " -> " << *(sdownRToGluinoBr->getProduct()->get_name())
  << " has BR " << sdownRToGluinoBr->getBr()
  << std::endl
  << std::endl;
  slhaData.read_file( "SPS1a_spectrum.out" );
  testReadier.readyObserversForNewPoint();
  std::cout
  << std::endl
  << "exclusive BR checker: gluino -> sdown_L"
  << std::endl << "pointer = "
  << gluinoToSdownLBr
  << std::endl << *(gluinoToSdownLBr->getDecayer()->get_name())
  << " -> " << *(gluinoToSdownLBr->getProduct()->get_name())
  << " has BR " << gluinoToSdownLBr->getBr()
  << std::endl
  << "exclusive BR checker: sdown_L -> gluino"
  << std::endl << "pointer = "
  << sdownLToGluinoBr
  << std::endl << *(sdownLToGluinoBr->getDecayer()->get_name())
  << " -> " << *(sdownLToGluinoBr->getProduct()->get_name())
  << " has BR " << sdownLToGluinoBr->getBr()
  << std::endl
  << "gluino -> sdown_R"
  << std::endl << "pointer = "
  << gluinoToSdownRBr
  << std::endl << *(gluinoToSdownRBr->getDecayer()->get_name())
  << " -> " << *(gluinoToSdownRBr->getProduct()->get_name())
  << " has BR " << gluinoToSdownRBr->getBr()
  << std::endl
  << "exclusive BR checker: sdown_R -> gluino"
  << std::endl << "pointer = "
  << sdownLToGluinoBr
  << std::endl << *(sdownRToGluinoBr->getDecayer()->get_name())
  << " -> " << *(sdownRToGluinoBr->getProduct()->get_name())
  << " has BR " << sdownRToGluinoBr->getBr()
  << std::endl
  << std::endl
  << "gluino pointer = " << testInputHandler.getGluino()
  << std::endl
  << "sdown_L pointer = " << testInputHandler.getSdownL()
  << std::endl
  << "sdown_R pointer = " << testInputHandler.getSdownR()
  << std::endl;

  LHC_FASER::readierForNewPoint*
  testReadier2( new LHC_FASER::readierForNewPoint() );
  LHC_FASER::exclusiveBrCalculator*
  sdownLCharginoOneMinus( new LHC_FASER::exclusiveBrCalculator(
                                                  testInputHandler.getSdownL(),
                                             testInputHandler.getCharginoOne(),
                                                                false,
                                              testInputHandler.getNotInJets4(),
                                                              testReadier2 ) );
  LHC_FASER::exclusiveBrCalculator*
  sdownLCharginoOneMinus2( new LHC_FASER::exclusiveBrCalculator(
                                                  testInputHandler.getSdownL(),
                                             testInputHandler.getCharginoOne(),
                                                                false,
                                              testInputHandler.getNotInJets4(),
                                                              testReadier2 ) );
  LHC_FASER::exclusiveBrCalculator*
  sdownLCharginoOnePlus( new LHC_FASER::exclusiveBrCalculator(
                                                  testInputHandler.getSdownL(),
                                             testInputHandler.getCharginoOne(),
                                                               true,
                                              testInputHandler.getNotInJets4(),
                                                              testReadier2 ) );
  std::cout
  << std::endl << "sdownLCharginoOneMinus = "
  << sdownLCharginoOneMinus << ", BR = "
  << sdownLCharginoOneMinus->getBr()
  << std::endl << "sdownLCharginoOneMinus2 = "
  << sdownLCharginoOneMinus2 << ", BR = "
  << sdownLCharginoOneMinus2->getBr()
  << std::endl << "sdownLCharginoOnePlus BR = "
  << sdownLCharginoOnePlus->getBr();
  testReadier2->readyObserversForNewPoint();
  delete sdownLCharginoOnePlus;
  std::cout
  << std::endl << "sdownLCharginoOneMinus BR = "
  << sdownLCharginoOneMinus->getBr();
  testReadier2->readyObserversForNewPoint();
  delete testReadier2;
  delete sdownLCharginoOneMinus2;
  std::cout
  << std::endl << "sdownLCharginoOneMinus BR = "
  << sdownLCharginoOneMinus->getBr();
  std::cout
  << std::endl;
  delete sdownLCharginoOneMinus;
  **/


  /* testing the cross-section grid stuff:
   * (works)*//**

  std::cout
  << std::endl
  << "pathToGrids = "
  << *(testInputHandler.getPathToGrids())
  << std::endl;

  std::string
  pathToTestGrid( *(testInputHandler.getPathToCrossSectionGrids()) );
  pathToTestGrid.append( "/7TeV/gluino+gluino_LO.dat" );
  LHC_FASER::squareGrid testGrid( &pathToTestGrid,
                                  &bluh,
                                  NULL,
                                  &testInputHandler );
  LHC_FASER::crossSectionTable gluinoGluinoTable( &testGrid,
                 testInputHandler.getScoloredProductionCombinations()->at( 0 ),
                                                  1.0,
                                                  &testInputHandler );
  std::cout
  << std::endl << "x: " << testGrid.getLowestX() << " - "
  << testGrid.getHighestX();
  std::cout
  << std::endl << "y: " << testGrid.getLowestY() << " - "
  << testGrid.getHighestY();
  std::cout
  << std::endl << "value at ( 200.0, 200.0 ) = " << testGrid.valueAt( 200.0,
                                                                      200.0 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at ( 200.0, 205.0 ) = " << testGrid.valueAt( 200.0,
                                                                      205.0 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at ( 205.0, 200.0 ) = " << testGrid.valueAt( 205.0,
                                                                      200.0 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at ( 201.0, 201.0 ) = " << testGrid.valueAt( 201.0,
                                                                      201.0 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at ( 210.0, 210.0 ) = " << testGrid.valueAt( 210.0,
                                                                      210.0 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at ( 199.5, 201.0 ) = " << testGrid.valueAt( 199.5,
                                                                      201.0 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at ( 199.5, 199.5 ) = " << testGrid.valueAt( 199.5,
                                                                      199.5 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at ( 2000.0, 200.5 ) = " << testGrid.valueAt( 2000.0,
                                                                       200.5 );
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout
  << std::endl << "value at ( 2000.5, 200.5 ) = " << testGrid.valueAt( 2000.5,
                                                                       200.5 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at ( 2020.0, 200.5 ) = " << testGrid.valueAt( 2020.0,
                                                                       200.5 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at ( 2222.2, 200.5 ) = " << testGrid.valueAt( 2222.2,
                                                                       200.5 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at ( 200.5, 2000.0 ) = " << testGrid.valueAt( 200.5,
                                                                      2000.0 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at ( 200.5, 2020.0 ) = " << testGrid.valueAt( 200.5,
                                                                      2020.0 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at ( 200.5, 2222.2 ) = " << testGrid.valueAt( 200.5,
                                                                      2222.2 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at ( 1999.5, 1999.5 ) = " << testGrid.valueAt( 1999.5,
                                                                      1999.5 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at ( 2000.0, 2000.0 ) = " << testGrid.valueAt( 2000.0,
                                                                      2000.0 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at ( 2000.5, 2000.5 ) = " << testGrid.valueAt( 2000.5,
                                                                      2000.5 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at ( 2020.0, 2020.0 ) = " << testGrid.valueAt( 2020.0,
                                                                      2020.0 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at ( 2200.0, 2200.0 ) = " << testGrid.valueAt( 2200.0,
                                                                      2200.0 );
  std::cout << std::endl;
  std::cout
  << std::endl << "value at ( 2200.5, 2200.2 ) = " << testGrid.valueAt( 2200.5,
                                                                      2200.2 );
  std::cout << std::endl;
  std::cout
  << std::endl << "gluinoGluinoTable value = " << gluinoGluinoTable.getValue();

  pathToTestGrid.assign( *(testInputHandler.getPathToCrossSectionGrids()) );
  pathToTestGrid.append( "/7TeV/gluino+gluino_K.dat" );
  LHC_FASER::squareGrid testKFactorGrid( &pathToTestGrid,
                                         &bluh,
                                         NULL,
                                         &testInputHandler );
  std::cout
  << std::endl << "K-factor value at ( 200.5, 200.5 ) = "
  << testKFactorGrid.valueAt( 200.5,
                              200.5 );
  std::cout << std::endl;


  LHC_FASER::crossSectionHandler
  crossSectionHandlerTester( &testInputHandler );
  LHC_FASER::crossSectionTableSet*
  testCrossSectionSet( crossSectionHandlerTester.getTableSet( 7 ) );
  std::cout
  << std::endl << "testCrossSectionSet = " << testCrossSectionSet;
  //N.B.: for the next bit to work, crossSectionHandler::getTable() has to be
  //moved back to public.
  LHC_FASER::crossSectionTable*
  sdGoSigma( testCrossSectionSet->getTable(
             testInputHandler.getScoloredProductionCombinations()->at( 1 ) ) );
  LHC_FASER::crossSectionTable*
  suGoSigma( testCrossSectionSet->getTable(
            testInputHandler.getScoloredProductionCombinations()->at( 10 ) ) );
  LHC_FASER::crossSectionTableSet*
  testCrossSectionSet2( crossSectionHandlerTester.getTableSet( 7 ) );
  std::cout
  << std::endl << "testCrossSectionSet2 = " << testCrossSectionSet2;
  LHC_FASER::crossSectionTable*
  sdGoSigma2( testCrossSectionSet2->getTable(
             testInputHandler.getScoloredProductionCombinations()->at( 1 ) ) );
  std::cout
  << std::endl << "no extra squareGrid loading?";

  pathToTestGrid.assign( *(testInputHandler.getPathToCrossSectionGrids()) );
  pathToTestGrid.append( "/7TeV/sdown_L+gluino_LO.dat" );
  LHC_FASER::squareGrid testSdGoLo( &pathToTestGrid,
                                    &bluh,
                                    NULL,
                                    &testInputHandler );
  pathToTestGrid.assign( *(testInputHandler.getPathToCrossSectionGrids()) );
  pathToTestGrid.append( "/7TeV/sdown_L+gluino_K.dat" );
  LHC_FASER::squareGrid testSdGoK( &pathToTestGrid,
                                   &bluh,
                                   NULL,
                                   &testInputHandler );
  double sdMass( testInputHandler.getScoloredProductionCombinations()->at( 1
                                  )->getFirstParticle()->get_absolute_mass() );
  double goMass( testInputHandler.getScoloredProductionCombinations()->at( 1
                                 )->getSecondParticle()->get_absolute_mass() );
  std::cout << std::endl;
  std::cout
  << std::endl << "sd mass = " << sdMass
  << std::endl << "go mass = " << goMass
  << std::endl << "testSdGoLo value = " << testSdGoLo.valueAt( sdMass,
                                                               goMass )
  << std::endl << "testSdGoK value = " << testSdGoK.valueAt( sdMass,
                                                             goMass )
  << std::endl << "hence NLO value should be "
  << ( testSdGoLo.valueAt( sdMass,
                           goMass )
       * testSdGoK.valueAt( sdMass,
                            goMass ) )
  << std::endl << "sdGoSigma = " << sdGoSigma
  << std::endl << "sdGoSigma value = " << sdGoSigma->getValue()
  << std::endl << "sdGoSigma2 = " << sdGoSigma2
  << std::endl << "sdGoSigma2 value = " << sdGoSigma->getValue()
  << std::endl << "suGoSigma = " << suGoSigma
  << std::endl << "suGoSigma value = " << suGoSigma->getValue();
  **/


  /* testing the kinematics grid stuff:
   * (works)*//**

  std::string
  acceptanceGridPathToFile( *(testInputHandler.getPathToKinematicsGrids()) );
  std::cout
  << std::endl << "path to kinematics grids = " << acceptanceGridPathToFile;
  std::cout << std::endl;
  acceptanceGridPathToFile.append(
                        "/7TeV/jets/Atlas4jMET/squark+gluino_acceptance.dat" );
  LHC_FASER::acceptanceGrid
  jetsTestAcceptance( &acceptanceGridPathToFile );
  std::cout
  << std::endl << "jetsTestAcceptance.valueAt( "
  << jetsTestAcceptance.getLowestSquarkMass() << ", "
  << jetsTestAcceptance.getLowestGluinoMass() << ", 0.0, 0.0, 0 ) = "
  << jetsTestAcceptance.valueAt( jetsTestAcceptance.getLowestSquarkMass(),
                                 jetsTestAcceptance.getLowestGluinoMass(),
                                 0.0,
                                 0.0,
                                 0,
                                 false,
                                 false )
  << std::endl;
  std::cout
  << "jetsTestAcceptance.valueAt( 123.4, 567.8, 12.3, 45.6, 4 ) = "
  << jetsTestAcceptance.valueAt( 123.4,
                                 567.8,
                                 12.3,
                                 45.6,
                                 4,
                                 false,
                                 false )
  << std::endl;
  std::cout
  << "jetsTestAcceptance.valueAt( 234.5, 678.9, 12.3, 45.6, 4 ) = "
  << jetsTestAcceptance.valueAt( 234.5,
                                 678.9,
                                 12.3,
                                 45.6,
                                 4,
                                 false,
                                 false )
  << std::endl;
  std::cout
  << "jetsTestAcceptance.valueAt( 876.5, 432.1, 1.0, 1.0, 4 ) = "
  << jetsTestAcceptance.valueAt( 876.5,
                                 432.1,
                                 1.0,
                                 1.0,
                                 4,
                                 false,
                                 false )
  << std::endl;
  std::cout
  << "jetsTestAcceptance.valueAt( 501.0, 501.0, 0.0, 0.0, 4 ) = "
  << jetsTestAcceptance.valueAt( 501.0,
                                 501.0,
                                 0.0,
                                 0.0,
                                 4,
                                 false,
                                 false );
  std::cout << std::endl;
  acceptanceGridPathToFile.assign( *(testInputHandler.getPathToKinematicsGrids()) );
  std::cout
  << std::endl << "path to kinematics grids = " << acceptanceGridPathToFile;
  std::cout << std::endl;
  acceptanceGridPathToFile.append( "/7TeV/leptons/squark+gluino_acceptance.dat" );
  LHC_FASER::acceptanceGrid*
  leptonsTestAcceptance
  = new LHC_FASER::acceptanceGrid( &acceptanceGridPathToFile );
  std::cout
  << std::endl << "leptonsTestAcceptance->valueAt( "
  << leptonsTestAcceptance->getLowestSquarkMass() << ", "
  << leptonsTestAcceptance->getLowestGluinoMass() << ", 0.0, 0.0, 0 ) = "
  << leptonsTestAcceptance->valueAt(
                                  leptonsTestAcceptance->getLowestSquarkMass(),
                                  leptonsTestAcceptance->getLowestGluinoMass(),
                                     0.0,
                                     0.0,
                                     0,
                                     true,
                                     false )
  << std::endl
  << "leptonsTestAcceptance->valueAt( 123.4, 567.8, 12.3, 45.6, 4 ) = "
  << leptonsTestAcceptance->valueAt( 123.4,
                                     567.8,
                                     12.3,
                                     45.6,
                                     4,
                                     false,
                                     false )
  << std::endl
  << "leptonsTestAcceptance->valueAt( 234.5, 678.9, 12.3, 45.6, 4 ) = "
  << leptonsTestAcceptance->valueAt( 234.5,
                                     678.9,
                                     12.3,
                                     45.6,
                                     4,
                                     false,
                                     false )
  << std::endl
  << "leptonsTestAcceptance->valueAt( 876.5, 432.1, 1.0, 1.0, 4 ) = "
  << leptonsTestAcceptance->valueAt( 876.5,
                                     432.1,
                                     1.0,
                                     1.0,
                                     4,
                                     false,
                                     false )
  << std::endl
  << "leptonsTestAcceptance->valueAt( 501.0, 501.0, 0.0, 0.0, 4 ) = "
  << leptonsTestAcceptance->valueAt( 501.0,
                                     501.0,
                                     0.0,
                                     0.0,
                                     4,
                                     false,
                                     false );
  std::cout << std::endl;
  delete leptonsTestAcceptance;
  **/


  /* testing the lepton kinematics grid stuff:
   * (works)*//**/

  std::string
  acceptanceGridPathToLeptons(
                              *(testInputHandler.getPathToKinematicsGrids()) );
  std::cout
  << std::endl << "path to kinematics grids = " << acceptanceGridPathToLeptons;
  std::cout << std::endl;
  acceptanceGridPathToLeptons.append( "/7TeV/leptons" );
  LHC_FASER::leptonAcceptanceGrid
  testLeptonAcceptanceGrid( &acceptanceGridPathToLeptons,
                            &testInputHandler );
  std::cout
  << std::endl
  << "testLeptonAcceptanceGrid.getSquarkValue("
  << " sdown_R, neutralino_2, 4 => 6 GeV ) = "
  << testLeptonAcceptanceGrid.getSquarkValue( testInputHandler.getSdownR(),
                                           testInputHandler.getNeutralinoTwo(),
                                              4 );
  std::cout << std::endl;
  std::cout
  << std::endl
  << "testLeptonAcceptanceGrid.getGluinoValue("
  << " neutralino_2, 4 => 6 GeV ) = "
  << testLeptonAcceptanceGrid.getGluinoValue(
                                           testInputHandler.getNeutralinoTwo(),
                                              4,
                                              false,
                                              true );
  std::cout << std::endl;
  std::cout
  << std::endl
  << "testLeptonAcceptanceGrid.getGluinoValue("
  << " neutralino_2, 14 => 26 GeV ) = "
  << testLeptonAcceptanceGrid.getGluinoValue(
                                           testInputHandler.getNeutralinoTwo(),
                                              14,
                                              false,
                                              true );
  std::cout << std::endl;
  std::cout
  << std::endl
  << "testLeptonAcceptanceGrid.getGluinoValue("
  << " neutralino_2, 84 => 166 GeV ) = "
  << testLeptonAcceptanceGrid.getGluinoValue(
                                           testInputHandler.getNeutralinoTwo(),
                                              84,
                                              false,
                                              true );
  std::cout << std::endl;

  LHC_FASER::leptonAcceptanceTable*
  testSquarkBasedLeptonAcceptanceTable
  = new LHC_FASER::squarkBasedLeptonAcceptanceTable( &testLeptonAcceptanceGrid,
                                                testInputHandler.getSdownR(),
                                                     &testInputHandler );
  LHC_FASER::leptonAcceptanceParameterSet
  testLeptonAcceptanceParameterSet( &testInputHandler,
                                    testSquarkBasedLeptonAcceptanceTable,
                                    testInputHandler.getSdownR(),
                                    testInputHandler.getNeutralinoTwo() );
  std::cout
  << std::endl
  << "testLeptonAcceptanceParameterSet.getEffectiveSquarkMass() = "
  << testLeptonAcceptanceParameterSet.getEffectiveSquarkMass();
  std::cout
  << std::endl
  << "testLeptonAcceptanceParameterSet.acceptanceAt( 6.4, 10.0 ) = "
  << testLeptonAcceptanceParameterSet.acceptanceAt( 6.4,
                                                    10.0 );
  std::cout
  << std::endl
  << "testLeptonAcceptanceParameterSet.acceptanceAt( 7.4, 12.3 ) = "
  << testLeptonAcceptanceParameterSet.acceptanceAt( 7.4,
                                                    12.3 );
  std::cout
  << std::endl
  << "testLeptonAcceptanceParameterSet.acceptanceAt( 1007.4, 12.3 ) = "
  << testLeptonAcceptanceParameterSet.acceptanceAt( 1007.4,
                                                    12.3 );
  std::cout << std::endl;
  delete testSquarkBasedLeptonAcceptanceTable;

  LHC_FASER::leptonAcceptanceHandler
  testLeptonAcceptanceHandler( &testInputHandler,
                               testInputHandler.getPathToKinematicsGrids() );
  LHC_FASER::leptonAcceptanceParameterSet*
  testLeptonAcceptanceParameterSetPointer(
                testLeptonAcceptanceHandler.getLeptonAcceptanceParameterSet( 7,
                                                 testInputHandler.getScharmR(),
                                             testInputHandler.getCharginoOne(),
                                                                           2.0,
                                                                       9.0 ) );
  std::cout
  << std::endl
  << "testLeptonAcceptanceParameterSetPointer uses scharm_R rather than"
  << " sdown_L, and importantly that the grid's lepton cut was 9.0 GeV rather"
  << " than 10.0 GeV.";
  std::cout
  << std::endl
  << "testLeptonAcceptanceParameterSetPointer->getEffectiveSquarkMass() = "
  << testLeptonAcceptanceParameterSetPointer->getEffectiveSquarkMass();
  std::cout
  << std::endl
  << "testLeptonAcceptanceParameterSetPointer->acceptanceAt( 6.4, 10.0 ) = "
  << testLeptonAcceptanceParameterSetPointer->acceptanceAt( 6.4,
                                                            10.0 )
  << std::endl
  << "testLeptonAcceptanceParameterSetPointer->acceptanceAt( 7.4, 12.3 ) = "
  << testLeptonAcceptanceParameterSetPointer->acceptanceAt( 7.4,
                                                            12.3 );
  std::cout << std::endl;
   /**/


  /* testing the jet kinematics grid stuff:
   * (needs work)*//**
  acceptanceGridPath.assign( *(testInputHandler.getPathToKinematicsGrids()) );
  acceptanceGridPath.append( "/7TeV/jets/Atlas4jMET" );
  std::cout
  << std::endl << "path to jets kinematics grids = " << acceptanceGridPath;
  std::cout << std::endl;
  LHC_FASER::jetAcceptanceTable testJetKinematicsTable( &acceptanceGridPath,
                                                        &bluh,
                                                        4,
                                                        &testInputHandler );
  for( std::vector< LHC_FASER::signedParticleShortcutPair* >::const_iterator
       pairIterator(
               testInputHandler.getScoloredProductionCombinations()->begin() );
       testInputHandler.getScoloredProductionCombinations()->end()
       > pairIterator;
       ++pairIterator )
  {
    std::cout
    << std::endl
    << "testJetKinematicsTable.getAcceptance( [ "
    << *((*pairIterator)->getFirstParticle()->get_name_or_antiname(
                                   (*pairIterator)->firstIsNotAntiparticle() ))
    << ", "
    << *((*pairIterator)->getSecondParticle()->get_name_or_antiname(
                                  (*pairIterator)->secondIsNotAntiparticle() ))
    << " ], "
    << std::endl
    << " [ "
    << ?? gotta think of something here. silly to look at all possibilities
    << "] )->get_acceptance() = ";


  }

  acceptanceGridPath.assign(
                        *(testInputHandler.getPathToKinematicsGrids()) );
  acceptanceGridPath.append( "/7TeV/leptons" );
  std::cout
  << std::endl
  << "path to leptons kinematics grids = " << acceptanceGridPath;
  std::cout << std::endl;
  LHC_FASER::lepton_acceptance_table test_lepton_kinematics_table(
                                                         &acceptanceGridPath,
                                                                   7,
                                                         &testInputHandler );
  std::cout
  << std::endl << "test_lepton_kinematics_table.get_energy() = "
  << test_lepton_kinematics_table.get_energy()
  << std::endl << "test_lepton_kinematics_table.get_acceptance( "
  << "[ ";
  if( testInputHandler.getScoloredProductionCombinations()->at( sdL_go
                                                       )->is_first_particle() )
    {
      std::cout
      << *(testInputHandler.getScoloredProductionCombinations()->at(
                                    sdL_go )->getFirstParticle()->get_name());
    }
  else
    {
      std::cout
      << *(testInputHandler.getScoloredProductionCombinations()->at(
                                sdL_go )->getFirstParticle()->get_antiname());
    }
  std::cout << " + ";
  if( testInputHandler.getScoloredProductionCombinations()->at( sdL_go
                                                      )->is_second_particle() )
    {
      std::cout
      << *(testInputHandler.getScoloredProductionCombinations()->at(
                                   sdL_go )->getSecondParticle()->get_name());
    }
  else
    {
      std::cout
      << *(testInputHandler.getScoloredProductionCombinations()->at(
                               sdL_go )->getSecondParticle()->get_antiname());
    }
  std::cout << " ]," << std::endl;
  std::cout
  << "[ "
  << *(gluino_cascade3->get_colored_cascades()->at( 1
                                                  )->getDecayer()->get_name())
  << " -> " << *(gluino_cascade3->get_colored_cascades()->at( 1
                                                  )->getProduct()->get_name())
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
                                                  )->getDecayer()->get_name())
  << " -> " << *(sdown_L_cascade->get_colored_cascades()->at( 1
                                                  )->getProduct()->get_name())
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
      testInputHandler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                        )->acceptance_at( 23.4,
                                                                        20.0 );
  std::cout << std::endl;
  std::cout
  << "...->get_acceptance( -1.0, 10.0 )"
  << std::endl << " = "
  << test_lepton_kinematics_table.get_acceptance(
      testInputHandler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                        )->acceptance_at( -1.0,
                                                                        10.0 );
  std::cout << std::endl;
  std::cout
  << "...->get_acceptance( 1.0, 10.0 )"
  << std::endl << " = "
  << test_lepton_kinematics_table.get_acceptance(
      testInputHandler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                        )->acceptance_at( 1.0,
                                                                        10.0 );
  std::cout << std::endl;
  std::cout
  << "...->get_acceptance( 10.0, 10.0 )"
  << std::endl << " = "
  << test_lepton_kinematics_table.get_acceptance(
      testInputHandler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                        )->acceptance_at( 10.0,
                                                                        10.0 );
  std::cout << std::endl;
  std::cout
  << "...->get_acceptance( 12.3, 10.0 )"
  << std::endl << " = "
  << test_lepton_kinematics_table.get_acceptance(
      testInputHandler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                        )->acceptance_at( 12.3,
                                                                        10.0 );
  std::cout << std::endl;
  std::cout
  << "...->get_acceptance( 20.0, 10.0 )"
  << std::endl << " = "
  << test_lepton_kinematics_table.get_acceptance(
      testInputHandler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                        )->acceptance_at( 20.0,
                                                                        10.0 );
  std::cout << std::endl;
  std::cout
  << "...->get_acceptance( 79.9, 10.0 )"
  << std::endl << " = "
  << test_lepton_kinematics_table.get_acceptance(
      testInputHandler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                        )->acceptance_at( 79.9,
                                                                        10.0 );
  std::cout << std::endl;
  std::cout
  << "...->get_acceptance( 80.0, 10.0 )"
  << std::endl << " = "
  << test_lepton_kinematics_table.get_acceptance(
      testInputHandler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                        )->acceptance_at( 80.0,
                                                                        10.0 );
  std::cout << std::endl;
  std::cout
  << "...->get_acceptance( 100.0, 10.0 )"
  << std::endl << " = "
  << test_lepton_kinematics_table.get_acceptance(
      testInputHandler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                       )->acceptance_at( 100.0,
                                                                        10.0 );
  std::cout << std::endl;
  std::cout
  << std::endl << "test_lepton_kinematics_table.get_acceptance( "
  << "[ ";
  if( testInputHandler.getScoloredProductionCombinations()->at( sdL_go
                                                       )->is_first_particle() )
    {
      std::cout
      << *(testInputHandler.getScoloredProductionCombinations()->at(
                                    sdL_go )->getFirstParticle()->get_name());
    }
  else
    {
      std::cout
      << *(testInputHandler.getScoloredProductionCombinations()->at(
                                sdL_go )->getFirstParticle()->get_antiname());
    }
  std::cout << " + ";
  if( testInputHandler.getScoloredProductionCombinations()->at( sdL_go
                                                      )->is_second_particle() )
    {
      std::cout
      << *(testInputHandler.getScoloredProductionCombinations()->at(
                                   sdL_go )->getSecondParticle()->get_name());
    }
  else
    {
      std::cout
      << *(testInputHandler.getScoloredProductionCombinations()->at(
                               sdL_go )->getSecondParticle()->get_antiname());
    }
  std::cout << " ]," << std::endl;
  std::cout
  << "[ "
  << *(gluino_cascade3->get_colored_cascades()->back(
                                                  )->getDecayer()->get_name())
  << " -> " << *(gluino_cascade3->get_colored_cascades()->back(
                                                  )->getProduct()->get_name())
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
                                                  )->getDecayer()->get_name())
  << " -> " << *(sdown_L_cascade->get_colored_cascades()->at( 1
                                                  )->getProduct()->get_name())
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
      testInputHandler.getScoloredProductionCombinations()->at( sdL_go ),
                               gluino_cascade3->get_colored_cascades()->back(),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                        )->acceptance_at( 23.4,
                                                                        20.0 );

  LHC_FASER::kinematics_table
  kinematics_table_tester( &testJetKinematicsTable,
                           4,
                           &test_lepton_kinematics_table );
  std::cout
  << std::endl << "kinematics_table_tester.get_column() = "
  << kinematics_table_tester.get_column()
  << std::endl << "kinematics_table_tester.get_jet_acceptance( "
  << "[ ";
  if( testInputHandler.getScoloredProductionCombinations()->at( sdL_go
                                                       )->is_first_particle() )
    {
      std::cout
      << *(testInputHandler.getScoloredProductionCombinations()->at(
                                    sdL_go )->getFirstParticle()->get_name());
    }
  else
    {
      std::cout
      << *(testInputHandler.getScoloredProductionCombinations()->at(
                                sdL_go )->getFirstParticle()->get_antiname());
    }
  std::cout << " + ";
  if( testInputHandler.getScoloredProductionCombinations()->at( sdL_go
                                                      )->is_second_particle() )
    {
      std::cout
      << *(testInputHandler.getScoloredProductionCombinations()->at(
                                   sdL_go )->getSecondParticle()->get_name());
    }
  else
    {
      std::cout
      << *(testInputHandler.getScoloredProductionCombinations()->at(
                               sdL_go )->getSecondParticle()->get_antiname());
    }
  std::cout << " ]," << std::endl;
  std::cout
  << "[ "
  << *(gluino_cascade3->get_colored_cascades()->at( 1
                                                  )->getDecayer()->get_name())
  << " -> " << *(gluino_cascade3->get_colored_cascades()->at( 1
                                                  )->getProduct()->get_name())
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
                                                  )->getDecayer()->get_name())
  << " -> " << *(sdown_L_cascade->get_colored_cascades()->at( 1
                                                  )->getProduct()->get_name())
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
      testInputHandler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
          sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance();
  std::cout << std::endl;
  std::cout
  << "kinematics_table_tester.get_lepton_acceptance(...)->get_acceptance("
  << " 12.34, 8.0 )"
  << std::endl << " = "
  << kinematics_table_tester.get_lepton_acceptance(
      testInputHandler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                       )->acceptance_at( 12.34,
                                                                         8.0 );

  std::string
  kinematics_at_energy_name(
                        *(testInputHandler.getPathToKinematicsGrids()) );
  kinematics_at_energy_name.append( "/7TeV" );
  std::string jet_signal_name( "Atlas4jMET" );
  LHC_FASER::kinematics_table_set
  kinematics_table_set_tester( 7,
                               &kinematics_at_energy_name,
                               &jet_signal_name,
                               &test_lepton_kinematics_table,
                               &testInputHandler );
  std::cout << std::endl;
  std::cout
  << "kinematics_table_set_tester.getTable( 4"
  << " )->get_jet_acceptance(...)->get_acceptance() = "
  << kinematics_table_set_tester.getTable( 4 )->get_jet_acceptance(
      testInputHandler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
          sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance();
  std::cout << std::endl;
  std::cout
  << "kinematics_table_set_tester.getTable( 4"
  << " )->get_lepton_acceptance(...)->get_acceptance( 12.34, 8.0 ) = "
  << kinematics_table_set_tester.getTable( 4 )->get_lepton_acceptance(
      testInputHandler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                       )->acceptance_at( 12.34,
                                                                         8.0 );
  std::cout << std::endl;

  LHC_FASER::kinematics_handler kinematics_handler_tester(
                                                         &testInputHandler );
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
      testInputHandler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
           sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance()
  << std::endl
  << "->get_lepton_acceptance(...)->get_acceptance( 12.34, 8.0 ) = "
  << kinematics_handler_table_one->get_lepton_acceptance(
      testInputHandler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                       )->acceptance_at( 12.34,
                                                                         8.0 )
  << std::endl
  << "kinematics_handler_table_two = " << kinematics_handler_table_two
  << std::endl
  << "->get_jet_acceptance(...)->get_acceptance() = "
  << kinematics_handler_table_two->get_jet_acceptance(
      testInputHandler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
           sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance()
  << std::endl
  << "->get_lepton_acceptance(...)->get_acceptance( 12.34, 8.0 ) = "
  << kinematics_handler_table_two->get_lepton_acceptance(
      testInputHandler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                       )->acceptance_at( 12.34,
                                                                         8.0 )
  << std::endl
  << "kinematics_handler_table_three = " << kinematics_handler_table_three
  << std::endl
  << "->get_jet_acceptance(...)->get_acceptance() = "
  << kinematics_handler_table_three->get_jet_acceptance(
      testInputHandler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
           sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance()
  << std::endl
  << "->get_lepton_acceptance(...)->get_acceptance( 12.34, 8.0 ) = "
  << kinematics_handler_table_three->get_lepton_acceptance(
      testInputHandler.getScoloredProductionCombinations()->at( sdL_go ),
                              gluino_cascade3->get_colored_cascades()->at( 1 ),
            sdown_L_cascade->get_colored_cascades()->at( 1 ) )->get_acceptance(
                                                       )->acceptance_at( 12.34,
                                                                         8.0 )
  << std::endl;
  **/


  /* testing the cascade calculation stuff:
   * (needs work)*//**

  LHC_FASER::kinematics_handler test_kinematics( &testInputHandler );
  std::string test_signal_name( "Atlas4jMET" );
  LHC_FASER::lepton_acceptance_value*
  test_lepton_acceptance_value
  = test_kinematics.getTable( 7,
                               &test_signal_name,
                               4 )->get_lepton_acceptance(
           testInputHandler.getScoloredProductionCombinations()->at( 1 ),
                               testInputHandler.get_colored_cascade_set(
            testInputHandler.getGluino() )->get_colored_cascades()->at( 1 ),
                               testInputHandler.get_colored_cascade_set(
         testInputHandler.getSdownL() )->get_colored_cascades()->at( 1 ) );

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
                                              testInputHandler.getSdownL(),
                                                                   true,
                                      testInputHandler.getNeutralinoFour(),
                                                                   true,
                                              testInputHandler.getSmuonR(),
                                                           &testInputHandler,
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
                                              testInputHandler.getSdownL(),
                                                                  true,
                                      testInputHandler.getNeutralinoFour(),
                                                                  true,
                                             testInputHandler.getStauOne(),
                                                           &testInputHandler,
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
                                              testInputHandler.getSdownL(),
                                                            true,
                                      testInputHandler.getNeutralinoFour(),
                                                            true,
                                                    testInputHandler.getZ(),
                                                           &testInputHandler,
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
  = testInputHandler.getExclusiveBrCalculator(
                                      testInputHandler.getNeutralinoFour(),
                                              testInputHandler.getSmuonR(),
                                              true,
                               testInputHandler.getEmptyList() )->getBr();
  double
  smuon_R_to_neutralino_one
  = testInputHandler.getExclusiveBrCalculator(
                                              testInputHandler.getSmuonR(),
                                       testInputHandler.getNeutralinoOne(),
                                                    true,
                               testInputHandler.getEmptyList() )->getBr();
  double
  neutralino_four_to_stau_one
  = testInputHandler.getExclusiveBrCalculator(
                                      testInputHandler.getNeutralinoFour(),
                                             testInputHandler.getStauOne(),
                                                    true,
                               testInputHandler.getEmptyList() )->getBr();
  double
  stau_one_to_neutralino_one
  = testInputHandler.getExclusiveBrCalculator(
                                             testInputHandler.getStauOne(),
                                       testInputHandler.getNeutralinoOne(),
                                                    true,
                               testInputHandler.getEmptyList() )->getBr();
  double
  neutralino_four_to_Z
  = testInputHandler.getNeutralinoFour()->inspect_direct_decay_handler(
                                        )->get_branching_ratio_for_exact_match(
                                    testInputHandler.getZ()->get_PDG_code(),
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
                                              testInputHandler.getSdownL(),
                                                                   true,
                                      testInputHandler.getNeutralinoFour(),
                                                                   true,
                                              testInputHandler.getSmuonR(),
                                                           &testInputHandler,
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
                                              testInputHandler.getSdownL(),
                                                                    true,
                                      testInputHandler.getNeutralinoFour(),
                                                                    true,
                                             testInputHandler.getStauOne(),
                                                           &testInputHandler,
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
                                              testInputHandler.getSdownL(),
                                                            true,
                                      testInputHandler.getNeutralinoFour(),
                                                            true,
                                                    testInputHandler.getZ(),
                                                           &testInputHandler,
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
   = new LHC_FASER::cascade_handler( &testInputHandler );

  LHC_FASER::charge_summed_cascade_calculator*
  charge_summed_test_10_10_40_a
  = cascade_handler_tester->get_charge_summed_calculator( 10.0,
                                                          10.0,
                                                          40.0,
                                                  test_lepton_acceptance_value,
                                              testInputHandler.getSdownL(),
                                                          true,
                                      testInputHandler.getNeutralinoFour(),
                                                          true );
  LHC_FASER::charge_summed_cascade_calculator*
  charge_summed_test_20_10_40
  = cascade_handler_tester->get_charge_summed_calculator( 20.0,
                                                          10.0,
                                                          40.0,
                                                  test_lepton_acceptance_value,
                                              testInputHandler.getSdownL(),
                                                          true,
                                      testInputHandler.getNeutralinoFour(),
                                                          true );
  LHC_FASER::charge_summed_cascade_calculator*
  charge_summed_test_10_10_40_b
  = cascade_handler_tester->get_charge_summed_calculator( 10.0,
                                                          10.0,
                                                          40.0,
                                                  test_lepton_acceptance_value,
                                              testInputHandler.getSdownL(),
                                                          true,
                                      testInputHandler.getNeutralinoFour(),
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

	

  /* testing the signal stuff:
   * (needs work)*//**

  LHC_FASER::colored_cascade_set*
  gluino_BR_sum_tester
  = testInputHandler.get_colored_cascade_set(
                                             testInputHandler.getGluino() );
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
          = testInputHandler.getExclusiveBrCalculator(
                                               testInputHandler.getGluino(),
                                                   (*BR_iterator)->get_EWino(),
                                                            true,
                                         testInputHandler.getEmptyList() );

          if( NULL != BR_sum_helper )
            {

              gluino_BR_sum += BR_sum_helper->getBr();

            }

          if( !((*BR_iterator)->get_EWino()->counts_as_self_conjugate()) )
            {

              BR_sum_helper
              = testInputHandler.getExclusiveBrCalculator(
                                               testInputHandler.getGluino(),
                                                   (*BR_iterator)->get_EWino(),
                                                                false,
                                         testInputHandler.getEmptyList() );

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
          = testInputHandler.getExclusiveBrCalculator(
                                               testInputHandler.getGluino(),
                                                 (*BR_iterator)->getProduct(),
                                                            true,
                                         testInputHandler.getEmptyList() );
          if( NULL != BR_sum_helper )
            {

              BR_sum_helper_double = BR_sum_helper->getBr();

              BR_sum_helper
              = testInputHandler.getExclusiveBrCalculator(
                                                 (*BR_iterator)->getProduct(),
                                                   (*BR_iterator)->get_EWino(),
                                                                true,
                                         testInputHandler.getEmptyList() );

              if( NULL != BR_sum_helper )
                {

                  gluino_BR_sum
                  += ( 2.0 * BR_sum_helper_double * BR_sum_helper->getBr() );

                }

              if( !((*BR_iterator)->get_EWino()->counts_as_self_conjugate()) )
                {

                  BR_sum_helper
                  = testInputHandler.getExclusiveBrCalculator(
                                                 (*BR_iterator)->getProduct(),
                                                   (*BR_iterator)->get_EWino(),
                                                                    false,
                                         testInputHandler.getEmptyList() );

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
                                                         &testInputHandler );

  LHC_FASER::kinematics_handler test_kinematics_handler( &testInputHandler );

  LHC_FASER::cascade_handler test_cascade_handler( &testInputHandler );
  LHC_FASER::signalShortcuts test_signal_shortcuts( &testInputHandler,
                                                   &test_cross_section_handler,
                                                     &test_kinematics_handler,
                                                     &test_cascade_handler );

  LHC_FASER::signal_handler test_bad_signal( "testBadSignal",
                                             &test_signal_shortcuts );

  std::cout
  << std::endl
  << *(test_bad_signal.get_name()) << " has signal value "
  << test_bad_signal.getValue();
  std::cout << std::endl;

  LHC_FASER::signal_handler
  cross_section_breakdown_test( "sigmaBreakdownTest",
                                &test_signal_shortcuts );

  std::cout
  << std::endl
  << *(cross_section_breakdown_test.get_name()) << " has signal value "
  << cross_section_breakdown_test.getValue();
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
  << test_Atlas_signal.getValue();
  std::cout << std::endl;

  slhaData.read_file( "SPS2_spectrum.out" );
  std::cout << std::endl;
  std::cout
  << std::endl << "without updating for new point, "
  << *(cross_section_breakdown_test.get_name()) << " has signal value "
  << cross_section_breakdown_test.getValue()
  << std::endl
  << *(test_Atlas_signal.get_name()) << " has signal value "
  << test_Atlas_signal.getValue();
  std::cout << std::endl;

  testReadier.readyObserversForNewPoint();
  std::cout
  << std::endl << "after updating for new point, "
  << *(cross_section_breakdown_test.get_name()) << " has signal value "
  << cross_section_breakdown_test.getValue()
  << std::endl
  << *(test_Atlas_signal.get_name()) << " has signal value "
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
          = testInputHandler.getExclusiveBrCalculator(
                                               testInputHandler.getGluino(),
                                                   (*BR_iterator)->get_EWino(),
                                                            true,
                                         testInputHandler.getEmptyList() );

          if( NULL != BR_sum_helper )
            {

              gluino_BR_sum += BR_sum_helper->getBr();

            }

          if( !((*BR_iterator)->get_EWino()->counts_as_self_conjugate()) )
            {

              BR_sum_helper
              = testInputHandler.getExclusiveBrCalculator(
                                               testInputHandler.getGluino(),
                                                   (*BR_iterator)->get_EWino(),
                                                                false,
                                         testInputHandler.getEmptyList() );

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
          = testInputHandler.getExclusiveBrCalculator(
                                               testInputHandler.getGluino(),
                                                 (*BR_iterator)->getProduct(),
                                                            true,
                                         testInputHandler.getEmptyList() );
          if( NULL != BR_sum_helper )
            {

              BR_sum_helper_double = BR_sum_helper->getBr();

              BR_sum_helper
              = testInputHandler.getExclusiveBrCalculator(
                                                 (*BR_iterator)->getProduct(),
                                                   (*BR_iterator)->get_EWino(),
                                                                true,
                                         testInputHandler.getEmptyList() );

              if( NULL != BR_sum_helper )
                {

                  gluino_BR_sum
                  += ( 2.0 * BR_sum_helper_double * BR_sum_helper->getBr() );

                }

              if( !((*BR_iterator)->get_EWino()->counts_as_self_conjugate()) )
                {

                  BR_sum_helper
                  = testInputHandler.getExclusiveBrCalculator(
                                                 (*BR_iterator)->getProduct(),
                                                   (*BR_iterator)->get_EWino(),
                                                                    false,
                                         testInputHandler.getEmptyList() );

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

      slhaData.read_file( "SPS1a_spectrum.out" );
      testReadier.readyObserversForNewPoint();
      std::cout
      << std::endl << "SPS1a: "
      << *(test_Atlas_signal.get_name()) << " has signal value "
      << test_Atlas_signal.getValue();
      std::cout << std::endl;

      slhaData.read_file( "SPS2_spectrum.out" );
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

      testReadier.readyObserversForNewPoint();
      std::cout
      << std::endl << "SPS2: "
      << *(test_Atlas_signal.get_name()) << " has signal value "
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
  << *(sigma_SPS1a->get_name()) << " (sigma_SPS1a) has signal value "
  << sigma_SPS1a->getValue();
  std::cout << std::endl;
  std::cout
  << *(bad_SPS2->get_name()) << " (bad_SPS2) has signal value "
  << bad_SPS2->getValue();
  std::cout << std::endl;
  std::cout
  << *(Atlas_SPS1a->get_name()) << " (Atlas_SPS1a) has signal value "
  << Atlas_SPS1a->getValue();
  std::cout << std::endl;
  std::cout
  << *(Atlas_SPS2->get_name()) << " (Atlas_SPS2) has signal value "
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
              << *((*signal_iterator)->get_name()) << ": ";
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
