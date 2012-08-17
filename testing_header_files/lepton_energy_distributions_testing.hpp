/*
 * lepton_energy_distributions_testing.hpp
 *
 *  Created on: Oct 30, 2011
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

#ifndef LEPTON_ENERGY_DISTRIBUTIONS_TESTING_HPP_
#define LEPTON_ENERGY_DISTRIBUTIONS_TESTING_HPP_

#include "basic_testing_stuff.hpp"
#include "../LHC-FASER_derived_lepton_distributions.hpp"

namespace LHC_FASER
{
  class squarkMassReturner : public effectiveSquarkMassHolder
  {
  public:
    squarkMassReturner(
                  CppSLHA::particle_property_set const* const squarkPointer ) :
        squarkPointer( squarkPointer )
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
      return squarkPointer->get_absolute_mass();
    }

  protected:
    CppSLHA::particle_property_set const* const squarkPointer;
  };


  // this holds the set of distributions & handles adding them.
  class distributionSet
  {
  public:
    distributionSet( readierForNewPoint* const testReadier ) :
        testReadier( testReadier ),
        hardMuon( new hardMuonFromTau() ),
        softMuon( new softMuonFromTau() ),
        hardPion( new hardPionFromTau() ),
        softPion( new softPionFromTau() )
    {
      // just an initialization list.
    }
    ~distributionSet()
    {
      for( std::vector< leptonEnergyDistribution* >::iterator
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
    addDistributions( leptonEnergyDistribution* const inputDistribution,
                      std::string distributionName )
    {
      std::string* testName( new std::string( distributionName ) );
      testName->append( "_ll" );
      distributionsVector.push_back( inputDistribution );
      namesVector.push_back( testName );
      leptonEnergyDistribution*
      tauDistribution( new LHC_FASER::visibleTauDecayProduct( testReadier,
                                                             inputDistribution,
                                                              hardMuon ) );
      testName = new std::string( distributionName );
      testName->append( "_tLm" );
      distributionsVector.push_back( tauDistribution );
      namesVector.push_back( testName );
      tauDistribution = new LHC_FASER::visibleTauDecayProduct( testReadier,
                                                             inputDistribution,
                                                               softMuon );
      testName = new std::string( distributionName );
      testName->append( "_tRm" );
      distributionsVector.push_back( tauDistribution );
      namesVector.push_back( testName );
      tauDistribution = new LHC_FASER::visibleTauDecayProduct( testReadier,
                                                             inputDistribution,
                                                               softPion );
      testName = new std::string( distributionName );
      testName->append( "_tLp" );
      distributionsVector.push_back( tauDistribution );
      namesVector.push_back( testName );
      tauDistribution = new LHC_FASER::visibleTauDecayProduct( testReadier,
                                                             inputDistribution,
                                                               hardPion );
      testName = new std::string( distributionName );
      testName->append( "_tRp" );
      distributionsVector.push_back( tauDistribution );
      namesVector.push_back( testName );
    }

    std::vector< leptonEnergyDistribution* > const*
    getDistributionsVector()
    const
    {
      return &distributionsVector;
    }

    std::vector< std::string* > const*
    getNamesVector()
    const
    {
      return &namesVector;
    }


  protected:
    readierForNewPoint* testReadier;
    tauDecayCoefficient const* const hardMuon;
    tauDecayCoefficient const* const softMuon;
    tauDecayCoefficient const* const hardPion;
    tauDecayCoefficient const* const softPion;
    std::vector< leptonEnergyDistribution* >
    distributionsVector;
    std::vector< std::string* > namesVector;
  };


  class leptonEnergyDistributionsTestingClass
  {
  public:
    leptonEnergyDistributionsTestingClass( basicStuffTesting* basicStuff ) :
      slhaData( basicStuff->getSlha() ),
      testReadier( basicStuff->getReadier() ),
      testDistributionSet( testReadier ),
      testSpectrum( new CppSLHA::EW_scale_spectrum(
                        basicStuff->getSlha()->inspect_particle_spectrum() ) ),
      currentDistribution( NULL )
    {
      // just an initialization list.
    }
    ~leptonEnergyDistributionsTestingClass()
    {
      // does nothing.
    }


    void
    performTest()
    {
      CppSLHA::particle_property_set const* const
      sdownLPointer( testSpectrum->inspect_particle_property_set(
                                                CppSLHA::PDG_code::sdown_L ) );
      squarkMassReturner sdownLMass( sdownLPointer );
      CppSLHA::particle_property_set* const
      scharmLPointer( testSpectrum->get_particle_property_set(
                                               CppSLHA::PDG_code::scharm_L ) );
      CppSLHA::particle_property_set const* const
      stopOnePointer( testSpectrum->inspect_particle_property_set(
                                               CppSLHA::PDG_code::stop_one ) );
      CppSLHA::particle_property_set const* const
      stopTwoPointer( testSpectrum->inspect_particle_property_set(
                                               CppSLHA::PDG_code::stop_two ) );
      squarkMassReturner stopTwoMass( stopTwoPointer );
      CppSLHA::particle_property_set const* const
      neutralinoFourPointer( testSpectrum->inspect_particle_property_set(
                                        CppSLHA::PDG_code::neutralino_four ) );
      CppSLHA::particle_property_set* const
      neutralinoThreePointer( testSpectrum->get_particle_property_set(
                                       CppSLHA::PDG_code::neutralino_three ) );
      CppSLHA::particle_property_set const* const
      neutralinoTwoPointer( testSpectrum->inspect_particle_property_set(
                                         CppSLHA::PDG_code::neutralino_two ) );
      CppSLHA::particle_property_set const* const
      charginoTwoPointer( testSpectrum->inspect_particle_property_set(
                                           CppSLHA::PDG_code::chargino_two ) );
      CppSLHA::particle_property_set const* const
      charginoOnePointer( testSpectrum->inspect_particle_property_set(
                                           CppSLHA::PDG_code::chargino_one ) );
      CppSLHA::particle_property_set const* const
      smuonLPointer( testSpectrum->inspect_particle_property_set(
                                                CppSLHA::PDG_code::smuon_L ) );
      CppSLHA::particle_property_set const* const
      smuonRPointer( testSpectrum->inspect_particle_property_set(
                                                CppSLHA::PDG_code::smuon_R ) );
      CppSLHA::particle_property_set const* const
      muonSneutrinoLPointer( testSpectrum->inspect_particle_property_set(
                                       CppSLHA::PDG_code::muon_sneutrino_L ) );
      CppSLHA::particle_property_set const* const
      muonSneutrinoRPointer( testSpectrum->inspect_particle_property_set(
                                       CppSLHA::PDG_code::muon_sneutrino_R ) );
      CppSLHA::particle_property_set const* const
      lightNeutralEwsbScalarPointer(
                                   testSpectrum->inspect_particle_property_set(
                              CppSLHA::PDG_code::light_neutral_EWSB_scalar ) );
      CppSLHA::particle_property_set* const
      chargedEwsbScalarPointer( testSpectrum->get_particle_property_set(
                                    CppSLHA::PDG_code::charged_EWSB_scalar ) );
      CppSLHA::particle_property_set const* const
      zPointer( testSpectrum->inspect_particle_property_set(
                                                      CppSLHA::PDG_code::Z ) );
      CppSLHA::particle_property_set const* const
      wPointer( testSpectrum->inspect_particle_property_set(
                                                 CppSLHA::PDG_code::W_plus ) );
      CppSLHA::particle_property_set const* const
      neutralinoOnePointer( testSpectrum->inspect_particle_property_set(
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
      currentDistribution
      = new flatNearMuonPlusAntimuon( testReadier,
                                      slhaData,
                                      sdownLPointer,
                                      &sdownLMass,
                                      neutralinoTwoPointer,
                                      smuonRPointer );
      testDistributionSet.addDistributions( currentDistribution,
                                            "flat_near" );
      currentDistribution
      = new LHC_FASER::sameChiralityNearMuon( testReadier,
                                              slhaData,
                                              sdownLPointer,
                                              &sdownLMass,
                                              neutralinoTwoPointer,
                                              smuonRPointer );
      testDistributionSet.addDistributions( currentDistribution,
                                            "same_near" );
      currentDistribution
      = new LHC_FASER::oppositeChiralityNearMuon( testReadier,
                                                  slhaData,
                                                  sdownLPointer,
                                                  &sdownLMass,
                                                  neutralinoTwoPointer,
                                                  smuonRPointer );
      testDistributionSet.addDistributions( currentDistribution,
                                            "opp_near" );
      currentDistribution
      = new LHC_FASER::flatFarMuonPlusAntimuon( testReadier,
                                                slhaData,
                                                sdownLPointer,
                                                &sdownLMass,
                                                neutralinoTwoPointer,
                                                smuonRPointer,
                                                neutralinoOnePointer );
      testDistributionSet.addDistributions( currentDistribution,
                                            "flat_far" );
      currentDistribution
      = new LHC_FASER::sameChiralityFarMuon( testReadier,
                                             slhaData,
                                             sdownLPointer,
                                             &sdownLMass,
                                             neutralinoTwoPointer,
                                             smuonRPointer,
                                             neutralinoOnePointer );
      testDistributionSet.addDistributions( currentDistribution,
                                            "same_far" );
      currentDistribution
      = new LHC_FASER::oppositeChiralityFarMuon( testReadier,
                                                 slhaData,
                                                 sdownLPointer,
                                                 &sdownLMass,
                                                 neutralinoTwoPointer,
                                                 smuonRPointer,
                                                 neutralinoOnePointer );
      testDistributionSet.addDistributions( currentDistribution,
                                            "opp_far" );
      currentDistribution
      = new LHC_FASER::HiggsMuonPlusAntimuon( testReadier,
                                              slhaData,
                                              sdownLPointer,
                                              &sdownLMass,
                                              neutralinoFourPointer,
                                              lightNeutralEwsbScalarPointer,
                                              neutralinoOnePointer );
      testDistributionSet.addDistributions( currentDistribution,
                                            "hn4" );
      currentDistribution
      = new LHC_FASER::HiggsMuonPlusAntimuon( testReadier,
                                              slhaData,
                                              sdownLPointer,
                                              &sdownLMass,
                                              neutralinoThreePointer,
                                              lightNeutralEwsbScalarPointer,
                                              neutralinoOnePointer );
      testDistributionSet.addDistributions( currentDistribution,
                                            "hn3" );
      currentDistribution
      = new LHC_FASER::HiggsMuonPlusAntimuon( testReadier,
                                              slhaData,
                                              sdownLPointer,
                                              &sdownLMass,
                                              charginoTwoPointer,
                                              chargedEwsbScalarPointer,
                                              neutralinoOnePointer );
      testDistributionSet.addDistributions( currentDistribution,
                                            "H+" );
      currentDistribution
      = new LHC_FASER::zHandedMuon( testReadier,
                                    slhaData,
                                    sdownLPointer,
                                    &sdownLMass,
                                    neutralinoThreePointer,
                                    zPointer,
                                    neutralinoOnePointer,
                                    true,
                                    false );
      testDistributionSet.addDistributions( currentDistribution,
                                            "Zn3L" );
      currentDistribution
      = new LHC_FASER::zHandedMuon( testReadier,
                                    slhaData,
                                    sdownLPointer,
                                    &sdownLMass,
                                    neutralinoThreePointer,
                                    zPointer,
                                    neutralinoOnePointer,
                                    false,
                                    false );
      testDistributionSet.addDistributions( currentDistribution,
                                            "Zn3R" );
      currentDistribution
      = new LHC_FASER::zHandedMuon( testReadier,
                                    slhaData,
                                    sdownLPointer,
                                    &sdownLMass,
                                    neutralinoThreePointer,
                                    zPointer,
                                    neutralinoOnePointer,
                                    false,
                                    true );
      testDistributionSet.addDistributions( currentDistribution,
                                            "Zn3summed" );
      currentDistribution
      = new LHC_FASER::zHandedMuon( testReadier,
                                    slhaData,
                                    sdownLPointer,
                                    &sdownLMass,
                                    neutralinoFourPointer,
                                    zPointer,
                                    neutralinoOnePointer,
                                    true,
                                    false );
      testDistributionSet.addDistributions( currentDistribution,
                                            "Zn4L" );
      currentDistribution
      = new LHC_FASER::zHandedMuon( testReadier,
                                    slhaData,
                                    sdownLPointer,
                                    &sdownLMass,
                                    neutralinoFourPointer,
                                    zPointer,
                                    neutralinoOnePointer,
                                    false,
                                    false );
      testDistributionSet.addDistributions( currentDistribution,
                                            "Zn4R" );
      currentDistribution
      = new LHC_FASER::zHandedMuon( testReadier,
                                    slhaData,
                                    sdownLPointer,
                                    &sdownLMass,
                                    neutralinoFourPointer,
                                    zPointer,
                                    neutralinoOnePointer,
                                    false,
                                    true );
      testDistributionSet.addDistributions( currentDistribution,
                                            "Zn4summed" );
      currentDistribution
      = new LHC_FASER::vectorFromSquarkToMuon( testReadier,
                                               slhaData,
                                               sdownLPointer,
                                               &sdownLMass,
                                               scharmLPointer,
                                               wPointer );
      testDistributionSet.addDistributions( currentDistribution,
                                            "sdscW" );
      currentDistribution
      = new LHC_FASER::scalarFromSquarkToMuon( testReadier,
                                               slhaData,
                                               stopTwoPointer,
                                               &stopTwoMass,
                                               stopOnePointer,
                                               lightNeutralEwsbScalarPointer );
      testDistributionSet.addDistributions( currentDistribution,
                                            "st2st1h" );
      currentDistribution
      = new LHC_FASER::wMinusHandedMuon( testReadier,
                                         slhaData,
                                         sdownLPointer,
                                         &sdownLMass,
                                         charginoOnePointer,
                                         wPointer,
                                         neutralinoOnePointer,
                                         true );
      testDistributionSet.addDistributions( currentDistribution,
                                            "Wx1L" );
      currentDistribution
      = new LHC_FASER::wMinusHandedMuon( testReadier,
                                         slhaData,
                                         sdownLPointer,
                                         &sdownLMass,
                                         charginoOnePointer,
                                         wPointer,
                                         neutralinoOnePointer,
                                         false );
      testDistributionSet.addDistributions( currentDistribution,
                                            "Wx1R" );
      currentDistribution
      = new LHC_FASER::wMinusHandedMuon( testReadier,
                                         slhaData,
                                         sdownLPointer,
                                         &sdownLMass,
                                         charginoTwoPointer,
                                         wPointer,
                                         neutralinoOnePointer,
                                         true );
      testDistributionSet.addDistributions( currentDistribution,
                                            "Wx2L" );
      currentDistribution
      = new LHC_FASER::wMinusHandedMuon( testReadier,
                                         slhaData,
                                         sdownLPointer,
                                         &sdownLMass,
                                         charginoTwoPointer,
                                         wPointer,
                                         neutralinoOnePointer,
                                         false );
      testDistributionSet.addDistributions( currentDistribution,
                                            "Wx2R" );
      currentDistribution
      = new LHC_FASER::neutralinoThreeBodyDecay( testReadier,
                                                 slhaData,
                                                 sdownLPointer,
                                                 &sdownLMass,
                                                 neutralinoTwoPointer,
                                                 neutralinoOnePointer,
                                                 smuonLPointer,
                                                 smuonRPointer );
      testDistributionSet.addDistributions( currentDistribution,
                                            "nn3b" );
      currentDistribution
      = new LHC_FASER::charginoThreeBodyDecay( testReadier,
                                               slhaData,
                                               sdownLPointer,
                                               &sdownLMass,
                                               charginoOnePointer,
                                               neutralinoOnePointer,
                                               smuonLPointer,
                                               smuonRPointer,
                                               muonSneutrinoLPointer,
                                               muonSneutrinoRPointer );
      testDistributionSet.addDistributions( currentDistribution,
                                            "xx3b" );

      std::ofstream gnuplotCommandFile;
      gnuplotCommandFile.open( "./lepton_distributions/test_plotter.input" );
      if( gnuplotCommandFile.is_open() )
        // if the file was successfully opened...
      {
        gnuplotCommandFile
        << "set term postscript enhanced color solid" << std::endl;
        std::string gnuplotDataFileName;
        for( unsigned int graphCount( 0 );
             testDistributionSet.getDistributionsVector()->size() > graphCount;
             ++graphCount )
        {
          // open a file for writing out the data points.
          std::ofstream gnuplotDataFile;
          gnuplotDataFileName.assign( "./lepton_distributions/" );
          gnuplotDataFileName.append(
                   *(testDistributionSet.getNamesVector()->at( graphCount )) );
          gnuplotDataFileName.append( ".dat" );
          gnuplotDataFile.open( gnuplotDataFileName.c_str() );
          int const numberOfBins( 100 );
          double binSize;
          double binPosition;
          double plotValue;
          double areaSum;
          if( gnuplotDataFile.is_open() )
            // if the file was successfully opened...
          {
            gnuplotCommandFile
            << "set output \""
            << *(testDistributionSet.getNamesVector()->at( graphCount ))
            << ".eps\"" << std::endl
            << "set style line 1 lt rgb \"red\" lw 3" << std::endl
            << "plot '"
            << *(testDistributionSet.getNamesVector()->at( graphCount ))
            << ".dat' index 0 notitle with lines ls 1" << std::endl;

            std::cout
            << std::endl
            << *(testDistributionSet.getNamesVector()->at( graphCount ))
            << ": "
            << testDistributionSet.getDistributionsVector()->at( graphCount
                                                          )->getMinimumEnergy()
            << " to "
            << testDistributionSet.getDistributionsVector()->at( graphCount
                                                         )->getMaximumEnergy();
            std::cout << std::endl;

            binSize
            = ( ( testDistributionSet.getDistributionsVector()->at( graphCount
                                                          )->getMaximumEnergy()
                 - testDistributionSet.getDistributionsVector()->at( graphCount
                                                        )->getMinimumEnergy() )
                / numberOfBins );
            areaSum = 0.0;
            for( int binCounter( 0 );
                 numberOfBins >= binCounter;
                 ++binCounter )
            {
              binPosition
              = ( (double)binCounter * binSize
                  + testDistributionSet.getDistributionsVector()->at(
                                            graphCount )->getMinimumEnergy() );
              plotValue
              = testDistributionSet.getDistributionsVector()->at( graphCount
                                                     )->valueAt( binPosition );
              areaSum += ( plotValue * binSize );
              gnuplotDataFile << binPosition << " " << plotValue << std::endl;
            }
            gnuplotDataFile.close();
            std::cout << "area = " << areaSum;
            std::cout << std::endl;
          }
          else
          {
            std::cout
            << std::endl << "error!  output file could not be opened!";
            exit( EXIT_FAILURE );
          }
        }
        gnuplotCommandFile.close();
      }
      else
      {
        std::cout << std::endl << "error! gnuplot file could not be opened!";
        exit( EXIT_FAILURE );
      }
    }


  protected:
    CppSLHA::CppSLHA0 const* const slhaData;
    readierForNewPoint* testReadier;
    distributionSet testDistributionSet;
    CppSLHA::EW_scale_spectrum* testSpectrum;
    leptonEnergyDistribution* currentDistribution;
  };

}  // end of LHC_FASER namespace

#endif /* LEPTON_ENERGY_DISTRIBUTIONS_TESTING_HPP_ */
