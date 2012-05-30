/*
 * signal_timing_testing.hpp
 *
 *  Created on: Dec 03, 2011
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

#ifndef SIGNAL_TIMING_TESTING_HPP_
#define SIGNAL_TIMING_TESTING_HPP_

#include "basic_testing_stuff.hpp"
#include "input_testing.hpp"
#include "../LHC-FASER.hpp"

namespace LHC_FASER
{
  typedef std::pair< signalHandler*,
                     std::vector< double >* > signalAndTimes;

  class lhcFaserTimer
  {
  public:
    lhcFaserTimer( std::string const spectrumFileToLoad,
                   lhcFaser* theLhcFaser,
                   std::vector< signalHandler* >* const signalVector,
                   std::string const
                   nameOfFileLoadedByTheLhcFaser = "test_spectrum.out" ) :
      spectrumFileToLoad( spectrumFileToLoad ),
      nameOfFileLoadedByTheLhcFaser( nameOfFileLoadedByTheLhcFaser ),
      theLhcFaser( theLhcFaser ),
      signalsAndTheirTimes(),
      systemCallString( "cp " ),
      startTime(),
      endTime(),
      secondsTaken( 0.0 )
    {
      systemCallString.append( spectrumFileToLoad );
      systemCallString.append( " test_spectrum.out" );

      for( std::vector< signalHandler* >::iterator
           signalIterator( signalVector->begin() );
           signalVector->end() > signalIterator;
           ++signalIterator )
      {
        signalsAndTheirTimes.push_back( new signalAndTimes( *signalIterator,
                                                 new std::vector< double > ) );
      }
    }
    ~lhcFaserTimer()
    {
      for( std::vector< signalAndTimes* >::iterator
           deletionIterator( signalsAndTheirTimes.begin() );
           signalsAndTheirTimes.end() > deletionIterator;
           ++deletionIterator )
      {
        delete (*deletionIterator)->second;
        delete *deletionIterator;
      }
    }

    std::string const*
    getSpectrumFilename()
    const
    {
      return &spectrumFileToLoad;
    }

    void
    loadSpectrumAndCalculate()
    {
      std::cout << std::endl;
      std::cout << std::endl;
      gettimeofday( &startTime,
                    NULL );
      std::cout
      << std::endl << "copying " << spectrumFileToLoad;
      system( systemCallString.c_str() );
      gettimeofday( &endTime,
                    NULL );
      secondsTaken
      = ( (double)( endTime.tv_sec - startTime.tv_sec )
          + 0.000001 * (double)( endTime.tv_usec - startTime.tv_usec ) );
      std::cout
      << " took " << secondsTaken << " seconds.";
      std::cout << std::endl;

      gettimeofday( &startTime,
                    NULL );
      std::cout
      << std::endl << "reading it in & updating";
      theLhcFaser->updateForNewSlha();
      gettimeofday( &endTime,
                    NULL );
      secondsTaken
      = ( (double)( endTime.tv_sec - startTime.tv_sec )
          + 0.000001 * (double)( endTime.tv_usec - startTime.tv_usec ) );
      std::cout
      << " took " << secondsTaken << " seconds.";
      std::cout << std::endl;

      for( std::vector< signalAndTimes* >::iterator
           signalIterator( signalsAndTheirTimes.begin() );
           signalsAndTheirTimes.end() > signalIterator;
           ++signalIterator )
      {
        gettimeofday( &startTime,
                      NULL );
        std::cout
        << std::endl << "calculating the point's value for "
        << *((*signalIterator)->first->getName()) << ": ";
        std::cout << (*signalIterator)->first->getValue();
        gettimeofday( &endTime,
                      NULL );
        secondsTaken
        = ( (double)( endTime.tv_sec - startTime.tv_sec )
            + 0.000001 * (double)( endTime.tv_usec
                - startTime.tv_usec ) );
        std::cout
        << "; took " << secondsTaken << " seconds.";
        std::cout << std::endl;
        (*signalIterator)->second->push_back( secondsTaken );
      }
    }

    void
    summarizeTimes()
    {
      std::cout
      << std::endl
      << spectrumFileToLoad << ":";
      std::cout << std::endl;
      for( std::vector< signalAndTimes* >::iterator
           signalIterator( signalsAndTheirTimes.begin() );
           signalsAndTheirTimes.end() > signalIterator;
           ++signalIterator )
      {
        secondsTaken = 0.0;
        for( std::vector< double >::iterator
             timeIterator( (*signalIterator)->second->begin() );
             (*signalIterator)->second->end() > timeIterator;
             ++timeIterator )
        {
          secondsTaken += *timeIterator;
        }
        secondsTaken *= ( 1.0 / (double)((*signalIterator)->second->size()) );
        std::cout
        << std::endl
        << *((*signalIterator)->first->getName()) << ": " << secondsTaken;
        std::cout << std::endl;
      }
      std::cout << std::endl;
    }


  protected:
    std::string spectrumFileToLoad;
    std::string nameOfFileLoadedByTheLhcFaser;
    lhcFaser* theLhcFaser;
    std::vector< signalAndTimes* > signalsAndTheirTimes;
    std::string systemCallString;
    timeval startTime;
    timeval endTime;
    double secondsTaken;
  };


  class signalTimingTesting
  {
  public:
    signalTimingTesting( basicStuffTesting* basicStuff,
                         inputTesting* inputTester ) :
        basicStuff( basicStuff ),
        inputTester( inputTester ),
        testInputHandler( inputTester->getInput() )
    {
      // just an initialization list.
    }
    ~signalTimingTesting()
    {
      // does nothing.
    }


    void
    performTest( int const recalculationAmount )
    {
      timeval startTime;
      timeval endTime;
      double secondsTaken;

      std::vector< signalHandler* > signalVector;
      signalHandler* signalPointer( NULL );
      std::string signalName;

      //std::string gridPath( "/home/bol/BOL_work/projects/");
      //gridPath.append( "Eclipse_workspace/LHC-FASER_unit_tester/testing/" );
      //gridPath.append( "grids/" );
      std::string gridPath( "./grids/");

      system( "cp LM1_spectrum.out test_spectrum.out" );

      gettimeofday( &startTime,
                    NULL );
      std::cout
      << std::endl << "constructing testLhcFaser";
      lhcFaser testLhcFaser( "test_spectrum.out",
                             gridPath,
                             "fb" );
      testLhcFaser.setVerbosity( true );
      gettimeofday( &endTime,
                    NULL );
      secondsTaken = ( (double)( endTime.tv_sec - startTime.tv_sec )
                       + 0.000001
                         * (double)( endTime.tv_usec - startTime.tv_usec ) );
      std::cout
      << " took " << secondsTaken << " seconds.";
      std::cout << std::endl;
      double initialConstruction( secondsTaken );

      gettimeofday( &startTime,
                    NULL );
      signalName.assign( "Atlas4jMET_0l_7TeV" );
      std::cout
      << std::endl << "adding " << signalName;
      testLhcFaser.addSignal( signalName );
      signalPointer = testLhcFaser.getSignal( signalName );
      if( NULL != signalPointer )
      {
        signalVector.push_back( signalPointer );
      }
      gettimeofday( &endTime,
                    NULL );
      secondsTaken
      = ( (double)( endTime.tv_sec - startTime.tv_sec )
          + 0.000001 * (double)( endTime.tv_usec - startTime.tv_usec ) );
      std::cout
      << " took " << secondsTaken << " seconds.";
      std::cout << std::endl;

      gettimeofday( &startTime,
                    NULL );
      signalName.assign( "Atlas3jMET_1l_7TeV" );
      std::cout
      << std::endl << "adding " << signalName;
      signalPointer = testLhcFaser.addSignal( signalName );
      if( NULL != signalPointer )
      {
        signalVector.push_back( signalPointer );
      }
      gettimeofday( &endTime,
                    NULL );
      secondsTaken
      = ( (double)( endTime.tv_sec - startTime.tv_sec )
          + 0.000001 * (double)( endTime.tv_usec - startTime.tv_usec ) );
      std::cout
      << " took " << secondsTaken << " seconds.";
      std::cout << std::endl;

      gettimeofday( &startTime,
                    NULL );
      signalName.assign( "noJetOrMETCut_sssf_7TeV" );
      std::cout
      << std::endl << "adding " << signalName;
      signalPointer = testLhcFaser.addSignal( signalName );
      if( NULL != signalPointer )
      {
        signalVector.push_back( signalPointer );
      }
      gettimeofday( &endTime,
                    NULL );
      secondsTaken
      = ( (double)( endTime.tv_sec - startTime.tv_sec )
          + 0.000001 * (double)( endTime.tv_usec - startTime.tv_usec ) );
      std::cout
      << " took " << secondsTaken << " seconds.";
      std::cout << std::endl;

      std::vector< lhcFaserTimer* > signalTimers;
      signalTimers.push_back( new lhcFaserTimer( "SPS1a_spectrum.out",
                                                 &testLhcFaser,
                                                 &signalVector ) );
      signalTimers.push_back( new lhcFaserTimer( "SPS1b_spectrum.out",
                                                 &testLhcFaser,
                                                 &signalVector ) );
      signalTimers.push_back( new lhcFaserTimer( "SPS2_spectrum.out",
                                                 &testLhcFaser,
                                                 &signalVector ) );
      signalTimers.push_back( new lhcFaserTimer( "SPS3_spectrum.out",
                                                 &testLhcFaser,
                                                 &signalVector ) );
      /*signalTimers.push_back( new lhcFaserTimer( "SPS4_spectrum.out",
                                                 &testLhcFaser,
                                                 &signalVector ) );*/
      // SPS4 is borked.
      signalTimers.push_back( new lhcFaserTimer( "SPS5_spectrum.out",
                                                 &testLhcFaser,
                                                 &signalVector ) );
      signalTimers.push_back( new lhcFaserTimer( "SPS6_spectrum.out",
                                                 &testLhcFaser,
                                                 &signalVector ) );
      signalTimers.push_back( new lhcFaserTimer( "SPS7_spectrum.out",
                                                 &testLhcFaser,
                                                 &signalVector ) );
      signalTimers.push_back( new lhcFaserTimer( "SPS8_spectrum.out",
                                                 &testLhcFaser,
                                                 &signalVector ) );
      signalTimers.push_back( new lhcFaserTimer( "SPS9_spectrum.out",
                                                 &testLhcFaser,
                                                 &signalVector ) );
      signalTimers.push_back( new lhcFaserTimer( "LM1_spectrum.out",
                                                 &testLhcFaser,
                                                 &signalVector ) );
      signalTimers.push_back( new lhcFaserTimer( "LM8_spectrum.out",
                                                 &testLhcFaser,
                                                 &signalVector ) );
      signalTimers.push_back( new lhcFaserTimer( "SU4_spectrum.out",
                                                 &testLhcFaser,
                                                 &signalVector ) );
      signalTimers.push_back( new lhcFaserTimer( "convoluted_spectrum.out",
                                                 &testLhcFaser,
                                                 &signalVector ) );

      std::cout
      << std::endl << "timing sum of 1st calculations of each signal:";
      gettimeofday( &startTime,
                    NULL );
      std::string systemCallString( "" );
      std::string const* spectrumFile;
      for( std::vector< lhcFaserTimer* >::iterator
           spectrumIterator( signalTimers.begin() );
           signalTimers.end() > spectrumIterator;
           ++spectrumIterator )
      {
        spectrumFile = (*spectrumIterator)->getSpectrumFilename();
        systemCallString.assign( "cp " );
        systemCallString.append( *spectrumFile );
        systemCallString.append( " test_spectrum.out" );
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout
        << std::endl << "copying " << *spectrumFile;
        system( systemCallString.c_str() );
        std::cout << std::endl;
        std::cout
        << std::endl << "reading it in & updating";
        testLhcFaser.updateForNewSlha();
        for( std::vector< signalHandler* >::iterator
             signalIterator( signalVector.begin() );
             signalVector.end() > signalIterator;
             ++signalIterator )
        {
          std::cout
          << std::endl << "calculating "
          << *((*signalIterator)->getName()) << ": ";
          std::cout << (*signalIterator)->getValue();
        }
      }
      gettimeofday( &endTime,
                    NULL );
      secondsTaken
      = ( (double)( endTime.tv_sec - startTime.tv_sec )
          + 0.000001 * (double)( endTime.tv_usec
              - startTime.tv_usec ) );
      std::cout
      << "; it took " << secondsTaken << " seconds.";
      std::cout << std::endl;
      double totalForFirstCalculations( secondsTaken );

      for( int reloadCount( recalculationAmount );
           0 < reloadCount;
           --reloadCount )
      {
        std::cout
        << std::endl << "number of reloads to go: " << reloadCount;
        for( std::vector< lhcFaserTimer* >::iterator
             spectrumIterator( signalTimers.begin() );
             signalTimers.end() > spectrumIterator;
             ++spectrumIterator )
        {
          (*spectrumIterator)->loadSpectrumAndCalculate();
        }
      }
      // now we calculate the averages & display the summary:
      std::cout
      << std::endl
      << "summary:"
      << std::endl
      << "initial construction of the lhcFaser: "
      << initialConstruction << " seconds" << std::endl
      << "total time to run the full set of " << signalTimers.size()
      << " points once: "
      << totalForFirstCalculations << " seconds" << std::endl;
      for( std::vector< lhcFaserTimer* >::iterator
           spectrumIterator( signalTimers.begin() );
           signalTimers.end() > spectrumIterator;
           ++spectrumIterator )
      {
        (*spectrumIterator)->summarizeTimes();
      }
      std::cout << std::endl;
      std::cout
      << std::endl
      << "need to check that these event rates are reasonable."
      << std::endl;
    }


  protected:
    basicStuffTesting* basicStuff;
    inputTesting* inputTester;
    inputHandler const* const testInputHandler;
  };

}  // end of LHC_FASER namespace

#endif /* SIGNAL_TIMING_TESTING_HPP_ */
