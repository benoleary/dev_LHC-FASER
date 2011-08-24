/*
 * LHC-FASER_input_handling_stuff.cpp
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

#include "LHC-FASER_input_handling_stuff.hpp"

namespace LHC_FASER
{
  updateDependentGaugeInputs::updateDependentGaugeInputs(
                                           CppSLHA::SLHA_BLOCK const* smInputs,
                                            particlePointer const wPlusPointer,
                                                particlePointer const zPointer,
                                          readierForNewPoint* const readier ) :
    getsReadiedForNewPoint( readier ),
    smInputs( smInputs ),
    wPlusPointer( wPlusPointer ),
    zPointer( zPointer )
  {
    // just an initialization list.
  }

  updateDependentGaugeInputs::~updateDependentGaugeInputs()
  {
    // does nothing.
  }


  void
  updateDependentGaugeInputs::calculateGaugeCouplings()
  {
    // we work out the gauge couplings from the SMINPUTS BLOCK:
    if( NULL != smInputs )
    {
      electromagneticFineStructureConstant
      = sqrt( ( ( 4.0 * M_PI ) / smInputs->get_entry( 1 ) ) );
      weakCosine = ( wPlusPointer->get_absolute_mass()
                     / zPointer->get_absolute_mass() );
      weakSine = sqrt( ( 1.0 - weakCosine * weakCosine ) );
      hyperchargeGaugeCoupling
      = ( electromagneticFineStructureConstant / weakCosine );
      weakGaugeCoupling = ( electromagneticFineStructureConstant / weakSine );
      colorGaugeCoupling = sqrt( ( 4.0 * M_PI * smInputs->get_entry( 3 ) ) );
    }
    else
    {
      std::cout
      << std::endl
      << "LHC-FASER::warning! input_handler::input_handler() did not"
      <<" provide a valid pointer to SMINPUTS. gauge couplings are set to "
      << CppSLHA::CppSLHA_global::really_wrong_value_string;
      std::cout << std::endl;

      electromagneticFineStructureConstant
      = CppSLHA::CppSLHA_global::really_wrong_value;
      weakSine = CppSLHA::CppSLHA_global::really_wrong_value;
      weakCosine = CppSLHA::CppSLHA_global::really_wrong_value;
      hyperchargeGaugeCoupling = CppSLHA::CppSLHA_global::really_wrong_value;
      weakGaugeCoupling = CppSLHA::CppSLHA_global::really_wrong_value;
      colorGaugeCoupling = CppSLHA::CppSLHA_global::really_wrong_value;
    }
  }



  updateDependentAbsoluteMasses::updateDependentAbsoluteMasses(
                                           particlePointer const gluinoPointer,
                                          readierForNewPoint* const readier ) :
    getsReadiedForNewPoint( readier ),
    gluinoPointer( gluinoPointer )
  {
    // just an initialization list.
  }

  updateDependentAbsoluteMasses::~updateDependentAbsoluteMasses()
  {
    // does nothing.
  }



  updateDependentAverageMass::updateDependentAverageMass(
                                        particleVectorPointer const setPointer,
                                          readierForNewPoint* const readier ) :
    getsReadiedForNewPoint( readier ),
    setPointer( setPointer )
  {
    // just an initialization list.
  }

  updateDependentAverageMass::~updateDependentAverageMass()
  {
    // does nothing.
  }



  inputHandler::inputHandler( CppSLHA::CppSLHA0 const* const cppSlhaPointer,
                       CppSLHA::EW_scale_spectrum const* const spectrumPointer,
                              std::string const pathToGrids,
                              readierForNewPoint* const readier ) :
    readier( readier ),
    verbosity( false ),
    usingNloFlag( true ),
    cppSlhaPointer( cppSlhaPointer ),
    spectrumPointer( spectrumPointer ),
    pathToGrids( pathToGrids ),
    pathToCrossSectionGrids( pathToGrids ),
    pathToKinematicsGrids( pathToGrids ),
    downPointer( spectrumPointer->inspect_particle_property_set(
                                                   CppSLHA::PDG_code::down ) ),
    upPointer( spectrumPointer->inspect_particle_property_set(
                                                     CppSLHA::PDG_code::up ) ),
    strangePointer( spectrumPointer->inspect_particle_property_set(
                                                CppSLHA::PDG_code::strange ) ),
    charmPointer( spectrumPointer->inspect_particle_property_set(
                                                  CppSLHA::PDG_code::charm ) ),
    bottomPointer( spectrumPointer->inspect_particle_property_set(
                                                 CppSLHA::PDG_code::bottom ) ),
    topPointer( spectrumPointer->inspect_particle_property_set(
                                                    CppSLHA::PDG_code::top ) ),
    electronPointer( spectrumPointer->inspect_particle_property_set(
                                               CppSLHA::PDG_code::electron ) ),
    electronNeutrinoPointer( spectrumPointer->inspect_particle_property_set(
                                      CppSLHA::PDG_code::electron_neutrino ) ),
    muonPointer( spectrumPointer->inspect_particle_property_set(
                                                   CppSLHA::PDG_code::muon ) ),
    muonNeutrinoPointer( spectrumPointer->inspect_particle_property_set(
                                          CppSLHA::PDG_code::muon_neutrino ) ),
    tauLeptonPointer( spectrumPointer->inspect_particle_property_set(
                                             CppSLHA::PDG_code::tau_lepton ) ),
    tauNeutrinoPointer( spectrumPointer->inspect_particle_property_set(
                                           CppSLHA::PDG_code::tau_neutrino ) ),
    gluonPointer( spectrumPointer->inspect_particle_property_set(
                                                  CppSLHA::PDG_code::gluon ) ),
    photonPointer( spectrumPointer->inspect_particle_property_set(
                                                 CppSLHA::PDG_code::photon ) ),
    zPointer( spectrumPointer->inspect_particle_property_set(
                                                      CppSLHA::PDG_code::Z ) ),
    wPlusPointer( spectrumPointer->inspect_particle_property_set(
                                                 CppSLHA::PDG_code::W_plus ) ),
    lightNeutralEwsbScalarPointer(
                                spectrumPointer->inspect_particle_property_set(
                              CppSLHA::PDG_code::light_neutral_EWSB_scalar ) ),
    heavyNeutralEwsbScalarPointer(
                                spectrumPointer->inspect_particle_property_set(
                              CppSLHA::PDG_code::heavy_neutral_EWSB_scalar ) ),
    neutralEwsbPseudoscalarPointer(
                                spectrumPointer->inspect_particle_property_set(
                              CppSLHA::PDG_code::neutral_EWSB_pseudoscalar ) ),
    chargedEwsbScalarPointer( spectrumPointer->inspect_particle_property_set(
                                    CppSLHA::PDG_code::charged_EWSB_scalar ) ),
    gluinoPointer( spectrumPointer->inspect_particle_property_set(
                                                 CppSLHA::PDG_code::gluino ) ),
    sdownLPointer( spectrumPointer->inspect_particle_property_set(
                                                CppSLHA::PDG_code::sdown_L ) ),
    supLPointer( spectrumPointer->inspect_particle_property_set(
                                                  CppSLHA::PDG_code::sup_L ) ),
    sstrangeLPointer( spectrumPointer->inspect_particle_property_set(
                                             CppSLHA::PDG_code::sstrange_L ) ),
    scharmLPointer( spectrumPointer->inspect_particle_property_set(
                                               CppSLHA::PDG_code::scharm_L ) ),
    sbottomOnePointer( spectrumPointer->inspect_particle_property_set(
                                            CppSLHA::PDG_code::sbottom_one ) ),
    stopOnePointer( spectrumPointer->inspect_particle_property_set(
                                               CppSLHA::PDG_code::stop_one ) ),
    sdownRPointer( spectrumPointer->inspect_particle_property_set(
                                                CppSLHA::PDG_code::sdown_R ) ),
    supRPointer( spectrumPointer->inspect_particle_property_set(
                                                  CppSLHA::PDG_code::sup_R ) ),
    sstrangeRPointer( spectrumPointer->inspect_particle_property_set(
                                             CppSLHA::PDG_code::sstrange_R ) ),
    scharmRPointer( spectrumPointer->inspect_particle_property_set(
                                               CppSLHA::PDG_code::scharm_R ) ),
    sbottomTwoPointer( spectrumPointer->inspect_particle_property_set(
                                            CppSLHA::PDG_code::sbottom_two ) ),
    stopTwoPointer( spectrumPointer->inspect_particle_property_set(
                                               CppSLHA::PDG_code::stop_two ) ),
    selectronLPointer( spectrumPointer->inspect_particle_property_set(
                                            CppSLHA::PDG_code::selectron_L ) ),
    electronSneutrinoLPointer( spectrumPointer->inspect_particle_property_set(
                                   CppSLHA::PDG_code::electron_sneutrino_L ) ),
    smuonLPointer( spectrumPointer->inspect_particle_property_set(
                                                CppSLHA::PDG_code::smuon_L ) ),
    muonSneutrinoLPointer( spectrumPointer->inspect_particle_property_set(
                                       CppSLHA::PDG_code::muon_sneutrino_L ) ),
    stauOnePointer( spectrumPointer->inspect_particle_property_set(
                                               CppSLHA::PDG_code::stau_one ) ),
    tauSneutrinoLPointer( spectrumPointer->inspect_particle_property_set(
                                        CppSLHA::PDG_code::tau_sneutrino_L ) ),
    selectronRPointer( spectrumPointer->inspect_particle_property_set(
                                            CppSLHA::PDG_code::selectron_R ) ),
    electronSneutrinoRPointer( spectrumPointer->inspect_particle_property_set(
                                   CppSLHA::PDG_code::electron_sneutrino_R ) ),
    smuonRPointer( spectrumPointer->inspect_particle_property_set(
                                                CppSLHA::PDG_code::smuon_R ) ),
    muonSneutrinoRPointer( spectrumPointer->inspect_particle_property_set(
                                       CppSLHA::PDG_code::muon_sneutrino_R ) ),
    stauTwoPointer( spectrumPointer->inspect_particle_property_set(
                                               CppSLHA::PDG_code::stau_two ) ),
    tauSneutrinoRPointer( spectrumPointer->inspect_particle_property_set(
                                        CppSLHA::PDG_code::tau_sneutrino_R ) ),
    neutralinoOnePointer( spectrumPointer->inspect_particle_property_set(
                                         CppSLHA::PDG_code::neutralino_one ) ),
    neutralinoTwoPointer( spectrumPointer->inspect_particle_property_set(
                                         CppSLHA::PDG_code::neutralino_two ) ),
    neutralinoThreePointer( spectrumPointer->inspect_particle_property_set(
                                       CppSLHA::PDG_code::neutralino_three ) ),
    neutralinoFourPointer( spectrumPointer->inspect_particle_property_set(
                                        CppSLHA::PDG_code::neutralino_four ) ),
    charginoOnePointer( spectrumPointer->inspect_particle_property_set(
                                           CppSLHA::PDG_code::chargino_one ) ),
    charginoTwoPointer( spectrumPointer->inspect_particle_property_set(
                                           CppSLHA::PDG_code::chargino_two ) ),
    decayCheckers( new decayCheckerHandler( readier ) ),
    nmixPointer( cppSlhaPointer->inspect_BLOCK( "NMIX" ) ),
    umixPointer( cppSlhaPointer->inspect_BLOCK( "UMIX" ) ),
    vmixPointer( cppSlhaPointer->inspect_BLOCK( "VMIX" ) ),
    sbotmixPointer( cppSlhaPointer->inspect_BLOCK( "SBOTMIX" ) ),
    stopmixPointer( cppSlhaPointer->inspect_BLOCK( "STOPMIX" ) ),
    staumixPointer( cppSlhaPointer->inspect_BLOCK( "STAUMIX" ) ),
    ydPointer( cppSlhaPointer->inspect_BLOCK( "YD" ) ),
    yuPointer( cppSlhaPointer->inspect_BLOCK( "YU" ) ),
    yePointer( cppSlhaPointer->inspect_BLOCK( "YE" ) ),
    hardMuonFromTauFunction(),
    softMuonFromTauFunction(),
    hardPionFromTauFunction(),
    softPionFromTauFunction()
  {
    updateDependentInputValues
    = new updateDependentGaugeInputs( cppSlhaPointer->inspect_BLOCK(
                                                                  "SMINPUTS" ),
                                      wPlusPointer,
                                      zPointer,
                                      readier );

    updateDependentMasses = new updateDependentAbsoluteMasses( gluinoPointer,
                                                               readier );
    updateDependentAverageSquarks4Mass
    = new updateDependentAverageMass( squarks4,
                                      readier );

    // the default cross-sections are those generated with the MSTW2008 PDF
    // set, & the default kinematics grids are those from PYTHIA 8.0:
    pathToCrossSectionGrids.append( "/cross-sections/MSTW2008" );
    pathToKinematicsGrids.append( "/kinematics/PYTHIA8" );


    // find all the shortcut pointers to particles:

    coloredSparticles.push_back( gluinoPointer );
    sparticleProductionCombinations.push_back(
                                 new signedParticleShortcutPair( gluinoPointer,
                                                                          true,
                                                                 gluinoPointer,
                                                                      true ) );
    addSquarkAsColoredSparticleAndInProduction( sdownLPointer );
    squarks4.push_back( sdownLPointer );
    squarks5.push_back( sdownLPointer );
    sdownTypes.push_back( sdownLPointer );
    addSquarkAsColoredSparticleAndInProduction( sdownRPointer );
    squarks4.push_back( sdownRPointer );
    squarks5.push_back( sdownRPointer );
    sdownTypes.push_back( sdownRPointer );
    addSquarkAsColoredSparticleAndInProduction( supLPointer );
    squarks4.push_back( supLPointer );
    squarks5.push_back( supLPointer );
    supTypes.push_back( supLPointer );
    addSquarkAsColoredSparticleAndInProduction( supRPointer );
    squarks4.push_back( supRPointer );
    squarks5.push_back( supRPointer );
    supTypes.push_back( supRPointer );
    addSquarkAsColoredSparticleAndInProduction( sstrangeLPointer );
    squarks4.push_back( sstrangeLPointer );
    squarks5.push_back( sstrangeLPointer );
    sdownTypes.push_back( sstrangeLPointer );
    addSquarkAsColoredSparticleAndInProduction( sstrangeRPointer );
    squarks4.push_back( sstrangeRPointer );
    squarks5.push_back( sstrangeRPointer );
    sdownTypes.push_back( sstrangeRPointer );
    addSquarkAsColoredSparticleAndInProduction( scharmLPointer );
    squarks4.push_back( scharmLPointer );
    squarks5.push_back( scharmLPointer );
    supTypes.push_back( scharmLPointer );
    addSquarkAsColoredSparticleAndInProduction( scharmRPointer );
    squarks4.push_back( scharmRPointer );
    squarks5.push_back( scharmRPointer );
    supTypes.push_back( scharmRPointer );
    addSquarkAsColoredSparticleAndInProduction( sbottomOnePointer );
    squarks5.push_back( sbottomOnePointer );
    sdownTypes.push_back( sbottomOnePointer );
    addSquarkAsColoredSparticleAndInProduction( sbottomTwoPointer );
    squarks5.push_back( sbottomTwoPointer );
    sdownTypes.push_back( sbottomTwoPointer );
    addSquarkAsColoredSparticleAndInProduction( stopOnePointer );
    supTypes.push_back( stopOnePointer );
    addSquarkAsColoredSparticleAndInProduction( stopTwoPointer );
    supTypes.push_back( stopTwoPointer );

    sparticleProductionCombinations.push_back(
                                 new signedParticleShortcutPair( sdownLPointer,
                                                                 true,
                                                                 sdownRPointer,
                                                                 false ) );
    sparticleProductionCombinations.push_back(
                                 new signedParticleShortcutPair( sdownRPointer,
                                                                 true,
                                                                 sdownLPointer,
                                                                 false ) );
    sparticleProductionCombinations.push_back(
                                 new signedParticleShortcutPair( sdownLPointer,
                                                                 true,
                                                                 supRPointer,
                                                                 false ) );
    sparticleProductionCombinations.push_back(
                              new signedParticleShortcutPair( sdownRPointer,
                                                                 true,
                                                                 supLPointer,
                                                                 false ) );
    sparticleProductionCombinations.push_back(
                               new signedParticleShortcutPair( supLPointer,
                                                               true,
                                                               sdownRPointer,
                                                               false ) );
    sparticleProductionCombinations.push_back(
                                   new signedParticleShortcutPair( supRPointer,
                                                                   true,
                                                                 sdownLPointer,
                                                                   false ) );
    sparticleProductionCombinations.push_back(
                                   new signedParticleShortcutPair( supLPointer,
                                                                   true,
                                                                   supRPointer,
                                                                   false ) );
    sparticleProductionCombinations.push_back(
                                   new signedParticleShortcutPair( supRPointer,
                                                                   true,
                                                                   supLPointer,
                                                                   false ) );
    sparticleProductionCombinations.push_back(
                                 new signedParticleShortcutPair( sdownLPointer,
                                                                 true,
                                                                 supLPointer,
                                                                 true ) );
    sparticleProductionCombinations.push_back(
                                 new signedParticleShortcutPair( sdownRPointer,
                                                                 true,
                                                                 supRPointer,
                                                                 true ) );

    exclusiveBrs = new exclusiveBrHandler( &sdownTypes,
                                           &supTypes,
                                           readier );


    // make all the shortcut sparticle lists:

    electroweakinos.push_back( neutralinoOnePointer );
    neutralinos.push_back( neutralinoOnePointer );
    electroweakinos.push_back( neutralinoTwoPointer );
    unstableElectroweakinos.push_back( neutralinoTwoPointer );
    neutralinos.push_back( neutralinoTwoPointer );
    unstableNeutralinos.push_back( neutralinoTwoPointer );
    electroweakinos.push_back( neutralinoThreePointer );
    unstableElectroweakinos.push_back( neutralinoThreePointer );
    neutralinos.push_back( neutralinoThreePointer );
    unstableNeutralinos.push_back( neutralinoThreePointer );
    electroweakinos.push_back( neutralinoFourPointer );
    unstableElectroweakinos.push_back( neutralinoFourPointer );
    neutralinos.push_back( neutralinoFourPointer );
    unstableNeutralinos.push_back( neutralinoFourPointer );
    electroweakinos.push_back( charginoOnePointer );
    unstableElectroweakinos.push_back( charginoOnePointer );
    charginos.push_back( charginoOnePointer );
    electroweakinos.push_back( charginoTwoPointer );
    unstableElectroweakinos.push_back( charginoTwoPointer );
    charginos.push_back( charginoTwoPointer );

    chargedSleptons.push_back( selectronLPointer );
    chargedSleptons.push_back( selectronRPointer );
    chargedSleptons.push_back( smuonLPointer );
    chargedSleptons.push_back( smuonRPointer );
    chargedSleptons.push_back( stauOnePointer );
    chargedSleptons.push_back( stauTwoPointer );

    sneutrinos.push_back( electronSneutrinoLPointer );
    sneutrinos.push_back( electronSneutrinoRPointer );
    sneutrinos.push_back( muonSneutrinoLPointer );
    sneutrinos.push_back( muonSneutrinoRPointer );
    sneutrinos.push_back( tauSneutrinoLPointer );
    sneutrinos.push_back( tauSneutrinoRPointer );

    neutralEwsbScalarsAndPseudoscalars.push_back(
                                               lightNeutralEwsbScalarPointer );
    neutralEwsbScalarsAndPseudoscalars.push_back(
                                               heavyNeutralEwsbScalarPointer );
    neutralEwsbScalarsAndPseudoscalars.push_back(
                                              neutralEwsbPseudoscalarPointer );

    ewsbBosonsAndMassiveVectorBosons.push_back( wPlusPointer );
    ewsbBosonsAndMassiveVectorBosons.push_back( zPointer );
    ewsbBosonsAndMassiveVectorBosons.push_back(
                                               lightNeutralEwsbScalarPointer );
    ewsbBosonsAndMassiveVectorBosons.push_back(
                                               heavyNeutralEwsbScalarPointer );
    ewsbBosonsAndMassiveVectorBosons.push_back(
                                              neutralEwsbPseudoscalarPointer );
    ewsbBosonsAndMassiveVectorBosons.push_back( chargedEwsbScalarPointer );

    // make the jet particle code lists:
    for( std::vector< CppSLHA::particle_property_set* >::const_iterator
         particleIterator
         = spectrumPointer->get_particle_property_sets()->begin();
         spectrumPointer->get_particle_property_sets()->end()
         > particleIterator;
         ++particleIterator )
      // go through each particle in the spectrum...
    {
      if( (*particleIterator)->counts_as_jet() )
        // if we find a jet...
      {
        jets.push_back( (*particleIterator)->get_PDG_code() );
        // note the jet particle.
        if( !((*particleIterator)->counts_as_self_conjugate()) )
        {
          jets.push_back( -((*particleIterator)->get_PDG_code()) );
        }
      }
    }

    notInJets5.push_back( CppSLHA::PDG_code::top );
    notInJets5.push_back( -(CppSLHA::PDG_code::top) );
    notInJets4.assign( notInJets5.begin(),
                       notInJets5.end() );
    notInJets4.push_back( CppSLHA::PDG_code::bottom );
    notInJets4.push_back( -(CppSLHA::PDG_code::bottom) );

    bool notInNotInJets;
    for( std::list< int >::iterator jetIterator = jets.begin();
         jets.end() != jetIterator;
         ++jetIterator )
    {
      // assume that this jet is not in notInJets5:
      notInNotInJets = true;
      for( std::list< int >::const_iterator
           notInJetsIterator = notInJets5.begin();
           ( ( notInJets5.end() != notInJetsIterator )
             &&
             notInNotInJets );
           ++notInJetsIterator )
      {
        if( *notInJetsIterator == *jetIterator )
        {
          notInNotInJets = false;
        }
      }
      if( notInNotInJets )
      {
        jets5.push_back( *jetIterator );
      }
      // assume that this jet is not in notInJets4:
      notInNotInJets = true;
      for( std::list< int >::const_iterator
           notInJetsIterator = notInJets4.begin();
           ( ( notInJets4.end() != notInJetsIterator )
             &&
             notInNotInJets );
           ++notInJetsIterator )
      {
        if( *notInJetsIterator == *jetIterator )
        {
          notInNotInJets = false;
        }
      }
      if( notInNotInJets )
      {
        jets4.push_back( *jetIterator );
      }
    }
  }

  inputHandler::~inputHandler()
  {
    for( std::vector< signedParticleShortcutPair* >::iterator
         deletionIterator = sparticleProductionCombinations.begin();
         sparticleProductionCombinations.end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
    for( std::vector< onshellSquarkMassHolder* >::iterator
         deletionIterator = onshellSquarkEffectiveMasses.begin();
         onshellSquarkEffectiveMasses.end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }

    delete decayCheckers;
    delete exclusiveBrs;
    delete updateDependentInputValues;
    delete updateDependentMasses;
    delete updateDependentAverageSquarks4Mass;
  }


  double
  inputHandler::quarkOrLeptonLeftHandedness( int const sfermionCode,
                                              int const ewinoCode )
  const
  /* this returns the square of the coupling associated with the left-handed
   * projection operator of the Feynman rule for the vertex of the given
   * sfermion with the given electroweakino & the appropriate quark, divided
   * by the sum of the squares of the couplings for the left- & right-handed
   * projection operators, based on Rosiek's Feynman rules (hence always for
   * incoming quarks & leptons). this has to be changed if general squark
   * flavor mixing is incorporated. also, it ignores the possibility of
   * right-handed neutrinos.
   */
  {
    /* I might as well keep some notes here:
     * Rosiek's paper with all the MSSM Feynman rules is great, but it has
     * a rather annoying convention on the mixing matrices, in that he uses
     * a different convention to most others, most notably the convention
     * used in the SLHA format.
     * neutralinos:
     * Rosiek's Z_N(j,k) is NMIX(k,j) if neutralino_k has positive mass in
     * the SLHA file. if it's negative in the SLHA file, then Z_N(j,k) is
     * instead -i * NMIX(k,j).
     * Rosiek's Z_+(j,k) is VMIX(k,j), & Rosiek's Z_-(j,k) is UMIX(k,j).
     * in the absence of CP violation, the U & V mixing matrices can be
     * chosen to be real.
     */

    // by the way, if we get a gluino, we assume it has the same handedness as
    // if we got a sbottom_1.

    if( CppSLHA::PDG_code::gluino == sfermionCode )
    {
      return quarkOrLeptonLeftHandedness( CppSLHA::PDG_code::sbottom_one,
                                          ewinoCode );
    }
    else if( ( CppSLHA::PDG_code::chargino_one == ewinoCode )
             ||
             ( CppSLHA::PDG_code::chargino_one == -ewinoCode )
             ||
             ( CppSLHA::PDG_code::chargino_two == ewinoCode )
             ||
             ( CppSLHA::PDG_code::chargino_two == -ewinoCode ) )
    {
      // if it's an up-type squark, we assume that it's a positive chargino;
      // if it's a down-type squark, we assume that it's a negative chargino:
      if( ( CppSLHA::PDG_code::sup_R == sfermionCode )
          ||
          ( CppSLHA::PDG_code::scharm_R == sfermionCode )
          ||
          ( CppSLHA::PDG_code::sdown_R == sfermionCode )
          ||
          ( CppSLHA::PDG_code::sstrange_R == sfermionCode )
          ||
          ( CppSLHA::PDG_code::selectron_L == sfermionCode )
          ||
          ( CppSLHA::PDG_code::smuon_L == sfermionCode )
          ||
          ( CppSLHA::PDG_code::stau_one == sfermionCode )
          ||
          ( CppSLHA::PDG_code::selectron_R == sfermionCode )
          ||
          ( CppSLHA::PDG_code::smuon_R == sfermionCode )
          ||
          ( CppSLHA::PDG_code::stau_two == sfermionCode )
          ||
          ( CppSLHA::PDG_code::electron_sneutrino_R == sfermionCode )
          ||
          ( CppSLHA::PDG_code::muon_sneutrino_R == sfermionCode )
          ||
          ( CppSLHA::PDG_code::tau_sneutrino_R == sfermionCode ) )
      {
        return 1.0;
        /* the right-chiral squarks only interact with the Higgsino component
         * of charginos, so produce only *left*-handed quarks, while charged
         * sleptons produce only neutrinos in their interactions with
         * charginos, which are assumed to always be left-handed. (right-chiral
         * sneutrinos also only would interact with the Higgsino components.)
         */
      }
      else
      {
        double leftCoupling;
        double rightCoupling;
        double leftSquared;
        int whichChargino = 1;
        if( CppSLHA::PDG_code::chargino_two == ewinoCode )
        {
          whichChargino = 2;
        }
        if( ( NULL == umixPointer )
            ||
            ( NULL == vmixPointer )
            ||
            ( NULL == ydPointer )
            ||
            ( NULL == yuPointer )
            ||
            ( NULL == sbotmixPointer )
            ||
            ( NULL == stopmixPointer )
            ||
            ( NULL == staumixPointer ) )
        {
          if( isVerbose() )
          {
            std::cout
            << std::endl << "LHC-FASER::warning!"
            << " input_handler::quark_or_lepton_left_handedness"
            << " couldn't find at least one of UMIX, VMIX, YD, YU,"
            << " SBOTMIX, STOPMIX or STAUMIX. returning "
            << CppSLHA::CppSLHA_global::really_wrong_value_string;
            std::cout << std::endl;
          }
          return CppSLHA::CppSLHA_global::really_wrong_value;
        }
        else
        {
          if( CppSLHA::PDG_code::sup_L == sfermionCode )
          {
            leftCoupling
            = ( getWeakGaugeCoupling() * vmixPointer->get_entry( whichChargino,
                                                      1 ) );
            // actually the Feynman rule has an additional minus sign, but it
            // doesn't matter because this is getting squared.
            rightCoupling = ( ydPointer->get_entry( 1,
                                                    1 )
                              * umixPointer->get_entry( whichChargino,
                                                        2 ) );
            // actually, the Feynman rule has an additional minus sign, but it
            // doesn't matter because this is getting squared.
          }
          else if( CppSLHA::PDG_code::scharm_L == sfermionCode )
          {
            leftCoupling
            = ( getWeakGaugeCoupling() * vmixPointer->get_entry( whichChargino,
                                                      1 ) );
            // actually, the Feynman rule has an additional minus sign, but it
            // doesn't matter because this is getting squared.
            rightCoupling = ( ydPointer->get_entry( 2,
                                                    2 )
                              * umixPointer->get_entry( whichChargino,
                                                        2 ) );
            // actually, the Feynman rule has an additional minus sign, but it
            // doesn't matter because this is getting squared.
          }
          else if( CppSLHA::PDG_code::stop_one == sfermionCode )
          {
            leftCoupling = ( yuPointer->get_entry( 3,
                                                   3 )
                             * stopmixPointer->get_entry( 1,
                                                          2 )
                             * vmixPointer->get_entry( whichChargino,
                                                       2 )
                             - getWeakGaugeCoupling()
                               * stopmixPointer->get_entry( 1,
                                                            1 )
                               * vmixPointer->get_entry( whichChargino,
                                                         1 ) );
            // the (relative) signs of the terms do matter in this case.
            rightCoupling = ( ydPointer->get_entry( 3,
                                                    3 )
                              * stopmixPointer->get_entry( 1,
                                                           1 )
                              * umixPointer->get_entry( whichChargino,
                                                        2 ) );
            // actually, the Feynman rule has an additional minus sign, but it
            // doesn't matter because this is getting squared.
          }
          else if( CppSLHA::PDG_code::stop_two == sfermionCode )
          {
            leftCoupling = ( yuPointer->get_entry( 3,
                                                   3 )
                             * stopmixPointer->get_entry( 2,
                                                          2 )
                             * vmixPointer->get_entry( whichChargino,
                                                       2 )
                             - getWeakGaugeCoupling()
                               * stopmixPointer->get_entry( 2,
                                                            1 )
                               * vmixPointer->get_entry( whichChargino,
                                                         1 ) );
            // the (relative) signs of the terms do matter in this case.
            rightCoupling = ( ydPointer->get_entry( 3,
                                                    3 )
                              * stopmixPointer->get_entry( 2,
                                                           1 )
                              * umixPointer->get_entry( whichChargino,
                                                        2 ) );
            // actually, the Feynman rule has an additional minus sign, but it
            // doesn't matter because this is getting squared.
          }
          else if( CppSLHA::PDG_code::sdown_L == sfermionCode )
          {
            leftCoupling
            = ( getWeakGaugeCoupling() * umixPointer->get_entry( whichChargino,
                                                      1 ) );
            // actually the Feynman rule has an additional minus sign, but it
            // doesn't matter because this is getting squared.
            rightCoupling = ( yuPointer->get_entry( 1,
                                                    1 )
                              * vmixPointer->get_entry( whichChargino,
                                                        2 ) );
            // actually, the Feynman rule has an additional minus sign, but it
            // doesn't matter because this is getting squared.
          }
          else if( CppSLHA::PDG_code::sstrange_L == sfermionCode )
          {
            leftCoupling
            = ( getWeakGaugeCoupling() * umixPointer->get_entry( whichChargino,
                                                      1 ) );
            // actually, the Feynman rule has an additional minus sign, but it
            // doesn't matter because this is getting squared.
            rightCoupling = ( yuPointer->get_entry( 2,
                                                    2 )
                              * vmixPointer->get_entry( whichChargino,
                                                        2 ) );
            // actually, the Feynman rule has an additional minus sign, but it
            // doesn't matter because this is getting squared.
          }
          else if( CppSLHA::PDG_code::sbottom_one == sfermionCode )
          {
            leftCoupling = ( ydPointer->get_entry( 3,
                                                   3 )
                             * sbotmixPointer->get_entry( 1,
                                                          2 )
                             * umixPointer->get_entry( whichChargino,
                                                       2 )
                             - getWeakGaugeCoupling()
                               * sbotmixPointer->get_entry( 1,
                                                            1 )
                               * umixPointer->get_entry( whichChargino,
                                                         1 ) );
            // the (relative) signs of the terms do matter in this case.
            rightCoupling = ( yuPointer->get_entry( 3,
                                                    3 )
                              * sbotmixPointer->get_entry( 1,
                                                           1 )
                              * vmixPointer->get_entry( whichChargino,
                                                        2 ) );
            // actually, the Feynman rule has an additional minus sign, but it
            // doesn't matter because this is getting squared.
          }
          else if( CppSLHA::PDG_code::sbottom_two == sfermionCode )
          {
            leftCoupling = ( ydPointer->get_entry( 3,
                                                   3 )
                             * sbotmixPointer->get_entry( 2,
                                                          2 )
                             * umixPointer->get_entry( whichChargino,
                                                       2 )
                             - getWeakGaugeCoupling()
                               * sbotmixPointer->get_entry( 2,
                                                            1 )
                               * umixPointer->get_entry( whichChargino,
                                                         1 ) );
            // the (relative) signs of the terms do matter in this case.
            rightCoupling = ( yuPointer->get_entry( 3,
                                                    3 )
                              * sbotmixPointer->get_entry( 2,
                                                           1 )
                              * vmixPointer->get_entry( whichChargino,
                                                        2 ) );
            // actually, the Feynman rule has an additional minus sign, but it
            // doesn't matter because this is getting squared.
          }
          else if( CppSLHA::PDG_code::electron_sneutrino_L == sfermionCode )
          {
            leftCoupling
            = ( getWeakGaugeCoupling() * vmixPointer->get_entry( whichChargino,
                                                      1 ) );
            // actually the Feynman rule has an additional minus sign, but it
            // doesn't matter because this is getting squared.
            rightCoupling = ( yePointer->get_entry( 1,
                                                    1 )
                              * umixPointer->get_entry( whichChargino,
                                                        2 ) );
            // actually, the Feynman rule has an additional minus sign, but it
            // doesn't matter because this is getting squared.
          }
          else if( CppSLHA::PDG_code::muon_sneutrino_L == sfermionCode )
          {
            leftCoupling
            = ( getWeakGaugeCoupling() * vmixPointer->get_entry( whichChargino,
                                                      1 ) );
            // actually, the Feynman rule has an additional minus sign, but it
            // doesn't matter because this is getting squared.
            rightCoupling = ( yePointer->get_entry( 2,
                                                    2 )
                              * umixPointer->get_entry( whichChargino,
                                                        2 ) );
            // actually, the Feynman rule has an additional minus sign, but it
            // doesn't matter because this is getting squared.
          }
          else if( CppSLHA::PDG_code::tau_sneutrino_L == sfermionCode )
          {
            leftCoupling
            = ( getWeakGaugeCoupling() * vmixPointer->get_entry( whichChargino,
                                                      1 ) );
            // actually, the Feynman rule has an additional minus sign, but it
            // doesn't matter because this is getting squared.
            rightCoupling = ( yePointer->get_entry( 3,
                                                    3 )
                              * umixPointer->get_entry( whichChargino,
                                                        2 ) );
            // actually, the Feynman rule has an additional minus sign, but it
            // doesn't matter because this is getting squared.
          }
          else
          {
            if( isVerbose() )
            {
              std::cout
              << std::endl << "LHC-FASER::warning!"
              << " input_handler::quark_or_lepton_left_handedness"
              << " didn't recognize the sfermion PDG code \""
              << sfermionCode << "\". returning "
              << CppSLHA::CppSLHA_global::really_wrong_value_string;
              std::cout << std::endl;
            }
            leftSquared = CppSLHA::CppSLHA_global::really_wrong_value;
            rightSquared = CppSLHA::CppSLHA_global::really_wrong_value;
          }
          leftSquared = ( leftCoupling * leftCoupling );
          return ( leftSquared
                   / ( leftSquared + rightCoupling * rightCoupling ) );
        }
      }  // end of if-else over the types of sfermion.
    }  // end of if it's a chargino.
    else
      // otherwise we assume that it's a neutralino...
    {
      int whichNeutralino = 1;
      // we start by assuming that it's neutralino_1.
      if( CppSLHA::PDG_code::neutralino_two == ewinoCode )
      {
        whichNeutralino = 2;
      }
      else if( CppSLHA::PDG_code::neutralino_three == ewinoCode )
      {
        whichNeutralino = 3;
      }
      else if( CppSLHA::PDG_code::neutralino_three == ewinoCode )
      {
        whichNeutralino = 4;
      }
      double leftCoupling;
      double rightCoupling;
      double leftSquared;
      if( ( NULL == nmixPointer )
          ||
          ( NULL == ydPointer )
          ||
          ( NULL == yuPointer )
          ||
          ( NULL == sbotmixPointer )
          ||
          ( NULL == stopmixPointer )
          ||
          ( NULL == staumixPointer ) )
      {
        if( isVerbose() )
        {
          std::cout
          << std::endl << "LHC-FASER::warning!"
          << " input_handler::quark_or_lepton_left_handedness"
          << " couldn't find at least one of NMIX, YD, YU,"
          << " SBOTMIX, STOPMIX or STAUMIX. returning "
          << CppSLHA::CppSLHA_global::really_wrong_value_string;
          std::cout << std::endl;
        }
        return CppSLHA::CppSLHA_global::really_wrong_value;
      }
      else
      {
        if( CppSLHA::PDG_code::sdown_L == sfermionCode )
        {
          leftCoupling
          = ( M_SQRT1_2 * ( getWeakGaugeCoupling()
                            * nmixPointer->get_entry( whichNeutralino,
                                                      2 )
                            - ( getHyperchargeGaugeCoupling() / 3.0 )
                                * nmixPointer->get_entry( whichNeutralino,
                                                          1 ) ) );
          // by the way, M_SQRT1_2 is 1 over the square root of 2, defined in
          // math.h.
          rightCoupling = ( ydPointer->get_entry( 1,
                                                  1 )
                            * nmixPointer->get_entry( whichNeutralino,
                                                      3 ) );
        }
        else if( CppSLHA::PDG_code::sdown_R == sfermionCode )
        {
          leftCoupling = ( ydPointer->get_entry( 1,
                                                 1 )
                           * nmixPointer->get_entry( whichNeutralino,
                                                     3 ) );
          rightCoupling = ( M_SQRT1_2 * ( getHyperchargeGaugeCoupling() / 1.5 )
                            * nmixPointer->get_entry( whichNeutralino,
                                                      1 ) );
          // actually, the Feynman rule has an additional minus sign, but it
          // doesn't matter because this is getting squared.
        }
        else if( CppSLHA::PDG_code::sup_L == sfermionCode )
        {
          leftCoupling
          = ( M_SQRT1_2 * ( getWeakGaugeCoupling()
                            * nmixPointer->get_entry( whichNeutralino,
                                                      2 )
                            - ( getHyperchargeGaugeCoupling() / 3.0 )
                              * nmixPointer->get_entry( whichNeutralino,
                                                        1 ) ) );
          rightCoupling = ( yuPointer->get_entry( 1,
                                                  1 )
                            * nmixPointer->get_entry( whichNeutralino,
                                                      4 ) );
          // actually, the Feynman rule has an additional minus sign, but it
          // doesn't matter because this is getting squared.
        }
        else if( CppSLHA::PDG_code::sup_R == sfermionCode )
        {
          leftCoupling = ( yuPointer->get_entry( 1,
                                                 1 )
                           * nmixPointer->get_entry( whichNeutralino,
                                                     4 ) );
          // actually, the Feynman rule has an additional minus sign, but it
          // doesn't matter because this is getting squared.
          rightCoupling = ( M_SQRT1_2 * ( getHyperchargeGaugeCoupling() / 0.75 )
                            * nmixPointer->get_entry( whichNeutralino,
                                                      1 ) );
        }
        else if( CppSLHA::PDG_code::sstrange_L == sfermionCode )
        {
          leftCoupling
          = ( M_SQRT1_2 * ( getWeakGaugeCoupling()
                            * nmixPointer->get_entry( whichNeutralino,
                                                      2 )
                            - ( getHyperchargeGaugeCoupling() / 3.0 )
                              * nmixPointer->get_entry( whichNeutralino,
                                                        1 ) ) );
          rightCoupling = ( ydPointer->get_entry( 2,
                                                  2 )
                            * nmixPointer->get_entry( whichNeutralino,
                                                      3 ) );
        }
        else if( CppSLHA::PDG_code::sstrange_R == sfermionCode )
        {
          leftCoupling = ( ydPointer->get_entry( 2,
                                                 2 )
                           * nmixPointer->get_entry( whichNeutralino,
                                                     3 ) );
          rightCoupling = ( M_SQRT1_2 * ( getHyperchargeGaugeCoupling() / 1.5 )
                            * nmixPointer->get_entry( whichNeutralino,
                                                      1 ) );
          // actually, the Feynman rule has an additional minus sign, but it
          // doesn't matter because this is getting squared.
        }
        else if( CppSLHA::PDG_code::scharm_L == sfermionCode )
        {
          leftCoupling
          = ( M_SQRT1_2 * ( getWeakGaugeCoupling()
                            * nmixPointer->get_entry( whichNeutralino,
                                                      2 )
                            - ( getHyperchargeGaugeCoupling() / 3.0 )
                              * nmixPointer->get_entry( whichNeutralino,
                                                        1 ) ) );
          rightCoupling = ( yuPointer->get_entry( 2,
                                                  2 )
                            * nmixPointer->get_entry( whichNeutralino,
                                                      4 ) );
          // actually, the Feynman rule has an additional minus sign, but it
          // doesn't matter because this is getting squared.
        }
        else if( CppSLHA::PDG_code::scharm_R == sfermionCode )
        {
          leftCoupling = ( yuPointer->get_entry( 2,
                                                 2 )
                           * nmixPointer->get_entry( whichNeutralino,
                                                     4 ) );
          // actually, the Feynman rule has an additional minus sign, but it
          // doesn't matter because this is getting squared.
          rightCoupling = ( M_SQRT1_2 * ( getHyperchargeGaugeCoupling() / 0.75 )
                            * nmixPointer->get_entry( whichNeutralino,
                                                      1 ) );
        }
        else if( CppSLHA::PDG_code::sbottom_one == sfermionCode )
        {
          leftCoupling
          = ( ( M_SQRT1_2
                * ( getWeakGaugeCoupling() * nmixPointer->get_entry( whichNeutralino,
                                                          2 )
                    - ( getHyperchargeGaugeCoupling() / 3.0 )
                      * nmixPointer->get_entry( whichNeutralino,
                                                1 ) ) )
                      * sbotmixPointer->get_entry( 1,
                                                   1 )
                    + ( ydPointer->get_entry( 3,
                                              3 )
                        * nmixPointer->get_entry( whichNeutralino,
                                                  3 ) )
                        * sbotmixPointer->get_entry( 1,
                                                     2 ) );
          // the (relative) signs of the terms do matter in this case.
          rightCoupling = ( ( ydPointer->get_entry( 3,
                                                    3 )
                              * nmixPointer->get_entry( whichNeutralino,
                                                        3 ) )
                              * sbotmixPointer->get_entry( 1,
                                                           1 )
                              - ( M_SQRT1_2 * ( getHyperchargeGaugeCoupling() / 1.5 )
                                  * nmixPointer->get_entry( whichNeutralino,
                                                            1 ) )
                                  * sbotmixPointer->get_entry( 1,
                                                               2 ) );
          // the (relative) signs of the terms do matter in this case.
        }
        else if( CppSLHA::PDG_code::sbottom_two == sfermionCode )
        {
          leftCoupling
          = ( ( M_SQRT1_2
                * ( getWeakGaugeCoupling() * nmixPointer->get_entry( whichNeutralino,
                                                          2 )
                    - ( getHyperchargeGaugeCoupling()  / 3.0 )
                      * nmixPointer->get_entry( whichNeutralino,
                                                1 ) ) )
              * sbotmixPointer->get_entry( 2,
                                           1 )
              + ( ydPointer->get_entry( 3,
                                        3 )
                  * nmixPointer->get_entry( whichNeutralino,
                                            3 ) )
                  * sbotmixPointer->get_entry( 2,
                                               2 ) );
          // the (relative) signs of the terms do matter in this case.
          rightCoupling = ( ( ydPointer->get_entry( 3,
                                                    3 )
                              * nmixPointer->get_entry( whichNeutralino,
                                                        3 ) )
                              * sbotmixPointer->get_entry( 2,
                                                           1 )
                              - ( M_SQRT1_2 * ( getHyperchargeGaugeCoupling() / 1.5 )
                                  * nmixPointer->get_entry( whichNeutralino,
                                                            1 ) )
                                  * sbotmixPointer->get_entry( 2,
                                                               2 ) );
          // the (relative) signs of the terms do matter in this case.
        }
        else if( CppSLHA::PDG_code::stop_one == sfermionCode )
        {
          leftCoupling
          = ( ( M_SQRT1_2
                * ( getWeakGaugeCoupling() * nmixPointer->get_entry( whichNeutralino,
                                                          2 )
                    - ( getHyperchargeGaugeCoupling() / 3.0 )
                      * nmixPointer->get_entry( whichNeutralino,
                                                1 ) ) )
              * stopmixPointer->get_entry( 1,
                                           1 )
              - ( yuPointer->get_entry( 3,
                                        3 )
                  * nmixPointer->get_entry( whichNeutralino,
                                            4 ) )
                  * stopmixPointer->get_entry( 1,
                                               2 ) );
          // the (relative) signs of the terms do matter in this case.
          rightCoupling = ( ( M_SQRT1_2 * ( getHyperchargeGaugeCoupling() / 0.75 )
                              * nmixPointer->get_entry( whichNeutralino,
                                                        1 ) )
                              * stopmixPointer->get_entry( 1,
                                                           2 )
                              - ( yuPointer->get_entry( 3,
                                                        3 )
                                  * nmixPointer->get_entry( whichNeutralino,
                                                            4 ) )
                                  * stopmixPointer->get_entry( 1,
                                                               1 ) );
          // the (relative) signs of the terms do matter in this case.
        }
        else if( CppSLHA::PDG_code::stop_two == sfermionCode )
        {
          leftCoupling
          = ( ( M_SQRT1_2
                * ( getWeakGaugeCoupling() * nmixPointer->get_entry( whichNeutralino,
                                                          2 )
                    - ( getHyperchargeGaugeCoupling() / 3.0 )
                      * nmixPointer->get_entry( whichNeutralino,
                                                1 ) ) )
                * stopmixPointer->get_entry( 2,
                                             1 )
                - ( yuPointer->get_entry( 3,
                                          3 )
                    * nmixPointer->get_entry( whichNeutralino,
                                              4 ) )
                    * stopmixPointer->get_entry( 2,
                                                 2 ) );
          // the (relative) signs of the terms do matter in this case.
          rightCoupling = ( ( M_SQRT1_2 * ( getHyperchargeGaugeCoupling() / 0.75 )
                              * nmixPointer->get_entry( whichNeutralino,
                                                        1 ) )
                              * stopmixPointer->get_entry( 2,
                                                           2 )
                              - ( yuPointer->get_entry( 3,
                                                        3 )
                                  * nmixPointer->get_entry( whichNeutralino,
                                                            4 ) )
                                  * stopmixPointer->get_entry( 2,
                                                               1 ) );
          // the (relative) signs of the terms do matter in this case.
        }
        else if( CppSLHA::PDG_code::selectron_L == sfermionCode )
        {
          leftCoupling
          = ( M_SQRT1_2 * ( getWeakGaugeCoupling()
                            * nmixPointer->get_entry( whichNeutralino,
                                                      2 )
                            + getHyperchargeGaugeCoupling()
                              * nmixPointer->get_entry( whichNeutralino,
                                                        1 ) ) );
          rightCoupling = ( yePointer->get_entry( 1,
                                                  1 )
                            * nmixPointer->get_entry( whichNeutralino,
                                                      3 ) );
        }
        else if( CppSLHA::PDG_code::selectron_R == sfermionCode )
        {
          leftCoupling = ( yePointer->get_entry( 1,
                                                 1 )
                           * nmixPointer->get_entry( whichNeutralino,
                                                     3 ) );
          rightCoupling = ( M_SQRT1_2 * ( getHyperchargeGaugeCoupling() * 2.0 )
                            * nmixPointer->get_entry( whichNeutralino,
                                                          1 ) );
          // actually, the Feynman rule has an additional minus sign, but it
          // doesn't matter because this is getting squared.
        }
        else if( CppSLHA::PDG_code::smuon_L == sfermionCode )
        {
          leftCoupling
          = ( M_SQRT1_2 * ( getWeakGaugeCoupling()
                            * nmixPointer->get_entry( whichNeutralino,
                                                      2 )
                            + getHyperchargeGaugeCoupling()
                              * nmixPointer->get_entry( whichNeutralino,
                                                        1 ) ) );
          rightCoupling = ( yePointer->get_entry( 2,
                                                  2 )
                            * nmixPointer->get_entry( whichNeutralino,
                                                      3 ) );
        }
        else if( CppSLHA::PDG_code::smuon_R == sfermionCode )
        {
          leftCoupling = ( yePointer->get_entry( 2,
                                                 2 )
                           * nmixPointer->get_entry( whichNeutralino,
                                                     3 ) );
          rightCoupling = ( M_SQRT1_2 * ( getHyperchargeGaugeCoupling() * 2.0 )
                            * nmixPointer->get_entry( whichNeutralino,
                                                      1 ) );
          // actually, the Feynman rule has an additional minus sign, but it
          // doesn't matter because this is getting squared.
        }
        else if( CppSLHA::PDG_code::stau_one == sfermionCode )
        {
          leftCoupling
          = ( ( M_SQRT1_2
                * ( getWeakGaugeCoupling() * nmixPointer->get_entry( whichNeutralino,
                                                          2 )
                    + getHyperchargeGaugeCoupling()
                      * nmixPointer->get_entry( whichNeutralino,
                                                1 ) ) )
              * staumixPointer->get_entry( 1,
                                           1 )
              + ( yePointer->get_entry( 3,
                                        3 )
                  * nmixPointer->get_entry( whichNeutralino,
                                            3 ) )
                  * staumixPointer->get_entry( 1,
                                               2 ) );
          // the (relative) signs of the terms do matter in this case.
          rightCoupling = ( ( yePointer->get_entry( 3,
                                                    3 )
                              * nmixPointer->get_entry( whichNeutralino,
                                                        3 ) )
                              * staumixPointer->get_entry( 1,
                                                           1 )
                            - ( M_SQRT1_2 * ( getHyperchargeGaugeCoupling() * 2.0 )
                                * nmixPointer->get_entry( whichNeutralino,
                                                          1 ) )
                                * staumixPointer->get_entry( 1,
                                                             2 ) );
          // the (relative) signs of the terms do matter in this case.
        }
        else if( CppSLHA::PDG_code::stau_two == sfermionCode )
        {
          leftCoupling
          = ( ( M_SQRT1_2
                * ( getWeakGaugeCoupling() * nmixPointer->get_entry( whichNeutralino,
                                                          2 )
                    + getHyperchargeGaugeCoupling()
                      * nmixPointer->get_entry( whichNeutralino,
                                                1 ) ) )
              * staumixPointer->get_entry( 2,
                                           1 )
              + ( yePointer->get_entry( 3,
                                        3 )
                  * nmixPointer->get_entry( whichNeutralino,
                                            3 ) )
                  * staumixPointer->get_entry( 2,
                                               2 ) );
          // the (relative) signs of the terms do matter in this case.
          rightCoupling = ( ( yePointer->get_entry( 3,
                                                    3 )
                              * nmixPointer->get_entry( whichNeutralino,
                                                        3 ) )
                            * staumixPointer->get_entry( 2,
                                                         1 )
                            - ( M_SQRT1_2 * ( getHyperchargeGaugeCoupling() * 2.0 )
                                * nmixPointer->get_entry( whichNeutralino,
                                                          1 ) )
                              * staumixPointer->get_entry( 2,
                                                           2 ) );
          // the (relative) signs of the terms do matter in this case.
        }
        else
        {
          if( isVerbose() )
          {
            std::cout
            << std::endl << "LHC-FASER::warning!"
            << " input_handler::quark_or_lepton_left_handedness"
            << " didn't recognize the sfermion PDG code \""
            << sfermionCode << "\". returning "
            << CppSLHA::CppSLHA_global::really_wrong_value_string;
            std::cout << std::endl;
          }
          leftSquared = CppSLHA::CppSLHA_global::really_wrong_value;
          rightSquared = CppSLHA::CppSLHA_global::really_wrong_value;
        }
      }
      leftSquared = ( leftCoupling * leftCoupling );
      return ( leftSquared
               / ( leftSquared + rightCoupling * rightCoupling ) );
    }  // end of if it's a neutralino.
  }
}  // end of LHC_FASER namespace.

