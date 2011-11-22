/*
 * LHC-FASER_input_handling_stuff.hpp
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

#ifndef LHC_FASER_INPUT_HANDLING_STUFF_HPP_
#define LHC_FASER_INPUT_HANDLING_STUFF_HPP_

#include "LHC-FASER_global_stuff.hpp"
#include "LHC-FASER_base_lepton_distribution_stuff.hpp"
#include "LHC-FASER_sparticle_decay_stuff.hpp"

namespace LHC_FASER
{
  /* this class holds the various gauge couplings which are re-calculated for
   * every new point. it is kept separate from the inputHandler class for ease
   * of (abusing) const correctness.
   */
  class updateDependentGaugeInputs : public getsReadiedForNewPoint
  {
  public:
    updateDependentGaugeInputs( CppSLHA::SLHA_BLOCK const* smInputs,
                                particlePointer const wPlusPointer,
                                particlePointer const zPointer,
                                readierForNewPoint* const readierPointer );
    ~updateDependentGaugeInputs();

    double
    getElectromagneticFineStructureConstant();
    double
    getHyperchargeGaugeCoupling();
    double
    getWeakGaugeCoupling();
    double
    getColorGaugeCoupling();
    double
    getWeakSine();
    double
    getWeakCosine();

  protected:
    CppSLHA::SLHA_BLOCK const* smInputs;
    particlePointer const wPlusPointer;
    particlePointer const zPointer;
    // the gauge couplings:
    double electromagneticFineStructureConstant;
    double hyperchargeGaugeCoupling;
    // the gauge coupling for U(1) at m_Z.
    double weakGaugeCoupling;
    // the gauge coupling for SU(2) at m_Z.
    double colorGaugeCoupling;
    // the gauge coupling for SU(3) at m_Z.
    double weakSine;
    // the sine of the weak mixing angle.
    double weakCosine;
    // the cosine of the weak mixing angle.

    void
    calculateGaugeCouplings();
  };


  /* this class is to hold frequently-accessed absolute masses which are
   * re-calculated for every new point. it is kept separate from the
   * inputHandler class for ease of (abusing) const correctness.
   */
  class updateDependentAbsoluteMasses : public getsReadiedForNewPoint
  {
  public:
    updateDependentAbsoluteMasses( particlePointer const gluinoPointer,
                                   readierForNewPoint* const readierPointer );
    ~updateDependentAbsoluteMasses();

    double
    getGluinoMass();

  protected:
    particlePointer const gluinoPointer;
    double gluinoMass;
  };


  /* this class is to hold average absolute masses which are re-calculated for
   * every new point. it is kept separate from the inputHandler class for ease
   * of (abusing) const correctness.
   */
  class updateDependentAverageMass : public getsReadiedForNewPoint
  {
  public:
    updateDependentAverageMass( particleVectorPointer const setPointer,
                                readierForNewPoint* const readierPointer );
    ~updateDependentAverageMass();

    double
    getAverageMass();

  protected:
    particleVectorPointer const setPointer;
    double averageMass;
  };


  // this class holds a set of pointers to various objects for ease of
  // reference to particles.
  class inputHandler
  {
  public:
    inputHandler( CppSLHA::CppSLHA0 const* const cppSlhaPointer,
                  CppSLHA::EW_scale_spectrum const* const spectrumPointer,
                  std::string const pathToGrids,
                  readierForNewPoint* const readierPointer );
    ~inputHandler();

    bool
    isVerbose()
    const;
    /* if this is true, lots of warnings will be printed. if false, the code
     * should stay quiet without warning of the impending segmentation
     * faults...
     */
    void
    setVerbosity( bool isVerbose );
    bool
    usingNlo()
    const;
    void
    setUsingNloFlag( bool usingNloFlag );
    std::string const*
    getPathToGrids()
    const;
    std::string*
    getPathToCrossSectionGrids();
    std::string const*
    inspectPathToCrossSectionGrids()
    const;
    void
    setPathToCrossSectionGrids(
                            std::string const* const pathToCrossSectionGrids );
    std::string*
    getPathToKinematicsGrids();
    std::string const*
    inspectPathToKinematicsGrids()
    const;
    void
    setPathToKinematicsGrids( std::string const* const pathToKinematicsGrids );
    CppSLHA::CppSLHA0 const*
    getCppSlha()
    const;
    CppSLHA::EW_scale_spectrum const*
    getSpectrum()
    const;
    double
    getElectromagneticFineStructureConstant()
    const;
    double
    getHyperchargeGaugeCoupling()
    const;
    double
    getWeakGaugeCoupling()
    const;
    double
    getColorGaugeCoupling()
    const;
    double
    getWeakSine()
    const;
    double
    getWeakCosine()
    const;
    readierForNewPoint*
    getReadier()
    const;
    particlePointer
    getDown()
    const;
    particlePointer
    getUp()
    const;
    particlePointer
    getStrange()
    const;
    particlePointer
    getCharm()
    const;
    particlePointer
    getBottom()
    const;
    particlePointer
    getTop()
    const;
    particlePointer
    getElectron()
    const;
    particlePointer
    getElectronNeutrino()
    const;
    particlePointer
    getMuon()
    const;
    particlePointer
    getMuonNeutrino()
    const;
    particlePointer
    getTauLepton()
    const;
    particlePointer
    getTauNeutrino()
    const;
    particlePointer
    getGluon()
    const;
    particlePointer
    getPhoton()
    const;
    particlePointer
    getZ()
    const;
    particlePointer
    getWPlus()
    const;
    particlePointer
    getLightNeutralEwsbScalar()
    const;
    particlePointer
    getHeavyNeutralEwsbScalar()
    const;
    particlePointer
    getNeutralEwsbPseudoscalar()
    const;
    particlePointer
    getChargedEwsbScalar()
    const;
    particlePointer
    getGluino()
    const;
    double
    getGluinoMass()
    const;
    particlePointer
    getSdownL()
    const;
    particlePointer
    getSupL()
    const;
    particlePointer
    getSstrangeL()
    const;
    particlePointer
    getScharmL()
    const;
    particlePointer
    getSbottomOne()
    const;
    particlePointer
    getStopOne()
    const;
    particlePointer
    getSdownR()
    const;
    particlePointer
    getSupR()
    const;
    particlePointer
    getSstrangeR()
    const;
    particlePointer
    getScharmR()
    const;
    particlePointer
    getSbottomTwo()
    const;
    particlePointer
    getStopTwo()
    const;
    particlePointer
    getSelectronL()
    const;
    particlePointer
    getElectronSneutrinoL()
    const;
    particlePointer
    getSmuonL()
    const;
    particlePointer
    getMuonSneutrinoL()
    const;
    particlePointer
    getStauOne()
    const;
    particlePointer
    getTauSneutrinoL()
    const;
    particlePointer
    getSelectronR()
    const;
    particlePointer
    getElectronSneutrinoR()
    const;
    particlePointer
    getSmuonR()
    const;
    particlePointer
    getMuonSneutrinoR()
    const;
    particlePointer
    getStauTwo()
    const;
    particlePointer
    getTauSneutrinoR()
    const;
    particlePointer
    getNeutralinoOne()
    const;
    particlePointer
    getNeutralinoTwo()
    const;
    particlePointer
    getNeutralinoThree()
    const;
    particlePointer
    getNeutralinoFour()
    const;
    particlePointer
    getCharginoOne()
    const;
    particlePointer
    getCharginoTwo()
    const;
    bool
    isIn( int const pdgCode,
          particleVectorPointer const particleSet )
    const;
    /* this checks to see if the particle with the absolute value of the given
     * PDG particle code is in the given collection. (hence it checks to see if
     * the given particle or its antiparticle are in the collection.)
     */
    exclusiveBrCalculator*
    getExclusiveBrCalculator( particlePointer const decayerParticle,
                              particlePointer const productParticle,
                              bool const productIsNotAntiparticle,
                              std::list< int > const* const exclusionList )
    const;
    std::list< int > const*
    getEmptyList()
    const;
    std::list< int > const*
    getJets()
    const;
    std::list< int > const*
    getJets5()
    const;
    std::list< int > const*
    getJets4()
    const;
    std::list< int > const*
    getNotInJets5()
    const;
    std::list< int > const*
    getNotInJets4()
    const;
    particleVectorPointer
    getColoredSparticles()
    const;
    particleVectorPointer
    getSquarks()
    const;
    particleVectorPointer
    getSquarks5()
    const;
    particleVectorPointer
    getSquarks4()
    const;
    particleVectorPointer
    getSdownTypes()
    const;
    particleVectorPointer
    getSupTypes()
    const;
    std::vector< signedParticleShortcutPair* > const*
    getScoloredProductionCombinations()
    const;
    particleVectorPointer
    getElectroweakinos()
    const;
    particleVectorPointer
    getUnstableElectroweakinos()
    const;
    particleVectorPointer
    getNeutralinos()
    const;
    particleVectorPointer
    getUnstableNeutralinos()
    const;
    particleVectorPointer
    getCharginos()
    const;
    particleVectorPointer
    getChargedSleptons()
    const;
    particleVectorPointer
    getSneutrinos()
    const;
    particleVectorPointer
    getNeutralEwsbScalarsAndPseudoscalars()
    const;
    particleVectorPointer
    getEwsbBosonsAndMassiveVectorBosons()
    const;
    particleVectorPointer
    getNeutralEwsbBosonsAndMassiveVectorBosons()
    const;
    particleVectorPointer
    getChargedEwsbBosonsAndMassiveVectorBosons()
    const;
    hardMuonFromTau const*
    getHardMuonFromTau()
    const;
    softMuonFromTau const*
    getSoftMuonFromTau()
    const;
    hardPionFromTau const*
    getHardPionFromTau()
    const;
    softPionFromTau const*
    getSoftPionFromTau()
    const;
    double
    getAverageSquarks4Mass()
    const;
    double
    quarkOrLeptonLeftHandedness( int const sfermionCode,
                                 int const ewinoCode );
    const
    /* this returns the square of the coupling associated with the left-handed
     * projection operator of the Feynman rule for the vertex of the given
     * sfermion with the given electroweakino & the appropriate quark, divided
     * by the sum of the squares of the couplings for the left- & right-handed
     * projection operators, based on Rosiek's Feynman rules. this has to be
     * changed if general squark flavor mixing is incorporated. also, it
     * ignores the possibility of right-handed neutrinos.
     */
    effectiveSquarkMassHolder*
    getOnShellEffectiveSquarkMass( particlePointer const  onShellSquark )
    const;
    effectiveSquarkMassHolder*
    getSquarkMinusBosonEffectiveMass( particlePointer const onShellSquark,
                                      particlePointer const onShellBoson,
                                      particlePointer const onShellEwino )
    const;
    effectiveSquarkMassHolder*
    getSquarkPlusBosonEffectiveMass( particlePointer const onShellSquark,
                                     particlePointer const onShellBoson,
                                     particlePointer const onShellEwino )
    const;

  protected:
    bool isVerbose;
    bool usingNloFlag;
    readierForNewPoint* const readierPointer;
    CppSLHA::CppSLHA0 const* const cppSlhaPointer;
    CppSLHA::EW_scale_spectrum const* const spectrumPointer;
    std::string pathToGrids;
    std::string pathToCrossSectionGrids;
    std::string pathToKinematicsGrids;

    particlePointer const downPointer;
    particlePointer const upPointer;
    particlePointer const strangePointer;
    particlePointer const charmPointer;
    particlePointer const bottomPointer;
    particlePointer const topPointer;
    particlePointer const electronPointer;
    particlePointer const electronNeutrinoPointer;
    particlePointer const muonPointer;
    particlePointer const muonNeutrinoPointer;
    particlePointer const tauLeptonPointer;
    particlePointer const tauNeutrinoPointer;
    particlePointer const gluonPointer;
    particlePointer const photonPointer;
    particlePointer const zPointer;
    particlePointer const wPlusPointer;

    particlePointer const lightNeutralEwsbScalarPointer;
    particlePointer const heavyNeutralEwsbScalarPointer;
    particlePointer const neutralEwsbPseudoscalarPointer;
    particlePointer const chargedEwsbScalarPointer;
    particlePointer const gluinoPointer;
    particlePointer const sdownLPointer;
    particlePointer const supLPointer;
    particlePointer const sstrangeLPointer;
    particlePointer const scharmLPointer;
    particlePointer const sbottomOnePointer;
    particlePointer const stopOnePointer;
    particlePointer const sdownRPointer;
    particlePointer const supRPointer;
    particlePointer const sstrangeRPointer;
    particlePointer const scharmRPointer;
    particlePointer const sbottomTwoPointer;
    particlePointer const stopTwoPointer;
    particlePointer const selectronLPointer;
    particlePointer const electronSneutrinoLPointer;
    particlePointer const smuonLPointer;
    particlePointer const muonSneutrinoLPointer;
    particlePointer const stauOnePointer;
    particlePointer const tauSneutrinoLPointer;
    particlePointer const selectronRPointer;
    particlePointer const electronSneutrinoRPointer;
    particlePointer const smuonRPointer;
    particlePointer const muonSneutrinoRPointer;
    particlePointer const stauTwoPointer;
    particlePointer const tauSneutrinoRPointer;
    particlePointer const neutralinoOnePointer;
    particlePointer const neutralinoTwoPointer;
    particlePointer const neutralinoThreePointer;
    particlePointer const neutralinoFourPointer;
    particlePointer const charginoOnePointer;
    particlePointer const charginoTwoPointer;

    // a few useful SLHA BLOCK shortcuts:
    CppSLHA::SLHA_BLOCK const* nmixPointer;
    CppSLHA::SLHA_BLOCK const* umixPointer;
    CppSLHA::SLHA_BLOCK const* vmixPointer;
    CppSLHA::SLHA_BLOCK const* sbotmixPointer;
    CppSLHA::SLHA_BLOCK const* stopmixPointer;
    CppSLHA::SLHA_BLOCK const* staumixPointer;
    CppSLHA::SLHA_BLOCK const* ydPointer;
    CppSLHA::SLHA_BLOCK const* yuPointer;
    CppSLHA::SLHA_BLOCK const* yePointer;

    updateDependentGaugeInputs* updateDependentInputValues;
    updateDependentAbsoluteMasses* updateDependentMasses;
    updateDependentAverageMass* updateDependentAverageSquarks4Mass;

    std::list< int > emptyList;
    // emptyList is there for when we don't care about excluding particles.
    /* the following std::list< int > are lists of PDG codes for particles:
     * jets refers to all the particles in the spectrum which return true for
     * counts_as_jet(), as well as their antiparticles, unless they are
     * self-conjugate.
     * notInFiveFlavorJets is just top & antitop, & notInFourFlavorJets is
     * notInFiveFlavorJets plus bottom & antibottom.
     * fiveFlavorJets is jets minus everything in notInFiveFlavorJets, likewise
     * fourFlavorJets is jets minus everything in fourFlavorJets.
     * other exotic heavy quarks should either be labelled as not counting as
     * a jet in the spectrum, or put in notInFiveFlavorJets in the same manner
     * as top.
     * hence, without additional light jet particles, fiveFlavorJets is
     * { gluon + the 5 light quarks (down, up, strange, charm, bottom)
     * + the 5 light antiquarks }, & fourFlavorJets is fiveFlavorJets with
     * bottom & antibottom removed.
     */
    std::list< int > sixFlavorJets;
    std::list< int > fiveFlavorJets;
    std::list< int > fourFlavorJets;
    std::list< int > notInFiveFlavorJets;
    std::list< int > notInFourFlavorJets;

    std::vector< particlePointer > coloredSparticles;
    /* fourFlavorSquarks refers to { sdown, sup, sstrange, scharm }, while
     * fiveFlavorSquarks is fourFlavorSquarks plus sbottom, & squarks is
     * fiveFlavorSquarks plus stop.
     */
    std::vector< particlePointer > sixFlavorSquarks;
    std::vector< particlePointer > fiveFlavorSquarks;
    std::vector< particlePointer > fourFlavorSquarks;
    std::vector< particlePointer > sdownTypes;
    std::vector< particlePointer > supTypes;
    std::vector< signedParticleShortcutPair* >
    sparticleProductionCombinations;
    std::vector< onShellSquarkMassHolder* > onShellSquarkEffectiveMasses;
    std::vector< squarkMinusBosonMassHolder* > squarkMinusBosonEffectiveMasses;
    std::vector< squarkPlusBosonMassHolder* > squarkPlusBosonEffectiveMasses;

    std::vector< particlePointer > allElectroweakinos;
    std::vector< particlePointer > unstableElectroweakinos;
    std::vector< particlePointer > allNeutralinos;
    std::vector< particlePointer > unstableNeutralinos;
    std::vector< particlePointer > allCharginos;
    std::vector< particlePointer > allChargedSleptons;
    std::vector< particlePointer > allSneutrinos;

    std::vector< particlePointer > neutralEwsbScalarsAndPseudoscalars;
    std::vector< particlePointer > ewsbBosonsAndMassiveVectorBosons;
    std::vector< particlePointer > neutralEwsbBosonsAndMassiveVectorBosons;
    std::vector< particlePointer > chargedEwsbBosonsAndMassiveVectorBosons;

    exclusiveBrHandler* exclusiveBrs;

    hardMuonFromTau hardMuonFromTauFunction;
    softMuonFromTau softMuonFromTauFunction;
    hardPionFromTau hardPionFromTauFunction;
    softPionFromTau softPionFromTauFunction;

    void
    addSquarkAsColoredSparticleAndInProduction(
                                           particlePointer const inputSquark );
    /* this adds the given squark to the list of squarks & the list of colored
     * sparticles, & also adds its relevant production combinations (with
     * a gluino, with its antisquark, with itself) to the list.
     */
  };




  // inline functions:

  inline double
  updateDependentGaugeInputs::getElectromagneticFineStructureConstant()
  {
    if( needsToPrepareForThisPoint() )
    {
      calculateGaugeCouplings();
    }
    return electromagneticFineStructureConstant;
  }

  inline double
  updateDependentGaugeInputs::getHyperchargeGaugeCoupling()
  {
    if( needsToPrepareForThisPoint() )
    {
      calculateGaugeCouplings();
    }
    return hyperchargeGaugeCoupling;
  }

  inline double
  updateDependentGaugeInputs::getWeakGaugeCoupling()
  {
    if( needsToPrepareForThisPoint() )
    {
      calculateGaugeCouplings();
    }
    return weakGaugeCoupling;
  }

  inline double
  updateDependentGaugeInputs::getColorGaugeCoupling()
  {
    if( needsToPrepareForThisPoint() )
    {
      calculateGaugeCouplings();
    }
    return colorGaugeCoupling;
  }

  inline double
  updateDependentGaugeInputs::getWeakSine()
  {
    if( needsToPrepareForThisPoint() )
    {
      calculateGaugeCouplings();
    }
    return weakSine;
  }

  inline double
  updateDependentGaugeInputs::getWeakCosine()
  {
    if( needsToPrepareForThisPoint() )
    {
      calculateGaugeCouplings();
    }
    return weakCosine;
  }



  inline double
  updateDependentAbsoluteMasses::getGluinoMass()
  {
    if( needsToPrepareForThisPoint() )
      {
        gluinoMass = gluinoPointer->get_absolute_mass();
      }
    return gluinoMass;
  }



  inline double
  updateDependentAverageMass::getAverageMass()
  {
    if( needsToPrepareForThisPoint() )
    {
      averageMass = 0.0;
      for( std::vector< particlePointer >::const_iterator
           particleIterator = setPointer->begin();
           setPointer->end() > particleIterator;
           ++particleIterator )
      {
        averageMass += (*particleIterator)->get_absolute_mass();
      }
      averageMass *= ( 1.0 / (double)setPointer->size() );
    }
    return averageMass;
  }



  inline bool
  inputHandler::isVerbose()
  const
  // if this is true, lots of warnings will be printed. if false, the code
  // should stay quiet without warning of the impending segmentation faults...
  {
    return isVerbose;
  }

  inline void
  inputHandler::setVerbosity( bool isVerbose )
  {
    this->isVerbose = isVerbose;
  }

  inline bool
  inputHandler::usingNlo()
  const
  {
    return usingNloFlag;
  }

  inline void
  inputHandler::setUsingNloFlag( bool usingNloFlag )
  {
    this->usingNloFlag = usingNloFlag;
  }

  inline std::string const*
  inputHandler::getPathToGrids()
  const
  {
    return &pathToGrids;
  }

  inline std::string*
  inputHandler::getPathToCrossSectionGrids()
  {
    return &pathToCrossSectionGrids;
  }

  inline std::string const*
  inputHandler::inspectPathToCrossSectionGrids()
  const
  {
    return &pathToCrossSectionGrids;
  }

  inline void
  inputHandler::setPathToCrossSectionGrids(
                             std::string const* const pathToCrossSectionGrids )
  {
    this->pathToCrossSectionGrids.assign( *pathToCrossSectionGrids );
  }

  inline std::string*
  inputHandler::getPathToKinematicsGrids()
  {
    return &pathToKinematicsGrids;
  }

  inline std::string const*
  inputHandler::inspectPathToKinematicsGrids()
  const
  {
    return &pathToKinematicsGrids;
  }

  inline void
  inputHandler::setPathToKinematicsGrids(
                               std::string const* const pathToKinematicsGrids )
  {
    this->pathToKinematicsGrids.assign( *pathToKinematicsGrids );
  }

  inline CppSLHA::CppSLHA0 const*
  inputHandler::getCppSlha()
  const
  {
    return cppSlhaPointer;
  }

  inline CppSLHA::EW_scale_spectrum const*
  inputHandler::getSpectrum()
  const
  {
    return spectrumPointer;
  }

  inline double
  inputHandler::getElectromagneticFineStructureConstant()
  const
  {
    return
    updateDependentInputValues->getElectromagneticFineStructureConstant();
  }

  inline double
  inputHandler::getHyperchargeGaugeCoupling()
  const
  {
    return updateDependentInputValues->getHyperchargeGaugeCoupling();
  }

  inline double
  inputHandler::getWeakGaugeCoupling()
  const
  {
    return updateDependentInputValues->getWeakGaugeCoupling();
  }

  inline double
  inputHandler::getColorGaugeCoupling()
  const
  {
    return updateDependentInputValues->getColorGaugeCoupling();
  }

  inline double
  inputHandler::getWeakSine()
  const
  {
    return updateDependentInputValues->getWeakSine();
  }

  inline double
  inputHandler::getWeakCosine()
  const
  {
    return updateDependentInputValues->getWeakCosine();
  }

  inline readierForNewPoint*
  inputHandler::getReadier()
  const
  {
    return readierPointer;
  }

  inline particlePointer
  inputHandler::getDown()
  const
  {
    return downPointer;
  }

  inline particlePointer
  inputHandler::getUp()
  const
  {
    return upPointer;
  }

  inline particlePointer
  inputHandler::getStrange()
  const
  {
    return strangePointer;
  }

  inline particlePointer
  inputHandler::getCharm()
  const
  {
    return charmPointer;
  }

  inline particlePointer
  inputHandler::getBottom()
  const
  {
    return bottomPointer;
  }

  inline particlePointer
  inputHandler::getTop()
  const
  {
    return topPointer;
  }

  inline particlePointer
  inputHandler::getElectron()
  const
  {
    return electronPointer;
  }

  inline particlePointer
  inputHandler::getElectronNeutrino()
  const
  {
    return electronNeutrinoPointer;
  }

  inline particlePointer
  inputHandler::getMuon()
  const
  {
    return muonPointer;
  }

  inline particlePointer
  inputHandler::getMuonNeutrino()
  const
  {
    return muonNeutrinoPointer;
  }

  inline particlePointer
  inputHandler::getTauLepton()
  const
  {
    return tauLeptonPointer;
  }

  inline particlePointer
  inputHandler::getTauNeutrino()
  const
  {
    return tauNeutrinoPointer;
  }

  inline particlePointer
  inputHandler::getGluon()
  const
  {
    return gluonPointer;
  }

  inline particlePointer
  inputHandler::getPhoton()
  const
  {
    return photonPointer;
  }

  inline particlePointer
  inputHandler::getZ()
  const
  {
    return zPointer;
  }

  inline particlePointer
  inputHandler::getWPlus()
  const
  {
    return wPlusPointer;
  }

  inline particlePointer
  inputHandler::getLightNeutralEwsbScalar()
  const
  {
    return lightNeutralEwsbScalarPointer;
  }

  inline particlePointer
  inputHandler::getHeavyNeutralEwsbScalar()
  const
  {
    return heavyNeutralEwsbScalarPointer;
  }

  inline particlePointer
  inputHandler::getNeutralEwsbPseudoscalar()
  const
  {
    return neutralEwsbPseudoscalarPointer;
  }

  inline particlePointer
  inputHandler::getChargedEwsbScalar()
  const
  {
    return chargedEwsbScalarPointer;
  }

  inline particlePointer
  inputHandler::getGluino()
  const
  {
    return gluinoPointer;
  }

  inline double
  inputHandler::getGluinoMass()
  const
  {
    return updateDependentMasses->getGluinoMass();
  }

  inline particlePointer
  inputHandler::getSdownL()
  const
  {
    return sdownLPointer;
  }

  inline particlePointer
  inputHandler::getSupL()
  const
  {
    return supLPointer;
  }

  inline particlePointer
  inputHandler::getSstrangeL()
  const
  {
    return sstrangeLPointer;
  }

  inline particlePointer
  inputHandler::getScharmL()
  const
  {
    return scharmLPointer;
  }

  inline particlePointer
  inputHandler::getSbottomOne()
  const
  {
    return sbottomOnePointer;
  }

  inline particlePointer
  inputHandler::getStopOne()
  const
  {
    return stopOnePointer;
  }

  inline particlePointer
  inputHandler::getSdownR()
  const
  {
    return sdownRPointer;
  }

  inline particlePointer
  inputHandler::getSupR()
  const
  {
    return supRPointer;
  }

  inline particlePointer
  inputHandler::getSstrangeR()
  const
  {
    return sstrangeRPointer;
  }

  inline particlePointer
  inputHandler::getScharmR()
  const
  {
    return scharmRPointer;
  }

  inline particlePointer
  inputHandler::getSbottomTwo()
  const
  {
    return sbottomTwoPointer;
  }

  inline particlePointer
  inputHandler::getStopTwo()
  const
  {
    return stopTwoPointer;
  }

  inline particlePointer
  inputHandler::getSelectronL()
  const
  {
    return selectronLPointer;
  }

  inline particlePointer
  inputHandler::getElectronSneutrinoL()
  const
  {
    return electronSneutrinoLPointer;
  }

  inline particlePointer
  inputHandler::getSmuonL()
  const
  {
    return smuonLPointer;
  }

  inline particlePointer
  inputHandler::getMuonSneutrinoL()
  const
  {
    return muonSneutrinoLPointer;
  }

  inline particlePointer
  inputHandler::getStauOne()
  const
  {
    return stauOnePointer;
  }

  inline particlePointer
  inputHandler::getTauSneutrinoL()
  const
  {
    return tauSneutrinoLPointer;
  }

  inline particlePointer
  inputHandler::getSelectronR()
  const
  {
    return selectronRPointer;
  }

  inline particlePointer
  inputHandler::getElectronSneutrinoR()
  const
  {
    return electronSneutrinoRPointer;
  }

  inline particlePointer
  inputHandler::getSmuonR()
  const
  {
    return smuonRPointer;
  }

  inline particlePointer
  inputHandler::getMuonSneutrinoR()
  const
  {
    return muonSneutrinoRPointer;
  }

  inline particlePointer
  inputHandler::getStauTwo()
  const
  {
    return stauTwoPointer;
  }

  inline particlePointer
  inputHandler::getTauSneutrinoR()
  const
  {
    return tauSneutrinoRPointer;
  }

  inline particlePointer
  inputHandler::getNeutralinoOne()
  const
  {
    return neutralinoOnePointer;
  }

  inline particlePointer
  inputHandler::getNeutralinoTwo()
  const
  {
    return neutralinoTwoPointer;
  }

  inline particlePointer
  inputHandler::getNeutralinoThree()
  const
  {
    return neutralinoThreePointer;
  }

  inline particlePointer
  inputHandler::getNeutralinoFour()
  const
  {
    return neutralinoFourPointer;
  }

  inline particlePointer
  inputHandler::getCharginoOne()
  const
  {
    return charginoOnePointer;
  }

  inline particlePointer
  inputHandler::getCharginoTwo()
  const
  {
    return charginoTwoPointer;
  }

  inline bool
  inputHandler::isIn( int const pdgCode,
                      particleVectorPointer const particleSet )
  const
  /* this checks to see if the particle with the absolute value of the given
   * PDG particle code is in the given collection. (hence it checks to see if
   * the given particle or its antiparticle are in the collection.)
   */
  {
    bool returnBool( false );
    for( std::vector< particlePointer >::const_iterator
         setIterator( particleSet->begin() );
         particleSet->end() > setIterator;
         ++setIterator )
    {
      if( ( pdgCode == (*setIterator)->get_PDG_code() )
          ||
          ( pdgCode == -(*setIterator)->get_PDG_code() ) )
      {
        returnBool = true;
      }
    }
    return returnBool;
  }

  inline exclusiveBrCalculator*
  inputHandler::getExclusiveBrCalculator(
                                        particlePointer const decayingParticle,
                                         particlePointer const productParticle,
                                          bool const productIsNotAntiparticle,
                                  std::list< int > const* const exclusionList )
  const
  {
    return exclusiveBrs->getBrCalculator( decayingParticle,
                                          productParticle,
                                          productIsNotAntiparticle,
                                          exclusionList );
  }

  inline std::list< int > const*
  inputHandler::getEmptyList()
  const
  {
    return &emptyList;
  }

  inline std::list< int > const*
  inputHandler::getJets()
  const
  {
    return &sixFlavorJets;
  }

  inline std::list< int > const*
  inputHandler::getJets5()
  const
  {
    return &fiveFlavorJets;
  }

  inline std::list< int > const*
  inputHandler::getJets4()
  const
  {
    return &fourFlavorJets;
  }

  inline std::list< int > const*
  inputHandler::getNotInJets5()
  const
  {
    return &notInFiveFlavorJets;
  }

  inline std::list< int > const*
  inputHandler::getNotInJets4()
  const
  {
    return &notInFourFlavorJets;
  }

  inline particleVectorPointer
  inputHandler::getColoredSparticles()
  const
  {
    return &coloredSparticles;
  }

  inline particleVectorPointer
  inputHandler::getSquarks()
  const
  {
    return &sixFlavorSquarks;
  }

  inline particleVectorPointer
  inputHandler::getSquarks5()
  const
  {
    return &fiveFlavorSquarks;
  }

  inline particleVectorPointer
  inputHandler::getSquarks4()
  const
  {
    return &fourFlavorSquarks;
  }

  inline particleVectorPointer
  inputHandler::getSdownTypes()
  const
  {
    return &sdownTypes;
  }

  inline particleVectorPointer
  inputHandler::getSupTypes()
  const
  {
    return &supTypes;
  }

  inline std::vector< signedParticleShortcutPair* > const*
  inputHandler::getScoloredProductionCombinations()
  const
  {
    return &sparticleProductionCombinations;
  }

  inline particleVectorPointer
  inputHandler::getElectroweakinos()
  const
  {
    return &allElectroweakinos;
  }

  inline particleVectorPointer
  inputHandler::getUnstableElectroweakinos()
  const
  {
    return &unstableElectroweakinos;
  }

  inline particleVectorPointer
  inputHandler::getNeutralinos()
  const
  {
    return &allNeutralinos;
  }

  inline particleVectorPointer
  inputHandler::getUnstableNeutralinos()
  const
  {
    return &unstableNeutralinos;
  }

  inline particleVectorPointer
  inputHandler::getCharginos()
  const
  {
    return &allCharginos;
  }

  inline particleVectorPointer
  inputHandler::getChargedSleptons()
  const
  {
    return &allChargedSleptons;
  }

  inline particleVectorPointer
  inputHandler::getSneutrinos()
  const
  {
    return &allSneutrinos;
  }

  inline particleVectorPointer
  inputHandler::getNeutralEwsbScalarsAndPseudoscalars()
  const
  {
    return &neutralEwsbScalarsAndPseudoscalars;
  }

  inline particleVectorPointer
  inputHandler::getEwsbBosonsAndMassiveVectorBosons()
  const
  {
    return &ewsbBosonsAndMassiveVectorBosons;
  }

  inline particleVectorPointer
  inputHandler::getNeutralEwsbBosonsAndMassiveVectorBosons()
  const
  {
    return &neutralEwsbBosonsAndMassiveVectorBosons;
  }

  inline particleVectorPointer
  inputHandler::getChargedEwsbBosonsAndMassiveVectorBosons()
  const
  {
    return &chargedEwsbBosonsAndMassiveVectorBosons;
  }

  inline hardMuonFromTau const*
  inputHandler::getHardMuonFromTau()
  const
  {
    return &hardMuonFromTauFunction;
  }

  inline softMuonFromTau const*
  inputHandler::getSoftMuonFromTau()
  const
  {
    return &softMuonFromTauFunction;
  }

  inline hardPionFromTau const*
  inputHandler::getHardPionFromTau()
  const
  {
    return &hardPionFromTauFunction;
  }

  inline softPionFromTau const*
  inputHandler::getSoftPionFromTau()
  const
  {
    return &softPionFromTauFunction;
  }

  inline double
  inputHandler::getAverageSquarks4Mass()
  const
  {
    return updateDependentAverageSquarks4Mass->getAverageMass();
  }

  inline void
  inputHandler::addSquarkAsColoredSparticleAndInProduction(
                                            particlePointer const inputSquark )
  /* this adds the given squark to the list of squarks & the list of colored
   * sparticles, & also adds its relevant production combinations (with a
   * gluino, with its antisquark, with itself) to the list.
   */
  {
    sixFlavorSquarks.push_back( inputSquark );
    coloredSparticles.push_back( inputSquark );
    sparticleProductionCombinations.push_back(
                                   new signedParticleShortcutPair( inputSquark,
                                                                   true,
                                                                 gluinoPointer,
                                                                   true ) );
    sparticleProductionCombinations.push_back(
                                   new signedParticleShortcutPair( inputSquark,
                                                                   true,
                                                                   inputSquark,
                                                                   false ) );
    sparticleProductionCombinations.push_back(
                                   new signedParticleShortcutPair( inputSquark,
                                                                   true,
                                                                   inputSquark,
                                                                   true ) );
    onShellSquarkEffectiveMasses.push_back(
                                  new onShellSquarkMassHolder( inputSquark ) );
  }

  inline effectiveSquarkMassHolder*
  inputHandler::getOnShellEffectiveSquarkMass( particlePointer onShellSquark )
  const
  {
    effectiveSquarkMassHolder* returnPointer( NULL );
    for( std::vector< onShellSquarkMassHolder* >::const_iterator
         squarkIterator( onShellSquarkEffectiveMasses.begin() );
         onShellSquarkEffectiveMasses.end() > squarkIterator;
         ++squarkIterator )
    {
      if( (*squarkIterator)->getOnShellSquark() == onShellSquark )
      {
        returnPointer = *squarkIterator;
        squarkIterator = onShellSquarkEffectiveMasses.end();
      }
    }
    return returnPointer;
  }

  inline effectiveSquarkMassHolder*
  inputHandler::getSquarkMinusBosonEffectiveMass(
                                           particlePointer const onShellSquark,
                                            particlePointer const onShellBoson,
                                           particlePointer const onShellEwino )
  const
  {
    effectiveSquarkMassHolder* returnPointer( NULL );
    for( std::vector< squarkMinusBosonMassHolder* >::const_iterator
         massIterator( squarkMinusBosonEffectiveMasses.begin() );
         squarkMinusBosonEffectiveMasses.end() > massIterator;
         ++massIterator )
    {
      if( (*massIterator)->isEquivalent( onShellSquark,
                                         onShellBoson,
                                         onShellEwino ) )
      {
        returnPointer = *massIterator;
        massIterator = squarkMinusBosonEffectiveMasses.end();
      }
    }
    return returnPointer;
  }

  inline effectiveSquarkMassHolder*
  inputHandler::getSquarkPlusBosonEffectiveMass(
                                           particlePointer const onShellSquark,
                                            particlePointer const onShellBoson,
                                           particlePointer const onShellEwino )
  const
  {
    effectiveSquarkMassHolder* returnPointer( NULL );
    for( std::vector< squarkPlusBosonMassHolder* >::const_iterator
         massIterator( squarkPlusBosonEffectiveMasses.begin() );
         squarkPlusBosonEffectiveMasses.end() > massIterator;
         ++massIterator )
    {
      if( (*massIterator)->isEquivalent( onShellSquark,
                                         onShellBoson,
                                         onShellEwino ) )
      {
        returnPointer = *massIterator;
        massIterator = squarkMinusBosonEffectiveMasses.end();
      }
    }
    return returnPointer;
  }

}

#endif /* LHC_FASER_INPUT_HANDLING_STUFF_HPP_ */
