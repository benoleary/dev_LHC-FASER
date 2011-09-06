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
#include "LHC-FASER_derived_lepton_distributions.hpp"
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
                                readierForNewPoint* const readier )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~updateDependentGaugeInputs()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    getElectromagneticFineStructureConstant()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getHyperchargeGaugeCoupling()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getWeakGaugeCoupling()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getColorGaugeCoupling()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getWeakSine()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getWeakCosine()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

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
    calculateGaugeCouplings()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };


  /* this class is to hold frequently-accessed absolute masses which are
   * re-calculated for every new point. it is kept separate from the
   * inputHandler class for ease of (abusing) const correctness.
   */
  class updateDependentAbsoluteMasses : public getsReadiedForNewPoint
  {
  public:
    updateDependentAbsoluteMasses( particlePointer const gluinoPointer,
                                   readierForNewPoint* const readier )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~updateDependentAbsoluteMasses()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    getGluinoMass()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

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
                                readierForNewPoint* const readier )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~updateDependentAverageMass()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    getAverageMass()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

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
                  readierForNewPoint* const readier )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    ~inputHandler()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    bool
    isVerbose()
    const
    /* if this is true, lots of warnings will be printed. if false, the code
     * should stay quiet without warning of the impending segmentation
     * faults...
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setVerbosity( bool verbosity )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    bool
    usingNlo()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    void
    setUsingNloFlag( bool usingNloFlag )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    std::string const*
    getPathToGrids()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    std::string*
    getPathToCrossSectionGrids()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    std::string const*
    inspectPathToCrossSectionGrids()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    std::string*
    getPathToKinematicsGrids()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    std::string const*
    inspectPathToKinematicsGrids()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    CppSLHA::CppSLHA0 const*
    getCppSlha()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    CppSLHA::EW_scale_spectrum const*
    getSpectrum()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getElectromagneticFineStructureConstant()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getHyperchargeGaugeCoupling()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getWeakGaugeCoupling()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getColorGaugeCoupling()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getWeakSine()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getWeakCosine()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    readierForNewPoint*
    getReadier()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getDown()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getUp()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getStrange()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getCharm()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getBottom()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getTop()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getElectron()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getElectronNeutrino()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getMuon()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getMuonNeutrino()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getTauLepton()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getTauNeutrino()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getGluon()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getPhoton()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getZ()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getWPlus()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getLightNeutralEwsbScalar()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getHeavyNeutralEwsbScalar()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getNeutralEwsbPseudoscalar()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getChargedEwsbScalar()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getGluino()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getGluinoMass()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getSdownL()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getSupL()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getSstrangeL()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getScharmL()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getSbottomOne()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getStopOne()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getSdownR()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getSupR()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getSstrangeR()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getScharmR()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getSbottomTwo()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getStopTwo()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getSelectronL()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getElectronSneutrinoL()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getSmuonL()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getMuonSneutrinoL()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getStauOne()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getTauSneutrinoL()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getSelectronR()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getElectronSneutrinoR()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getSmuonR()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getMuonSneutrinoR()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getStauTwo()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getTauSneutrinoR()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getNeutralinoOne()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getNeutralinoTwo()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getNeutralinoThree()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getNeutralinoFour()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getCharginoOne()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particlePointer
    getCharginoTwo()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    bool
    isIn( int const pdgCode,
          particleVectorPointer const particleSet )
    const
    /* this checks to see if the particle with the absolute value of the given
     * PDG particle code is in the given collection. (hence it checks to see if
     * the given particle or its antiparticle are in the collection.)
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    /*directDecayChecker*
    getDirectDecayChecker( particlePointer const decayerParticle,
                           particleVectorPointer const productSet )
    const
    *//* code after the classes in this .hpp file, or in the .cpp file. *///;
    exclusiveBrCalculator*
    getExclusiveBrCalculator( particlePointer const decayerParticle,
                              particlePointer const productParticle,
                              bool const productIsNotAntiparticle,
                              std::list< int > const* const exclusionList )
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    std::list< int > const*
    getEmptyList()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    std::list< int > const*
    getJets()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    std::list< int > const*
    getJets5()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    std::list< int > const*
    getJets4()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    std::list< int > const*
    getNotInJets5()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    std::list< int > const*
    getNotInJets4()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particleVectorPointer
    getColoredSparticles()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particleVectorPointer
    getSquarks()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particleVectorPointer
    getSquarks5()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particleVectorPointer
    getSquarks4()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particleVectorPointer
    getSdownTypes()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particleVectorPointer
    getSupTypes()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    std::vector< signedParticleShortcutPair* > const*
    getScoloredProductionCombinations()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particleVectorPointer
    getElectroweakinos()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particleVectorPointer
    getUnstableElectroweakinos()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particleVectorPointer
    getNeutralinos()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particleVectorPointer
    getUnstableNeutralinos()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particleVectorPointer
    getCharginos()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particleVectorPointer
    getChargedSleptons()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particleVectorPointer
    getSneutrinos()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particleVectorPointer
    getNeutralEwsbScalarsAndPseudoscalars()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    particleVectorPointer
    getEwsbBosonsAndMassiveVectorBosons()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    hardMuonFromTau const*
    getHardMuonFromTau()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    softMuonFromTau const*
    getSoftMuonFromTau()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    hardPionFromTau const*
    getHardPionFromTau()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    softPionFromTau const*
    getSoftPionFromTau()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getAverageSquarks4Mass()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    quarkOrLeptonLeftHandedness( int const sfermionCode,
                                 int const ewinoCode )
    const
    /* this returns the square of the coupling associated with the left-handed
     * projection operator of the Feynman rule for the vertex of the given
     * sfermion with the given electroweakino & the appropriate quark, divided
     * by the sum of the squares of the couplings for the left- & right-handed
     * projection operators, based on Rosiek's Feynman rules. this has to be
     * changed if general squark flavor mixing is incorporated. also, it
     * ignores the possibility of right-handed neutrinos.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    effectiveSquarkMassHolder*
    getOnShellEffectiveSquarkMass( particlePointer onshellSquark )
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    bool verbosity;
    bool usingNloFlag;
    readierForNewPoint* const readier;
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
     * notInJets5 is just top & antitop, & notInJets4 is notInJets5 plus
     * bottom & antibottom.
     * jets5 is jets minus everything in not_int_jets5, likewise jets4 is jets
     * minus everything in not_int_jets4.
     * other exotic heavy quarks should either be labelled as not counting as
     * a jet in the spectrum, or put in notInJets5 in the same manner as top.
     * hence, without additional light jet particles, jets5 is
     * { gluon + the 5 light quarks (down, up, strange, charm, bottom)
     * + the 5 light antiquarks }, & jets4 is jets5 with bottom & antibottom
     * removed.
     */
    std::list< int > jets;
    std::list< int > jets5;
    std::list< int > jets4;
    std::list< int > notInJets5;
    std::list< int > notInJets4;

    std::vector< particlePointer > coloredSparticles;
    // squarks4 refers to { sdown, sup, sstrange, scharm }, while squarks5
    // is squarks4 plus sbottom, & squarks is squarks5 plus stop.
    std::vector< particlePointer > squarks;
    std::vector< particlePointer > squarks5;
    std::vector< particlePointer > squarks4;
    std::vector< particlePointer > sdownTypes;
    std::vector< particlePointer > supTypes;
    std::vector< signedParticleShortcutPair* >
    sparticleProductionCombinations;
    std::vector< onshellSquarkMassHolder* > onshellSquarkEffectiveMasses;

    std::vector< particlePointer > electroweakinos;
    std::vector< particlePointer > unstableElectroweakinos;
    std::vector< particlePointer > neutralinos;
    std::vector< particlePointer > unstableNeutralinos;
    std::vector< particlePointer > charginos;
    std::vector< particlePointer > chargedSleptons;
    std::vector< particlePointer > sneutrinos;

    std::vector< particlePointer > neutralEwsbScalarsAndPseudoscalars;
    std::vector< particlePointer > ewsbBosonsAndMassiveVectorBosons;

    //decayCheckerHandler* decayCheckers;
    exclusiveBrHandler* exclusiveBrs;

    hardMuonFromTau hardMuonFromTauFunction;
    softMuonFromTau softMuonFromTauFunction;
    hardPionFromTau hardPionFromTauFunction;
    softPionFromTau softPionFromTauFunction;

    void
    addSquarkAsColoredSparticleAndInProduction(
                                            particlePointer const inputSquark )
    /* this adds the given squark to the list of squarks & the list of colored
     * sparticles, & also adds its relevant production combinations (with
     * a gluino, with its antisquark, with itself) to the list.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };  // end of inputHandler class.




  // inline functions:

  inline double
  updateDependentGaugeInputs::getElectromagneticFineStructureConstant()
  {
    if( needsToPrepareForThisPoint() )
    {
      calculateGaugeCouplings();
      finishPreparingForThisPoint();
    }
    return electromagneticFineStructureConstant;
  }

  inline double
  updateDependentGaugeInputs::getHyperchargeGaugeCoupling()
  {
    if( needsToPrepareForThisPoint() )
    {
      calculateGaugeCouplings();
      finishPreparingForThisPoint();
    }
    return hyperchargeGaugeCoupling;
  }

  inline double
  updateDependentGaugeInputs::getWeakGaugeCoupling()
  {
    if( needsToPrepareForThisPoint() )
    {
      calculateGaugeCouplings();
      finishPreparingForThisPoint();
    }
    return weakGaugeCoupling;
  }

  inline double
  updateDependentGaugeInputs::getColorGaugeCoupling()
  {
    if( needsToPrepareForThisPoint() )
    {
      calculateGaugeCouplings();
      finishPreparingForThisPoint();
    }
    return colorGaugeCoupling;
  }

  inline double
  updateDependentGaugeInputs::getWeakSine()
  {
    if( needsToPrepareForThisPoint() )
    {
      calculateGaugeCouplings();
      finishPreparingForThisPoint();
    }
    return weakSine;
  }

  inline double
  updateDependentGaugeInputs::getWeakCosine()
  {
    if( needsToPrepareForThisPoint() )
    {
      calculateGaugeCouplings();
      finishPreparingForThisPoint();
    }
    return weakCosine;
  }



  inline double
  updateDependentAbsoluteMasses::getGluinoMass()
  {
    if( needsToPrepareForThisPoint() )
      {
        gluinoMass = gluinoPointer->get_absolute_mass();
        finishPreparingForThisPoint();
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
      finishPreparingForThisPoint();
    }
    return averageMass;
  }



  inline bool
  inputHandler::isVerbose()
  const
  // if this is true, lots of warnings will be printed. if false, the code
  // should stay quiet without warning of the impending segmentation faults...
  {
    return verbosity;
  }

  inline void
  inputHandler::setVerbosity( bool verbosity )
  {
    this->verbosity = verbosity;
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

  inline std::string*
  inputHandler::getPathToKinematicsGrids()
  {
    return &pathToKinematicsGrids;
  }

  inline std::string const*
  inputHandler::inspectPathToCrossSectionGrids()
  const
  {
    return &pathToCrossSectionGrids;
  }

  inline std::string const*
  inputHandler::inspectPathToKinematicsGrids()
  const
  {
    return &pathToKinematicsGrids;
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
    return readier;
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
    bool returnBool = false;
    for( std::vector< particlePointer >::const_iterator
         setIterator = particleSet->begin();
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

  /*inline directDecayChecker*
  inputHandler:: getDirectDecayChecker( particlePointer const decayingParticle,
                                       particleVectorPointer const productSet )
  const
  {
    return decayCheckers->getDecayChecker( decayingParticle,
                                           productSet );
  }*/

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
    return &jets;
  }

  inline std::list< int > const*
  inputHandler::getJets5()
  const
  {
    return &jets5;
  }

  inline std::list< int > const*
  inputHandler::getJets4()
  const
  {
    return &jets4;
  }

  inline std::list< int > const*
  inputHandler::getNotInJets5()
  const
  {
    return &notInJets5;
  }

  inline std::list< int > const*
  inputHandler::getNotInJets4()
  const
  {
    return &notInJets4;
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
    return &squarks;
  }

  inline particleVectorPointer
  inputHandler::getSquarks5()
  const
  {
    return &squarks5;
  }

  inline particleVectorPointer
  inputHandler::getSquarks4()
  const
  {
    return &squarks4;
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
    return &electroweakinos;
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
    return &neutralinos;
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
    return &charginos;
  }

  inline particleVectorPointer
  inputHandler::getChargedSleptons()
  const
  {
    return &chargedSleptons;
  }

  inline particleVectorPointer
  inputHandler::getSneutrinos()
  const
  {
    return &sneutrinos;
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
    squarks.push_back( inputSquark );
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
    onshellSquarkEffectiveMasses.push_back(
                                  new onshellSquarkMassHolder( inputSquark ) );
  }

  inline effectiveSquarkMassHolder*
  inputHandler::getOnShellEffectiveSquarkMass( particlePointer onshellSquark )
  const
  {
    effectiveSquarkMassHolder* returnPointer( NULL );
    for( std::vector< onshellSquarkMassHolder* >::const_iterator
         squarkIterator( onshellSquarkEffectiveMasses.begin() );
         onshellSquarkEffectiveMasses.end() > squarkIterator;
         ++squarkIterator )
    {
      if( (*squarkIterator)->getOnshellSquark() == onshellSquark )
      {
        returnPointer = *squarkIterator;
      }
    }
    return returnPointer;
  }

}

#endif /* LHC_FASER_INPUT_HANDLING_STUFF_HPP_ */
