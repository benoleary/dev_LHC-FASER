/*
 * LHC-FASER_full_cascade_stuff.hpp
 *
 *  Created on: 30 Mar 2011
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

#ifndef LHC_FASER_FULL_CASCADE_STUFF_HPP_
#define LHC_FASER_FULL_CASCADE_STUFF_HPP_

#include "LHC-FASER_input_handling_stuff.hpp"
#include "LHC-FASER_sparticle_decay_stuff.hpp"
#include "LHC-FASER_lepton_kinematics_stuff.hpp"
#include "LHC-FASER_electroweak_cascade_collection_stuff.hpp"

namespace LHC_FASER
{
  /* this class is for keeping track of a particular cascade decay of a colored
   * sparticle down to the last electroweak decay, such that it keeps track of
   * all the sparticles of a cascade decay. it is an abstract base class.
   *
   * try to remember that there's 1 of these per beam energy! thus any given
   * scolored could have e.g. a set of fullCascades for 7 TeV & another set
   * for 14 TeV (with the same cascades because the hierarchy is the same, BUT
   * different ewinoCascades & bosonCascades because each looks up its own
   * effective squark mass, which depends on the beam energy (though really
   * should also depend on which events are selected by jet+MET cuts, but it'd
   * be too impractical to implement that)).
   *
   * within the approximations made by this code, instances of this class cover
   * cascade decays to *both* charge of chargino, if the electroweakino is
   * such. getBrToEwino returns the sum of the branching ratios to both of the
   * charges (in the case of an initial squark, which does not decay through
   * a gluino, this will be just the branching ratio for whichever charge is
   * allowed, while if there a gluino in the cascade, for the part of the
   * cascade coming from the decay of the gluino (so, not including for example
   * the decay of a squark to another squark via a vector or spin-0 boson,
   * should this happen *earlier* in the cascade than the decay to & from the
   * gluino), it returns the sum of BRs which is twice the BR to a single
   * charge, because of the approximations used. getAcceptance should take care
   * to avoid double-counting for derived classes with cascades which start
   * with gluinos, insofar as if no asymmetric charge is requested (only
   * applicable to the leptons requested), it returns the acceptance for a
   * single charge (since the charge-conjugate has the same acceptance), but if
   * a combination with an asymmetric charge is requested, it returns half the
   * acceptance, to account for only half the BR returned actually
   * corresponding to cascades that can give the requested signal.
   *
   * because fullCascadeSet stores its fullCascades in a
   * minimalAllocationVector, the constructor doesn't set any of the instance's
   * pointers or values: this is left to setProperties. well, since the base
   * class is abstract, its constructor takes an argument which is unique to
   * each derived class, so it can be set safely since it should never to reset
   * to a different type. (typeOfColorfulCascade is used because I don't want
   * to get into the messy habit of dynamic casting.)
   */
  class fullCascade
  {
  public:
    typedef std::pair< particlePointer,
                       int > particleWithInt;
    enum colorfulCascadeType
    {
      // the longer the cascade, the less accurate it is likely to be...
      notSet /* marked as not yet a proper cascade */,
      squarkToEwino /* direct squark to electroweakino */,
      gluinoToEwino /* direct gluino to electroweakino */,
      squarkToSquarkThenMore
      /* squark to boson plus other squark then shorter cascade */
      /* treated by jets+MET grid as squarkToEwino */,
      squarkToGauginoThenMore /* squark to gaugino then shorter cascade */
      /* treated by jets+MET grid as squarkToEwino or squarkToGluinoToEwino */,
      gluinoOrElectroweakinoToSquarkThenMore
      /* gluino or electroweakino to squark then shorter cascade */
      /* treated by jets+MET grid as gluinoToSquarkToEwino */,
      sizeOfEnumPlusOne /* possibly useful */
    };

    fullCascade( colorfulCascadeType const typeOfColorfulCascade,
                 int const firstDecayBodyNumber,
                 int const maximumSmFermionsFromElectroweakCascades = 0,
                 inputHandler const* const inputShortcut = NULL,
                 particlePointer const initialScolored = NULL,
       double const beamEnergy = CppSLHA::CppSLHA_global::really_wrong_value );
    virtual
    ~fullCascade();

    particlePointer
    getInitialSparticle()
    const;
    virtual particlePointer
    getElectroweakinoAtEndOfScoloredness()
    const;
    void
    resetCachedBranchingRatio();
    /* this sets branchingRatioNeedsToBeRecalculated to true, so that
     * cachedBranchingRatio will be re-calculated next time the branching ratio
     * is requested. I'd prefer to do this with inheritance from
     * getsReadiedForNewPoint, but I'd have to write a variant of
     * minimalAllocationVector to allow for arguments to be given to
     * constructors so that the appropriate readierForNewPoint can be given
     * to the fullCascade constructor.
     */
    virtual double
    getBrToEwino( std::list< int > const* excludedSmParticles );
    /* this works out the branching ratio for the decays of colored sparticles
     * down to the 1st electroweakino in the cascade (including any decays to
     * EWSB bosons + squarks), using the given list of excluded SM particles
     * to exclude unwanted parts of the BR. by default, it multiplies the BR
     * for the initial decay by subcascadePointer->getBrToEwino.
     */
    double
    getTotalBrToEwino();
    /* this works out the branching ratio for the decays of colored
     * sparticles down to the 1st electroweakino in the cascade (including
     * any decays to bosons + squarks), without excluding any particles.
     */
    std::vector< particleWithInt* > const*
    getCascadeDefiner()
    const;
    colorfulCascadeType
    getColorfulCascadeType()
    const;
    virtual double
    getAcceptance( bool const initialSparticleIsNotAntiparticle,
                   acceptanceCutSet* const acceptanceCuts,
                   int const numberOfAdditionalJets,
                   int numberOfNegativeElectrons,
                   int numberOfPositiveElectrons,
                   int numberOfNegativeMuons,
                   int numberOfPositiveMuons )
    = 0;
    /* this calls the appropriate functions on ewinoCascade &, if not NULL,
     * bosonCascades, to build the required acceptance, taking into account
     * whether the charges should be swapped if scoloredIsNotAntiparticle is
     * false.
     */
    double
    getAcceptanceWithMinimumJets( bool const initialSparticleIsNotAntiparticle,
                                  acceptanceCutSet* const acceptanceCuts,
                                  int const minimumNumberOfAdditionalJets,
                                  int const numberOfNegativeElectrons,
                                  int const numberOfPositiveElectrons,
                                  int const numberOfNegativeMuons,
                                  int const numberOfPositiveMuons );
    /* this sums up getAcceptance for
     * numberOfAdditionalJets = minimumNumberOfAdditionalJets
     * to ( numberOfSmFermionsFromElectroweakDecaysPerFullCascade
     *      - numberOfNegativeElectrons
     *      - numberOfPositiveElectrons
     *      - numberOfNegativeMuons
     *      - numberOfPositiveMuons )
     */
    double
    leptonFlavorSummedWithMinimumJets(
                                  bool const initialSparticleIsNotAntiparticle,
                                        acceptanceCutSet* const acceptanceCuts,
                                       int const minimumNumberOfAdditionalJets,
                                             int const numberOfNegativeLeptons,
                                           int const numberOfPositiveLeptons );
    /* this sums up getAcceptanceWithMinimumJets for numberOfNegativeLeptons &
     * numberOfPositiveLeptons split between electrons & muons in every
     * combination.
     */
    virtual double
    specifiedJetsOneOssfMinusOsdfPair( acceptanceCutSet* const acceptanceCuts,
                                       int const numberOfAdditionalJets )
    = 0;
    // this should add up all combinations of jets while looking only for
    // one OSSF-OSDF pair.
    virtual double
    unspecifiedJetsSpecifiedOssfMinusOsdfPairs(
                                        acceptanceCutSet* const acceptanceCuts,
                                               int const numberOfLeptonPairs )
    = 0;
    // this should add up all appropriate OSSF-OSDF pair combinations.
    // the following just call various combinations of getAcceptance(...):
    double
    specifiedJetsSpecifiedChargeSummedLeptons(
                                        acceptanceCutSet* const acceptanceCuts,
                                              int const numberOfAdditionalJets,
                                               int const numberOfLeptons );
    double
    unspecifiedJetsSpecifiedChargeSummedLeptons(
                                        acceptanceCutSet* const acceptanceCuts,
                                                 int const numberOfLeptons );
    std::string
    getAsString()
    const;
    // this is for error messages and testing output.
    int
    getMaximumSmFermionsFromElectroweakCascades()
    const;


  protected:
    static int const maximumJetsFromEwCascade;
    static int const maximumLeptonsOfEachTypeFromEwCascade;

    inputHandler const* inputShortcut;
    particlePointer initialSparticle;
    double beamEnergy;
    colorfulCascadeType const typeOfColorfulCascade;
    int maximumSmFermionsFromElectroweakCascades;
    int const firstDecayBodyNumber;
    // this is the number of particles which are produced by the decay of
    // initialSparticle.
    fullCascade* subcascadePointer;
    minimalAllocationVector< particleWithInt > cascadeDefiner;
    particleWithInt* cascadeSegment;
    // this is just used for setting up cascadeDefiner properly.
    std::list< int > soughtDecayProductList;
    /* by default, there is only 1 particle to be sought in initialSparticle's
     * decays. the constructor sets up soughtDecayProductList to have just a
     * single entry, & setProperties just sets that entry, & nothing else makes
     * any changes to soughtDecayProductList. derived classes involving decays
     * to bosons add entries to their own instances.
     */
    double cachedBranchingRatio;
    bool branchingRatioNeedsToBeRecalculated;

    void
    buildOn( fullCascade* const copySource );
    /* this copies the basic stuff (inputShortcut, beamEnergy, cascadeDefiner)
     * from copySource, though not typeOfColorfulCascade (which is always set
     * by constructors) or initialSparticle (which depends on the type of
     * derived class), & also it sets up soughtDecayProductList to look for
     * copySource->initialSparticle, & adds copySource->initialSparticle at the
     * end of cascadeDefiner, with firstDecayBodyNumber.
     */
  };


  namespace fullCascadeType
  {
    // this is derived class for direct squark to electroweakino fullCascades.
    // it should be further specified to be sdown-type or sup-type.
    class squarkDirectlyToElectroweak : public fullCascade
    {
    public:
      squarkDirectlyToElectroweak( inputHandler const* const inputShortcut,
                                   particlePointer const initialSquark,
                                   double const beamEnergy,
                              electroweakCascadeSet* const directEwinoCascades,
                          int const maximumSmFermionsFromElectroweakCascades );
      virtual
      ~squarkDirectlyToElectroweak();

      virtual bool
      isOpen()
      = 0;
      /* this returns true if the squark is heavy enough to decay into the
       * electroweakino, false otherwise. it also sorts out whether it should
       * be using the decays involving a W boson.
       */
      virtual particlePointer
      getElectroweakinoAtEndOfScoloredness()
      const;
      virtual double
      getBrToEwino( std::list< int > const* excludedSmParticles );
      /* this works out the branching ratio for the decays of colored
       * sparticles down to the 1st electroweakino in the cascade (including
       * any decays to bosons + squarks), using the given list of excluded SM
       * particles to exclude unwanted parts of the BR.
       */
      virtual double
      unspecifiedJetsSpecifiedOssfMinusOsdfPairs(
                                        acceptanceCutSet* const acceptanceCuts,
                                               int const numberOfLeptonPairs );
      // this should add up all appropriate OSSF-OSDF pair combinations.


    protected:
      electroweakCascadeSet* directEwinoCascades;
    };


    namespace squarkDirectlyToElectroweakType
    {
      // this is derived class for direct sdown-type to electroweakino
      // fullCascades.
      class sdownType : public squarkDirectlyToElectroweak
      {
      public:
        sdownType( inputHandler const* const inputShortcut,
                   particlePointer const initialSquark,
                   double const beamEnergy,
                   electroweakCascadeSet* const directEwinoCascades );
        virtual
        ~sdownType();

        virtual bool
        isOpen();
        /* this returns true if the squark is heavy enough to decay into the
         * electroweakino, false otherwise. it also sorts out whether it should
         * be using the decays involving a W boson.
         */
        virtual double
        getAcceptance( bool const initialSparticleIsNotAntiparticle,
                       acceptanceCutSet* const acceptanceCuts,
                       int const numberOfAdditionalJets,
                       int numberOfNegativeElectrons,
                       int numberOfPositiveElectrons,
                       int numberOfNegativeMuons,
                       int numberOfPositiveMuons );
        /* this calls the appropriate functions on ewinoCascades to build the
         * required acceptance, taking into account whether the charges should
         * be swapped if scoloredIsNotAntiparticle is false.
         */
        virtual double
        specifiedJetsOneOssfMinusOsdfPair(
                                        acceptanceCutSet* const acceptanceCuts,
                                           int const numberOfAdditionalJets );
        // this should add up all combinations of jets while looking only for
        // one OSSF-OSDF pair.


      protected:
        bool ewinoFlipsCharge;
      };


      /* this is derived class for direct sup-type to electroweakino
       * fullCascades, which also covers if the squark decays to an additional
       * W.
       */
      class supType : public squarkDirectlyToElectroweak
      {
      public:
        supType( inputHandler const* const inputShortcut,
                 particlePointer const initialSquark,
                 double const beamEnergy,
                 electroweakCascadeSet* const directEwinoCascades,
                 electroweakCascadeSet* const ewinoWithWCascades,
                 electroweakCascadeSet* const bosonCascades );
        virtual
        ~supType();

        virtual bool
        isOpen();
        /* this returns true if the squark is heavy enough to decay into the
         * electroweakino, false otherwise. it also sorts out whether it should
         * be using the decays involving a W boson.
         */
        virtual double
        getAcceptance( bool const initialSparticleIsNotAntiparticle,
                       acceptanceCutSet* const acceptanceCuts,
                       int const numberOfAdditionalJets,
                       int numberOfNegativeElectrons,
                       int numberOfPositiveElectrons,
                       int numberOfNegativeMuons,
                       int numberOfPositiveMuons );
        /* this calls the appropriate functions on directEwinoCascades to build
         * the required acceptance or ewinoWithWCascades & bosonCascades if the
         * decay with a boson is not negligible, taking into account whether
         * the charges should be swapped if scoloredIsNotAntiparticle is false.
         */
        virtual double
        specifiedJetsOneOssfMinusOsdfPair(
                                        acceptanceCutSet* const acceptanceCuts,
                                           int const numberOfAdditionalJets );
        // this should add up all combinations of jets while looking only for
        // one OSSF-OSDF pair.


      protected:
        electroweakCascadeSet* ewinoWithWCascades;
        electroweakCascadeSet* bosonCascades;
        particlePointer wBoson;
        double ewinoMass;
        bool shouldUseDecaysWithW;
        double directFraction;
        double wFraction;
        std::list< int > decayProductListIncludingW;

        bool
        decayWithWIsNotNegligible();
        // this returns true if the decay involving a W boson is not
        // negligible, also setting up the relevant fractions, false otherwise.
        double
        getCombinedAcceptance( acceptanceCutSet* const acceptanceCuts,
                               int const numberOfAdditionalJets,
                               int const numberOfNegativeElectrons,
                               int const numberOfPositiveElectrons,
                               int const numberOfNegativeMuons,
                               int const numberOfPositiveMuons );
      };

    }  // end of squarkDirectlyToElectroweakType namespace


    // this is derived class for direct gluino to electroweakino
    // fullCascades..
    class gluinoDirectlyToElectroweak : public fullCascade
    {
    public:
      gluinoDirectlyToElectroweak( inputHandler const* const inputShortcut,
                                   double const beamEnergy,
                                  electroweakCascadeSet* const ewinoCascades );
      virtual
      ~gluinoDirectlyToElectroweak();

      bool
      isOpen();
      /* this returns true if the squark is heavy enough to decay into the
       * electroweakino, false otherwise. it also sorts out whether it should
       * be using the decays involving a W boson.
       */
      virtual double
      getBrToEwino( std::list< int > const* excludedSmParticles );
      /* this works out the branching ratio for the decays of colored
       * sparticles down to the 1st electroweakino in the cascade (including
       * any decays to bosons + squarks), using the given list of excluded SM
       * particles to exclude unwanted parts of the BR.
       */
      virtual particlePointer
      getElectroweakinoAtEndOfScoloredness()
      const;
      virtual double
      getAcceptance( bool const initialSparticleIsNotAntiparticle,
                     acceptanceCutSet* const acceptanceCuts,
                     int const numberOfAdditionalJets,
                     int numberOfNegativeElectrons,
                     int numberOfPositiveElectrons,
                     int numberOfNegativeMuons,
                     int numberOfPositiveMuons );
      /* this calls the appropriate functions on ewinoCascade, summing over the
       * 2 charge versions if the electroweakino is a chargino. it ignores the
       * value of scoloredIsNotAntiparticle.
       */
      virtual double
      specifiedJetsOneOssfMinusOsdfPair(
                                        acceptanceCutSet* const acceptanceCuts,
                                         int const numberOfAdditionalJets );
      // this should add up all combinations of jets while looking only for
      // one OSSF-OSDF pair.
      virtual double
      unspecifiedJetsSpecifiedOssfMinusOsdfPairs(
                                        acceptanceCutSet* const acceptanceCuts,
                                               int const numberOfLeptonPairs );
      // this should add up all appropriate OSSF-OSDF pair combinations.


    protected:
      electroweakCascadeSet* ewinoCascades;
    };


    /* this is derived class for compound fullCascades beginning with a squark,
     * which decays to a boson, a jet, & the initial sparticle of
     * subcascadePointer.
     */
    class squarkByBosonToCompound : public fullCascade
    {
    public:
      squarkByBosonToCompound();
      virtual
      ~squarkByBosonToCompound();

      virtual fullCascade*
      setProperties( particlePointer const initialSquark,
                     electroweakCascadeSet* const bosonCascades,
                     fullCascade* const subcascadePointer );
      virtual double
      getAcceptance( bool const initialSparticleIsNotAntiparticle,
                     acceptanceCutSet* const acceptanceCuts,
                     int const numberOfAdditionalJets,
                     int numberOfNegativeElectrons,
                     int numberOfPositiveElectrons,
                     int numberOfNegativeMuons,
                     int numberOfPositiveMuons );
      /* this calls the appropriate functions on subcascadePointer &
       * bosonCascades to build the required acceptance, taking into account
       * whether the charges should be swapped if scoloredIsNotAntiparticle is
       * false.
       */
      virtual double
      specifiedJetsOneOssfMinusOsdfPair(
                                        acceptanceCutSet* const acceptanceCuts,
                                         int const numberOfAdditionalJets );
      // this should add up all combinations of jets while looking only for
      // one OSSF-OSDF pair.
      virtual double
      unspecifiedJetsSpecifiedOssfMinusOsdfPairs(
                                        acceptanceCutSet* const acceptanceCuts,
                                               int const numberOfLeptonPairs );
      // this should add up all appropriate OSSF-OSDF pair combinations.


    protected:
      electroweakCascadeSet* bosonCascades;
      bool bosonFlipsCharge;
    };


    /* this is derived class for compound fullCascades beginning with a squark,
     * which decays to a jet & the initial sparticle of subcascadePointer. it
     * still needs specializing to either a sdown-type or sup-type.
     */
    class squarkByJetToCompound : public fullCascade
    {
    public:
      squarkByJetToCompound();
      virtual
      ~squarkByJetToCompound();

      virtual fullCascade*
      setProperties( particlePointer const initialSquark,
                     fullCascade* const subcascadePointer,
          electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource )
      = 0;
      virtual double
      getAcceptance( bool const initialSparticleIsNotAntiparticle,
                     acceptanceCutSet* const acceptanceCuts,
                     int const numberOfAdditionalJets,
                     int numberOfNegativeElectrons,
                     int numberOfPositiveElectrons,
                     int numberOfNegativeMuons,
                     int numberOfPositiveMuons )
      = 0;
      /* this calls the appropriate functions on subcascadePointer to build the
       * required acceptance, taking into account whether the charges should be
       * swapped if scoloredIsNotAntiparticle is false.
       */
      double
      unspecifiedJetsSpecifiedOssfMinusOsdfPairs(
                                      acceptanceCutSet* const acceptanceCuts,
                                               int const numberOfLeptonPairs );
      // this should add up all appropriate OSSF-OSDF pair combinations.


    //protected:
      // nothing.
    };


    namespace squarkByJetToCompoundType
    {
      /* this is derived class for compound fullCascades beginning with a
       * sdown-type, which decays to a jet & the initial sparticle of
       * subcascadePointer.
       */
      class sdownType : public squarkByJetToCompound
      {
      public:
        sdownType();
        virtual
        ~sdownType();

        virtual fullCascade*
        setProperties( particlePointer const initialSquark,
                       fullCascade* const subcascadePointer,
                       electroweakCascadesForOneBeamEnergy* const
                       electroweakCascadeSource = NULL );
        virtual double
        getAcceptance( bool const initialSparticleIsNotAntiparticle,
                       acceptanceCutSet* const acceptanceCuts,
                       int const numberOfAdditionalJets,
                       int numberOfNegativeElectrons,
                       int numberOfPositiveElectrons,
                       int numberOfNegativeMuons,
                       int numberOfPositiveMuons );
        /* this calls the appropriate functions on subcascadePointer to build
         * the required acceptance, taking into account whether the charges
         * should be swapped if scoloredIsNotAntiparticle is false.
         */
        double
        specifiedJetsOneOssfMinusOsdfPair(
                                        acceptanceCutSet* const acceptanceCuts,
                                           int const numberOfAdditionalJets );
        // this should add up all combinations of jets while looking only for
        // one OSSF-OSDF pair.


      //protected:
        // nothing.
      };


      /* this is derived class for compound fullCascades beginning with a
       * sup-type, which decays to a jet, maybe a boson, & the initial
       * sparticle of subcascadePointer.
       */
      class supType : public squarkByJetToCompound
      {
      public:
        supType();
        virtual
        ~supType();

        virtual fullCascade*
        setProperties( particlePointer const initialSquark,
                       fullCascade* const subcascadePointer,
         electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource );
        virtual double
        getAcceptance( bool const initialSparticleIsNotAntiparticle,
                       acceptanceCutSet* const acceptanceCuts,
                       int const numberOfAdditionalJets,
                       int numberOfNegativeElectrons,
                       int numberOfPositiveElectrons,
                       int numberOfNegativeMuons,
                       int numberOfPositiveMuons );
        /* this calls the appropriate functions on subcascadePointer to build
         * the required acceptance, & also with bosonCascades if the decay with
         * a boson is not negligible, taking into account whether the charges
         * should be swapped if scoloredIsNotAntiparticle is false.
         */
        double
        specifiedJetsOneOssfMinusOsdfPair(
                                        acceptanceCutSet* const acceptanceCuts,
                                           int const numberOfAdditionalJets );
        // this should add up all combinations of jets while looking only for
        // one OSSF-OSDF pair.


      protected:
        electroweakCascadesForOneBeamEnergy* electroweakCascadeSource;
        electroweakCascadeSet* bosonCascades;
        particlePointer wBoson;
        double ewinoMass;
        bool shouldUseDecaysWithW;
        double directFraction;
        double wFraction;
        std::list< int > decayProductListIncludingW;

        bool
        decayWithWIsNotNegligible();
        // this returns true if the decay involving a W boson is not
        // negligible, also setting up the relevant fractions, false otherwise.
        double
        getCombinedAcceptance( acceptanceCutSet* const acceptanceCuts,
                               int const numberOfAdditionalJets,
                               int const numberOfNegativeElectrons,
                               int const numberOfPositiveElectrons,
                               int const numberOfNegativeMuons,
                               int const numberOfPositiveMuons );
      };

    }  // end of squarkByJetToCompoundType namespace
    typedef squarkByJetToCompoundType::sdownType sdownByJetToCompound;
    typedef squarkByJetToCompoundType::supType supByJetToCompound;


    /* this is derived class for compound fullCascades beginning with a
     * gluino or electroweakino, which decays to a jet, maybe a boson, & the
     * initial sparticle of subcascadePointer.
     */
    class gluinoOrElectroweakinoToCompound : public fullCascade
    {
    public:
      gluinoOrElectroweakinoToCompound();
      virtual
      ~gluinoOrElectroweakinoToCompound();

      virtual double
      specifiedJetsOneOssfMinusOsdfPair(
                                        acceptanceCutSet* const acceptanceCuts,
                                         int const numberOfAdditionalJets );
      // this should add up all combinations of jets while looking only for
      // one OSSF-OSDF pair.
      virtual double
      unspecifiedJetsSpecifiedOssfMinusOsdfPairs(
                                        acceptanceCutSet* const acceptanceCuts,
                                               int const numberOfLeptonPairs );
      // this should add up all appropriate OSSF-OSDF pair combinations.


    protected:
      electroweakCascadesForOneBeamEnergy* electroweakCascadeSource;
      electroweakCascadeSet* bosonCascades;
      particlePointer wBoson;
      double ewinoMass;
      bool shouldUseDecaysWithW;
      double directFraction;
      double wFraction;
      std::list< int > decayProductListIncludingW;
    };


    namespace gluinoOrElectroweakinoToCompoundType
    {
      /* this is derived class for compound fullCascades beginning with a
       * gluino or neutralino, which decays to a jet, maybe a boson & the
       * initial sparticle of subcascadePointer.
       */
      class gluinoOrNeutralinoSet : public gluinoOrElectroweakinoToCompound
      {
      public:
        gluinoOrNeutralinoSet();
        virtual
        ~gluinoOrNeutralinoSet();

        virtual fullCascade*
        setProperties( particlePointer const initialSparticle,
                       fullCascade* const subcascadePointer,
         electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource );
        virtual double
        getAcceptance( bool const initialSparticleIsNotAntiparticle,
                       acceptanceCutSet* const acceptanceCuts,
                       int const numberOfAdditionalJets,
                       int numberOfNegativeElectrons,
                       int numberOfPositiveElectrons,
                       int numberOfNegativeMuons,
                       int numberOfPositiveMuons );
        /* this calls the appropriate functions on subcascadePointer to build
         * the required acceptance, & also with bosonCascades if the decay with
         * a boson is not negligible, taking into account whether the charges
         * should be swapped if scoloredIsNotAntiparticle is false.
         */


      protected:
        bool
        decayWithWIsNotNegligible();
        // this returns true if the decay involving a W boson is not
        // negligible, also setting up the relevant fractions, false otherwise.
        double
        getCombinedAcceptance( acceptanceCutSet* const acceptanceCuts,
                               int const numberOfAdditionalJets,
                               int const numberOfNegativeElectrons,
                               int const numberOfPositiveElectrons,
                               int const numberOfNegativeMuons,
                               int const numberOfPositiveMuons );
      };


      /* this is derived class for compound fullCascades beginning with a
       * chargino, which decays to a jet, maybe a boson, & the initial
       * sparticle of subcascadePointer.
       */
      class charginoSet : public gluinoOrElectroweakinoToCompound
      {
      public:
        charginoSet();
        virtual
        ~charginoSet();

        virtual fullCascade*
        setProperties( particlePointer const initialSparticle,
                       fullCascade* const subcascadePointer,
         electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource );
        virtual double
        getAcceptance( bool const initialSparticleIsNotAntiparticle,
                       acceptanceCutSet* const acceptanceCuts,
                       int const numberOfAdditionalJets,
                       int numberOfNegativeElectrons,
                       int numberOfPositiveElectrons,
                       int numberOfNegativeMuons,
                       int numberOfPositiveMuons );
        /* this calls the appropriate functions on subcascadePointer to build
         * the required acceptance, & also with bosonCascades if the decay with
         * a boson is not negligible, taking into account whether the charges
         * should be swapped if scoloredIsNotAntiparticle is false.
         */


      protected:
        double squarkWithWFraction;
        double antisquarkWithWFraction;
        bool decayingToSupType;

        bool
        decayWithWIsNotNegligible();
        // this returns true if the decay involving a W boson is not
        // negligible, also setting up the relevant fractions, false otherwise.
        double
        getCombinedAcceptance( acceptanceCutSet* const acceptanceCuts,
                               int const numberOfAdditionalJets,
                               int const numberOfNegativeElectrons,
                               int const numberOfPositiveElectrons,
                               int const numberOfNegativeMuons,
                               int const numberOfPositiveMuons );
        // if this is called, it is the decay of a positive chargino.
      };

    }  // end of gluinoOrElectroweakinoType namespace
    typedef gluinoOrElectroweakinoToCompoundType::gluinoOrNeutralinoSet
    gluinoOrNeutralinoToCompound;
    typedef gluinoOrElectroweakinoToCompoundType::charginoSet
    charginoToCompound;

  }
  typedef
  fullCascadeType::gluinoOrNeutralinoToCompound
  compoundGluinoOrNeutralinoFullCascade;
  typedef
  fullCascadeType::squarkDirectlyToElectroweak
  squarkDirectlyToElectroweakCascade;
  typedef
  fullCascadeType::gluinoDirectlyToElectroweak
  gluinoDirectlyToElectroweakCascade;


  /* this class holds the possible cascades for a given scolored, as well as a
   * pointer to a std::list of the other fullCascadeSet instances, which get
   * ordered by mass of initial colored sparticle by whichever fullCascadeSet
   * is called 1st.
   */
  class fullCascadeSet : public getsReadiedForNewPoint
  {
  public:
    static bool
    massOrdered( fullCascadeSet* firstSet,
                 fullCascadeSet* secondSet );

    fullCascadeSet( inputHandler const* const inputShortcut,
                    particlePointer const initialSparticle,
                electroweakCascadesForOneBeamEnergy* const electroweakCascades,
                    double const beamEnergy );
    virtual
    ~fullCascadeSet();

    std::vector< fullCascade* >*
    getOpenCascades();
    // this calls setUpCascades() if it needs to be readied for this point.
    particlePointer
    getInitialSparticle()
    const;


  protected:
    inputHandler const* const inputShortcut;
    particlePointer const initialSparticle;
    double initialSparticleMass;
    electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource;
    std::vector< fullCascade* > openCascades;
    // this holds pointers to all the open cascades.
    std::list< fullCascadeSet* >* orderedCascadeSets;
    std::list< fullCascadeSet* >::iterator setIterator;
    std::vector< fullCascade* >* potentialSubcascades;
    double subcascadeBranchingRatio;
    double const beamEnergy;
    std::list< int > singleSpecifiedDecayProductList;
    int singleSpecifiedDecayProduct;
    std::vector< particlePointer > const* potentialDecayProducts;
    fullCascadeSet* potentialSubcascadeSet;

    virtual void
    setUpCascades()
    = 0;
    // this should set up all the open cascades.
  };


  class fullCascadeSetOrderer : public getsReadiedForNewPoint
  {
  public:
    fullCascadeSetOrderer( inputHandler const* const inputShortcut );
    ~fullCascadeSetOrderer();

    std::list< fullCascadeSet* >*
    getNeutralinoColoredCascades();
    std::list< fullCascadeSet* >*
    getCharginoColoredCascades();
    std::list< fullCascadeSet* >*
    getSdownTypeCascades();
    std::list< fullCascadeSet* >*
    getSupTypeCascades();
    std::list< fullCascadeSet* >*
    getNeutralinoColoredCascadesWithoutOrdering();
    std::list< fullCascadeSet* >*
    getCharginoColoredCascadesWithoutOrdering();
    std::list< fullCascadeSet* >*
    getSdownTypeCascadesWithoutOrdering();
    std::list< fullCascadeSet* >*
    getSupTypeCascadesWithoutOrdering();


  protected:
    inputHandler const* const inputShortcut;
    std::list< fullCascadeSet* > sdownTypeColoredCascades;
    std::list< fullCascadeSet* > supTypeColoredCascades;
    std::list< fullCascadeSet* > allNeutralinoColoredCascades;
    std::list< fullCascadeSet* > allCharginoColoredCascades;
    std::list< fullCascadeSet* > relevantNeutralinoColoredCascades;
    std::list< fullCascadeSet* > relevantCharginoColoredCascades;
    double lightestSquarkMass;

    void
    orderCascades();
  };


  namespace fullCascadeSetType
  {
    /* this class specializes fullCascade for the case of cascades beginning
     * with a squark, but is still abstract, since it has to be further
     * specialized to sdown-type or sup-type.
     */
    class squarkSet : public fullCascadeSet
    {
    public:
      squarkSet( inputHandler const* const inputShortcut,
                 particlePointer const initialScolored,
           electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource,
                 fullCascadeSetOrderer* const setOrderer,
                 fullCascadeSet* const gluinoFullCascade,
                 double const beamEnergy,
                 bool const ewinoCodeIsAlwaysPositive = true,
                 bool const bosonCodeIsAlwaysPositive = true );
      virtual
      ~squarkSet();


    protected:
      fullCascadeSetOrderer* const setOrderer;
      fullCascadeSet* const gluinoFullCascade;
      std::vector< squarkDirectlyToElectroweakCascade* > directToEwinoCascades;
      minimalAllocationVector< fullCascadeType::squarkByBosonToCompound >
      compoundByBosonCascades;
      bool ewinoCodeIsAlwaysPositive;
      bool bosonCodeIsAlwaysPositive;

      virtual void
      setUpCascades();
      // this should set up all the open cascades.
      virtual void
      buildSquarkCompoundCascades()
      = 0;
      // this does the job of finding the right squark + EW boson combinations.
      void
      addEwinoCompoundCascade( bool ewinoCodeIsPositive );
      // this relies on orderedCascadeSets having already been set correctly.
      void
      addByJetCascade();
      // this relies on singleSpecifiedDecayProductList.front() &
      // potentialSubcascadeSet having already been set correctly.
      void
      addSquarkCompoundCascade( bool bosonCodeIsPositive );
      // this relies on orderedCascadeSets & potentialDecayProducts having
      // already been set correctly.
      virtual void
      clearCompoundByJetCascades()
      = 0;
      // this should clear the minimalAllocationVector in the derived class.
      virtual fullCascade*
      getNewCompoundByJetCascade( fullCascade* const subcascadePointer )
      = 0;
      // this should add the appropriate new cascade to the
      // minimalAllocationVector in the derived class.
    };


    namespace squarkSetType
    {
      // this class specializes fullCascade for the case of cascades beginning
      // with a down-type squark.
      class sdownType : public squarkSet
      {
      public:
        sdownType( inputHandler const* const inputShortcut,
                   particlePointer const initialScolored,
           electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource,
                   fullCascadeSetOrderer* const setOrderer,
                   fullCascadeSet* const gluinoFullCascade,
                   double const beamEnergy );
        virtual
        ~sdownType();


      protected:
        minimalAllocationVector< fullCascadeType::sdownByJetToCompound >
        compoundByJetCascades;

        virtual void
        findOpenDirectCascades();
        // this puts all open direct cascades into openCascades.
        virtual void
        buildSquarkCompoundCascades();
        // this does the job of finding the right squark + EW boson
        // combinations.
        virtual void
        clearCompoundByJetCascades();
        // this should clear the minimalAllocationVector in the derived class.
        virtual fullCascade*
        getNewCompoundByJetCascade( fullCascade* const subcascadePointer );
        // this should add the appropriate new cascade to the
        // minimalAllocationVector in the derived class.
      };


      // this class specializes fullCascade for the case of cascades beginning
      // with an up-type squark.
      class supType : public squarkSet
      {
      public:
        supType( inputHandler const* const inputShortcut,
                 particlePointer const initialScolored,
           electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource,
                 fullCascadeSetOrderer* const setOrderer,
                 fullCascadeSet* const gluinoFullCascade,
                 double const beamEnergy );
        virtual
        ~supType();


      protected:
        minimalAllocationVector< fullCascadeType::supByJetToCompound >
        compoundByJetCascades;
        std::list< int > twoSpecifiedDecayProductsList;
        particlePointer appropriateSdownForWDecay;

        virtual void
        findOpenDirectCascades();
        // this puts all open direct cascades into openCascades.
        virtual void
        buildSquarkCompoundCascades();
        // this does the job of finding the right squark + EW boson
        // combinations.
        virtual void
        clearCompoundByJetCascades();
        // this should clear the minimalAllocationVector in the derived class.
        virtual fullCascade*
        getNewCompoundByJetCascade( fullCascade* const subcascadePointer );
        // this should add the appropriate new cascade to the
        // minimalAllocationVector in the derived class.
      };

    }  // end of squarkSetType namespace


    // this class specializes fullCascade for the case of cascades beginning
    // with a gluino or an electroweakino.
    class gluinoOrElectroweakinoSet : public fullCascadeSet
    {
    public:
      gluinoOrElectroweakinoSet( inputHandler const* const inputShortcut,
                                 particlePointer const initialSparticle,
           electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource,
                                 fullCascadeSetOrderer* const setOrderer,
                                 double const beamEnergy );
      virtual
      ~gluinoOrElectroweakinoSet();


    protected:
      fullCascadeSetOrderer* const setOrderer;

      virtual void
      setUpCascades();
      // this should set up all the open cascades.
      virtual void
      clearCompoundCascades()
      = 0;
      // this should clear the minimalAllocationVector in the derived class.
      virtual fullCascade*
      getNewCompoundCascade( fullCascade* const subcascadePointer )
      = 0;
      // this should add the appropriate new cascade to the
      // minimalAllocationVector in the derived class.
      void
      buildSquarkCompoundCascades();
      // this does the job of finding the right squark subcascades.
    };


    namespace gluinoOrElectroweakinoSetType
    {
      // this class specializes fullCascade for the case of cascades beginning
      // with a gluino.
      class gluinoSet : public gluinoOrElectroweakinoSet
      {
      public:
        gluinoSet( inputHandler const* const inputShortcut,
           electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource,
                   fullCascadeSetOrderer* const setOrderer,
                   double const beamEnergy );
        virtual
        ~gluinoSet();


      protected:
        std::vector< gluinoDirectlyToElectroweakCascade* >
        directToEwinoCascades;
        minimalAllocationVector< compoundGluinoOrNeutralinoFullCascade >
        compoundCascades;

        virtual void
        setUpCascades();
        // this should set up all the open cascades.
        virtual void
        clearCompoundCascades();
        // this should clear the minimalAllocationVector in the derived class.
        virtual fullCascade*
        getNewCompoundCascade( fullCascade* const subcascadePointer );
        // this should add the appropriate new cascade to the
        // minimalAllocationVector in the derived class.
      };


      // this class specializes fullCascade for the case of cascades beginning
      // with a neutralino.
      class neutralinoSet : public gluinoOrElectroweakinoSet
      {
      public:
        neutralinoSet( inputHandler const* const inputShortcut,
                       particlePointer const initialSparticle,
           electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource,
                       fullCascadeSetOrderer* const setOrderer,
                       double const beamEnergy );
        virtual
        ~neutralinoSet();


      protected:
        minimalAllocationVector< compoundGluinoOrNeutralinoFullCascade >
        compoundCascades;

        virtual void
        clearCompoundCascades();
        // this should clear the minimalAllocationVector in the derived class.
        virtual fullCascade*
        getNewCompoundCascade( fullCascade* const subcascadePointer );
        // this should add the appropriate new cascade to the
        // minimalAllocationVector in the derived class.
      };


      // this class specializes fullCascade for the case of cascades beginning
      // with a chargino.
      class charginoSet : public gluinoOrElectroweakinoSet
      {
      public:
        charginoSet( inputHandler const* const inputShortcut,
                     particlePointer const initialSparticle,
           electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource,
                       fullCascadeSetOrderer* const setOrderer,
                       double const beamEnergy );
        virtual
        ~charginoSet();


      protected:
        minimalAllocationVector< fullCascadeType::charginoToCompound >
        compoundCascades;

        virtual void
        clearCompoundCascades();
        // this should clear the minimalAllocationVector in the derived class.
        virtual fullCascade*
        getNewCompoundCascade( fullCascade* const subcascadePointer );
        // this should add the appropriate new cascade to the
        // minimalAllocationVector in the derived class.
      };

    }  // end of gluinoOrElectroweakinoSetType namespace

  }  // end of fullCascadeSetType namespace


  // this class holds all the fullCascadeSets for a given beam energy & gives
  // out pointers to them.
  class fullCascadeSetsForOneBeamEnergy
  {
  public:
    fullCascadeSetsForOneBeamEnergy( inputHandler const* const inputShortcut,
           electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource,
                                     double const beamEnergy );
    ~fullCascadeSetsForOneBeamEnergy();

    fullCascadeSet*
    getFullCascadeSet( particlePointer const initialScolored );
    /* this returns the fullCascadeSet for the requested colored sparticle, or
     * NULL if we were asked for a sparticle that is not the gluino or in
     * inputShortcut->getSquarks() or ->getElectroweakinos() (though not set up
     * for using compound decays of electroweakinos directly (used just by the
     * cascades themselves for building compound cascades, with the intention
     * of only cascades beginning with colored sparticles being requested,
     * there seems to be no problem in letting them be given out).
     */
    int
    getBeamEnergy()
    const;


  protected:
    inputHandler const* const inputShortcut;
    electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource;
    double const beamEnergy;
    fullCascadeSetOrderer cascadeOrderer;
    fullCascadeSetType::gluinoOrElectroweakinoSetType::gluinoSet
    gluinoCascadeSet;
    std::list< fullCascadeSet* >* cascadeSetList;
  };


  // this class holds all the fullCascadeSetsForOneBeamEnergy instances & gives
  // out pointers to fullCascadeSets.
  class fullCascadeSetFactory
  {
  public:
    fullCascadeSetFactory( inputHandler const* const inputShortcut,
                   electroweakCascadeHandler* const electroweakCascadeSource );
    ~fullCascadeSetFactory();

    fullCascadeSetsForOneBeamEnergy*
    getFullCascadeSetsForOneBeamEnergy( int const beamEnergy );
    /* this looks to see if there is an existing
     * fullCascadeSetsForOneBeamEnergy with the requested values, & if not,
     * makes 1, & returns the pointer.
     */
    fullCascadeSet*
    getFullCascadeSet( particlePointer const initialScolored,
                       int const beamEnergy );


  protected:
    inputHandler const* const inputShortcut;
    electroweakCascadeHandler* const electroweakCascadeSource;
    std::vector< fullCascadeSetsForOneBeamEnergy* > cascadeSetsPerEnergy;
  };





  // inline functions:


  inline particlePointer
  fullCascade::getInitialSparticle()
  const
  {
    return initialSparticle;
  }

  inline particlePointer
  fullCascade::getElectroweakinoAtEndOfScoloredness()
  const
  {
    return subcascadePointer->getElectroweakinoAtEndOfScoloredness();
  }

  inline void
  fullCascade::resetCachedBranchingRatio()
  /* this sets branchingRatioNeedsToBeRecalculated to true, so that
   * cachedBranchingRatio will be re-calculated next time the branching ratio
   * is requested. I'd prefer to do this with inheritance from
   * getsReadiedForNewPoint, but I'd have to write a variant of
   * minimalAllocationVector to allow for arguments to be given to
   * constructors so that the appropriate readierForNewPoint can be given
   * to the fullCascade constructor.
   */
  {
    branchingRatioNeedsToBeRecalculated = true;
  }

  inline double
  fullCascade::getTotalBrToEwino()
  /* this works out the branching ratio for the decays of colored
   * sparticles down to the 1st electroweakino in the cascade (including
   * any decays to bosons + squarks), without excluding any particles.
   */
  {
    if( branchingRatioNeedsToBeRecalculated )
    {
      cachedBranchingRatio
      = initialSparticle->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                                       &soughtDecayProductList,
                                               inputShortcut->getEmptyList() );
      if( NULL != subcascadePointer )
      {
        cachedBranchingRatio *= subcascadePointer->getTotalBrToEwino();
      }
      branchingRatioNeedsToBeRecalculated = false;
    }
    return cachedBranchingRatio;
  }

  inline double
  fullCascade::getBrToEwino( std::list< int > const* excludedSmParticles )
  /* this works out the branching ratio for the decays of colored sparticles
   * down to the 1st electroweakino in the cascade (including any decays to
   * EWSB bosons + squarks), using the given list of excluded SM particles
   * to exclude unwanted parts of the BR. by default, it multiplies the BR
   * for the initial decay by subcascadePointer->getBrToEwino.
   */
  {
    double
    returnValue( initialSparticle->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                                       &soughtDecayProductList,
                                                       excludedSmParticles ) );
    if( NULL != subcascadePointer )
    {
      returnValue *= subcascadePointer->getBrToEwino( excludedSmParticles );
    }
    if( initialSparticle->counts_as_self_conjugate() )
    {
      returnValue += returnValue;
    }
    return returnValue;
  }

  inline std::vector< fullCascade::particleWithInt* > const*
  fullCascade::getCascadeDefiner()
  const
  {
    return cascadeDefiner.inspectVector();
  }

  inline fullCascade::colorfulCascadeType
  fullCascade::getColorfulCascadeType()
  const
  {
    return typeOfColorfulCascade;
  }

  inline double
  fullCascade::getAcceptanceWithMinimumJets(
                                  bool const initialSparticleIsNotAntiparticle,
                                        acceptanceCutSet* const acceptanceCuts,
                                       int const minimumNumberOfAdditionalJets,
                                           int const numberOfNegativeElectrons,
                                           int const numberOfPositiveElectrons,
                                             int const numberOfNegativeMuons,
                                             int const numberOfPositiveMuons )
  /* this sums up getAcceptance for
   * numberOfAdditionalJets = minimumNumberOfAdditionalJets
   * to ( numberOfSmFermionsFromElectroweakDecaysPerFullCascade
   *      - numberOfNegativeElectrons
   *      - numberOfPositiveElectrons
   *      - numberOfNegativeMuons
   *      - numberOfPositiveMuons )
   */
  {
    double returnDouble( 0.0 );
    for( int numberOfJets( ( maximumSmFermionsFromElectroweakCascades
                             - numberOfNegativeElectrons
                             - numberOfPositiveElectrons
                             - numberOfNegativeMuons
                             - numberOfPositiveMuons ) );
         minimumNumberOfAdditionalJets <= numberOfJets;
         --numberOfJets )
    {
      returnDouble += getAcceptance( initialSparticleIsNotAntiparticle,
                                     acceptanceCuts,
                                     numberOfJets,
                                     numberOfNegativeElectrons,
                                     numberOfPositiveElectrons,
                                     numberOfNegativeMuons,
                                     numberOfPositiveMuons );
    }
    return returnDouble;
  }

  inline double
  fullCascade::leptonFlavorSummedWithMinimumJets(
                                  bool const initialSparticleIsNotAntiparticle,
                                        acceptanceCutSet* const acceptanceCuts,
                                       int const minimumNumberOfAdditionalJets,
                                             int const numberOfNegativeLeptons,
                                            int const numberOfPositiveLeptons )
  /* this sums up getAcceptanceWithMinimumJets for numberOfNegativeLeptons &
   * numberOfPositiveLeptons split between electrons & muons in every
   * combination.
   */
  {
    double returnDouble( 0.0 );
    for( int numberOfNegativeElectrons( numberOfNegativeLeptons );
         0 <= numberOfNegativeElectrons;
         --numberOfNegativeElectrons )
    {
      for( int numberOfPositiveElectrons( numberOfPositiveLeptons );
           0 <= numberOfPositiveElectrons;
           --numberOfPositiveElectrons )
      {
      returnDouble
      += getAcceptanceWithMinimumJets( initialSparticleIsNotAntiparticle,
                                       acceptanceCuts,
                                       minimumNumberOfAdditionalJets,
                                       numberOfNegativeElectrons,
                                       numberOfPositiveElectrons,
                       ( numberOfNegativeLeptons - numberOfNegativeElectrons ),
                     ( numberOfPositiveLeptons - numberOfPositiveElectrons ) );
      }
    }
    return returnDouble;
  }

  inline void
  fullCascade::buildOn( fullCascade* const copySource )
  /* this copies the basic stuff (inputShortcut, beamEnergy, cascadeDefiner)
   * from copySource, though not typeOfColorfulCascade (which is always set by
   * constructors) or initialSparticle (which depends on the type of derived
   * class), & also it sets up soughtDecayProductList to look for
   * copySource->initialSparticle, & adds copySource->initialSparticle at the
   * end of cascadeDefiner, with firstDecayBodyNumber.
   */
  {
    this->inputShortcut = copySource->inputShortcut;
    this->beamEnergy = copySource->beamEnergy;
    this->subcascadePointer = copySource;
    this->soughtDecayProductList.front()
    = copySource->initialSparticle->get_PDG_code();
    this->cascadeDefiner.clearEntries();
    for( int cascadeFiller( 0 );
         copySource->cascadeDefiner.getSize() > cascadeFiller;
         ++cascadeFiller )
    {
      this->cascadeSegment = this->cascadeDefiner.addNewAtEnd();
      this->cascadeSegment->first
      = copySource->cascadeDefiner.getPointer( cascadeFiller )->first;
      this->cascadeSegment->second
      = copySource->cascadeDefiner.getPointer( cascadeFiller )->second;
    }
    this->cascadeSegment = cascadeDefiner.addNewAtEnd();
    this->cascadeSegment->first = copySource->initialSparticle;
    this->cascadeSegment->second = firstDecayBodyNumber;
  }

  inline std::string
  fullCascade::getAsString()
  const
  // this is for error messages and testing output.
  {
    std::stringstream stringBuilder;
    stringBuilder
    << *(initialSparticle->get_name());
    for( std::vector< particleWithInt* >::const_reverse_iterator
         cascadeParticleIterator( cascadeDefiner.inspectVector()->rbegin() );
         cascadeDefiner.inspectVector()->rend() > cascadeParticleIterator;
         ++cascadeParticleIterator )
    {
      stringBuilder
      << " => "
      << *((*cascadeParticleIterator)->first->get_name());
    }
    return stringBuilder.str();
  }

  inline int
  fullCascade::getMaximumSmFermionsFromElectroweakCascades()
  const
  {
    return maximumSmFermionsFromElectroweakCascades;
  }



  namespace fullCascadeType
  {
    inline double
    squarkDirectlyToElectroweak::getBrToEwino(
                                  std::list< int > const* excludedSmParticles )
    /* this works out the branching ratio for the decays of colored
     * sparticles down to the 1st electroweakino in the cascade (including
     * any decays to bosons + squarks), using the given list of excluded SM
     * particles to exclude unwanted parts of the BR.
     */
    {
      return initialSparticle->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                                       &soughtDecayProductList,
                                                         excludedSmParticles );
    }

    inline particlePointer
    squarkDirectlyToElectroweak::getElectroweakinoAtEndOfScoloredness()
    const
    {
      return directEwinoCascades->getElectroweakDecayer();
    }

    inline double
    squarkDirectlyToElectroweak::unspecifiedJetsSpecifiedOssfMinusOsdfPairs(
                                        acceptanceCutSet* const acceptanceCuts,
                                                int const numberOfLeptonPairs )
    // this should add up all appropriate OSSF-OSDF pair combinations.
    {
      if( 0 <= numberOfLeptonPairs )
      {
        return 1.0;
      }
      else if( 1 == numberOfLeptonPairs )
      {
        return directEwinoCascades->getOssfMinusOsdf( acceptanceCuts );
      }
      else
      {
        return 0.0;
      }
    }



    namespace squarkDirectlyToElectroweakType
    {
      inline bool
      sdownType::isOpen()
      // this returns true if the squark is heavy enough to decay into the
      // electroweakino, false otherwise.
      {
        if( ( initialSparticle->get_absolute_mass()
              > directEwinoCascades->getElectroweakDecayer(
                                                       )->get_absolute_mass() )
            &&
            ( lhcFaserGlobal::negligibleBr < getTotalBrToEwino() ) )
        {
          return true;
        }
        else
        {
          return false;
        }
      }

      inline double
      sdownType::specifiedJetsOneOssfMinusOsdfPair(
                                        acceptanceCutSet* const acceptanceCuts,
                                             int const numberOfAdditionalJets )
      // this should add up all combinations of jets while looking only for
      // one OSSF-OSDF pair.
      {
        /* if there is no possibility of extra jets from a boson decay, we can
         * only get a single OSSF-OSDF pair from an electroweakino decay, which
         * cannot produce additional jets.
         */
        if( 0 <= numberOfAdditionalJets )
        {
          return directEwinoCascades->getOssfMinusOsdf( acceptanceCuts );
        }
        else
        {
          return 0.0;
        }
      }



      inline bool
      supType::decayWithWIsNotNegligible()
      // this returns true if the decay involving a W boson is not negligible,
      // also setting up the relevant fractions, false otherwise.
      {
        wFraction
        = ( initialSparticle->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                                  &decayProductListIncludingW )
            / initialSparticle->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                                   &soughtDecayProductList ) );
        if( lhcFaserGlobal::negligibleBr < wFraction )
        {
          directFraction = ( 1.0 - wFraction );
          return true;
        }
        else
        {
          return false;
        }
      }

    }  // end of squarkDirectlyToElectroweakType namespace



    inline double
    gluinoDirectlyToElectroweak::getBrToEwino(
                                  std::list< int > const* excludedSmParticles )
    /* this works out the branching ratio for the decays of colored
     * sparticles down to the 1st electroweakino in the cascade (including
     * any decays to bosons + squarks), using the given list of excluded SM
     * particles to exclude unwanted parts of the BR.
     */
    {
      double
      returnValue( initialSparticle->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                                       &soughtDecayProductList,
                                                       excludedSmParticles ) );
      if( !(ewinoCascades->getElectroweakDecayer(
                                               )->counts_as_self_conjugate()) )
        // if we need to add the charge-conjugate BR too...
      {
        returnValue += returnValue;
      }
      return returnValue;
    }

    inline particlePointer
    gluinoDirectlyToElectroweak::getElectroweakinoAtEndOfScoloredness()
    const
    {
      return ewinoCascades->getElectroweakDecayer();
    }

    inline double
    gluinoDirectlyToElectroweak::unspecifiedJetsSpecifiedOssfMinusOsdfPairs(
                                        acceptanceCutSet* const acceptanceCuts,
                                                int const numberOfLeptonPairs )
    // this should add up all appropriate OSSF-OSDF pair combinations.
    {
      if( 0 <= numberOfLeptonPairs )
      {
        return 1.0;
      }
      else if( 1 == numberOfLeptonPairs )
      {
        return ewinoCascades->getOssfMinusOsdf( acceptanceCuts );
      }
      else
      {
        return 0.0;
      }
    }

    inline double
    gluinoDirectlyToElectroweak::specifiedJetsOneOssfMinusOsdfPair(
                                        acceptanceCutSet* const acceptanceCuts,
                                             int const numberOfAdditionalJets )
    // this should add up all combinations of jets while looking only for
    // one OSSF-OSDF pair.
    {
      /* if there is no possibility of extra jets from a boson decay, we can
       * only get a single OSSF-OSDF pair from an electroweakino decay, which
       * cannot produce additional jets.
       */
      if( 0 <= numberOfAdditionalJets )
      {
        return ewinoCascades->getOssfMinusOsdf( acceptanceCuts );
      }
      else
      {
        return 0.0;
      }
    }

    inline bool
    gluinoDirectlyToElectroweak::isOpen()
    // this returns true if the squark is heavy enough to decay into the
    // electroweakino, false otherwise.
    {
      if( ( initialSparticle->get_absolute_mass()
            > ewinoCascades->getElectroweakDecayer()->get_absolute_mass() )
          &&
          ( lhcFaserGlobal::negligibleBr < getTotalBrToEwino() ) )
      {
        return true;
      }
      else
      {
        return false;
      }
    }



    inline double
    squarkByJetToCompound::unspecifiedJetsSpecifiedOssfMinusOsdfPairs(
                                        acceptanceCutSet* const acceptanceCuts,
                                                int const numberOfLeptonPairs )
    // this should add up all appropriate OSSF-OSDF pair combinations.
    {
      return subcascadePointer->unspecifiedJetsSpecifiedOssfMinusOsdfPairs(
                                                                acceptanceCuts,
                                                         numberOfLeptonPairs );
    }



    namespace squarkByJetToCompoundType
    {
      inline fullCascade*
      sdownType::setProperties( particlePointer const initialSquark,
                                fullCascade* const subcascadePointer,
          electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource )
      {
        initialSparticle = initialSquark;
        buildOn( subcascadePointer );
        this->maximumSmFermionsFromElectroweakCascades
        = subcascadePointer->getMaximumSmFermionsFromElectroweakCascades();
        resetCachedBranchingRatio();
        return this;
      }

      inline double
      sdownType::getAcceptance( bool const initialSparticleIsNotAntiparticle,
                                acceptanceCutSet* const acceptanceCuts,
                                int const numberOfAdditionalJets,
                                int numberOfNegativeElectrons,
                                int numberOfPositiveElectrons,
                                int numberOfNegativeMuons,
                                int numberOfPositiveMuons )
      /* this calls the appropriate functions on subcascadePointer to build the
       * required acceptance, taking into account whether the charges should be
       * swapped if scoloredIsNotAntiparticle is false.
       */
      {
        return
        subcascadePointer->getAcceptance( !initialSparticleIsNotAntiparticle,
                                          acceptanceCuts,
                                          numberOfAdditionalJets,
                                          numberOfNegativeElectrons,
                                          numberOfPositiveElectrons,
                                          numberOfNegativeMuons,
                                          numberOfPositiveMuons );
        /* !initialSparticleIsNotAntiparticle is used because a subsequent
         * gluino or neutralino doesn't care, while a subsquent chargino needs
         * to know that we want the case of a negatively-charged chargino.
         */
      }

      inline double
      sdownType::specifiedJetsOneOssfMinusOsdfPair(
                                        acceptanceCutSet* const acceptanceCuts,
                                             int const numberOfAdditionalJets )
      // this should add up all combinations of jets while looking only for
      // one OSSF-OSDF pair.
      {
        /* if there is no possibility of extra jets from a boson decay, we can
         * only get a single OSSF-OSDF pair from an electroweakino decay, which
         * cannot produce additional jets.
         */
        if( 0 <= numberOfAdditionalJets )
        {
          return subcascadePointer->specifiedJetsOneOssfMinusOsdfPair(
                                                                acceptanceCuts,
                                                      numberOfAdditionalJets );
        }
        else
        {
          return 0.0;
        }
      }



      inline fullCascade*
      supType::setProperties( particlePointer const initialSquark,
                              fullCascade* const subcascadePointer,
          electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource )
      {
        initialSparticle = initialSquark;
        this->electroweakCascadeSource = electroweakCascadeSource;
        buildOn( subcascadePointer );
        maximumSmFermionsFromElectroweakCascades
        = ( subcascadePointer->getMaximumSmFermionsFromElectroweakCascades()
            + 2 );
        resetCachedBranchingRatio();
        decayProductListIncludingW.front()
        = subcascadePointer->getInitialSparticle()->get_PDG_code();
        shouldUseDecaysWithW = decayWithWIsNotNegligible();
        return this;
      }

    }  // end of squarkByJetToCompoundType namespace



    inline double
  gluinoOrElectroweakinoToCompound::unspecifiedJetsSpecifiedOssfMinusOsdfPairs(
                                        acceptanceCutSet* const acceptanceCuts,
                                                int const numberOfLeptonPairs )
    // this should add up all appropriate OSSF-OSDF pair combinations.
    {
      return subcascadePointer->unspecifiedJetsSpecifiedOssfMinusOsdfPairs(
                                                                acceptanceCuts,
                                                         numberOfLeptonPairs );
    }

  }  // end of fullCascadeType namespace



  inline bool
  fullCascadeSet::massOrdered( fullCascadeSet* firstSet,
                               fullCascadeSet* secondSet )
  // this returns true if firstSet's initialSparticle has equal or lighter mass
  // than secondSet, false otherwise.
  {
    if( firstSet->initialSparticle->get_absolute_mass()
        <= secondSet->initialSparticle->get_absolute_mass() )
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  inline std::vector< fullCascade* >*
  fullCascadeSet::getOpenCascades()
  // this calls setUpCascades() if it needs to be readied for this point.
  {
    if( needsToPrepareForThisPoint() )
    {
      openCascades.clear();
      initialSparticleMass = initialSparticle->get_absolute_mass();
      setUpCascades();
    }
    return &openCascades;
  }

  inline particlePointer
  fullCascadeSet::getInitialSparticle()
  const
  {
    return initialSparticle;
  }



  inline std::list< fullCascadeSet* >*
  fullCascadeSetOrderer::getNeutralinoColoredCascades()
  {
    if( needsToPrepareForThisPoint() )
    {
      orderCascades();
    }
    return &relevantNeutralinoColoredCascades;
  }

  inline std::list< fullCascadeSet* >*
  fullCascadeSetOrderer::getCharginoColoredCascades()
  {
    if( needsToPrepareForThisPoint() )
    {
      orderCascades();
    }
    return &relevantCharginoColoredCascades;
  }

  inline std::list< fullCascadeSet* >*
  fullCascadeSetOrderer::getSdownTypeCascades()
  {
    if( needsToPrepareForThisPoint() )
    {
      orderCascades();
    }
    return &sdownTypeColoredCascades;
  }

  inline std::list< fullCascadeSet* >*
  fullCascadeSetOrderer::getSupTypeCascades()
  {
    if( needsToPrepareForThisPoint() )
    {
      orderCascades();
    }
    return &supTypeColoredCascades;
  }

  inline std::list< fullCascadeSet* >*
  fullCascadeSetOrderer::getNeutralinoColoredCascadesWithoutOrdering()
  {
    return &allNeutralinoColoredCascades;
  }

  inline std::list< fullCascadeSet* >*
  fullCascadeSetOrderer::getCharginoColoredCascadesWithoutOrdering()
  {
    return &allCharginoColoredCascades;
  }

  inline std::list< fullCascadeSet* >*
  fullCascadeSetOrderer::getSdownTypeCascadesWithoutOrdering()
  {
    return &sdownTypeColoredCascades;
  }

  inline std::list< fullCascadeSet* >*
  fullCascadeSetOrderer::getSupTypeCascadesWithoutOrdering()
  {
    return &supTypeColoredCascades;
  }



  namespace fullCascadeSetType
  {
    inline void
    squarkSet::addEwinoCompoundCascade( bool ewinoCodeIsPositive )
    // this relies on orderedCascadeSets having already been set correctly.
    {
      setIterator = orderedCascadeSets->begin();
      while( ( setIterator != orderedCascadeSets->end() )
             &&
             ( initialSparticleMass
               > (*setIterator)->getInitialSparticle()->get_absolute_mass() ) )
        /* we go through the list of given electroweakinos heavier than the
         * lightest squark until we get to either the end of the list or a
         * neutralino heavier than this down-type squark:
         */
      {
        if( ewinoCodeIsPositive )
        {
          singleSpecifiedDecayProductList.front()
          = (*setIterator)->getInitialSparticle()->get_PDG_code();
        }
        else
        {
          singleSpecifiedDecayProductList.front()
          = -((*setIterator)->getInitialSparticle()->get_PDG_code());
        }
        potentialSubcascadeSet = *setIterator;
        addByJetCascade();
        // now we move on to the next potential electroweakino cascade:
        ++setIterator;
      }
    }



    namespace squarkSetType
    {
      inline void
      sdownType::buildSquarkCompoundCascades()
      {
        // now we look at compound cascades from down-type squarks:
        orderedCascadeSets = setOrderer->getSdownTypeCascades();
        potentialDecayProducts
        = inputShortcut->getNeutralEwsbBosonsAndMassiveVectorBosons();
        addSquarkCompoundCascade( true );

        // now we look at compound cascades from up-type squarks:
        orderedCascadeSets = setOrderer->getSupTypeCascades();
        potentialDecayProducts
        = inputShortcut->getChargedEwsbBosonsAndMassiveVectorBosons();
        addSquarkCompoundCascade( false );
      }

      inline void
      sdownType::clearCompoundByJetCascades()
      // this should clear the minimalAllocationVector in the derived class.
      {
        compoundByJetCascades.clearEntries();
      }

      inline fullCascade*
      sdownType::getNewCompoundByJetCascade(
                                         fullCascade* const subcascadePointer )
      // this should add the appropriate new cascade to the
      // minimalAllocationVector in the derived class.
      {
        return
        compoundByJetCascades.addNewAtEnd()->setProperties( initialSparticle,
                                                            subcascadePointer,
                                                    electroweakCascadeSource );
      }

      inline void
      sdownType::findOpenDirectCascades()
      // this puts all open direct cascades into openCascades.
      {
        for( std::vector< squarkDirectlyToElectroweakCascade* >::iterator
             cascadeIterator( directToEwinoCascades.begin() );
             directToEwinoCascades.end() > cascadeIterator;
             ++cascadeIterator )
        {
          if( (*cascadeIterator)->isOpen() )
          {
            openCascades.push_back( *cascadeIterator );
          }
        }
      }



      inline void
      supType::buildSquarkCompoundCascades()
      {
        // now we look at compound cascades from up-type squarks:
        orderedCascadeSets = setOrderer->getSupTypeCascades();
        potentialDecayProducts
        = inputShortcut->getNeutralEwsbBosonsAndMassiveVectorBosons();
        addSquarkCompoundCascade( true );

        // now we look at compound cascades from down-type squarks:
        orderedCascadeSets = setOrderer->getSdownTypeCascades();
        potentialDecayProducts
        = inputShortcut->getChargedEwsbBosonsAndMassiveVectorBosons();
        addSquarkCompoundCascade( true );
      }

      inline void
      supType::clearCompoundByJetCascades()
      // this should clear the minimalAllocationVector in the derived class.
      {
        compoundByJetCascades.clearEntries();
      }

      inline fullCascade*
      supType::getNewCompoundByJetCascade(
                                         fullCascade* const subcascadePointer )
      // this should add the appropriate new cascade to the
      // minimalAllocationVector in the derived class.
      {
        return
        compoundByJetCascades.addNewAtEnd()->setProperties( initialSparticle,
                                                            subcascadePointer,
                                                    electroweakCascadeSource );
      }

      inline void
      supType::findOpenDirectCascades()
      // this puts all open direct cascades into openCascades.
      {
        for( std::vector< squarkDirectlyToElectroweakCascade* >::iterator
             cascadeIterator( directToEwinoCascades.begin() );
             directToEwinoCascades.end() > cascadeIterator;
             ++cascadeIterator )
        {
          if( (*cascadeIterator)->isOpen() )
          {
            openCascades.push_back( *cascadeIterator );
          }
        }
      }

    }  // end of squarkSetType namespace



    inline void
    gluinoOrElectroweakinoSet::setUpCascades()
    {
      // 1st we clear the compound cascades:
      clearCompoundCascades();

      // now we look at compound cascades from squarks:
      orderedCascadeSets = setOrderer->getSdownTypeCascades();
      buildSquarkCompoundCascades();
      orderedCascadeSets = setOrderer->getSupTypeCascades();
      buildSquarkCompoundCascades();
    }



    namespace gluinoOrElectroweakinoSetType
    {
      inline fullCascade*
      gluinoSet::getNewCompoundCascade( fullCascade* const subcascadePointer )
      // this should add the appropriate new cascade to the
      // minimalAllocationVector in the derived class.
      {
        return compoundCascades.addNewAtEnd()->setProperties( initialSparticle,
                                                             subcascadePointer,
                                                    electroweakCascadeSource );
      }


      inline void
      gluinoSet::clearCompoundCascades()
      // this should clear the minimalAllocationVector in the derived class.
      {
        compoundCascades.clearEntries();
      }



      inline fullCascade*
      neutralinoSet::getNewCompoundCascade(
                                         fullCascade* const subcascadePointer )
      // this should add the appropriate new cascade to the
      // minimalAllocationVector in the derived class.
      {
        return compoundCascades.addNewAtEnd()->setProperties( initialSparticle,
                                                             subcascadePointer,
                                                    electroweakCascadeSource );
      }


      inline void
      neutralinoSet::clearCompoundCascades()
      // this should clear the minimalAllocationVector in the derived class.
      {
        compoundCascades.clearEntries();
      }



      inline fullCascade*
      charginoSet::getNewCompoundCascade(
                                         fullCascade* const subcascadePointer )
      // this should add the appropriate new cascade to the
      // minimalAllocationVector in the derived class.
      {
        return compoundCascades.addNewAtEnd()->setProperties( initialSparticle,
                                                             subcascadePointer,
                                                    electroweakCascadeSource );
      }


      inline void
      charginoSet::clearCompoundCascades()
      // this should clear the minimalAllocationVector in the derived class.
      {
        compoundCascades.clearEntries();
      }

    }  // end of gluinoOrElectroweakinoSetType namespace

  }  // end of fullCascadeSetType namespace



  inline int
  fullCascadeSetsForOneBeamEnergy::getBeamEnergy()
  const
  {
    return beamEnergy;
  }



  inline fullCascadeSet*
  fullCascadeSetFactory::getFullCascadeSet(
                                         particlePointer const initialScolored,
                                            int const beamEnergy )
  {
    return getFullCascadeSetsForOneBeamEnergy( beamEnergy
                                       )->getFullCascadeSet( initialScolored );
  }

}  // end of LHC_FASER namespace.

#endif /* LHC_FASER_FULL_CASCADE_STUFF_HPP_ */
