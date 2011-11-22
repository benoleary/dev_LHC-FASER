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

#ifndef LHC_FASER_FULL_CASCADE_STUFF_HPP_
#define LHC_FASER_FULL_CASCADE_STUFF_HPP_

#include "LHC-FASER_input_handling_stuff.hpp"
#include "LHC-FASER_sparticle_decay_stuff.hpp"
#include "LHC-FASER_lepton_kinematics_stuff.hpp"
#include "LHC-FASER_electroweak_cascade_collection_stuff.hpp"

namespace LHC_FASER
{
  /* new plan:
   * basic squarkToEwino, gluinoToEwino
   * => sxFullCascade, gxFullCascade, & sups have supxFullCascade,
   * which has an electroweakCascadeSet for decays without a W, & a pair of
   * electroweakCascadeSets for decays with a W.
   * gjm, sjm, gbm, sbm (m for more) build on above + also on each other.
   * es for intermediate electroweakinos.
   * put all sq & all unstable ewinos into list, order by mass, pop out
   * ewinos lighter than all sq & all ewinos heavier than all sq.
   * while( m_sq < m_go ), makeSquarksLighterThanGluino, then makeGluino,
   * then makeSquarksHeavierThanGluino, making es objects as encountered.
   *
   * squarkToEwino, gluinoToEwino both have basic getAcceptance. sm & gm build
   * on getAcceptance of subcascadePointer, sbm & gbm build on own
   * bosonCascades * getAcceptance of subcascadePointer.
   *
   * need to write supDecayToWPlusAndNeutralino class derived from
   * electroweakCascadeSet class. should just have pointers to 2
   * electroweakCascades but giving virtualSdown effectiveSquarkMasses.
   * (reorganize electroweakCascadeSet so that constructor uses initialization
   * functions that can be used by supDecayToWPlusAndNeutralino.)
   */


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
      gluinoToSquarkThenMore /* gaugino to squark then shorter cascade */
      /* treated by jets+MET grid as gluinoToSquarkToEwino */,
      ewinoToMore /* electroweakino to shorter cascade */
      /* shouldn't be seen by jets+MET grid */,
      sizeOfEnumPlusOne /* possibly useful */
    };

    fullCascade( colorfulCascadeType const typeOfColorfulCascade,
                 int const firstDecayBodyNumber,
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
    getBrToEwino( std::list< int > const* excludedSmParticles )
    = 0;
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
    specifiedJetsOneOssfMinusOsdfPair( acceptanceCutSet* const acceptanceCuts,
                                       int const numberOfAdditionalJets );
    /* this calls
     * ewinoCascades->getOssfMinusOsdf( initialSparticleIsNotAntiparticle,
     *                                  acceptanceCuts ),
     * & if bosonCascades is not NULL, it calls
     * ewinoCascades->getAcceptance( [ { 0, 1, 2 } jets
     *                                 + 0 leptons for acceptanceCuts ] )
     * & the same for bosonCascades, as well as
     * ->getOssfMinusOsdf( acceptanceCuts ).
     */
    double
    unspecifiedJetsSpecifiedOssfMinusOsdfPairs(
                                        acceptanceCutSet* const acceptanceCuts,
                                               int const numberOfLeptonPairs );
    /* this only calls ewinoCascades->getOssfMinusOsdf( acceptanceCuts ),
     * & if bosonCascades is not NULL,
     * bosonCascades->getOssfMinusOsdf( acceptanceCuts ).
     */

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

  protected:
    static int const numberOfSmFermionsFromElectroweakDecaysPerPartialCascade;
    static int const numberOfSmFermionsFromElectroweakDecaysPerFullCascade;

    inputHandler const* inputShortcut;
    particlePointer initialSparticle;
    double beamEnergy;
    colorfulCascadeType const typeOfColorfulCascade;
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
    // this is derived class for direct sdown-type to electroweakino
    // fullCascades.
    class sdownDirectlyToElectroweak : public fullCascade
    {
    public:
      sdownDirectlyToElectroweak( inputHandler const* const inputShortcut,
                                  particlePointer const initialSquark,
                                  double const beamEnergy,
                                  electroweakCascadeSet* const ewinoCascades );
      virtual
      ~sdownDirectlyToElectroweak();

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
      /* this calls the appropriate functions on ewinoCascades to build the
       * required acceptance, taking into account whether the charges should be
       * swapped if scoloredIsNotAntiparticle is false.
       */

    protected:
      bool ewinoFlipsCharge;
      electroweakCascadeSet* ewinoCascades;
    };


    // this is derived class for direct sup-type to electroweakino
    // fullCascades, which also covers if the squark decays to an additional W.
    class supDirectlyToElectroweak : public fullCascade
    {
    public:
      supDirectlyToElectroweak( inputHandler const* const inputShortcut,
                                particlePointer const initialSquark,
                                double const beamEnergy,
                              electroweakCascadeSet* const directEwinoCascades,
                               electroweakCascadeSet* const ewinoWithWCascades,
                                electroweakCascadeSet* const bosonCascades );
      virtual
      ~supDirectlyToElectroweak();

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
      /* this calls the appropriate functions on directEwinoCascades to build
       * the required acceptance or ewinoWithWCascades & bosonCascades if the
       * decay with a boson is not negligible, taking into account whether the
       * charges should be swapped if scoloredIsNotAntiparticle is false.
       */

    protected:
      electroweakCascadeSet* directEwinoCascades;
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
      // this returns true if the decay involving a W boson is not negligible,
      // also setting up the relevant fractions, false otherwise.
      double
      getCombinedAcceptance( acceptanceCutSet* const acceptanceCuts,
                             int const numberOfAdditionalJets,
                             int const numberOfNegativeElectrons,
                             int const numberOfPositiveElectrons,
                             int const numberOfNegativeMuons,
                             int const numberOfPositiveMuons );
    };


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

    protected:
      electroweakCascadeSet* bosonCascades;
      bool bosonFlipsCharge;
    };


    /* this is derived class for compound fullCascades beginning with a
     * sdown-type, which decays to a jet & the initial sparticle of
     * subcascadePointer.
     */
    class sdownByJetToCompound : public fullCascade
    {
    public:
      sdownByJetToCompound();
      virtual
      ~sdownByJetToCompound();

      virtual fullCascade*
      setProperties( particlePointer const initialSquark,
                     fullCascade* const subcascadePointer );
      virtual double
      getAcceptance( bool const initialSparticleIsNotAntiparticle,
                     acceptanceCutSet* const acceptanceCuts,
                     int const numberOfAdditionalJets,
                     int numberOfNegativeElectrons,
                     int numberOfPositiveElectrons,
                     int numberOfNegativeMuons,
                     int numberOfPositiveMuons );
      /* this calls the appropriate functions on subcascadePointer to build the
       * required acceptance, taking into account whether the charges should be
       * swapped if scoloredIsNotAntiparticle is false.
       */

    //protected:
      // nothing.
    };


    /* this is derived class for compound fullCascades beginning with a
     * sup-type, which decays to a jet, maybe a boson, & the initial sparticle
     * of subcascadePointer.
     */
    class supByJetToCompound : public fullCascade
    {
    public:
      supByJetToCompound();
      virtual
      ~supByJetToCompound();

      virtual fullCascade*
      setProperties( particlePointer const initialSquark,
                     fullCascade* const subcascadePointer );
      virtual double
      getAcceptance( bool const initialSparticleIsNotAntiparticle,
                     acceptanceCutSet* const acceptanceCuts,
                     int const numberOfAdditionalJets,
                     int numberOfNegativeElectrons,
                     int numberOfPositiveElectrons,
                     int numberOfNegativeMuons,
                     int numberOfPositiveMuons );
      /* this calls the appropriate functions on subcascadePointer to build the
       * required acceptance, & also with bosonCascades if the decay with a
       * boson is not negligible, taking into account whether the charges
       * should be swapped if scoloredIsNotAntiparticle is false.
       */

    protected:
      electroweakCascadeSet* const bosonCascades;
      particlePointer wBoson;
      double ewinoMass;
      bool shouldUseDecaysWithW;
      double directFraction;
      double wFraction;
      std::list< int > decayProductListIncludingW;

      bool
      decayWithWIsNotNegligible();
      // this returns true if the decay involving a W boson is not negligible,
      // also setting up the relevant fractions, false otherwise.
      double
      getCombinedAcceptance( acceptanceCutSet* const acceptanceCuts,
                             int const numberOfAdditionalJets,
                             int const numberOfNegativeElectrons,
                             int const numberOfPositiveElectrons,
                             int const numberOfNegativeMuons,
                             int const numberOfPositiveMuons );
    };


    /* this is derived class for compound fullCascades beginning with a
     * gluino or neutralino, which decays to a jet, maybe a boson & the initial
     * sparticle of subcascadePointer.
     */
    class gluinoOrNeutralinoToCompound : public fullCascade
    {
    public:
      gluinoOrNeutralinoToCompound();
      virtual
      ~gluinoOrNeutralinoToCompound();

      virtual fullCascade*
      setProperties( particlePointer const initialSparticle,
                     fullCascade* const subcascadePointer );
      virtual double
      getAcceptance( bool const initialSparticleIsNotAntiparticle,
                     acceptanceCutSet* const acceptanceCuts,
                     int const numberOfAdditionalJets,
                     int numberOfNegativeElectrons,
                     int numberOfPositiveElectrons,
                     int numberOfNegativeMuons,
                     int numberOfPositiveMuons );
      /* this calls the appropriate functions on subcascadePointer to build the
       * required acceptance, & also with bosonCascades if the decay with a
       * boson is not negligible, taking into account whether the charges
       * should be swapped if scoloredIsNotAntiparticle is false.
       */

    protected:
      electroweakCascadeSet* const bosonCascades;
      particlePointer wBoson;
      double ewinoMass;
      bool shouldUseDecaysWithW;
      double directFraction;
      double wFraction;
      std::list< int > decayProductListIncludingW;

      bool
      decayWithWIsNotNegligible();
      // this returns true if the decay involving a W boson is not negligible,
      // also setting up the relevant fractions, false otherwise.
      double
      getCombinedAcceptance( acceptanceCutSet* const acceptanceCuts,
                             int const numberOfAdditionalJets,
                             int const numberOfNegativeElectrons,
                             int const numberOfPositiveElectrons,
                             int const numberOfNegativeMuons,
                             int const numberOfPositiveMuons );
    };


    /* this is derived class for compound fullCascades beginning with a
     * chargino, which decays to a jet, maybe a boson & the initial sparticle
     * of subcascadePointer.
     */
    class charginoToCompound : public fullCascade
    {
    public:
      charginoToCompound();
      virtual
      ~charginoToCompound();

      virtual fullCascade*
      setProperties( particlePointer const initialSparticle,
                     fullCascade* const subcascadePointer );
      virtual double
      getAcceptance( bool const initialSparticleIsNotAntiparticle,
                     acceptanceCutSet* const acceptanceCuts,
                     int const numberOfAdditionalJets,
                     int numberOfNegativeElectrons,
                     int numberOfPositiveElectrons,
                     int numberOfNegativeMuons,
                     int numberOfPositiveMuons );
      /* this calls the appropriate functions on subcascadePointer to build the
       * required acceptance, & also with bosonCascades if the decay with a
       * boson is not negligible, taking into account whether the charges
       * should be swapped if scoloredIsNotAntiparticle is false.
       */

    protected:
      electroweakCascadeSet* const bosonCascades;
      particlePointer wBoson;
      double ewinoMass;
      bool shouldUseDecaysWithW;
      double directFraction;
      double wFraction;
      std::list< int > decayProductListIncludingW;
      double squarkWithWFraction;
      double antisquarkWithWFraction;
      bool decayingToSupType;

      bool
      decayWithWIsNotNegligible();
      // this returns true if the decay involving a W boson is not negligible,
      // also setting up the relevant fractions, false otherwise.
      double
      getCombinedAcceptance( acceptanceCutSet* const acceptanceCuts,
                             int const numberOfAdditionalJets,
                             int const numberOfNegativeElectrons,
                             int const numberOfPositiveElectrons,
                             int const numberOfNegativeMuons,
                             int const numberOfPositiveMuons );
      // if this is called, it is the decay of a positive chargino.
    };

  }


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
    fullCascadeSet* const gluinoFullCascade;
    std::list< fullCascadeSet* >* orderedCascadeSets;
    std::vector< fullCascade* >* potentialSubcascades;
    double subcascadeBranchingRatio;

    double const beamEnergy;
    std::list< int > singleSpecifiedDecayProductList;
    std::list< int > twoSpecifiedDecayProductsList;
    int singleSpecifiedDecayProduct;
    std::vector< particlePointer >* potentialDecayProducts;

    void
    ensureCascadesAreUpToDate();
    // this calls setUpCascades() if needsToPrepareForThisPoint() is true.
    void
    setUpCascades();
    // this clears openCascades, then sets it to be filled with all open
    // fullCascades that initialSparticle has for this point.
    virtual void
    buildLongerCascades()
    = 0;
    // this should set up all the cascades that involve taking subcascades from
    // the other fullCascadeSets.
    void
    addOpenCascade( fullCascade* cascadeToAdd );
    // this adds the cascade to openCascades & also resets it so that it'll
    // recalculate its branching ratio when next requested.
  };



  class fullCascadeSetOrderer : public getsReadiedForNewPoint
  {
  public:
    fullCascadeSetOrderer();
    ~fullCascadeSetOrderer();

    std::list< fullCascadeSet* > const*
    getOrderedFullCascadeSets();
    void
    addFullCascadeSet( fullCascadeSet* const fullCascadeSetPointer );

  protected:
    std::list< fullCascadeSet* > allFullCascadeSets;
    std::list< fullCascadeSet* > orderedAndParedList;
  };


  namespace fullCascadeSetType
  {
    // this class specializes fullCascade for the case of cascades beginning
    // with a down-type squark.
    class sdownTypeSet : public fullCascadeSet
    {
    public:
      sdownTypeSet( inputHandler const* const inputShortcut,
                   electroweakCascadesForOneBeamEnergy* const ewCascadeHandler,
                    particlePointer const initialScolored,
                    fullCascadeSetOrderer* const setOrderer,
                    double const beamEnergy );
      virtual
      ~sdownTypeSet();

    protected:
      fullCascadeSetOrderer* const setOrderer;
      std::vector< fullCascadeType::sdownDirectlyToElectroweak* >
      directToEwinoCascades;
      minimalAllocationVector< fullCascadeType::squarkByBosonToCompound >
      compoundByBosonCascades;
      minimalAllocationVector< fullCascadeType::sdownByJetToCompound >
      compoundByJetCascades;

      virtual void
      buildLongerCascades();
      // this should set up all the cascades that involve taking subcascades
      // from the other fullCascadeSets.
    };


    // this class specializes fullCascade for the case of cascades beginning
    // with an up-type squark.
    class supTypeSet : public fullCascadeSet
    {
    public:
      supTypeSet( inputHandler const* const inputShortcut,
                 electroweakCascadesForOneBeamEnergy* const ewCascadeHandler,
                 particlePointer const initialScolored,
                 readiableCascadeSetList* const squarkCascadeSetList,
                 fullCascadeSet* const gluinoCascadeSet,
                 double const beamEnergy );
      virtual
      ~supTypeSet();

    protected:
      fullCascadeSetOrderer* const setOrderer;
      std::vector< fullCascadeType::supDirectlyToElectroweak* >
      directToEwinoCascades;
      minimalAllocationVector< fullCascadeType::squarkByBosonToCompound >
      compoundByBosonCascades;
      minimalAllocationVector< fullCascadeType::supByJetToCompound >
      compoundByJetCascades;

      virtual void
      buildLongerCascades();
      // this should set up all the cascades that involve taking subcascades
      // from the other fullCascadeSets.
    };


    // this class specializes fullCascade for the case of cascades beginning
    // with a gluino.
    class gluinoSet : public fullCascadeSet
    {
    public:
      gluinoSet( inputHandler const* const inputShortcut,
                 electroweakCascadesForOneBeamEnergy* const ewCascadeHandler,
                 readiableCascadeSetList* const squarkCascadeSetList,
                 double const beamEnergy );
      virtual
      ~gluinoSet();

    protected:
      virtual void
      buildLongerCascades();
      // this should set up all the cascades that involve taking subcascades
      // from the other fullCascadeSets.
    };

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
     * inputShortcut->getSquarks().
     */
    int
    getBeamEnergy()
    const;

  protected:
    inputHandler const* const inputShortcut;
    electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource;
    double const beamEnergy;
    std::vector< sdownTypeSet* > squarkCascadeSets;
    gluinoSet* gluinoCascadeSet;
    fullCascadeSetOrderer cascadeOrderer;
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
    for( int
         numberOfJets( ( numberOfSmFermionsFromElectroweakDecaysPerFullCascade
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
         cascadeParticleIterator( cascadeDefiner->rbegin() );
         cascadeDefiner->rend() > cascadeParticleIterator;
         ++cascadeParticleIterator )
    {
      stringBuilder
      << " => "
      << *((*cascadeParticleIterator)->first->get_name());
    }
    return stringBuilder.str();
  }



  namespace fullCascadeType
  {
    inline bool
    sdownDirectlyToElectroweak::isOpen()
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
    sdownDirectlyToElectroweak::getBrToEwino(
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
    sdownDirectlyToElectroweak::getElectroweakinoAtEndOfScoloredness()
    const
    {
      return ewinoCascades->getElectroweakDecayer();
    }

    inline double
    sdownDirectlyToElectroweak::getAcceptance(
                                          bool const initialSparticleIsNotAntiparticle,
                                        acceptanceCutSet* const acceptanceCuts,
                                              int const numberOfAdditionalJets,
                                              int numberOfNegativeElectrons,
                                              int numberOfPositiveElectrons,
                                              int numberOfNegativeMuons,
                                              int numberOfPositiveMuons )
    /* this calls the appropriate functions on ewinoCascades to build the
     * required acceptance, taking into account whether the charges should be
     * swapped if scoloredIsNotAntiparticle is false.
     */
    {
      if( ( initialSparticleIsNotAntiparticle
            &&
            ewinoFlipsCharge )
          ||
          ( !initialSparticleIsNotAntiparticle
            &&
            !ewinoFlipsCharge ) )
        // if the electroweakino is an antiparticle, we swap the charges.
      {
        int swappingInt( numberOfPositiveElectrons );
        numberOfPositiveElectrons = numberOfNegativeElectrons;
        numberOfNegativeElectrons = swappingInt;
        swappingInt = numberOfPositiveMuons;
        numberOfPositiveMuons = numberOfNegativeMuons;
        numberOfNegativeMuons = swappingInt;
      }
      return ewinoCascades->getAcceptance( acceptanceCuts,
                                           numberOfAdditionalJets,
                                           numberOfNegativeElectrons,
                                           numberOfPositiveElectrons,
                                           numberOfNegativeMuons,
                                           numberOfPositiveMuons );
    }



    inline bool
    supDirectlyToElectroweak::decayWithWIsNotNegligible()
    // this returns true if the decay involving a W boson is not negligible,
    // also setting up the relevant fractions, false otherwise.
    {
      wFraction = ( initialSparticle->inspect_direct_decay_handler(
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

    inline bool
    supDirectlyToElectroweak::isOpen()
    /* this returns true if the squark is heavy enough to decay into the
     * electroweakino, false otherwise. it also sorts out whether it should
     * be using the decays involving a W boson.
     */
    {
      ewinoMass
      = directEwinoCascades->getElectroweakDecayer()->get_absolute_mass();
      if( ( initialSparticle->get_absolute_mass()
            > ewinoCascades->getElectroweakDecayer()->get_absolute_mass() )
          &&
          ( lhcFaserGlobal::negligibleBr < getTotalBrToEwino() ) )
      {
        if( initialSparticle->get_absolute_mass()
            > ( ewinoMass + wBoson->get_absolute_mass() ) )
        {
          shouldUseDecaysWithW = decayWithWIsNotNegligible();
        }
        else
        {
          shouldUseDecaysWithW = false;
        }
        return true;
      }
      else
      {
        return false;
      }
    }

    inline double
    supDirectlyToElectroweak::getBrToEwino(
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
    supDirectlyToElectroweak::getElectroweakinoAtEndOfScoloredness()
    const
    {
      return ewinoCascades->getElectroweakDecayer();
    }



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



    inline fullCascade*
    squarkByBosonToCompound::setProperties(
                                           particlePointer const initialSquark,
                                    electroweakCascadeSet* const bosonCascades,
                                         fullCascade* const subcascadePointer )
    {
      initialSparticle = initialSquark;
      buildOn( subcascadePointer );
      resetCachedBranchingRatio();
      this->bosonCascades = bosonCascades;
      if( ( CppSLHA::PDG_code::W_plus
            == bosonCascades->getElectroweakDecayer()->get_PDG_code() )
          &&
          ( inputShortcut->isIn( initialSquark->get_PDG_code(),
                                 inputShortcut->getSdownTypes() ) ) )
        // if we have to worry about which sign of PDG code to use...
      {
        bosonFlipsCharge = true;
        soughtDecayProductList.back() = -(CppSLHA::PDG_code::W_plus);
      }
      else
      {
        bosonFlipsCharge = false;
        soughtDecayProductList.back()
        = bosonCascades->getElectroweakDecayer()->get_PDG_code();
      }
      return this;
    }



    inline void
    sdownByJetToCompound::setProperties( particlePointer const initialSquark,
                                         fullCascade* const subcascadePointer )
    {
      initialSparticle = initialSquark;
      buildOn( subcascadePointer );
      resetCachedBranchingRatio();
    }

    inline double
    sdownByJetToCompound::getAcceptance(
                                  bool const initialSparticleIsNotAntiparticle,
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
      /* !initialSparticleIsNotAntiparticle is used because a subsequent gluino
       * or neutralino doesn't care, while a subsquent chargino needs to know
       * that we want the case of a negatively-charged chargino.
       */
    }



    inline fullCascade*
    supByJetToCompound::setProperties( particlePointer const initialSquark,
                                       fullCascade* const subcascadePointer,
                                   electroweakCascadeSet* const bosonCascades )
    {
      initialSparticle = initialSquark;
      this->bosonCascades = bosonCascades;
      buildOn( subcascadePointer );
      resetCachedBranchingRatio();
      decayProductListIncludingW.front()
      = subcascadePointer->getInitialSparticle()->get_PDG_code();
      decayProductListIncludingW.back()
      = bosonCascades->getElectroweakDecayer()->get_PDG_code();
      shouldUseDecaysWithW = decayWithWIsNotNegligible();
      return this;
    }



    inline bool
    gluinoOrNeutralinoToCompound::decayWithWIsNotNegligible()
    // this returns true if the decay involving a W boson is not negligible,
    // also setting up the relevant fractions, false otherwise.
    {
      wFraction = ( initialSparticle->inspect_direct_decay_handler(
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

    inline fullCascade*
    gluinoOrNeutralinoToCompound::setProperties(
                                        particlePointer const initialSparticle,
                                          fullCascade* const subcascadePointer,
                                   electroweakCascadeSet* const bosonCascades )
    {
      this->initialSparticle = initialSparticle;
      this->bosonCascades = bosonCascades;
      buildOn( subcascadePointer );
      resetCachedBranchingRatio();
      decayProductListIncludingW.front()
      = subcascadePointer->getInitialSparticle()->get_PDG_code();
      // only sup-types are going to be produced with virtual antitops by
      // gluinos or neutralinos:
      if( inputShortcut->isIn( decayProductListIncludingW.front(),
                               inputShortcut->getSupTypes() ) )
      {
        decayProductListIncludingW.back()
        = -(bosonCascades->getElectroweakDecayer()->get_PDG_code());
        shouldUseDecaysWithW = decayWithWIsNotNegligible();
      }
      else
      {
        shouldUseDecaysWithW = false;
      }
      return this;
    }



    inline bool
    charginoToCompound::decayWithWIsNotNegligible()
    // this returns true if the decay involving a W boson is not negligible,
    // also setting up the relevant fractions, false otherwise.
    {
      // 1st we look for chargino -> squark + antiquark + W:
      wFraction = ( initialSparticle->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                                  &decayProductListIncludingW )
                    / initialSparticle->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                                   &soughtDecayProductList ) );
      squarkWithWFraction = wFraction;
      // now we look for chargino -> antisquark + quark + W:
      decayProductListIncludingW.front()
      = -(decayProductListIncludingW.front());
      wFraction += ( initialSparticle->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                                  &decayProductListIncludingW )
                    / initialSparticle->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                                   &soughtDecayProductList ) );
      if( lhcFaserGlobal::negligibleBr < wFraction )
      {
        directFraction = ( 1.0 - wFraction );
        antisquarkWithWFraction = ( wFraction - squarkWithWFraction );
        return true;
      }
      else
      {
        return false;
      }
    }

    inline fullCascade*
    charginoToCompound::setProperties( particlePointer const initialSparticle,
                                       fullCascade* const subcascadePointer,
                                   electroweakCascadeSet* const bosonCascades )
    {
      this->initialSparticle = initialSparticle;
      this->bosonCascades = bosonCascades;
      buildOn( subcascadePointer );
      resetCachedBranchingRatio();
      decayProductListIncludingW.front()
      = subcascadePointer->getInitialSparticle()->get_PDG_code();
      decayProductListIncludingW.back()
      = bosonCascades->getElectroweakDecayer()->get_PDG_code();
      // positive charginos decay to positive W bosons.
      if( inputShortcut->isIn( decayProductListIncludingW.front(),
                               inputShortcut->getSdownTypes() ) )
      {
        decayingToSupType = false;
      }
      else
      {
        decayProductListIncludingW.front()
        = subcascadePointer->getInitialSparticle()->get_PDG_code();
        decayingToSupType = true;
      }
      shouldUseDecaysWithW = decayWithWIsNotNegligible();
      return this;
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
    ensureCascadesAreUpToDate();
    return &openCascades;
  }

  inline std::vector< sxFullCascade* > const*
  fullCascadeSet::getSxCascades()
  {
    ensureCascadesAreUpToDate();
    return &openSxCascades;
  }

  inline std::vector< gxFullCascade* > const*
  fullCascadeSet::getGxCascades()
  {
    ensureCascadesAreUpToDate();
    return &openGxCascades;
  }

  inline std::vector< sjgxFullCascade* > const*
  fullCascadeSet::getSjgxCascades()
  {
    ensureCascadesAreUpToDate();
    return sjgxCascades.inspectVector();
  }

  inline std::vector< gjsxFullCascade* > const*
  fullCascadeSet::getGjsxCascades()
  {
    ensureCascadesAreUpToDate();
    return gjsxCascades.inspectVector();
  }

  inline std::vector< sjgjsxFullCascade* > const*
  fullCascadeSet::getSjgjsxCascades()
  {
    ensureCascadesAreUpToDate();
    return sjgjsxCascades.inspectVector();
  }

  inline std::vector< sbsxFullCascade* > const*
  fullCascadeSet::getSbsxCascades()
  {
    ensureCascadesAreUpToDate();
    return sbsxCascades.inspectVector();
  }

  inline std::vector< gbsxFullCascade* > const*
  fullCascadeSet::getGbsxCascades()
  {
    ensureCascadesAreUpToDate();
    return gbsxCascades.inspectVector();
  }

  inline std::vector< gjsbsxFullCascade* > const*
  fullCascadeSet::getGjsbsxCascades()
  {
    ensureCascadesAreUpToDate();
    return gjsbsxCascades.inspectVector();
  }

  inline std::vector< sbgxFullCascade* > const*
  fullCascadeSet::getSbgxCascades()
  {
    ensureCascadesAreUpToDate();
    return sbgxCascades.inspectVector();
  }

  inline std::vector< sbsjgxFullCascade* > const*
  fullCascadeSet::getSbsjgxCascades()
  {
    ensureCascadesAreUpToDate();
    return sbsjgxCascades.inspectVector();
  }

  inline std::vector< sbsjgjsxFullCascade* > const*
  fullCascadeSet::getSbsjgjsxCascades()
  {
    ensureCascadesAreUpToDate();
    return sbsjgjsxCascades.inspectVector();
  }

  inline std::vector< sbgjsxFullCascade* > const*
  fullCascadeSet::getSbgjsxCascades()
  {
    ensureCascadesAreUpToDate();
    return sbgjsxCascades.inspectVector();
  }

  inline std::vector< sjgbsxFullCascade* > const*
  fullCascadeSet::getSjgbsxCascades()
  {
    ensureCascadesAreUpToDate();
    return sjgbsxCascades.inspectVector();
  }

  inline std::vector< sjgjsbsxFullCascade* > const*
  fullCascadeSet::getSjgjsbsxCascades()
  {
    ensureCascadesAreUpToDate();
    return sjgjsbsxCascades.inspectVector();
  }

  inline particlePointer
  fullCascadeSet::getInitialSparticle()
  const
  {
    return initialSparticle;
  }

  inline void
  fullCascadeSet::ensureCascadesAreUpToDate()
  // this calls setUpCascades() if needsToPrepareForThisPoint() is true.
  {
    if( needsToPrepareForThisPoint() )
    {
      setUpCascades();
    }
  }

  inline void
  fullCascadeSet::addOpenCascade( fullCascade* cascadeToAdd )
  // this adds the cascade to openCascades & also resets it so that it'll
  // recalculate its branching ratio when next requested.
  {
    openCascades.push_back( cascadeToAdd );
    cascadeToAdd->resetCachedBranchingRatio();
  }



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
