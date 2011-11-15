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
 *      LHC-FASER also requires grids of lookup acceptanceValues. These should also be
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
   * basic sx, gx, sbx, gbx
   * gm, sm (m for more) build on above + also on each other
   * es for intermediate electroweakinos
   * put all sq & all unstable ewinos into list, order by mass, pop out
   * ewinos lighter than all sq & all ewinos heavier than all sq.
   * while( m_sq < m_go ), makeSquarksLighterThanGluino, then makeGluino,
   * then makeSquarksHeavierThanGluino, making es objects as encountered.
   *
   * sx, gx, sbx, & gbx all have basic getAcceptance. sm & gm build on own (if
   * not NULL) bosonCascades * getAcceptance of subcascadePointer.
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
      /* we do not take into account multiple decays to vector or spin-0 bosons
       * plus other colored sparticles, because the cascades would be too long,
       * with too low branching ratios, as well as being quite difficult to get
       * right...
       */
      unset /* marked as not yet a proper cascade */,
      sx /* direct squark to electroweakino */,
      gx /* direct gluino to electroweakino */,
      sjgx /* squark to gluino to electroweakino */,
      gjsx /* gluino to squark to electroweakino */,
      sjgjsx
      /* squark to gluino to squark to electroweakino */
      /* treated by jets+MET grid as gsx */,
      sbsx /* squark to boson plus other squark to electroweakino */
      /* treated by jets+MET grid as sx */,
      gbsx /* gluino to boson plus squark to electroweakino */
      /* treated by jets+MET grid as sx */,
      gjsbsx
      /* gluino to squark to boson plus other squark to electroweakino */
      /* treated by jets+MET grid as gsx */,
      sbgx /* squark to boson plus gluino to electroweakino */
      /* treated by jets+MET grid as sgx */,
      sbsjgx
      /* squark to vector boson plus other squark to gluino to
         electroweakino */
      /* treated by jets+MET grid as sgx */,
      sbsjgjsx
      /* squark to boson plus other squark to gluino to yet another squark to
         electroweakino */
      /* treated by jets+MET grid as gsx */,
      sbgjsx
      /* squark to boson plus gluino to other squark to electroweakino */
      /* treated by jets+MET grid as gsx */,
      sjgbsx
      /* squark to gluino to boson plus other squark to electroweakino */
      /* treated by jets+MET grid as gsx */,
      sjgjsbsx
      /* squark to gluino to other squark to boson plus yet another squark to
         electroweakino */
      /* treated by jets+MET grid as gsx */,
      sizeOfEnumPlusOne /* possibly useful */
    };

    fullCascade( colorfulCascadeType const typeOfColorfulCascade,
                 int const firstDecayBodyNumber,
                 double const chargeConjugateSumFactor );
    virtual
    ~fullCascade();

    particlePointer
    getInitialScolored()
    const;
    particlePointer
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
    std::vector< particleWithInt* > const*
    getCascadeDefiner()
    const;
    colorfulCascadeType
    getColorfulCascadeType()
    const;
    virtual double
    getAcceptance( bool const scoloredIsNotAntiparticle,
                   acceptanceCutSet* const acceptanceCuts,
                   int const numberOfAdditionalJets,
                   int const numberOfNegativeElectrons,
                   int const numberOfPositiveElectrons,
                   int const numberOfNegativeMuons,
                   int const numberOfPositiveMuons )
    = 0;
    /* this calls the appropriate functions on ewinoCascade &, if not NULL,
     * bosonCascades, to build the required acceptance, taking into account
     * whether the charges should be swapped if scoloredIsNotAntiparticle is
     * false.
     */
    double
    getAcceptanceWithMinimumJets( bool const scoloredIsNotAntiparticle,
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
     * ewinoCascades->getOssfMinusOsdf( scoloredIsNotAntiparticle,
     *                                 acceptanceCuts ),
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

  protected:
    static int const numberOfSmFermionsFromElectroweakDecaysPerPartialCascade;
    static int const numberOfSmFermionsFromElectroweakDecaysPerFullCascade;

    inputHandler const* inputShortcut;
    particlePointer initialScolored;
    double beamEnergy;
    colorfulCascadeType const typeOfColorfulCascade;
    int const firstDecayBodyNumber;
    // this is the number of particles which are produced by the decay of
    // initialScolored.
    double const chargeConjugateSumFactor;
    /* this should be 1.0 for a cascade beginning with a squark & 2.0 for a
     * cascade beginning with a gluino, to cover decays through
     * charge-conjugate versions.
     */
    bool ewinoFlipsCharge;
    bool bosonFlipsCharge;
    fullCascade* subcascadePointer;
    minimalAllocationVector< particleWithInt > cascadeDefiner;
    particleWithInt* cascadeSegment;
    // this is just used for setting up cascadeDefiner properly.
    electroweakCascadeSet* ewinoCascades;
    electroweakCascadeSet* bosonCascades;
    std::list< int > soughtDecayProductList;
    /* by default, there is only 1 particle to be sought in initialScolored's
     * decays. the constructor sets up soughtDecayProductList to have just a
     * single entry, & setProperties just sets that entry, & nothing else makes
     * any changes to soughtDecayProductList. derived classes involving decays
     * to bosons add entries to their own instances.
     */
    double cachedBranchingRatio;
    bool branchingRatioNeedsToBeRecalculated;

    void
    buildOn( fullCascade* const copySource );
    /* this copies the basic stuff (inputShortcut, beamEnergy, cascadeDefiner,
     * ewinoCascade, bosonCascade) from copySource, though not
     * typeOfColorfulCascade (which is always set by constructors) or
     * initialScolored (which depends on the type of derived class), & also it
     * sets up soughtDecayProductList to look for copySource->initialScolored,
     * & adds copySource->initialScolored at the end of cascadeDefiner, with
     * firstDecayBodyNumber.
     */
    double
    getCombinedAcceptance( bool const bosonChargeNotFlipped,
                           bool const ewinoChargeNotFlipped,
                           acceptanceCutSet* const acceptanceCuts,
                           int const numberOfAdditionalJets,
                           int const numberOfNegativeElectrons,
                           int const numberOfPositiveElectrons,
                           int const numberOfNegativeMuons,
                           int const numberOfPositiveMuons );
    /* this combines the acceptancesPerCutSet from ewinoCascades &
     * bosonCascades. Whether the relevant scoloreds are particles or
     * antiparticles should be decided by the derived class using this
     * function.
     */
  };


  // this is derived class for direct squark to electroweakino fullCascades.
  class sxFullCascade : public fullCascade
  {
  public:
    sxFullCascade();
    virtual
    ~sxFullCascade();

    virtual void
    setProperties( inputHandler const* const inputShortcut,
                   particlePointer const initialSquark,
                   double const beamEnergy,
                   electroweakCascadeSet* const ewinoCascades );
    bool
    isOpen();
    // this returns true if the squark is heavy enough to decay into the
    // electroweakino, false otherwise.
    virtual double
    getBrToEwino( std::list< int > const* excludedSmParticles );
    /* this works out the branching ratio for the decays of colored sparticles
     * down to the 1st electroweakino in the cascade (including any decays to
     * bosons + squarks), using the given list of excluded SM particles to
     * exclude unwanted parts of the BR.
     */
    virtual double
    getAcceptance( bool const scoloredIsNotAntiparticle,
                   acceptanceCutSet* const acceptanceCuts,
                   int const numberOfAdditionalJets,
                   int const numberOfNegativeElectrons,
                   int const numberOfPositiveElectrons,
                   int const numberOfNegativeMuons,
                   int const numberOfPositiveMuons );
    /* this calls the appropriate functions on ewinoCascades &, if not NULL,
     * bosonCascades, to build the required acceptance, taking into account
     * whether the charges should be swapped if scoloredIsNotAntiparticle is
     * false.
     */

  //protected:
    //nothing.
  };


  // this is derived class for direct gluino to electroweakino fullCascades.
  class gxFullCascade : public fullCascade
  {
  public:
    gxFullCascade();
    virtual
    ~gxFullCascade();

    virtual void
    setProperties( inputHandler const* const inputShortcut,
                   double const beamEnergy,
                   electroweakCascadeSet* const ewinoCascades );
    bool
    isOpen();
    // this returns true if the gluino is heavy enough to decay into the
    // electroweakino, false otherwise.
    virtual double
    getBrToEwino( std::list< int > const* excludedSmParticles );
    /* this calls the appropriate functions on ewinoCascades &, if not NULL,
     * bosonCascades, to build the required acceptance, taking into account
     * whether the charges should be swapped if scoloredIsNotAntiparticle is
     * false.
     */
    virtual double
    getAcceptance( bool const scoloredIsNotAntiparticle,
                   acceptanceCutSet* const acceptanceCuts,
                   int const numberOfAdditionalJets,
                   int const numberOfNegativeElectrons,
                   int const numberOfPositiveElectrons,
                   int const numberOfNegativeMuons,
                   int const numberOfPositiveMuons );
    /* this calls the appropriate functions on ewinoCascades, summing over the
     * 2 charge versions if the electroweakino is a chargino. it ignores the
     * value of scoloredIsNotAntiparticle.
     */

  //protected:
    // nothing.
  };


  // this is derived class for squark to gluino to electroweakino fullCascades.
  class sjgxFullCascade : public fullCascade
  {
  public:
    sjgxFullCascade();
    virtual
    ~sjgxFullCascade();

    virtual void
    setProperties( particlePointer const initialSquark,
                   fullCascade* const gxCascade );
    virtual double
    getAcceptance( bool const scoloredIsNotAntiparticle,
                   acceptanceCutSet* const acceptanceCuts,
                   int const numberOfAdditionalJets,
                   int const numberOfNegativeElectrons,
                   int const numberOfPositiveElectrons,
                   int const numberOfNegativeMuons,
                   int const numberOfPositiveMuons );
    /* this calls the appropriate functions on ewinoCascades &, if not NULL,
     * bosonCascades, to build the required acceptance, taking into account
     * whether the charges should be swapped if scoloredIsNotAntiparticle is
     * false.
     */

  //protected:
    // nothing.
  };


  // this is derived class for gluino to squark to electroweakino fullCascades.
  class gjsxFullCascade : public fullCascade
  {
  public:
    gjsxFullCascade();
    virtual
    ~gjsxFullCascade();

    virtual void
    setProperties( fullCascade* const sxCascade );
    virtual double
    getAcceptance( bool const scoloredIsNotAntiparticle,
                   acceptanceCutSet* const acceptanceCuts,
                   int const numberOfAdditionalJets,
                   int const numberOfNegativeElectrons,
                   int const numberOfPositiveElectrons,
                   int const numberOfNegativeMuons,
                   int const numberOfPositiveMuons );
    /* this calls the appropriate functions on ewinoCascades &, if not NULL,
     * bosonCascades, to build the required acceptance, taking into account
     * whether the charges should be swapped if scoloredIsNotAntiparticle is
     * false.
     */

  //protected:
    // nothing.
  };


  // this is derived class for squark to gluino to squark to electroweakino
  // fullCascades, treated by jets+MET grid as gsx.
  class sjgjsxFullCascade : public fullCascade
  {
  public:
    sjgjsxFullCascade();
    virtual
    ~sjgjsxFullCascade();

    virtual void
    setProperties( particlePointer const initialSquark,
                   fullCascade* const gjsxCascade );
    virtual double
    getAcceptance( bool const scoloredIsNotAntiparticle,
                   acceptanceCutSet* const acceptanceCuts,
                   int const numberOfAdditionalJets,
                   int const numberOfNegativeElectrons,
                   int const numberOfPositiveElectrons,
                   int const numberOfNegativeMuons,
                   int const numberOfPositiveMuons );
    /* this calls the appropriate functions on ewinoCascades &, if not NULL,
     * bosonCascades, to build the required acceptance, taking into account
     * whether the charges should be swapped if scoloredIsNotAntiparticle is
     * false.
     */

  //protected:
    // nothing.
  };


  // this is derived class for squark to boson plus other squark to
  // electroweakino fullCascades, treated by jets+MET grid as sx.
  class sbsxFullCascade : public fullCascade
  {
  public:
    sbsxFullCascade();
    virtual
    ~sbsxFullCascade();

    virtual void
    setProperties( particlePointer const initialSquark,
                   electroweakCascadeSet* const bosonCascades,
                   fullCascade* const sxCascade );
    virtual double
    getAcceptance( bool const scoloredIsNotAntiparticle,
                   acceptanceCutSet* const acceptanceCuts,
                   int const numberOfAdditionalJets,
                   int const numberOfNegativeElectrons,
                   int const numberOfPositiveElectrons,
                   int const numberOfNegativeMuons,
                   int const numberOfPositiveMuons );
    /* this calls the appropriate functions on ewinoCascades &, if not NULL,
     * bosonCascades, to build the required acceptance, taking into account
     * whether the charges should be swapped if scoloredIsNotAntiparticle is
     * false.
     */

  //protected:
    // nothing.
  };


  // this is derived class for gluino to boson plus squark to electroweakino
  // fullCascades, treated by jets+MET grid as sx.
  class gbsxFullCascade : public fullCascade
  {
  public:
    gbsxFullCascade();
    virtual
    ~gbsxFullCascade();

    virtual void
    setProperties( electroweakCascadeSet* const bosonCascades,
                   fullCascade* const sxCascade );
    virtual double
    getAcceptance( bool const scoloredIsNotAntiparticle,
                   acceptanceCutSet* const acceptanceCuts,
                   int const numberOfAdditionalJets,
                   int const numberOfNegativeElectrons,
                   int const numberOfPositiveElectrons,
                   int const numberOfNegativeMuons,
                   int const numberOfPositiveMuons );
    /* this calls the appropriate functions on ewinoCascades &, if not NULL,
     * bosonCascades, to build the required acceptance, taking into account
     * whether the charges should be swapped if scoloredIsNotAntiparticle is
     * false.
     */

  //protected:
    // nothing.
  };


  // this is derived class for gluino to squark to boson plus other squark to
  // electroweakino fullCascades, treated by jets+MET grid as gsx.
  class gjsbsxFullCascade : public fullCascade
  {
  public:
    gjsbsxFullCascade();
    virtual
    ~gjsbsxFullCascade();

    virtual void
    setProperties( fullCascade* const sbsxCascade );
    virtual double
    getAcceptance( bool const scoloredIsNotAntiparticle,
                   acceptanceCutSet* const acceptanceCuts,
                   int const numberOfAdditionalJets,
                   int const numberOfNegativeElectrons,
                   int const numberOfPositiveElectrons,
                   int const numberOfNegativeMuons,
                   int const numberOfPositiveMuons );
    /* the assumptions made by this code mean that the branching ratios to
     * electroweakinos multiplied by acceptances through squark & antisquark
     * factorize. if the electroweakino is a chargino, it is assumed that only
     * either the squark or the antisquark has a non-zero BR to the chargino,
     * hence the compound BR from the gluino to the positively-charged chargino
     * is given by getBrToEwino, which will be the BR through the appropriate
     * charge of squark.
     */

  //protected:
    // nothing.
  };


  // this is derived class for squark to boson plus gluino to electroweakino
  // fullCascades, treated by jets+MET grid as sgx.
  class sbgxFullCascade : public fullCascade
  {
  public:
    sbgxFullCascade();
    virtual
    ~sbgxFullCascade();

    virtual void
    setProperties( particlePointer const initialSquark,
                   electroweakCascadeSet* const bosonCascades,
                   fullCascade* const gxCascade );
    virtual double
    getAcceptance( bool const scoloredIsNotAntiparticle,
                   acceptanceCutSet* const acceptanceCuts,
                   int const numberOfAdditionalJets,
                   int const numberOfNegativeElectrons,
                   int const numberOfPositiveElectrons,
                   int const numberOfNegativeMuons,
                   int const numberOfPositiveMuons );
    /* this calls the appropriate functions on ewinoCascades &, if not NULL,
     * bosonCascades, to build the required acceptance, taking into account
     * whether the charges should be swapped if scoloredIsNotAntiparticle is
     * false.
     */

  //protected:
    // nothing.
  };


  // this is derived class for squark to boson plus other squark to gluino to
  // electroweakino fullCascades, treated by jets+MET grid as sgx.
  class sbsjgxFullCascade : public fullCascade
  {
  public:
    sbsjgxFullCascade();
    virtual
    ~sbsjgxFullCascade();

    virtual void
    setProperties( particlePointer const initialSquark,
                   electroweakCascadeSet* const bosonCascades,
                   fullCascade* const sjgxCascade );
    virtual double
    getAcceptance( bool const scoloredIsNotAntiparticle,
                   acceptanceCutSet* const acceptanceCuts,
                   int const numberOfAdditionalJets,
                   int const numberOfNegativeElectrons,
                   int const numberOfPositiveElectrons,
                   int const numberOfNegativeMuons,
                   int const numberOfPositiveMuons );
    /* this calls the appropriate functions on ewinoCascades &, if not NULL,
     * bosonCascades, to build the required acceptance, taking into account
     * whether the charges should be swapped if scoloredIsNotAntiparticle is
     * false.
     */

  //protected:
    // nothing.
  };


  /* this is derived class for squark to boson plus other squark to gluino to
   * yet another squark to electroweakino fullCascades, treated by jets+MET
   * grid as gsx.
   */
  class sbsjgjsxFullCascade : public fullCascade
  {
  public:
    sbsjgjsxFullCascade();
    virtual
    ~sbsjgjsxFullCascade();

    virtual void
    setProperties( particlePointer const initialSquark,
                   electroweakCascadeSet* const bosonCascades,
                   fullCascade* const sjgjsxCascade );
    virtual double
    getAcceptance( bool const scoloredIsNotAntiparticle,
                   acceptanceCutSet* const acceptanceCuts,
                   int const numberOfAdditionalJets,
                   int const numberOfNegativeElectrons,
                   int const numberOfPositiveElectrons,
                   int const numberOfNegativeMuons,
                   int const numberOfPositiveMuons );
    /* this calls the appropriate functions on ewinoCascades &, if not NULL,
     * bosonCascades, to build the required acceptance, taking into account
     * whether the charges should be swapped if scoloredIsNotAntiparticle is
     * false.
     */

  //protected:
    // nothing.
  };


  // this is derived class for squark to boson plus gluino to other squark to
  // electroweakino fullCascades, treated by jets+MET grid as gsx.
  class sbgjsxFullCascade : public fullCascade
  {
  public:
    sbgjsxFullCascade();
    virtual
    ~sbgjsxFullCascade();

    virtual void
    setProperties( particlePointer const initialSquark,
                   electroweakCascadeSet* const bosonCascades,
                   fullCascade* const gjsxCascade );
    virtual double
    getAcceptance( bool const scoloredIsNotAntiparticle,
                   acceptanceCutSet* const acceptanceCuts,
                   int const numberOfAdditionalJets,
                   int const numberOfNegativeElectrons,
                   int const numberOfPositiveElectrons,
                   int const numberOfNegativeMuons,
                   int const numberOfPositiveMuons );
    /* this calls the appropriate functions on ewinoCascades &, if not NULL,
     * bosonCascades, to build the required acceptance, taking into account
     * whether the charges should be swapped if scoloredIsNotAntiparticle is
     * false.
     */

  //protected:
    // nothing.
  };


  // this is derived class for squark to gluino to boson plus other squark to
  // electroweakino fullCascades, treated by jets+MET grid as gsx.
  class sjgbsxFullCascade : public fullCascade
  {
  public:
    sjgbsxFullCascade();
    virtual
    ~sjgbsxFullCascade();

    virtual void
    setProperties( particlePointer const initialSquark,
                   fullCascade* const gbsxCascade );
    virtual double
    getAcceptance( bool const scoloredIsNotAntiparticle,
                   acceptanceCutSet* const acceptanceCuts,
                   int const numberOfAdditionalJets,
                   int const numberOfNegativeElectrons,
                   int const numberOfPositiveElectrons,
                   int const numberOfNegativeMuons,
                   int const numberOfPositiveMuons );
    /* this calls the appropriate functions on ewinoCascades &, if not NULL,
     * bosonCascades, to build the required acceptance, taking into account
     * whether the charges should be swapped if scoloredIsNotAntiparticle is
     * false.
     */

  //protected:
    // nothing.
  };


  /* this is derived class for squark to gluino to other squark to boson plus
   * yet another squark to electroweakino fullCascades, treated by jets+MET
   * grid as gsx.
   */
  class sjgjsbsxFullCascade : public fullCascade
  {
  public:
    sjgjsbsxFullCascade();
    virtual
    ~sjgjsbsxFullCascade();

    virtual void
    setProperties( particlePointer const initialSquark,
                   fullCascade* const gjsvsxCascade );
    virtual double
    getAcceptance( bool const scoloredIsNotAntiparticle,
                   acceptanceCutSet* const acceptanceCuts,
                   int const numberOfAdditionalJets,
                   int const numberOfNegativeElectrons,
                   int const numberOfPositiveElectrons,
                   int const numberOfNegativeMuons,
                   int const numberOfPositiveMuons );
    /* this calls the appropriate functions on ewinoCascades &, if not NULL,
     * bosonCascades, to build the required acceptance, taking into account
     * whether the charges should be swapped if scoloredIsNotAntiparticle is
     * false.
     */

  //protected:
    // nothing.
  };


  /* this class holds the possible cascades for a given scolored, as well as a
   * pointer to a std::list of the other fullCascadeSet instances, which get
   * ordered by mass of initial colored sparticle by whichever fullCascadeSet
   * is called 1st.
   */
  class fullCascadeSet : public getsReadiedForNewPoint
  {
  public:
    typedef std::pair< std::list< fullCascadeSet* >*,
                       publicGetsReadiedForNewPoint* > readiableCascadeSetList;

    static bool
    massOrdered( fullCascadeSet* firstSet,
                 fullCascadeSet* secondSet );

    fullCascadeSet( inputHandler const* const inputShortcut,
                    particlePointer const initialScolored,
                electroweakCascadesForOneBeamEnergy* const electroweakCascades,
                    readiableCascadeSetList* const squarkCascadeSetList,
                    double const beamEnergy );
    virtual
    ~fullCascadeSet();

    std::vector< fullCascade* >*
    getOpenCascades();
    // this calls setUpCascades() if it needs to be readied for this point.
    std::vector< sxFullCascade* > const*
    getSxCascades();
    std::vector< gxFullCascade* > const*
    getGxCascades();
    std::vector< sjgxFullCascade* > const*
    getSjgxCascades();
    std::vector< gjsxFullCascade* > const*
    getGjsxCascades();
    std::vector< sjgjsxFullCascade* > const*
    getSjgjsxCascades();
    std::vector< sbsxFullCascade* > const*
    getSbsxCascades();
    std::vector< gbsxFullCascade* > const*
    getGbsxCascades();
    std::vector< gjsbsxFullCascade* > const*
    getGjsbsxCascades();
    std::vector< sbgxFullCascade* > const*
    getSbgxCascades();
    std::vector< sbsjgxFullCascade* > const*
    getSbsjgxCascades();
    std::vector< sbsjgjsxFullCascade* > const*
    getSbsjgjsxCascades();
    std::vector< sbgjsxFullCascade* > const*
    getSbgjsxCascades();
    std::vector< sjgbsxFullCascade* > const*
    getSjgbsxCascades();
    std::vector< sjgjsbsxFullCascade* > const*
    getSjgjsbsxCascades();
    particlePointer
    getInitialScolored()
    const;

  protected:
    inputHandler const* const inputShortcut;
    particlePointer const initialScolored;
    electroweakCascadesForOneBeamEnergy* const electroweakCascadeSource;
    std::vector< fullCascade* > openCascades;
    // this holds pointers to all the open cascades.

    // 1 of the sxCascades & gxCascades will be empty, the other has all the
    // direct decays of initialScolored to electroweakinos:
    minimalAllocationVector< sxFullCascade > sxCascades;
    minimalAllocationVector< gxFullCascade > gxCascades;
    std::vector< sxFullCascade* > openSxCascades;
    std::vector< gxFullCascade* > openGxCascades;

    /* all of these minimalAllocationVectors hold the cascades built upon
     * cascades from the other fullCascadeSets, & the bools note if they have
     * been worked out for this point or not:
     */
    minimalAllocationVector< sjgxFullCascade > sjgxCascades;
    minimalAllocationVector< gjsxFullCascade > gjsxCascades;
    minimalAllocationVector< sjgjsxFullCascade > sjgjsxCascades;
    minimalAllocationVector< sbsxFullCascade > sbsxCascades;
    minimalAllocationVector< gbsxFullCascade > gbsxCascades;
    minimalAllocationVector< gjsbsxFullCascade > gjsbsxCascades;
    minimalAllocationVector< sbgxFullCascade > sbgxCascades;
    minimalAllocationVector< sbsjgxFullCascade > sbsjgxCascades;
    minimalAllocationVector< sbsjgjsxFullCascade > sbsjgjsxCascades;
    minimalAllocationVector< sbgjsxFullCascade > sbgjsxCascades;
    minimalAllocationVector< sjgbsxFullCascade > sjgbsxCascades;
    minimalAllocationVector< sjgjsbsxFullCascade > sjgjsbsxCascades;

    readiableCascadeSetList* const squarkCascadeSetList;
    double const beamEnergy;
    std::list< int > soughtPositivePdgCodeList;
    std::list< int > soughtNegativePdgCodeList;

    void
    ensureCascadesAreUpToDate();
    // this calls setUpCascades() if needsToPrepareForThisPoint() is true.
    void
    setUpCascades();
    // this clears openCascades, then sets it to be filled with all open
    // fullCascades that initialScolored has for this point.
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


  // this class specializes fullCascade for the case of cascades beginning with
  // a squark.
  class squarkFullCascadeSet : public fullCascadeSet
  {
  public:
    squarkFullCascadeSet( inputHandler const* const inputShortcut,
                   electroweakCascadesForOneBeamEnergy* const ewCascadeHandler,
                          particlePointer const initialScolored,
                          readiableCascadeSetList* const squarkCascadeSetList,
                          fullCascadeSet* const gluinoCascadeSet,
                          double const beamEnergy );
    virtual
    ~squarkFullCascadeSet();

  protected:
    fullCascadeSet* const gluinoCascadeSet;

    virtual void
    buildLongerCascades();
    // this should set up all the cascades that involve taking subcascades from
    // the other fullCascadeSets.
  };


  // this class specializes fullCascade for the case of cascades beginning with
  // a gluino.
  class gluinoFullCascadeSet : public fullCascadeSet
  {
  public:
    gluinoFullCascadeSet( inputHandler const* const inputShortcut,
                   electroweakCascadesForOneBeamEnergy* const ewCascadeHandler,
                          readiableCascadeSetList* const squarkCascadeSetList,
                          double const beamEnergy );
    virtual
    ~gluinoFullCascadeSet();

  protected:
    virtual void
    buildLongerCascades();
    // this should set up all the cascades that involve taking subcascades from
    // the other fullCascadeSets.
  };


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
    std::vector< squarkFullCascadeSet* > squarkCascadeSets;
    gluinoFullCascadeSet* gluinoCascadeSet;
    fullCascadeSet::readiableCascadeSetList squarkCascadeSetList;
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
  fullCascade::getInitialScolored()
  const
  {
    return initialScolored;
  }

  inline particlePointer
  fullCascade::getElectroweakinoAtEndOfScoloredness()
  const
  {
    return ewinoCascades->getElectroweakDecayer();
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
  fullCascade::getBrToEwino( std::list< int > const* excludedSmParticles )
  /* this works out the branching ratio for the decays of colored sparticles
   * down to the 1st electroweakino in the cascade (including any decays to
   * bosons + squarks), using the given list of excluded SM particles to
   * exclude unwanted parts of the BR. by default, it multiplies the BR for the
   * initial decay by subcascadePointer->getBrToEwino.
   */
  {
    /* if initialScolored is a gluino & the electroweakino is a chargino, we
     * return the sum of BRs to both charges of chargino, provided that
     * chargeConjugateSumFactor has been set correctly.
     */
    if( branchingRatioNeedsToBeRecalculated )
    {
      cachedBranchingRatio =
      ( initialScolored->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                                       &soughtDecayProductList,
                                                          excludedSmParticles )
        * chargeConjugateSumFactor
        * subcascadePointer->getBrToEwino( excludedSmParticles ) );
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
                                          bool const scoloredIsNotAntiparticle,
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
      returnDouble += getAcceptance( scoloredIsNotAntiparticle,
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
  /* this copies the basic stuff (inputShortcut, beamEnergy, cascadeDefiner,
   * ewinoCascade, vectorCascade) from copySource, though not
   * typeOfColorfulCascade (which is always set by constructors) or
   * initialScolored (which depends on the type of derived class), & also it
   * sets up soughtDecayProductList to look for copySource->initialScolored,
   * & adds copySource->initialScolored at the end of cascadeDefiner, with
   * firstDecayBodyNumber.
   */
  {
    inputShortcut = copySource->inputShortcut;
    beamEnergy = copySource->beamEnergy;
    ewinoCascades = copySource->ewinoCascades;
    bosonCascades = copySource->bosonCascades;
    subcascadePointer = copySource;
    soughtDecayProductList.front()
    = copySource->initialScolored->get_PDG_code();
    cascadeDefiner.clearEntries();
    for( int cascadeFiller( 0 );
         copySource->cascadeDefiner.getSize() > cascadeFiller;
         ++cascadeFiller )
    {
      cascadeSegment = cascadeDefiner.addNewAtEnd();
      cascadeSegment->first
      = copySource->cascadeDefiner.getPointer( cascadeFiller )->first;
      cascadeSegment->second
      = copySource->cascadeDefiner.getPointer( cascadeFiller )->second;
    }
    cascadeSegment = cascadeDefiner.addNewAtEnd();
    cascadeSegment->first = copySource->initialScolored;
    cascadeSegment->second = firstDecayBodyNumber;
    ewinoFlipsCharge = copySource->ewinoFlipsCharge;
  }



  inline bool
  sxFullCascade::isOpen()
  // this returns true if the squark is heavy enough to decay into the
  // electroweakino, false otherwise.
  {
    if( initialScolored->get_absolute_mass()
        > ewinoCascades->getElectroweakDecayer()->get_absolute_mass() )
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  inline double
  sxFullCascade::getBrToEwino( std::list< int > const* excludedSmParticles )
  /* this works out the branching ratio for the decays of colored sparticles
   * down to the 1st electroweakino in the cascade (including any decays to
   * bosons + squarks), using the given list of excluded SM particles to
   * exclude unwanted parts of the BR.
   */
  {
    if( branchingRatioNeedsToBeRecalculated )
    {
      cachedBranchingRatio
      = initialScolored->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                                       &soughtDecayProductList,
                                                         excludedSmParticles );
      branchingRatioNeedsToBeRecalculated = false;
    }
    return cachedBranchingRatio;
  }

  inline double
  sxFullCascade::getAcceptance( bool const scoloredIsNotAntiparticle,
                                acceptanceCutSet* const acceptanceCuts,
                                int const numberOfAdditionalJets,
                                int const numberOfNegativeElectrons,
                                int const numberOfPositiveElectrons,
                                int const numberOfNegativeMuons,
                                int const numberOfPositiveMuons )
  /* this calls the appropriate functions on ewinoCascades &, if not NULL,
   * bosonCascades, to build the required acceptance, taking into account
   * whether the charges should be swapped if scoloredIsNotAntiparticle is
   * false.
   */
  {
    if( ( scoloredIsNotAntiparticle
          &&
          !ewinoFlipsCharge )
        ||
        ( !scoloredIsNotAntiparticle
          &&
          ewinoFlipsCharge ) )
      // if the electroweakino is a particle, we return the acceptance for the
      // particle.
    {
      return ewinoCascades->getAcceptance( acceptanceCuts,
                                           numberOfAdditionalJets,
                                           numberOfNegativeElectrons,
                                           numberOfPositiveElectrons,
                                           numberOfNegativeMuons,
                                           numberOfPositiveMuons );
    }
    else
      // if the electroweakino an antiparticle, we swap the charges.
    {
      return ewinoCascades->getAcceptance( acceptanceCuts,
                                           numberOfAdditionalJets,
                                           numberOfPositiveElectrons,
                                           numberOfNegativeElectrons,
                                           numberOfPositiveMuons,
                                           numberOfNegativeMuons );
    }
  }



  inline void
  gxFullCascade::setProperties( inputHandler const* const inputShortcut,
                                double const beamEnergy,
                                electroweakCascadeSet* const ewinoCascade )
  {
    this->inputShortcut = inputShortcut;
    this->initialScolored = inputShortcut->getGluino();
    this->beamEnergy = beamEnergy;
    this->ewinoCascades = ewinoCascade;
    soughtDecayProductList.front() =
                         ewinoCascade->getElectroweakDecayer()->get_PDG_code();
    resetCachedBranchingRatio();
    cascadeDefiner.clearEntries();
    // reset cascadeDefiner.
    cascadeSegment = cascadeDefiner.addNewAtEnd();
    cascadeSegment->first = ewinoCascade->getElectroweakDecayer();
    // gx means only 1 decay to be recorded.
    cascadeSegment->second = firstDecayBodyNumber;
    // gx also means that the decay is 2-body.
  }

  inline bool
  gxFullCascade::isOpen()
  // this returns true if the gluino is heavy enough to decay into the
  // electroweakino, false otherwise.
  {
    if( inputShortcut->getGluinoMass()
        > ewinoCascades->getElectroweakDecayer()->get_absolute_mass() )
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  inline double
  gxFullCascade::getBrToEwino( std::list< int > const* excludedSmParticles )
  /* this works out the branching ratio for the decays of colored sparticles
   * down to the 1st electroweakino in the cascade (including any decays to
   * vector bosons + squarks), using the given list of excluded SM particles
   * to exclude unwanted parts of the BR.
   */
  {
    if( branchingRatioNeedsToBeRecalculated )
    {
      if( ewinoCascades->getElectroweakDecayer()->counts_as_self_conjugate() )
        // if the gluino has only 1 charge version of this decay...
      {
        cachedBranchingRatio
        = initialScolored->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                                       &soughtDecayProductList,
                                                         excludedSmParticles );
      }
      else
      {
        cachedBranchingRatio
        = ( 2.0 * initialScolored->inspect_direct_decay_handler(
                                             )->get_branching_ratio_for_subset(
                                                       &soughtDecayProductList,
                                                       excludedSmParticles ) );
      }
      branchingRatioNeedsToBeRecalculated = false;
    }
    return cachedBranchingRatio;
  }



  inline void
  sjgxFullCascade::setProperties( particlePointer const initialSquark,
                                  fullCascade* const gxCascade )
  {
    initialScolored = initialSquark;
    buildOn( gxCascade );
    resetCachedBranchingRatio();
  }

  inline double
  sjgxFullCascade::getAcceptance( bool const scoloredIsNotAntiparticle,
                                  acceptanceCutSet* const acceptanceCuts,
                                  int const numberOfAdditionalJets,
                                  int const numberOfNegativeElectrons,
                                  int const numberOfPositiveElectrons,
                                  int const numberOfNegativeMuons,
                                  int const numberOfPositiveMuons )
  /* this calls the appropriate functions on ewinoCascades &, if not NULL,
   * bosonCascades, to build the required acceptance, taking into account
   * whether the charges should be swapped if scoloredIsNotAntiparticle is
   * false.
   */
  {
    return subcascadePointer->getAcceptance( scoloredIsNotAntiparticle,
                                             acceptanceCuts,
                                             numberOfAdditionalJets,
                                             numberOfNegativeElectrons,
                                             numberOfPositiveElectrons,
                                             numberOfNegativeMuons,
                                             numberOfPositiveMuons );
  }



  inline void
  gjsxFullCascade::setProperties( fullCascade* const sxCascade )
  {
    buildOn( sxCascade );
    initialScolored = inputShortcut->getGluino();
    resetCachedBranchingRatio();
  }

  inline double
  gjsxFullCascade::getAcceptance( bool const scoloredIsNotAntiparticle,
                                  acceptanceCutSet* const acceptanceCuts,
                                  int const numberOfAdditionalJets,
                                  int const numberOfNegativeElectrons,
                                  int const numberOfPositiveElectrons,
                                  int const numberOfNegativeMuons,
                                  int const numberOfPositiveMuons )
  /* the assumptions made by this code mean that the branching ratios to
   * electroweakinos multiplied by acceptances through squark & antisquark
   * factorize. if the electroweakino is a chargino, it is assumed that only
   * either the squark or the antisquark has a non-zero BR to the chargino,
   * hence the compound BR from the gluino to the positively-charged chargino
   * is given by getBrToEwino, which will be the BR through the appropriate
   * charge of squark.
   */
  {
    return ( 0.5 * ( subcascadePointer->getAcceptance( true,
                                                       acceptanceCuts,
                                                       numberOfAdditionalJets,
                                                     numberOfNegativeElectrons,
                                                     numberOfPositiveElectrons,
                                                       numberOfNegativeMuons,
                                                       numberOfPositiveMuons )
                     + subcascadePointer->getAcceptance( false,
                                                         acceptanceCuts,
                                                        numberOfAdditionalJets,
                                                     numberOfNegativeElectrons,
                                                     numberOfPositiveElectrons,
                                                         numberOfNegativeMuons,
                                                   numberOfPositiveMuons ) ) );
    // return the weighted acceptances of the squark & the antisquark versions.
  }



  inline void
  sjgjsxFullCascade::setProperties( particlePointer const initialSquark,
                                    fullCascade* const gjsxCascade )
  {
    initialScolored = initialSquark;
    buildOn( gjsxCascade );
    resetCachedBranchingRatio();
  }

  inline double
  sjgjsxFullCascade::getAcceptance( bool const scoloredIsNotAntiparticle,
                                    acceptanceCutSet* const acceptanceCuts,
                                    int const numberOfAdditionalJets,
                                    int const numberOfNegativeElectrons,
                                    int const numberOfPositiveElectrons,
                                    int const numberOfNegativeMuons,
                                    int const numberOfPositiveMuons )
  {
    return subcascadePointer->getAcceptance( scoloredIsNotAntiparticle,
                                             acceptanceCuts,
                                             numberOfAdditionalJets,
                                             numberOfNegativeElectrons,
                                             numberOfPositiveElectrons,
                                             numberOfNegativeMuons,
                                             numberOfPositiveMuons );
    // subascade takes care of counting over different charginos if needed.
  }



  inline void
  sbsxFullCascade::setProperties( particlePointer const initialSquark,
                                  electroweakCascadeSet* const bosonCascades,
                                  fullCascade* const sxCascade )
  {
    initialScolored = initialSquark;
    buildOn( sxCascade );
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
  }

  inline double
  sbsxFullCascade::getAcceptance( bool const scoloredIsNotAntiparticle,
                                  acceptanceCutSet* const acceptanceCuts,
                                  int const numberOfAdditionalJets,
                                  int const numberOfNegativeElectrons,
                                  int const numberOfPositiveElectrons,
                                  int const numberOfNegativeMuons,
                                  int const numberOfPositiveMuons )
  /* this calls the appropriate functions on ewinoCascades &, if not NULL,
   * bosonCascades, to build the required acceptance, taking into account
   * whether the charges should be swapped if scoloredIsNotAntiparticle is
   * false.
   */
  {
    // perversely, we assume an antiparticle 1st. also, squark decays via
    // bosons do not produce antisquarks.
    bool bosonChargeNotFlipped( bosonFlipsCharge );
    bool ewinoChargeNotFlipped( ewinoFlipsCharge );
    if( scoloredIsNotAntiparticle )
    {
      bosonChargeNotFlipped = !bosonFlipsCharge;
      ewinoChargeNotFlipped = !ewinoFlipsCharge;
    }
    return getCombinedAcceptance( bosonChargeNotFlipped,
                                  ewinoChargeNotFlipped,
                                  acceptanceCuts,
                                  numberOfAdditionalJets,
                                  numberOfNegativeElectrons,
                                  numberOfPositiveElectrons,
                                  numberOfNegativeMuons,
                                  numberOfPositiveMuons );
    // in the used conventions, an incoming squark decays to a vector boson
    // plus a squark, hence the same bool is used for both.
  }



  inline void
  gbsxFullCascade::setProperties( electroweakCascadeSet* const bosonCascades,
                                  fullCascade* const sxCascade )
  {
    buildOn( sxCascade );
    resetCachedBranchingRatio();
    initialScolored = inputShortcut->getGluino();
    this->bosonCascades = bosonCascades;
    if( ( CppSLHA::PDG_code::W_plus
          == bosonCascades->getElectroweakDecayer()->get_PDG_code() )
        &&
        ( inputShortcut->isIn( soughtDecayProductList.front(),
                               inputShortcut->getSupTypes() ) ) )
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
  }

  inline double
  gbsxFullCascade::getAcceptance( bool const scoloredIsNotAntiparticle,
                                  acceptanceCutSet* const acceptanceCuts,
                                  int const numberOfAdditionalJets,
                                  int const numberOfNegativeElectrons,
                                  int const numberOfPositiveElectrons,
                                  int const numberOfNegativeMuons,
                                  int const numberOfPositiveMuons )
  /* this calls the appropriate functions on ewinoCascades &, if not NULL,
   * bosonCascades, to build the required acceptance, taking into account
   * whether the charges should be swapped if scoloredIsNotAntiparticle is
   * false.
   */
  {
    return ( 0.5 * ( getCombinedAcceptance( !bosonFlipsCharge,
                                            !ewinoFlipsCharge,
                                            acceptanceCuts,
                                            numberOfAdditionalJets,
                                            numberOfNegativeElectrons,
                                            numberOfPositiveElectrons,
                                            numberOfNegativeMuons,
                                            numberOfPositiveMuons )
                    + getCombinedAcceptance( bosonFlipsCharge,
                                             ewinoFlipsCharge,
                                             acceptanceCuts,
                                             numberOfAdditionalJets,
                                             numberOfNegativeElectrons,
                                             numberOfPositiveElectrons,
                                             numberOfNegativeMuons,
                                             numberOfPositiveMuons ) ) );
    // getBrToEwino gives the sum of the decays via squark & antisquark, so the
    // acceptances factorize this way.
  }



  inline void
  gjsbsxFullCascade::setProperties( fullCascade* const sbsxCascade )
  {
    buildOn( sbsxCascade );
    resetCachedBranchingRatio();
    initialScolored = inputShortcut->getGluino();
  }

  inline double
  gjsbsxFullCascade::getAcceptance( bool const scoloredIsNotAntiparticle,
                                    acceptanceCutSet* const acceptanceCuts,
                                    int const numberOfAdditionalJets,
                                    int const numberOfNegativeElectrons,
                                    int const numberOfPositiveElectrons,
                                    int const numberOfNegativeMuons,
                                    int const numberOfPositiveMuons )
  /* the assumptions made by this code mean that the branching ratios to
   * electroweakinos multiplied by acceptances through squark & antisquark
   * factorize. if the electroweakino is a chargino, it is assumed that only
   * either the squark or the antisquark has a non-zero BR to the chargino,
   * hence the compound BR from the gluino to the positively-charged chargino
   * is given by getBrToEwino, which will be the BR through the appropriate
   * charge of squark.
   */
  {
    return ( 0.5 * ( subcascadePointer->getAcceptance( true,
                                                       acceptanceCuts,
                                                       numberOfAdditionalJets,
                                                     numberOfNegativeElectrons,
                                                     numberOfPositiveElectrons,
                                                       numberOfNegativeMuons,
                                                       numberOfPositiveMuons )
                     + subcascadePointer->getAcceptance( false,
                                                         acceptanceCuts,
                                                        numberOfAdditionalJets,
                                                     numberOfNegativeElectrons,
                                                     numberOfPositiveElectrons,
                                                         numberOfNegativeMuons,
                                                   numberOfPositiveMuons ) ) );
    // return the weighted acceptances of the squark & the antisquark versions.
  }



  inline void
  sbgxFullCascade::setProperties( particlePointer const initialSquark,
                                  electroweakCascadeSet* const bosonCascades,
                                  fullCascade* const gxCascade )
  {
    initialScolored = initialSquark;
    buildOn( gxCascade );
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
  }

  inline double
  sbgxFullCascade::getAcceptance( bool const scoloredIsNotAntiparticle,
                                  acceptanceCutSet* const acceptanceCuts,
                                  int const numberOfAdditionalJets,
                                  int const numberOfNegativeElectrons,
                                  int const numberOfPositiveElectrons,
                                  int const numberOfNegativeMuons,
                                  int const numberOfPositiveMuons )
  /* this calls the appropriate functions on ewinoCascades &, if not NULL,
   * bosonCascades, to build the required acceptance, taking into account
   * whether the charges should be swapped if scoloredIsNotAntiparticle is
   * false.
   */
  {
    // perversely, we assume an antiparticle 1st.
    bool bosonChargeNotFlipped( bosonFlipsCharge );
    if( scoloredIsNotAntiparticle )
    {
      bosonChargeNotFlipped = !bosonFlipsCharge;
    }
    return getCombinedAcceptance( bosonChargeNotFlipped,
                                  true /* gluinos don't care. */,
                                  acceptanceCuts,
                                  numberOfAdditionalJets,
                                  numberOfNegativeElectrons,
                                  numberOfPositiveElectrons,
                                  numberOfNegativeMuons,
                                  numberOfPositiveMuons );
  }



  inline void
  sbsjgxFullCascade::setProperties( particlePointer const initialSquark,
                                    electroweakCascadeSet* const bosonCascades,
                                    fullCascade* const sjgxCascade )
  {
    initialScolored = initialSquark;
    buildOn( sjgxCascade );
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
  }

  inline double
  sbsjgxFullCascade::getAcceptance( bool const scoloredIsNotAntiparticle,
                                    acceptanceCutSet* const acceptanceCuts,
                                    int const numberOfAdditionalJets,
                                    int const numberOfNegativeElectrons,
                                    int const numberOfPositiveElectrons,
                                    int const numberOfNegativeMuons,
                                    int const numberOfPositiveMuons )
  /* this calls the appropriate functions on ewinoCascades &, if not NULL,
   * bosonCascades, to build the required acceptance, taking into account
   * whether the charges should be swapped if scoloredIsNotAntiparticle is
   * false.
   */
  {
    // perversely, we assume an antiparticle 1st.
    bool bosonChargeNotFlipped( bosonFlipsCharge );
    if( scoloredIsNotAntiparticle )
    {
      bosonChargeNotFlipped = !bosonFlipsCharge;
    }
    return getCombinedAcceptance( bosonChargeNotFlipped,
                                  true /* gluinos don't care. */,
                                  acceptanceCuts,
                                  numberOfAdditionalJets,
                                  numberOfNegativeElectrons,
                                  numberOfPositiveElectrons,
                                  numberOfNegativeMuons,
                                  numberOfPositiveMuons );
  }



  inline void
  sbsjgjsxFullCascade::setProperties( particlePointer const initialSquark,
                                    electroweakCascadeSet* const vectorCascade,
                                      fullCascade* const sjgjsxCascade )
  {
    initialScolored = initialSquark;
    buildOn( sjgjsxCascade );
    resetCachedBranchingRatio();
    this->bosonCascades = vectorCascade;
    if( ( CppSLHA::PDG_code::W_plus
          == vectorCascade->getElectroweakDecayer()->get_PDG_code() )
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
  }

  inline double
  sbsjgjsxFullCascade::getAcceptance( bool const scoloredIsNotAntiparticle,
                                      acceptanceCutSet* const acceptanceCuts,
                                      int const numberOfAdditionalJets,
                                      int const numberOfNegativeElectrons,
                                      int const numberOfPositiveElectrons,
                                      int const numberOfNegativeMuons,
                                      int const numberOfPositiveMuons )
  /* this calls the appropriate functions on ewinoCascades &, if not NULL,
   * bosonCascades, to build the required acceptance, taking into account
   * whether the charges should be swapped if scoloredIsNotAntiparticle is
   * false.
   */
  {
    // perversely, we assume an antiparticle 1st.
    bool bosonChargeNotFlipped( bosonFlipsCharge );
    if( scoloredIsNotAntiparticle )
    {
      bosonChargeNotFlipped = !bosonFlipsCharge;
    }
    return getCombinedAcceptance( bosonChargeNotFlipped,
                                  true /* gluinos don't care. */,
                                  acceptanceCuts,
                                  numberOfAdditionalJets,
                                  numberOfNegativeElectrons,
                                  numberOfPositiveElectrons,
                                  numberOfNegativeMuons,
                                  numberOfPositiveMuons );
  }



  inline void
  sbgjsxFullCascade::setProperties( particlePointer const initialSquark,
                                    electroweakCascadeSet* const vectorCascade,
                                    fullCascade* const gjsxCascade )
  {
    initialScolored = initialSquark;
    buildOn( gjsxCascade );
    resetCachedBranchingRatio();
    this->bosonCascades = vectorCascade;
    if( ( CppSLHA::PDG_code::W_plus
          == vectorCascade->getElectroweakDecayer()->get_PDG_code() )
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
  }

  inline double
  sbgjsxFullCascade::getAcceptance( bool const scoloredIsNotAntiparticle,
                                    acceptanceCutSet* const acceptanceCuts,
                                    int const numberOfAdditionalJets,
                                    int const numberOfNegativeElectrons,
                                    int const numberOfPositiveElectrons,
                                    int const numberOfNegativeMuons,
                                    int const numberOfPositiveMuons )
  /* this calls the appropriate functions on ewinoCascades &, if not NULL,
   * bosonCascades, to build the required acceptance, taking into account
   * whether the charges should be swapped if scoloredIsNotAntiparticle is
   * false.
   */
  {
    // perversely, we assume an antiparticle 1st.
    bool bosonChargeNotFlipped( bosonFlipsCharge );
    if( scoloredIsNotAntiparticle )
    {
      bosonChargeNotFlipped = !bosonFlipsCharge;
    }
    return getCombinedAcceptance( bosonChargeNotFlipped,
                                  true /* gluinos don't care. */,
                                  acceptanceCuts,
                                  numberOfAdditionalJets,
                                  numberOfNegativeElectrons,
                                  numberOfPositiveElectrons,
                                  numberOfNegativeMuons,
                                  numberOfPositiveMuons );
  }



  inline void
  sjgbsxFullCascade::setProperties( particlePointer const initialSquark,
                                    fullCascade* const gbsxCascade )
  {
    initialScolored = initialSquark;
    buildOn( gbsxCascade );
    resetCachedBranchingRatio();
  }

  inline double
  sjgbsxFullCascade::getAcceptance( bool const scoloredIsNotAntiparticle,
                                    acceptanceCutSet* const acceptanceCuts,
                                    int const numberOfAdditionalJets,
                                    int const numberOfNegativeElectrons,
                                    int const numberOfPositiveElectrons,
                                    int const numberOfNegativeMuons,
                                    int const numberOfPositiveMuons )
  /* this calls the appropriate functions on ewinoCascades &, if not NULL,
   * bosonCascades, to build the required acceptance, taking into account
   * whether the charges should be swapped if scoloredIsNotAntiparticle is
   * false.
   */
  {
    return subcascadePointer->getAcceptance( scoloredIsNotAntiparticle,
                                             acceptanceCuts,
                                             numberOfAdditionalJets,
                                             numberOfNegativeElectrons,
                                             numberOfPositiveElectrons,
                                             numberOfNegativeMuons,
                                             numberOfPositiveMuons );
    // subcascadePointer takes care of counting over different charginos if
    // needed.
  }



  inline void
  sjgjsbsxFullCascade::setProperties( particlePointer const initialSquark,
                                      fullCascade* const gjsvsxCascade )
  {
    initialScolored = initialSquark;
    buildOn( gjsvsxCascade );
    resetCachedBranchingRatio();
  }

  inline double
  sjgjsbsxFullCascade::getAcceptance( bool const scoloredIsNotAntiparticle,
                                      acceptanceCutSet* const acceptanceCuts,
                                      int const numberOfAdditionalJets,
                                      int const numberOfNegativeElectrons,
                                      int const numberOfPositiveElectrons,
                                      int const numberOfNegativeMuons,
                                      int const numberOfPositiveMuons )
  /* this calls the appropriate functions on ewinoCascades &, if not NULL,
   * bosonCascades, to build the required acceptance, taking into account
   * whether the charges should be swapped if scoloredIsNotAntiparticle is
   * false.
   */
  {
    return subcascadePointer->getAcceptance( scoloredIsNotAntiparticle,
                                             acceptanceCuts,
                                             numberOfAdditionalJets,
                                             numberOfNegativeElectrons,
                                             numberOfPositiveElectrons,
                                             numberOfNegativeMuons,
                                             numberOfPositiveMuons );
    // subcascadePointer takes care of counting over different charginos if
    // needed.
  }



  inline bool
  fullCascadeSet::massOrdered( fullCascadeSet* firstSet,
                               fullCascadeSet* secondSet )
  // this returns true if firstSet's initialScolored has equal or lighter mass
  // than secondSet, false otherwise.
  {
    if( firstSet->initialScolored->get_absolute_mass()
        <= secondSet->initialScolored->get_absolute_mass() )
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
  fullCascadeSet::getInitialScolored()
  const
  {
    return initialScolored;
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
