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
#include "LHC-FASER_kinematics_stuff.hpp"
#include "LHC-FASER_derived_electroweak_cascade_stuff.hpp"

namespace LHC_FASER
{

  /* this class is for keeping all the stuff relevant to the cascade of the
   * EW_decayer of a full_cascade. it returns acceptances for individual
   * charges of leptons; N.B: the charges assume that the initial colored
   * sparticle is either a gluino or a positive-baryon-number (which happens to
   * correspond with positive PDG code) squark! further N.B.: the initial
   * colored sparticle is *not* known by the EW_decayer_values -
   * decaying_scolored is the *last* colored sparticle in the cascade decay.
   */
  // gone!



  /* this class keeps track of the acceptances for the cascade decay of a
   * colored sparticle directly to an on-shell electroweakino, which then
   * decays into 2 SM fermions plus the LSP (unless it already is the LSP) (the
   * 2 SM fermions could be tau leptons, & we actually allow each tau lepton
   * to decay into a visible product).
   *
   * each channel for the decay of the electroweakino caches its acceptance for
   * a given set of cuts.
   */
  // gone!



  /* this class is for keeping track of a particular cascade decay of a colored
   * sparticle down to the last electroweak decay, such that it keeps track of
   * all the sparticles of a cascade decay. it is an abstract base class.
   *
   * try to remember that there's 1 of these per beam energy! thus any given
   * scolored could have e.g. a set of full_cascades for 7 TeV & another set
   * for 14 TeV (with the same cascades because the hierarchy is the same, BUT
   * different scolored_directly_to_EWino_cascade &
   * scolored_to_squark_plus_vector_cascade instances because each looks up its
   * own effective squark mass, which depends on the beam energy (though really
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
   * the decay of a squark to another squark via a vector boson, should this
   * happen *earlier* in the cascade than the decay to & from the gluino), it
   * returns the sum of BRs which is twice the BR to a single charge, because
   * of the approximations used. getAcceptance should take care to avoid
   * double-counting for derived classes with cascades which start with
   * gluinos, insofar as if no asymmetric charge is requested (only applicable
   * to the leptons requested), it returns the acceptance for a single charge
   * (since the charge-conjugate has the same acceptance), but if a combination
   * with an asymmetric charge is requested, it returns half the acceptance, to
   * account for only half the BR returned actually corresponding to cascades
   * that can give the requested signal.
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
    enum colorfulCascadeType
    {
      /* we do not take into account multiple decays to vector bosons plus
       * other colored sparticles, because the cascades would be too long, with
       * too low branching ratios, as well as being quite difficult to get
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
      svsx /* squark to vector boson plus other squark to electroweakino */
      /* treated by jets+MET grid as sx */,
      gvsx /* gluino to vector boson plus squark to electroweakino */
      /* treated by jets+MET grid as sx */,
      gjsvsx
      /* gluino to squark to vector boson plus other squark to
         electroweakino */
      /* treated by jets+MET grid as gsx */,
      svgx /* squark to vector boson plus gluino to electroweakino */
      /* treated by jets+MET grid as sgx */,
      svsjgx
      /* squark to vector boson plus other squark to gluino to
         electroweakino */
      /* treated by jets+MET grid as sgx */,
      svsjgjsx
      /* squark to vector boson plus other squark to gluino to yet another
         squark to electroweakino */
      /* treated by jets+MET grid as gsx */,
      svgjsx
      /* squark to vector boson plus gluino to other squark to
         electroweakino */
      /* treated by jets+MET grid as gsx */,
      sjgvsx
      /* squark to gluino to vector boson plus other squark to
         electroweakino */
      /* treated by jets+MET grid as gsx */,
      sjgjsvsx
      /* squark to gluino to other squark to vector boson plus yet another
         squark to electroweakino */
      /* treated by jets+MET grid as gsx */
    };

    fullCascade( colorfulCascadeType const typeOfColorfulCascade )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~fullCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    static fullCascade*
    getNewFullCascade()
    // this just returns a new fullCascade using the default constructor.
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    std::pair< CppSLHA::particle_property_set const*,
                                              int >*
    fullCascade::getNewPair()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    input_handler const*
    getShortcut()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    getBrToEwino( std::list< int > const* excludedSmParticles )
    /* this works out the branching ratio for the decays of colored sparticles
     * down to the 1st electroweakino in the cascade (including any decays to
     * vector bosons + squarks), using the given list of excluded SM particles
     * to exclude unwanted parts of the BR.
     */
    = 0;

    minimalAllocationVector< std::pair< CppSLHA::particle_property_set const*,
                                        int > >*
    getCascadeDefiner()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual double
    getAcceptance( bool const scoloredIsNotAntiparticle,
                   acceptanceCutSet* const cuts,
                   int const numberOfAdditionalJets,
                   int const numberOfNegativeElectrons,
                   int const numberOfPositiveElectrons,
                   int const numberOfNegativeMuons,
                   int const numberOfPositiveMuons )
    /* this calls the appropriate functions on ewinoCascade &, if not NULL,
     * vectorCascade, to build the required acceptance, taking into account
     * whether the charges should be swapped if scoloredIsNotAntiparticle is
     * false.
     */
    = 0;

    double
    getAcceptanceWithMinimumJets( bool const scoloredIsNotAntiparticle,
                                  acceptanceCutSet* const cuts,
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
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    specifiedJetsOneOssfMinusOsdfPair( acceptanceCutSet* const cuts,
                                       int const numberOfAdditionalJets )
    /* this calls
     * ewinoCascade->getOssfMinusOsdf( scoloredIsNotAntiparticle,
     *                                 cuts ),
     * & if vectorCascade is not NULL, it calls
     * ewinoCascade->getAcceptance( [ { 0, 1, 2 } jets + 0 leptons for cuts ] )
     * & the same for vectorCascade, as well as
     * ->getOssfMinusOsdf( cuts ).
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    unspecifiedJetsSpecifiedOssfMinusOsdfPairs( acceptanceCutSet* const cuts,
                                                int const numberOfLeptonPairs )
    // this only calls ewinoCascade->getOssfMinusOsdf( cuts ),
    // & if vectorCascade is not NULL, vectorCascade->getOssfMinusOsdf( cuts ).
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    // the following just call various combinations of getAcceptance(...):

    double
    specifiedJetsSpecifiedChargeSummedLeptons( acceptanceCutSet* const cuts,
                                              int const numberOfAdditionalJets,
                                               int const numberOfLeptons )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    unspecifiedJetsSpecifiedChargeSummedLeptons( acceptanceCutSet* const cuts,
                                                 int const numberOfLeptons )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    static int const
    numberOfSmFermionsFromElectroweakDecaysPerPartialCascade = 2;
    static int const numberOfSmFermionsFromElectroweakDecaysPerFullCascade = 4;

    input_handler const* shortcut;
    CppSLHA::particle_property_set const* initialScolored;
    double beamEnergy;
    colorfulCascadeType const typeOfColorfulCascade;
    minimalAllocationVector< std::pair< CppSLHA::particle_property_set const*,
                                        int > > cascadeDefiner;
    std::pair< CppSLHA::particle_property_set const*,
               int >* cascadeSegment;
    // this is just used for setting up cascadeDefiner properly.
    electroweakCascade* ewinoCascade;
    electroweakCascade* vectorCascade;
    std::list< int > soughtDecayProductList;
    /* by default, there is only 1 particle to be sought in initialScolored's
     * decays. the constructor sets up soughtDecayProductList to have just a
     * single entry, & setProperties just sets that entry, & nothing else makes
     * any changes to soughtDecayProductList. derived classes involving decays
     * to vector bosons add entries to their own instances.
     */
    //double innerReturnDouble;
    //double outerReturnDouble;
    // these are just used for the various functions that add up various
    // combinations of getAcceptance.
    // actually, that's just asking for bugs, & it's almost certainly not
    // faster than declaring the doubles in functions when they're needed.

    void
    copyCascadeDefiner(
     minimalAllocationVector< std::pair< CppSLHA::particle_property_set const*,
                                         int > > cascadeDefiner )
    // this copies the values in the pairs of the given vector into
    // this->cascadeDefiner.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    getCombinedAcceptance( bool const vectorScoloredIsNotAntiparticle,
                           bool const ewinoScoloredIsNotAntiparticle,
                           acceptanceCutSet* const cuts,
                           int const numberOfAdditionalJets,
                           int const numberOfNegativeElectrons,
                           int const numberOfPositiveElectrons,
                           int const numberOfNegativeMuons,
                           int const numberOfPositiveMuons )
    /* this combines the acceptances from ewinoCascade & vectorCascade. Whether
     * the relevant scoloreds are particles or antiparticles should be decided
     * by the derived class using this function.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    bool
    chargeAsymmetricCharginoSignal( int const numberOfNegativeElectrons,
                                    int const numberOfPositiveElectrons,
                                    int const numberOfNegativeMuons,
                                    int const numberOfPositiveMuons )
    /* this returns true if the acceptance should be halved because only 1
     * charge of chargino should be used. (it returns false if the
     * electroweakino is a neutralino.)
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };


  // this is derived class for direct squark to electroweakino fullCascades.
  class sxFullCascade : public fullCascade
  {
  public:
    sxFullCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~sxFullCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    setProperties( input_handler const* const shortcut,
                   CppSLHA::particle_property_set const* const initialSquark,
                   double const beamEnergy,
                   electroweakCascade* const ewinoCascade )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    getBrToEwino( std::list< int > const* excludedSmParticles )
    /* this works out the branching ratio for the decays of colored sparticles
     * down to the 1st electroweakino in the cascade (including any decays to
     * vector bosons + squarks), using the given list of excluded SM particles
     * to exclude unwanted parts of the BR.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual double
    getAcceptance( bool const scoloredIsNotAntiparticle,
                   acceptanceCutSet* const cuts,
                   int const numberOfAdditionalJets,
                   int const numberOfNegativeElectrons,
                   int const numberOfPositiveElectrons,
                   int const numberOfNegativeMuons,
                   int const numberOfPositiveMuons )
    /* this calls the appropriate functions on ewinoCascade &, if not NULL,
     * vectorCascade, to build the required acceptance, taking into account
     * whether the charges should be swapped if scoloredIsNotAntiparticle is
     * false.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    //protected:
    // nothing.
  };

  // this is derived class for direct gluino to electroweakino fullCascades.
  class gxFullCascade : public fullCascade
  {
  public:
    gxFullCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~gxFullCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    setProperties( input_handler const* const shortcut,
                   double const beamEnergy,
                   electroweakCascade* const ewinoCascade )
    /* code after the classes in this .hpp file, or in the .cpp file. */;


    double
    getBrToEwino( std::list< int > const* excludedSmParticles )
    /* this works out the branching ratio for the decays of colored sparticles
     * down to the 1st electroweakino in the cascade (including any decays to
     * vector bosons + squarks), using the given list of excluded SM particles
     * to exclude unwanted parts of the BR.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual double
    getAcceptance( bool const scoloredIsNotAntiparticle,
                   acceptanceCutSet* const cuts,
                   int const numberOfAdditionalJets,
                   int const numberOfNegativeElectrons,
                   int const numberOfPositiveElectrons,
                   int const numberOfNegativeMuons,
                   int const numberOfPositiveMuons )
    /* this calls the appropriate functions on ewinoCascade, summing over the 2
     * charge versions if the electroweakino is a chargino. it ignores the value
     * of scoloredIsNotAntiparticle.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  //protected:
    // nothing.
  };

  // this is derived class for squark to gluino to electroweakino fullCascades.
  class sjgxFullCascade : public fullCascade
  {
  public:
    sjgxFullCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~sjgxFullCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    setProperties( CppSLHA::particle_property_set const* const initialSquark,
                   gxFullCascade* const gxCascade )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    getBrToEwino( std::list< int > const* excludedSmParticles )
    /* this works out the branching ratio for the decays of colored sparticles
     * down to the 1st electroweakino in the cascade (including any decays to
     * vector bosons + squarks), using the given list of excluded SM particles
     * to exclude unwanted parts of the BR.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual double
    getAcceptance( bool const scoloredIsNotAntiparticle,
                   acceptanceCutSet* const cuts,
                   int const numberOfAdditionalJets,
                   int const numberOfNegativeElectrons,
                   int const numberOfPositiveElectrons,
                   int const numberOfNegativeMuons,
                   int const numberOfPositiveMuons )
    /* this calls the appropriate functions on ewinoCascade &, if not NULL,
     * vectorCascade, to build the required acceptance, taking into account
     * whether the charges should be swapped if scoloredIsNotAntiparticle is
     * false.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    gxFullCascade* gxCascade;
  };

  // this is derived class for gluino to squark to electroweakino fullCascades.
  class gjsxFullCascade : public fullCascade
  {
  public:
    gjsxFullCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~gjsxFullCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    setProperties( sxFullCascade* const sxCascade )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    getBrToEwino( std::list< int > const* excludedSmParticles )
    /* this works out the branching ratio for the decays of colored sparticles
     * down to the 1st electroweakino in the cascade (including any decays to
     * vector bosons + squarks), using the given list of excluded SM particles
     * to exclude unwanted parts of the BR.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual double
    getAcceptance( bool const scoloredIsNotAntiparticle,
                   acceptanceCutSet* const cuts,
                   int const numberOfAdditionalJets,
                   int const numberOfNegativeElectrons,
                   int const numberOfPositiveElectrons,
                   int const numberOfNegativeMuons,
                   int const numberOfPositiveMuons )
    /* this calls the appropriate functions on ewinoCascade &, if not NULL,
     * vectorCascade, to build the required acceptance, taking into account
     * whether the charges should be swapped if scoloredIsNotAntiparticle is
     * false.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    sxFullCascade* sxCascade;
  };

  // this is derived class for squark to gluino to squark to electroweakino
  // fullCascades, treated by jets+MET grid as gsx.
  class sjgjsxFullCascade : public fullCascade
  {
  public:
    sjgjsxFullCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~sjgjsxFullCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    setProperties( CppSLHA::particle_property_set const* const initialSquark,
                   gjsxFullCascade* const gjsxCascade )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    getBrToEwino( std::list< int > const* excludedSmParticles )
    /* this works out the branching ratio for the decays of colored sparticles
     * down to the 1st electroweakino in the cascade (including any decays to
     * vector bosons + squarks), using the given list of excluded SM particles
     * to exclude unwanted parts of the BR.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual double
    getAcceptance( bool const scoloredIsNotAntiparticle,
                   acceptanceCutSet* const cuts,
                   int const numberOfAdditionalJets,
                   int const numberOfNegativeElectrons,
                   int const numberOfPositiveElectrons,
                   int const numberOfNegativeMuons,
                   int const numberOfPositiveMuons )
    /* this calls the appropriate functions on ewinoCascade &, if not NULL,
     * vectorCascade, to build the required acceptance, taking into account
     * whether the charges should be swapped if scoloredIsNotAntiparticle is
     * false.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    gjsxFullCascade* gjsxCascade;
  };

  // this is derived class for squark to vector boson plus other squark to
  // electroweakino fullCascades, treated by jets+MET grid as sx.
  class svsxFullCascade : public fullCascade
  {
  public:
    svsxFullCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~svsxFullCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    setProperties( CppSLHA::particle_property_set const* const initialSquark,
                   electroweakCascade* const vectorCascade,
                   sxFullCascade* const sxCascade )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    getBrToEwino( std::list< int > const* excludedSmParticles )
    /* this works out the branching ratio for the decays of colored sparticles
     * down to the 1st electroweakino in the cascade (including any decays to
     * vector bosons + squarks), using the given list of excluded SM particles
     * to exclude unwanted parts of the BR.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual double
    getAcceptance( bool const scoloredIsNotAntiparticle,
                   acceptanceCutSet* const cuts,
                   int const numberOfAdditionalJets,
                   int const numberOfNegativeElectrons,
                   int const numberOfPositiveElectrons,
                   int const numberOfNegativeMuons,
                   int const numberOfPositiveMuons )
    /* this calls the appropriate functions on ewinoCascade &, if not NULL,
     * vectorCascade, to build the required acceptance, taking into account
     * whether the charges should be swapped if scoloredIsNotAntiparticle is
     * false.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    sxFullCascade* sxCascade;
  };

  // this is derived class for gluino to vector boson plus squark to
  // electroweakino fullCascades, treated by jets+MET grid as sx.
  class gvsxFullCascade : public fullCascade
  {
  public:
    gvsxFullCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~gvsxFullCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    setProperties( electroweakCascade* const vectorCascade,
                   sxFullCascade* const sxCascade )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    getBrToEwino( std::list< int > const* excludedSmParticles )
    /* this works out the branching ratio for the decays of colored sparticles
     * down to the 1st electroweakino in the cascade (including any decays to
     * vector bosons + squarks), using the given list of excluded SM particles
     * to exclude unwanted parts of the BR.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual double
    getAcceptance( bool const scoloredIsNotAntiparticle,
                   acceptanceCutSet* const cuts,
                   int const numberOfAdditionalJets,
                   int const numberOfNegativeElectrons,
                   int const numberOfPositiveElectrons,
                   int const numberOfNegativeMuons,
                   int const numberOfPositiveMuons )
    /* this calls the appropriate functions on ewinoCascade &, if not NULL,
     * vectorCascade, to build the required acceptance, taking into account
     * whether the charges should be swapped if scoloredIsNotAntiparticle is
     * false.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    sxFullCascade* sxCascade;
  };

  // this is derived class for gluino to squark to vector boson plus other
  // squark to electroweakino fullCascades, treated by jets+MET grid as gsx.
  class gjsvsxFullCascade : public fullCascade
  {
  public:
    gjsvsxFullCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~gjsvsxFullCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    setProperties( svsxFullCascade* const svsxCascade )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    getBrToEwino( std::list< int > const* excludedSmParticles )
    /* this works out the branching ratio for the decays of colored sparticles
     * down to the 1st electroweakino in the cascade (including any decays to
     * vector bosons + squarks), using the given list of excluded SM particles
     * to exclude unwanted parts of the BR.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual double
    getAcceptance( bool const scoloredIsNotAntiparticle,
                   acceptanceCutSet* const cuts,
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
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    svsxFullCascade* svsxCascade;
  };

  // this is derived class for squark to vector boson plus gluino to
  // electroweakino fullCascades, treated by jets+MET grid as sgx.
  class svgxFullCascade : public fullCascade
  {
  public:
    svgxFullCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~svgxFullCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    setProperties( CppSLHA::particle_property_set const* const initialSquark,
                   electroweakCascade* const vectorCascade,
                   gxFullCascade* const gxCascade )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    getBrToEwino( std::list< int > const* excludedSmParticles )
    /* this works out the branching ratio for the decays of colored sparticles
     * down to the 1st electroweakino in the cascade (including any decays to
     * vector bosons + squarks), using the given list of excluded SM particles
     * to exclude unwanted parts of the BR.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual double
    getAcceptance( bool const scoloredIsNotAntiparticle,
                   acceptanceCutSet* const cuts,
                   int const numberOfAdditionalJets,
                   int const numberOfNegativeElectrons,
                   int const numberOfPositiveElectrons,
                   int const numberOfNegativeMuons,
                   int const numberOfPositiveMuons )
    /* this calls the appropriate functions on ewinoCascade &, if not NULL,
     * vectorCascade, to build the required acceptance, taking into account
     * whether the charges should be swapped if scoloredIsNotAntiparticle is
     * false.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    gxFullCascade* gxCascade;
  };

  // this is derived class for squark to vector boson plus other squark to
  // gluino to electroweakino fullCascades, treated by jets+MET grid as sgx.
  class svsjgxFullCascade : public fullCascade
  {
  public:
    svsjgxFullCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~svsjgxFullCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    setProperties( CppSLHA::particle_property_set const* const initialSquark,
                   electroweakCascade* const vectorCascade,
                   sjgxFullCascade* const sjgxCascade )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    getBrToEwino( std::list< int > const* excludedSmParticles )
    /* this works out the branching ratio for the decays of colored sparticles
     * down to the 1st electroweakino in the cascade (including any decays to
     * vector bosons + squarks), using the given list of excluded SM particles
     * to exclude unwanted parts of the BR.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual double
    getAcceptance( bool const scoloredIsNotAntiparticle,
                   acceptanceCutSet* const cuts,
                   int const numberOfAdditionalJets,
                   int const numberOfNegativeElectrons,
                   int const numberOfPositiveElectrons,
                   int const numberOfNegativeMuons,
                   int const numberOfPositiveMuons )
    /* this calls the appropriate functions on ewinoCascade &, if not NULL,
     * vectorCascade, to build the required acceptance, taking into account
     * whether the charges should be swapped if scoloredIsNotAntiparticle is
     * false.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    sjgxFullCascade* sjgxCascade;
  };

  /* this is derived class for squark to vector boson plus other squark to
   * gluino to yet another squark to electroweakino fullCascades, treated by
   * jets+MET grid as gsx.
   */
  class svsjgjsxFullCascade : public fullCascade
  {
  public:
    svsjgjsxFullCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~svsjgjsxFullCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    setProperties( CppSLHA::particle_property_set const* const initialSquark,
                   electroweakCascade* const vectorCascade,
                   sjgjsxFullCascade* const sjgjsxCascade )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    getBrToEwino( std::list< int > const* excludedSmParticles )
    /* this works out the branching ratio for the decays of colored sparticles
     * down to the 1st electroweakino in the cascade (including any decays to
     * vector bosons + squarks), using the given list of excluded SM particles
     * to exclude unwanted parts of the BR.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual double
    getAcceptance( bool const scoloredIsNotAntiparticle,
                   acceptanceCutSet* const cuts,
                   int const numberOfAdditionalJets,
                   int const numberOfNegativeElectrons,
                   int const numberOfPositiveElectrons,
                   int const numberOfNegativeMuons,
                   int const numberOfPositiveMuons )
    /* this calls the appropriate functions on ewinoCascade &, if not NULL,
     * vectorCascade, to build the required acceptance, taking into account
     * whether the charges should be swapped if scoloredIsNotAntiparticle is
     * false.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    sjgjsxFullCascade* sjgjsxCascade;
  };

  // this is derived class for squark to vector boson plus gluino to other
  // squark to electroweakino fullCascades, treated by jets+MET grid as gsx.
  class svgjsxFullCascade : public fullCascade
  {
  public:
    svgjsxFullCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~svgjsxFullCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    setProperties( CppSLHA::particle_property_set const* const initialSquark,
                   electroweakCascade* const vectorCascade,
                   gjsxFullCascade* const gjsxCascade )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    getBrToEwino( std::list< int > const* excludedSmParticles )
    /* this works out the branching ratio for the decays of colored sparticles
     * down to the 1st electroweakino in the cascade (including any decays to
     * vector bosons + squarks), using the given list of excluded SM particles
     * to exclude unwanted parts of the BR.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual double
    getAcceptance( bool const scoloredIsNotAntiparticle,
                   acceptanceCutSet* const cuts,
                   int const numberOfAdditionalJets,
                   int const numberOfNegativeElectrons,
                   int const numberOfPositiveElectrons,
                   int const numberOfNegativeMuons,
                   int const numberOfPositiveMuons )
    /* this calls the appropriate functions on ewinoCascade &, if not NULL,
     * vectorCascade, to build the required acceptance, taking into account
     * whether the charges should be swapped if scoloredIsNotAntiparticle is
     * false.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    gjsxFullCascade* gjsxCascade;
  };

  // this is derived class for squark to gluino to vector boson plus other
  // squark to electroweakino fullCascades, treated by jets+MET grid as gsx.
  class sjgvsxFullCascade : public fullCascade
  {
  public:
    sjgvsxFullCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~sjgvsxFullCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    setProperties( CppSLHA::particle_property_set const* const initialSquark,
                   gvsxFullCascade* const gvsxCascade )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    getBrToEwino( std::list< int > const* excludedSmParticles )
    /* this works out the branching ratio for the decays of colored sparticles
     * down to the 1st electroweakino in the cascade (including any decays to
     * vector bosons + squarks), using the given list of excluded SM particles
     * to exclude unwanted parts of the BR.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual double
    getAcceptance( bool const scoloredIsNotAntiparticle,
                   acceptanceCutSet* const cuts,
                   int const numberOfAdditionalJets,
                   int const numberOfNegativeElectrons,
                   int const numberOfPositiveElectrons,
                   int const numberOfNegativeMuons,
                   int const numberOfPositiveMuons )
    /* this calls the appropriate functions on ewinoCascade &, if not NULL,
     * vectorCascade, to build the required acceptance, taking into account
     * whether the charges should be swapped if scoloredIsNotAntiparticle is
     * false.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    gvsxFullCascade* gvsxCascade;
  };

  /* this is derived class for squark to gluino to other squark to vector boson
   * plus yet another squark to electroweakino fullCascades, treated by
   * jets+MET grid as gsx.
   */
  class sjgjsvsxFullCascade : public fullCascade
  {
  public:
    sjgjsvsxFullCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~sjgjsvsxFullCascade()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    void
    setProperties( CppSLHA::particle_property_set const* const initialSquark,
                        gjsvsxFullCascade* const gjsvsxCascade )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    double
    getBrToEwino( std::list< int > const* excludedSmParticles )
    /* this works out the branching ratio for the decays of colored sparticles
     * down to the 1st electroweakino in the cascade (including any decays to
     * vector bosons + squarks), using the given list of excluded SM particles
     * to exclude unwanted parts of the BR.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual double
    getAcceptance( bool const scoloredIsNotAntiparticle,
                   acceptanceCutSet* const cuts,
                   int const numberOfAdditionalJets,
                   int const numberOfNegativeElectrons,
                   int const numberOfPositiveElectrons,
                   int const numberOfNegativeMuons,
                   int const numberOfPositiveMuons )
    /* this calls the appropriate functions on ewinoCascade &, if not NULL,
     * vectorCascade, to build the required acceptance, taking into account
     * whether the charges should be swapped if scoloredIsNotAntiparticle is
     * false.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    gjsvsxFullCascade* gjsvsxCascade;
  };


  /* this class holds the possible cascades for a given scolored, as well as a
   * pointer to a std::list of the other fullCascadeSet instances, which get
   * ordered by mass of initial colored sparticle by whichever fullCascadeSet
   * is called 1st.
   */
  class fullCascadeSet : public readied_for_new_point
  {
  public:
    fullCascadeSet( input_handler const* const shortcut,
                   CppSLHA::particle_property_set const* const initialScolored,
                    std::list< fullCascadeSet* >* const cascadeSetList,
                    double const beamEnergy )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~fullCascadeSet()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    std::vector< fullCascade* >*
    getOpenCascades()
    // this should call setUpCascades() if it needs to be readied for this
    // point.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    minimalAllocationVector< fullCascade > openCascades;
    input_handler const* const shortcut;
    std::list< fullCascadeSet* >* const cascadeSetList;
    double const beamEnergy;

    void
    setUpCascades()
    // this clears openCascades, then sets it to be filled with all open
    // fullCascades that initialScolored has for this point.
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };



  // this class holds all the full_cascades & gives out pointers to a vector of
  // valid full_cascades which is reset with each updated point.
  class full_cascade_handler : public readied_for_new_point
  {

  public:

    full_cascade_handler()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~full_cascade_handler()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    std::vector< full_cascade* > const*
    get_valid_cascades()
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:

    std::vector< full_cascade* > valid_cascades;

    // collections of the full_cascades also.

    void
    find_all_valid_cascades()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  };





  // inline functions:

  inline fullCascade*
  fullCascade::getNewFullCascade()
  // this just returns a new fullCascade using the default constructor.
  {
    return new fullCascade();
  }

  inline std::pair< CppSLHA::particle_property_set const*,
                                                   int >*
  fullCascade::getNewPair()
  {
    return new std::pair< CppSLHA::particle_property_set const*,
                          int >( NULL,
                                 CppSLHA::CppSLHA_global::really_wrong_value );
  }

  input_handler const*
  fullCascade::getShortcut()
  const
  {
    return shortcut;
  }

  inline minimalAllocationVector< std::pair<
                                         CppSLHA::particle_property_set const*,
                                             int > >*
  fullCascade::getCascadeDefiner()
  {
    return &cascadeDefiner;
  }

  inline double
  fullCascade::getAcceptanceWithMinimumJets(
                                          bool const scoloredIsNotAntiparticle,
                                             acceptanceCutSet* const cuts,
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
    double returnDouble = 0.0;
    for( int numberOfJets
         = ( numberOfSmFermionsFromElectroweakDecaysPerFullCascade
             - numberOfNegativeElectrons - numberOfPositiveElectrons
             - numberOfNegativeMuons - numberOfPositiveMuons );
         minimumNumberOfAdditionalJets <= numberOfJets;
         --numberOfJets )
    {
      returnDouble += getAcceptance( scoloredIsNotAntiparticle,
                                          cuts,
                                          numberOfJets,
                                          numberOfNegativeElectrons,
                                          numberOfPositiveElectrons,
                                          numberOfNegativeMuons,
                                          numberOfPositiveMuons );
    }
    return returnDouble;
  }

  inline bool
  fullCascade::chargeAsymmetricCharginoSignal(
                                           int const numberOfNegativeElectrons,
                                           int const numberOfPositiveElectrons,
                                               int const numberOfNegativeMuons,
                                              int const numberOfPositiveMuons )
  /* this returns true if the acceptance should be halved because only 1
   * charge of chargino should be used. (it returns false if the
   * electroweakino is a neutralino.)
   */
  {
    if( ewinoCascade->getElectroweakDecayer()->counts_as_self_conjugate()
        ||
        ( ( numberOfNegativeElectrons == numberOfPositiveElectrons )
          &&
          ( numberOfNegativeMuons == numberOfPositiveMuons ) ) )
    {
      return false;
    }
    else
    {
      return true;
    }
  }

  inline void
  fullCascade::copyCascadeDefiner(
     minimalAllocationVector< std::pair< CppSLHA::particle_property_set const*,
                                         int > > cascadeDefiner )
  {
    this->cascadeDefiner.clear();
    // reset this->cascadeDefiner.
    for( int cascadeFiller = 0;
         cascadeDefiner->getSize() > cascadeFiller;
         ++cascadeFiller )
    {
      cascadeSegment = this->cascadeDefiner.addNewAtEnd();
      cascadeSegment->first
      = cascadeDefiner->getPointer( cascadeFiller )->first;
      cascadeSegment->second
      = cascadeDefiner->getPointer( cascadeFiller )->second;
    }
  }

  inline void
  sxFullCascade::setProperties( input_handler const* const shortcut,
                     CppSLHA::particle_property_set const* const initialSquark,
                                double const beamEnergy,
                                electroweakCascade* const ewinoCascade )
  {
    this->shortcut = shortcut;
    this->initialScolored = initialSquark;
    this->beamEnergy = beamEnergy;
    this->ewinoCascade = ewinoCascade;
    soughtDecayProductList.front() =
                         ewinoCascade->getElectroweakDecayer()->get_PDG_code();
    cascadeDefiner.clear();
    // reset cascadeDefiner.
    cascadeSegment = cascadeDefiner.addNewAtEnd();
    cascadeSegment->first = ewinoCascade->getElectroweakDecayer();
    // sx means only 1 decay to be recorded.
    cascadeSegment->second = 2;
    // sx also means that the decay is 2-body.
  }

  inline double
  sxFullCascade::getBrToEwino( std::list< int > const* excludedSmParticles )
  /* this works out the branching ratio for the decays of colored sparticles
   * down to the 1st electroweakino in the cascade (including any decays to
   * vector bosons + squarks), using the given list of excluded SM particles
   * to exclude unwanted parts of the BR.
   */
  {
    return initialScolored->inspect_direct_decay_handler(
                    )->get_branching_ratio_for_subset( &soughtDecayProductList,
                                                       excludedSmParticles );
  }



  inline void
  gxFullCascade::setProperties( input_handler const* const shortcut,
                                double const beamEnergy,
                                electroweakCascade* const ewinoCascade )
  {
    this->shortcut = shortcut;
    this->initialScolored = shortcut->get_gluino();
    this->beamEnergy = beamEnergy;
    this->ewinoCascade = ewinoCascade;
    soughtDecayProductList.front() =
                         ewinoCascade->getElectroweakDecayer()->get_PDG_code();
    cascadeDefiner.clear();
    // reset cascadeDefiner.
    cascadeSegment = cascadeDefiner.addNewAtEnd();
    cascadeSegment->first = ewinoCascade->getElectroweakDecayer();
    // gx means only 1 decay to be recorded.
    cascadeSegment->second = 3;
    // gx also means that the decay is 2-body.
  }

  inline double
  gxFullCascade::getBrToEwino( std::list< int > const* excludedSmParticles )
  /* this works out the branching ratio for the decays of colored sparticles
   * down to the 1st electroweakino in the cascade (including any decays to
   * vector bosons + squarks), using the given list of excluded SM particles
   * to exclude unwanted parts of the BR.
   */
  {
    if( ewinoCascade->getElectroweakDecayer()->counts_as_self_conjugate() )
      // if the gluino has only 1 charge version of this decay...
    {
      return initialScolored->inspect_direct_decay_handler(
                    )->get_branching_ratio_for_subset( &soughtDecayProductList,
                                                       excludedSmParticles );
    }
    else
    {
      return ( 2.0 * initialScolored->inspect_direct_decay_handler(
                    )->get_branching_ratio_for_subset( &soughtDecayProductList,
                                                       excludedSmParticles ) );
    }
  }



  inline void
  sjgxFullCascade::setProperties(
                     CppSLHA::particle_property_set const* const initialSquark,
                                  gxFullCascade* const gxCascade )
  {
    this->gxCascade = gxCascade;
    this->shortcut = gxCascade->shortcut;
    this->initialScolored = initialSquark;
    this->beamEnergy = gxCascade->beamEnergy;
    this->ewinoCascade = gxCascade->ewinoCascade;
    soughtDecayProductList.front() = CppSLHA::PDG_code::gluino;
    copyCascadeDefiner( gxCascade->getCascadeDefiner() );
    // reset cascadeDefiner & fill it with a copy of gxCascade's version.
    cascadeSegment = cascadeDefiner.addNewAtEnd();
    cascadeSegment->first = gxCascade->getInitialScolored();
    // sjgx means 1 more decay to be recorded on top of those of the given
    // cascade.
    cascadeSegment->second = 2;
    // sjgx also means that the initial decay is 2-body.
  }

  inline double
  sjgxFullCascade::getBrToEwino( std::list< int > const* excludedSmParticles )
  /* this works out the branching ratio for the decays of colored sparticles
   * down to the 1st electroweakino in the cascade (including any decays to
   * vector bosons + squarks), using the given list of excluded SM particles
   * to exclude unwanted parts of the BR.
   */
  {
    return ( initialScolored->inspect_direct_decay_handler(
                    )->get_branching_ratio_for_subset( &soughtDecayProductList,
                                                       excludedSmParticles )
             * gxCascade->getBrToEwino( excludedSmParticles ) );
  }

  inline double
  sjgxFullCascade::getAcceptance( bool const scoloredIsNotAntiparticle,
                                  acceptanceCutSet* const cuts,
                                  int const numberOfAdditionalJets,
                                  int const numberOfNegativeElectrons,
                                  int const numberOfPositiveElectrons,
                                  int const numberOfNegativeMuons,
                                  int const numberOfPositiveMuons )
  /* this calls the appropriate functions on ewinoCascade &, if not NULL,
   * vectorCascade, to build the required acceptance, taking into account
   * whether the charges should be swapped if scoloredIsNotAntiparticle is
   * false.
   */
  {
    return gxCascade->getAcceptance( scoloredIsNotAntiparticle,
                                     cuts,
                                     numberOfAdditionalJets,
                                     numberOfNegativeElectrons,
                                     numberOfPositiveElectrons,
                                     numberOfNegativeMuons,
                                     numberOfPositiveMuons );
  }



  inline std::vector< fullCascade* >*
  fullCascadeSet::getOpenCascades()
  // this should call setUpCascades() if it needs to be readied for this
  // point.
  {
    if( needs_to_prepare_for_this_point() )
    {
      setUpCascades();
    }
    return cascadeSet.getOpenCascades();
  }
}  // end of LHC_FASER namespace.

#endif /* LHC_FASER_FULL_CASCADE_STUFF_HPP_ */
