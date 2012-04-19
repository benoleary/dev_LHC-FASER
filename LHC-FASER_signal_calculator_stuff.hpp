/*
 * LHC-FASER_signal_calculator_stuff.hpp
 *
 *  Created on: 17 Nov 2010
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

#ifndef LHC_FASER_SIGNAL_CALCULATOR_STUFF_HPP_
#define LHC_FASER_SIGNAL_CALCULATOR_STUFF_HPP_

#include "LHC-FASER_global_stuff.hpp"
#include "LHC-FASER_signal_data_collection_stuff.hpp"


/* NOTES FOR WRITING NEW SIGNALS:
 *
 * I have tried to make it as straightforward as possible to write new signals,
 * but inevitably there is going to have to be some new code written for each
 * new signal. Until I get around to writing a better version of these notes,
 * you'll have to go by the existing signals as examples.
 *
 * Things to remember:
 * - valueForCurrentCascades(...) has to be written for any class derived from
 *   signalCalculator, even if it doesn't use it because it over-rides
 *   calculateValue(...) itself. however, most signals should just have their
 *   own valueForCurrentCascades(...), unless they do something strange.
 * - make sure that exclusiveBrHandler::alwaysNeglectedDecays doesn't
 *   interfere with the new signal. exclusiveBrHandler returns a NULL pointer
 *   if the decay is in alwaysNeglectedDecays.
 * - fire hot, hot hurt.
 */

namespace LHC_FASER
{
  // this abstract base class encapsulates the general lepton acceptance
  // separate from the jet acceptance.
  class leptonAcceptanceForCascadePair
  {
  public:
    leptonAcceptanceForCascadePair();
    virtual
    ~leptonAcceptanceForCascadePair();

    virtual double
    withExactlyNJets( signalDefinitionSet const* const signalDefinitions,
                      int const exactNumberOfJets,
                      fullCascade* const firstCascade,
                      bool const firstIsNotAntiparticle,
                      fullCascade* const secondCascade,
                      bool const secondIsNotAntiparticle ) const = 0;
    virtual double
    withAtLeastNJets( signalDefinitionSet const* const signalDefinitions,
                      int const minimumNumberOfJets,
                      fullCascade* const firstCascade,
                      bool const firstIsNotAntiparticle,
                      fullCascade* const secondCascade,
                      bool const secondIsNotAntiparticle ) const;


  protected:
    static int const maximumNumberOfJets;
  };


  namespace leptonAcceptanceStyle
  {
    // this class derived from leptonAcceptanceForCascadePair is the most
    // general I could think of.
    class fullySpecified : public leptonAcceptanceForCascadePair
    {
    public:
      fullySpecified( int const numberOfNegativeElectrons,
                      int const numberOfPositiveElectrons,
                      int const numberOfNegativeMuons,
                      int const numberOfPositiveMuons );
      virtual
      ~fullySpecified();

      virtual double
      withExactlyNJets( signalDefinitionSet const* const signalDefinitions,
                        int const exactNumberOfJets,
                        fullCascade* const firstCascade,
                        bool const firstIsNotAntiparticle,
                        fullCascade* const secondCascade,
                        bool const secondIsNotAntiparticle ) const;


    protected:
      int const numberOfNegativeElectrons;
      int const numberOfPositiveElectrons;
      int const numberOfNegativeMuons;
      int const numberOfPositiveMuons;
    };


    // this class derived from leptonAcceptanceForCascadePair just returns
    // 100% lepton / additional jet acceptance.
    class noLeptonCutNorExtraJetCut : public leptonAcceptanceForCascadePair
    {
    public:
      noLeptonCutNorExtraJetCut();
      virtual
      ~noLeptonCutNorExtraJetCut();

      virtual double
      withExactlyNJets( signalDefinitionSet const* const signalDefinitions,
                        int const minimumNumberOfJets,
                        fullCascade* const firstCascade,
                        bool const firstIsNotAntiparticle,
                        fullCascade* const secondCascade,
                        bool const secondIsNotAntiparticle ) const;
      virtual double
      withAtLeastNJets( signalDefinitionSet const* const signalDefinitions,
                        int const minimumNumberOfJets,
                        fullCascade* const firstCascade,
                        bool firstIsNotAntiparticle,
                        fullCascade* const secondCascade,
                        bool const secondIsNotAntiparticle ) const;
    };


    // this class derived from leptonAcceptanceForCascadePair sums up over
    // all combinations of electrons & muons & charges.
    class chargeAndFlavorSummed : public leptonAcceptanceForCascadePair
    {
    public:
      chargeAndFlavorSummed( int const numberOfLeptons );
      virtual
      ~chargeAndFlavorSummed();

      virtual double
      withExactlyNJets( signalDefinitionSet const* const signalDefinitions,
                        int const exactNumberOfJets,
                        fullCascade* const firstCascade,
                        bool const firstIsNotAntiparticle,
                        fullCascade* const secondCascade,
                        bool const secondIsNotAntiparticle ) const;


    protected:
      int const numberOfLeptons;
    };


    // this class derived from leptonAcceptanceForCascadePair sums up over
    // all combinations of charges, distinguishing electrons from muons.
    class chargeSummed : public leptonAcceptanceForCascadePair
    {
    public:
      chargeSummed( int const numberOfElectrons,
                    int const numberOfMuons );
      virtual
      ~chargeSummed();

      virtual double
      withExactlyNJets( signalDefinitionSet const* const signalDefinitions,
                        int const exactNumberOfJets,
                        fullCascade* const firstCascade,
                        bool const firstIsNotAntiparticle,
                        fullCascade* const secondCascade,
                        bool const secondIsNotAntiparticle ) const;


    protected:
      int const numberOfElectrons;
      int const numberOfMuons;
    };


    // this class derived from leptonAcceptanceForCascadePair sums up over
    // all combinations of electrons & muons charges, distinguishing charges.
    class flavorSummed : public leptonAcceptanceForCascadePair
    {
    public:
      flavorSummed( int const numberOfNegativeLeptons,
                    int const numberOfPositiveLeptons );
      virtual
      ~flavorSummed();

      virtual double
      withExactlyNJets( signalDefinitionSet const* const signalDefinitions,
                        int const exactNumberOfJets,
                        fullCascade* const firstCascade,
                        bool const firstIsNotAntiparticle,
                        fullCascade* const secondCascade,
                        bool const secondIsNotAntiparticle ) const;


    protected:
      int const numberOfNegativeLeptons;
      int const numberOfPositiveLeptons;
    };


    // this class derived from leptonAcceptanceForCascadePair returns the
    // acceptance for *exactly one* OSSF-OSDF pair.
    class ossfMinusOsdf : public leptonAcceptanceForCascadePair
    {
    public:
      ossfMinusOsdf();
      virtual
      ~ossfMinusOsdf();

      virtual double
      withExactlyNJets( signalDefinitionSet const* const signalDefinitions,
                        int const exactNumberOfJets,
                        fullCascade* const firstCascade,
                        bool const firstIsNotAntiparticle,
                        fullCascade* const secondCascade,
                        bool const secondIsNotAntiparticle ) const;
    };


    // this class derived from leptonAcceptanceForCascadePair returns the
    // acceptance for *exactly one* same-sign, same-flavor pair.
    class sameSignSameFlavor : public leptonAcceptanceForCascadePair
    {
    public:
      sameSignSameFlavor();
      virtual
      ~sameSignSameFlavor();

      virtual double
      withExactlyNJets( signalDefinitionSet const* const signalDefinitions,
                        int const exactNumberOfJets,
                        fullCascade* const firstCascade,
                        bool const firstIsNotAntiparticle,
                        fullCascade* const secondCascade,
                        bool const secondIsNotAntiparticle ) const;
    };


  }  // end of leptonAcceptanceStyle namespace


  /* this is an abstract base class to do the actual calculation of the value
   * of the event rate for a signal. the main differences in derived classes
   * are what acceptances they access from the handlers, how they put them
   * together, & how they estimate the uncertainty factor.
   */
  class signalCalculator
  {
  public:
    signalCalculator( signalDefinitionSet const* const signalDefinitions );
    virtual
    ~signalCalculator();

    /* actually, this doesn't work, for various reasons. currently C++ does not
     * allow virtual static functions. I was only going to use it to force me
     * to remember to write the appropriate static functions in derived classes
     * anyway.
     */
    /*virtual static
    signalCalculator*
    getCalculator( std::string const* const argumentString,
                   signalDefinitionSet const* const signalDefinitions )
    = 0;*/
    /* this should be over-written in each derived class to construct a new
     * instance of the derived class based on arguments, & return a pointer to
     * the new instance of the derived class, or NULL if the string could not
     * be parsed as a correct key for the derived class.
     */
    virtual bool
    calculateValue( double* signalValue,
                    double* uncertaintyFactor );
    /* this calculates the event rate for the signal & puts its value in
     * signalValue, & puts an estimate for the uncertainty into
     * uncertaintyFactor, & returns true if it did all this successfully.
     *
     * by default, this calls goThroughCascadesNormally( signalValue,
     *                                                   uncertaintyFactor )
     * which calls valueForCurrentCascades(...) for each pair of cascades.
     * it is expected that valueForCurrentCascades(...) for each derived class
     * - checks the cross-section times BRs for the cascade pair - if it's high
     *   enough, it proceeds
     * - checks the lepton acceptances for the cascade pair - if they're high
     *   enough in combination with the BRs & the cross-section, it proceeds
     * - obtains its jet acceptance from its kinematics table, checks
     *   cross-section * BRs * acceptances, & if high enough, returns the
     *   product, or 0.0 otherwise.
     */


  protected:
    static bool
    parseBeamEnergy( std::string& argumentString,
                     signalDefinitionSet* const signalDefinitions );
    /* this looks for "_7TeV", "_07TeV", "_10TeV", or "_14TeV", removes it from
     * argumentString & sets the beam energy appropriately. returns false if it
     * could not find any (& does not modify argumentString in this case).
     */
    static leptonAcceptanceForCascadePair*
    parseLeptonAcceptance( std::string& argumentString,
                           signalDefinitionSet* const signalDefinitions );
    /* this looks for strings encoding the type of lepton cuts to use. the
     * strings are, where # stands for any string representing an integer, in the
     * order in which they are checked:
     * "_ossf" : ossfMinusOsdf
     * "_sssf" : sameSignSameFlavor
     * "_noExtraCut" : noLeptonCutNorExtraJetCut
     * "_#l" : chargeAndFlavorSummed
     * "_#lm#lp" : flavorSummed
     * "_#epm#mpm" : chargeSummed
     * "_#em#ep#mm#mp" : fullySpecified
     * parseLeptonTransverseMomentumCuts is then called on the remainder of
     * argumentString. finally, it creates a new leptonAcceptanceForCascadePair
     * & returns a pointer to it. NULL is returned if argumentString could not
     * be properly interpretted.
     */
    static bool
    separateOutIntegerSubstring( int& destinationInt,
                                 std::string& argumentString );
    /* this checks if argumentString begins with a numeric character,
     * & if so destinationInt is set to the value of the int made from the
     * unbroken sequences of numeric characters from the start of
     * argumentString up to the first non-numeric character, & these characters
     * are removed from the start of argumentString, & true is returned.
     * otherwise argumentString is left as it is & false is returned.
     */
    static bool
    separateOutDecimalSubstring( double& destinationDouble,
                                 std::string& argumentString );
    /* this uses separateOutIntegerSubstring to check if argumentString begins
     * with a numeric character or '.'. destinationDouble is then set to be
     * equal to the int that separateOutIntegerSubstring found (or 0 if
     * argumentString begins with '.'), & then if there is an int after the '.'
     * (if there is a '.'), that becomes the part of destinationDouble after
     * the decimal point.
     */
    static bool
    parseLeptonTransverseMomentumCuts( std::string& argumentString,
                                signalDefinitionSet* const signalDefinitions );
    /* this looks for "_pTl" then a double then "GeV", then a subsequent
     * double, & a subsequent "GeV", interpretted as the primary & secondary
     * lepton cuts respectively. it returns true if there is nothing more in
     * argumentString, false otherwise.
     */
    static bool
    findAndRemoveSubstring( std::string const& soughtSubstring,
                            std::string& searchedString );
    /* this removes soughtSubstring from searchedString if found as a
     * substring, & returns true. if it is not found, searchedString is left
     * untouched & false is returned.
     */

    signalDefinitionSet const* const signalDefinitions;
    inputHandler const* const inputShortcut;
    std::vector< productionChannelPointerSet* > productionChannels;
    //fullCascade* firstCascade;
    //fullCascade* secondCascade;
    std::vector< fullCascade* >* firstCascades;
    std::vector< fullCascade* >* secondCascades;
    std::list< int > excludedFinalStateParticles;
    double firstCascadeBrToEwino;
    double secondCascadeBrToEwino;
    signedParticleShortcutPair const* initialPair;
    bool firstSparticleIsNotAntiparticle;
    bool secondSparticleIsNotAntiparticle;
    double subchannelCrossSectionTimesBrToEwinos;

    virtual double
    valueForCurrentCascades( fullCascade* firstCascade,
                             fullCascade* secondCascade )
    = 0;
    /* this function should assume that firstCascadeBrToEwino,
     * secondCascadeBrToEwino, initialPair, firstSparticleIsNotAntiparticle,
     * secondCascadeBrToEwino, & subchannelCrossSectionTimesBrToEwinos have all
     * been set correctly.
     * it is expected that valueForCurrentCascades(...) for each derived class
     * - checks the cross-section times BRs for the cascade pair - if it's high
     *   enough, it proceeds
     * - checks the lepton acceptances for the cascade pair - if they're high
     *   enough in combination with the BRs & the cross-section, it proceeds
     * - obtains its jet acceptance from its kinematics table, checks
     *   cross-section * BRs * acceptances, & if high enough, returns the
     *   product, or 0.0 otherwise.
     */
    bool
    goThroughCascadesNormally( double* signalValue,
                               double* uncertaintyFactor );
    // this calls valueForCurrentCascades(...) for each pair of cascades for
    // each production channel.
  };


  namespace signalClasses
  {
    /* this is a derived class which just always returns
     * CppSLHA::CppSLHA_global::really_wrong_value as its value. it is used for
     * malformed signal names.
     */
    class reallyWrongCalculator : public signalCalculator
    {
    public:
      reallyWrongCalculator(
                          signalDefinitionSet const* const signalDefinitions );
      virtual
      ~reallyWrongCalculator();

      bool
      calculateValue( double* signalValue,
                      double* uncertaintyFactor );
      // this always returns false, & always sets signalValue &
      // uncertaintyFactor to CppSLHA::CppSLHA_global::really_wrong_value.


    protected:
      virtual double
      valueForCurrentCascades( fullCascade* firstCascade,
                               fullCascade* secondCascade );
      // this shouldn't be called.
    };


    /* this is a derived class to test whether I break down the production
     * channels & recombine them correctly. all acceptances are set to 1.0 so
     * that I should just recover the total cross-section.
     */
    class sigmaBreakdownTest : public signalCalculator
    {
    public:
      static signalCalculator*
      getCalculator( std::string const& argumentString,
                     signalDefinitionSet* const signalDefinitions );
      /* this either returns a pointer to a new sigmaBreakdownTest with cuts
       * decided by the given string, or NULL if the string could not be parsed
       * properly.
       */

      ~sigmaBreakdownTest();

      bool
      calculateValue( double* const signalValue,
                      double* const uncertaintyFactor );
      // this over-rides the base version because it does things a bit
      // differently in a few of the steps.


    protected:
      double channelBrTotal;

      sigmaBreakdownTest( signalDefinitionSet const* const signalDefinitions );

      virtual double
      valueForCurrentCascades( fullCascade* firstCascade,
                               fullCascade* secondCascade );
      // this shouldn't be called.
    };


    /* this is a derived class to calculateValue the "Atlas 4 jets plus missing
     * transverse momentum plus some lepton cut" signal.
     * it takes the kinematics from the Atlas4jMET grid & combines them with
     * cascade decays leading to however many leptons passing whichever cut.
     */
    class atlasFourJetMetPlusGivenLeptonCuts : public signalCalculator
    {
    public:
      static int const jetAcceptanceGridTableColumn;
      // this is dependent on the format of the grids.
      static double const defaultExtraJetCut;
      // this is the standard cut for the jets beyond the hardest cut for this
      // signal as implemented in this code.
      static double const defaultPrimaryLeptonCut;
      static double const defaultSecondaryLeptonCut;
      // the default Atlas lepton transverse momentum cuts are 20.0 GeV
      // (for a single lepton to *pass*) & 10.0 GeV (for all others to *fail*).

      static signalCalculator*
      getCalculator( std::string const& argumentString,
                     signalDefinitionSet* const signalDefinitions );
      /* this either returns a pointer to a new
       * atlasFourJetMetPlusGivenLeptonCuts with cuts decided by the given
       * string, or NULL if the string could not be parsed properly.
       */

      ~atlasFourJetMetPlusGivenLeptonCuts();


    protected:
      leptonAcceptanceForCascadePair const* const leptonAcceptanceCalculator;
      jetAcceptanceTable* fourJetKinematics;
      jetAcceptanceTable* threeJetKinematics;
      jetAcceptanceTable* twoJetKinematics;
      /* the 4-jet signal is complicated enough that we also consider only 3 of
       * the 4 required hard jets coming from (showered) decays to
       * electroweakinos, with the 4th coming from a cascade decay of the
       * electroweakinos. hence we need an extra kinematics set (though it is
       * only used for direct jet acceptance). we also consider only 2 of the 4
       * hard jets coming from the decays to electroweakinos & the other 2
       * coming from the decays of the electroweakinos.
       */
      double fourJetAcceptance;
      double threeJetAcceptance;
      // these are used as each pairing of cascades from each production
      // channel is calculated:
      double subchannelValue;
      double subchannelZeroOrMoreJets;
      double subchannelOneOrMoreJets;
      double subchannelTwoOrMoreJets;

      atlasFourJetMetPlusGivenLeptonCuts(
                            signalDefinitionSet const* const signalDefinitions,
      leptonAcceptanceForCascadePair const* const leptonAcceptanceCalculator );

      virtual double
      valueForCurrentCascades( fullCascade* firstCascade,
                               fullCascade* secondCascade );
      // see base version's description.
    };


    /* this is a derived class to calculateValue the "Atlas 3 jets plus missing
     * transverse momentum plus some lepton cut" signal.
     * it takes the kinematics from the Atlas3jMET grid & combines them with
     * cascade decays leading to however many leptons passing whichever cut.
     */
    class atlasThreeJetMetPlusGivenLeptonCuts : public signalCalculator
    {
    public:
      static int const jetAcceptanceGridTableColumn;
      // this is dependent on the format of the grids.
      static double const defaultExtraJetCut;
      // this is the standard cut for the jets beyond the hardest cut for this
      // signal as implemented in this code.
      static double const defaultPrimaryLeptonCut;
      static double const defaultSecondaryLeptonCut;
      // the default Atlas lepton transverse momentum cuts are 20.0 GeV
      // (for a single lepton to *pass*) & 10.0 GeV (for all others to *fail*).

      static signalCalculator*
      getCalculator( std::string const& argumentString,
                     signalDefinitionSet* const signalDefinitions );
      /* this either returns a pointer to a new
       * atlasThreeJetMetPlusGivenLeptonCuts with cuts decided by the given
       * string, or NULL if the string could not be parsed properly.
       */

      ~atlasThreeJetMetPlusGivenLeptonCuts();


    protected:
      leptonAcceptanceForCascadePair const* const leptonAcceptanceCalculator;
      jetAcceptanceTable* threeJetKinematics;
      jetAcceptanceTable* twoJetKinematics;
      /* the 3-jet signal is complicated enough that we also consider only 2 of
       * the 3 required hard jets coming from (showered) decays to
       * electroweakinos, with the 3rd coming from a cascade decay of the
       * electroweakinos. hence we need an extra kinematics set (though it is
       * only used for direct jet acceptance).
       */
      double threeJetAcceptance;
      // these are used as each pairing of cascades from each production
      // channel is calculated:
      double subchannelValue;
      double subchannelZeroOrMoreJets;
      double subchannelOneOrMoreJets;

      atlasThreeJetMetPlusGivenLeptonCuts(
                            signalDefinitionSet const* const signalDefinitions,
      leptonAcceptanceForCascadePair const* const leptonAcceptanceCalculator );

      virtual double
      valueForCurrentCascades( fullCascade* firstCascade,
                               fullCascade* secondCascade );
      // see base version's description.
    };


    /* this is a derived class to calculateValue the "no cut on jets or missing
     * transverse momentum, but some lepton cut" signal.
     * it only takes cascade decays leading to however many leptons passing
     * whichever cut.
     */
    class noJetMetButGivenLeptonCuts : public signalCalculator
    {
    public:
      static double const defaultPrimaryLeptonCut;
      static double const defaultSecondaryLeptonCut;
      // the default Atlas lepton transverse momentum cuts are 20.0 GeV
      // (for a single lepton to *pass*) & 10.0 GeV (for all others to *fail*).

      static signalCalculator*
      getCalculator( std::string const& argumentString,
                     signalDefinitionSet* const signalDefinitions );
      /* this either returns a pointer to a new noJetMetButGivenLeptonCuts with
       * cuts decided by the given string, or NULL if the string could not be
       * parsed properly.
       */

      ~noJetMetButGivenLeptonCuts();


    protected:
      leptonAcceptanceForCascadePair const* const leptonAcceptanceCalculator;
      // these are used as each pairing of cascades from each production
      // channel is calculated:
      double subchannelValue;

      noJetMetButGivenLeptonCuts(
                            signalDefinitionSet const* const signalDefinitions,
      leptonAcceptanceForCascadePair const* const leptonAcceptanceCalculator );

      virtual double
      valueForCurrentCascades( fullCascade* firstCascade,
                               fullCascade* secondCascade );
      // see base version's description.
    };

  }  // end of signalClasses namespace.


  /* this is a class to handle each individual signal to be calculated.
   * it accesses numbers common to different signals through the
   * crossSectionHandler, kinematicsHandler & cascadeHandler classes.
   * it takes a string encoding what it should calculate, & constructs a
   * signalCalculator object to actually calculate the signal with the
   * given handlers. it also takes care of updating the signal & storing its
   * value & estimated uncertainty.
   */
  class signalHandler : public getsReadiedForNewPoint
  {
  public:
    typedef signalCalculator* (*signalCalculatorCreator)( std::string const&,
                                                  signalDefinitionSet* const );

    signalHandler( std::string const signalName,
                   double const crossSectionUnitFactor,
                   signalDefinitionSet const* const signalDefinitions );
    ~signalHandler();

    std::string const*
    getName()
    const;
    double
    getValue();
    double
    getUpperUncertainty();
    double
    getLowerUncertainty();
    double
    getUpperValue();
    double
    getLowerValue();
    double
    getLogUncertainty();


  protected:
    std::string signalName;
    signalCalculator* rateCalculator;
    signalShortcuts const* const inputShortcut;
    signalDefinitionSet signalPreparationDefinitions;
    double signalValue;
    double uncertaintyFactor;
    double crossSectionUnitFactor;
    // this is to allow for the user to specify event rates in fb, pb or nb.
  };





  // inline functions:


  inline double
  leptonAcceptanceForCascadePair::withAtLeastNJets(
                            signalDefinitionSet const* const signalDefinitions,
                                                 int const minimumNumberOfJets,
                                               fullCascade* const firstCascade,
                                             bool const firstIsNotAntiparticle,
                                              fullCascade* const secondCascade,
                                     bool const secondIsNotAntiparticle ) const
  {
    double returnValue( 0.0 );
    for( int numberOfJets( minimumNumberOfJets );
         maximumNumberOfJets >= numberOfJets;
         ++numberOfJets )
    {
      returnValue += withExactlyNJets( signalDefinitions,
                                       numberOfJets,
                                       firstCascade,
                                       firstIsNotAntiparticle,
                                       secondCascade,
                                       secondIsNotAntiparticle );
    }
    return returnValue;
  }


  namespace leptonAcceptanceStyle
  {
    inline double
    noLeptonCutNorExtraJetCut::withExactlyNJets(
                            signalDefinitionSet const* const signalDefinitions,
                                                 int const exactNumberOfJets,
                                               fullCascade* const firstCascade,
                                             bool const firstIsNotAntiparticle,
                                              fullCascade* const secondCascade,
                                     bool const secondIsNotAntiparticle ) const
    {
      return 1.0;
    }

    inline double
    noLeptonCutNorExtraJetCut::withAtLeastNJets(
                            signalDefinitionSet const* const signalDefinitions,
                                                 int const minimumNumberOfJets,
                                               fullCascade* const firstCascade,
                                             bool const firstIsNotAntiparticle,
                                              fullCascade* const secondCascade,
                                     bool const secondIsNotAntiparticle ) const
    {
      return 1.0;
    }

  }  // end of leptonAcceptanceStyle namespace



  inline bool
  signalCalculator::calculateValue( double* signalValue,
                                    double* uncertaintyFactor )
  /* this calculates the event rate for the signal & puts its value in
   * signalValue, & puts an estimate for the uncertainty into
   * uncertaintyFactor, & returns true if it did all this successfully.
   *
   * by default, this calls goThroughCascadesNormally( signalValue,
   *                                                   uncertaintyFactor )
   * which calls valueForCurrentCascades(...) for each pair of cascades.
   * it is expected that valueForCurrentCascades(...) for each derived class
   * - checks the cross-section times BRs for the cascade pair - if it's high
   *   enough, it proceeds
   * - checks the lepton acceptances for the cascade pair - if they're high
   *   enough in combination with the BRs & the cross-section, it proceeds
   * - obtains its jet acceptance from its kinematics table, checks
   *   cross-section * BRs * acceptances, & if high enough, returns the
   *   product, or 0.0 otherwise.
   */
  {
    return goThroughCascadesNormally( signalValue,
                                      uncertaintyFactor );
  }



  inline std::string const*
  signalHandler::getName()
  const
  {
   return &signalName;
  }

  inline double
  signalHandler::getValue()
  {
    // if the signal has not been calculated for this point, update
    // signalValue before returning it:
    if( needsToPrepareForThisPoint() )
      {
        rateCalculator->calculateValue( &signalValue,
                                        &uncertaintyFactor );
        signalValue *= crossSectionUnitFactor;
      }
    return signalValue;
  }

  inline double
  signalHandler::getUpperUncertainty()
  {
    return ( getValue() * ( uncertaintyFactor - 1.0 ) );
  }

  inline double
  signalHandler::getLowerUncertainty()
  {
    return ( getValue() * ( 1.0 - ( 1.0 / uncertaintyFactor ) ) );
  }

  inline double
  signalHandler::getUpperValue()
  {
    return ( getValue() * uncertaintyFactor );
  }

  inline double
  signalHandler::getLowerValue()
  {
    return ( getValue() / uncertaintyFactor );
  }

  inline double
  signalHandler::getLogUncertainty()
  {
    // if the signal has not been calculated for this point, update
    // uncertaintyFactor before returning it:
    if( needsToPrepareForThisPoint() )
      {
        rateCalculator->calculateValue( &signalValue,
                                        &uncertaintyFactor );
      }
    return log( uncertaintyFactor );
  }



  namespace signalClasses
  {
    inline bool
    reallyWrongCalculator::calculateValue( double* signalValue,
                                           double* uncertaintyFactor )
    // this always returns false, & always sets signalValue &
    // uncertaintyFactor to CppSLHA::CppSLHA_global::really_wrong_value.
    {
      *signalValue = CppSLHA::CppSLHA_global::really_wrong_value;
      *uncertaintyFactor = CppSLHA::CppSLHA_global::really_wrong_value;
      return false;
    }

    inline double
    reallyWrongCalculator::valueForCurrentCascades( fullCascade* firstCascade,
                                                   fullCascade* secondCascade )
    // this shouldn't be called.
    {
      std::cout
      << std::endl
      << "LHC-FASER::error! if this message is being displayed, somehow"
      << "reallyWrongCalculator::valueForCurrentCascades( " << firstCascade
      << ", " << secondCascade << " ) has been called.";
      std::cout << std::endl;

      return 0.0;
    }



    inline signalCalculator*
    sigmaBreakdownTest::getCalculator( std::string const& argumentString,
                                 signalDefinitionSet* const signalDefinitions )
    /* this either returns a pointer to a new sigmaBreakdownTest with cuts
     * decided by the given string, or NULL if the string could not be parsed
     * properly.
     */
    {
      if( 0 == argumentString.compare( 0,
                                       18,
                                       "sigmaBreakdownTest" ) )
      {
        signalDefinitions->setExcludedStandardModelProducts(
                          signalDefinitions->getShortcuts()->getInputShortcuts(
                                                          )->getNotInJets5() );
        return new sigmaBreakdownTest( signalDefinitions );
      }
      else
      {
        return NULL;
      }
    }

    inline double
    sigmaBreakdownTest::valueForCurrentCascades( fullCascade* firstCascade,
                                                 fullCascade* secondCascade )
    // this shouldn't be called.
    {
      std::cout
      << std::endl
      << "LHC-FASER::error! if this message is being displayed, somehow"
      << "sigmaBreakdownTest::valueForCurrentCascades( " << firstCascade
      << ", " << secondCascade << " ) has been called.";
      std::cout << std::endl;

      return 0.0;
    }



    inline double
    noJetMetButGivenLeptonCuts::valueForCurrentCascades(
                                                     fullCascade* firstCascade,
                                                   fullCascade* secondCascade )
    // see base version's description.
    {
      return ( subchannelCrossSectionTimesBrToEwinos
               * leptonAcceptanceCalculator->withAtLeastNJets(
                                                             signalDefinitions,
                                                               0,
                                                               firstCascade,
                                               firstSparticleIsNotAntiparticle,
                                                               secondCascade,
                                          secondSparticleIsNotAntiparticle ) );
    }

  }  // end of signalClasses namespace

}  // end of LHC_FASER namespace.

#endif /* LHC_FASER_SIGNAL_CALCULATOR_STUFF_HPP_ */
