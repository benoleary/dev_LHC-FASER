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
  /* this is an abstract base class to do the actual calculation of the value
   * of the event rate for a signal. the main differences in derived classes
   * are what acceptances they access from the handlers, how they put them
   * together, & how they estimate the uncertainty factor.
   */
  class signalCalculator
  {
  public:
    signalCalculator( signalDefinitionSet* signalDefinitions );
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
                   signalDefinitionSet* const signalDefinitions )
    = 0;*/
    /* this should be over-written in each derived class to construct a new
     * instance of the derived class based on arguments, & return a pointer to
     * the new instance of the derived class.
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
    static int
    parseNumberOfLeptons( std::string const& argumentString );
    // this looks for "_0l", "_1l", or "_2l", & returns the appropriate int.
    // returns -1 if it could not find any.
    static bool
    parseBeamEnergy( std::string const& argumentString,
                     signalDefinitionSet* const signalDefinitions );
    // this looks for "_7TeV", "_07TeV", "_10TeV", or "_14TeV", & sets the beam
    // energy appropriately. returns false if it could not find any.
    static void
    parseLeptonTransverseMomentumCuts( std::string const& argumentString,
                                signalDefinitionSet* const signalDefinitions,
                                       double const defaultPrimaryLeptonCut,
                                      double const defaultSecondaryLeptonCut );
    /* this looks for "_pTl" then a double then "GeV", then a subsequent
     * "_pTl", a subsequent double, & a subsequent "GeV", interpretted as the
     * primary & secondary lepton cuts respectively.
     */

    signalDefinitionSet signalDefinitions;
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
    signalHandler( std::string const signalName,
                   double const crossSectionUnitFactor,
                   signalDefinitionSet* const signalPreparationDefinitions );
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
    signalShortcuts* const inputShortcut;
    double signalValue;
    double uncertaintyFactor;
    double crossSectionUnitFactor;
    // this is to allow for the user to specify event rates in fb, pb or nb.
  };


  namespace signalCalculatorClasses
  {
    /* this is a derived class which just always returns
     * CppSLHA::CppSLHA_global::really_wrong_value as its value. it is used for
     * malformed signal names.
     */
    class reallyWrongCalculator : public signalCalculator
    {
    public:
      reallyWrongCalculator( signalDefinitionSet* const signalDefinitions );
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
      sigmaBreakdownTest( signalDefinitionSet* const signalDefinitions );
      ~sigmaBreakdownTest();

      bool
      calculateValue( double* const signalValue,
                      double* const uncertaintyFactor );
      // this over-rides the base version because it does things a bit
      // differently in a few of the steps.


    protected:
      double channelBrTotal;

      virtual double
      valueForCurrentCascades( fullCascade* firstCascade,
                               fullCascade* secondCascade );
      // this shouldn't be called.
    };


    class leptonAcceptanceForCascadePair
    {
    public:
      leptonAcceptanceForCascadePair(
                                signalDefinitionSet* const signalDefinitions );
      virtual
      ~leptonAcceptanceForCascadePair();

      virtual double
      withExactlyNJets( int const exactNumberOfJets,
                        fullCascade* firstCascade,
                        bool firstIsNotAntiparticle,
                        fullCascade* secondCascade,
                        bool secondIsNotAntiparticle ) = 0;
      virtual double
      withAtLeastNJets( int const minimumNumberOfJets,
                        fullCascade* firstCascade,
                        bool firstIsNotAntiparticle,
                        fullCascade* secondCascade,
                        bool secondIsNotAntiparticle );


    protected:
      static int const maximumNumberOfJets;

      signalDefinitionSet* const signalDefinitions;
    };


    namespace leptonAcceptanceStyle
    {
      class fullySpecified : public leptonAcceptanceForCascadePair
      {
      public:
        fullySpecified( signalDefinitionSet* const signalDefinitions,
                        int const numberOfNegativeElectrons,
                        int const numberOfPositiveElectrons,
                        int const numberOfNegativeMuons,
                        int const numberOfPositiveMuons );
        virtual
        ~fullySpecified();

        virtual double
        withExactlyNJets( int const exactNumberOfJets,
                          fullCascade* firstCascade,
                          bool firstIsNotAntiparticle,
                          fullCascade* secondCascade,
                          bool secondIsNotAntiparticle );
        bool
        isSameAs( fullySpecified const* const comparisonPointer );


      protected:
        int const numberOfNegativeElectrons;
        int const numberOfPositiveElectrons;
        int const numberOfNegativeMuons;
        int const numberOfPositiveMuons;
      };


      class noLeptonCutNorExtraJetCut : public leptonAcceptanceForCascadePair
      {
      public:
        noLeptonCutNorExtraJetCut(
                                signalDefinitionSet* const signalDefinitions );
        virtual
        ~noLeptonCutNorExtraJetCut();

        virtual double
        withExactlyNJets( int const exactNumberOfJets,
                          fullCascade* firstCascade,
                          bool firstIsNotAntiparticle,
                          fullCascade* secondCascade,
                          bool secondIsNotAntiparticle );
        virtual double
        withAtLeastNJets( int const minimumNumberOfJets,
                          fullCascade* firstCascade,
                          bool firstIsNotAntiparticle,
                          fullCascade* secondCascade,
                          bool secondIsNotAntiparticle );


      //protected:
        // nothing.
      };


      class chargeAndFlavorSummed : public leptonAcceptanceForCascadePair
      {
      public:
        chargeAndFlavorSummed( signalDefinitionSet* const signalDefinitions,
                               int const numberOfLeptons );
        virtual
        ~chargeAndFlavorSummed();

        virtual double
        withAtLeastNJets( int const minimumNumberOfJets,
                          fullCascade* firstCascade,
                          bool firstIsNotAntiparticle,
                          fullCascade* secondCascade,
                          bool secondIsNotAntiparticle );
        virtual double
        withExactlyNJets( int const exactNumberOfJets,
                          fullCascade* firstCascade,
                          bool firstIsNotAntiparticle,
                          fullCascade* secondCascade,
                          bool secondIsNotAntiparticle );
        bool
        isSameAs( chargeAndFlavorSummed const* const comparisonPointer );


      protected:
        int const numberOfLeptons;
      };


      class chargeSummed : public leptonAcceptanceForCascadePair
      {
      public:
        chargeSummed( signalDefinitionSet* const signalDefinitions,
                      int const numberOfElectrons,
                      int const numberOfMuons );
        virtual
        ~chargeSummed();

        virtual double
        withAtLeastNJets( int const minimumNumberOfJets,
                          fullCascade* firstCascade,
                          bool firstIsNotAntiparticle,
                          fullCascade* secondCascade,
                          bool secondIsNotAntiparticle );
        virtual double
        withExactlyNJets( int const exactNumberOfJets,
                          fullCascade* firstCascade,
                          bool firstIsNotAntiparticle,
                          fullCascade* secondCascade,
                          bool secondIsNotAntiparticle );
        bool
        isSameAs( chargeSummed const* const comparisonPointer );


      protected:
        int const numberOfElectrons;
        int const numberOfMuons;
      };


      class flavorSummed : public leptonAcceptanceForCascadePair
      {
      public:
        flavorSummed( signalDefinitionSet* const signalDefinitions,
                      int const numberOfNegativeLeptons,
                      int const numberOfPositiveLeptons );
        virtual
        ~flavorSummed();

        virtual double
        withAtLeastNJets( int const minimumNumberOfJets,
                          fullCascade* firstCascade,
                          bool firstIsNotAntiparticle,
                          fullCascade* secondCascade,
                          bool secondIsNotAntiparticle );
        virtual double
        withExactlyNJets( int const exactNumberOfJets,
                          fullCascade* firstCascade,
                          bool firstIsNotAntiparticle,
                          fullCascade* secondCascade,
                          bool secondIsNotAntiparticle );
        bool
        isSameAs( flavorSummed const* const comparisonPointer );


      protected:
        int const numberOfNegativeLeptons;
        int const numberOfPositiveLeptons;
      };


      class ossfMinusOsdf : public leptonAcceptanceForCascadePair
      {
      public:
        ossfMinusOsdf( signalDefinitionSet* const signalDefinitions,
                       bool const restrictedToOnePair );
        virtual
        ~ossfMinusOsdf();

        virtual double
        withAtLeastNJets( int const minimumNumberOfJets,
                          fullCascade* firstCascade,
                          bool firstIsNotAntiparticle,
                          fullCascade* secondCascade,
                          bool secondIsNotAntiparticle );
        virtual double
        withExactlyNJets( int const exactNumberOfJets,
                          fullCascade* firstCascade,
                          bool firstIsNotAntiparticle,
                          fullCascade* secondCascade,
                          bool secondIsNotAntiparticle );
        bool
        isSameAs( ossfMinusOsdf const* const comparisonPointer );


      protected:
        bool const restrictedToOnePair;
      };

    }  // end of leptonAcceptanceStyle namespace


    class leptonAcceptanceForCascadePairFactory
    {
    public:
      leptonAcceptanceForCascadePairFactory();
      ~leptonAcceptanceForCascadePairFactory();

      leptonAcceptanceForCascadePair*
      getFullySpecified( signalDefinitionSet* const signalDefinitions,
                         int const numberOfNegativeElectrons,
                         int const numberOfPositiveElectrons,
                         int const numberOfNegativeMuons,
                         int const numberOfPositiveMuons );
      leptonAcceptanceForCascadePair*
      getNoLeptonCut( signalDefinitionSet* const signalDefinitions );
      leptonAcceptanceForCascadePair*
      getChargeAndFlavorSummed( signalDefinitionSet* const signalDefinitions,
                                int const numberOfLeptons );
      leptonAcceptanceForCascadePair*
      getChargeSummed( signalDefinitionSet* const signalDefinitions,
                       int const numberOfElectrons,
                       int const numberOfMuons );
      leptonAcceptanceForCascadePair*
      getFlavorSummed( signalDefinitionSet* const signalDefinitions,
                       int const numberOfNegativeLeptons,
                       int const numberOfPositiveLeptons );
      leptonAcceptanceForCascadePair*
      getOssfMinusOsdf( signalDefinitionSet* const signalDefinitions,
                        bool const restrictedToOnePair );


    protected:
      std::vector< leptonAcceptanceStyle::fullySpecified* > fullySpecifiedSet;
      leptonAcceptanceStyle::noLeptonCutNorExtraJetCut* noLeptonOrExtraJetCuts;
      std::vector< leptonAcceptanceStyle::chargeAndFlavorSummed* >
      chargeAndFlavorSummedSet;
      std::vector< leptonAcceptanceStyle::chargeSummed* > chargeSummedSet;
      std::vector< leptonAcceptanceStyle::flavorSummed* > flavorSummedSet;
      std::vector< leptonAcceptanceStyle::ossfMinusOsdf* > ossfMinusOsdfSet;
    };


    /* this is a derived class to calculateValue the "Atlas 4 jets plus missing
     * transverse momentum plus N leptons" signal.
     * it takes the kinematics from the Atlas4jMET grid & combines them with
     * cascade decays leading to N leptons passing the cut.
     */
    class atlasFourJetMetNLeptons : public signalCalculator
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
      getCalculator( std::string const* const argumentString,
                     signalDefinitionSet* const signalDefinitions );
      /* this either returns a pointer to a new atlasFourJetMetNLeptons with
       * cuts decided by the given string, or a pointer to a
       * reallyWrongCalculator.
       */

      ~atlasFourJetMetNLeptons();


    protected:
      int const numberOfLeptons;
      jetAcceptanceTable* fourJetKinematics;
      jetAcceptanceTable* threeJetKinematics;
      /* the 4-jet signal is complicated enough that we also consider only 3 of
       * the 4 required hard jets coming from (showered) decays to
       * electroweakinos, with the 4th coming from a cascade decay of the
       * electroweakinos. hence we need an extra kinematics set (though it is
       * only used for direct jet acceptance).
       */
      double fourJetAcceptance;
      // these are used as each pairing of cascades from each production
      // channel is calculated:
      double subchannelValue;
      double subchannelZeroOrMoreJetsNLeptons;
      double subchannelOneOrMoreJetsNLeptons;

      atlasFourJetMetNLeptons( signalDefinitionSet* const signalDefinitions,
                               int const numberOfLeptons );

      virtual double
      valueForCurrentCascades( fullCascade* firstCascade,
                               fullCascade* secondCascade );
      // see base version's description.
    };


    /* this is a derived class to calculate the "Atlas 3 jets plus missing
     * transverse momentum plus exactly 1 lepton" signal.
     * it takes the kinematics from the Atlas3jMET grid & combines them with
     * cascade decays leading to 1 lepton passing the cut.
     */
    class atlasThreeJetMetNLeptons : public signalCalculator
    {
    public:
      static int const jetAcceptanceGridTableColumn;
      // this is dependent on the format of the grids.
      static double const defaultExtraJetCut;
      // this is the standard cut for the jets beyond the hardest cut for this
      // signal as implemented in this code.
      static double const defaultPrimaryLeptonCut;
      static double const defaultSecondaryLeptonCut;
      // the default Atlas3jMET1l lepton transverse momentum cuts are 20.0 GeV
      // (for a single lepton to *pass*) & 10.0 GeV (for all others to *fail*).

      static signalCalculator*
      getCalculator( std::string const* const argumentString,
                     signalDefinitionSet* const signalDefinitions );
      // this either returns a pointer to a new atlasThreeJetMetNLeptons with
      // cuts decided by the given string, or a reallyWrongCalculator pointer.

      ~atlasThreeJetMetNLeptons();


    protected:
      int const numberOfLeptons;
      // this is used as each pairing of cascades from each production
      // channel is calculated.
      double subchannelZeroOrMoreJetsNLeptons;

      atlasThreeJetMetNLeptons( signalDefinitionSet* const signalDefinitions,
                                int const numberOfLeptons );

      virtual double
      valueForCurrentCascades( fullCascade* firstCascade,
                               fullCascade* secondCascade );
      // see base version's description.
    };


    /* this is a derived class to calculate  the "same sign dilepton" signal.
     * it assumes no jet or MET cuts, so combines acceptances just from lepton
     * cuts combined with cross-sections & branching ratios.
     */
    class sameSignDilepton : public signalCalculator
    {
    public:
      static double const defaultLeptonCut;
      // the default same-sign dilepton lepton transverse momentum cut is
      // 20.0 GeV.

      static signalCalculator*
      getCalculator( std::string const* const argumentString,
                     signalDefinitionSet* const signalDefinitions );
      // this either returns a pointer to a new sameSignDilepton with cuts
      // decided by the given string, or a reallyWrongCalculator pointer.

      ~sameSignDilepton();


    protected:
      // this is used as each pairing of cascades from each production
      // channel is calculated.
      double subchannelDileptonAcceptance;

      sameSignDilepton( signalDefinitionSet* const signalDefinitions );

      double
      getDileptonFromSingleCascade( fullCascade* dileptonCascade,
                                  bool const dileptonCascaderIsNotAntiparticle,
                                    fullCascade* otherCascade,
                                    bool const otherIsNotAntiparticle );
      double
      getOneLeptonFromEachCascade( fullCascade* firstCascade,
                                   bool const firstIsNotAntiparticle,
                                   fullCascade* secondCascade,
                                   bool const secondIsNotAntiparticle );
      virtual double
      valueForCurrentCascades( fullCascade* firstCascade,
                               fullCascade* secondCascade );
      // see base version's description.
    };

  }  // end of signalCalculatorClasses namespace.





  // inline functions:


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




  namespace signalCalculatorClasses
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
    leptonAcceptanceForCascadePair::withAtLeastNJets(
                                                 int const minimumNumberOfJets,
                                                     fullCascade* firstCascade,
                                                   bool firstIsNotAntiparticle,
                                                    fullCascade* secondCascade,
                                                 bool secondIsNotAntiparticle )
    {
      double returnValue( 0.0 );
      for( int numberOfJets( minimumNumberOfJets );
           maximumNumberOfJets >= numberOfJets;
           ++numberOfJets )
      {
        returnValue += withExactlyNJets( numberOfJets,
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
      noLeptonCutNorExtraJetCut::withExactlyNJets( int const exactNumberOfJets,
                                                   fullCascade* firstCascade,
                                                   bool firstIsNotAntiparticle,
                                                   fullCascade* secondCascade,
                                                 bool secondIsNotAntiparticle )
      {
        return 1.0;
      }

      inline double
      noLeptonCutNorExtraJetCut::withAtLeastNJets( int const exactNumberOfJets,
                                                   fullCascade* firstCascade,
                                                   bool firstIsNotAntiparticle,
                                                   fullCascade* secondCascade,
                                                 bool secondIsNotAntiparticle )
      {
        return 1.0;
      }


      class chargeAndFlavorSummed : public leptonAcceptanceForCascadePair
      {
      public:
        chargeAndFlavorSummed( signalDefinitionSet* const signalDefinitions,
                               int const numberOfLeptons );
        virtual
        ~chargeAndFlavorSummed();

        virtual double
        withAtLeastNJets( int const minimumNumberOfJets,
                          fullCascade* firstCascade,
                          fullCascade* secondCascade );
        virtual double
        withExactlyNJets( int const exactNumberOfJets,
                          fullCascade* firstCascade,
                          fullCascade* secondCascade );
        bool
        isSameAs( chargeAndFlavorSummed const* const comparisonPointer );


      protected:
        int const numberOfLeptons;
      };


      class chargeSummed : public leptonAcceptanceForCascadePair
      {
      public:
        chargeSummed( signalDefinitionSet* const signalDefinitions,
                      int const numberOfElectrons,
                      int const numberOfMuons );
        virtual
        ~chargeSummed();

        virtual double
        withAtLeastNJets( int const minimumNumberOfJets,
                          fullCascade* firstCascade,
                          fullCascade* secondCascade );
        virtual double
        withExactlyNJets( int const exactNumberOfJets,
                          fullCascade* firstCascade,
                          fullCascade* secondCascade );
        bool
        isSameAs( chargeSummed const* const comparisonPointer );


      protected:
        int const numberOfElectrons;
        int const numberOfMuons;
      };


      class flavorSummed : public leptonAcceptanceForCascadePair
      {
      public:
        flavorSummed( signalDefinitionSet* const signalDefinitions,
                      int const numberOfNegativeLeptons,
                      int const numberOfPositiveLeptons );
        virtual
        ~flavorSummed();

        virtual double
        withAtLeastNJets( int const minimumNumberOfJets,
                          fullCascade* firstCascade,
                          fullCascade* secondCascade );
        virtual double
        withExactlyNJets( int const exactNumberOfJets,
                          fullCascade* firstCascade,
                          fullCascade* secondCascade );
        bool
        isSameAs( flavorSummed const* const comparisonPointer );


      protected:
        int const numberOfNegativeLeptons;
        int const numberOfPositiveLeptons;
      };


      class ossfMinusOsdf : public leptonAcceptanceForCascadePair
      {
      public:
        ossfMinusOsdf( signalDefinitionSet* const signalDefinitions,
                       bool const restrictedToOnePair );
        virtual
        ~ossfMinusOsdf();

        virtual double
        withAtLeastNJets( int const minimumNumberOfJets,
                          fullCascade* firstCascade,
                          fullCascade* secondCascade );
        virtual double
        withExactlyNJets( int const exactNumberOfJets,
                          fullCascade* firstCascade,
                          fullCascade* secondCascade );
        bool
        isSameAs( ossfMinusOsdf const* const comparisonPointer );


      protected:
        bool const restrictedToOnePair;
      };

    }  // end of leptonAcceptanceStyle namespace


    inline double
    sameSignDilepton::getDileptonFromSingleCascade(
                                                  fullCascade* dileptonCascade,
                                  bool const dileptonCascaderIsNotAntiparticle,
                                                    fullCascade* otherCascade,
                                            bool const otherIsNotAntiparticle )
    {
      return ( ( dileptonCascade->leptonFlavorSummedWithMinimumJets(
                                             dileptonCascaderIsNotAntiparticle,
                                                            &signalDefinitions,
                                                                     0,
                                                                     2,
                                                                     0 )
                 + dileptonCascade->leptonFlavorSummedWithMinimumJets(
                                             dileptonCascaderIsNotAntiparticle,
                                                            &signalDefinitions,
                                                                       0,
                                                                       0,
                                                                       2 ) )
               * otherCascade->getAcceptanceWithMinimumJets(
                                                        otherIsNotAntiparticle,
                                                            &signalDefinitions,
                                                             0,
                                                             0,
                                                             0,
                                                             0,
                                                             0 ) );
    }

  }  // end of signalCalculatorClasses namespace

}  // end of LHC_FASER namespace.

#endif /* LHC_FASER_SIGNAL_CALCULATOR_STUFF_HPP_ */
