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
 * - make sure that exclusive_BR_handler::always_neglected_decays doesn't
 *   interfere with the new signal. exclusive_BR_handler returns a NULL pointer
 *   if the decay is in always_neglected_decays.
 * - remember to include antiparticles in jets_to_tag & jets_to_neglect.
 * - fire hot, hot hurt.
 */

namespace LHC_FASER
{

  /* this is an abstract base class to do the actual calculation of the value
   * of the event rate for a signal. the main differences in derived classes
   * are what values they access from the handlers, how they put them
   * together, & how they estimate the uncertainty factor.
   */
  class signalCalculator
  {
  public:
    signalCalculator( signalDefinitionSet* signalDefinitions )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~signalCalculator()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    /* actually, this doesn't work, for various reasons. currently C++ does not
     * allow virtual static functions. I was only going to use it to force me
     * to remember to write the appropriate static functions in derived classes
     * anyway.
     */
    /*virtual static
    signalCalculator*
    getCalculator( std::string const* const arguments,
                   signalDefinitionSet* const signalDefinitions )
    = 0;*/
    /* this should be over-written in each derived class to construct a new
     * instance of the derived class based on arguments, & return a pointer to
     * the new instance of the derived class.
     */

    virtual bool
    calculate( double* signalValue,
               double* uncertaintyFactor )
    = 0;
    /* this calculates the event rate for the signal & puts its value in
     * signal_value, & puts an estimate for the uncertainty into
     * uncertainty_factor, & returns true if it did all this successfully.
     *
     *
     * the flow for signal_calculator::calculate() is:
     * sets up a set of channels (once per calculator instance)
     * - each of these channels associates with a cross-section table
     *
     * once per new point, each channel:
     * - checks its cross-section - if it's high enough, it proceeds
     * - obtains the sQCD_to_EWinos for its scoloreds
     * - checks the scew BRs - if they're high enough, it proceeds
     * - obtains its jet acceptance from its kinematics table, checks
     *   cross-section * scew BRs * jet acceptance, if high enough, proceeds
     * - obtains lepton acceptance from kinematics table, passes it with scew
     *   pairs to cascade handler, now has lepton & additional jet acceptances
     * - puts it all together
     * then all the channel totals are added together
     */

  protected:
    signalDefinitionSet signalDefinitions;
    std::vector< productionChannelPointerSet* > channels;
    //full_cascade* first_cascade;
    //full_cascade* second_cascade;
    std::vector< fullCascade* >* firstCascades;
    std::vector< fullCascade* >* secondCascades;
    std::list< int > excludedFinalStateParticles;
    double firstCascadeBrToEwino;
    double secondCascadeBrToEwino;
  };  // end of signal_calculator class.


  /* this is a class to handle each individual signal to be calculated.
   * it accesses numbers common to different signals through the
   * crossSectionHandler, kinematics_handler & cascade_handler classes.
   * it takes a string encoding what it should calculate, & constructs a
   * signalCalculator object to actually calculate the signal with the
   * given handlers. it also takes care of updating the signal & storing its
   * value & estimated uncertainty.
   */
  class signalHandler : public readied_for_new_point
  {
  public:
    signalHandler( std::string const signalName,
                   signal_shortcuts* const shortcut,
                   double const crossSectionUnitFactor )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~signalHandler()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    std::string const*
    getName()
    const
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getValue()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getUpperUncertainty()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getLowerUncertainty()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getUpperValue()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getLowerValue()
    /* code after the classes in this .hpp file, or in the .cpp file. */;
    double
    getLogUncertainty()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    std::string signalName;
    signalCalculator* rateCalculator;
    signal_shortcuts* const shortcut;
    double signalValue;
    double uncertaintyFactor;
    double crossSectionUnitFactor;
    // this is to allow for the user to specify event rates in fb, pb or nb.
  };  // end of signal_handler class.


  /* this is a derived class which just always returns
   * CppSLHA::CppSLHA_global::really_wrong_value as its value. it is used for
   * malformed signal names.
   */
  class reallyWrongCalculator : public signalCalculator
  {

  public:

    reallyWrongCalculator( signal_shortcuts* const shortcut )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    virtual
    ~reallyWrongCalculator()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    bool
    calculate( double* signalValue,
               double* uncertaintyFactor )
    // this always returns false, & always sets signal_value &
    // uncertainty_factor to CppSLHA::CppSLHA_global::really_wrong_value.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  //protected:
    //nothing.
  };

  /* this is a derived class to test whether I break down the production
   * channels & recombine them correctly. all acceptances are set to 1.0 so
   * that I should just recover the total cross-section.
   */
  class sigmaBreakdownTestCalculator : public signalCalculator
  {
  public:
    sigmaBreakdownTestCalculator(
                                 signalDefinitionSet* const signalDefinitions )
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~sigmaBreakdownTestCalculator()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    bool
    calculate( double* const signalValue,
               double* const uncertaintyFactor )
    /* this calculates the event rate for the signal & puts its value in
     * signal_value, & puts an estimate for the uncertainty into
     * uncertainty_factor, & returns true if it did all this successfully.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    double channelBrTotal;
  };  // end of sigmaBreakdownTestCalculator class.

  /* this is a derived class to calculate the "Atlas 4 jets plus missing
   * transverse momentum plus 0 leptons" signal.
   * it takes the kinematics from the Atlas4jMET grid & combines them with
   * cascade decays leading to 0 leptons passing the cut.
   */
  class Atlas4jMET0l_calculator : public signalCalculator
  {

  public:
    static int const jetAcceptanceGridTableColumn = 4;
    // this is dependent on the format of the grids.
    static double const defaultExtraJetCut = 40.0;
    // this is the standard cut for the jets beyond the hardest cut for this
    // signal as implemented in this code.
    static double const defaultLeptonCut = 10.0;
    // the default Atlas4jMET0l lepton transverse momentum cut is 10.0 GeV.

    static signalCalculator*
    getCalculator( std::string const* const arguments,
                   signalDefinitionSet* const signalDefinitions )
    // this either returns a pointer to a new Atlas4jMET0l_calculator with cuts
    // decided by the given string, or a pointer to a really_wrong_calculator.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~Atlas4jMET0l_calculator()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    bool
    calculate( double* const signalValue,
               double* const uncertaintyFactor )
    /* this calculates the event rate for the signal & puts its value in
     * signal_value, & puts an estimate for the uncertainty into
     * uncertainty_factor, & returns true if it did all this successfully.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;

  protected:
    kinematics_table* fourJetKinematics;
    kinematics_table* threeJetKinematics;
    /* the 4-jet signal is complicated enough that we also consider only 3 of
     * the 4 required hard jets coming from (showered) decays to
     * electroweakinos, with the 4th coming from a cascade decay of the
     * electroweakinos. hence we need an extra kinematics set (though it is
     * only used for direct jet acceptance).
     */
    double fourJetAcceptance;
    //double threeJetAcceptance;

    // these are used as each pairing of cascades from each production channel
    // is calculated.
    double subchannelCrossSectionTimesBrToEwinos;
    double subchannelValue;
    double subchannelZeroOrMoreJetsZeroLeptons;
    double subchannelOneOrMoreJetsZeroLeptons;

    Atlas4jMET0l_calculator( signalDefinitionSet* const signalDefinitions )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };  // end of Atlas4jMET0l_calculator class.

  /* this is a derived class to calculate the "Atlas 3 jets plus missing
   * transverse momentum plus exactly 1 lepton" signal.
   * it takes the kinematics from the Atlas3jMET grid & combines them with
   * cascade decays leading to 1 lepton passing the cut.
   */
  class Atlas3jMET1l_calculator : public signalCalculator
  {

  public:
    static int const jetAcceptanceGridTableColumn = 4;
    // this is dependent on the format of the grids.
    static double const defaultExtraJetCut = 40.0;
    // this is the standard cut for the jets beyond the hardest cut for this
    // signal as implemented in this code.
    static double const defaultPrimaryLeptonCut = 20.0;
    static double const defaultSecondaryLeptonCut = 10.0;
    // the default Atlas4jMET0l lepton transverse momentum cuts are 20.0 GeV
    // (for a single lepton to *pass*) & 10.0 GeV (for all others to *fail*).

    static signalCalculator*
    getCalculator( std::string const* const arguments,
                   signalDefinitionSet* const signalDefiner )
    // this either returns a pointer to a new Atlas3jMET1l_calculator with cuts
    // decided by the given string, or a pointer to a really_wrong_calculator.
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    ~Atlas3jMET1l_calculator()
    /* code after the classes in this .hpp file, or in the .cpp file. */;

    bool
    calculate( double* const signalValue,
               double* const uncertaintyFactor )
    /* this calculates the event rate for the signal & puts its value in
     * signal_value, & puts an estimate for the uncertainty into
     * uncertainty_factor, & returns true if it did all this successfully.
     */
    /* code after the classes in this .hpp file, or in the .cpp file. */;


  protected:
    // these are used as each pairing of cascades from each production channel
    // is calculated.
    double subchannelCrossSectionTimesBrToEwinos;
    double subchannelZeroOrMoreJetsOneLepton;

    Atlas3jMET1l_calculator( signalDefinitionSet* const signalDefinitions )
    /* code after the classes in this .hpp file, or in the .cpp file. */;
  };  // end of Atlas3jMET1l_calculator class.





  // inline functions:

  inline std::string const*
  signalHandler::getName()
  const
  {
   return &signalName;
  }

  inline double
  signalHandler::getValue()
  {
    // debugging:
    /**std::cout
    << std::endl
    << "debugging: signalHandler::getValue() called."
    << " needs_to_prepare_for_this_point() = "
    << needs_to_prepare_for_this_point();
    std::cout << std::endl;**/

    // if the signal has not been calculated for this point, update
    // signal_value before returning it:
    if( needs_to_prepare_for_this_point() )
      {
        rateCalculator->calculate( &signalValue,
                                   &uncertaintyFactor );
        signalValue *= crossSectionUnitFactor;
        finish_preparing_for_this_point();
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
    // uncertainty_factor before returning it:
    if( needs_to_prepare_for_this_point() )
      {
        rateCalculator->calculate( &signalValue,
                                   &uncertaintyFactor );
        finish_preparing_for_this_point();
      }
    return log( uncertaintyFactor );
  }



  inline bool
  reallyWrongCalculator::calculate( double* signalValue,
                                    double* uncertaintyFactor )
  // this always returns false, & always sets signal_value &
  // uncertainty_factor to CppSLHA::CppSLHA_global::really_wrong_value.
  {
    // debugging:
    /**std::cout
    << std::endl
    << "debugging: reallyWrongCalculator::calculate( "
    << signalValue << ", " << uncertaintyFactor << " ) called.";
    std::cout << std::endl;**/
    *signalValue = CppSLHA::CppSLHA_global::really_wrong_value;
    *uncertaintyFactor = CppSLHA::CppSLHA_global::really_wrong_value;
    return false;
  }
}  // end of LHC_FASER namespace.

#endif /* LHC_FASER_SIGNAL_CALCULATOR_STUFF_HPP_ */
