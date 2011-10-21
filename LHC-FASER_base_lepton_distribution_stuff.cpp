/*
 * LHC-FASER_base_lepton_distribution_stuff.cpp
 *
 *  Created on: 17 Mar 2010
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
 *      LHC-FASER_base_lepton_distribution_stuff.hpp
 *      LHC-FASER_base_lepton_distribution_stuff.cpp
 *      LHC-FASER_derived_lepton_distributions.hpp
 *      LHC-FASER_derived_lepton_distributions.cpp
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


#include "LHC-FASER_base_lepton_distribution_stuff.hpp"


namespace LHC_FASER
{
  effectiveSquarkMassHolder::effectiveSquarkMassHolder()
  {
    // does nothing.
  }

  effectiveSquarkMassHolder::~effectiveSquarkMassHolder()
  {
    // does nothing.
  }



  onshellSquarkMassHolder::onshellSquarkMassHolder(
                  CppSLHA::particle_property_set const* const onshellSquark ) :
    onshellSquark( onshellSquark )
  {
    // just an initialization list.
  }

  onshellSquarkMassHolder::~onshellSquarkMassHolder()
  {
    // does nothing.
  }



  leptonDistributionExpansionTerm::leptonDistributionExpansionTerm(
                                                       int const powerOfEnergy,
                                                    int const powerOfLogarithm,
                                                 double const coefficientValue,
                       leptonDistributionExpansionTerm* const referenceTerm ) :
    powerOfEnergy( powerOfEnergy ),
    powerOfLogarithm( powerOfLogarithm ),
    coefficientValue( coefficientValue ),
    referenceTerm( referenceTerm )
  {
    // just an initialization list.
  }


  leptonDistributionExpansionTerm::~leptonDistributionExpansionTerm()
  {
    // does nothing.
  }


  double
  leptonDistributionExpansionTerm::evaluate( double const inputEnergy )
  const
  // this evaluates this term for the given input energy.
  {
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "leptonDistributionExpansionTerm::evaluate( " << inputEnergy
    << " ) called. coefficientValue = " << coefficientValue;
    std::cout << std::endl;**/

    if( 0.0 != coefficientValue )
    {
      double returnValue( coefficientValue );

      for( int powerCount( powerOfEnergy );
           0 < powerCount;
           --powerCount )
      {
        returnValue *= inputEnergy;
      }

      double logarithmOfEnergy( log( inputEnergy ) );
      for( int logarithmCounter( powerOfLogarithm );
           0 < logarithmCounter;
           --logarithmCounter )
      {
        returnValue *= logarithmOfEnergy;
      }

      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "returning " << returnValue;
      std::cout << std::endl;**/
      return returnValue;
    }
    else
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "returning 0.0";
      std::cout << std::endl;**/
      return 0.0;
    }
  }

  double
  leptonDistributionExpansionTerm::getArea( double startEnergy,
                                            double endEnergy )
  const
  // this gives the definite integral of the term from startEnergy to
  // endEnergy.
  {
    double returnValue( 0.0 );

    // if this version of the virtual function is being used,
    // powerOfEnergy >= 0 & powerOfLogarithm >= 1.
    if( 0.0 != coefficientValue )
    {
      if( minimumInputEnergy > startEnergy )
      {
        startEnergy = minimumInputEnergy;
      }

      if( endEnergy > startEnergy )
      {
        int const energyPowerPlusOne( ( 1 + powerOfEnergy ) );
        int predivisorOfPowersOfEnergyPowerPlusOne( -energyPowerPlusOne
                                                    * energyPowerPlusOne );
        double const minusEnergyPowerMinusOne( -(double)energyPowerPlusOne );
        // the terms in the sum alternate sign, so I put it here.
        double integralToStart( coefficientValue );
        double integralToEnd( coefficientValue );

        for( int powerCount( powerOfEnergy );
             0 <= powerCount;
             --powerCount )
        {
          integralToStart *= startEnergy;
          integralToEnd *= endEnergy;
        }
        /* now integralToStart is
         * coefficientValue * startEnergy^( 1 + powerOfEnergy )
         * & similar for integralToEnd.
         */

        double logarithmOfStartEnergy( log( startEnergy ) );
        double logarithmOfEndEnergy( log( endEnergy ) );

        int logPowerFactorial( 1 );
        double commonSumFactor;
        double
        startSumTerm( minusEnergyPowerMinusOne * logarithmOfStartEnergy );
        double endSumTerm( minusEnergyPowerMinusOne * logarithmOfEndEnergy );
        double startSum( 1.0 + startSumTerm );
        double endSum( 1.0 + endSumTerm );
        /* now startSum is the sum of the terms corresponding to the lowest 2
         * powers of ln( startEnergy ) (i.e. 0 & 1, so the constant term & the
         * term linear in ln( startEnergy )), similar for endSum.
         * these need to be summed with the terms with powers of
         * ln( startEnergy ) from 2 to powerOfLogarithm.
         * logPowerFactorial = 1!, &
         * predivisorOfPowersOfEnergyPowerPlusOne already has its 1st 2 factors
         * of energyPowerPlusOneAsInt. logPowerFactorial needs multiplying by
         * the integers 2 to powerOfLogarithm (to end up as powerOfLogarithm!)
         * & predivisorOfPowersOfEnergyPowerPlusOne needs to be multiplied by
         * energyPowerPlusOneAsDouble ( powerOfLogarithm - 1 ) more times (to
         * end up as energyPowerPlusOneAsDouble^powerOfLogarithm).
         * the following loop goes over the integers 2 to powerOfLogarithm,
         * hence ( powerOfLogarithm - 1 ) iterations.
         */
        for( int logarithmCounter( 2 );
             powerOfLogarithm >= logarithmCounter;
             ++logarithmCounter )
        {
          logPowerFactorial *= logarithmCounter;
          predivisorOfPowersOfEnergyPowerPlusOne *= -energyPowerPlusOne;
          commonSumFactor
          = ( minusEnergyPowerMinusOne / (double)logarithmCounter );

          // debugging:
          /**std::cout << std::endl << "debugging:"
          << std::endl
          << "logarithmCounter = " << logarithmCounter
          << std::endl
          << "predivisorOfPowersOfEnergyPowerPlusOne = "
          << predivisorOfPowersOfEnergyPowerPlusOne
          << std::endl
          << "commonSumFactor: " << commonSumFactor
          << " * logarithmOfStartEnergy: " << logarithmOfStartEnergy
          << " * startSumTerm: " << startSumTerm;
          std::cout << std::endl;**/

          startSumTerm *= ( logarithmOfStartEnergy * commonSumFactor );
          startSum += startSumTerm;
          endSumTerm *= ( logarithmOfEndEnergy * commonSumFactor );
          endSum += endSumTerm;

          // debugging:
          /**std::cout
          << " = startSumTerm = " << startSumTerm;
          std::cout << std::endl;**/
        }
        integralToStart *= startSum;
        integralToEnd *= endSum;

        returnValue
        = ( ( (double)logPowerFactorial * ( integralToEnd - integralToStart ) )
            / (double)predivisorOfPowersOfEnergyPowerPlusOne );
      }
    }
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "leptonDistributionExpansionTerm::getArea( " << startEnergy << ", "
    << endEnergy << " ) for " << coefficientValue << " E^" << powerOfEnergy
    << " L^" << powerOfLogarithm << " returning " << returnValue;
    std::cout << std::endl;**/

    return returnValue;
  }

  double const leptonDistributionExpansionTerm::minimumInputEnergy( 0.001 );
   // 1 MeV should be sufficient as the energy resolution.



  leptonDistributionInverseTerm::leptonDistributionInverseTerm(
                                                 double const coefficientValue,
                       leptonDistributionExpansionTerm* const referenceTerm ) :
    leptonDistributionExpansionTerm( -1,
                                     0,
                                     coefficientValue,
                                     referenceTerm )
  {
    // just an initialization.
  }

  leptonDistributionInverseTerm::~leptonDistributionInverseTerm()
  {
    // does nothing.
  }



  leptonDistributionConstantTerm::leptonDistributionConstantTerm(
                                                 double const coefficientValue,
                       leptonDistributionExpansionTerm* const referenceTerm ) :
    leptonDistributionExpansionTerm( 0,
                                     0,
                                     coefficientValue,
                                     referenceTerm )
  {
    // just an initialization.
  }

  leptonDistributionConstantTerm::~leptonDistributionConstantTerm()
  {
    // does nothing.
  }


  leptonDistributionPowerTerm::leptonDistributionPowerTerm(
                                                       int const powerOfEnergy,
                                                 double const coefficientValue,
                       leptonDistributionExpansionTerm* const referenceTerm ) :
    leptonDistributionExpansionTerm( powerOfEnergy,
                                     0,
                                     coefficientValue,
                                     referenceTerm )
  {
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "leptonDistributionPowerTerm::leptonDistributionPowerTerm( "
    << powerOfEnergy << ", " << coefficientValue << ", " << referenceTerm
    << " ) called.";
    std::cout << std::endl;**/

    // just an initialization.
  }

  leptonDistributionPowerTerm::~leptonDistributionPowerTerm()
  {
    // does nothing.
  }


  double
  leptonDistributionPowerTerm::getArea( double startEnergy,
                                        double endEnergy )
  const
  // this gives the definite integral of the term from startEnergy to
  // endEnergy.
  {
    double returnValue( 0.0 );

    // if this version of the virtual function is being used,
    // powerOfEnergy >= 1.
    if( 0.0 != coefficientValue )
    {
      if( minimumInputEnergy > startEnergy )
      {
        startEnergy = minimumInputEnergy;
      }

      if( endEnergy > startEnergy )
      {
        double energyPowerPlusOne( ( 1.0 + (double)powerOfEnergy ) );
        double integralToStart( startEnergy );
        double integralToEnd( endEnergy );

        for( int powerCount( powerOfEnergy );
             0 < powerCount;
             --powerCount )
        {
          integralToStart *= startEnergy;
          integralToEnd *= endEnergy;
        }
        /* now integralToStart is startEnergy^( 1 + powerOfEnergy )
         * & similar for integralToEnd.
         */
        returnValue
        = ( ( ( integralToEnd - integralToStart ) * coefficientValue )
            / energyPowerPlusOne );
      }
    }
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "leptonDistributionPowerTerm::getArea( " << startEnergy << ", "
    << endEnergy << " ) returning " << returnValue;
    std::cout << std::endl;**/
    return returnValue;
  }



  segmentTermSet::segmentTermSet() :
    segmentStart( CppSLHA::CppSLHA_global::really_wrong_value ),
    segmentEnd( CppSLHA::CppSLHA_global::really_wrong_value )
  {
    // just an initialization list.
  }

  segmentTermSet::~segmentTermSet()
  {
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "segmentTerms.size() = " << segmentTerms.size();
    std::cout << std::endl;**/

    for( std::vector< leptonDistributionExpansionTerm* >::iterator
         deletionIterator( segmentTerms.begin() );
         segmentTerms.end() > deletionIterator;
         ++deletionIterator )
    {
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "~segmentTermSet() deleting a leptonDistributionExpansionTerm "
      << *deletionIterator << ", "
      << (*deletionIterator)->getCoefficient()
      << " E^" << (*deletionIterator)->getPowerOfEnergy()
      << " L^" << (*deletionIterator)->getPowerOfLogarithm();
      std::cout << std::endl;**/
      delete *deletionIterator;
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "successfully deleted.";
      std::cout << std::endl;**/
    }
  }


  leptonDistributionExpansionTerm*
  segmentTermSet::addTerm( int const powerOfEnergy,
                           int const powerOfLogarithm,
                          leptonDistributionExpansionTerm* const referenceTerm,
                           double const coefficientValue )
  {
    if( 0 == powerOfLogarithm )
    {
      if( -1 == powerOfEnergy )
      {
        termAdder = new leptonDistributionInverseTerm( coefficientValue,
                                                       referenceTerm );
      }
      else if( 0 == powerOfEnergy )
      {
        termAdder = new leptonDistributionConstantTerm( coefficientValue,
                                                        referenceTerm );
      }
      else
      {
        termAdder = new leptonDistributionPowerTerm( powerOfEnergy,
                                                     coefficientValue,
                                                     referenceTerm );
      }
    }
    else
    {
      termAdder = new leptonDistributionExpansionTerm( powerOfEnergy,
                                                       powerOfLogarithm,
                                                       coefficientValue,
                                                       referenceTerm );
    }
    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "segmentTermSet::addTerm adding " << termAdder;
    std::cout << std::endl;**/

    segmentTerms.push_back( termAdder );
    return termAdder;
  }

  double
  segmentTermSet::evaluate( double const inputEnergy )
  const
  // this evaluates the sum of the set of terms for the given input energy.
  {
    // debugging:
    /**std::cout << std::endl << "debugging: "
    << "segmentTermSet::evaluate( " << inputEnergy
    << " ) called";
    std::cout << std::endl;**/

    double returnValue( 0.0 );

    // only evaluate if the input is within the segment's range:
    if( ( inputEnergy >= segmentStart )
        &&
        ( 0.0 < inputEnergy )
        &&
        ( inputEnergy <= segmentEnd ) )
      {
        // debugging:
        /**int segmentTermCounter( 0 );
        std::cout << std::endl << "debugging: "
        << "inputEnergy within range, segmentTerms.size() = "
        << segmentTerms.size();
        std::cout << std::endl;**/

        for( std::vector< leptonDistributionExpansionTerm* >::const_iterator
             summationIterator = segmentTerms.begin();
             summationIterator < segmentTerms.end();
             ++summationIterator )
          {
            // debugging:
            /**std::cout << std::endl << "debugging: "
            << "segmentTermCounter = " << ++segmentTermCounter;
            std::cout << std::endl;**/

            returnValue += (*summationIterator)->evaluate( inputEnergy );
          }
      }
    return returnValue;
  }



  tauDecayCoefficient::tauDecayCoefficient()
  {
    // just an initialization.
  }

  tauDecayCoefficient::~tauDecayCoefficient()
  {
    // does nothing.
  }



  hardMuonFromTau::hardMuonFromTau() :
    tauDecayCoefficient()
  {
    // just an initialization.
  }

  hardMuonFromTau::~hardMuonFromTau()
  {
    // does nothing.
  }


  double
  hardMuonFromTau::operator()( int const visibleProductEnergyPower
                            /* the power of the muon's energy for the term. */,
                               int const visibleProductLogPower
           /* the power of the logarithm of the muon's energy for the term. */,
                               int const tauLeptonEnergyPower
                      /* the power of the tau lepton's energy for the term. */,
                               int const tauLeptonLogPower
     /* the power of the logarithm of the tau lepton's energy for the term. */,
                               double tauMinEnergy
         /* the minimum of the range of the tau lepton energies considered. */,
                               double const tauMaxEnergy
         /* the maximum of the range of the tau lepton energies considered. */,
                               bool const isInsideRange )
  const
  /* this returns the coefficientValue for the muon distribution which would come
   * from a left-handed tau lepton distribution of the given power of the tau
   * lepton's energy & power of logarithm thereof, for requested integer powers
   * of the muon's energy or logarithm thereof, which also depends on whether
   * the muon's energy is inside the range of the tau lepton's energy for this
   * segment.
   */
  {
    // debugging:
    /**std::cout << std::endl << "debugging: "
    << "hardMuonFromTau::operator()( " << visibleProductEnergyPower << ", "
    << visibleProductLogPower << ", "
    << tauLeptonEnergyPower << ", "
    << tauLeptonLogPower << ", "
    << tauMinEnergy << ", "
    << tauMaxEnergy << ", "
    << isInsideRange
    << " ) called";
    std::cout << std::endl;**/

    if( leptonDistributionExpansionTerm::minimumInputEnergy > tauMinEnergy )
    {
      tauMinEnergy = leptonDistributionExpansionTerm::minimumInputEnergy;
    }

    double returnValue( 0.0 );

    if( ( -1 == tauLeptonEnergyPower )
        &&
        ( 0 == tauLeptonLogPower ) )
    {
      if( ( -1 == visibleProductEnergyPower )
          &&
          ( 0 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = 0.75;
        }
      }
      else if( ( 0 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( -1.0 / tauMaxEnergy );
        if( !isInsideRange )
        {
          returnValue += ( 1.0 / tauMinEnergy );
        }
      }
      else if( ( 3 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( 0.25 / ( tauMaxEnergy * tauMaxEnergy
                                 * tauMaxEnergy * tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue -= ( 0.25 / ( tauMinEnergy * tauMinEnergy
                                    * tauMinEnergy * tauMinEnergy ) );
        }
      }
    }  // end of if( -1 == tau_energy_power ) && ( 0 == tau_log_power )
    else if( ( 0 == tauLeptonEnergyPower )
             &&
             ( 0 == tauLeptonLogPower ) )
    {
      if( ( 0 == visibleProductEnergyPower )
          &&
          ( 0 == visibleProductLogPower ) )
      {
        returnValue = log( tauMaxEnergy );
        if( isInsideRange )
        {
          returnValue -= ( 1.0 / 3.0 );
        }
        else
        {
          returnValue -= log( tauMinEnergy );
        }
      }
      else if( ( 0 == visibleProductEnergyPower )
               &&
               ( 1 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = -1.0;
        }
      }
      else if( ( 3 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue
        = ( 1.0 / ( 3.0 * tauMaxEnergy * tauMaxEnergy * tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue
          -= ( 1.0 / ( 3.0 * tauMinEnergy * tauMinEnergy * tauMinEnergy ) );
        }
      }
    }  // end of if( 0 == tau_energy_power ) && ( 0 == tau_log_power )
    else if( ( 0 == tauLeptonEnergyPower )
             &&
             ( 1 == tauLeptonLogPower ) )
    {
      if( ( 0 == visibleProductEnergyPower )
          &&
          ( 0 == visibleProductLogPower ) )
      {
        double logOfTauMaxEnergy( log( tauMaxEnergy ) );
        returnValue
        = ( 0.5 * logOfTauMaxEnergy * logOfTauMaxEnergy );
        if( isInsideRange )
        {
          returnValue -= ( 1.0 / 9.0 );
        }
        else
        {
          double logOfTauMinEnergy( log( tauMinEnergy ) );
          returnValue
          -= ( 0.5 * logOfTauMinEnergy * logOfTauMinEnergy );
        }
      }
      else if( ( 0 == visibleProductEnergyPower )
              &&
              ( 1 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = ( -1.0 / 3.0 );
        }
      }
      else if( ( 0 == visibleProductEnergyPower )
               &&
               ( 2 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = -0.5;
        }
      }
      else if( ( 3 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue
        = ( ( 1.0 + 3.0 * log( tauMaxEnergy ) )
            / ( 9.0 * tauMaxEnergy * tauMaxEnergy * tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue
          -= ( ( 1.0 + 3.0 * log( tauMinEnergy ) )
               / ( 9.0 * tauMinEnergy * tauMinEnergy * tauMinEnergy ) );
        }
      }
    }  // end of if( 0 == tau_energy_power ) && ( 1 == tau_log_power )
    else if( ( 1 == tauLeptonEnergyPower )
             &&
             ( 0 == tauLeptonLogPower ) )
    {
      if( ( 0 == visibleProductEnergyPower )
          &&
          ( 0 == visibleProductLogPower ) )
      {
        returnValue = tauMaxEnergy;
        if( !isInsideRange )
        {
          returnValue -= tauMinEnergy;
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = -1.5;
        }
      }
      else if( ( 3 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( 0.5 / ( tauMaxEnergy * tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue -= ( 0.5 / ( tauMinEnergy * tauMinEnergy ) );
        }
      }
    }  // end of if( 1 == tau_energy_power ) && ( 0 == tau_log_power )
    else if( ( 1 == tauLeptonEnergyPower )
             &&
             ( 1 == tauLeptonLogPower ) )
    {
      if( ( 0 == visibleProductEnergyPower )
          &&
          ( 0 == visibleProductLogPower ) )
      {
        returnValue
        = ( tauMaxEnergy * ( log( tauMaxEnergy ) - 1.0 ) );
        if( !isInsideRange )
        {
          returnValue
          += ( tauMinEnergy * ( 1.0 - log( tauMinEnergy ) ) );
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = 0.75;
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 1 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = -1.5;
        }
      }
      else if( ( 3 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( ( 0.25 + 0.5 * log( tauMaxEnergy ) )
                        / ( tauMaxEnergy * tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue
          -= ( ( 0.25 + 0.5 * log( tauMinEnergy ) )
               / ( tauMinEnergy * tauMinEnergy ) );
        }
      }
    }  // end of if( 1 == tau_energy_power ) && ( 1 == tau_log_power )
    else if( ( 2 == tauLeptonEnergyPower )
             &&
             ( 0 == tauLeptonLogPower ) )
    {
      if( ( 0 == visibleProductEnergyPower )
          &&
          ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( 0.5 * tauMaxEnergy * tauMaxEnergy );
        if( !isInsideRange )
        {
          returnValue -= ( 0.5 * tauMinEnergy * tauMinEnergy );
        }
      }
      else if( ( 2 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = -1.5;
        }
      }
      else if( ( 3 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( 1.0 / tauMaxEnergy );
        if( !isInsideRange )
        {
          returnValue -= ( 1.0 / tauMinEnergy );
        }
      }
    }  // end of if( 2 == tau_energy_power ) && ( 0 == tau_log_power )
    return returnValue;
  }



  softMuonFromTau::softMuonFromTau() :
    tauDecayCoefficient()
  {
    // just an initialization.
  }

  softMuonFromTau::~softMuonFromTau()
  {
    // does nothing.
  }


  double
  softMuonFromTau::operator()( int const visibleProductEnergyPower
                            /* the power of the muon's energy for the term. */,
                               int const visibleProductLogPower
           /* the power of the logarithm of the muon's energy for the term. */,
                               int const tauLeptonEnergyPower
                       /* the power of the tau lepton's energyfor the term. */,
                               int const tauLeptonLogPower
     /* the power of the logarithm of the tau lepton's energy for the term. */,
                               double tauMinEnergy
         /* the minimum of the range of the tau lepton energies considered. */,
                               double const tauMaxEnergy
         /* the maximum of the range of the tau lepton energies considered. */,
                               bool const isInsideRange )
  const
  /* this returns the coefficientValue for the muon distribution which would come
   * from a right-handed tau lepton distribution of the given power of the tau
   * lepton's energy & power of logarithm thereof, for requested integer powers
   * of the muon's energy or logarithm thereof, which also depends on whether
   * the muon's energy is inside the range of the tau lepton's energy for this
   * segment.
   */
  {
    if( leptonDistributionExpansionTerm::minimumInputEnergy > tauMinEnergy )
    {
      tauMinEnergy = leptonDistributionExpansionTerm::minimumInputEnergy;
    }
    double returnValue( 0.0 );
    if( ( -1 == tauLeptonEnergyPower )
        &&
        ( 0 == tauLeptonLogPower ) )
    {
      if( ( -1 == visibleProductEnergyPower )
          &&
          ( 0 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = 0.5;
        }
      }
      else if( ( 0 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( -1.0 / tauMaxEnergy );
        if( !isInsideRange )
        {
          returnValue += ( 1.0 / tauMinEnergy );
        }
      }
      else if( ( 2 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( 1.0 / ( tauMaxEnergy * tauMaxEnergy * tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue
          -= ( 1.0 / ( tauMinEnergy * tauMinEnergy * tauMinEnergy ) );
        }
      }
      else if( ( 3 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( -0.5 / ( tauMaxEnergy * tauMaxEnergy
                                 * tauMaxEnergy * tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue += ( 0.5 / ( tauMinEnergy * tauMinEnergy
                                   * tauMinEnergy * tauMinEnergy ) );
        }
      }
    }  // end of if( -1 == tau_energy_power ) && ( 0 == tau_log_power )
    else if( ( 0 == tauLeptonEnergyPower )
             &&
             ( 0 == tauLeptonLogPower ) )
    {
      if( ( 0 == visibleProductEnergyPower )
          &&
          ( 0 == visibleProductLogPower ) )
      {
        returnValue = log( tauMaxEnergy );
        if( isInsideRange )
        {
          returnValue -= ( 5.0 / 6.0 );
        }
        else
        {
          returnValue -= log( tauMinEnergy );
        }
      }
      else if( ( 0 == visibleProductEnergyPower )
               &&
               ( 1 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = -1.0;
        }
      }
      else if( ( 2 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( 1.5 / ( tauMaxEnergy * tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue -= ( 1.5 / ( tauMinEnergy * tauMinEnergy ) );
        }
      }
      else if( ( 3 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue
        = ( -2.0 / ( 3.0 * tauMaxEnergy * tauMaxEnergy * tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue
          += ( 2.0 / ( 3.0 * tauMinEnergy * tauMinEnergy * tauMinEnergy ) );
        }
      }
    }  // end of if( 0 == tau_energy_power ) && ( 0 == tau_log_power )
    else if( ( 0 == tauLeptonEnergyPower )
             &&
             ( 1 == tauLeptonLogPower ) )
    {
      if( ( 0 == visibleProductEnergyPower )
          &&
          ( 0 == visibleProductLogPower ) )
      {
        double logOfTauMaxEnergy( log( tauMaxEnergy ) );
        returnValue = ( 0.5 * logOfTauMaxEnergy * logOfTauMaxEnergy );
        if( isInsideRange )
        {
          returnValue -= ( 19.0 / 36.0 );
        }
        else
        {
          double logOfTauMinEnergy = log( tauMinEnergy );
          returnValue -= ( 0.5 * logOfTauMinEnergy * logOfTauMinEnergy );
        }
      }
      else if( ( 0 == visibleProductEnergyPower )
               &&
               ( 1 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = ( -5.0 / 6.0 );
        }
      }
      else if( ( 0 == visibleProductEnergyPower )
               &&
               ( 2 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = -0.5;
        }
      }
      else if( ( 2 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( ( 0.75 + 1.5 * log( tauMaxEnergy ) )
                        / ( tauMaxEnergy * tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue -= ( ( 0.75 + 1.5 * log( tauMinEnergy ) )
                           / ( tauMinEnergy * tauMinEnergy ) );
        }
      }
      else if( ( 3 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( ( -2.0 * ( 1.0 + 3.0 * log( tauMaxEnergy ) ) )
                      / ( 9.0 * tauMaxEnergy * tauMaxEnergy * tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue += ( ( 2.0 * ( 1.0 + 3.0 * log( tauMinEnergy ) ) )
                      / ( 9.0 * tauMinEnergy * tauMinEnergy * tauMinEnergy ) );
        }
      }
    }  // end of if( 0 == tau_energy_power ) && ( 1 == tau_log_power )
    else if( ( 1 == tauLeptonEnergyPower )
             &&
             ( 0 == tauLeptonLogPower ) )
    {
      if( ( 0 == visibleProductEnergyPower )
          &&
          ( 0 == visibleProductLogPower ) )
      {
        returnValue = tauMaxEnergy;
        if( !isInsideRange )
        {
          returnValue -= tauMinEnergy;
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = -3.0;
        }
      }
      else if( ( 2 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( 3.0 / tauMaxEnergy );
        if( !isInsideRange )
        {
          returnValue -= ( 3.0 / tauMinEnergy );
        }
      }
      else if( ( 3 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( -1.0 / ( tauMaxEnergy * tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue += ( 1.0 / ( tauMinEnergy * tauMinEnergy ) );
        }
      }
    }  // end of if( 1 == tau_energy_power ) && ( 0 == tau_log_power )
    else if( ( 1 == tauLeptonEnergyPower )
             &&
             ( 1 == tauLeptonLogPower ) )
    {
      if( ( 0 == visibleProductEnergyPower )
          &&
          ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( tauMaxEnergy * ( log( tauMaxEnergy ) - 1.0 ) );
        if( !isInsideRange )
        {
          returnValue
          += ( tauMinEnergy * ( 1.0 - log( tauMinEnergy ) ) );
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = -1.5;
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 1 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = -3.0;
        }
      }
      else if( ( 2 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue
        = ( ( 3.0 * ( 1.0 + log( tauMaxEnergy ) ) ) / tauMaxEnergy );
        if( !isInsideRange )
        {
          returnValue
          -= ( ( 3.0 * ( 1.0 + log( tauMinEnergy ) ) ) / tauMinEnergy );
        }
      }
      else if( ( 3 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( ( -0.5 - log( tauMaxEnergy ) )
                        / ( tauMaxEnergy * tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue += ( ( 0.5 + log( tauMinEnergy ) )
                           / ( tauMinEnergy * tauMinEnergy ) );
        }
      }
    }  // end of if( 1 == tau_energy_power ) && ( 1 == tau_log_power )
    else if( ( 2 == tauLeptonEnergyPower )
             &&
             ( 0 == tauLeptonLogPower ) )
    {
      if( ( 0 == visibleProductEnergyPower )
          &&
          ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( 0.5 * tauMaxEnergy * tauMaxEnergy );
        if( !isInsideRange )
        {
          returnValue -= ( 0.5 * tauMinEnergy * tauMinEnergy );
        }
      }
      else if( ( 2 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( -3.0 * log( tauMaxEnergy ) );
        if( isInsideRange )
        {
          returnValue += 1.5;
        }
        else
        {
          returnValue += ( 3.0 * log( tauMinEnergy ) );
        }
      }
      else if( ( 2 == visibleProductEnergyPower )
               &&
               ( 1 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = 3.0;
        }
      }
      else if( ( 3 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( -2.0 / ( tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue += ( 2.0 / ( tauMinEnergy ) );
        }
      }
    }  // end of if( 2 == tau_energy_power ) && ( 0 == tau_log_power )
    return returnValue;
  }



  hardPionFromTau::hardPionFromTau() :
    tauDecayCoefficient()
  {
    // just an initialization.
  }

  hardPionFromTau::~hardPionFromTau()
  {
    // does nothing.
  }


  double
  hardPionFromTau::operator()( int const visibleProductEnergyPower
                            /* the power of the pion's energy for the term. */,
                               int const visibleProductLogPower
           /* the power of the logarithm of the pion's energy for the term. */,
                               int const tauLeptonEnergyPower
                      /* the power of the tau lepton's energy for the term. */,
                               int const tauLeptonLogPower
     /* the power of the logarithm of the tau lepton's energy for the term. */,
                               double tauMinEnergy
         /* the minimum of the range of the tau lepton energies considered. */,
                               double const tauMaxEnergy
         /* the maximum of the range of the tau lepton energies considered. */,
                               bool const isInsideRange )
  const
  /* this returns the coefficientValue for the pion distribution which would come
   * from a right-handed tau lepton distribution of the given power of the
   * tau lepton's energy & power of logarithm thereof, for requested integer
   * powers of the pion's energy or logarithm thereof, which also depends on
   * whether the pion's energy is inside the range of the tau lepton's energy
   * for this segment.
   */
  {
    if( leptonDistributionExpansionTerm::minimumInputEnergy > tauMinEnergy )
    {
      tauMinEnergy = leptonDistributionExpansionTerm::minimumInputEnergy;
    }
    double returnValue( 0.0 );
    if( ( -1 == tauLeptonEnergyPower )
        &&
        ( 0 == tauLeptonLogPower ) )
    {
      if( ( -1 == visibleProductEnergyPower )
          &&
          ( 0 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = 0.25;
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( ( -0.25 - 0.5 * log( tauMaxEnergy ) )
                        / ( tauMaxEnergy * tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue += ( ( 0.25 + 0.5 * log( tauMinEnergy ) )
                           / ( tauMinEnergy * tauMinEnergy ) );
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 1 == visibleProductLogPower ) )
      {
        returnValue = ( 0.5 / ( tauMaxEnergy * tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue -= ( 0.5 / ( tauMinEnergy * tauMinEnergy ) );
        }
      }
    }  // end of if( -1 == tau_energy_power ) && ( 0 == tau_log_power )
    else if( ( 0 == tauLeptonEnergyPower )
             &&
             ( 0 == tauLeptonLogPower ) )
    {
      if( ( 0 == visibleProductEnergyPower )
          &&
          ( 0 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = 1.0;
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( ( -1.0 - log( tauMaxEnergy ) ) / tauMaxEnergy );
        if( !isInsideRange )
        {
          returnValue += ( ( 1.0 + log( tauMinEnergy ) )
                           / tauMinEnergy );
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 1 == visibleProductLogPower ) )
      {
        returnValue = ( 1.0 / tauMaxEnergy );
        if( !isInsideRange )
        {
          returnValue -= ( 1.0 / tauMinEnergy );
        }
      }
    }  // end of if( 0 == tau_energy_power ) && ( 0 == tau_log_power )
    else if( ( 0 == tauLeptonEnergyPower )
             &&
             ( 1 == tauLeptonLogPower ) )
    {
      if( ( 0 == visibleProductEnergyPower )
          &&
          ( 0 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = 2.0;
        }
      }
      else if( ( 0 == visibleProductEnergyPower )
               &&
               ( 1 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = 1.0;
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        double logOfTauMaxEnergy( log( tauMaxEnergy ) );
        returnValue
        = ( ( -2.0 - logOfTauMaxEnergy * ( 2.0 + logOfTauMaxEnergy ) )
            / tauMaxEnergy );
        if( !isInsideRange )
        {
          double logOfTauMinEnergy( log( tauMinEnergy ) );
          returnValue
          += ( ( 2.0 + logOfTauMinEnergy * ( 2.0 + logOfTauMinEnergy ) )
               / tauMinEnergy );
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 1 == visibleProductLogPower ) )
      {
        returnValue = ( ( 1.0 + log( tauMaxEnergy ) ) / tauMaxEnergy );
        if( !isInsideRange )
        {
          returnValue -= ( ( 1.0 + log( tauMinEnergy ) ) / tauMinEnergy );
        }
      }
    }  // end of if( 0 == tau_energy_power ) && ( 1 == tau_log_power )
    else if( ( 1 == tauLeptonEnergyPower )
             &&
             ( 0 == tauLeptonLogPower ) )
    {
      if( ( 1 == visibleProductEnergyPower )
          &&
          ( 0 == visibleProductLogPower ) )
      {
        double logOfTauMaxEnergy( log( tauMaxEnergy ) );
        returnValue = ( 0.5 * logOfTauMaxEnergy * logOfTauMaxEnergy );
        if( !isInsideRange )
        {
          double logOfTauMinEnergy( log( tauMinEnergy ) );
          returnValue -= ( 0.5 * logOfTauMinEnergy * logOfTauMinEnergy );
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 1 == visibleProductLogPower ) )
      {
        returnValue = -log( tauMaxEnergy );
        if( !isInsideRange )
        {
          returnValue += log( tauMinEnergy );
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 2 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = 0.5;
        }
      }
    }  // end of if( 1 == tau_energy_power ) && ( 0 == tau_log_power )
    else if( ( 1 == tauLeptonEnergyPower )
             &&
             ( 1 == tauLeptonLogPower ) )
    {
      if( ( 1 == visibleProductEnergyPower )
          &&
          ( 0 == visibleProductLogPower ) )
      {
        double logOfTauMaxEnergy( log( tauMaxEnergy ) );
        returnValue
        = ( ( logOfTauMaxEnergy * logOfTauMaxEnergy * logOfTauMaxEnergy )
            / 3.0 );
        if( !isInsideRange )
        {
          double logOfTauMinEnergy( log( tauMinEnergy ) );
          returnValue -= ( ( logOfTauMinEnergy
                             * logOfTauMinEnergy * logOfTauMinEnergy ) / 3.0 );
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 1 == visibleProductLogPower ) )
      {
        double logOfTauMaxEnergy( log( tauMaxEnergy ) );
        returnValue = ( -0.5 * logOfTauMaxEnergy * logOfTauMaxEnergy );
        if( !isInsideRange )
        {
          double logOfTauMinEnergy( log( tauMinEnergy ) );
          returnValue += ( 0.5 * logOfTauMinEnergy * logOfTauMinEnergy );
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 3 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = ( 1.0 / 6.0 );
        }
      }
    }  // end of if( 1 == tau_energy_power ) && ( 1 == tau_log_power )
    else if( ( 2 == tauLeptonEnergyPower )
             &&
             ( 0 == tauLeptonLogPower ) )
    {
      if( ( 1 == visibleProductEnergyPower )
          &&
          ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( tauMaxEnergy * ( log( tauMaxEnergy ) - 1.0 ) );
        if( !isInsideRange )
        {
          returnValue += ( tauMinEnergy * ( 1.0 - log( tauMinEnergy ) ) );
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 1 == visibleProductLogPower ) )
      {
        returnValue = -tauMaxEnergy;
        if( !isInsideRange )
        {
          returnValue += tauMinEnergy;
        }
      }
      else if( ( 2 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = 1.0;
        }
      }
    }  // end of if( 2 == tau_energy_power ) && ( 0 == tau_log_power )
    return returnValue;
  }



  softPionFromTau::softPionFromTau() :
    tauDecayCoefficient()
  {
    // just an initialization.
  }

  softPionFromTau::~softPionFromTau()
  {
    // does nothing.
  }


  double
  softPionFromTau::operator()( int const visibleProductEnergyPower
                            /* the power of the pion's energy for the term. */,
                               int const visibleProductLogPower
           /* the power of the logarithm of the pion's energy for the term. */,
                               int const tauLeptonEnergyPower
                       /* the power of the tau lepton's energyfor the term. */,
                               int const tauLeptonLogPower
     /* the power of the logarithm of the tau lepton's energy for the term. */,
                               double tauMinEnergy
         /* the minimum of the range of the tau lepton energies considered. */,
                               double const tauMaxEnergy
         /* the maximum of the range of the tau lepton energies considered. */,
                               bool const isInsideRange )
  const
  /* this returns the coefficientValue for the pion distribution which would come
   * from a left-handed tau lepton distribution of the given power of the tau
   * lepton's energy & power of logarithm thereof, for requested integer
   * powers of the pion's energy or logarithm thereof, which also depends on
   * whether the pion's energy is inside the range of the tau lepton's energy
   * for this segment.
   */
  {
    if( leptonDistributionExpansionTerm::minimumInputEnergy > tauMinEnergy )
    {
      tauMinEnergy = leptonDistributionExpansionTerm::minimumInputEnergy;
    }
    double returnValue = 0.0;
    if( ( -1 == tauLeptonEnergyPower )
        &&
        ( 0 == tauLeptonLogPower ) )
    {
      if( ( -1 == visibleProductEnergyPower )
          &&
          ( 0 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = 0.25;
        }
      }
      else if( ( 0 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( -1.0 / tauMaxEnergy );
        if( !isInsideRange )
        {
          returnValue += ( 1.0 / tauMinEnergy );
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( ( 0.5 * log( tauMaxEnergy ) + 0.75 )
                        / ( tauMaxEnergy * tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue -= ( ( 0.5 * log( tauMinEnergy ) + 0.75 )
                           / ( tauMinEnergy * tauMinEnergy ) );
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 1 == visibleProductLogPower ) )
      {
        returnValue = ( -0.5 / ( tauMaxEnergy * tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue += ( 0.5 / ( tauMinEnergy * tauMinEnergy ) );
        }
      }
    }  // end of if( -1 == tau_energy_power ) && ( 0 == tau_log_power )
    else if( ( 0 == tauLeptonEnergyPower )
             &&
             ( 0 == tauLeptonLogPower ) )
    {
      if( ( 0 == visibleProductEnergyPower )
          &&
          ( 0 == visibleProductLogPower ) )
      {
        returnValue = log( tauMaxEnergy );
        if( isInsideRange )
        {
          returnValue -= 2.0;
        }
        else
        {
          returnValue -= log( tauMinEnergy );
        }
      }
      else if( ( 0 == visibleProductEnergyPower )
               &&
               ( 1 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = -1.0;
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( ( 2.0 + log( tauMaxEnergy ) ) / tauMaxEnergy );
        if( !isInsideRange )
        {
          returnValue -= ( ( 2.0 + log( tauMinEnergy ) ) / tauMinEnergy );
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 1 == visibleProductLogPower ) )
      {
        returnValue = ( -1.0 / tauMaxEnergy );
        if( !isInsideRange )
        {
          returnValue += ( 1.0 / tauMinEnergy );
        }
      }
    }  // end of if( 0 == tau_energy_power ) && ( 0 == tau_log_power )
    else if( ( 0 == tauLeptonEnergyPower )
             &&
             ( 1 == tauLeptonLogPower ) )
    {
      if( ( 0 == visibleProductEnergyPower )
          &&
          ( 0 == visibleProductLogPower ) )
      {
        double logOfTauMaxEnergy( log( tauMaxEnergy ) );
        returnValue = ( 0.5 * logOfTauMaxEnergy * logOfTauMaxEnergy );
        if( isInsideRange )
        {
          returnValue -= 3.0;
        }
        else
        {
          double logOfTauMinEnergy( log( tauMinEnergy ) );
          returnValue -= ( 0.5 * logOfTauMinEnergy * logOfTauMinEnergy );
        }
      }
      else if( ( 0 == visibleProductEnergyPower )
               &&
               ( 1 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = -2.0;
        }
      }
      else if( ( 0 == visibleProductEnergyPower )
               &&
               ( 2 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = -0.5;
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        double logOfTauMaxEnergy( log( tauMaxEnergy ) );
        returnValue
        = ( ( 3.0 + logOfTauMaxEnergy * ( 3.0 + logOfTauMaxEnergy ) )
            / tauMaxEnergy );
        if( !isInsideRange )
        {
          double logOfTauMinEnergy( log( tauMinEnergy ) );
          returnValue
          -= ( ( 3.0 + logOfTauMinEnergy * ( 3.0 + logOfTauMinEnergy ) )
               / tauMinEnergy );
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 1 == visibleProductLogPower ) )
      {
        returnValue = ( ( -1.0 - log( tauMaxEnergy ) ) / tauMaxEnergy );
        if( !isInsideRange )
        {
          returnValue += ( ( 1.0 + log( tauMinEnergy ) ) / tauMinEnergy );
        }
      }
    }  // end of if( 0 == tau_energy_power ) && ( 1 == tau_log_power )
    else if( ( 1 == tauLeptonEnergyPower )
             &&
             ( 0 == tauLeptonLogPower ) )
    {
      if( ( 0 == visibleProductEnergyPower )
          &&
          ( 0 == visibleProductLogPower ) )
      {
        returnValue = tauMaxEnergy;
        if( !isInsideRange )
        {
          returnValue -= tauMinEnergy;
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        double logOfTauMaxEnergy( log( tauMaxEnergy ) );
        returnValue
        = ( -logOfTauMaxEnergy * ( 0.5 * logOfTauMaxEnergy + 1.0 ) );
        if( isInsideRange )
        {
          returnValue -= 1.0;
        }
        else
        {
          double logOfTauMinEnergy( log( tauMinEnergy ) );
          returnValue
          += ( logOfTauMinEnergy * ( 0.5 * logOfTauMinEnergy + 1.0 ) );
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 1 == visibleProductLogPower ) )
      {
        returnValue = log( tauMaxEnergy );
        if( isInsideRange )
        {
          returnValue += 1.0;
        }
        else
        {
          returnValue -= log( tauMinEnergy );
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 2 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = -0.5;
        }
      }
    }  // end of if( 1 == tau_energy_power ) && ( 0 == tau_log_power )
    else if( ( 1 == tauLeptonEnergyPower )
             &&
             ( 1 == tauLeptonLogPower ) )
    {
      if( ( 0 == visibleProductEnergyPower )
          &&
          ( 0 == visibleProductLogPower ) )
      {
        returnValue
        = ( tauMaxEnergy * ( log( tauMaxEnergy ) - 1.0 ) );
        if( !isInsideRange )
        {
          returnValue
          += ( tauMinEnergy * ( 1.0 - log( tauMinEnergy ) ) );
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        double logOfTauMaxEnergy( log( tauMaxEnergy ) );
        returnValue = ( ( -logOfTauMaxEnergy * logOfTauMaxEnergy
                          * ( 1.5 + logOfTauMaxEnergy ) ) / 3.0 );
        if( isInsideRange )
        {
          returnValue += 1.0;
        }
        else
        {
          double logOfTauMinEnergy( log( tauMinEnergy ) );
          returnValue += ( ( logOfTauMinEnergy * logOfTauMinEnergy
                             * ( 1.5 + logOfTauMinEnergy ) ) / 3.0 );
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 1 == visibleProductLogPower ) )
      {
        double logOfTauMaxEnergy( log( tauMaxEnergy ) );
        returnValue = ( 0.5 * logOfTauMaxEnergy * logOfTauMaxEnergy );
        if( isInsideRange )
        {
          returnValue -= 1.0;
        }
        else
        {
          double logOfTauMinEnergy( log( tauMinEnergy ) );
          returnValue -= ( 0.5 * logOfTauMinEnergy * logOfTauMinEnergy );
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 2 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = 0.5;
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 3 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = ( -1 / 6.0 );
        }
      }
    }  // end of if( 1 == tau_energy_power ) && ( 1 == tau_log_power )
    else if( ( 2 == tauLeptonEnergyPower )
             &&
             ( 0 == tauLeptonLogPower ) )
    {
      if( ( 0 == visibleProductEnergyPower )
          &&
          ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( 0.5 * tauMaxEnergy * tauMaxEnergy );
        if( !isInsideRange )
        {
          returnValue -= ( 0.5 * tauMinEnergy * tauMinEnergy );
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( -tauMaxEnergy * log( tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue += ( tauMinEnergy * log( tauMinEnergy ) );
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
                &&
                ( 1 == visibleProductLogPower ) )
      {
        returnValue = tauMaxEnergy;
        if( !isInsideRange )
        {
          returnValue -= tauMinEnergy;
        }
      }
      else if( ( 2 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = -0.5;
        }
      }
    }  // end of if( 2 == tau_energy_power ) && ( 0 == tau_log_power )
    return returnValue;
  }



  tauSegmentTriple::tauSegmentTriple(
                          std::vector< segmentTermSet* >* segmentSetToPopulate,
                                      segmentTermSet* referenceSegment,
                                  tauDecayCoefficient const* const tauDecay ) :
    referenceSegment( referenceSegment ),
    tauDecay( tauDecay )
  {
    double coefficientChecker;
    lowSegment = new segmentTermSet();
    lowSegment->setSegmentStart(
                         leptonDistributionExpansionTerm::minimumInputEnergy );
    highSegment = new segmentTermSet();
    for( std::vector< leptonDistributionExpansionTerm* >::const_iterator
         termIterator( referenceSegment->getTerms()->begin() );
         referenceSegment->getTerms()->end() > termIterator;
         ++termIterator )
    {
      for( int energyPowerCounter( -1 );
           3 >= energyPowerCounter;
           ++energyPowerCounter )
      {
        for( int logPowerCounter( 0 );
             3 >= logPowerCounter;
             ++logPowerCounter )
        {
          /* each term is checked to see if it's non-zero. "3.1" & "4.3"
           * were chosen since there should be no accidental
           * cancellations with these acceptanceValues, since there are just
           * logarithms, powers of the end energies, & some rational
           * numbers, but none are multiples of 3.1 or 4.3 so there can
           * be no cancellation.
           */
          coefficientChecker = (*tauDecay)( energyPowerCounter,
                                            logPowerCounter,
                                           (*termIterator)->getPowerOfEnergy(),
                                        (*termIterator)->getPowerOfLogarithm(),
                                            3.1,
                                            4.3,
                                            false );
          if( 0.0 != coefficientChecker )
          {
            // debugging:
            /**std::cout << std::endl << "debugging:"
            << std::endl
            << "tauSegmentTriple::tauSegmentTriple(...) adding E^"
            << energyPowerCounter << " L^" << logPowerCounter << " to low.";
            std::cout << std::endl;**/

            lowSegment->addTerm( energyPowerCounter,
                                 logPowerCounter,
                                 *termIterator );
          }
          coefficientChecker = (*tauDecay)( energyPowerCounter,
                                            logPowerCounter,
                                           (*termIterator)->getPowerOfEnergy(),
                                        (*termIterator)->getPowerOfLogarithm(),
                                            3.1,
                                            4.3,
                                            true );
          if( 0.0 != coefficientChecker )
          {
            // debugging:
            /**std::cout << std::endl << "debugging:"
            << std::endl
            << "tauSegmentTriple::tauSegmentTriple(...) adding E^"
            << energyPowerCounter << " L^" << logPowerCounter << " to high.";
            std::cout << std::endl;**/
            highSegment->addTerm( energyPowerCounter,
                                  logPowerCounter,
                                  *termIterator );
          }
        }  // end of loop over powers of log( energy ).
      }  // end of loop over powers of energy.
    }  // end of loop over leptonDistributionExpansionTerm pointers.

    // now that all the lowSegment & highSegment terms have been added to the
    // segments, we add them to the vector of segmentTermSets.
    segmentSetToPopulate->push_back( lowSegment );
    segmentSetToPopulate->push_back( highSegment );
  }

  tauSegmentTriple::~tauSegmentTriple()
  {

    // does nothing.

  }


  void
  tauSegmentTriple::updateSegments()
  {
    lowSegment->setSegmentEnd( referenceSegment->getSegmentStart() );
    highSegment->setSegmentRange( referenceSegment->getSegmentStart(),
                                  referenceSegment->getSegmentEnd() );

    // debugging:
    /**std::cout << std::endl << "debugging:"
    << std::endl
    << "tauSegmentTriple::updateSegments() called. lowSegment = "
    << lowSegment->getSegmentStart() << " to " << lowSegment->getSegmentEnd()
    << ", highSegment = " << highSegment->getSegmentStart() << " to "
    << highSegment->getSegmentEnd();
    std::cout << std::endl;**/

    for( std::vector< leptonDistributionExpansionTerm* >::const_iterator
         term_iterator( lowSegment->getTerms()->begin() );
         lowSegment->getTerms()->end() > term_iterator;
         ++term_iterator )
    {
      (*term_iterator)->setCoefficient(
                         (*term_iterator)->getReferenceTerm()->getCoefficient()
                           * (*tauDecay)( (*term_iterator)->getPowerOfEnergy(),
                                       (*term_iterator)->getPowerOfLogarithm(),
                      (*term_iterator)->getReferenceTerm()->getPowerOfEnergy(),
                   (*term_iterator)->getReferenceTerm()->getPowerOfLogarithm(),
                                        referenceSegment->getSegmentStart(),
                                        referenceSegment->getSegmentEnd(),
                                        false ) );
      // debugging:
      /**std::cout << std::endl << "debugging:"
      << std::endl
      << "(lowSegment) reference term: "
      << (*term_iterator)->getReferenceTerm()->getCoefficient()
      << " * E^"
      << (*term_iterator)->getReferenceTerm()->getPowerOfEnergy()
      << " * L^"
      << (*term_iterator)->getReferenceTerm()->getPowerOfLogarithm()
      << std::endl
      << "=> " << (*term_iterator)->getCoefficient()
      << " * E^" << (*term_iterator)->getPowerOfEnergy()
      << " * L^" << (*term_iterator)->getPowerOfLogarithm();
      std::cout << std::endl;**/
    }
    for( std::vector< leptonDistributionExpansionTerm* >::const_iterator
         term_iterator( highSegment->getTerms()->begin() );
         highSegment->getTerms()->end() > term_iterator;
         ++term_iterator )
    {
      (*term_iterator)->setCoefficient(
                         (*term_iterator)->getReferenceTerm()->getCoefficient()
                           * (*tauDecay)( (*term_iterator)->getPowerOfEnergy(),
                                       (*term_iterator)->getPowerOfLogarithm(),
                      (*term_iterator)->getReferenceTerm()->getPowerOfEnergy(),
                   (*term_iterator)->getReferenceTerm()->getPowerOfLogarithm(),
                                       referenceSegment->getSegmentStart(),
                                       referenceSegment->getSegmentEnd(),
                                       true ) );
    }
  }



  leptonEnergyDistribution::leptonEnergyDistribution(
                                          readierForNewPoint* const readier,
                                       CppSLHA::CppSLHA0 const* const spectrum,
                     CppSLHA::particle_property_set const* const firstParticle,
                          effectiveSquarkMassHolder* const effectiveSquarkMass,
                    CppSLHA::particle_property_set const* const secondParticle,
                     CppSLHA::particle_property_set const* const thirdParticle,
                 CppSLHA::particle_property_set const* const fourthParticle ) :
    getsReadiedForNewPoint( readier ),
    spectrum( spectrum ),
    firstParticle( firstParticle ),
    effectiveSquarkMass( effectiveSquarkMass ),
    secondParticle( secondParticle ),
    thirdParticle( thirdParticle ),
    fourthParticle( fourthParticle )
  {
    // just an initialization list.
  }

  leptonEnergyDistribution::~leptonEnergyDistribution()
  {
    // does nothing, because the segmentTermSets are not dynamically allocated.
  }


  void
  leptonEnergyDistribution::reset()
  {
    if( NULL != effectiveSquarkMass )
    {
      firstMass = effectiveSquarkMass->getEffectiveSquarkMass();
    }
    else
    {
      firstMass = CppSLHA::CppSLHA_global::really_wrong_value;
    }
    if( NULL != secondParticle )
    {
      secondMass = secondParticle->get_absolute_mass();
    }
    else
    {
      secondMass = CppSLHA::CppSLHA_global::really_wrong_value;
    }
    if( NULL != thirdParticle )
    {
      thirdMass = thirdParticle->get_absolute_mass();
    }
    else
    {
      thirdMass = CppSLHA::CppSLHA_global::really_wrong_value;
    }
    if( NULL != fourthParticle )
    {
      fourthMass = fourthParticle->get_absolute_mass();
    }
    else
    {
      fourthMass = CppSLHA::CppSLHA_global::really_wrong_value;
    }
    calculateCoefficients();
    normalizeCoefficients();
  }



  visibleTauDecayProduct::visibleTauDecayProduct(
                                    readierForNewPoint* const given_readier,
                      leptonEnergyDistribution* const given_tau_distribution,
                               tauDecayCoefficient const* const tau_decay ) :
    leptonEnergyDistribution( given_readier,
                              NULL,
                              NULL,
                              NULL,
                              NULL,
                              NULL,
                              NULL ),
    tauDistribution( given_tau_distribution )
  {
    // debugging:
    /**std::cout
    << std::endl
    << "debugging: "
    << "visibleTauDecayProduct( "
    << tauDistribution << ", " << tauDecay << " )";
    std::cout
    << std::endl;**/
    minimumEnergy = leptonDistributionExpansionTerm::minimumInputEnergy;

    for( std::vector< segmentTermSet* >::iterator
         segmentIterator( given_tau_distribution->getSegments()->begin() );
         given_tau_distribution->getSegments()->end() > segmentIterator;
         ++segmentIterator )
      {

        tauTriples.push_back( new tauSegmentTriple( &segments,
                                                    *segmentIterator,
                                                    tau_decay ) );

      }  // end of loop through segments.

  }

  visibleTauDecayProduct::~visibleTauDecayProduct()
  {

    for( std::vector< tauSegmentTriple* >::iterator
         deletion_iterator = tauTriples.begin();
         tauTriples.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

  }


  void
  visibleTauDecayProduct::calculateCoefficients()
  {
    // debugging:
    /**std::cout
    << std::endl
    << "debugging: "
    << "visibleTauDecayProduct::calculateCoefficients()";
    std::cout
    << std::endl;**/

    maximumEnergy = tauDistribution->getMaximumEnergy();
    // now all the underlying stuff should have been updated for the point.
    for( std::vector< tauSegmentTriple* >::iterator
         tripleIterator( tauTriples.begin() );
         tauTriples.end() > tripleIterator;
         ++tripleIterator )
      {
        (*tripleIterator)->updateSegments();
      }  // end of loop through segments.
  }

}  // end of LHC_FASER namespace

