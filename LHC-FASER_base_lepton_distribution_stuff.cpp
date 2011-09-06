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
 *      LHC-FASER also requires grids of lookup values. These should also be
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
        return_value *= logarithmOfEnergy;
      }

      return returnValue;
    }
    else
    {
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
        double energyPowerPlusOne( ( 1.0 + (double)powerOfEnergy ) );
        double integralToStart( ( coefficientValue / energyPowerPlusOne ) );
        double integralToEnd( integralToStart );

        for( int powerCount( powerOfEnergy );
             0 <= powerCount;
             --powerCount )
        {
          integralToStart *= startEnergy;
          integralToEnd *= endEnergy;
        }
        /* now integralToStart is
         * coefficientValue * startEnergy^( 1 + powerOfEnergy )
         *                  / ( 1 + powerOfEnergy )
         * & similar for integralToEnd.
         */

        energyPowerPlusOne *= -1.0;
        // the terms in the sum alternate sign, so I put it here.
        double
        predivisorOfInversePowersOfEnergyPowerPlusOne( energyPowerPlusOne );
        double logarithmOfStartEnergy( log( startEnergy ) );
        double logarithmOfEndEnergy( log( endEnergy ) );

        int logPowerFactorial( 1 );
        double commonSumFactor( energyPowerPlusOne );
        double startSumTerm( energyPowerPlusOne * logarithmOfStartEnergy );
        double endSumTerm( energyPowerPlusOne * logarithmOfEndEnergy );
        double startSum( 1.0 + startSumTerm );
        double endSum( 1.0 + endSumTerm );
        /* now startSum is the sum of the terms corresponding to the lowest 2
         * powers of ln( startEnergy ) (i.e. 0 & 1, so the constant term & the
         * term linear in ln( startEnergy )), similar for endSum.
         * these need to be summed with the terms with powers of
         * ln( startEnergy ) from 2 to powerOfLogarithm.
         * logPowerFactorial = 1!, &
         * predivisorOfInversePowersOfEnergyPowerPlusOne already has its 1st
         * factor of energyPowerPlusOne. logPowerFactorial needs multiplying by
         * the integers 2 to powerOfLogarithm (to end up as powerOfLogarithm!)
         * & predivisorOfInversePowersOfEnergyPowerPlusOne needs to be
         * multiplied by energyPowerPlusOne ( powerOfLogarithm - 1 ) more
         * times (to end up as energyPowerPlusOne^powerOfLogarithm).
         * the following loop goes over the integers 2 to powerOfLogarithm,
         * hence ( powerOfLogarithm - 1 ) iterations.
         */
        for( int logarithmCounter( 2 );
             powerOfLogarithm >= logarithmCounter;
             ++logarithmCounter )
        {
          logPowerFactorial *= ( logPowerFactorial + 1 );
          predivisorOfInversePowersOfEnergyPowerPlusOne *= energyPowerPlusOne;

          commonSumFactor *= ( energyPowerPlusOne / (double)powerOfLogarithm );
          startSumTerm *= ( logarithmOfStartEnergy * commonSumFactor );
          startSum += startSumTerm;
          endSumTerm *= ( logarithmOfEndEnergy * commonSumFactor );
          endSum += endSumTerm;
        }
        integralToStart
        *= ( ( (double)logPowerFactorial
               / predivisorOfInversePowersOfEnergyPowerPlusOne )
             * startSum );
        integralToEnd
        *= ( ( (double)logPowerFactorial
               / predivisorOfInversePowersOfEnergyPowerPlusOne )
             * endSum );
        returnValue = ( integralToEnd - integralToStart );
      }
    }
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
    // does nothing.
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
                               double const tauMinEnergy
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
    << "hardMuonFromTau::operator()( " << muon_energy_power << ", "
    << muon_log_power << ", "
    << tau_energy_power << ", "
    << tau_log_power << ", "
    << Emin << ", "
    << Emax << ", "
    << is_inside_range
    << " ) called";
    std::cout << std::endl;**/

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
          returnValue = 1.0;
        }
      }
      else if( ( 0 == visibleProductEnergyPower )
          &&
          ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( -4.0 / ( 3.0 * tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue += ( 4.0 / ( 3.0 * tauMinEnergy ) );
        }
      }
      else if( ( 3 == visibleProductEnergyPower )
          &&
          ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( 1.0 / ( 3.0 * tauMaxEnergy * tauMaxEnergy
                                    * tauMaxEnergy * tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue -= ( 1.0 / ( 3.0 * tauMinEnergy * tauMinEnergy
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
        returnValue = ( ( 4.0 * log( tauMaxEnergy ) ) / 3.0 );
        if( isInsideRange )
        {
          returnValue -= ( 4.0 / 9.0 );
        }
        else
        {
          returnValue -= ( ( 4.0 * log( tauMinEnergy ) ) / 3.0 );
        }
      }
      else if( ( 0 == visibleProductEnergyPower )
               &&
               ( 1 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = ( -4.0 / 3.0 );
        }
      }
      else if( ( 3 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue
        = ( 4.0 / ( 9.0 * tauMaxEnergy * tauMaxEnergy * tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue
          -= ( 4.0 / ( 9.0 * tauMinEnergy * tauMinEnergy * tauMinEnergy ) );
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
        = ( ( 2.0 * logOfTauMaxEnergy * logOfTauMaxEnergy ) / 3.0 );
        if( isInsideRange )
        {
          returnValue -= ( 4.0 / 27.0 );
        }
        else
        {
          double logOfTauMinEnergy( log( tauMinEnergy ) );
          returnValue
          -= ( ( 2.0 * logOfTauMinEnergy * logOfTauMinEnergy ) / 3.0 );
        }
      }
      else if( ( 0 == visibleProductEnergyPower )
              &&
              ( 1 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = ( -4.0 / 9.0 );
        }
      }
      else if( ( 0 == visibleProductEnergyPower )
               &&
               ( 2 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = ( -2.0 / 3.0 );
        }
      }
      else if( ( 3 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue
        = ( ( 4.0 * ( 1.0 + 3.0 * log( tauMaxEnergy ) ) )
            / ( 27.0 * tauMaxEnergy * tauMaxEnergy * tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue
          -= ( ( 4.0 * ( 1.0 + 3.0 * log( tauMinEnergy ) ) )
               / ( 27.0 * tauMinEnergy * tauMinEnergy * tauMinEnergy ) );
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
        returnValue = ( ( 4.0 * tauMaxEnergy ) / 3.0 );
        if( !isInsideRange )
        {
          returnValue -= ( ( 4.0 * tauMinEnergy ) / 3.0 );
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = -2.0;
        }
      }
      else if( ( 3 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( 2.0 / ( 3.0 * tauMaxEnergy * tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue -= ( 2.0 / ( 3.0 * tauMinEnergy * tauMinEnergy ) );
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
        = ( ( 4.0 * tauMaxEnergy * ( log( tauMaxEnergy ) - 1.0 ) ) / 3.0 );
        if( !isInsideRange )
        {
          returnValue
          += ( ( 4.0 * tauMinEnergy * ( 1.0 - log( tauMinEnergy ) ) ) / 3.0 );
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
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
               ( 1 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = -2.0;
        }
      }
      else if( ( 3 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( ( 1.0 + 2.0 * log( tauMaxEnergy ) )
                        / ( 3.0 * tauMaxEnergy * tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue
          -= ( ( 1.0 + 2.0 * log( tauMinEnergy ) )
               / ( 3.0 * tauMinEnergy * tauMinEnergy ) );
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
        returnValue = ( ( 2.0 * tauMaxEnergy * tauMaxEnergy ) / 3.0 );
        if( !isInsideRange )
        {
          returnValue -= ( ( 2.0 * tauMinEnergy * tauMinEnergy ) / 3.0 );
        }
      }
      else if( ( 2 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = -2.0;
        }
      }
      else if( ( 3 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( 4.0 / ( 3.0 * tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue -= ( 4.0 / ( 3.0 * tauMinEnergy ) );
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
                               double const tauMinEnergy
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
          returnValue = 1.0;
        }
      }
      else if( ( 0 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( -2.0 / tauMaxEnergy );
        if( !isInsideRange )
        {
          returnValue += ( 2.0 / tauMinEnergy );
        }
      }
      else if( ( 2 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( 2.0 / ( tauMaxEnergy * tauMaxEnergy * tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue
          -= ( 2.0 / ( tauMinEnergy * tauMinEnergy * tauMinEnergy ) );
        }
      }
      else if( ( 3 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( -1.0 / ( tauMaxEnergy * tauMaxEnergy
                                 * tauMaxEnergy * tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue += ( 1.0 / ( tauMinEnergy * tauMinEnergy
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
        returnValue = ( 2.0 * log( tauMaxEnergy ) );
        if( isInsideRange )
        {
          returnValue -= ( 5.0 / 3.0 );
        }
        else
        {
          returnValue -= ( 2.0 * log( tauMinEnergy ) );
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
      else if( ( 2 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( 3.0 / ( tauMaxEnergy * tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue -= ( 3.0 / ( tauMinEnergy * tauMinEnergy ) );
        }
      }
      else if( ( 3 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue
        = ( -4.0 / ( 3.0 * tauMaxEnergy * tauMaxEnergy * tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue
          += ( 4.0 / ( 3.0 * tauMinEnergy * tauMinEnergy * tauMinEnergy ) );
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
        returnValue = ( logOfTauMaxEnergy * logOfTauMaxEnergy );
        if( isInsideRange )
        {
          returnValue -= ( 19.0 / 18.0 );
        }
        else
        {
          double logOfTauMinEnergy = log( tauMinEnergy );
          returnValue -= ( logOfTauMinEnergy * logOfTauMinEnergy );
        }
      }
      else if( ( 0 == visibleProductEnergyPower )
               &&
               ( 1 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = ( -5.0 / 3.0 );
        }
      }
      else if( ( 0 == visibleProductEnergyPower )
               &&
               ( 2 == visibleProductLogPower ) )
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
        returnValue = ( ( 1.5 + 3.0 * log( tauMaxEnergy ) )
                        / ( tauMaxEnergy * tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue -= ( ( 1.5 + 3.0 * log( tauMinEnergy ) )
                           / ( tauMinEnergy * tauMinEnergy ) );
        }
      }
      else if( ( 3 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( ( -4.0 * ( 1.0 + 3.0 * log( tauMaxEnergy ) ) )
                      / ( 9.0 * tauMaxEnergy * tauMaxEnergy * tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue += ( ( 4.0 * ( 1.0 + 3.0 * log( tauMinEnergy ) ) )
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
        returnValue = ( 2.0 * tauMaxEnergy );
        if( !isInsideRange )
        {
          returnValue -= ( 2.0 * tauMinEnergy );
        }
      }
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = -6.0;
        }
      }
      else if( ( 2 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( 6.0 / tauMaxEnergy );
        if( !isInsideRange )
        {
          returnValue -= ( 6.0 / tauMinEnergy );
        }
      }
      else if( ( 3 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( -2.0 / ( tauMaxEnergy * tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue += ( 2.0 / ( tauMinEnergy * tauMinEnergy ) );
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
        returnValue = ( 2.0 * tauMaxEnergy * ( log( tauMaxEnergy ) - 1.0 ) );
        if( !isInsideRange )
        {
          returnValue
          += ( 2.0 * tauMinEnergy * ( 1.0 - log( tauMinEnergy ) ) );
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
      else if( ( 1 == visibleProductEnergyPower )
               &&
               ( 1 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = -6.0;
        }
      }
      else if( ( 2 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue
        = ( ( 6.0 * ( 1.0 + log( tauMaxEnergy ) ) ) / tauMaxEnergy );
        if( !isInsideRange )
        {
          returnValue
          -= ( ( 6.0 * ( 1.0 + log( tauMinEnergy ) ) ) / tauMinEnergy );
        }
      }
      else if( ( 3 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( ( -1.0 * ( 1.0 + 2.0 * log( tauMaxEnergy ) ) )
                        / ( tauMaxEnergy * tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue += ( ( 1.0 + 2.0 * log( tauMinEnergy ) )
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
        returnValue = ( tauMaxEnergy * tauMaxEnergy );
        if( !isInsideRange )
        {
          returnValue -= ( tauMinEnergy * tauMinEnergy );
        }
      }
      else if( ( 2 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( -6.0 * log( tauMaxEnergy ) );
        if( isInsideRange )
        {
          returnValue += 3.0;
        }
        else
        {
          returnValue += ( 6.0 * log( tauMinEnergy ) );
        }
      }
      else if( ( 2 == visibleProductEnergyPower )
               &&
               ( 1 == visibleProductLogPower ) )
      {
        if( isInsideRange )
        {
          returnValue = 6.0;
        }
      }
      else if( ( 3 == visibleProductEnergyPower )
               &&
               ( 0 == visibleProductLogPower ) )
      {
        returnValue = ( -4.0 / ( tauMaxEnergy ) );
        if( !isInsideRange )
        {
          returnValue += ( 4.0 / ( tauMinEnergy ) );
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
  hardPionFromTau::operator()( int const pion_energy_power
                            /* the power of the pion's energy for the term. */,
                                   int const pion_log_power
           /* the power of the logarithm of the pion's energy for the term. */,
                                   int const tau_energy_power
                      /* the power of the tau lepton's energy for the term. */,
                                   int const tau_log_power
     /* the power of the logarithm of the tau lepton's energy for the term. */,
                                   double const Emin
         /* the minimum of the range of the tau lepton energies considered. */,
                                   double const Emax
         /* the maximum of the range of the tau lepton energies considered. */,
                                   bool const is_inside_range )
  const
  /* this returns the coefficientValue for the pion distribution which would come
   * from a right-handed tau lepton distribution of the given power of the
   * tau lepton's energy & power of logarithm thereof, for requested integer
   * powers of the pion's energy or logarithm thereof, which also depends on
   * whether the pion's energy is inside the range of the tau lepton's energy
   * for this segment.
   */
  {

    double return_value = 0.0;

    if( ( -1 == tau_energy_power )
        &&
        ( 0 == tau_log_power ) )
      {

        if( ( -1 == pion_energy_power )
            &&
            ( 0 == pion_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = 1.0;

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 0 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value
            = ( ( -1.0 * ( 1.0 + 2.0 * log_max ) ) / ( Emax * Emax ) );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value
                += ( ( ( 1.0 + 2.0 * log_min ) ) / ( Emin * Emin ) );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 1 == pion_log_power ) )
          {

            return_value = ( 2.0 / ( Emax * Emax ) );

            if( !is_inside_range )
              {

                return_value -= ( 2.0 / ( Emin * Emin ) );

              }

          }

      }  // end of if( -1 == tau_energy_power ) && ( 0 == tau_log_power )
    else if( ( 0 == tau_energy_power )
             &&
             ( 0 == tau_log_power ) )
      {

        if( ( 0 == pion_energy_power )
            &&
            ( 0 == pion_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = 4.0;

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 0 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( ( -4.0 * ( 1.0 + log_max ) ) / Emax );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value += ( ( 4.0 * ( 1.0 + log_min ) ) / Emin );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 1 == pion_log_power ) )
          {

            return_value = ( 4.0 / Emax );

            if( !is_inside_range )
              {

                return_value -= ( 4.0 / Emin );

              }

          }

      }  // end of if( 0 == tau_energy_power ) && ( 0 == tau_log_power )
    else if( ( 0 == tau_energy_power )
             &&
             ( 1 == tau_log_power ) )
      {

        if( ( 0 == pion_energy_power )
            &&
            ( 0 == pion_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = 8.0;

              }

          }
        else if( ( 0 == pion_energy_power )
                 &&
                 ( 1 == pion_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = 4.0;

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 0 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value
            = ( ( -4.0 * ( 2.0 + 2.0 * log_max + log_max * log_max ) )
                / Emax );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value
                += ( ( 4.0 * ( 2.0 + 2.0 * log_min + log_min * log_min ) )
                     / Emin );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 1 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( ( 4.0 * ( 1.0 + log_max ) ) / Emax );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value -= ( ( 4.0 * ( 1.0 + log_min ) ) / Emin );

              }

          }

      }  // end of if( 0 == tau_energy_power ) && ( 1 == tau_log_power )
    else if( ( 1 == tau_energy_power )
             &&
             ( 0 == tau_log_power ) )
      {

        if( ( 1 == pion_energy_power )
            &&
            ( 0 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( 2.0 * log_max * log_max );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value -= ( 2.0 * log_min * log_min );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 1 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( -4.0 * log_max );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value += ( 4.0 * log_min );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 2 == pion_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = 2.0;

              }

          }

      }  // end of if( 1 == tau_energy_power ) && ( 0 == tau_log_power )
    else if( ( 1 == tau_energy_power )
             &&
             ( 1 == tau_log_power ) )
      {

        if( ( 1 == pion_energy_power )
            &&
            ( 0 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( ( 4.0 * log_max * log_max * log_max ) / 3.0 );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value
                -= ( ( 4.0 * log_min * log_min * log_min ) / 3.0 );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 1 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( -2.0 * log_max * log_max );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value += ( 2.0 * log_min * log_min );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 3 == pion_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = ( 2.0 / 3.0 );

              }

          }

      }  // end of if( 1 == tau_energy_power ) && ( 1 == tau_log_power )
    else if( ( 2 == tau_energy_power )
             &&
             ( 0 == tau_log_power ) )
      {

        if( ( 1 == pion_energy_power )
            &&
            ( 0 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( 4.0 * Emax * ( log_max - 1.0 ) );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value += ( 4.0 * Emin * ( 1.0 - log_min ) );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 1 == pion_log_power ) )
          {

            return_value = ( -4.0 * Emax );

            if( !is_inside_range )
              {

                return_value += ( 4.0 * Emin );

              }

          }
        else if( ( 2 == pion_energy_power )
                 &&
                 ( 0 == pion_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = 4.0;

              }

          }

      }  // end of if( 2 == tau_energy_power ) && ( 0 == tau_log_power )

    return return_value;

  }


  soft_pion_from_tau::soft_pion_from_tau() :
    tauDecayCoefficient()
  {

    // just an initialization.

  }

  soft_pion_from_tau::~soft_pion_from_tau()
  {

    // does nothing.

  }


  double
  soft_pion_from_tau::operator()( int const pion_energy_power
                            /* the power of the pion's energy for the term. */,
                                   int const pion_log_power
           /* the power of the logarithm of the pion's energy for the term. */,
                                   int const tau_energy_power
                       /* the power of the tau lepton's energyfor the term. */,
                                   int const tau_log_power
     /* the power of the logarithm of the tau lepton's energy for the term. */,
                                   double const Emin
         /* the minimum of the range of the tau lepton energies considered. */,
                                   double const Emax
         /* the maximum of the range of the tau lepton energies considered. */,
                                   bool const is_inside_range )
  const
  /* this returns the coefficientValue for the pion distribution which would come
   * from a left-handed tau lepton distribution of the given power of the tau
   * lepton's energy & power of logarithm thereof, for requested integer
   * powers of the pion's energy or logarithm thereof, which also depends on
   * whether the pion's energy is inside the range of the tau lepton's energy
   * for this segment.
   */
  {

    double return_value = 0.0;

    if( ( -1 == tau_energy_power )
        &&
        ( 0 == tau_log_power ) )
      {

        if( ( -1 == pion_energy_power )
            &&
            ( 0 == pion_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = 1.0;

              }

          }
        else if( ( 0 == pion_energy_power )
                 &&
                 ( 0 == pion_log_power ) )
          {

            return_value = ( -4.0 / Emax );

            if( !is_inside_range )
              {

                return_value += ( 4.0 / Emin );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 0 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( ( 2.0 * log_max + 3.0 ) / ( Emax * Emax ) );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value -= ( ( 2.0 * log_min + 3.0 ) / ( Emin * Emin ) );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 1 == pion_log_power ) )
          {

            return_value = ( -2.0 / ( Emax * Emax ) );

            if( !is_inside_range )
              {

                return_value += ( 2.0 / ( Emin * Emin ) );

              }

          }

      }  // end of if( -1 == tau_energy_power ) && ( 0 == tau_log_power )
    else if( ( 0 == tau_energy_power )
             &&
             ( 0 == tau_log_power ) )
      {

        if( ( 0 == pion_energy_power )
            &&
            ( 0 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( 4.0 * log_max );

            if( is_inside_range )
              {

                return_value -= 8.0;

              }
            else
              {
                double log_min = log( Emin );

                return_value -= ( 4.0 * log_min );

              }

          }
        else if( ( 0 == pion_energy_power )
                 &&
                 ( 1 == pion_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = -4.0;

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 0 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( ( 4.0 * ( 2.0 + log_max ) ) / Emax );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value -= ( ( 4.0 * ( 2.0 + log_min ) ) / Emin );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 1 == pion_log_power ) )
          {

            return_value = ( -4.0 / Emax );

            if( !is_inside_range )
              {

                return_value += ( 4.0 / Emin );

              }

          }

      }  // end of if( 0 == tau_energy_power ) && ( 0 == tau_log_power )
    else if( ( 0 == tau_energy_power )
             &&
             ( 1 == tau_log_power ) )
      {

        if( ( 0 == pion_energy_power )
            &&
            ( 0 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( 2.0 * log_max * log_max );

            if( is_inside_range )
              {

                return_value -= 12.0;

              }
            else
              {

                double log_min = log( Emin );

                return_value -= ( 2.0 * log_min * log_min );

              }

          }
        else if( ( 0 == pion_energy_power )
                 &&
                 ( 1 == pion_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = -8.0;

              }

          }
        else if( ( 0 == pion_energy_power )
                 &&
                 ( 2 == pion_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = -2.0;

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 0 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value
            = ( ( 4.0 * ( 3.0 + 3.0 * log_max + log_max * log_max ) ) / Emax );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value
                -= ( ( 4.0 * ( 3.0 + 3.0 * log_min + log_min * log_min ) )
                     / Emin );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 1 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( ( -4.0 * ( 1.0 + log_max ) ) / Emax );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value += ( ( 4.0 * ( 1.0 + log_min ) ) / Emin );

              }

          }

      }  // end of if( 0 == tau_energy_power ) && ( 1 == tau_log_power )
    else if( ( 1 == tau_energy_power )
             &&
             ( 0 == tau_log_power ) )
      {

        if( ( 0 == pion_energy_power )
            &&
            ( 0 == pion_log_power ) )
          {

            return_value = ( 4.0 * Emax );

            if( !is_inside_range )
              {

                return_value -= ( 4.0 * Emin );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 0 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( -2.0 * log_max * ( log_max + 2.0 ) );

            if( is_inside_range )
              {

                return_value -= 4.0;

              }
            else
              {

                double log_min = log( Emin );

                return_value += ( 2.0 * log_min * ( log_min + 2.0 ) );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 1 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( 4.0 * log_max );

            if( is_inside_range )
              {

                return_value += 4.0;

              }
            else
              {

                double log_min = log( Emin );

                return_value -= ( 4.0 * log_min );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 2 == pion_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = -2.0;

              }

          }

      }  // end of if( 1 == tau_energy_power ) && ( 0 == tau_log_power )
    else if( ( 1 == tau_energy_power )
             &&
             ( 1 == tau_log_power ) )
      {

        if( ( 0 == pion_energy_power )
            &&
            ( 0 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( 4.0 * Emax * ( log_max - 1.0 ) );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value += ( 4.0 * Emin * ( 1.0 - log_min ) );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 0 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value
            = ( ( -2.0 * log_max * log_max * ( 3.0 + 2.0 * log_max ) ) / 3.0 );

            if( is_inside_range )
              {

                return_value += 4.0;

              }
            else
              {

                double log_min = log( Emin );

                return_value
                += ( ( 2.0 * log_min * log_min * ( 3.0 + 2.0 * log_min ) )
                     / 3.0 );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 1 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( 2.0 * log_max * log_max );

            if( is_inside_range )
              {

                return_value -= 4.0;

              }
            else
              {

                double log_min = log( Emin );

                return_value -= ( 2.0 * log_min * log_min );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 2 == pion_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = 2.0;

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 3 == pion_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = ( -2.0 / 3.0 );

              }

          }

      }  // end of if( 1 == tau_energy_power ) && ( 1 == tau_log_power )
    else if( ( 2 == tau_energy_power )
             &&
             ( 0 == tau_log_power ) )
      {

        if( ( 0 == pion_energy_power )
            &&
            ( 0 == pion_log_power ) )
          {

            return_value = ( 2.0 * Emax * Emax );

            if( !is_inside_range )
              {

                return_value -= ( 2.0 * Emin * Emin );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 0 == pion_log_power ) )
          {

            double log_max = log( Emax );

            return_value = ( -4.0 * Emax * log_max );

            if( !is_inside_range )
              {

                double log_min = log( Emin );

                return_value += ( 4.0 * Emin * log_min );

              }

          }
        else if( ( 1 == pion_energy_power )
                 &&
                 ( 1 == pion_log_power ) )
          {

            return_value = ( 4.0 * Emax );

            if( !is_inside_range )
              {

                return_value -= ( 4.0 * Emin );

              }

          }
        else if( ( 2 == pion_energy_power )
                 &&
                 ( 0 == pion_log_power ) )
          {

            if( is_inside_range )
              {

                return_value = -2.0;

              }

          }

      }  // end of if( 2 == tau_energy_power ) && ( 0 == tau_log_power )

    return return_value;

  }



  tau_segment_triple::tau_segment_triple(
                           std::vector< segmentTermSet* >* given_segment_set,
                                          segmentTermSet* given_segment,
                         tauDecayCoefficient const* const given_tau_decay ) :
    reference_segment( given_segment ),
    tau_decay( given_tau_decay )
  {

    double coefficient_checker;

    low_segment = new segmentTermSet();
    low_segment->setSegmentStart( 0.0 );

    high_segment = new segmentTermSet();

    for( std::vector< leptonDistributionExpansionTerm* >::const_iterator
         term_iterator = given_segment->getTerms()->begin();
         given_segment->getTerms()->end() > term_iterator;
         ++term_iterator )
      {

        for( int energy_power_counter = -1;
             3 >= energy_power_counter;
             ++energy_power_counter )
          {

            for( int log_power_counter = 0;
                 3 >= log_power_counter;
                 ++log_power_counter )
              {

                /* each term is checked to see if it's non-zero. "3.1" & "4.3"
                 * were chosen since there should be no accidental
                 * cancellations with these values, since there are just
                 * logarithms, powers of the end energies, & some rational
                 * numbers, but none are multiples of 3.1 or 4.3 so there can
                 * be no cancellation.
                 */

                coefficient_checker
                = (*tau_decay)( energy_power_counter,
                                log_power_counter,
                                (*term_iterator)->getPowerOfEnergy(),
                                (*term_iterator)->getPowerOfLogarithm(),
                                3.1,
                                4.3,
                                false );

                if( 0.0 != coefficient_checker )
                  {

                    low_segment->addTerm(
                                        new leptonDistributionExpansionTerm(
                                                          energy_power_counter,
                                                             log_power_counter,
                                   CppSLHA::CppSLHA_global::really_wrong_value,
                                                            *term_iterator ) );

                  }

                coefficient_checker
                = (*tau_decay)( energy_power_counter,
                                log_power_counter,
                                (*term_iterator)->getPowerOfEnergy(),
                                (*term_iterator)->getPowerOfLogarithm(),
                                3.1,
                                4.3,
                                true );

                if( 0.0 != coefficient_checker )
                  {

                    high_segment->addTerm(
                                        new leptonDistributionExpansionTerm(
                                                          energy_power_counter,
                                                             log_power_counter,
                                   CppSLHA::CppSLHA_global::really_wrong_value,
                                                            *term_iterator ) );

                  }

              }  // end of loop over powers of log( energy ).

          }  // end of loop over powers of energy.

      }  // end of loop over leptonDistributionExpansionTerm pointers.

    // now that all the low_segment & high_segment terms have been added to the
    // segments, we add them to the vector of segment_term_sets.
    given_segment_set->push_back( low_segment );
    given_segment_set->push_back( high_segment );

  }

  tau_segment_triple::~tau_segment_triple()
  {

    // does nothing.

  }


  void
  tau_segment_triple::update_segments()
  {

    low_segment->setSegmentEnd( reference_segment->getSegmentStart() );
    high_segment->setSegmentRange( reference_segment->getSegmentStart(),
                                     reference_segment->getSegmentEnd() );

    for( std::vector< leptonDistributionExpansionTerm* >::const_iterator
         term_iterator( low_segment->getTerms()->begin() );
         low_segment->getTerms()->end() > term_iterator;
         ++term_iterator )
      {

        (*term_iterator)->setCoefficient( reference_segment->getCoefficient(
                                       (*term_iterator)->getPowerOfEnergy(),
                                   (*term_iterator)->getPowerOfLogarithm() )
                       * (*tau_decay)( (*term_iterator)->getPowerOfEnergy(),
                                    (*term_iterator)->getPowerOfLogarithm(),
                                       (*term_iterator)->getReferenceTerm(
                                                      )->getPowerOfEnergy(),
                                       (*term_iterator)->getReferenceTerm(
                                                   )->getPowerOfLogarithm(),
                                       0.0,
                                       low_segment->getSegmentEnd(),
                                       false ) );

      }
    for( std::vector< leptonDistributionExpansionTerm* >::const_iterator
         term_iterator( high_segment->getTerms()->begin() );
         high_segment->getTerms()->end() > term_iterator;
         ++term_iterator )
      {

        (*term_iterator)->setCoefficient( reference_segment->getCoefficient(
                                       (*term_iterator)->getPowerOfEnergy(),
                                   (*term_iterator)->getPowerOfLogarithm() )
                       * (*tau_decay)( (*term_iterator)->getPowerOfEnergy(),
                                    (*term_iterator)->getPowerOfLogarithm(),
                                       (*term_iterator)->getReferenceTerm(
                                                      )->getPowerOfEnergy(),
                                       (*term_iterator)->getReferenceTerm(
                                                   )->getPowerOfLogarithm(),
                                       high_segment->getSegmentStart(),
                                       high_segment->getSegmentEnd(),
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
    for( std::vector< segmentTermSet* >::iterator
         deletionIterator = segments.begin();
         segments.end() > deletionIterator;
         ++deletionIterator )
      {
        delete *deletionIterator;
      }
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
  }



  visible_tau_decay_product::visible_tau_decay_product(
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
    tau_distribution( given_tau_distribution )
  {

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: "
    << "visible_tau_decay_product( "
    << tau_distribution << ", " << tau_decay << " )";
    std::cout
    << std::endl;**/

    minimumEnergy = 0.0;

    for( std::vector< segmentTermSet* >::iterator
         segment_iterator = given_tau_distribution->getSegments()->begin();
         given_tau_distribution->getSegments()->end() > segment_iterator;
         ++segment_iterator )
      {

        tau_triples.push_back( new tau_segment_triple( &segments,
                                                       *segment_iterator,
                                                       tau_decay ) );

      }  // end of loop through segments.

  }

  visible_tau_decay_product::~visible_tau_decay_product()
  {

    for( std::vector< tau_segment_triple* >::iterator
         deletion_iterator = tau_triples.begin();
         tau_triples.end() > deletion_iterator;
         ++deletion_iterator )
      {

        delete *deletion_iterator;

      }

  }


  void
  visible_tau_decay_product::calculateCoefficients()
  {

    // debugging:
    /**std::cout
    << std::endl
    << "debugging: "
    << "visible_tau_decay_product::calculate_coefficients()";
    std::cout
    << std::endl;**/

    maximumEnergy = tau_distribution->getMaximumEnergy();
    // now all the underlying stuff should have been updated for the point.

    for( std::vector< tau_segment_triple* >::iterator
         triple_iterator( tau_triples.begin() );
         tau_triples.end() > triple_iterator;
         ++triple_iterator )
      {

        (*triple_iterator)->update_segments();

      }  // end of loop through segments.

  }

}  // end of LHC_FASER namespace

