/*
 * LHC-FASER_kinematics_stuff.hpp
 *
 *  Created on: 02 Mar 2011
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

#include "LHC-FASER_jet_kinematics_stuff.hpp"

namespace LHC_FASER
{
  squarkMassForGridDecider::squarkMassForGridDecider(
                                         inputHandler const* const shortcut ) :
    shortcut( shortcut )
  {
    // just an initialization list.
  }

  squarkMassForGridDecider::~squarkMassForGridDecider()
  {
    // does nothing.
  }


  heavierThanGluinoSquarkMassForGrid::heavierThanGluinoSquarkMassForGrid(
                                         inputHandler const* const shortcut ) :
    squarkMassForGridDecider( shortcut )
  {
    // just an initialization list.
  }

  heavierThanGluinoSquarkMassForGrid::~heavierThanGluinoSquarkMassForGrid()
  {
    // does nothing.
  }


  firstMassForGrid::firstMassForGrid( inputHandler const* const shortcut ) :
    squarkMassForGridDecider( shortcut )
  {
    // just an initialization list.
  }

  firstMassForGrid::~firstMassForGrid()
  {
    // does nothing.
  }


  secondMassForGrid::secondMassForGrid( inputHandler const* const shortcut ) :
    squarkMassForGridDecider( shortcut )
  {
    // just an initialization list.
  }

  secondMassForGrid::~secondMassForGrid()
  {
    // does nothing.
  }


  averageSquarkMassForGrid::averageSquarkMassForGrid(
                                         inputHandler const* const shortcut ) :
    squarkMassForGridDecider( shortcut )
  {
    // just an initialization list.
  }

  averageSquarkMassForGrid::~averageSquarkMassForGrid()
  {
    // does nothing.
  }



  jetAcceptanceTable::usedCascades const
  jetAcceptanceTable::lastEnumElement( jetAcceptanceTable::sizeOfEnum );
  int const
  jetAcceptanceTable::sizeOfUsedCascades(
                                    (int)jetAcceptanceTable::lastEnumElement );

  jetAcceptanceTable::jetAcceptanceTable(
                                    std::string const* const gridFilesLocation,
                                          std::string const* const jetCutName,
                                          int const acceptanceColumn,
                                         inputHandler const* const shortcut ) :
    jetCutName( *jetCutName ),
    acceptanceColumn( acceptanceColumn ),
    shortcut( shortcut ),
    heavierThanGluinoSquarkMass( shortcut ),
    useFirstMass( shortcut ),
    useSecondMass( shortcut ),
    useAverageMass( shortcut )
  {
    std::string gridFileBaseName( *gridFilesLocation );
    gridFileBaseName.append( "/" );
    gridFileBaseName.append( *jetCutName );
    std::string gridFileName( gridFileBaseName );
    gridFileName.append( "/gluino+gluino_acceptance.dat" );
    gluinoGluinoGrid = new acceptanceGrid( &gridFileName );
    gridFileName.assign( gridFileBaseName );
    gridFileName.append( "/squark+gluino_acceptance.dat" );
    squarkGluinoGrid = new acceptanceGrid( &gridFileName );
    gridFileName.assign( gridFileBaseName );
    gridFileName.append( "/squark+antisquark_acceptance.dat" );
    squarkAntisquarkGrid = new acceptanceGrid( &gridFileName );
    gridFileName.assign( gridFileBaseName );
    gridFileName.append( "/squark+squark_acceptance.dat" );
    squarkSquarkGrid = new acceptanceGrid( &gridFileName );
    for( int vectorAdder( sizeOfUsedCascades );
         0 < vectorAdder;
         --vectorAdder )
    {
      gridMatrixRow = gridsMatrix.addNewAtEnd();
      for( int pairAdder( sizeOfUsedCascades );
           0 < pairAdder;
           --pairAdder )
      {
        gridMatrixRow->addNewAtEnd();
      }
    }

    // now we specify the various cases:
    gridMatrixRow = gridsMatrix.getPointer( (int)gx );
    gridMatrixElement = gridMatrixRow->getPointer( (int)gx );
    gridMatrixElement->first = gluinoGluinoGrid;
    gridMatrixElement->second = &heavierThanGluinoSquarkMass;
    gridMatrixElement = gridMatrixRow->getPointer( (int)sx );
    gridMatrixElement->first = squarkGluinoGrid;
    gridMatrixElement->second = &useSecondMass;
    gridMatrixElement = gridMatrixRow->getPointer( (int)gsx );
    gridMatrixElement->first = gluinoGluinoGrid;
    gridMatrixElement->second = &useSecondMass;
    gridMatrixElement = gridMatrixRow->getPointer( (int)sgx );
    gridMatrixElement->first = squarkGluinoGrid;
    gridMatrixElement->second = &useSecondMass;

    gridMatrixRow = gridsMatrix.getPointer( (int)sx );
    gridMatrixElement = gridMatrixRow->getPointer( (int)gx );
    gridMatrixElement->first = squarkGluinoGrid;
    gridMatrixElement->second = &useFirstMass;
    gridMatrixElement = gridMatrixRow->getPointer( (int)sx );
    gridMatrixElement->first = squarkAntisquarkGrid;
    gridMatrixElement->second = &useAverageMass;
    gridMatrixElement = gridMatrixRow->getPointer( (int)gsx );
    gridMatrixElement->first = squarkGluinoGrid;
    gridMatrixElement->second = &useAverageMass;
    gridMatrixElement = gridMatrixRow->getPointer( (int)sgx );
    gridMatrixElement->first = squarkAntisquarkGrid;
    gridMatrixElement->second = &useFirstMass;

    gridMatrixRow = gridsMatrix.getPointer( (int)gsx );
    gridMatrixElement = gridMatrixRow->getPointer( (int)gx );
    gridMatrixElement->first = gluinoGluinoGrid;
    gridMatrixElement->second = &useFirstMass;
    gridMatrixElement = gridMatrixRow->getPointer( (int)sx );
    gridMatrixElement->first = squarkGluinoGrid;
    gridMatrixElement->second = &useAverageMass;
    gridMatrixElement = gridMatrixRow->getPointer( (int)gsx );
    gridMatrixElement->first = gluinoGluinoGrid;
    gridMatrixElement->second = &useAverageMass;
    gridMatrixElement = gridMatrixRow->getPointer( (int)sgx );
    gridMatrixElement->first = gluinoGluinoGrid;
    gridMatrixElement->second = &useFirstMass;

    gridMatrixRow = gridsMatrix.getPointer( (int)sgx );
    gridMatrixElement = gridMatrixRow->getPointer( (int)gx );
    gridMatrixElement->first = squarkGluinoGrid;
    gridMatrixElement->second = &useFirstMass;
    gridMatrixElement = gridMatrixRow->getPointer( (int)sx );
    gridMatrixElement->first = squarkAntisquarkGrid;
    gridMatrixElement->second = &useSecondMass;
    gridMatrixElement = gridMatrixRow->getPointer( (int)gsx );
    gridMatrixElement->first = gluinoGluinoGrid;
    gridMatrixElement->second = &useSecondMass;
    gridMatrixElement = gridMatrixRow->getPointer( (int)sgx );
    gridMatrixElement->first = squarkAntisquarkGrid;
    gridMatrixElement->second = &useAverageMass;
  }

  jetAcceptanceTable::~jetAcceptanceTable()
  {
    delete gluinoGluinoGrid;
    delete squarkGluinoGrid;
    delete squarkAntisquarkGrid;
    delete squarkSquarkGrid;
  }

  int
  jetAcceptanceTable::getIntForCascadeType(
                                         fullCascade const* const givenCascade,
                                     particlePointer* const squarkFromCascade )
  {
    usedCascades typeToCountAs;
    if( fullCascade::gx == givenCascade->getColofulCascadeType() )
    {
      typeToCountAs = gx;
      *squarkFromCascade = NULL;
    }
    else if( ( fullCascade::sx == givenCascade->getColofulCascadeType() )
             ||
             ( fullCascade::svsx == givenCascade->getColofulCascadeType() ) )
    {
      typeToCountAs = sx;
      *squarkFromCascade = givenCascade->getCascadeDefiner()->back()->first;
    }
    else if( ( fullCascade::sjgx == givenCascade->getColofulCascadeType() )
             ||
             ( fullCascade::svgx == givenCascade->getColofulCascadeType() )
             ||
             ( fullCascade::svsjgx == givenCascade->getColofulCascadeType() ) )
    {
      typeToCountAs = sgx;
      *squarkFromCascade = givenCascade->getCascadeDefiner()->at( 2 )->first;
    }
    else
    {
      typeToCountAs = gsx;
      *squarkFromCascade = givenCascade->getCascadeDefiner()->at( 1 )->first;
    }
    return (int)typeToCountAs;
  }

  double
  jetAcceptanceTable::getAcceptance(
                        signedParticleShortcutPair const* const initialPair,
                                     fullCascade const* const firstCascade,
                                     fullCascade const* const secondCascade )
  {
    gridMatrixRow = gridsMatrix.getPointer( getIntForCascadeType( firstCascade,
                                                       &firstCascadeSquark ) );
    gridMatrixElement = gridMatrixRow->getPointer( getIntForCascadeType(
                                                                 secondCascade,
                                                      &secondCascadeSquark ) );
    /* gridMatrixElement is now a pair of pointers, 1 to the appropriate grid,
     * the other to the appropriate functor for obtaining the squark mass to
     * use. now we check to see if we need to redirect to squarkSquarkGrid:
     */
    gridToUse = gridMatrixElement->first;
    if( ( squarkAntisquarkGrid == gridToUse )
        &&
        ( initialPair->firstIsNotAntiparticle()
          == initialPair->secondIsNotAntiparticle() ) )
    {
      gridToUse = squarkSquarkGrid;
    }
    // now we just use the general valueAt(...) function of acceptanceGrid:
    return
    gridToUse->valueAt( (*(gridMatrixElement->second))( firstCascadeSquark,
                                                        secondCascadeSquark ),
                        shortcut->getGluinoMass(),
                        firstCascade->getCascadeDefiner()->front(
                                                 )->first->get_absolute_mass(),
                        secondCascade->getCascadeDefiner()->front(
                                                 )->first->get_absolute_mass(),
                        acceptanceColumn,
                        false,
                        false );
  }



  jetAcceptanceTablesForOneBeamEnergy::jetAcceptanceTablesForOneBeamEnergy(
                                           inputHandler const* const shortcut,
                                                          int const beamEnergy,
                               std::string const* const gridFileSetLocation ) :
    shortcut( shortcut ),
    beamEnergy( beamEnergy ),
    gridFileSetLocation( *gridFileSetLocation )
  {
    std::stringstream energyStream( "" );
    energyStream << "/" << beamEnergy << "TeV/jets/";
    this->gridFileSetLocation.append( energyStream.str() );
  }

  jetAcceptanceTablesForOneBeamEnergy::~jetAcceptanceTablesForOneBeamEnergy()
  {
    for( std::vector< jetAcceptanceTable* >::iterator
         deletionIterator( jetTables.begin() );
         jetTables.end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
  }

  jetAcceptanceTable*
  jetAcceptanceTablesForOneBeamEnergy::getTable(
                                           std::string const* const jetCutName,
                                                 int const acceptanceColumn )
  // this returns the jetAcceptanceTable for the requested signal & column.
  {
    jetAcceptanceTable* returnPointer( NULL );
    // we look to see if we already have a jetAcceptanceTable for this jet+MET
    // cut set:
    for( std::vector< jetAcceptanceTable* >::iterator
         searchIterator( jetTables.begin() );
         jetTables.end() > searchIterator;
         ++searchIterator )
    {
      if((*searchIterator)->isRequested( jetCutName,
                                         acceptanceColumn ) )
      {
        returnPointer = *searchIterator;
        searchIterator = jetTables.end();
      }
    }
    if( NULL == returnPointer )
      // if we do not already have a jetAcceptanceTable for this jet+MET cut
      // set, we make a new instance:
    {
      returnPointer = new jetAcceptanceTable( &gridFileSetLocation,
                                              jetCutName,
                                              acceptanceColumn,
                                              shortcut );
      jetTables.push_back( returnPointer );
    }
    return returnPointer;
  }



  jetPlusMetAcceptanceHandler::jetPlusMetAcceptanceHandler(
                                           inputHandler const* const shortcut,
                               std::string const* const gridFileSetLocation ) :
    shortcut( shortcut ),
    gridFileSetLocation( *gridFileSetLocation )
  {
    // just an initialization list.
  }

  jetPlusMetAcceptanceHandler::~jetPlusMetAcceptanceHandler()
  {
    for( std::vector< jetAcceptanceTablesForOneBeamEnergy* >::iterator
         deletionIterator( acceptanceTables.begin() );
         acceptanceTables.end() > deletionIterator;
         ++deletionIterator )
    {
      delete *deletionIterator;
    }
  }

  jetAcceptanceTablesForOneBeamEnergy*
  jetPlusMetAcceptanceHandler::getJetAcceptanceTablesForOneBeamEnergy(
                                                         int const beamEnergy )
  /* this looks to see if there is an existing
   * jetAcceptanceTablesForOneBeamEnergy with the requested beamEnergy, & if
   * not, makes 1, & returns the pointer.
   */
  {
    jetAcceptanceTablesForOneBeamEnergy* returnPointer( NULL );
    // we look to see if we already have a leptonAcceptancesForOneBeamEnergy
    // for these values:
    for( std::vector< jetAcceptanceTablesForOneBeamEnergy* >::iterator
         searchIterator( acceptanceTables.begin() );
         acceptanceTables.end() > searchIterator;
         ++searchIterator )
    {
      if( beamEnergy == (*searchIterator)->getBeamEnergy() )
      {
        returnPointer = *searchIterator;
        searchIterator = acceptanceTables.end();
      }
    }
    if( NULL == returnPointer )
      // if we do not already have a leptonAcceptanceParameterSet for this
      // electroweakino, we make a new instance:
    {
      returnPointer = new jetAcceptanceTablesForOneBeamEnergy( shortcut,
                                                               beamEnergy,
                                                        &gridFileSetLocation );
      acceptanceTables.push_back( returnPointer );
    }
    return returnPointer;
  }

}  // end of LHC_FASER namespace.

