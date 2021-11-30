/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016,2017 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "PathProjectionCalculator.h"
#include "core/ActionWithValue.h"
#include "core/ActionRegister.h"
#include "core/ActionSet.h"
#include <string.h>

namespace PLMD {
namespace mapping {

void PathProjectionCalculator::registerKeywords(Keywords& keys) {
  keys.add("compulsory","METRIC","the method to use for computing the displacement vectors between the reference frames");
  keys.add("compulsory","REFERENCE","labels for actions that contain reference coordinates for each point on the path");
}

PathProjectionCalculator::PathProjectionCalculator( Action* act ):
mypath_obj(NULL)
{
  ActionWithArguments* aarg=dynamic_cast<ActionWithArguments*>( act );
  if( aarg ) {
      mypath_obj = aarg->getPntrToArgument(0);
      // Check that we have only one argument as input
      if( aarg->getNumberOfArguments()!=1 ) act->error("should only have one argument to this function");
  }
  // Ensure that values are stored in base calculation and that PLUMED doesn't try to calculate this in the stream
  if( mypath_obj ) mypath_obj->buildDataStore( act->getLabel() );
  // Check that the input is a matrix
  if( mypath_obj ) if( mypath_obj->getRank()!=2 ) act->error("the input to this action should be a matrix");
  // Get the labels for the reference points
  std::vector<std::string> reference_data; act->parseVector("REFERENCE", reference_data);
  for(unsigned i=0; i<reference_data.size(); ++i) {
      ActionWithValue* av = act->plumed.getActionSet().selectWithLabel<ActionWithValue*>( reference_data[i] );
      if( !av || av->getName()!="CONSTANT_VALUE" ) act->error("input " + reference_data[i] + " is not a CONSTANT_VALUE action");
      Value* myval = av->copyOutput(0); refargs.push_back( myval );
      if( myval->getRank()==2 && reference_data.size()!=1 ) act->error("should only be one matrix in input to path projection object");
      if( myval->getRank()>0 && myval->getShape()[0]!=refargs[0]->getShape()[0] ) act->error("mismatch in number of reference frames in input to reference_data"); 
  }
  // Create a plumed main object to compute distances between reference configurations
  int s=sizeof(double);
  metric.cmd("setRealPrecision",&s);
  metric.cmd("setNoVirial");
  metric.cmd("setMDEngine","plumed");
  int nat=0; metric.cmd("setNatoms",&nat);
  unsigned nargs=refargs.size(); if( refargs[0]->getRank()==2 ) nargs = refargs[0]->getShape()[1];
  std::string str_nargs; Tools::convert( nargs, str_nargs ); std::string period_str=" PERIODIC=NO";
  if( mypath_obj && mypath_obj->isPeriodic() ) { std::string min, max; mypath_obj->getDomain( min, max ); period_str=" PERIODIC=" + min + "," + max; }
  metric.cmd("createValue arg1: PUT SHAPE=" + str_nargs + period_str);
  metric.cmd("createValue arg2: PUT SHAPE=" + str_nargs + period_str);
  double tstep=1.0; metric.cmd("setTimestep",&tstep);
  std::string inp; act->parse("METRIC",inp); inp += " ARG1=arg2 ARG2=arg1"; const char* cinp=inp.c_str();
  std::vector<std::string> input=Tools::getWords(inp);
  if( input.size()==1 && !actionRegister().check(input[0]) ) {
      metric.cmd("setPlumedDat",cinp); metric.cmd("init");
  } else {
      metric.cmd("init"); metric.cmd("readInputLine",cinp);
  }
  // Now setup stuff to retrieve the final displacement
  unsigned aind = metric.getActionSet().size()-1;
  while( true ) {
     const ActionShortcut* as=dynamic_cast<const ActionShortcut*>( metric.getActionSet()[aind].get() );
     if( !as ) break ; aind = aind - 1; plumed_assert( aind>=0 );
  } 
  ActionWithValue* fav = dynamic_cast<ActionWithValue*>( metric.getActionSet()[aind].get() );
  if( !fav ) act->error("final value should calculate relevant value that you want as reference");
  std::string name = (fav->copyOutput(0))->getName(); long rank; metric.cmd("getDataRank " + name, &rank );
  if( rank==0 ) rank=1;
  std::vector<long> ishape( rank ); metric.cmd("getDataShape " + name, &ishape[0] );
  unsigned nvals=1; for(unsigned i=0;i<ishape.size();++i) nvals *= ishape[i];
  data.resize( nvals ); metric.cmd("setMemoryForData " + name, &data[0] );

}

unsigned PathProjectionCalculator::getNumberOfFrames() const {
  return refargs[0]->getShape()[0];
}

void PathProjectionCalculator::computeVectorBetweenFrames( const unsigned& ifrom, const unsigned& ito ) {
  int step = 1; metric.cmd("setStep",&step); 
  std::vector<double> valdata1( data.size() ), valdata2( data.size() );
  getReferenceConfiguration( ito, valdata2 ); getReferenceConfiguration( ifrom, valdata1 );
  metric.cmd("setValue arg1", &valdata1[0] ); 
  metric.cmd("setValue arg2", &valdata2[0] );
  metric.cmd("calc");
}

void PathProjectionCalculator::getDisplaceVector( const unsigned& ifrom, const unsigned& ito, std::vector<double>& displace ) {
  if( displace.size()!=data.size() ) displace.resize( data.size() );
  computeVectorBetweenFrames( ifrom, ito ); for(unsigned i=0;i<data.size();++i) displace[i] = data[i]; 
}

void PathProjectionCalculator::getReferenceConfiguration( const unsigned& iframe, std::vector<double>& refpos ) const {
  if( refpos.size()!=data.size() ) refpos.resize( data.size() );
  if( refargs[0]->getRank()==2 ) {
      for(unsigned i=0; i<refpos.size(); ++i) refpos[i] = refargs[0]->get( iframe*refpos.size() + i ); 
  } else {
      for(unsigned i=0; i<refpos.size(); ++i) refpos[i] = refargs[i]->get(iframe);
  }
} 

void PathProjectionCalculator::setReferenceConfiguration( const unsigned& iframe, std::vector<double>& refpos ) {
  plumed_dbg_assert( refpos.size()==data.size() );
  if( refargs[0]->getRank()==2 ) {
      for(unsigned i=0; i<refpos.size(); ++i) refargs[0]->set( iframe*refpos.size() + i, refpos[i] ); 
  } else {
      for(unsigned i=0; i<refpos.size(); ++i) refargs[i]->set( iframe, refpos[i] );
  }
} 

}
}

