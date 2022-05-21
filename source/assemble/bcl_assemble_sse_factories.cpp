// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "assemble/bcl_assemble_sse_factories.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse_factory_conformation.h"
#include "assemble/bcl_assemble_sse_factory_mc.h"
#include "sspred/bcl_sspred_sse_factory_highest.h"
#include "sspred/bcl_sspred_sse_factory_threshold.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SSEFactories::SSEFactories() :
      e_Conformation(        AddEnum( "Conformation"   , util::ShPtr< assemble::SSEFactoryInterface>( new         SSEFactoryConformation()))),
      e_PredictionMC(        AddEnum( "SSPredMC"       , util::ShPtr< assemble::SSEFactoryInterface>( new         SSEFactoryMC( sspred::GetMethods().e_Undefined, 0.0)))),
      e_PredictionThreshold( AddEnum( "SSPredThreshold", util::ShPtr< assemble::SSEFactoryInterface>( new sspred::SSEFactoryThreshold( sspred::GetMethods().e_Undefined)))),
      e_PredictionHighest(   AddEnum( "SSPredHighest"  , util::ShPtr< assemble::SSEFactoryInterface>( new sspred::SSEFactoryHighest( sspred::GetMethods().e_Undefined))))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSEFactories::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief create SSEFactoryInterface from enum, threshold and method
    //! @param FACTORY enum of the factory to be used
    //! @param SS_METHOD ss_pred to be used
    //! @param SSTYPE_THRESHOLDS sstype thresholds
    //! @return ShPtr to assemble::SSEFactoryInterface
    util::ShPtr< SSEFactoryInterface> SSEFactories::Create
    (
      const SSEFactory &FACTORY,
      const sspred::Method &SS_METHOD,
      const storage::Map< biol::SSType, double> &SSTYPE_THRESHOLDS
    ) const
    {
      // create factory
      util::ShPtr< SSEFactoryInterface> sp_factory( ( *FACTORY)->Clone());

      // set the member
      sp_factory->SetMethod( SS_METHOD);
      sp_factory->SetThresholds( SSTYPE_THRESHOLDS);

      // end
      return sp_factory;
    }

    //! @brief construct on access function for all SSEFactories
    //! @return reference to only instances of SSEFactories
    const SSEFactories &GetSSEFactories()
    {
      return SSEFactories::GetEnums();
    }

  } // namespace assemble

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< assemble::SSEFactoryInterface>, assemble::SSEFactories>;

  } // namespace util
} // namespace bcl
