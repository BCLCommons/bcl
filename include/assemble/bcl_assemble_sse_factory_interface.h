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

#ifndef BCL_ASSEMBLE_SSE_FACTORY_INTERFACE_H_
#define BCL_ASSEMBLE_SSE_FACTORY_INTERFACE_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "sspred/bcl_sspred.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEFactoryInterface
    //! @brief Interface class to be used for deriving SSEFactories,for generating set of SSEs from given AASequences.
    //! @details SSEFactoryInterface provides the interface for SSEFactory classes that use various predictions, constraints,
    //! etc. to generate/predict set of SSEs from the given amino acid sequences.
    //!
    //! @remarks example unnecessary
    //! @author karakam
    //! @date 22.02.2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEFactoryInterface :
      public math::FunctionInterfaceSerializable< biol::AASequence, SSEPool>
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief virtual copy constructor
      virtual SSEFactoryInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief set the ss method used
      //! @param SS_METHOD sspred method
      virtual void SetMethod( const sspred::Method &SS_METHOD) = 0;

      //! @brief set the thresholds to use
      //! @brief SSTYPE_THRESHOLDS the thresholds to use for the sstypes desired
      virtual void SetThresholds( const storage::Map< biol::SSType, double> &SSTYPE_THRESHOLDS) = 0;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that returns a set of SSEs for the given AASequence
      //! @brief SEQUENCE AASequence from which the SSEPool is going to be built
      //! @return SSEPool built from provided SEQUENCE
      virtual SSEPool
      operator()( const biol::AASequence &SEQUENCE) const = 0;

      //! @brief operator that returns a set of SSEs for the given vector of AASequences
      //! @brief SEQUENCES Vector of AASequences from which the SSEPool is going to be built
      //! @return SSEPool built from provided SEQUENCES
      SSEPool
      operator()( const util::SiPtrVector< biol::AASequence> &SEQUENCES) const;

    }; // class SSEFactoryInterface

  } // namespace assemble
} // namespace bcl

#endif //BCL_ASSEMBLE_SSE_FACTORY_INTERFACE_H_
