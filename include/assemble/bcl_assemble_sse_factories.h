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

#ifndef BCL_ASSEMBLE_SSE_FACTORIES_H_
#define BCL_ASSEMBLE_SSE_FACTORIES_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_assemble_sse_factory_interface.h"
#include "util/bcl_util_enumerate.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEFactories
    //! @brief TODO: add a brief comment
    //! @details TODO: add an detailed description to this class
    //!
    //! @see @link example_assemble_sse_factories.cpp @endlink
    //! @author woetzen
    //! @date Jul 2, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEFactories :
      public util::Enumerate< util::ShPtr< SSEFactoryInterface>, SSEFactories>
    {
        friend class util::Enumerate< util::ShPtr< SSEFactoryInterface>, SSEFactories>;
    public:

    //////////
    // data //
    //////////

      const SSEFactory e_Conformation;        //!< from backbone conformation
      const SSEFactory e_PredictionMC;        //!< from prediction optimized with MCM
      const SSEFactory e_PredictionThreshold; //!< from prediction applying a threshold for each sstype
      const SSEFactory e_PredictionHighest;   //!< from prediction - highest prediction and if above threshold

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SSEFactories();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief create SSEFactoryInterface from enum, threshold and method
      //! @param FACTORY enum of the factory to be used
      //! @param SS_METHOD ss_pred to be used
      //! @param SSTYPE_THRESHOLDS sstype thresholds
      //! @return ShPtr to assemble::SSEFactoryInterface
      util::ShPtr< SSEFactoryInterface> Create
      (
        const SSEFactory &FACTORY,
        const sspred::Method &SS_METHOD,
        const storage::Map< biol::SSType, double> &SSTYPE_THRESHOLDS
      ) const;

    }; // class SSEFactories

    //! @brief construct on access function for all SSEFactories
    //! @return reference to only instances of SSEFactories
    BCL_API
    const SSEFactories &GetSSEFactories();

  } // namespace assemble

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< assemble::SSEFactoryInterface>, assemble::SSEFactories>;

  } // namespace util
} // namespace bcl

#endif // BCL_ASSEMBLE_SSE_FACTORIES_H_ 
