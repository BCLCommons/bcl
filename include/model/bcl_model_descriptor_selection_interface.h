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

#ifndef BCL_MODEL_DESCRIPTOR_SELECTION_INTERFACE_H_
#define BCL_MODEL_DESCRIPTOR_SELECTION_INTERFACE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DescriptorSelectionInterface
    //! @brief interface to generalize all descriptor selection schemata
    //!
    //! @see @link example_namespace_descriptor_selection_interface.cpp @endlink
    //! @author butkiem1
    //! @date May 31, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DescriptorSelectionInterface :
      public util::SerializableInterface
    {

    private:

    //////////
    // data //
    //////////

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new Interface
      virtual DescriptorSelectionInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief get default initialized descriptor set
      //! @param TOTAL descriptor set with all available descriptor groups
      //! @return initialized descriptor set
      virtual util::ObjectDataLabel GetInitialDescriptorSet( const util::ObjectDataLabel &TOTAL) const = 0;

    ///////////////
    // operators //
    ///////////////

      //! @brief assemble all descriptor combinations based on an initial and a total descriptor set as an object label
      //! @param INITIAL initial descriptor set as an object label
      //! @param TOTAL all available descriptor groups in descriptor selection process
      //! @return container with all possible descriptor combinations based on initial descriptor set
      virtual const storage::Vector< util::ObjectDataLabel> operator()
      (
        const util::ObjectDataLabel &INITIAL,
        const util::ObjectDataLabel &TOTAL
      ) const = 0;

    }; // class DescriptorSelectionInterface

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_DESCRIPTOR_SELECTION_INTERFACE_H_
