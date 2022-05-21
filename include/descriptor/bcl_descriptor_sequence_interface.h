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

#ifndef BCL_DESCRIPTOR_SEQUENCE_INTERFACE_H_
#define BCL_DESCRIPTOR_SEQUENCE_INTERFACE_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_has_cache.h"
#include "iterate/bcl_iterate_generic.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SequenceInterface
    //! @brief objects that will be described by descriptors in this namespace need to inherit from this interface
    //!
    //! @tparam t_DataType type of the underlying sequence
    //!
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date Dec 04, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class SequenceInterface :
      public HasCache< util::ObjectInterface>
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new SequenceInterface
      virtual SequenceInterface< t_DataType> *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief return the length of the sequence in question
      //! @return the length of the sequence in question
      virtual size_t GetSize() const
      {
        return GetIterator().GetSize();
      }

      //! @brief get the iterator for the sequence
      //! @return the iterator for the sequence
      virtual iterate::Generic< const t_DataType> GetIterator() const = 0;

    }; // template class SequenceInterface

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_SEQUENCE_INTERFACE_H_
