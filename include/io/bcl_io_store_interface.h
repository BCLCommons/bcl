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

#ifndef BCL_IO_STORE_INTERFACE_H_
#define BCL_IO_STORE_INTERFACE_H_

// include the namespace header
#include "bcl_io.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace io
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class StoreInterface
    //! @brief class defines an interface for storing objects of templated type
    //!
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date Jan 09, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_StoredType>
    class StoreInterface :
      public virtual util::SerializableInterface
    {

    public:

      //! @brief Clone the StoreInterface
      //! @return pointer to new StoreInterface
      virtual StoreInterface *Clone() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief store an item
      //! @param ITEM item to store
      //! @return key associated with the stored item
      virtual std::string Store( const t_StoredType &ITEM) = 0;

    }; // class StoreInterface

  } // namespace io
} // namespace bcl

#endif // BCL_IO_STORE_INTERFACE_H_
