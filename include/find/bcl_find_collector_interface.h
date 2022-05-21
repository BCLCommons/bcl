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

#ifndef BCL_FIND_COLLECTOR_INTERFACE_H_
#define BCL_FIND_COLLECTOR_INTERFACE_H_

// include the namespace header
#include "bcl_find.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_class_descriptor.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace find
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CollectorInterface
    //! @brief Class is used for getting a SiPtrList of substituents (e.x. all CA atoms) from an
    //! argument (e.x. protein model).
    //!
    //! @tparam t_ReturnType is the type of substituents which will be collected
    //! @tparam t_ArgumenType is the type of object from which the substituents will be collected
    //!
    //! @remarks example unnecessary
    //! @author alexanns
    //! @date 03/276/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ReturnType, typename t_ArgumentType>
    class CollectorInterface :
      virtual public util::SerializableInterface
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! @brief Clone function.
      //! @return pointer to a new CollectorInterface
      virtual CollectorInterface *Clone() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! Collect the t_ReturnType objects in t_ArgumentType
      //! @param ARGUMENT entity that contains a t_ReturnType
      //! @return returns Group of the collected t_ReturnType objects
      virtual t_ReturnType Collect( const t_ArgumentType &ARGUMENT) const = 0;

    }; // template class CollectorInterface

  } // namespace find
} // namespace bcl

#endif //BCL_FIND_COLLECTOR_INTERFACE_H_
