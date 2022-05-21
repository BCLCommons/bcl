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

#ifndef BCL_FIND_LOCATOR_CRITERIA_INTERFACE_H_
#define BCL_FIND_LOCATOR_CRITERIA_INTERFACE_H_

// include the namespace header
#include "bcl_find.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace find
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LocatorCriteriaInterface
    //! @brief Interface to all locator criteria classes, that are used to identify and
    //! return a t_ReturnType from an object of t_ArgumentType
    //!
    //! @remarks example unnecessary
    //! @author alexanns
    //! @date 01/16/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ReturnType, typename t_ArgumentType, typename t_CriteriaType>
    class LocatorCriteriaInterface :
      public virtual util::SerializableInterface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief clone function
      //! @return pointer to a new LocatorCriteriaInterface
      virtual LocatorCriteriaInterface *Clone() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief locate the t_ReturnType in t_ArgumentType
      //! @param ARGUMENT entity that contains a t_ReturnType
      //! @param CRITERIA type of criteria
      //! @return returns the located t_ReturnType
      virtual t_ReturnType Locate( const t_ArgumentType &ARGUMENT, const t_CriteriaType &CRITERIA) const = 0;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      virtual io::Serializer GetSerializer() const
      {
        io::Serializer serializer;
        serializer.SetClassDescription( "This class needs a real implementation of GetSerializer().");

        return serializer;
      }

    }; // template class LocatorCriteriaInterface

  } // namespace find
} // namespace bcl

#endif //BCL_FIND_LOCATOR_CRITERIA_INTERFACE_H_
