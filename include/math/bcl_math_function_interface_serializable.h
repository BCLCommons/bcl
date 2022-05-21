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

#ifndef BCL_MATH_FUNCTION_INTERFACE_SERIALIZABLE_H_
#define BCL_MATH_FUNCTION_INTERFACE_SERIALIZABLE_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_function_interface.h"
#include "io/bcl_io_serializer.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FunctionInterfaceSerializable
    //! @brief just like FunctionInterface, except ultimately derived from SerializableInterface
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date Nov 01, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class FunctionInterfaceSerializable :
      public FunctionInterface< t_ArgumentType, t_ResultType>,
      virtual public util::SerializableInterface
    {
    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief virtual copy constructor
      virtual FunctionInterfaceSerializable *Clone() const = 0;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      virtual const std::string &GetAlias() const
      {
        return this->GetScheme();
      }

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      virtual io::Serializer GetSerializer() const
      {
        io::Serializer serializer;
        serializer.SetClassDescription( "This class needs a real implementation of GetSerializer().");

        return serializer;
      }

    }; // template class FunctionInterfaceSerializable

  } // namespace math
} // namespace bcl

#endif //BCL_MATH_FUNCTION_INTERFACE_SERIALIZABLE_H_
