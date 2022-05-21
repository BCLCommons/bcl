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

#ifndef BCL_OPTI_PROCESSOR_INTERFACE_H_
#define BCL_OPTI_PROCESSOR_INTERFACE_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProcessorInterface
    //! @brief Interface for algorithms processing optimization targets.
    //! @details This interface should be implemented by classes used during the pre- and post-processing of the
    //! optimization target within the optimization framework. For example this interface can be used for
    //! printers to write results or for initial filtering of the argument.
    //!
    //! @remarks example unnecessary
    //! @author fischea
    //! @date Apr 17, 2017
    //!
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_ArgumentType>
    class ProcessorInterface :
      public util::SerializableInterface
    {

    ////////////////////////////////////
    //  construction and destruction  //
    ////////////////////////////////////

    public:

      //! @brief virtual copy constructor
      //! @return pointer to a new ProcessorInterface
      virtual ProcessorInterface *Clone() const = 0;

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      virtual const std::string &GetAlias() const = 0;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      virtual io::Serializer GetSerializer() const = 0;

    ///////////////
    // operators //
    ///////////////

      //! @brief optimizes a given argument
      //! @param ARGUMENT data to be optimized
      virtual void operator()( t_ArgumentType &ARGUMENT) const = 0;

    }; // class ProcessorInterface

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_PROCESSOR_INTERFACE_H_
