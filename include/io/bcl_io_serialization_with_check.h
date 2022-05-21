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

#ifndef BCL_IO_SERIALIZATION_WITH_CHECK_H_
#define BCL_IO_SERIALIZATION_WITH_CHECK_H_

// include the namespace header
#include "bcl_io.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_io_serialization_builtin.h"
#include "command/bcl_command_parameter_check_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace io
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SerializationWithCheck
    //! @brief data labels for string
    //!
    //! @see @link example_io_serialization_with_check.cpp @endlink
    //! @author mendenjl
    //! @date Nov 01, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_DataType>
    class SerializationWithCheck :
      public SerializationBuiltin< t_DataType>
    {
    //////////
    // data //
    //////////

      //! parameter check
      util::ShPtr< command::ParameterCheckInterface> m_Check;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct from base class
      //! @param PARAMETER_CHECK any util::ParameterCheckInterface derived class
      //! @param DATA reference to the associated member variable
      SerializationWithCheck
      (
        const command::ParameterCheckInterface &PARAMETER_CHECK,
        const t_DataType *DATA
      );

      //! @brief Clone function
      //! @return pointer to new SerializationWithCheck
      SerializationWithCheck *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief set the given object using the label
      //! @param OBJECT object to read
      //! @param LABEL label that is used to set the string
      //! @param ERROR_STREAM stream to write errors to
      //! @return bool - if the parameter was set (parameters are not set if they are not allowed parameters)
      bool TryReadObject( t_DataType &OBJECT, const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM) const;

      //! @brief writes the help for the label
      //! @param OSTREAM the stream to which the help is written to
      //! @param INDENT the amount of indent
      //! @return the given stream to which the help was written to
      std::ostream &WriteHelp( std::ostream &OSTREAM, const size_t INDENT = 0) const;

    }; // class SerializationWithCheck

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API SerializationWithCheck< std::string>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SerializationWithCheck< util::ObjectDataLabel>;

  } // namespace io
} // namespace bcl

#endif // BCL_IO_SERIALIZATION_WITH_CHECK_H_

