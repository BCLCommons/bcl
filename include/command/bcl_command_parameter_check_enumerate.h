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

#ifndef BCL_COMMAND_PARAMETER_CHECK_ENUMERATE_H_
#define BCL_COMMAND_PARAMETER_CHECK_ENUMERATE_H_

// include the namespace header
#include "bcl_command.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_command_parameter_check_interface.h"
#include "util/bcl_util_object_data_label.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ParameterCheckEnumerate
    //! @brief ParameterCheckInterface implementation for util::Enumerate derived classes
    //! @details each Enumerate derived class has enums with different names, that can be passed over the command line
    //! but needs the be validated as a valid enum. "Undefined" will be accepted as valid parameter value, but will not
    //! show in the command line as a value. This is supposed to be used in places, where enumerators are dynamically
    //! extended by static initialization using the AddEnum function on the enumerator.
    //!
    //! @tparam t_Enumerate util::Enumerate derived class
    //!
    //! @see @link example_command_parameter_check_enumerate.cpp @endlink
    //! @author woetzen
    //! @date Sep 12, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Enumerate>
    class ParameterCheckEnumerate :
      public ParameterCheckInterface
    {

    public:

      static const util::SiPtr< const util::ObjectInterface> s_Instance; //!< single instance of that class

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new ParameterCheckEnumerate< t_Enumerate>
      ParameterCheckEnumerate< t_Enumerate> *Clone() const
      {
        return new ParameterCheckEnumerate< t_Enumerate>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief returns if PARAMETER is allowed, i.e. is in the enum list
      //! @param PARAMETER the parameter to check
      //! @param PARAMETER_NAME the name of the parameter being checked
      //! @param ERROR_STREAM the stream to which errors are written
      //! @return returns true if the parameter is allowed, false otherwise
      bool IsAllowedParameter
      (
        const std::string &PARAMETER,
        const std::string &PARAMETER_NAME,
        std::ostream &ERROR_STREAM
      ) const
      {
        typename t_Enumerate::EnumType try_enum;
        return try_enum.TryRead( util::ObjectDataLabel( "", PARAMETER), ERROR_STREAM);
      }

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief writes the help for the commandline
      //! @param OSTREAM output stream to be written to
      //! @param INDENT the amount to indent each new line after the first
      //! @return ostream which was written to
      std::ostream &WriteHelp( std::ostream &OSTREAM, const size_t INDENT = 0) const
      {
        typename t_Enumerate::EnumType try_enum;
        // end
        return try_enum.WriteHelp( OSTREAM, INDENT);
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // return the stream
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      //! @see ParameterCheckInterface::Write
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // return the stream
        return OSTREAM;
      }

    }; // template class ParameterCheckEnumerate

    // instantiate s_Instance
    template< typename t_Enumerate>
    const util::SiPtr< const util::ObjectInterface> ParameterCheckEnumerate< t_Enumerate>::s_Instance
    (
      GetObjectInstances().AddInstance( new ParameterCheckEnumerate< t_Enumerate>())
    );

  } // namespace command
} // namespace bcl

#endif // BCL_COMMAND_PARAMETER_CHECK_ENUMERATE_H_
