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

#ifndef BCL_COMMAND_PARAMETER_INTERFACE_H_
#define BCL_COMMAND_PARAMETER_INTERFACE_H_

// include the namespace header
#include "bcl_command.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ParameterInterface
    //! @brief is a command line helper class that defines a parameter in the commandline.
    //! it provides get and set functions as well as output helper
    //!
    //! @remarks example unnecessary
    //! @author woetzen, heinzes1, karakam
    //! @date 5/5/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ParameterInterface :
      public util::ObjectInterface
    {

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns Name
      //! @return std::string - name
      virtual const std::string &GetName() const = 0;

      //! @brief set the name
      //! @param NAME name for the parameter
      virtual void SetName( const std::string &NAME) = 0;

      //! @brief returns description
      //! @return std::string-description
      virtual const std::string &GetDescription() const = 0;

      //! @brief set the parameter from a given string VALUE
      //! @param VALUE the value of the parameter as given in the commandline as an argument
      //! @param ERROR_STREAM stream to write an error to
      //! @return true on success and parameter was compatible
      virtual bool SetParameter( const std::string &VALUE, std::ostream &ERROR_STREAM) = 0;

      //! @brief returns m_ParameterInterface
      //! @return value of the parameter
      virtual const std::string &GetValue() const = 0;

      //! @brief returns m_DefaultParameterInterface
      //! @return value of the default parameter
      virtual const std::string &GetDefaultValue() const = 0;

      //! @brief set the default parameter
      //! @param PARAMETER default value for the parameter
      virtual void SetDefaultParameter( const std::string &PARAMETER) = 0;

      //! @brief return parameter converted to a numerical value
      //! @return template< typename T1> - parameter converted to a numerical value
      template< typename T1>
      T1 GetNumericalValue() const
      {
        return util::ConvertStringToNumericalValue< T1>( GetValue());
      }

      //! @brief returns if parameter was set in the command line
      //! @return m_WasSetInCommandline
      virtual bool GetWasSetInCommandLine() const = 0;

      //! @brief return whether default value was given
      //! @return m_DefaultGiven
      virtual bool GetWasDefaultGiven() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief returns if PARAMETER is allowed, i.e. the object behind the ParameterCheckInterface allows it
      //! @param PARAMETER the parameter to check
      //! @param ERROR_STREAM the stream to which errors are written
      //! @return returns true if the parameter is allowed, false otherwise
      virtual bool IsAllowedParameter( const std::string &VALUE, std::ostream &ERROR_STREAM) const = 0;

      //! @brief reset into original state
      virtual void Reset() = 0;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief writes the help for the commandline
      //! @param OSTREAM the stream to which the help is written to
      //! @param INDENT the amount of indent to write before the help
      //! @return the given stream to which the help was written to
      virtual std::ostream &WriteHelp( std::ostream &OSTREAM, const size_t INDENT = 0) const = 0;

      //! @brief writes the user provided commandline
      //! @param OSTREAM the stream to which the commandline is written
      //! @return given ostream to which the commandline was written
      virtual std::ostream &WriteUserCommand( std::ostream &OSTREAM) const = 0;
      
    }; // class ParameterInterface
  } // namespace command
} // namespace bcl

#endif //BCL_COMMAND_PARAMETER_INTERFACE_H_
