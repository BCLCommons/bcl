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

#ifndef BCL_COMMAND_PARAMETER_H_
#define BCL_COMMAND_PARAMETER_H_

// include the namespace header
#include "bcl_command.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_command_parameter_check_interface.h"
#include "bcl_command_parameter_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Parameter
    //! @brief This class is a command line helper class and contains a parameter and you may pass a function\n
    //! derived from ParameterCheckInterface that will check the parameter, if you want to set it from\n
    //! command line.
    //!
    //! @see @link example_command_parameter.cpp @endlink
    //! @author heinzes1, woetzen, karakam
    //! @date 10/12/2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Parameter :
      public ParameterInterface
    {

    private:

    //////////
    // data //
    //////////

      std::string m_Name; //!< name of parameter

      std::string m_Description; //!< describes the parameter

      util::ShPtr< ParameterCheckInterface> m_ParameterCheck; //!< function that checks if parameter is valid

      std::string m_Parameter; //!< this is the actual parameter

      bool m_WasSetInCommandline; //!< this flag will be true if the parameter was set in the commandline

      bool m_DefaultGiven; //!< this flag is set to true if a default is given

      std::string m_DefaultParameter; //!< stores the default

    public:

      static const util::SiPtr< const util::ObjectInterface> s_Instance; //!< single instance of that class

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Parameter();

      //! @brief construct from description string, parameter check function
      //! @param NAME name of parameter
      //! @param DESCRIPTION description string of parameter
      Parameter( const std::string &NAME, const std::string &DESCRIPTION);

      //! @brief construct from description string, parameter check function
      //! @param NAME name of parameter
      //! @param DESCRIPTION description string of parameter
      //! @param PARAMETER_CHECK ParameterCheckInterface derived class to be used to check parameter
      Parameter
      (
        const std::string &NAME,
        const std::string &DESCRIPTION,
        const ParameterCheckInterface &PARAMETER_CHECK
      );

      //! @brief construct from name, description string and default parameter string
      //! @param NAME name of parameter
      //! @param DESCRIPTION description string of parameter
      //! @param DEFAULT_PARAMETER Default string value of the parameter
      Parameter
      (
        const std::string &NAME,
        const std::string &DESCRIPTION,
        const std::string &DEFAULT_PARAMETER
      );

      //! @brief construct from description string, parameter check, and default parameter string
      //! @param NAME name of parameter
      //! @param DESCRIPTION description string of parameter
      //! @param PARAMETER_CHECK ParameterCheckInterface derived class to be used to check parameter
      //! @param DEFAULT_PARAMETER Default string value of the parameter
      Parameter
      (
        const std::string &NAME,
        const std::string &DESCRIPTION,
        const ParameterCheckInterface &PARAMETER_CHECK,
        const std::string &DEFAULT_PARAMETER
      );

      //! @brief copy constructor
      //! @return pointer to new Parameter object
      Parameter *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return std::string - the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns Name
      //! @return std::string - name
      const std::string &GetName() const;

      //! @brief set the name
      //! @param NAME name for the parameter
      void SetName( const std::string &NAME);

      //! @brief returns description
      //! @return std::string-description
      const std::string &GetDescription() const;

      //! @brief checks if PARAMETER is allowed and returns
      //! @param PARAMETER - the parameter we want to set
      //! @param ERROR_STREAM - the stream to which errors should be written
      //! @return bool - if the parameter was set (parameters are not set if they are not allowed parameters)
      bool SetParameter( const std::string &PARAMETER, std::ostream &ERROR_STREAM);

      //! @brief returns m_Parameter
      //! @return value of the parameter
      const std::string &GetValue() const;

      //! @brief returns m_DefaultParameter
      //! @return value of the default parameter
      const std::string &GetDefaultValue() const;

      //! @brief set the default parameter
      //! @param PARAMETER default value for the parameter
      void SetDefaultParameter( const std::string &PARAMETER);

      //! @brief return parameter converted to a numerical value
      //! @return template< typename T1> - parameter converted to a numerical value
      template< typename T1>
      T1 GetNumericalValue() const
      {
        return util::ConvertStringToNumericalValue< T1>( m_Parameter);
      }

      //! @brief returns if parameter was set in the command line
      //! @return m_WasSetInCommandline
      bool GetWasSetInCommandLine() const;

      //! @brief return whether default value was given
      //! @return m_DefaultGiven
      bool GetWasDefaultGiven() const;

      //! @brief access to the parameter check
      //! @return the ParameterCheck
      const ParameterCheckInterface &GetParameterCheck() const;

      //! @brief set the parameter check
      //! @param CHECK the new parameter check
      void SetParameterCheck( const ParameterCheckInterface &CHECK);

    ////////////////
    // operations //
    ////////////////

      //! @brief returns if PARAMETER is allowed, i.e. the object behind the ParameterCheckInterface allows it
      //! @param PARAMETER the parameter to check
      //! @param ERROR_STREAM the stream to which errors are written
      //! @return returns true if the parameter is allowed, false otherwise
      bool IsAllowedParameter( const std::string &PARAMETER, std::ostream &ERROR_STREAM) const;

      //! @brief reset into original state
      void Reset();

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief writes the help for the commandline
      //! @param OSTREAM the stream to which the help is written to
      //! @param INDENT number of indentations
      //! @return the given stream to which the help was written to
      std::ostream &WriteHelp( std::ostream &OSTREAM, const size_t INDENT = 0) const;

      //! @brief writes the user provided commandline
      //! @param OSTREAM the stream to which the commandline is written
      //! @return given ostream to which the commandline was written
      std::ostream &WriteUserCommand( std::ostream &OSTREAM) const;

    protected:

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

    }; // class Parameter

    //! @brief StringVectorFromFilenameParameter gives all the strings in a file as a Vector
    //! @param PARAMETER the name of the file
    //! @return return a Vector which has all of the strings contained within the file denoted by "FLAG"
    BCL_API storage::Vector< std::string> StringVectorFromFilenameParameter( const ParameterInterface &PARAMETER);

  } // namespace command
} // namespace bcl

#endif //BCL_COMMAND_PARAMETER_H_
