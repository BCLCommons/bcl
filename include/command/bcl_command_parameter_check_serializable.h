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

#ifndef BCL_COMMAND_PARAMETER_CHECK_SERIALIZABLE_H_
#define BCL_COMMAND_PARAMETER_CHECK_SERIALIZABLE_H_

// include the namespace header
#include "bcl_command.h"

// include other forward headers - sorted alphabetically
#include "io/bcl_io.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_command_parameter_check_interface.h"
#include "util/bcl_util_own_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ParameterCheckSerializable
    //! @brief This class is a command line helper class derived from ParameterCheckInterface that checks, whether a given\n
    //! parameter is within a given list.
    //!
    //! @see @link example_command_parameter_check_serializable.cpp @endlink
    //! @author mendenjl
    //! @date Jan 19, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ParameterCheckSerializable :
      public ParameterCheckInterface
    {

    private:

    //////////
    // data //
    //////////

      //! Serializable object that can be used to check whether something is a valid value
      util::OwnPtr< io::SerializationInterface> m_SerializableInstance;

    public:

      static const util::SiPtr< const util::ObjectInterface> s_Instance; //!< single instance of that class

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ParameterCheckSerializable();

      //! @brief construct from serializable object
      //! @param OBJECT vector of allowed strings
      ParameterCheckSerializable( const io::SerializationInterface &OBJECT);

      //! @brief virtual copy constructor
      //! @return pointer to new ParameterCheckSerializable object
      ParameterCheckSerializable *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return const ref std::string - the class name
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief returns if PARAMETER is allowed, i.e. is in the list of allowed strings
      //! @param PARAMETER the parameter to check
      //! @param PARAMETER_NAME the name of the parameter being checked
      //! @param ERROR_STREAM the stream to which errors are written
      //! @return returns true if the parameter is allowed, false otherwise
      bool IsAllowedParameter
      (
        const std::string &PARAMETER,
        const std::string &PARAMETER_NAME,
        std::ostream &ERROR_STREAM
      ) const;

      //! @brief check the parameter, returning CheckResult
      //! @param PARAMETER the parameter to check
      //! @param PARAMETER_NAME the name of the parameter being checked
      //! @param ERROR_STREAM the stream to which errors are written
      //! @return returns CheckResult
      io::ValidationResult Check
      (
        const std::string &PARAMETER,
        const std::string &PARAMETER_NAME,
        std::ostream &ERROR_STREAM
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief writes the help for the commandline
      //! @param OSTREAM - the stream to which you should write
      //! @param INDENT - amount to indent each new line after the first
      //! @return std::ostream &OSTREAM - return the stream after you write to it
      //! @see WriteList
      std::ostream &WriteHelp( std::ostream &OSTREAM, const size_t INDENT = 0) const;

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      //! @see ParameterCheckInterface::Write
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // ParameterAllowed

  } // namespace command
} // namespace bcl

#endif //BCL_COMMAND_PARAMETER_CHECK_SERIALIZABLE_H_
