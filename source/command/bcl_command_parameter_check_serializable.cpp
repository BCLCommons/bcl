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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "command/bcl_command_parameter_check_serializable.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_guesser.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ParameterCheckSerializable::s_Instance
    (
      GetObjectInstances().AddInstance( new ParameterCheckSerializable())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ParameterCheckSerializable::ParameterCheckSerializable()
    {
    }

    //! @brief construct from serializable object
    //! @param OBJECT vector of allowed strings
    ParameterCheckSerializable::ParameterCheckSerializable( const io::SerializationInterface &OBJECT) :
      m_SerializableInstance( OBJECT.Clone())
    {
    }

    //! @brief copy constructor
    //! @return pointer to new ParameterCheckSerializable object
    ParameterCheckSerializable *ParameterCheckSerializable::Clone() const
    {
      return new ParameterCheckSerializable( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return const ref std::string - the class name
    const std::string &ParameterCheckSerializable::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns if PARAMETER is allowed, i.e. is in the list of allowed strings
    //! @param PARAMETER the parameter to check
    //! @param PARAMETER_NAME the name of the parameter being checked
    //! @param ERROR_STREAM the stream to which errors are written
    //! @return returns true if the parameter is allowed, false otherwise
    bool ParameterCheckSerializable::IsAllowedParameter
    (
      const std::string &PARAMETER,
      const std::string &PARAMETER_NAME,
      std::ostream &ERROR_STREAM
    ) const
    {
      // try to read the parameter into an object data label
      util::ObjectDataLabel internal;
      if( !internal.TryAssign( PARAMETER, ERROR_STREAM))
      {
        return false;
      }

      // if no object was given, the parameter is allowed
      if( !m_SerializableInstance.IsDefined())
      {
        return true;
      }

      // copy the pointer / clone the object, to allow for read attempt
      util::OwnPtr< io::SerializationInterface> fresh_instance( m_SerializableInstance);

      // try to read the variable; return true on success
      return fresh_instance->TryRead( internal, ERROR_STREAM);
    }

    //! @brief check the parameter, returning CheckResult
    //! @param PARAMETER the parameter to check
    //! @param PARAMETER_NAME the name of the parameter being checked
    //! @param ERROR_STREAM the stream to which errors are written
    //! @return returns CheckResult
    io::ValidationResult ParameterCheckSerializable::Check
    (
      const std::string &PARAMETER,
      const std::string &PARAMETER_NAME,
      std::ostream &ERROR_STREAM
    ) const
    {
      // by default, check whether the parameter is equal to the global help string, if so, write the parameter's help
      // and return e_Help
      if( PARAMETER == io::ValidationResult::GetHelpString())
      {
        ERROR_STREAM << "Help for <" << PARAMETER_NAME << ">:\n";
        WriteHelp( ERROR_STREAM);
        return io::e_Help;
      }
      // create an object data label out of the parameter
      util::ObjectDataLabel label;
      if( !label.TryAssign( PARAMETER, ERROR_STREAM))
      {
        return io::e_Invalid;
      }

      // copy the pointer / clone the object, to allow for read attempt
      util::OwnPtr< io::SerializationInterface> fresh_instance( m_SerializableInstance);
      return fresh_instance->ValidateRead( label, ERROR_STREAM);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief writes the help for the commandline
    //! @param OSTREAM - the stream to which you should write
    //! @param INDENT - amount to indent each new line after the first
    //! @return std::ostream &OSTREAM - return the stream after you write to it
    //! @see WriteList
    std::ostream &ParameterCheckSerializable::WriteHelp( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // call write help o the internal object
      m_SerializableInstance->WriteHelp( OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    //! @see ParameterCheckInterface::Read
    std::istream &ParameterCheckSerializable::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SerializableInstance, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    //! @see ParameterCheckInterface::Write
    std::ostream &ParameterCheckSerializable::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SerializableInstance, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace command
} // namespace bcl
