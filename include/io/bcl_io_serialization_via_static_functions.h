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

#ifndef BCL_IO_SERIALIZATION_VIA_STATIC_FUNCTIONS_H_
#define BCL_IO_SERIALIZATION_VIA_STATIC_FUNCTIONS_H_

// include the namespace header
#include "bcl_io.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_io_serialization_base.h"
#include "command/bcl_command_parameter_check_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace io
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SerializationViaStaticFunctions
    //! @brief handles labels for cases when data is retrieved through functions that are static or outside class scope
    //!
    //! @see @link example_io_serialization_via_static_functions.cpp @endlink
    //! @author mendenjl
    //! @date Nov 01, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_DataType>
    class SerializationViaStaticFunctions :
      public SerializationBase< t_DataType>
    {
    //////////
    // data //
    //////////

      std::string( *m_Getter)( const t_DataType &); //!< Function that gets a string representation of the member

       //! Function that gets a t_DataType using a particular string and error stream
      bool ( *m_Setter)( t_DataType &, const std::string &, std::ostream &);

      //! parameter check
      util::ShPtr< command::ParameterCheckInterface> m_Check;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SerializationViaStaticFunctions()
      {
      }

      //! @brief construct from base class
      //! @param PARAMETER parameter to copy
      //! @param GETTER function to get the string
      //! @param SETTER function to set the string
      //! @param DATA reference to the associated member variable
      SerializationViaStaticFunctions
      (
        std::string ( *GETTER)( const t_DataType &),
        bool ( *SETTER)( t_DataType &, const std::string &, std::ostream &),
        const t_DataType *DATA = NULL
      ) :
        SerializationBase< t_DataType>( DATA),
        m_Getter( GETTER),
        m_Setter( SETTER)
      {
      }

      //! @brief construct from base class
      //! @param PARAMETER_CHECK any util::ParameterCheckInterface derived class
      //! @param GETTER function that gets a string representation of the value from the class
      //! @param SETTER function that sets the value of the class using a string and error stream
      //! @param DATA reference to the associated member variable
      SerializationViaStaticFunctions
      (
        const command::ParameterCheckInterface &PARAMETER_CHECK,
        std::string ( *GETTER)( const t_DataType &),
        bool ( *SETTER)( t_DataType &, const std::string &, std::ostream &),
        const t_DataType *DATA = NULL
      ) :
        SerializationBase< t_DataType>( DATA),
        m_Getter( GETTER),
        m_Setter( SETTER),
        m_Check( PARAMETER_CHECK.Clone())
      {
      }

      //! @brief Clone function
      //! @return pointer to new SerializationViaStaticFunctions
      SerializationViaStaticFunctions *Clone() const
      {
        return new SerializationViaStaticFunctions( *this);
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

      //! @brief Get the label for an object of the given type
      //! @param OBJECT the object to get the label for
      //! @param WITH_DATA whether to include any data members, else, only include initialization members
      util::ObjectDataLabel GetLabelForObject( const t_DataType &OBJECT, const bool &WITH_DATA) const
      {
        return util::ObjectDataLabel( "", ( *m_Getter)( OBJECT));
      }

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief set the given object using the label
      //! @param OBJECT object to read
      //! @param LABEL label that is used to set the string
      //! @param ERR_STREAM stream to write errors to
      //! @return bool - if the parameter was set (parameters are not set if they are not allowed parameters)
      bool TryReadObject( t_DataType &OBJECT, const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM) const
      {
        // get the original length of the error stream, usually 0
        if( !m_Check.IsDefined() || m_Check->IsAllowedParameter( LABEL.GetValue(), "", ERR_STREAM))
        {
          return ( *m_Setter)( OBJECT, LABEL.GetValue(), ERR_STREAM);
        }

        return false;
      }

      //! @brief writes the help for the label
      //! @param OSTREAM the stream to which the help is written to
      //! @param INDENT the amount of indent
      //! @return the given stream to which the help was written to
      std::ostream &WriteHelp( std::ostream &OSTREAM, const size_t INDENT = 0) const
      {
        if( m_Check.IsDefined())
        {
          m_Check->WriteHelp( OSTREAM, INDENT);
        }
        return OSTREAM;
      }

    }; // class SerializationViaStaticFunctions

  } // namespace io
} // namespace bcl

#endif // BCL_IO_SERIALIZATION_VIA_STATIC_FUNCTIONS_H_

