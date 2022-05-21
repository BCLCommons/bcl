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
#include "descriptor/bcl_descriptor_example_string_sequence.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_base.h"
#include "descriptor/bcl_descriptor_iterator.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // enumerate this class for use in examples
    util::SiPtr< const util::ObjectInterface> StringSequence::s_Instance
    (
      util::Enumerated< StringSequence>::AddInstance( new StringSequence())
    );

    //! @brief constructor from a string
    StringSequence::StringSequence( const std::string &STRING) :
      m_String( STRING)
    {
    }

    //! @brief virtual copy constructor
    StringSequence *StringSequence::Clone() const
    {
      return new StringSequence( *this);
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &StringSequence::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the iterator for the sequence
    //! @return the iterator for the sequence
    iterate::Generic< const char> StringSequence::GetIterator() const
    {
      return iterate::Generic< const char>( m_String.begin(), m_String.end());
    }

    //! @brief get the alias for this class over the command line
    const std::string &StringSequence::GetAlias() const
    {
      static const std::string s_name( "String");
      return s_name;
    }

    //! @brief get the internally-held string
    const std::string &StringSequence::GetString() const
    {
      return m_String;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief get the iteration result as a string (semi-colons between entries)
    //! @param DESCRIPTOR the descriptor to use
    //! @param STRING the string of interest
    //! @param PRECISION the fixed float point precision # of digits to use
    //! @return the resulting string
    std::string StringSequence::WriteIterations
    (
      const util::Implementation< Base< char, float> > &DESCRIPTOR,
      const std::string &STRING,
      const size_t &PRECISION
    )
    {
      Iterator< char> itr( DESCRIPTOR->GetEffectiveType());
      StringSequence string_seq( STRING);
      util::Implementation< Base< char, float> > descriptor_copy( DESCRIPTOR);
      itr.SetObject( string_seq);
      descriptor_copy->SetObject( string_seq);
      std::ostringstream strstrm;
      util::Format formatter;
      formatter.FFP( PRECISION);
      for( ; itr.NotAtEnd(); ++itr)
      {
        const linal::VectorConstReference< float> result( ( *descriptor_copy)( itr));
        for
        (
          linal::VectorConstReference< float>::const_iterator
            itr_flt( result.Begin()), itr_flt_end( result.End());
          itr_flt != itr_flt_end;
          ++itr_flt
        )
        {
          strstrm << formatter( *itr_flt) << ' ';
        }
        strstrm << "; ";
      }

      return strstrm.str();
    }

    //! @brief get the iteration result as a string (semi-colons between entries)
    //! @param DESCRIPTOR the descriptor to use
    //! @param STRING the string of interest
    //! @return the resulting string
    std::string StringSequence::WriteIterations
    (
      const util::Implementation< Base< char, char> > &DESCRIPTOR,
      const std::string &STRING
    )
    {
      Iterator< char> itr( DESCRIPTOR->GetEffectiveType());
      StringSequence string_seq( STRING);
      util::Implementation< Base< char, char> > descriptor_copy( DESCRIPTOR);
      itr.SetObject( string_seq);
      descriptor_copy->SetObject( string_seq);
      std::ostringstream strstrm;
      for( ; itr.NotAtEnd(); ++itr)
      {
        const linal::VectorConstReference< char> result( ( *descriptor_copy)( itr));
        for
        (
          linal::VectorConstReference< char>::const_iterator
            itr_char( result.Begin()), itr_char_end( result.End());
          itr_char != itr_char_end;
          ++itr_char
        )
        {
          if( isprint( *itr_char))
          {
            strstrm << *itr_char;
          }
          else
          {
            strstrm << ' ';
          }
        }
        strstrm << "; ";
      }

      return strstrm.str();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer StringSequence::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Stores a string.");
      serializer.AddInitializer
      (
         "",
         "string to iterate over",
         io::Serialization::GetAgent( &m_String)
      );

      return serializer;
    }

  } // namespace descriptor
} // namespace bcl
