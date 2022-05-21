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
#include "descriptor/bcl_descriptor_molecule_cached_string.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // single instance of that class
    const util::SiPtr< const util::ObjectInterface> MoleculeCachedString::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, char> >::AddInstance
      (
        new MoleculeCachedString()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MoleculeCachedString::MoleculeCachedString() :
      m_CachedPropertyName( "example"),
      m_NumberCharacters( 12)
    {
    }

    //! virtual copy constructor
    MoleculeCachedString *MoleculeCachedString::Clone() const
    {
      return new MoleculeCachedString( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MoleculeCachedString::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MoleculeCachedString::GetAlias() const
    {
      static const std::string s_name( "Cached");
      return s_name;
    }

    //! @brief calculate the descriptors for a given sequence iterator position, ignoring cache
    //! @param ITR sequence iterator pointing to the tuple of sequence elements of interest; will be for the native dimension of the class
    //! @param STORAGE storage for the descriptor
    //! The inheriting class can assume that ITR is of native dimension and that STORAGE is zeroed and of the proper
    //! dimension
    void MoleculeCachedString::RecalculateImpl
    (
      const Iterator< chemistry::AtomConformationalInterface> &ITR,
      linal::VectorReference< char> &STORAGE
    )
    {
      util::SiPtr< const chemistry::ConformationInterface> si_molecule( this->GetCurrentObject());
      const chemistry::ConformationInterface &molecule( *si_molecule);

      // look for the property directly on the molecule first; this avoids
      // problems with conversion size_t <> float

      util::SiPtr< const std::string> si_ptr_string_found;

      // check whether the property is scalar, in which case its value may have been stored
      std::string property_in_cache;
      if( m_CachedPropertyName.IsScalar())
      {
        if( molecule.IsPropertyStored( m_CachedPropertyName.GetValue()))
        {
          // return the cached property as a string
          property_in_cache = util::TrimString( molecule.GetMDLProperty( m_CachedPropertyName.GetValue()));
          si_ptr_string_found = util::SiPtr< const std::string>( property_in_cache);
        }
      }
      else
      {
        const std::string property_string( m_CachedPropertyName.ToString());
        if( molecule.IsPropertyStored( property_string))
        {
          // return the cached property as a string
          property_in_cache = util::TrimString( property_string);
          si_ptr_string_found = util::SiPtr< const std::string>( property_in_cache);
        }
      }

      if( !si_ptr_string_found.IsDefined())
      {
        // find the property in the cache
        util::SiPtr< const CacheMap::value_type> prop_ptr
        (
          molecule.FindInCache( m_CachedPropertyName)
        );

        // convert vector into string; use io::Serialize so that nan's are output as "nan" independent of machine
        if( prop_ptr->IsEmpty())
        {
          // check for empty vector
          return;
        }

        std::stringstream value;
        linal::Vector< float>::const_iterator itr( prop_ptr->Begin()), itr_end( prop_ptr->End());
        // write the first value
        io::Serialize::Write( *itr, value);
        for( ++itr; itr != itr_end; ++itr)
        {
          // write a space and then the next value
          io::Serialize::Write( *itr, value << ' ');
        }
        property_in_cache = util::TrimString( value.str());
        si_ptr_string_found = util::SiPtr< const std::string>( property_in_cache);
      }
      const std::string &s( *si_ptr_string_found);
      if( s.size() < m_NumberCharacters)
      {
        std::copy( s.begin(), s.end(), STORAGE.Begin());
        std::fill( STORAGE.Begin() + s.size(), STORAGE.End(), ' ');
      }
      else
      {
        std::copy( s.begin(), s.begin() + m_NumberCharacters, STORAGE.Begin());
      }
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoleculeCachedString::GetSerializer() const
    {
      io::Serializer init;
      init.SetClassDescription( "Retrieve a string already cached on the molecule");
      init.AddInitializer
      (
        "",
        "the name of a property stored on the molecule",
        io::Serialization::GetAgent( &m_CachedPropertyName),
        ""
      );
      init.AddInitializer
      (
        "size",
        "number of characters expected.  If the string is longer than this, it will be truncated to this length."
        "If it is shorter, it will be padded with spaces to this number of characters, provided that it exists on the "
        "molecule (otherwise it will be blank)",
        io::Serialization::GetAgent( &m_NumberCharacters),
        "12"
      );

      return init;
    }

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    bool MoleculeCachedString::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERROR_STREAM
    )
    {
      if( LABEL.GetValue() == GetAlias() && LABEL.GetNumberArguments() == size_t( 1))
      {
        // use the argument of the label
        m_CachedPropertyName = LABEL.GetArgument( 0);
        m_NumberCharacters = 12;
      }
      else if( LABEL.GetValue() == GetAlias() && LABEL.GetNumberArguments() == size_t( 2))
      {
        if( LABEL.GetArgument( 0).GetName() == "size")
        {
          m_CachedPropertyName = LABEL.GetArgument( 1);
          m_NumberCharacters = util::ConvertStringToNumericalValue< size_t>( LABEL.GetArgument( 0).GetValue());
        }
        else
        {
          m_CachedPropertyName = LABEL.GetArgument( 0);
          m_NumberCharacters = util::ConvertStringToNumericalValue< size_t>( LABEL.GetArgument( 1).GetValue());
        }
      }
      else
      {
        // look for the property with the name given in the label
        m_CachedPropertyName = LABEL;
      }
      return true;
    }
  } // namespace descriptor
} // namespace bcl
