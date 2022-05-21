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
#include "chemistry/bcl_chemistry_small_molecule_string_properties_cached.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // single instance of that class
    const util::SiPtr< const util::ObjectInterface> SmallMoleculeStringPropertiesCached::s_Instance
    (
      util::Enumerated< StringPropertyInterface>::AddInstance
      (
        new SmallMoleculeStringPropertiesCached()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SmallMoleculeStringPropertiesCached::SmallMoleculeStringPropertiesCached() :
      m_CachedPropertyName( "example")
    {
    }

    //! virtual copy constructor
    SmallMoleculeStringPropertiesCached *SmallMoleculeStringPropertiesCached::Clone() const
    {
      return new SmallMoleculeStringPropertiesCached( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SmallMoleculeStringPropertiesCached::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &SmallMoleculeStringPropertiesCached::GetAlias() const
    {
      static const std::string s_name( "Cached");
      return s_name;
    }

    //! @brief operator the implements the assignment operation on the two arguments returning a result
    //! @param MOLECULE the molecule to calculate the property for
    //! @return the property as a string
    std::string SmallMoleculeStringPropertiesCached::operator()( const ConformationInterface &MOLECULE) const
    {
      // look for the property directly on the molecule first; this avoids
      // problems with conversion size_t <> float 

      // check whether the property is scalar, in which case its value may have been stored
      if( m_CachedPropertyName.IsScalar())
      {
        // return the cached property as a string
        return util::TrimString( MOLECULE.GetMDLProperty( m_CachedPropertyName.GetValue()));
      }
      else
      {
        // return the cached property as a string
        return util::TrimString( MOLECULE.GetMDLProperty( m_CachedPropertyName.ToString()));
      }

      // find the property in the cache
      util::SiPtr< const descriptor::CacheMap::value_type> prop_ptr
      (
        MOLECULE.FindInCache( m_CachedPropertyName)
      );

      // convert vector into string; use io::Serialize so that nan's are output as "nan" independent of machine
      if( prop_ptr->IsEmpty())
      {
        // check for empty vector
        return std::string();
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

      // return the property as a string
      return value.str();
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SmallMoleculeStringPropertiesCached::GetSerializer() const
    {
      io::Serializer init;
      init.SetClassDescription( "Retrieve a property already cached on the molecule");
      init.AddInitializer
      (
        "",
        "the name of a property stored on the molecule",
        io::Serialization::GetAgent( &m_CachedPropertyName)
      );

      return init;
    }

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    bool SmallMoleculeStringPropertiesCached::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERROR_STREAM
    )
    {
      if( LABEL.GetValue() == GetAlias() && LABEL.GetNumberArguments() == size_t( 1))
      {
        // use the argument of the label
        m_CachedPropertyName = LABEL.GetArgument( 0);
      }
      else
      {
        // look for the property with the name given in the label
        m_CachedPropertyName = LABEL;
      }
      return true;
    }
  } // namespace chemistry
} // namespace bcl
