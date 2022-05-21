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
#include "chemistry/bcl_chemistry_small_molecule_string_properties_numeric.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_iterator.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // single instance of this class
    const util::SiPtr< const util::ObjectInterface> SmallMoleculeStringPropertiesNumeric::s_Instance
    (
      util::Enumerated< StringPropertyInterface>::AddInstance
      (
        new SmallMoleculeStringPropertiesNumeric()
      )
    );

    // Instance of this class that serves as the default instance
    const util::SiPtr< const util::ObjectInterface> SmallMoleculeStringPropertiesNumeric::s_DefaultInstance
    (
      util::Enumerated< StringPropertyInterface>::AddInstance
      (
        new SmallMoleculeStringPropertiesNumeric( "")
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor from alias
    SmallMoleculeStringPropertiesNumeric::SmallMoleculeStringPropertiesNumeric( const std::string &ALIAS) :
      m_Alias( ALIAS)
    {
    }

    //! virtual copy constructor
    SmallMoleculeStringPropertiesNumeric *SmallMoleculeStringPropertiesNumeric::Clone() const
    {
      return new SmallMoleculeStringPropertiesNumeric( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SmallMoleculeStringPropertiesNumeric::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &SmallMoleculeStringPropertiesNumeric::GetAlias() const
    {
      return m_Alias;
    }

    //! @brief operator the implements the assignment operation on the two arguments returning a result
    //! @param MOLECULE the molecule to calculate the property for
    //! @return the property as a string
    std::string SmallMoleculeStringPropertiesNumeric::operator()( const ConformationInterface &MOLECULE) const
    {
      if( !m_Property.IsDefined())
      {
        return std::string();
      }

      m_Property->SetObject( MOLECULE);

      // convert vector into string; use io::Serialize so that nan's are output as "nan" independent of machine
      std::stringstream value;
      for
      (
        descriptor::Iterator< AtomConformationalInterface> itr_desc( m_Property->GetType(), MOLECULE);
        itr_desc.NotAtEnd();
        ++itr_desc
      )
      {
        linal::VectorConstReference< float> values( m_Property->operator()( itr_desc));
        for
        (
          linal::Vector< float>::const_iterator itr( values.Begin()), itr_end( values.End());
          itr != itr_end;
          ++itr
        )
        {
          io::Serialize::Write( *itr, value) << ' ';
        }
      }

      // return the string
      return util::TrimString( value.str());
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SmallMoleculeStringPropertiesNumeric::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Calculates a numeric molecular or atomic property from the molecule and converts it into a string"
      );

      parameters.AddOptionalInitializer
      (
        "",
        "data label for an atom or small molecule property",
        io::Serialization::GetAgent( &m_Property)
      );

      return parameters;
    }

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    //! @return result of any validation performed internally
    io::ValidationResult SmallMoleculeStringPropertiesNumeric::PreReadHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERROR_STREAM
    )
    {
      if( LABEL.GetValue() == GetAlias() && LABEL.GetNumberArguments() == 1)
      {
        io::ValidationResult could_read( m_Property.ValidateRead( LABEL.GetArgument( 0), ERROR_STREAM));
        if( !could_read)
        {
          return could_read;
        }
      }
      else
      {
        io::ValidationResult could_read( m_Property.ValidateRead( LABEL, ERROR_STREAM));
        if( !could_read)
        {
          return could_read;
        }
      }

      // this object is already completely serialized at this point
      return io::ValidationResult( io::e_Complete);
    }

  } // namespace chemistry
} // namespace bcl
