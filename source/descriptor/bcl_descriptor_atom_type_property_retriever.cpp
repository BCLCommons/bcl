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
#include "descriptor/bcl_descriptor_atom_type_property_retriever.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AtomTypePropertyRetriever::AtomTypePropertyRetriever() :
      m_Alias( "UnknownAtomProperty"),
      m_Property( chemistry::AtomTypeData::s_NumberOfProperties)
    {
    }

    //! @brief constructor from a property
    //! @param PROPERTY the property to retrieve
    AtomTypePropertyRetriever::AtomTypePropertyRetriever
    (
      const chemistry::AtomTypeData::Properties &PROPERTY
    ) :
      m_Alias( std::string( "Atom_") + chemistry::AtomTypeData::PropertyEnum( PROPERTY).GetString()),
      m_Property( PROPERTY)
    {
    }

    //! @brief virtual copy constructor
    AtomTypePropertyRetriever *AtomTypePropertyRetriever::Clone() const
    {
      return new AtomTypePropertyRetriever( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &AtomTypePropertyRetriever::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &AtomTypePropertyRetriever::GetAlias() const
    {
      return m_Alias;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void AtomTypePropertyRetriever::Calculate
    (
      const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      STORAGE( 0) = ELEMENT->GetAtomType()->GetAtomTypeProperty( m_Property);
    } // Recalculate

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AtomTypePropertyRetriever::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Retrieves the " + chemistry::AtomTypeData::PropertyEnum( m_Property).GetString() + " of desired atom"
      );
      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
