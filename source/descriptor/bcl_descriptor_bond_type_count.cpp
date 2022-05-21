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
#include "descriptor/bcl_descriptor_bond_type_count.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> BondTypeCount::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new BondTypeCount()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    BondTypeCount::BondTypeCount() :
      m_BondProperty(),
      m_DesiredValue( 1)
    {
    }

    //! @brief constructor from number of steps, and mapped atom property
    BondTypeCount::BondTypeCount
    (
      const chemistry::ConfigurationalBondTypeData::Data &PROPERTY,
      const size_t &PARAMETER
    ) :
      m_BondProperty( PROPERTY),
      m_DesiredValue( PARAMETER)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to new BondTypeCount
    BondTypeCount *BondTypeCount::Clone() const
    {
      return new BondTypeCount( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &BondTypeCount::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &BondTypeCount::GetAlias() const
    {
      static const std::string s_name( "BondTypeCount");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void BondTypeCount::Calculate
    (
      const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      STORAGE = ELEMENT->CountNonValenceBondsWithProperty( m_BondProperty, m_DesiredValue);
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer BondTypeCount::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Counts the number of bonds that satisfy a condition (property == value)"
      );

      parameters.AddInitializer
      (
        "property",
        "bond property to query",
        io::Serialization::GetAgent( &m_BondProperty),
        "Identity"
      );
      parameters.AddInitializer
      (
        "value",
        "value to calculate number of bonds satisfying property == value",
        io::Serialization::GetAgent( &m_DesiredValue),
        "1"
      );
      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
