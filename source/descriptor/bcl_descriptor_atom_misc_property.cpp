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
#include "descriptor/bcl_descriptor_atom_misc_property.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_interface.h"
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> AtomMiscProperty::s_Instance
    (
      util::Enumerated< BaseElement< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new AtomMiscProperty()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AtomMiscProperty::AtomMiscProperty() :
      m_MiscPropertyString(),
      m_FeaturesPerAtom( 0)
    {
    }

    //! @brief constructor from number of steps, and mapped atom property
    AtomMiscProperty::AtomMiscProperty
    (
      const std::string &NAME,
      const size_t &PROPERTIES_PER_ATOM
    ) :
      m_MiscPropertyString( NAME),
      m_FeaturesPerAtom( PROPERTIES_PER_ATOM)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to new AtomMiscProperty
    AtomMiscProperty *AtomMiscProperty::Clone() const
    {
      return new AtomMiscProperty( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &AtomMiscProperty::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &AtomMiscProperty::GetAlias() const
    {
      static const std::string s_name( "Atom_MiscProperty");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void AtomMiscProperty::Calculate
    (
      const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      util::SiPtr< const chemistry::ConformationInterface> si_molecule( this->GetCurrentObject());
      const chemistry::ConformationInterface &molecule( *si_molecule);

      // check for the property in the cache
      if( molecule.IsCached( m_MiscPropertyString))
      {
        STORAGE = molecule.GetFromCache( m_MiscPropertyString)( ELEMENT.GetPosition());
      }

      // check for the property in the storage
      if( molecule.IsPropertyStored( m_MiscPropertyString))
      {
        STORAGE = molecule.GetMDLPropertyAsVector( m_MiscPropertyString)( ELEMENT.GetPosition());
      }
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AtomMiscProperty::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "retrieves an atom property that cannot be calculated directly by the bcl, e.g. spectra"
      );

      parameters.AddInitializer
      (
        "",
        "name of the miscellaneous property",
        io::Serialization::GetAgent( &m_MiscPropertyString)
      );
      parameters.AddInitializer
      (
        "values per atom",
        "expected number of values per atom",
        io::Serialization::GetAgentWithRange( &m_FeaturesPerAtom, 1, 100000),
        "1"
      );

      return parameters;
    }
  } // namespace descriptor
} // namespace bcl
