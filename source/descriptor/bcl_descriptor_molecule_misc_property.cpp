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
#include "descriptor/bcl_descriptor_molecule_misc_property.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_interface.h"
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> MoleculeMiscProperty::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new MoleculeMiscProperty()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MoleculeMiscProperty::MoleculeMiscProperty() :
      m_MiscPropertyString(),
      m_FeaturesPerMolecule( 0)
    {
    }

    //! @brief constructor from a property data label, which must include a miscellaneous property
    //! @param NAME the property name to use
    //! @param PROPERTIES_PER_MOLECULE # of properties to expect to be given per small molecule
    MoleculeMiscProperty::MoleculeMiscProperty
    (
      const std::string &NAME,
      const size_t &PROPERTIES_PER_MOLECULE
    ) :
      m_MiscPropertyString( NAME),
      m_FeaturesPerMolecule( PROPERTIES_PER_MOLECULE)
    {
      BCL_Assert
      (
        PROPERTIES_PER_MOLECULE > 0 && util::IsDefined( PROPERTIES_PER_MOLECULE),
        "invalid # of properties per misc small molecule property named "
        + NAME + ", # properties was " + util::Format()( PROPERTIES_PER_MOLECULE)
      );
    }

    //! @brief virtual copy constructor
    //! @return pointer to new MoleculeMiscProperty
    MoleculeMiscProperty *MoleculeMiscProperty::Clone() const
    {
      return new MoleculeMiscProperty( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &MoleculeMiscProperty::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &MoleculeMiscProperty::GetAlias() const
    {
      static const std::string s_name( "MiscProperty");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void MoleculeMiscProperty::Calculate( linal::VectorReference< float> &STORAGE)
    {
      util::SiPtr< const chemistry::ConformationInterface> si_molecule( this->GetCurrentObject());
      const chemistry::ConformationInterface &molecule( *si_molecule);

      // check for the property in the storage
      if( molecule.IsPropertyStored( m_MiscPropertyString))
      {
        STORAGE.CopyValues( molecule.GetMDLPropertyAsVector( m_MiscPropertyString));
      }
      // check for the property in the cache
      else if( molecule.IsCached( m_MiscPropertyString))
      {
        STORAGE.CopyValues( molecule.GetFromCache( m_MiscPropertyString));
      }
      else
      {
        STORAGE = util::GetUndefined< float>();
      }

    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoleculeMiscProperty::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "retrieves a molecule property that cannot be calculated directly by the bcl, e.g. biological data"
      );

      parameters.AddInitializer
      (
        "",
        "name of the miscellaneous property",
        io::Serialization::GetAgent( &m_MiscPropertyString)
      );
      parameters.AddInitializer
      (
        "values per molecule",
        "expected number of values per molecule",
        io::Serialization::GetAgentWithRange( &m_FeaturesPerMolecule, 1, 100000),
        "1"
      );

      return parameters;
    }
  } // namespace descriptor
} // namespace bcl
