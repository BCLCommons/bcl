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
#include "descriptor/bcl_descriptor_molecule_storage_id.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_molecule_storage_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> MoleculeStorageId::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, char> >::AddInstance
      (
        new MoleculeStorageId( new chemistry::MoleculeStorageFile( true))
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MoleculeStorageId::MoleculeStorageId()
    {
    }

    //! @brief constructor from storage created with new
    MoleculeStorageId::MoleculeStorageId
    (
      io::StoreInterface< chemistry::ConformationInterface> *const STORAGE
    ) :
      m_Alias( STORAGE->GetAlias() + "ID"),
      m_ConformationStorage( STORAGE)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to new MoleculeStorageId
    MoleculeStorageId *MoleculeStorageId::Clone() const
    {
      return new MoleculeStorageId( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &MoleculeStorageId::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &MoleculeStorageId::GetAlias() const
    {
      return m_Alias;
    }

    //! @brief get the feature siz  e under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t MoleculeStorageId::GetNormalSizeOfFeatures() const
    {
      return 14;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void MoleculeStorageId::Calculate( linal::VectorReference< char> &STORAGE)
    {
      util::SiPtr< const chemistry::ConformationInterface> sp_mol( GetCurrentObject());
      const chemistry::ConformationInterface &molecule( *sp_mol);

      // get the id of this molecule
      const std::string id
      (
        util::Format().W( GetNormalSizeOfFeatures() - 1)
        (
          m_ConformationStorage.IsDefined()
          ? m_ConformationStorage->Store( molecule)
          : std::string( "0")
        )
      );
      // always prepend with a space for readability
      STORAGE( 0) = ' ';
      std::copy( id.begin(), id.end(), STORAGE.Begin() + 1);
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoleculeStorageId::GetSerializer() const
    {
      io::Serializer serializer( m_ConformationStorage.GetSerializer());
      serializer.SetClassDescription( "ID of this molecule from a " + m_ConformationStorage.GetAlias());
      return serializer;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool MoleculeStorageId::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
    {
      return m_ConformationStorage->ReadInitializerSuccessHook( LABEL, ERR_STREAM);
    }

  } // namespace descriptor
} // namespace bcl
