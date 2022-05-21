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
#include "descriptor/bcl_descriptor_sdf_file_id.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    // separate, anonymous namespace to prevent exporting symbols
    namespace
    {
      // add each of the possible instances to the enumerated instances
      util::ObjectInterface *AddInstances()
      {
        // keep a pointer to the last created instance
        util::ObjectInterface *last_instance( NULL);
        for( size_t level_type( 0); level_type < chemistry::MoleculeStorageFile::s_NumberOfFindTypes; ++level_type)
        {
          last_instance =
            util::Enumerated< Base< chemistry::AtomConformationalInterface, char> >::AddInstance
            (
              new SdfFileId( static_cast< chemistry::MoleculeStorageFile::FindType>( level_type))
            );
        }
        return last_instance;
      }
    }
    const util::SiPtr< const util::ObjectInterface> SdfFileId::s_Instance( AddInstances());

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from storage created with new
    SdfFileId::SdfFileId
    (
      chemistry::MoleculeStorageFile::FindTypeEnum TYPE
    ) :
      m_Alias( TYPE.GetString() + "IndexInSdfFile"),
      m_Level( TYPE),
      m_MoleculeStorage()
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to new SdfFileId
    SdfFileId *SdfFileId::Clone() const
    {
      return new SdfFileId( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &SdfFileId::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &SdfFileId::GetAlias() const
    {
      return m_Alias;
    }

    //! @brief get the feature siz  e under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t SdfFileId::GetNormalSizeOfFeatures() const
    {
      return 14;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void SdfFileId::Calculate( linal::VectorReference< char> &STORAGE)
    {
      util::SiPtr< const chemistry::ConformationInterface> sp_mol( GetCurrentObject());
      const chemistry::ConformationInterface &molecule( *sp_mol);

      std::pair< size_t, bool> find_result( 0, false);
      for
      (
        storage::Vector< chemistry::MoleculeStorageFile>::const_iterator
          itr( m_MoleculeStorage.Begin()), itr_end( m_MoleculeStorage.End());
        itr != itr_end;
        ++itr
      )
      {
        std::pair< size_t, bool> find_results_inner( itr->Find( molecule, m_Level));
        find_result.first += find_results_inner.first;
        find_result.second = find_results_inner.second;
        if( find_result.second)
        {
          break;
        }
      }
      if( !find_result.second)
      {
        BCL_MessageStd( "Molecule with name " + molecule.GetName() + " not found, returning number of molecules in file");
      }

      // get the id of this molecule
      const std::string id( util::Format().W( GetNormalSizeOfFeatures() - 1)( find_result.first));
      // always prepend with a space for readability
      STORAGE( 0) = ' ';
      std::copy( id.begin(), id.end(), STORAGE.Begin() + 1);
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SdfFileId::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Index of this molecule in an sdf file, searching at the " + m_Level.GetString() + " level");
      serializer.AddInitializer
      (
        "",
        "Filenames for all sdf files to consider; returned index indicates the index within the files as if they were "
        "concatenated",
        io::Serialization::GetAgent( &m_MoleculeStorage)
      );

      return serializer;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool SdfFileId::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
    {
      for
      (
        storage::Vector< chemistry::MoleculeStorageFile>::iterator
          itr( m_MoleculeStorage.Begin()), itr_end( m_MoleculeStorage.End());
        itr != itr_end;
        ++itr
      )
      {
        itr->ReadInitializerSuccessHook( LABEL, ERR_STREAM);
      }
      return true;
    }

  } // namespace descriptor
} // namespace bcl
