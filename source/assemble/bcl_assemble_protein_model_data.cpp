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
#include "assemble/bcl_assemble_protein_model_data.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  ///////////
  // enums //
  ///////////

    //! @brief conversion to a string from a Type
    //! @param TYPE the type to get a string for
    //! @return a string representing that type
    const std::string &ProteinModelData::GetTypeName( const Type &TYPE)
    {
      static const std::string s_descriptors[ s_NumberTypes + 1] =
      {
        "pool",
        "native_model",
        "native_filtered_model",
        "sasa",
        "multiplier",
        "loop_domain_locators",
        "pdb_file",
        "identification",
        "membrane",
        "slicelist_manager",
        "undefined",
        GetStaticClassName< ProteinModelData>()
      };
      return s_descriptors[ size_t( TYPE)];
    }

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ProteinModelData::s_Instance
    (
      GetObjectInstances().AddInstance( new ProteinModelData())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProteinModelData::ProteinModelData() :
      m_DataMap()
    {
      m_DataMap[ e_Undefined] = util::ShPtr< util::ObjectInterface>();
    }

    //! @brief Clone function
    //! @return pointer to new ProteinModelData
    ProteinModelData *ProteinModelData::Clone() const
    {
      return new ProteinModelData( *this);
    }

    //! @brief hardcopy all member shared pointers
    //! @return pointer to hardcopied data
    util::ShPtr< ProteinModelData> ProteinModelData::HardCopy() const
    {
      // initialize new pmd
      util::ShPtr< ProteinModelData> new_pmd( new ProteinModelData());

      // iterate over the map
      for
      (
        storage::Map< TypeEnum, util::ShPtr< util::ObjectInterface> >::const_iterator itr( m_DataMap.Begin()),
          itr_end( m_DataMap.End());
        itr != itr_end; ++itr
      )
      {
        // add the hardcopied data
        if( itr->second.IsDefined())
        {
          new_pmd->Insert( itr->first, itr->second.HardCopy());
        }
      }

      // end
      return new_pmd;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinModelData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief insert data to the map
    //! @param KEY the identifier for the data
    //! @param SP_DATA ShPtr to the ObjectInterface derived data
    //! @return true, if data with that key did not exist and was inserted, false otherwise
    bool ProteinModelData::Insert( const Type &KEY, const util::ShPtr< util::ObjectInterface> &SP_DATA)
    {
      // check if key exists
      if( m_DataMap.Find( KEY) != m_DataMap.End())
      {
        return false;
      }

      // set the data
      m_DataMap[ KEY] = SP_DATA;

      // success
      return true;
    }

    //! @brief insert data from another object into the map
    //! @param PROTEIN_MODEL_DATA object containing data to be added
    //! @param REPLACE if true, data of the same key is replaced
    void ProteinModelData::Insert( const ProteinModelData &PROTEIN_MODEL_DATA, const bool REPLACE)
    {
      // iterate over the map
      for
      (
        storage::Map< TypeEnum, util::ShPtr< util::ObjectInterface> >::const_iterator
          itr( PROTEIN_MODEL_DATA.m_DataMap.Begin()),
          itr_end( PROTEIN_MODEL_DATA.m_DataMap.End());
        itr != itr_end; ++itr
      )
      {
        // if the insert fails and it should be replaced
        if( !Insert( itr->first, itr->second) && REPLACE)
        {
          Replace( itr->first, itr->second);
        }
      }
    }

    //! @brief replace data in the map
    //! @param KEY the identifier for the data
    //! @param SP_DATA ShPtr to the ObjectInterface derived data
    //! @return true if data with that key existed and was replaced, false otherwise
    bool ProteinModelData::Replace( const Type &KEY, const util::ShPtr< util::ObjectInterface> &SP_DATA)
    {
      // check if key exists
      if( m_DataMap.Find( KEY) == m_DataMap.End())
      {
        return false;
      }

      // set the data
      m_DataMap[ KEY] = SP_DATA;

      // success
      return true;
    }

    //! @brief access the data with a key
    //! @param KEY the key for the data
    //! @return ShPtr to the data, undefined ShPtr returned if there is no data stored with the given key
    const util::ShPtr< util::ObjectInterface> &ProteinModelData::GetData( const Type &KEY) const
    {
      // find entry for KEY
      const storage::Map< TypeEnum, util::ShPtr< util::ObjectInterface> >::const_iterator itr( m_DataMap.Find( KEY));

      // check if key exists
      if( itr == m_DataMap.End())
      {
        return m_DataMap.Find( e_Undefined)->second;
      }

      // return the data
      return itr->second;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProteinModelData::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_DataMap, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProteinModelData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_DataMap, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace assemble
} // namespace bcl
