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

#ifndef BCL_ASSEMBLE_PROTEIN_MODEL_DATA_H_
#define BCL_ASSEMBLE_PROTEIN_MODEL_DATA_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinModelData
    //! @brief storage for protein model related data
    //! @details data can be a multiplier function, restraints, membrane, SSEPool etc... The data is stored in a map of
    //!          generic ShPtr to ObectInterface - so that caller of GetData has to cast it to the appropriate data type
    //!
    //! @see @link example_assemble_protein_model_data.cpp @endlink
    //! @author woetzen, karakam
    //! @date Nov 23, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinModelData :
      public util::ObjectInterface
    {

    public:

    ///////////
    // enums //
    ///////////

      //! types of protein model data
      enum Type
      {
        e_Pool,
        e_NativeModel,
        e_NativeFilteredModel,
        e_Sasa,
        e_Multiplier,
        e_LoopDomainLocators,
        e_PDBFile,
        e_Identification,
        e_Membrane,          //!< biol::Membrane object associated with that protein model
        e_SlicelistManager,
        e_Undefined,
        s_NumberTypes
      };

      //! @brief conversion to a string from a Type
      //! @param TYPE the type to get a string for
      //! @return a string representing that type
      static const std::string &GetTypeName( const Type &TYPE);

      //! @brief enum class wrapper for Type
      typedef util::WrapperEnum< Type, &GetTypeName, s_NumberTypes> TypeEnum;

    private:

      //! the map that contains the data with a string key
      storage::Map< TypeEnum, util::ShPtr< util::ObjectInterface> > m_DataMap;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ProteinModelData();

      //! @brief Clone function
      //! @return pointer to new ProteinModelData
      ProteinModelData *Clone() const;

      //! @brief hardcopy all member shared pointers
      //! @return pointer to hardcopied data
      util::ShPtr< ProteinModelData> HardCopy() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns whether the data map is empty
      //! @return whether the data map is empty
      bool IsEmpty() const
      {
        return m_DataMap.IsEmpty();
      }

      //! @brief gives the keys of the data indicating the data types that are contained in the data
      //! @return vector of data types that are contained in the data
      storage::Set< TypeEnum> GetKeys() const
      {
        return m_DataMap.GetKeys();
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief insert data to the map
      //! @param KEY the identifier for the data
      //! @param SP_DATA ShPtr to the ObjectInterface derived data
      //! @return true, if data with that key did not exist and was inserted, false otherwise
      bool Insert( const Type &KEY, const util::ShPtr< util::ObjectInterface> &SP_DATA);

      //! @brief insert data from another object into the map
      //! @param PROTEIN_MODEL_DATA object containing data to be added
      //! @param REPLACE if true, data of the same key is replaced
      void Insert( const ProteinModelData &PROTEIN_MODEL_DATA, const bool REPLACE = false);

      //! @brief replace data in the map
      //! @param KEY the identifier for the data
      //! @param SP_DATA ShPtr to the ObjectInterface derived data
      //! @return true if data with that key existed and was replaced, false otherwise
      bool Replace( const Type &KEY, const util::ShPtr< util::ObjectInterface> &SP_DATA);

      //! @brief access the data with a key
      //! @param KEY the key for the data
      //! @return ShPtr to the data, undefined ShPtr returned if there is no data stored with the given key
      const util::ShPtr< util::ObjectInterface> &GetData( const Type &KEY) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class ProteinModelData

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_PROTEIN_MODEL_DATA_H_
